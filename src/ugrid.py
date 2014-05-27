import numpy as np
"""
wrapper around unstructured netcdf output with various
helper functions.  mostly just a place to collect relevant
bits of code.
"""
import netCDF4
import trigrid
from matplotlib.dates import date2num
import datetime
import pytz


def cf_to_datenums(nc_t_var):
    """ parse the 'units since epoch' style of time axis
    to python datenums
    """ 
    units,origin = nc_t_var.units.split(" since ")
    try:
        origin_date = datetime.datetime.strptime(origin,'%Y-%m-%d %H:%M:%S %Z')
    except ValueError:
        origin_date = datetime.datetime.strptime(origin,'%Y-%m-%d %H:%M:%S')
        
    if origin_date.tzinfo is not None:
        origin_date = origin_date.astimezone(pytz.utc)
    else:
        # Not sure why this happens...
        print "Failed to get timezone - assuming netcdf time is UTC"
    tzero = date2num( origin_date )

    div = dict(seconds=86400.,
               minutes=60*24,
               hours=24,
               days=1)[units]
    return tzero + nc_t_var[:] / div

class Ugrid(object):
    surface_dzmin = 2*0.001 # common value, but no guarantee that this matches suntans code.
    
    def __init__(self,nc):
        """
        nc: path to netcdf dataset, or open dataset
        """

        if isinstance(nc,str):
            self.nc_filename = nc
            self.nc = netCDF4.Dataset(self.nc_filename,'r')
        else:
            self.nc_filename = None
            self.nc = nc

        # for operations which require a mesh to be specified, this is the default:
        self.mesh_name = self.mesh_names()[0]

    def data_variable_names(self,mesh_name):
        """ return list of variables which appear to have real data (i.e. not just mesh
        geometry / topology)
        """
        mesh_name = mesh_name or self.mesh_name
        
        data_names = []

        prefix = mesh_name+'_'
        for vname in self.nc.variables.keys():
            if vname.startswith(prefix):
                if self.nc.dimensions.has_key(vname):
                    continue
                if hasattr(self.nc.variables[vname],'cf_role'):
                    continue
                data_names.append( vname[len(prefix):] )
        return data_names
        
    def mesh_names(self):
        """
        Find the meshes in the file, based on cf_role == 'mesh_topology'
        """
        meshes = []
        for vname in self.nc.variables.keys():
            try:
                if self.nc.variables[vname].cf_role == 'mesh_topology':
                    meshes.append(vname)
            except AttributeError:
                pass
        return meshes

    _node_cache = None
    def get_node_array(self,node_coordinates):
        """ given a 'node_x_coordinates node_y_coordinates'
        attribute from a mesh_topology variable, return an [Nnodes,3]
        array for the node locations.  The z coordinate is set to 0.
        Note that this array is cached, so z should not be modified
        in the return value.
        """
        if self._node_cache is None:
            self._node_cache = {}

        if not self._node_cache.has_key(node_coordinates):
            node_x_name,node_y_name = node_coordinates.split()
            node_x = self.nc.variables[node_x_name][:]
            node_y = self.nc.variables[node_y_name][:]
            node_z = 0.0 * node_x

            self._node_cache[node_coordinates] = np.vstack( (node_x,node_y,node_z) ).T
        return self._node_cache[node_coordinates]
    
    def Ncells(self,mesh_name=None):
        mesh_name = mesh_name or self.mesh_name
        mesh = self.nc.variables[mesh_name]
        # FIX: this makes a big assumption on the order of dimensions!
        return len( self.nc.dimensions[self.nc.variables[mesh.face_node_connectivity].dimensions[0]] )

    def Nkmax(self,mesh_name=None):
        mesh_name = mesh_name or self.mesh_name

        for dim_name in 'n%s_layers'%mesh_name, 'nMeshGlobal_layers':
            try:
                return len( self.nc.dimensions[dim_name] )
            except KeyError:
                pass
        raise Exception,"Failed to find vertical dimension"

    def get_cell_velocity(self,time_step,mesh_name=None):
        mesh_name = mesh_name or self.mesh_name

        U = np.zeros( (self.Ncells(mesh_name),self.Nkmax(),2), np.float64)
        U[:,:,0] = self.nc.variables[mesh_name + '_cell_east_velocity'][:,:,time_step]
        U[:,:,1] = self.nc.variables[mesh_name + '_cell_north_velocity'][:,:,time_step]
        return U

    def get_cell_scalar(self,label,time_step,mesh_name=None):
        mesh_name = mesh_name or self.mesh_name
        
        # totally untested!
        return self.nc.variables[mesh_name + '_' + label][:,:,time_step]
    
    def to_trigrid(self,mesh_name=None,skip_edges=False):
        nc = self.nc
        mesh_name = mesh_name or self.mesh_name
        mesh = nc.variables[mesh_name]

        node_x_name,node_y_name = mesh.node_coordinates.split()
        
        node_xy = np.array( [nc.variables[node_x_name][...],
                             nc.variables[node_y_name][...]]).T
        faces = nc.variables[mesh.face_node_connectivity][...]
        edges = nc.variables[mesh.edge_node_connectivity][...] # [N,2]
        g = trigrid.TriGrid(points=node_xy,cells=faces,edges=edges)
        
        if not skip_edges:
            # g.make_edges_from_cells() # this completely recreates the edges
            # instead, we need to discern the edge-cell connectivity from the
            # supplied edges
            pass
        return g

    def vertical_averaging_weights(self,time_step,ztop=None,zbottom=None,dz=None,mesh_name=None):
        """
        reimplementation of sunreader.Sunreader::averaging_weights
        
        Returns weights as array [Nk] to average over a cell-centered quantity
        for the range specified by ztop,zbottom, and dz.

        range is specified by 2 of the 3 of ztop, zbottom, dz, all non-negative.
        ztop: distance from freesurface
        zbottom: distance from bed
        dz: thickness

        if the result would be an empty region, return nans.

        this thing is slow! - lots of time in adjusting all_dz
        """
        mesh_name = mesh_name or self.mesh_name
        
        mesh = self.nc.variables[mesh_name]
        layer_var = 'nMeshGlobal_layers' # default
        for name in self.nc.dimensions.keys(): # but also try looking for it.
            try:
                if self.nc.variables[name].standard_name == 'ocean_zlevel_coordinate':
                    layer_var = name
                    break
            except KeyError:
                pass
            except AttributeError:
                pass
        # 
        layers = self.nc.variables[layer_var]
        layer_bounds = self.nc.variables[ layers.bounds ][...]

        layer_interfaces = np.concatenate( (layer_bounds[:,0],layer_bounds[-1:,1]) )

        bed = -self.nc.variables['Mesh_depth'][...]

        # FIX: this assumes an order to the dimensions...
        Ncells = self.Ncells()

        # this is a bit trickier, because there could be lumping.  for now, it should work okay
        # with 2-d, but won't be good for 3-d
        Nk = np.searchsorted(-layer_interfaces,-bed)

        one_dz = -np.diff(layer_interfaces)
        all_dz = one_dz[None,:].repeat(Ncells,axis=0)
        all_k = np.arange(len(one_dz))[None,:].repeat(Ncells,axis=0)

        z = layer_bounds.min(axis=1) # bottom of each cell

        h = self.nc.variables[mesh_name+'_surface'][:,time_step]

        # adjust bed and 
        # 3 choices here..
        # try to clip to reasonable values at the same time:
        if ztop is not None:
            if ztop != 0:
                h = h - ztop # don't modify h
                # don't allow h to go below the bed
                h[ h<bed ] = bed
            if dz is not None:
                # don't allow bed to be below the real bed.
                bed = np.maximum( h - dz, bed)
        if zbottom is not None:
            # no clipping checks for zbottom yet.
            if zbottom != 0:
                bed = bed + zbottom # don't modify bed!
            if dz is not None:
                h = bed + dz

        # so now h and bed are elevations bounding the integration region

        ctops = np.searchsorted(-z - self.surface_dzmin, -h)

        # default h_to_ctop will use the dzmin appropriate for the surface,
        # but at the bed, it goes the other way - safest just to say dzmin=0,
        # and also clamp to known Nk
        cbeds = np.searchsorted(-z,-bed) + 1 # it's an exclusive index
        cbeds[ cbeds>Nk ] = Nk[ cbeds>Nk ]

        # seems that there is a problem with how dry cells are handled -
        # for the exploratorium display this ending up with a number of cells with
        # salinity close to 1e6.
        # in the case of a dry cell, ctop==cbed==Nk[i]
        drymask = (all_k < ctops[:,None]) | (all_k>=cbeds[:,None])
        all_dz[drymask] = 0.0
        ii = np.arange(Ncells)
        all_dz[ii,ctops] = h - z[ctops]
        all_dz[ii,cbeds-1] -= bed - z[cbeds-1]


        # make those weighted averages
        # have to add extra axis to get broadcasting correct
        all_dz = all_dz / np.sum(all_dz,axis=1)[:,None]
        return all_dz

    def datenums(self):
        """ return datenums, referenced to UTC
        """
        return cf_to_datenums(self.nc.variables['time'])


