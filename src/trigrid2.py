# trigrid2.py: a temporary side branch of trigrid.py
#  see about having structured arrays for cells,points,edges
#  so that subclasses can add fields which are kept in sync with
#  the basic elements.
#  this also introduces the difference between points - which are x,y coordinate pairs,
#    and nodes, which are parts of the graph (and have a point, but also other data)

# To sort out:
# closest_point => closest_node

import cPickle as pickle
# import pickle

try:
    from safe_pylab import *
    from matplotlib.collections import *
    import plot_wkb
except:
    pass

from numpy import *
from numpy.linalg import norm

import wkb2shp,field

import sys,time

# qgis clashes with the Rtree library (because it includes its own local copy).
# fall back to a wrapper around the qgis spatial index if it looks like we're running
# under qgis.
from safe_rtree import Rtree

try:
    from shapely import geometry,wkb
    import shapely.predicates
except ImportError:
    print "Shapely is not available!"
    geometry = "unavailable"
import priority_queue as pq
import code
import os,types

import pdb

try:
    try:
        from osgeo import ogr,osr
    except ImportError:
        import ogr,osr
except ImportError:
    print "GDAL failed to load"
    ogr = "unavailable"
    osr = ogr
    
from array_append import array_append,array_concatenate,concatenate_safe_dtypes


# edge markers:
CUT_EDGE = 37 # the marker for a cut edge
OPEN_EDGE = 3
LAND_EDGE = 1
DELETED_EDGE = -1

# edge-cell markers ( the cell ids that reside in the edge array
BOUNDARY = -1 # cell marker for edge of domain
UNMESHED = -2 # cell marker for edges not yet meshed

xxyy = array([0,0,1,1])
xyxy = array([0,1,0,1])

# element statuses -
STAT_OK = 0 # b/c zero is the default
STAT_DEL = 1

# rotate the given vectors/points through the CCW angle in radians
def rot(angle,pnts):
    R = array( [[cos(angle),-sin(angle)],
                [sin(angle),cos(angle)]] )
    return tensordot(R,pnts,axes=(1,-1) ).transpose() # may have to tweak this for multiple points

def signed_area(points):
    i = arange(points.shape[0])
    ip1 = (i+1)%(points.shape[0])
    return 0.5*(points[i,0]*points[ip1,1] - points[ip1,0]*points[i,1]).sum()
    
def is_ccw(points):
    return signed_area(points) > 0
    
def ensure_ccw(points):
    if not is_ccw(points):
        # print "Hey - you gave me CW points.  I will reverse"
        points = points[::-1]
    return points

def ensure_cw(points):
    if is_ccw(points):
        # print "Hey - you gave me CCW points.  I will reverse"
        points = points[::-1]
    return points


def outermost_rings( poly_list ):
    """ given a list of Polygons, return indices for those that are not inside
    any other polygon
    """
    areas = array( [p.area for p in poly_list])
    order = argsort(-1 * areas) # large to small
    outer = []

    for i in range(len(order)):
        ri = order[i]
        # print "Checking to see if poly %d is an outer polygon"%ri

        is_exterior = 1
        # check polygon ri (the ith largest) against all polygons
        # larger than it.
        for j in range(i):
            rj = order[j]

            if poly_list[rj].contains( poly_list[ri] ):
                # print "%d contains %d"%(rj,ri)
                is_exterior = 0 # ri is contained by rj, so not exterior
                break

        if is_exterior:
            # print "%d is exterior"%ri
            outer.append(ri)
    return outer


def circumcenter(p1,p2,p3):
    ref = p1
    
    p1x = p1[...,0] - ref[...,0] # ==0.0
    p1y = p1[...,1] - ref[...,1] # ==0.0
    p2x = p2[...,0] - ref[...,0]
    p2y = p2[...,1] - ref[...,1]
    p3x = p3[...,0] - ref[...,0]
    p3y = p3[...,1] - ref[...,1]

    vc = zeros( p1.shape, float64)
    
    # taken from TRANSFORMER_gang.f90
    dd=2.0*((p1x-p2x)*(p1y-p3y) -(p1x-p3x)*(p1y-p2y))
    b1=p1x**2+p1y**2-p2x**2-p2y**2
    b2=p1x**2+p1y**2-p3x**2-p3y**2 
    vc[...,0]=(b1*(p1y-p3y)-b2*(p1y-p2y))/dd + ref[...,0]
    vc[...,1]=(b2*(p1x-p2x)-b1*(p1x-p3x))/dd + ref[...,1]
    
    return vc

def dist(p1,p2):
    return sqrt( sum((p1-p2)**2) )


class TriGridError(Exception):
    pass

class NoSuchEdgeError(TriGridError):
    pass

class NoSuchCellError(TriGridError):
    pass


# cache the results of reading points.dat files for suntans grid files
# maps filenames to point arrays - you should probably
# just copy the array, though, since there is the possibility
# of altering the points array
points_dat_cache = {}


class TriGrid(object):
    index = None
    edge_index = None
    _vcenters = None
    verbose = 0
    default_clip = None

    # Define the data stored for each point
    node_dtype = [ ('x',(float64,2)),
                   ('stat',int8) ]
    cell_dtype  = [ ('nodes',(int32,3)), # node indexes
                    ('vor_center',(float64,2)),  # voronoi center
                    ('centroid',(float64,2)),
                    ('stat',int8)]
    edge_dtype = [ ('nodes',(int32,2)),
                   ('marker',int32),
                   ('cells',(int32,2)),
                   ('center',(float64,2)),
                   ('stat',int8)]

    # CONSTANTS:
    CUT_EDGE = 37 # the marker for a cut edge
    OPEN_EDGE = 3
    LAND_EDGE = 1
    DELETED_EDGE = -1

    # edge-cell markers ( the cell ids that reside in the edge array
    BOUNDARY = -1 # cell marker for edge of domain
    UNMESHED = -2 # cell marker for edges not yet meshed

    # element statuses -
    STAT_OK = 0 # b/c zero is the default
    STAT_DEL = 1

    def __init__(self,sms_fname=None,
                 tri_basename=None,
                 suntans_path=None,processor=None,
                 suntans_reader=None,
                 tec_file=None,
                 gmsh_basename=None,
                 edges=None,points=None,cells=None,
                 extra_node_fields=[],
                 extra_cell_fields=[],
                 extra_edge_fields=[]):
        # Do this early, so any early allocations use the right set of fields
        # but be careful not to modify the class-wide defs - just the instance's.
        self.node_dtype = self.node_dtype + extra_node_fields
        self.cell_dtype = self.cell_dtype + extra_cell_fields
        self.edge_dtype = self.edge_dtype + extra_edge_fields
        
        self.sunreader = None
        self._node_to_cells = None

        self.init_listeners()
        
        if sms_fname:
            self.read_sms(sms_fname)
        elif tri_basename:
            self.read_triangle(tri_basename)
        elif gmsh_basename:
            self.read_gmsh(gmsh_basename)
        elif suntans_path:
            self.processor = processor
            self.suntans_path = suntans_path
            self.read_suntans()
        elif suntans_reader:
            self.processor = processor
            self.sunreader = suntans_reader
            self.read_suntans()
        elif tec_file:
            self.read_tecplot(tec_file)
        elif points is not None:
            self.from_data(points,edges,cells)
        else:
            # This will create zero-length arrays for everyone.
            self.from_data(None,None,None)

    def file_path(self,conf_name):
        if self.sunreader:
            return self.sunreader.file_path(conf_name,self.processor)
        else:
            if conf_name == 'points':
                basename = 'points.dat'
            elif conf_name == 'edges':
                basename = 'edges.dat'
            elif conf_name == 'cells':
                basename = 'cells.dat'
            else:
                raise Exception,"Unknown grid conf. name: "+conf_name
            if self.processor is not None and conf_name != 'points':
                basename = basename + ".%i"%self.processor
            return self.suntans_path + '/' + basename
    def from_data(self,points,edges,cells):
        """ points: either None, or an [N,2] array of xy pairs.
        """
        if points is None:
            self.nodes = zeros( 0, self.node_dtype )
        else:
            self.nodes = zeros( len(points), self.node_dtype)
            self.nodes['x'] = points[:,:2]  # discard any z's that come in
        
        if edges is None:
            self.edges = zeros(0,self.edge_dtype)
        else:
            ne = len(edges)
            self.edges = zeros(ne,self.edge_dtype)
            
            # incoming edges may just have connectivity
            if edges.shape[1] == 2:
                self.edges['nodes'] = edges

                # best guess at the edge markers
                self.edges['marker'] = LAND_EDGE 
                self.edges['cells'][:,0] = UNMESHED
                self.edges['cells'][:,1] = BOUNDARY
            elif edges.shape[1] == 5:
                self.edges['nodes'] = edges[:,:2]
                self.edges['marker'] = edges[:,2]
                self.edges['cells'] = edges[:,3:5]
            else:
                raise Exception,"Edges should have 2 or 5 entries per edge"

        if cells is None:
            self.cells = zeros( 0, self.cell_dtype)
        else:
            self.cells = zeros( len(cells), self.cell_dtype)
            if len(cells)>0:
                self.cells['nodes'] = cells
                self.cells['centroid'] = nan
                self.cells['vor_center'] = nan

    def copy(self,types_only=False):
        """ Return a new instance which is a deep copy of self
        If types_only is true, then only the datatypes are copied,
        not the actual data.  Useful when making a partial copy of a grid
        """
        g = TriGrid()
        g.node_dtype = self.node_dtype
        g.cell_dtype = self.cell_dtype
        g.edge_dtype = self.edge_dtype

        if types_only:
            g.from_data(None,None,None)
        else:
            g.nodes = self.nodes.copy()
            g.edges = self.edges.copy()
            g.cells = self.cells.copy()

        return g

    def extra_node_fields(self):
        return [nametype[0] for nametype in self.node_dtype[ len(TriGrid.node_dtype): ]]

    def extra_edge_fields(self):
        return [nametype[0] for nametype in self.edge_dtype[ len(TriGrid.edge_dtype): ]]

    def extra_cell_fields(self):
        return [nametype[0] for nametype in self.cell_dtype[ len(TriGrid.cell_dtype): ]]
        
    def refresh_metadata(self):
        """ Call this when the cells, edges and nodes may be out of sync with indices
        and the like.
        """
        self.index = None
        self.edge_index = None
        self._calc_edge_centers = False
        self._calc_cell_centers = False
        self._calc_vcenters = False


    def read_suntans(self,use_cache=1):
        self.read_from = "Suntans"

        # read the points:
        points_fn = os.path.abspath( self.file_path("points") )
        if use_cache and points_dat_cache.has_key(points_fn):
            points = points_dat_cache[points_fn]
        else:
            points_fp = open(points_fn)

            points = []
            for line in points_fp:
                coords = map(float,line.split())
                if len(coords) >= 2:
                    points.append(coords[:2])
            points = array(points)
            if use_cache:
                points_dat_cache[points_fn] = points

        # read the cells:
        cell_fname = self.file_path("cells")
            
        cells_fp = open(cell_fname)

        cells = []
        for line in cells_fp:
            line = line.split()
            if len(line)==8:
                # then three point indices:
                cells.append( map(int,line[2:5]) )

        cells = array(cells)
        
        # Edges!
        # Each line is endpoint_i endpoint_i  edge_marker  cell_i cell_i
        edge_fname = self.file_path('edges')

        # edges are stored just as in the data file:
        #  point_i, point_i, marker, cell_i, cell_i
        edges_fp = open(edge_fname,"rt")

        edges = []
        for line in edges_fp:
            line = line.split()
            if len(line) == 5:
                edges.append(map(int,line))
        edges = array(edges)
        
        self.from_data(points,edges,cells)
        
    def read_gmsh(self,gmsh_basename):
        """ reads output from gmsh  - gmsh_basename.{nod,ele}
        """
        self.fname = gmsh_basename
        self.read_from = "GMSH"
        
        points = loadtxt( self.fname +".nod")
        id_offset = int(points[0,0]) # probably one-based

        points = points[:,1:3]
        
        print "Reading cells"
        elements = loadtxt( self.fname +".ele")

        cells = elements[:,1:4].astype(int32) - id_offset

        self.from_data(points=points,
                       edges = None,
                       cells=cells)
        self.make_edges_from_cells()

    def read_triangle(self,tri_basename):
        """ reads output from triangle, tri_basename.{ele,poly,node}
        """
        self.fname = tri_basename
        self.read_from = "Triangle"
        
        points_fp = open(tri_basename + ".node")

        Npoints,point_dimension,npoint_attrs,npoint_markers = map(int,points_fp.readline().split())

        points = zeros((Npoints,2),float64)
        id_offset = 0
        for i in range(self.Npoints()):
            line = points_fp.readline().split()
            # pnt_id may be 0-based or 1-based.
            pnt_id = int(line[0])
            if i == 0:
                id_offset = pnt_id
                print "Index offset is ",id_offset

            points[i,:] = map(float,line[1:3])
        points_fp.close()

        print "Reading cells"
        elements_fp = open(tri_basename + ".ele")
        Ncells,node_per_tri,ncell_attrs = map(int,elements_fp.readline().split())

        if node_per_tri != 3:
            raise Exception,"Please - just use 3-point triangles!"
        
        cells = zeros((Ncells,3),int32) 

        pnt2cells = {} # hash for figuring out edge neighbors below
        
        for i in range(Ncells):
            parsed = map(int,elements_fp.readline().split())
            cell_id = parsed[0]
            cells[i,:] = array(parsed[1:]) - id_offset
            # update the mapping so we can fill in edge neighbors below:
            for n in cells[i,:]:
                pnt2cells[n] = i

        edges_fn = tri_basename + ".edge"
        if os.path.exists(edges_fn):
            edges_fp = open()

            Nedges,nedge_markers = map(int,edges_fp.readline().split())

            edges = zeros((Nedges,5),int32)

            # each edge is stored as:  (pnt_a, pnt_b, default_marker,node_1,node_2)
            for i in range(Nedges):
                idx,pnta,pntb,marker = map(int,edges_fp.readline().split())
                # and a bit of work to figure out which cells border this edge:
                pnta -= id_offset
                pntb -= id_offset
                cells_a = pnt2cells[pnta]
                cells_b = pnt2cells[pntb]

                adj_cells = intersect1d(cells_a,cells_b)
                neighbor1 = adj_cells[0]

                if len(adj_cells) == 1:
                    neighbor2 = -1
                else:
                    neighbor2 = adj_cells[1]

                edges[i] = [pnta,pntb,marker,neighbor1,neighbor2]
            print "WARNING: triangle reading code doesn't respect which side of the edge each cell is on"
            self.from_data(points=points,
                           cells=cells,
                           edges=edges)
        else:
            print "No edges - will recreate from cells"
            self.from_data(points=points,
                           cells=cells,
                           edges=None)
            self.make_edges_from_cells()
            
    def read_tecplot(self,fname):
        self.read_from = 'tecplot'
        self.fname = fname

        fp = open(fname)

        while 1:
            line = fp.readline()
            if line.find('ZONE') == 0:
                break

        import re
        m = re.search(r'\s+N=\s*(\d+)\s+E=\s*(\d+)\s',line)
        if not m:
            print "Failed to parse: "
            print line
            raise Exception,"Tecplot parsing error"

        
        # first non-blank line has number of cells and edges:
        Ncells = int( m.group(2) )
        Npoints = int( m.group(1) )

        points = zeros((Npoints,2),float64)
        for i in range(Npoints):
            points[i,:] = map(float,fp.readline().split())

        print "Reading cells"
        cells = zeros((Ncells,3),int32) # store zero-based indices

        # we might be reading in the output from ortho, in which
        # it reports the number of unique cells, not the real number
        # cells

        i=0
        cell_hash = {}
        for line in fp:
            pnt_ids = array( map(int,line.split()) )

            my_key = tuple(sort(pnt_ids))
            if not cell_hash.has_key(my_key):
                cell_hash[my_key] = i
                
                # store them as zero-based
                cells[i] = pnt_ids - 1
                i += 1

        if i != Ncells:
            print "Reading %i cells, but expected to get %i"%(i,self.Ncells)
            cells = self.cells[:i,:]

        self.from_data(points=points,cells=cells)
        
        # At this point we have enough info to create the edges
        self.make_edges_from_cells()

    # these are used in some gui code
    _calc_cell_centers = False # whether cell centers should be kept up to date
    def cell_centers(self):
        if not self._calc_cell_centers:
            self._calc_cell_centers = True
            self.cells['centroid'] = self.nodes['x'][self.cells['nodes']].mean(axis=1)
        return self.cells['centroid']
    
    _calc_edge_centers = False
    def edge_centers(self):
        if not self._calc_edge_centers:
            self._calc_edge_centers = True
            self.edges['center'] = self.nodes['x'][self.edges['nodes']].mean(axis=1)

        return self.edges['center']
    
    def update_edge_center(self,j):
        self.edges[j]['center'][:] = self.nodes['x'][self.edges[j]['nodes']].mean(axis=0)

    def ghost_cells(self):
        """ Return a bool array, with ghost cells marked True
        Ghost cells are determined as any cell with an edge that has marker 6
        """
        ghost_edges = self.edges['marker'] == 6
        ghost_cells = self.edges[ghost_edges]['cells'].ravel()

        bitmap = zeros( self.Ncells(), bool8 )
        bitmap[ ghost_cells ] = True
        return bitmap

    def delete_unused_nodes(self):
        """ any nodes which aren't in any cells or edges will be removed.
        """
        all_nodes = arange(self.Npoints())

        cell_nodes = unique(ravel(self.cells['nodes']))
        edge_nodes = unique(ravel(self.edges['nodes']))
        deleted_nodes = nonzero( self.nodes['stat']==STAT_DEL)[0]

        okay_nodes = unique( concatenate( (cell_nodes,edge_nodes,deleted_nodes) ) )

        unused = setdiff1d(all_nodes,okay_nodes)

        for n in unused:
            self.delete_node(n)
        
    def renumber(self):
        """
        removes duplicate cells, edges, and nodes that are not
        referenced by any cell, as well as cells that have been deleted (==-1)
        """
        ### NODES ### 
        # remove lonesome nodes - based on the edges
        valid_edges = (self.edges['stat'] == STAT_OK)
        
        active_nodes = unique(ravel(self.edges['nodes'][valid_edges]))
        
        if len(active_nodes)>0 and any(active_nodes<0):
            raise Exception,"renumber: Active nodes includes some negative indices"
        nodes_old_to_new = -ones(self.Npoints(),int32)

        if len(active_nodes) > 0:
            self.nodes = self.nodes[active_nodes]
            # need a mapping from active node to its index -
            # explicitly ask for int32 for consistency
            nodes_old_to_new[active_nodes] = arange(active_nodes.shape[0],dtype=int32)
        else:
            self.nodes = zeros(0,self.node_dtype)
            
        if any(isnan(self.nodes['x'])):
            raise Exception,"renumber: some points have NaNs!"

        ### CELLS ###
        cell_hash = {} # sorted tuples of vertices
        active_cells = [] # list of indexes into the old ones - active_cells[new_index] = old_index
        for i in range(self.Ncells()):
            my_key = tuple(sort(self.cells[i]['nodes']) )

            if not cell_hash.has_key(my_key) and self.cells[i]['stat']==STAT_OK:
                # we're original and not deleted
                cell_hash[my_key] = i # value is ignored...
                active_cells.append( i )
        if len(active_cells)>0:
            self.cells = self.cells[active_cells]
            # map node indices onto the new indices
            self.cells['nodes'] = nodes_old_to_new[self.cells['nodes']]
            
            # i.e. old_cells[old_index] = new_index
            cells_old_to_new = -ones(self.Ncells(),int32)
            cells_old_to_new[active_cells] = arange(len(active_cells))
        else:
            self.cells = zeros( 0, self.cell_dtype)
            cells_old_to_new = None

        if len(self.cells)>0 and any(self.cells['nodes']) < 0:
            raise Exception,"renumber: after remapping indices, have negative node index in cells"

        ### EDGES ###
        edge_hash = {} # sorted tuples of vertices
        active_edges = [] # list of indexes into the old ones
        for j in range(self.Nedges()):
            my_key = tuple(sort(self.edges['nodes'][j]) )

            if not edge_hash.has_key(my_key) and self.edges['stat'][j]==STAT_OK:
                # we're original and not deleted
                edge_hash[my_key] = j # value is ignored...
                active_edges.append( j )

        if len(active_edges) > 0:
            self.edges = self.edges[active_edges]
            # Map node and cell indices:
            self.edges['nodes'] = nodes_old_to_new[self.edges['nodes']]
            real_cells = (self.edges['cells'] >= 0)
            if any(real_cells):
                self.edges['cells'][real_cells] = cells_old_to_new[self.edges['cells'][real_cells]]
            # let's -1 and -2 edges stay as they are...
        else:
            self.edges = zeros(0,self.edge_dtype)
            

        # clear out stale data
        self._node_to_cells = None
        self.index = None
        self.edge_index = None
        self._node_to_edges = None


        # return the mappings so that subclasses can catch up
        return {'cell_map':cells_old_to_new,
                'node_map':nodes_old_to_new,
                'valid_nodes':active_nodes}

    def write_Triangle(self,basename,boundary_nodes=None):
        """  duplicate some of the output of the Triangle program -
        particularly the .node and .ele files

        note that node and cell indices are taken as 1-based.

        if boundary_nodes is supplied, it should be an integer valued array of length Npoints,
        and give the boundary marker for each node (usually 0 for internal, nonzero for boundary).
        this can be used to specify a subset of the boundary nodes for a BC in SWAN.

        if not specified, boundary markers will be 0 for internal, 1 for external nodes.
        """
        node_fp = open(basename + ".node",'wt')
        node_fp.write("%d 2 0 1\n"%(self.Npoints()))
        for n in range(self.Npoints()):
            if boundary_nodes is not None:
                bmark = boundary_nodes[n]
            else:
                # id x y boundary marker
                bmark = 0
                if self.boundary_angle(n) != 0:
                    bmark = 1
            node_fp.write("%d %f %f %d\n"%(n+1,self.nodes['x'][n,0],self.nodes['x'][n,1], bmark ) )
        node_fp.close()

        ele_fp = open(basename + ".ele",'wt')
        ele_fp.write("%d 3 0\n"%(self.Ncells()))
        for i in range(self.Ncells()):
            ele_fp.write("%d %d %d %d\n"%(i+1,
                                          self.cells['nodes'][i,0]+1,
                                          self.cells['nodes'][i,1]+1,
                                          self.cells['nodes'][i,2]+1))
        ele_fp.close()

    def write_obj(self,fname):
        """ Output to alias wavefront
         - scales points to fall within [0,10]
        """
        fp = open(fname,'wt')

        pmax = self.nodes['x'].max(axis=0)
        pmin = self.nodes['x'].min(axis=0)
        rng = (pmax-pmin).max()

        scaled_points = (self.nodes['x'] - pmin)*(10/rng)

        for i in range(self.Npoints()):
            fp.write("v %f %f 0.0\n"%(scaled_points[i,0],scaled_points[i,1]))


        for i in range(self.Ncells()):
            fp.write("f %d %d %d\n"%(self.cells['nodes'][i,0]+1,
                                     self.cells['nodes'][i,1]+1,
                                     self.cells['nodes'][i,2]+1))

        fp.close()

        
    def write_tulip(self,fname):
        """ Write a basic representation of the grid to a tulip
        compatible file
        """
        fp = file(fname,'wt')

        fp.write("(tlp \"2.0\"\n")
        
        fp.write("(nodes ")

        for i in range(self.Npoints()):
            if self.nodes['stat'][i]==STAT_OK:
                fp.write(" %i"%i )

        fp.write(")\n")

        for e in range(self.Nedges()):
            if self.edges[e]['stat']==STAT_OK:
                fp.write("(edge %i %i %i)\n"%(e,self.edges['nodes'][e,0],self.edges['nodes'][e,1]))

        # and the locations of the nodes
        fp.write("(property 0 layout \"viewLayout\" \n")
        for i in range(self.Npoints()):
            if self.nodes['stat'][i]==STAT_OK:
                fp.write("  (node %i \"(%f,%f,0)\")\n"%(i,self.nodes['x'][i,0],self.nodes['x'][i,1]))

        fp.write(")\n")
        
        fp.write(")\n")
        fp.close()
        
    def write_sms(self,fname):
        fp = open(fname,'wt')

        fp.write("\n") # seems to start with blank line.
        fp.write("%i  %i\n"%(self.Ncells(),self.Npoints()))
        
        # each point has three numbers, though the third is apparently
        # always 0

        for i in range(self.Npoints()):
            fp.write("%10i %.11f %.11f %.11f\n"%(i+1,
                                                 self.nodes['x'][i,0],
                                                 self.nodes['x'][i,1],
                                                 0.0 ))
            
        # everything is a triangle

        # compute area, positive means CCW
        #  - turns out SMS wants the order to be consistent, but it always *creates* CCW
        # triangles.  so best to create CCW triangles
        bad = self.areas() < 0
        n_bad = sum(bad)
        
        if n_bad > 0:
            print "Found %i CW triangles that will be reversed"%n_bad
            self.cells['nodes'][bad,:] = self.cells['nodes'][bad,::-1]

        for i in range(self.Ncells()):
            fp.write("%i 3 %i %i %i\n"%(i+1,
                                        self.cells['nodes'][i,0]+1,
                                        self.cells['nodes'][i,1]+1,
                                        self.cells['nodes'][i,2]+1) )

        # And then go back and switch the marker for some of the edges:
        print "SMS output: omitting boundary information"
        fp.write("0 = Number of open boundaries\n")
        fp.write("0 = Total number of open boundary nodes\n")
        fp.write("0 = Number of land boundaries\n")
        fp.write("0 = Total number of land boundary nodes\n")
        fp.close()

    def areas(self):
        """ returns signed area, CCW is positive"""
        i = array([0,1,2])
        ip = array([1,2,0])
        xi =  self.nodes['x'][self.cells['nodes'][:,i],0]
        yi =  self.nodes['x'][self.cells['nodes'][:,i],1]
        xip = self.nodes['x'][self.cells['nodes'][:,ip],0]
        yip = self.nodes['x'][self.cells['nodes'][:,ip],1]

        A = 0.5 * (xi*yip-xip*yi).sum(axis=1)

        return A

    def read_sms(self,fname):
        self.fname = fname
        self.read_from = "SMS"

        fp = open(fname)

        # skip leading blank lines
        while 1:
            line = fp.readline().strip()
            if line != "":
                break

        # first non-blank line has number of cells and edges:
        Ncells,Npoints = map(int,line.split())

        # each point has three numbers, though the third is apparently
        # always 0

        points = zeros((Npoints,2),float64)
        for i in range(Npoints):
            line = fp.readline().split()
            # pnt_id is 1-based
            pnt_id = int(line[0])

            points[pnt_id-1] = map(float,line[1:3])

        print "Reading cells"
        # store zero-based indices, and assume
        # everything is a triangle
        cells = zeros((Ncells,3),int32)
        
        for i in range(Ncells):
            parsed = map(int,fp.readline().split())
            cell_id = parsed[0]
            nvertices = parsed[1]
            pnt_ids = array(parsed[2:])

            if nvertices != 3:
                raise "Assumption of all triangles is not true!"
            # store them as zero-based
            self.cells['nodes'][cell_id-1] = pnt_ids - 1

        self.from_data(points=points,
                       cells=cells,
                       edges=None)
                       
        # At this point we have enough info to create the edges
        self.make_edges_from_cells()

        # And then go back and switch the marker for some of the edges:
        print "Reading boundaries"
        def read_first_int():
            return int(fp.readline().split()[0])

        for btype in ['open','land']:
            if btype == 'open':
                marker = 3 # open - not sure if this is 2 or 3...
            else:
                marker = 1 # closed
                
            n_boundaries = read_first_int()
            print "Number of %s boundaries: %d"%(btype,n_boundaries)
            tot_boundary_nodes = read_first_int() # who cares...

            for boundary_i in range(n_boundaries):
                print "Reading %s boundary %d"%(btype,boundary_i+1)
                n_nodes_this_boundary = read_first_int()
                for i in range(n_nodes_this_boundary):
                    node_i = read_first_int() - 1 # zero-based
                    if i>0:
                        # update the marker in edges
                        if node_i < last_node_i:
                            pa,pb = node_i,last_node_i
                        else:
                            pa,pb = last_node_i,node_i
                            
                        try:
                            edge_i = self.nodes_to_edge((pa,pb))
                            self.edges['marker'][edge_i] = marker
                        except NoSuchEdgeError:
                            print "Couldn't find edge",(pa,pb)
                            print self.nodes[ [pa,pb] ]
                            raise
                            
                    last_node_i = node_i

        print "Done"


    def pnt2cells(self,pnt_i):
        raise Exception,"pnt2cells is deprecated - call node_to_cells()"
        
    def node_to_cells(self,node_i):
        """ Given a node index, return a list of the cells that use that node.
        """
        if self._node_to_cells is None:
            # build hash table for point->cell lookup
            self._node_to_cells = {}
            for i in range(self.Ncells()):
                for n in self.cells['nodes'][i]:
                    if not self._node_to_cells.has_key(n):
                        self._node_to_cells[n] = []
                    self._node_to_cells[n].append(i)
        return self._node_to_cells[node_i]

    def Nedges(self):
        return len(self.edges)
    def Ncells(self):
        return len(self.cells)
    def Npoints(self):
        return self.Nnodes()
    def Nnodes(self):
        return len(self.nodes)
    
    def pnt2edges(self,pnt_i):
        raise Exception,"pnt2edges is deprecated - use node_to_edges"
    
    _node_to_edges = None
    def node_to_edges(self,node_i):
        if self._node_to_edges is None:
            n2e = {}
            for e in range(self.Nedges()):
                # skip deleted edges
                if self.edges['stat'][e] != STAT_OK:
                    continue
                
                for n in self.edges[e]['nodes']:
                    if not n2e.has_key(n):
                        n2e[n] = []
                    n2e[n].append(e)
            self._node_to_edges = n2e

        if self._node_to_edges.has_key(node_i):
            return self._node_to_edges[node_i]
        else:
            return []

    def boundary_angle(self,pnt_i):
        """ returns the interior angle in radians, formed by the
        boundary at the given point
        """
        edges = self.node_to_edges(pnt_i)

        # find the absolute angle of each edge, as an angle CCW from
        # east

        angle_right=None # the angle of the edge with the domain on the right
        angle_left =None # angle of the edge with the domain on the left

        for edge in edges:
            # only care about the edges on the boundary:
            if self.edges['cells'][edge,1] != BOUNDARY:
                continue
            segment = self.edges[edge]['nodes']
            seg_reversed = 0
            if segment[0] != pnt_i:
                segment = segment[::-1]
                seg_reversed = 1
                
            # sanity check
            if segment[0] != pnt_i:
                raise "Well, where is %d in %s"%(pnt_i,segment)

            delta = self.nodes['x'][segment[1]] - self.nodes['x'][segment[0]]

            angle = arctan2(delta[1],delta[0])
            # print "Edge %i to %i has angle %g degrees"%(edge,segment[1],180*angle/pi)

            # on which side of this edge is the domain?
            my_cell = self.edges['cells'][edge,0]

            if my_cell == UNMESHED:
                # the paver enforces that cell markers are 3=>left,4=>right
                # so with the stored order of the edge, the pretend cell center
                # is always to the left
                if not seg_reversed:
                    xprod = -1
                else:
                    xprod = 1
            else:
                my_cell_middle = mean( self.nodes['x'][ self.cells['nodes'][my_cell] ] , axis=0 )

                delta_middle = my_cell_middle - self.nodes['x'][pnt_i]

                # and cross-product:
                xprod = cross(delta_middle,delta)
                
            # print "Cross-product is: ",xprod
            if xprod > 0:
                # the cell center lies to the right of this edge,
                # print "Edge to %i has domain to the right"%segment[1]
                angle_right = angle
            else:
                # print "Edge to %i has domain to the left"%segment[1]
                angle_left = angle
        if angle_left is None and angle_right is None:
            # it's an interior node, so no boundary angle...
            return 0.0

        if angle_left is None:
            print "Angle from point %i with domain to left is None!"%pnt_i
        if angle_right is None:
            print "Angle from point %i with domain to right is None!"%pnt_i
        
        boundary_angle = (angle_right - angle_left) % (2*pi)
        return boundary_angle
        

    def plot_bad_bcs(self):
        bad_bcs = ((self.edges['marker'] == 0) != (self.edges['cells'][:,1] >= 0))
        self.plot(edge_mask = bad_bcs)

    def plot_nodes(self,ids=None):
        if ids is None:
            ids = arange(self.Npoints())
            if self.default_clip is not None:
                c = self.default_clip
                valid = (self.nodes['x'][:,0] > c[0]) & (self.nodes['x'][:,0]<c[1]) & \
                        (self.nodes['x'][:,1] > c[2]) & (self.nodes['x'][:,1]<c[3])
                ids= ids[valid]
            
        [annotate(str(i),self.nodes['x'][i]) for i in ids if self.nodes['stat'][i]==STAT_OK]
        
    def plot_edge_marks(self,edge_mask=None,clip=None):
        """ label edges with c[nc1]-j[j],mark-c[nc2],
        rotated so it reads in the correct orientation for nc1, nc2
        edge_mask should be a boolean array of size Nedges()
        clip can be a list like matplotlib axis() - [xmin,xmax,ymin,ymax]
        """
        if clip is None:
            clip = self.default_clip

        if edge_mask is None and clip:
            ec = self.edge_centers()
            edge_mask = self.edges[:,0] >=0 & ((ec[:,0] >= clip[0]) & (ec[:,0]<=clip[1]) \
                            & (ec[:,1] >= clip[2]) & (ec[:,1]<=clip[3]) )
        else:
            edge_mask = self.edges[:,0] >= 0

        for e in nonzero(edge_mask)[0]:
            delta = self.nodes['x'][ self.edges['nodes'][e,1]] - self.nodes['x'][self.edges['nodes'][e,0]]
            angle = arctan2(delta[1],delta[0])
            annotate("c%d-j%d,%d-c%d"%(self.edges['cells'][e,0],e,self.edges['marker'][e,2],self.edges['cells'][e,1]),
                     ec[e],rotation=angle*180/pi - 90,ha='center',va='center')

    def plot(self,voronoi=False,line_collection_args={},
             all_cells=True,edge_values=None,
             edge_mask=None,vmin=None,vmax=None,ax=None,clip=None):
        """ vmin: if nan, don't set an array at all for the edges
        edge_values: defaults to the edge marker.
        """
        if ax is None:
            ax = gca()

        if self.Ncells() == 0:
            voronoi = False
            
        if voronoi:
            self.vor_plot = ax.plot(self.vcenters()[:,0],self.vcenters()[:,1],".")

        if self.Nedges() == 0:
            return

        if edge_mask is None:
            if not all_cells:
                edge_mask = self.edges['cells'][:,1] < 0
            else:
                edge_mask = (self.edges['stat'] == STAT_OK)

        if sum(edge_mask) == 0:
            return
        
        # g.edges[:,:2] pulls out every edge, and just the endpoint
        #   indices.
        # indexing points by this maps the indices to points
        # which then has the z-values sliced out

        segments = self.nodes['x'][self.edges['nodes'][edge_mask]]

        clip = clip or self.default_clip
        if clip:
            # segments is Nedges * {a,b} * {x,y}
            points_visible = (segments[...,0] >= clip[0]) & (segments[...,0]<=clip[1]) \
                             & (segments[...,1] >= clip[2]) & (segments[...,1]<=clip[3])
            # so now clip is a bool array of length Nedges
            clip = any( points_visible, axis=1)
            segments = segments[clip,...]
        
        line_coll = LineCollection(segments,**line_collection_args)

        if vmin is not None and isnan(vmin):
            print "Skipping the edge array"
        else:
            # allow for coloring the edges
            if edge_values is None:
                edge_values = self.edges['marker'][edge_mask]

            if clip is not None:
                edge_values = edge_values[clip]

            line_coll.set_array(edge_values)

            if vmin:
                line_coll.norm.vmin = vmin
            if vmax:
                line_coll.norm.vmax = vmax
        
        ax.add_collection(line_coll)

        self.edge_collection = line_coll

        axis('equal')
        if not voronoi:
            # the collections themselves do not automatically set the
            # bounds of the axis
            ax.axis(self.bounds())

    def plot_scalar(self,scalar,pdata=None,clip=None,ax=None,norm=None,cmap=None):
        """ Plot the scalar assuming it sits at the center of the
        cells (i.e. use the voronoi centers)
        scalar should be a 1d array, with length the same as the
        number of cells
        """
        if ax is None:
            ax = gca()
        
        if not pdata:
            # create a numpy array for all of the segments:
            # each segment has 4 points so that it closes the triangle
            segments = zeros((self.Ncells(),4,2),float64)
            for i in range(self.Ncells()):
                for j in range(4):
                    segments[i,j,:] = self.nodes['x'][self.cells['nodes'][i,j%3]]

            if clip:
                good_points = (self.nodes['x'][:,0] > clip[0]) & \
                              (self.nodes['x'][:,0] < clip[1]) & \
                              (self.nodes['x'][:,1] > clip[2]) & \
                              (self.nodes['x'][:,1] < clip[3])
                # how to map that onto segments?
                good_verts = good_points[self.cells['nodes']]
                good_cells = good_verts.sum(axis=1) == 3
                segments = segments[good_cells]
                scalar = scalar[good_cells]
                if len(scalar) == 0:
                    return None
                
            mask = isnan(scalar)
            if any(mask):
                segments = segments[~mask]
                scalar = scalar[~mask]
                if len(scalar) == 0:
                    return None
                
            patch_coll = PolyCollection(segments,edgecolors='None',antialiaseds=0,norm=norm,cmap=cmap)
            # is this sufficient for coloring? YES
            patch_coll.set_array(scalar)
            pdata = patch_coll
            ax.add_collection(patch_coll)

            ax.axis('equal')

            ax.axis(self.bounds())
        else:
            pdata.set_array(scalar)
            draw()
        return pdata
    def animate_scalar(self,scalar_frames,post_proc=None):
        clf() # clear figure, to get rid of colorbar, too
        vmin = scalar_frames.min()
        vmax = scalar_frames.max()
        print "Max,min: ",vmax,vmin
        
        pdata = self.plot_scalar(scalar_frames[0])
        title("Step 0")
        pdata.norm.vmin = vmin
        pdata.norm.vmax = vmax
        colorbar(pdata)

        show()

        for i in range(1,scalar_frames.shape[0]):
            title("Step %d"%i)
            self.plot_scalar(scalar_frames[i],pdata)
            if post_proc:
                post_proc()
        
    def bounds(self):
        valid = self.nodes['stat']==STAT_OK
        return (self.nodes['x'][valid,0].min(),self.nodes['x'][valid,0].max(),
                self.nodes['x'][valid,1].min(),self.nodes['x'][valid,1].max() )

    _calc_vcenters = False # Are voronoi centers up to date?
    def vcenters(self):
        if not self._calc_vcenters:
            p1 = self.nodes['x'][self.cells['nodes'][:,0]]
            p2 = self.nodes['x'][self.cells['nodes'][:,1]]
            p3 = self.nodes['x'][self.cells['nodes'][:,2]]

            self.cells['vor_center'] = circumcenter(p1,p2,p3)
            self._calc_vcenters = True

        return self.cells['vor_center']
    
    def faces(self,i):
        # returns an 3 element array giving the edge indices for the
        # cell i
        # the 0th edge goes from the 0th vertex to the 1st.
        f = array([-1,-1,-1])
        for nf in range(3):
            f[nf] = self.nodes_to_edge( (self.cells['nodes'][i,nf],
                                         self.cells['nodes'][i,(nf+1)%3]) )
        return f

    def add_from_shp(self,shp,tolerance=0.0,on_edge=None,on_node=None):
        """ Add nodes and edges from a LineString shapefile
        tolerance: assume nodes which are at least this close to each other in the input
        are actually the same node.

        returns array of new node indices, array of new edge indices

        on_edge: None, or func(self,j,feat), to set additional edge attributes
        on_node: None, or func(self,n,feat), to set additional node attributes
        """
        ods = ogr.Open(shp)
        layer = ods.GetLayer(0)

        new_nodes = []
        new_edges = []
        
        while 1:
            feat = layer.GetNextFeature()
            if feat is None:
                break
            geom = wkb.loads( feat.GetGeometryRef().ExportToWkb() )
            points = array(geom)
            
            Npoints = len(points)
            
            last_node = None
            for i in range(Npoints):
                pnt = points[i]
                n = None

                if self.Npoints()>0:
                    n_best = self.closest_node(pnt)
                    if norm(pnt - self.nodes['x'][n_best]) <= tolerance:
                        n = n_best
                if n is None:
                    n = self.add_node(pnt)
                    if on_node:
                        on_node(self,n,feat)
                    new_nodes.append(n)

                if last_node is not None:
                    if n!=last_node:
                        j = self.add_edge(n,last_node)
                        if on_edge:
                            on_edge(self,j,feat)
                        new_edges.append(j)
                last_node = n
        return array(new_nodes), array(new_edges)

    def write_shp(self,shpname,only_boundaries=1,edge_mask=None,overwrite=0,srs_text='EPSG:26910'):
        """ Write some portion of the grid to a shapefile.
        If only_boundaries is specified, write out only the edges that have non-zero marker

        For starters, this writes every edge as a separate feature, but at some point it
        may make polygons out of the edges.
        """
        if edge_mask is None:
            if only_boundaries:
                edge_mask = (self.edges['stat']==STAT_OK) & (self.edges['marker'] != 0)
            else:
                edge_mask = (self.edges['stat']==STAT_OK)

        input_wkbs = []
        for j in nonzero(edge_mask)[0]:
            e = self.edges['nodes'][j]
            geo = geometry.LineString( [self.nodes['x'][e[0]], self.nodes['x'][e[1]]] )
            input_wkbs.append( geo )

        # Create the shapefile
        if len(input_wkbs) == 0:
            print "NO FEATURES TO WRITE!"
            return
        wkb2shp.wkb2shp(shpname,
                        input_wkbs,
                        fields = self.edges[edge_mask],
                        srs_text = srs_text,
                        overwrite=overwrite)
        

    def write_contours_shp(self,shpname,cell_depths,V,overwrite=False):
        """ like write_shp, but collects edges for each depth in V.
        """
        # because that's how suntans reads depth - no sign
        V = abs(V)
        cell_depths = abs(cell_depths)

        if overwrite and os.path.exists(shpname):
            os.unlink(shpname)

        # Create the shapefile
        drv = ogr.GetDriverByName('ESRI Shapefile')
        ods = drv.CreateDataSource(shpname)
        srs = osr.SpatialReference()
        srs.SetFromUserInput('EPSG:26910')

        olayer = ods.CreateLayer(shpname,
                                 srs=srs,
                                 geom_type=ogr.wkbLineString)

        # create some fields:
        olayer.CreateField(ogr.FieldDefn('depth',ogr.OFTReal))
        olayer.CreateField(ogr.FieldDefn('edge',ogr.OFTInteger))
        
        fdef = olayer.GetLayerDefn()

        internal = (self.edges[:,4] >= 0)

        for v in V:
            print "Finding contour edges for depth=%f"%v

            # These could be tweaked a little bit to get closed polygons
            on_contour = (cell_depths[self.edges[:,3]] <= v ) != (cell_depths[self.edges[:,4]] <= v)
            edge_mask = on_contour & internal

            for j in nonzero(edge_mask)[0]:
                e = self.edges['nodes'][j]
                geo = geometry.LineString( [self.nodes['x'][e[0]], self.nodes['x'][e[1]]] )

                new_feat_geom = ogr.CreateGeometryFromWkb( geo.wkb )

                feat = ogr.Feature(fdef)
                feat.SetGeometryDirectly(new_feat_geom)
                feat.SetField('depth',float(v))
                feat.SetField('edge',int(j))

                olayer.CreateFeature(feat)
            olayer.SyncToDisk()

    def carve_thalweg(self,depths,threshold,start,mode,max_count=None):
        """  Ensures that there is a path of cells from the given start edge
        to deep water with all cells of at least threshold depth.

        start: edge index
        
        depths and threshold should all be as *soundings* - i.e. positive 

        mode is 'cells' - cell-centered depths
           or 'edges' - edge-centered depths

        max_count: max number of cells/edges to deepen along the path (starting
          at start).

        Modifies depths in place.
        """
        c = self.edges['cells'][start,0]

        # approach: breadth-first search for a cell that is deep enough.

        # Track who's been visited -
        # this records the index of the cell from which this cell was visited.
        visitors = -1 * ones(self.Ncells(),int32)

        # Initialize the list of cells to visit
        stack = [c]
        visitors[c] = c # sentinel - visits itself

        gold = None
        
        try:
            while 1:
                new_stack = []

                for c in stack:
                    # find the neighbors of this cell:
                    edges = self.cell_to_edges(c)
                    for e in edges:
                        if mode == 'edges' and depths[e] > threshold:
                            gold = c
                            raise StopIteration

                        # find the neighbor cell
                        if self.edges['cells'][e,0] == c:
                            nc = self.edges['cells'][e,1]
                        else:
                            nc = self.edges['cells'][e,0]
                            
                        # have the neighbor, but should we visit it?
                        if nc < 0 or visitors[nc] >= 0:
                            continue
                        visitors[nc] = c
                        new_stack.append(nc)

                        if mode == 'cells' and depths[nc] > threshold:
                            gold = nc
                            raise StopIteration

                # everyone at this level has been visited and we haven't hit gold.
                # on to the next ring of neighbors:
                stack=new_stack
        except StopIteration:
            pass

        # then trace back and update all the depths that are too small
        c = gold

        along_the_path = []
        while c != visitors[c]:
            if mode == 'edges':
                e = self.cells_to_edge( [c,visitors[c]])
                along_the_path.append(e)
                #if depths[e] < threshold:
                #    depths[e] = threshold
                
            c=visitors[c]
            if mode == 'cells':
                along_the_path.append(c)
                #if depths[c] < threshold:
                #    depths[c] = threshold
        if max_count is None or max_count > len(along_the_path):
            max_count = len(along_the_path)
        for item in along_the_path[-max_count:]:
            if depths[item] < threshold:
                depths[item] = threshold
                    
        # Take care of starting edge
        if mode == 'edges' and depths[start] < threshold:
            depths[start] = threshold

    def write_suntans(self,pathname):
        """ create cells.dat, edges.dat and points.dat
        from the TriGrid instance, all in the directory
        specified by pathname
        """
        if not os.path.exists(pathname):
            print "Creating folder ",pathname
            os.makedirs(pathname)
        
        # check for missing BCs
        missing_bcs = (self.edges['marker']==0) & (self.edges['cells'][:,1]<0)
        n_missing = missing_bcs.sum()
        if n_missing > 0:
            print "WARNING: %d edges are on the boundary but have marker==0"%n_missing
            print "Assuming they are closed boundaries!"
            # make a copy so that somebody can plot the bad cells afterwards
            # with plot_missing_bcs()
            my_edges = self.edges.copy()
            my_edges['marker'][missing_bcs] = LAND_EDGE
        else:
            my_edges = self.edges
        
        cells_fp = open(pathname + "/cells.dat","w")
        edges_fp = open(pathname + "/edges.dat","w")
        points_fp= open(pathname + "/points.dat","w")

        for i in range(self.Npoints()):
            points_fp.write("%.5f %.5f 0\n"%(self.nodes['x'][i,0],self.nodes['x'][i,1]))
        points_fp.close()
        

        # probably this can be done via the edges array
        vc = self.vcenters()
        
        for i in range(self.Ncells()):
            # each line in the cell output is
            # x, y of voronoi center (I think)
            # zero-based point-indices x 3
            # zero-based ?cell? indices x 3, for neighbors?

            # find the neighbors:
            # the first neighbor: need another cell that has
            # both self.cells[i,0] and self.cells[i,1] in its
            # list.
            my_set = set([i])
            n = [-1,-1,-1]
            for j in 0,1,2:
                adj1 = self.node_to_cells(self.cells['nodes'][i,j])
                adj2 = self.node_to_cells(self.cells['nodes'][i,(j+1)%3])
                neighbor = adj1.intersection(adj2).difference(my_set)
                if len(neighbor) == 1:
                    n[j] = neighbor.pop()
            
            cells_fp.write("%.5f %.5f %i %i %i %i %i %i\n"%(
                    vc[i,0],vc[i,1],
                    self.cells['nodes'][i,0],self.cells['nodes'][i,1],self.cells['nodes'][i,2],
                    n[0],n[1],n[2]))
        cells_fp.close()

        for edge in my_edges:
            # point_id, point_id, edge_type, cell, cell
            edges_fp.write("%i %i %i %i %i\n"%(
                edge['nodes'][0],edge['nodes'][1],
                edge['marker'],
                edge['cells'][0],edge['cells'][1]))
            
        edges_fp.close()

    def nodes_to_edge(self,nodes):
        # this way is slow - most of the time in the array ops
        el0 = self.node_to_edges(nodes[0])
        el1 = self.node_to_edges(nodes[1])
        for e in el0:
            if e in el1:
                return e
        raise NoSuchEdgeError,str(nodes)

    def nodes_to_cell(self,nodes):
        """ return the cell (if any) that is made up of the given nodes
        depends on pnt2cells
        """
        try:
            cells_a = self.node_to_cells(nodes[0])
            cells_b = self.node_to_cells(nodes[1])
            cells_c = self.node_to_cells(nodes[2])

            c = cells_a.intersect1d(cells_b).intersect1d(cells_c)

            if len(c) == 0:
                raise NoSuchCellError()
            elif len(c) > 1:
                raise Exception,"Nodes %s mapped to cells %s"%(nodes,c)
            else:
                return c[0]
        except KeyError:
            raise NoSuchCellError()

    def node_neighbors(self,node_id):
        """ Return a list of nodes which share an edge with the given node
        """
        nbrs = []
        for j in self.node_to_edges(node_id):
            if self.edges['nodes'][j,0] == node_id:
                nbrs.append( self.edges['nodes'][j,1])
            else:
                nbrs.append( self.edges['nodes'][j,0])
        return nbrs
    
    def cell_neighbors(self,cell_id,adjacent_only=0):
        """ return array of cell_ids for neighbors of this
        cell.  here neighbors are defined by sharing a vertex,
        not just sharing an edge, unless adjacent_only is specified.
        (in which case it only returns cells sharing an edge)
        """
        if not adjacent_only:
            neighbors = [self.node_to_cells(p) for p in self.cells['nodes'][cell_id]]
            return unique(reduce(lambda x,y: x+y,neighbors))
        else:
            nbrs = []
            for nc1,nc2 in self.edges['cells'][self.cell_to_edges(cell_id)]:
                if nc1 != cell_id and nc1 >= 0:
                    nbrs.append(nc1)
                if nc2 != cell_id and nc2 >= 0:
                    nbrs.append(nc2)
            return array(nbrs)

    _components = None
    def calculate_connected_components(self):
        """ Returns an array mapping each node to an id of its connected component.
        Also saves it to self._components - note that this is not guaranteed to be
        fresh...
        """
        self._components = components = -1 * ones( self.Npoints(), int32)

        ncomponents = 0
        for i in self.valid_node_iter():
            # scan ahead to find the next unlabeled point
            if components[i] >= 0:
                continue

            # so i is the starting point for the connected components search
            stack = [i]
            while len(stack):
                visit = stack.pop()
                components[visit]=ncomponents

                for nbr in self.node_neighbors(visit):
                    if components[nbr] < 0:
                        stack.append(nbr)
            ncomponents+=1

        return components

    def subgrids_by_connected_components(self):
        """ returns an iterator over subgrids defined by the connected components of this
        grid.  Components with only a node, and no edges, will be ignored.
        also only handles nodes and edges - no attempt is made to transfer cells.
        """
        cc=0
        # track which nodes in the global array have been processed
        # 0: unvisited
        # 1: visited
        nodes_processed = zeros( self.Nnodes(), int8 )
        extra_node_fields = self.extra_node_fields()
        extra_edge_fields = self.extra_edge_fields()

        for n in self.valid_node_iter():
            if nodes_processed[n]:
                continue

            # we've found an unprocessed node - 
            subg = self.copy(types_only=True)
            g2l_nodes = {}

            stack = [n]
            # global edge indexes to add later:
            edges_to_add=set()
            while len(stack):
                v = stack.pop(-1) # visited node
                nodes_processed[v]=1

                ln = g2l_nodes[v] = subg.add_node( self.nodes['x'][v] )
                for enf in extra_node_fields:
                    subg.nodes[ln][enf] = self.nodes[v][enf]

                for j in self.node_to_edges(v):
                    edges_to_add.add(j)
                    for nbr in self.edges[j]['nodes']:
                        if nbr != v and not g2l_nodes.has_key(nbr):
                            stack.append(nbr)

            if subg.Nnodes() < 2:
                continue

            print "processing non-unit connected component %d"%cc

            for j in edges_to_add:
                gnodes = self.edges[j]['nodes']
                ln0 = g2l_nodes[gnodes[0]]
                ln1 = g2l_nodes[gnodes[1]]
                lj = subg.add_edge( ln0,ln1,
                                    cleft=-1,cright=-1 )
                for eef in extra_edge_fields:
                    subg.edges[lj][eef] = self.edges[j][eef]

            yield subg
            cc+=1

        
    def make_edges_from_cells(self):
        # iterate over cells, and for each cell, if it's index
        # is smaller than a neighbor or if no neighbor exists,
        # write an edge record
        edges = []
        default_marker = 0

        # this will get built on demand later.
        self._node_to_edges = None
        
        for i in range(self.Ncells()):
            # find the neighbors:
            # the first neighbor: need another cell that has
            # both self.cells[i,0] and self.cells[i,1] in its
            # list.
            my_set = [i]
            n = [-1,-1,-1]
            for j in 0,1,2:
                pnt_a = self.cells['nodes'][i,j]
                pnt_b = self.cells['nodes'][i,(j+1)%3]
                    
                adj1 = self.node_to_cells(pnt_a) # cells that use pnt_a
                adj2 = self.node_to_cells(pnt_b) # cells that use pnt_b

                # the intersection is us and our neighbor
                #  so difference out ourselves...
                neighbor = setdiff1d( intersect2d(adj1,adj2),my_set)
                # and maybe we ge a neighbor, maybe not (we're a boundary)
                if len(neighbor) == 1:
                    n = neighbor.pop()
                else:
                    n = -1
                    
                if n==-1 or i<n:
                    # we get to add the edge:
                    edges.append((pnt_a,
                                  pnt_b,
                                  default_marker,
                                  i,n))

        # self.edges = array(edges,int32)
        self.edges = zeros(len(edges), self.edge_dtype)
        self.edges['nodes'] = edges[:,0:2]
        self.edges['marker'] =  edges[:,2]
        self.edges['cells'] = edges[:,3:5]

    def verify_bc(self,do_plot=True):
        """ check to make sure that all grid boundaries have a BC set 
        """
        #  point_i, point_i, marker, cell_i, cell_i

        # marker: 0=> internal,1=> closed, 3=> open
        #  make sure that any internal boundary has a second cell index
        #  assumes that all edges have the first cell index != -1
        bad_edges = find( (self.edges['marker'][:,2]==0) & (self.edges['cells'][:,1]<0 ) )

        if do_plot:
            for e in bad_edges:
                bad_points = self.edges['nodes'][e]
            
                plot(self.nodes['x'][bad_points,0],
                     self.nodes['x'][bad_points,1],'r-o')
        
        if len(bad_edges) > 0:
            print "BAD: there are %d edges without BC that have only 1 cell"%len(bad_edges)
            return 0
        else:
            return 1

    def cell_to_edges(self,cell_i):
        if self.cells['stat'][cell_i] != STAT_OK:
            raise "cell %i has been deleted"%cell_i
        
        # return indices to the three edges for this cell:
        nodes = self.cells['nodes'][cell_i] # the three vertices

        edges = [ self.nodes_to_edge( (nodes[i], nodes[(i+1)%3]) ) for i in range(3) ]
        return edges

    def cells_to_edge(self,nc):
        e1 = self.cell_to_edges(nc[0])
        e2 = self.cell_to_edges(nc[1])
        for e in e1:
            if e in e2:
                return e
        raise Exception,"Cells %d and %d don't share an edge"%(nc[0],nc[1])

    def build_index(self):
        if self.index is None:
            # assemble points into list of (id, [x x y y], None)
            if self.verbose > 1:
                print "building point index"
            # old rtree required that stream inputs have non-interleaved coordinates,
            # but new rtree allows for interleaved coordinates all the time.
            # best solution probably to specify interleaved=False

            
            # tuples = [(i,self.nodes['x'][i,xxyy],None)
            #           for i in range(self.Npoints())
            #           if self.nodes['stat'][i]==STAT_OK ]

            # Can we make that a generator?
            def gen_tuples():
                for i in range(self.Nnodes()):
                    if self.verbose and i % 10000 == 0:
                        print "build_index: %d/%d"%(i,self.Nnodes())
                    if self.nodes['stat'][i] == STAT_OK:
                        yield i,self.nodes['x'][i,xxyy],None
                
            self.index = Rtree(gen_tuples(),interleaved=False)
            if self.verbose > 1:
                print "done"
            
    def build_edge_index(self):
        if self.edge_index is None:
            print "building edge index"
            ec = self.edge_centers()
            tuples = [(i,ec[i,xxyy],None) for i in range(self.Nedges())]
            self.edge_index = Rtree(tuples,interleaved=False)
            print "done"

    def closest_node(self,p,count=1,boundary=0):
        """ Returns the count closest nodes to p
        boundary=1: only choose nodes on the boundary.
        """
        if boundary:
            # print "Searching for nearby boundary point"
            # this is slow, but I'm too lazy to add any sort of index specific to
            # boundary nodes.  Note that this will include interprocessor boundary
            # nodes, too.
            boundary_nodes = unique( self.edges['nodes'][self.edges['marker']>0] )
            dists = sum( (p - self.nodes['x'][boundary_nodes])**2, axis=1)
            order = argsort(dists)
            closest = boundary_nodes[ order[:count] ]
            
            if count == 1:
                return closest[0]
            else:
                return closest
        else:
            if self.index is None:
                self.build_index()

            p = array(p)

            # returns the index of the grid point closest to the given point:
            hits = self.index.nearest( p[xxyy], count)

            # newer versions of rtree return a generator:
            if isinstance( hits, types.GeneratorType):
                # so translate that into a list like we used to get.
                h = []
                for i in range(count):
                    try:
                        h.append(hits.next())
                    except StopIteration:
                        break
                hits = h
                        
            if count > 1:
                return hits
            elif len(hits)>0:
                return hits[0]
            else:
                return None

    def closest_edge(self,p):
        """ Return the index of the edge whose center is nearest the point p.
        """
        if self.edge_index is None:
            self.build_edge_index()

        p = array(p)
        hits = self.edge_index.nearest( p[xxyy], 1)
        
        # newer versions of rtree return a generator:
        if isinstance( hits, types.GeneratorType):
            # so translate that into a list like we used to get.
            return hits.next()
        else:
            return hits[0]
    
    def closest_cell(self,p,full=0,inside=False):
        """
        full=0: return None if the closest *point* is not in a cell on this subdomain
        full=1: exhaustively search all cells, even if the nearest point is not on this subdomain

        inside: require that the returned cell contains p, otherwise return None
        """
        # rather than carry around another index, reuse the point index
        p = array(p)
        
        i = self.closest_node(p)
        try:
            cells = self.node_to_cells(i)
        except KeyError:
            if not full:
                return None
            else:
                print "This must be on a subdomain.  The best point wasn't in one of our cells"
                cells = range(self.Ncells())

        if inside:
            pnt = geometry.Point(p[0],p[1])
            for c in cells:
                tri = geometry.Polygon(self.nodes['x'][self.cells['nodes'][c]])
                if tri.contains(pnt):
                    return c
            return None
        else:
            vc = self.vcenters()
            dists = ((p-vc[cells])**2).sum(axis=1)
            chosen = cells[argmin(dists)]

            # dist = sqrt( ((p-vc[chosen])**2).sum() )
            # print "Closest cell was %f [m] away"%dist
        return chosen
        
    def set_edge_markers(self,pnt1,pnt2,marker):
        """ Find the nodes closest to each of the two points,
        Search for the shortest path between them on the boundary.
        Set all of those edges' markers to marker
        """
        n1 = self.closest_node(pnt1)
        n2 = self.closest_node(pnt2)

        path,length = self.shortest_path(n1,n2,boundary_only=1)

        for i in range(len(path)-1):
            e = self.nodes_to_edge( path[i:i+2] )
            self.edges['marker'][e] = marker

    def shortest_path(self,n1,n2,boundary_only=0,max_cost = inf,valid_edges = None):
        """ dijkstra on the edge graph from n1 to n2
        boundary_only: limit search to edges on the boundary (have
        a -1 for cell2)
        returns the list of node indices and the total length of the path, or None,None
        when max_cost was encountered or no path exists.

        if valid_edges is set, it is a boolean array with True for edges which are
          allowed for finding the shortest path.
        """
        if self._components is not None:
            if self._components[n1] != self._components[n2]:
                return None,None
        
        queue = pq.priorityDictionary()
        queue[n1] = 0

        done = {}

        if boundary_only:
            if valid_edges:
                raise Exception,"No support for valid_edges and boundary_only at the same time"
            valid_edges = (self.edges[:]['cells'][1] == BOUNDARY)

        while 1:
            # find the queue-member with the lowest cost:
            if len(queue)==0:
                return None,None # no way to get there from here.
            best = queue.smallest()
            best_cost = queue[best]
            if best_cost > max_cost:
                # print "Too far"
                return None,None
            
            del queue[best]

            done[best] = best_cost

            if best == n2:
                # print "Found the ending point"
                break

            # figure out its neighbors
            # Why did I do it this way??
            #cells = self.node_to_cells(best)
            #all_points = unique( self.cells['nodes'][cells] )

            # Better to use edges, since we may have edges, but no cells
            edges = self.node_to_edges(best)
            all_nodes = unique( self.edges['nodes'][edges] )
            
            for p in all_nodes:
                if done.has_key(p):
                    # both for best and for points that we've already done
                    continue
                
                if valid_edges is not None:
                    e = self.nodes_to_edge( (best,p) )
                    if not valid_edges[e]:
                        continue
                
                dist = sqrt( ((self.nodes['x'][p] - self.nodes['x'][best])**2).sum() )
                new_cost = best_cost + dist

                if not queue.has_key(p):
                    queue[p] = inf

                if queue[p] > new_cost:
                    queue[p] = new_cost

        # reconstruct the path:
        path = [n2]

        while 1:
            p = path[-1]
            if p == n1:
                break

            # figure out its neighbors
            # cells = self.node_to_cells(p)
            # all_points = unique( self.cells['nodes'][cells] )
            edges = self.node_to_edges(p)
            all_nodes = unique(self.edges['nodes'][edges])

            found_prev = 0
            for nbr in all_nodes:
                if nbr == p or not done.has_key(nbr):
                    continue

                dist = sqrt( ((self.nodes['x'][p] - self.nodes['x'][nbr])**2).sum() )

                if done[p] == done[nbr] + dist:
                    path.append(nbr)
                    found_prev = 1
                    break
            if not found_prev:
                return None,None

        return array( path[::-1] ), done[n2]

    def spanning_tree(self,nodes):
        """ Return an array of edges which link all of the given nodes
        no real attempt is made to make this optimal - in practice it really
        shouldn't matter since the nodes will always be limited in number and
        very close to each other.
        it does check for the simplest case of two nodes either joined by an
        existing edge
        """
        nodes = list(nodes)
        if len(nodes) < 2:
            return [] # nothing to span
        a = nodes[0]
        spanning_edges = set()
        for i in range(1,len(nodes)):
            b = nodes[i]
            # the most likely case:
            try:
                j = self.nodes_to_edge( [a,b] )
                spanning_edges.add(j)
            except NoSuchEdgeError:
                nlist,path_length = self.shortest_path(a,b)
                for c_d in zip( nlist[:-1],
                                nlist[1:] ):
                    spanning_edges.add( self.nodes_to_edge( c_d ) )
        return spanning_edges

    ### graph modification api calls

    def delete_node_and_merge(self,n):
        """ For a degree 2 node, remove it and make one edge out its two edges.
        this used to be in paver, but I don't think there is any reason it can't
        be here in trigrid.
        """
        edges = self.node_to_edges(n)
        
        if self.verbose > 1:
            print "Deleting node %d, with edges %s"%(n,edges)

        if len(edges) == 2:
            if self.verbose > 1:
                print "  deleting node %d, will merge edges %d and %d"%(n,edges[0],edges[1])
            e = self.merge_edges( edges[0], edges[1] )
        elif len(edges) != 0:
            print "Trying to delete node",n
            annotate("del",self.nodes['x'][n])
            print "Edges are:",self.edges[edges]
            
            raise Exception,"Can't delete node with %d edges"%len(edges)

        edges = self.node_to_edges(n)

        if len(edges) != 0:
            print "Should have removed all edges to node %d, but there are still some"%n
            
        self.delete_node(n)
        return e
    
    def unmerge_edges(self,e1,e2,e1data,e2data):
        self.edges[e1] = e1data
        self.edges[e2] = e2data

        # too lazy to do this right now, so to be safe just kill it
        self._node_to_edges = None
        
    def merge_edges(self,e1,e2):
        """ returns the id of the new edge, which for now will always be one of e1 and e2
        (and the other will have been deleted
        """
        if self.verbose > 1:
            print "Merging edges %d %d"%(e1,e2)
            print " edge %d: nodes %d %d"%(e1,self.edges['nodes'][e1,0],self.edges['nodes'][e1,1])
            print " edge %d: nodes %d %d"%(e2,self.edges['nodes'][e2,0],self.edges['nodes'][e2,1])

        
        B = intersect1d( self.edges['nodes'][e1], self.edges['nodes'][e2] )[0]

        # try to keep ordering the same (not sure if this is necessary)
        if self.edges['nodes'][e1,0] == B:
            e1,e2 = e2,e1

        # push the operation with the re-ordered edge nodes, so that we know (i.e.
        # live_dt knows) which of the edges is current, and which is being undeleted.
        self.push_op(self.unmerge_edges, e1, e2, self.edges[e1].copy(), self.edges[e2].copy() )
            
        # pick C from e2
        if self.edges['nodes'][e2,0] == B:
            C = self.edges['nodes'][e2,1]
        else:
            C = self.edges['nodes'][e2,0]

        if self.edges['nodes'][e1,0] == B:
            self.edges['nodes'][e1,0] = C
            A = self.edges['nodes'][e1,1]
        else:
            self.edges['nodes'][e1,1] = C
            A = self.edges['nodes'][e1,0]

        # this removes e2 from _pnt2edges for B & C
        # because of mucking with the edge data, better to handle the
        # entire rollback in merge_edges
        self.delete_edge(e2,rollback=0)
        
        # fix up edge lookup tables:
        if self._node_to_edges is not None:
            self._node_to_edges[C].append(e1)

            # B is still listed for e1
            b_edges = self._node_to_edges[B]
            if b_edges != [e1]:
                print "Merging edges.  Remaining node_to_edges[B=%d] = "%B,b_edges
                print "is not equal to e1 = ",[e1]
            self._node_to_edges[B] = []

        # and callbacks:
        self.updated_edge(e1)
        return e1
        
    def undelete_node(self,i,p):
        self.nodes[i] = p

        if self.index is not None:
            self.index.insert(i, self.nodes['x'][i,xxyy] )

    def delete_node(self,i,remove_edges=1):
        if self.verbose > 1:
            print "delete_node: %d, remove_edges=%s"%(i,remove_edges)
            
        if remove_edges:
            # make a copy so that as delete_edge modifies 
            # _node_to_edges we still have the original list
            nbr_edges = list(self.node_to_edges(i))

            for e in nbr_edges:
                self.delete_edge(e)

        self.push_op(self.undelete_node,i,self.nodes[i].copy())
            
        # nodes are marked as deleted by setting the x coordinate
        # to NaN, and remove from index
        if self.index is not None:
            coords = self.nodes['x'][i,xxyy]
            self.index.delete(i, coords )
            
        self.nodes['stat'][i] = STAT_DEL
        self.deleted_node(i)

    def undelete_cell(self,c,cdata,edge_updates):
        self.cells[c] = cdata
        self._calc_vcenters = False # Lazy

        for e,vals in edge_updates:
            self.edges[e] = vals
            
        if self._node_to_cells is not None: 
            for i in nodes:
                if not self._node_to_cells.has_key(i):
                    self._node_to_cells[i] = set()
                self._node_to_cells[i].append(c)
                
    def delete_cell(self,c,replace_with=-2,rollback=1):
        """ replace_with: the value to set on edges that used to reference
             this cell. -2 => leave an internal hole
                        -1 => create an 'island'
        """
        nA,nB,nC = self.cells['nodes'][c]
        ab = self.nodes_to_edge([nA,nB])
        bc = self.nodes_to_edge([nB,nC])
        ca = self.nodes_to_edge([nC,nA])

        edge_updates = [ [ab,self.edges[ab].copy()],
                         [bc,self.edges[bc].copy()],
                         [ca,self.edges[ca].copy()] ]
                     
        self.push_op(self.undelete_cell,c,self.cells[c].copy(),edge_updates)

        for e in [ab,bc,ca]:
            if self.edges['cells'][e,0] == c:
                check = 0
            elif self.edges['cells'][e,1] == c:
                check = 1
            else:
                print "Cell: %d  check on edge %d  with nbrs: %d %d"%(
                    c,e,self.edges[e]['cells'][0],self.edges[e]['cells'][1])
                
                raise Exception,"Deleting cell, but edge has no reference to it"
            self.edges['cells'][e,check] = replace_with
            
            # optional - update edge marker, and for now just assume it will
            # be a land edge (other BC types are generally handled later anyway)
            if replace_with == -1:
                # print "Deleting cell and replace_with is",replace_with
                if self.edges['marker'][e] == 0:
                    # print "Internal edge becoming a land edge"
                    self.edges['marker'][e] = LAND_EDGE
                    self.updated_edge(e)
        
        self.cells['stat'][c] = STAT_DEL
        self.cells['vor_center'][c] = nan
            
        if self._node_to_cells is not None:
            for n in [nA,nB,nC]:
                self._node_to_cells[n].remove(c)
            
        self.deleted_cell(c)

    def undelete_edge(self,e,e_data):
        self.edges[e] = e_data

        # fix up indexes:
        if self._node_to_edges is not None:
            for n in self.edges['nodes'][e]:
                if not self._node_to_edges.has_key(n):
                    self._node_to_edges[n] = []
                self._node_to_edges[n].append(e)

        if self.edge_index is not None:
            coords = self.edge_centers()[e][xxyy]
            self.edge_index.insert(e,coords)
        
    def delete_edge(self,e,rollback=1):
        """ for now, just make it into a degenerate edge
        specify rollback=0 to skip recording the undo information
        """
        if self.verbose > 1:
            print "Deleting edge %d:"%e

        # remove any neighboring cells first
        cell_nbrs = self.edges['cells'][e]

        if any(cell_nbrs == -1):
            replace_with = -1
        else:
            replace_with = -2
        
        for c in cell_nbrs:
            if c >= 0:
                self.delete_cell(c,replace_with=replace_with,rollback=rollback)

        # clear out indexes
        if self._node_to_edges is not None:
            self._node_to_edges[self.edges['nodes'][e,0]].remove(e)
            self._node_to_edges[self.edges['nodes'][e,1]].remove(e)

        if self.edge_index is not None:
            coords = self.edge_centers()[e][xxyy]
            self.edge_index.delete(e,coords)

        if rollback:
            self.push_op(self.undelete_edge,e,self.edges[e].copy())
        
        # mark edge deleted
        self.edges['stat'][e] = STAT_DEL
        # These aren't strictly needed, but can help in debugging sometimes
        self.edges['marker'][e] = DELETED_EDGE # DELETED
        self.edges['cells'][e] = -37

        # signal to anyone who cares
        self.deleted_edge(e)

    def valid_edges(self):
        """ returns an array of indices for valid edges - i.e. not deleted"""
        return nonzero(self.edges['stat']==STAT_OK)[0]

    # Iterators for valid elements:
    def valid_node_iter(self):
        for i in range(self.Npoints()):
            if self.nodes['stat'][i] == STAT_OK:
                yield i
                
    def valid_edge_iter(self):
        for j in range(self.Nedges()):
            if self.edges['stat'][j] == STAT_OK:
                yield j
    def valid_cell_iter(self):
        for n in range(self.Ncells()):
            if self.cells['stat'][n] == STAT_OK:
                yield n
                
    
    def split_edge(self,nodeA,nodeB,nodeC):
        """ take the existing edge AC and insert node B in the middle of it
        """
        e1 = self.nodes_to_edge([nodeA,nodeC])
                
        if any( self.edges['cells'][e1] >= 0 ):
            print "While trying to split the edge %d (%d-%d) with node %d"%(e1,nodeA,nodeC,nodeB)
            [annotate(str(i),self.nodes['x'][i]) for i in [nodeA,nodeC,nodeB]]
            print "The cell neighbors of the edge are:",self.edges['cells'][e1]
            raise Exception,"You can't split an edge that already has cells"

        # Push the specific operations instead of an aggregate split_edge op
        # self.push_op(self.unsplit_edge,e1,len(self.edges))
        # 2011-01-29: this used to be in the opp. order - but that implies
        #   an invalid state
        self.push_op(self.unmodify_edge,e1,self.edges[e1].copy())
        self.push_op(self.unadd_edge,self.Nedges())
        
        self.edges = array_append( self.edges, self.edges[e1] )
        e2 = self.Nedges() - 1
        
        # first make the old edge from AC to AB
        if self.edges['nodes'][e1,0] == nodeC:
            self.edges['nodes'][e1,0] = nodeB
            self.edges['nodes'][e2,1] = nodeB
        else:
            self.edges['nodes'][e1,1] = nodeB
            self.edges['nodes'][e2,0] = nodeB

        ## Update center
        self.update_edge_center(e1)
        self.update_edge_center(e2)
        

        # handle updates to indices
        #   update node_to_edges
        if self._node_to_edges is not None:
            # nodeA is fine.
            # nodeB has to get both edges:
            self._node_to_edges[nodeB] = [e1,e2]
            # nodeC 
            i = self._node_to_edges[nodeC].index(e1)
            self._node_to_edges[nodeC][i] = e2

        self.updated_edge(e1)
        self.created_edge(e2)

        return e2

    def unadd_edge(self,old_length):
        new_e = old_length

        if self._node_to_edges is not None:
            for n in self.edges['nodes'][new_e]:
                self._node_to_edges[n].remove(new_e)
        
        self.edges = self.edges[:old_length]

    def unmodify_edge(self, e, old_data):
        # print "unmodifying edge %d reverting to %s"%(e,old_data)
        if self._node_to_edges is not None:
            a,b = self.edges['nodes'][e]
            self._node_to_edges[a].remove(e)
            self._node_to_edges[b].remove(e)

            a,b = old_data['nodes']
            self._node_to_edges[a].append(e)
            self._node_to_edges[b].append(e)
            
        self.edges[e] = old_data

    def add_edge(self,nodeA,nodeB,marker=0,cleft=-2,cright=-2,coerce_boundary=None):
        """ returns the number of the edge
        for cells that are marked -2, this will check to see if a new cell can
        be made on that side with other unmeshed edges
        """
        try:
            e = self.nodes_to_edge([nodeA,nodeB])
            raise Exception,"edge betwen %d and %d already exists"%(nodeA,nodeB)
        except NoSuchEdgeError:
            pass

        # dynamic resizing for edges:
        self.push_op(self.unadd_edge,len(self.edges))
        new_edge = zeros( (),dtype=self.edge_dtype )
        new_edge['nodes'] = [nodeA,nodeB]
        new_edge['marker'] = marker
        new_edge['cells'] = [cleft,cright]
        self.edges = array_append( self.edges, new_edge)

        this_edge = self.Nedges()-1
        edge_ab = this_edge # for consistency in the mess of code below

        self.cells_from_last_new_edge = []
        
        if cleft == -2 or cright == -2:
            # First get any candidates, based just on connectivity
            edges_from_a = self.node_to_edges(nodeA) 
            edges_from_b = self.node_to_edges(nodeB) 

            neighbors_from_a = setdiff1d( self.edges['nodes'][edges_from_a].ravel(), [nodeA,nodeB] )
            neighbors_from_b = setdiff1d( self.edges['nodes'][edges_from_b].ravel(), [nodeA,nodeB] )

            # nodes that are connected to both a and b
            candidates = intersect1d( neighbors_from_a, neighbors_from_b )

            if len(candidates) > 0:
                # is there a candidate on our right?
                ab = self.nodes['x'][nodeB] - self.nodes['x'][nodeA]
                ab_left = rot(pi/2,ab)

                new_cells = []

                for c in candidates:
                    ac = self.nodes['x'][c] - self.nodes['x'][nodeA]
                     
                    if dot(ac,ab_left) < 0: # this one is on the right of AB
                        # make a stand-in A & B that are in CCW order in this cell
                        ccwA,ccwB = nodeB,nodeA
                        check_cell_ab = 1 # the relevant cell for the new edge
                    else:
                        ccwA,ccwB = nodeA,nodeB
                        check_cell_ab = 0
                        
                        
                    edge_ac = self.nodes_to_edge((ccwA,c))
                    if self.edges['nodes'][edge_ac,0] == ccwA:
                        # then the edge really is stored ac
                        check_cell_ac = 1
                    else:
                        check_cell_ac = 0

                    edge_bc = self.nodes_to_edge((ccwB,c))
                    if self.edges['nodes'][edge_bc,0] == ccwB:
                        check_cell_bc = 0
                    else:
                        check_cell_bc = 1

                    # so now we have edge_ab, edge_ac, edge_bc as edge ids for the
                    # edges that make up a new cell, and corresponding check_cell_ab
                    # check_cell_ac and check_cell_bc that index the adj. cell that is
                    # facing into this new cell.

                    ccw_edges = [edge_ab,edge_bc,edge_ac]
                    check_cells = [check_cell_ab, check_cell_bc, check_cell_ac]
                    
                    adj_ids = [ self.edges['cells'][e,check]
                                for e,check in zip(ccw_edges,check_cells) ]
                    adj_ids = array( adj_ids )

                    if any(adj_ids >= 0) and any( adj_ids != adj_ids[0]):
                        # bad.  one edge thinks there is already a cell here, but
                        # the others doesn't agree.
                        print "During call to add_edge(nodeA=%d nodeB=%d marker=%d cleft=%d cright=%d coerce=%s"%(nodeA,
                                                                                                                  nodeB,
                                                                                                                  marker,
                                                                                                                  cleft,
                                                                                                                  cright,
                                                                                                                  coerce_boundary)
                        raise Exception,"cell neighbor values for new cell using point %d are inconsistent: %s"%(c,adj_ids)
                    elif all(adj_ids == -1):
                        # leave them be, no new cell, all 3 edges are external
                        pass
                    elif coerce_boundary == -1:
                        # no new cell - everybody gets -1
                        self.edges['cells'][edge_ab,check_cell_ab] = -1
                        self.edges['cells'][edge_ac,check_cell_ac] = -1
                        self.edges['cells'][edge_bc,check_cell_bc] = -1
                    elif all( adj_ids == -2 ) or coerce_boundary == -2:
                        # make new cell, everybody gets that cell id
                        # Create the cell and get it's id:
                        new_cells.append( self.add_cell([ccwA,ccwB,c]) )

                        # update everybody's cell markers:
                        self.push_op(self.unmodify_edge, edge_ac, self.edges[edge_ac].copy() )
                        self.push_op(self.unmodify_edge, edge_bc, self.edges[edge_bc].copy() )
                        self.edges['cells'][edge_ac,check_cell_ac] = new_cells[-1]
                        self.edges['cells'][edge_bc,check_cell_bc] = new_cells[-1]
                        self.edges['cells'][edge_ab,check_cell_ab] = new_cells[-1]

                        # extend boundary - the fun one
                        # only when there was an external edge that now falls inside
                        # the new cell => mark the *other* side of the other edges to
                        # -1
                        if any(adj_ids==-1):
                            for i in range(3):
                                if adj_ids[i] == -2:
                                    # make its outside cell a -1
                                    # the 1-check gives us the outside cell nbr
                                    self.edges['cells'][ccw_edges[i],1-check_cells[i]] = -1
                                    # go ahead and set a closed edge, too
                                    self.edges['marker'][ccw_edges[i]] = LAND_EDGE
                                else:
                                    # as long as this edge wasn't originally -1,-1
                                    # (which ought to be illegal), it's safe to say
                                    # that it is now internal
                                    self.edges['marker'][ccw_edges[i]] = 0 # internal
                                # either way, let people know that markers have changed,
                                # but wait until later to signal on the new edge since
                                # it is not in the indices yet
                                if ccw_edges[i] != edge_ab:
                                    self.updated_edge(ccw_edges[i])

                self.cells_from_last_new_edge = new_cells

        # update node_to_edges
        if self._node_to_edges is not None:
            for n in [nodeA,nodeB]:
                if not self._node_to_edges.has_key(n):
                    self._node_to_edges[n] = []
                self._node_to_edges[n].append(this_edge)

        self.created_edge(this_edge)
                
        return this_edge

    def unadd_node(self,old_length):
        if self.index is not None:
            curr_len = len(self.nodes)
            for i in range(old_length,curr_len):
                coords = self.nodes['x'][i,xxyy]
                self.index.delete(i, coords )
        
        self.nodes = self.nodes[:old_length]
        
    def add_node(self,P):
        self.push_op(self.unadd_node,len(self.nodes))
        new_node = zeros( (), self.node_dtype)
        new_node['x'][:2] = P
        self.nodes = array_append( self.nodes, new_node )

        new_i = self.Npoints() - 1

        if self.index is not None:
            self.index.insert(new_i, self.nodes['x'][new_i,xxyy] )

        self.created_node(new_i)
        
        return new_i

    def unadd_cell(self,old_length):
        # remove entries from _node_to_cells
        #  the cell that was added is at the end:
        if self._node_to_cells is not None:
            new_c = old_length
            for n in self.cells['nodes'][new_c]:
                self._node_to_cells[n].remove(new_c)
        
        self.cells = self.cells[:old_length]
                                
    def add_cell(self,c):
        c = array(c)
        self.push_op(self.unadd_cell,len(self.cells))
        
        i = array([0,1,2])
        ip = array([1,2,0])
        xi =  self.nodes['x'][c[i],0]
        yi =  self.nodes['x'][c[i],1]
        xip = self.nodes['x'][c[ip],0]
        yip = self.nodes['x'][c[ip],1]
        
        A = 0.5 * (xi*yip-xip*yi).sum()
        if A < 0:
            print "WARNING: attempt to add CW cell.  Reversing"
            c = c[::-1]
            
        new_cell = zeros( (), self.cell_dtype)
        new_cell['nodes'] = c
        
        self.cells = array_append( self.cells, new_cell )
        self._calc_vcenters = False
        
        this_cell = self.Ncells() - 1

        if self._node_to_cells is not None: # could be smarter and actually update.
            for i in c:
                if not self._node_to_cells.has_key(i):
                    self._node_to_cells[i] = []
                self._node_to_cell[i].append(this_cell)

        self.created_cell(this_cell)
                
        return this_cell

    def edges_to_rings(self, edgemask=None, ccw=1):
        """ using only the edges for which edgemask is true,
        construct rings.  if edgemask is not given, use all of the
        current edges

        if ccw is 1, only non-intersecting ccw rings will be return
        if ccw is 0, only non-intersecting cw rings will be return
        """
        if edgemask is not None:
            edges = self.edges['nodes'][edgemask]
            masked_grid = TriGrid(points=self.nodes['x'],edges=edges['nodes'] )
            return masked_grid.edges_to_rings(edgemask=None)

        # remember which edges have already been assigned to a ring
        edges_used = zeros( self.Nedges(), int8 )

        rings = []
        
        for start_e in range(self.Nedges()):
            if edges_used[start_e]:
                continue

            # start tracing with the given edge -
            # it's hard to know beforehand which side of this edge is facing into
            # the domain, so start with the assumption that it obeys our convention
            # that going from edge[i,0] to edge[i,1] the interior is to the left
            # once a ring has been constructed, check to see if it has negative area
            # in which case we repeat the process with the opposite ordering.

            # one problem, though, is that interior rings really should have negative
            # area, since they will become holes.

            # at least the one with the largest area is correct.
            # Then any that are inside it should have negative areas...
            
            # what if we found all rings with positive area and all with negative
            # area.  then we'd have all the information ready for choosing who is
            # inside whom, and which orientation is correct?

            failed_edges_used1 = None
            failed_edges_used2 = None
            
            for flip in [0,1]:
                e = start_e
                edges_used[e] = 1 # tentatively used.
                a,b = self.edges['nodes'][e]
                if flip:
                    a,b = b,a

                if self.verbose > 1:
                    print "Starting ring trace with nodes ",a,b
                ring = [a,b] # stores node indices

                node_count = 1
                while 1:
                    node_count += 1

                    # used to be node_count > self.Npoints(), but since we step
                    # one extra bit around the circle, then go back and remove one
                    # node, I think it should be 1+.
                    if node_count > 1+self.Npoints():
                        # debug
                        self.plot()
                        pnts = self.nodes['x'][ring] 
                        plot(pnts[:,0],pnts[:,1],'ro')
                        # /debug
                        raise Exception,"Traced too far.  Something is wrong.  bailing"
                    
                    b_edges = self.node_to_edges(b)

                    if len(b_edges) == 2:
                        # easy case - one other edge leaves.
                        # new edge e
                        if b_edges[0] == e:
                            e = b_edges[1]
                        else:
                            e = b_edges[0]
                    else:
                        # calculate angles for all the edges, CCW relative to the
                        # x-axis (atan2 convention)
                        angles = []
                        for next_e in b_edges:
                            c = setdiff1d(self.edges['nodes'][next_e],[b])[0]
                            d = self.nodes['x'][c] - self.nodes['x'][b]
                            angles.append( arctan2(d[1],d[0]) )
                        angles = array(angles)
                        
                        e_idx = b_edges.index(e)
                        e_angle = angles[e_idx]
                        angles = (angles-e_angle) % (2*pi)
                        next_idx = argsort(angles)[-1]
                        e = b_edges[next_idx]
                        
                    if self.edges['nodes'][e,0] == b:
                        c = self.edges['nodes'][e,1]
                    else:
                        c = self.edges['nodes'][e,0]

                    # now we have a new edge e, and the next node in the ring c
                    if edges_used[e] == 0:
                        edges_used[e] = 1 # mark as tentatively used
                    else:
                        # we guessed wrong, and now we should just bail and flip the
                        # other way but it's not so slow just to keep going and figure it
                        # out later.
                        # print "Could be smarter and abort now."
                        pass

                    if len(ring) >= 2 and b==ring[0] and c == ring[1]:
                        #print " %d,%d == %d,%d we've come full circle.  well done"%(b,c,
                        #                                                            ring[0],ring[1])
                        break
                    ring.append(c)
                    a,b = b,c

                # remove that last one where we figured out that we were really all the way
                # around.
                ring = ring[:-1]
                points = self.nodes['x'][ring]
                if bool(ccw) == bool(is_ccw(points)):
                    # print "great, got correctly oriented ring (ccw=%s)"%ccw
                    edges_used[ edges_used==1 ] = 2 # really used
                    rings.append( array(ring) )
                    break # breaks out of the flip loop
                else:
                    # print "ring orientation wrong, wanted ccw=%s"%ccw
                    if flip:
                        area = signed_area(points)
                        
                        if self.verbose > 1:
                            print "Failed to get positive area either way:"
                            print "Ring area is ",area

                        if isnan(area):
                            print "Got nan area:"
                            print points
                            raise Exception,"NaN area trying to figure out rings"
                        
                        # raise Exception,"Failed to make positive area ring in either direction"
                        # I think this is actually valid - when Angel Island gets joined to
                        # Tiburon, if you start on Angel island either way you go you trace
                        # a region CCW.
                        
                        # however, nodes that were visited in both directions
                        # should 'probably' be marked so we don't visit them more.
                        # really not sure how this will fair with a multiple-bowtie
                        # issue...
                        failed_edges_used2 = where(edges_used==1)[0]
                        
                        edges_used[ intersect1d( failed_edges_used1,
                                                 failed_edges_used2 ) ] = 2
                    
                    # otherwise try again going the other direction,
                    # unmark edges, but remember them
                    failed_edges_used1 = where(edges_used==1)[0]
                    edges_used[ edges_used==1 ] = 0 # back into the pool

        if self.verbose > 0:
            print "Done creating rings: %d rings in total"%len(rings)
        return rings

    def edges_to_polygons(self,edgemask):
        """ use the edges (possibly masked by given edgemask) to create
        a shapely.geometry.Polygon() for each top-level polygon, ordered
        by decreasing area
        """
        rings_and_holes = self.edges_to_rings_and_holes(edgemask)

        polys = []
        for r,inner_rings in rings_and_holes:
            outer_points = self.nodes['x'][r]
            inner_points = [self.nodes['x'][ir] for ir in inner_rings]
            polys.append( geometry.Polygon( outer_points, inner_points) )

        areas = array([p.area for p in polys])
        order = argsort(-1 * areas)

        return [polys[i] for i in order]
        
    def edges_to_rings_and_holes(self,edgemask):
        """ using only the edges for which edgemask is true,
        construct polygons with holes.  if edgemask is not given, use all of the
        current edges

        This calls edges_to_rings to get both ccw rings and cw rings, and
        then determines which rings are inside which.

        returns a list [ [ outer_ring_nodes, [inner_ring1,inner_ring2,...]], ... ]
        """
        if edgemask is not None:
            edges = self.edges[edgemask,:2]
            masked_grid = TriGrid(points=self.nodes['x'],edges=edges['nodes'])
            return masked_grid.edges_to_rings_and_holes(edgemask=None)

        # print "calling edges_to_rings (ccw)"
        ccw_rings = self.edges_to_rings(ccw=1)
        # print "calling edges_to_rings (cw)"
        cw_rings  = self.edges_to_rings(ccw=0)

        # print "constructing polygons"
        # make single-ring polygons out of each:
        ccw_polys = [geometry.Polygon(self.nodes['x'][r]) for r in ccw_rings]
        cw_polys  = [geometry.Polygon(self.nodes['x'][r]) for r in cw_rings]

        # assume that the ccw poly with the largest area is the keeper.
        # technically we should consider all ccw polys that are not inside
        # any other poly
        ccw_areas = [p.area for p in ccw_polys]

        outer_rings = outermost_rings( ccw_polys )

        # Then for each outer_ring, search for cw_polys that fit inside it.
        outer_polys = []  # each outer polygon, followed by a list of its holes

        # print "finding the nesting order of polygons"
        for oi in outer_rings:
            outer_poly = ccw_polys[oi]

            # all cw_polys that are contained by this outer ring.
            # This is where the predicate error is happening -
            possible_children_i = []
            for i in range(len(cw_polys)):
                try:
                    if i!=oi and outer_poly.contains( cw_polys[i] ):
                        if not cw_polys[i].contains(outer_poly):
                            possible_children_i.append(i)
                        else:
                            print "Whoa - narrowly escaped disaster with a congruent CW poly"
                except shapely.predicates.PredicateError:
                    print "Failed while comparing rings - try negative buffering"
                    d = sqrt(cw_polys[i].area)
                    inner_poly = cw_polys[i].buffer(-d*0.00001,4)

                    if outer_poly.contains( inner_poly ):
                        possible_children_i.append(i)

            # the original list comprehension, but doesn't handle degenerate
            # case
            
            # possible_children_i = [i for i in range(len(cw_polys)) \
            #                        if outer_poly.contains( cw_polys[i] ) and i!=oi ]
                
            possible_children_poly = [cw_polys[i] for i in possible_children_i]

            # of the possible children, only the ones that are inside another child are
            # really ours.  outermost_rings will return indices into possible_children, so remap
            # those back to proper cw_poly indices to get children.
            children = [possible_children_i[j] for j in outermost_rings( possible_children_poly )]

            outer_polys.append( [ccw_rings[oi],
                                 [cw_rings[i] for i in children]] )
            
        return outer_polys

    def trim_to_left(self, path):
        """ Given a path, trim all cells to the left of it.
        """
        # mark the cut edges:
        for i in range(len(path)-1):
            e = self.nodes_to_edge( path[i:i+2] )

            if self.edges['marker'][e] in [0,CUT_EDGE]:
                # record at least ones that are really cut, in case some of
                # of the cut edges are actually on the boundary
                cut_edge = (path[i],path[i+1],e)

                self.edges[e]['marker'] = CUT_EDGE

        # choose the first cell, based on the last edge that was touched above:

        # the actual points:
        a = self.nodes['x'][cut_edge[0]]
        b = self.nodes['x'][cut_edge[1]]
        # the edge index
        edge = cut_edge[2]

        # the two cells that form this edge:
        cell1,cell2 = self.edges['cells'][edge]
        if cell1 >= 0:
            other_point1 = setdiff1d( self.cells['nodes'][cell1], cut_edge[:2] )[0]
        else:
            other_point1 = None
        if cell2 >= 0:
            other_point2 = setdiff1d( self.cells['nodes'][cell2], cut_edge[:2] )[0]
        else:
            other_point2 = None

        parallel = (b-a)
        # manually rotate 90deg CCW
        bad = array([ -parallel[1],parallel[0]] )

        if dot(self.nodes['x'][other_point1],bad) > dot(self.nodes['x'][other_point2],bad):
            bad_cell = cell1
        else:
            bad_cell = cell2

        print "Deleting"
        self.recursive_delete(bad_cell)
        print "Renumbering"
        self.renumber()

    def recursive_delete(self,c,renumber = 1):
        del_count = 0
        to_delete = [c]

        # things the queue have not been processed at all...

        while len(to_delete) > 0:
            # grab somebody:
            c = to_delete.pop()
            if self.cells[c,0] == -1:
                continue

            # get their edges
            nodea,nodeb,nodec = self.cells['nodes'][c]

            my_edges = [self.nodes_to_edge( (nodea,nodeb) ),
                        self.nodes_to_edge( (nodeb,nodec) ),
                        self.nodes_to_edge( (nodec,nodea) ) ]

            # mark it deleted - should this really be using
            # the self.delete_cell() API?
            self.cells['stat'][c] = STAT_DEL
            del_count += 1

            # add their neighbors to the queue to be processed:
            for e in my_edges:
                if self.edges['marker'][e] == 0:# only on non-cut, internal edges:
                    c1,c2 = self.edges['cells'][e]
                    if c1 == c:
                        nbr = c2
                    else:
                        nbr = c1

                    if nbr >= 0:
                        to_delete.append(nbr)
        print "Deleted %i cells"%del_count


    def leaf_nodes(self,use_envelope=False,tolerance=0.0,use_degree=True):
        """  return the indexes nodes which are deemed 'leaves', by a combination
        of
         use_envelope: node falls on the rectangular envelope, i.e. self.bounds()
         use_degree: node is degree 1
        """
        leaves = []
        
        if use_envelope:
            xmin,xmax,ymin,ymax = self.bounds()
            xmin+=tolerance
            xmax-=tolerance
            ymin+=tolerance
            ymax-=tolerance

            valid = (self.nodes['x'][:,0] < xmin) | (self.nodes['x'][:,0]>xmax) | \
                    (self.nodes['x'][:,1] < ymin) | (self.nodes['x'][:,1]>ymax)
            valid = valid & (self.nodes['stat']==STAT_OK)

            node_iterable = nonzero(valid)[0]
        else:
            node_iterable = self.valid_node_iter()
            
        for n in node_iterable:
            if (not use_degree) or len(self.node_to_edges(n)) == 1:
                leaves.append(n)
        return leaves

    def stitch_in_grid(self,gridB,tolerance=0.01,use_B_envelope=False):
        """ Stitch another grid into this grid.
        gridB: another grid, which will be modified, and joined to self.

        Currently joining is only performed at leaf nodes - nodes which have exactly one
        edge.  This will fail when the two grid share edges, or the two grids together would
        form a threeway (or more) intersection.

        tolerance: the max distance between two leaf nodes for them to be considered for joining.
        use_B_envelope: if true, only points on the boundary of the envelope of gridB will be considered.
          potentially speeds up stitches of tiled grids. (envelope = rectilinear bounds).  This test
          will use the same tolerance as leaf-to-leaf comparisons.
        """
        # Aleaves = self.leaf_nodes()
        Bleaves = gridB.leaf_nodes(use_envelope=use_B_envelope,tolerance=tolerance)
        
        ##
        to_join = [] # [(node from A,node from B), ...]

        #for nA in Aleaves:
        for nB in Bleaves:
            nA = self.closest_node( gridB.nodes['x'][nB] )
            d = dist( self.nodes['x'][nA],
                      gridB.nodes['x'][nB] )
            if d < tolerance:
                print "Joining nodes with separation of ",d
                to_join.append( (nA,nB) )

        # prune the common nodes from B, and remember the current B index of the next nodes,
        # so that we can come back and
        leaf_edgesB = [gridB.node_to_edges(nB)[0] for nA,nB in to_join]

        # save all of the data from these edge - don't forget those elevations!
        leaf_edgesBdata = gridB.edges[leaf_edgesB].copy()

        # delete
        for (nA,nB) in to_join:
            gridB.delete_node(nB,remove_edges=1)

        # merge grid arrays - have to update index changes for elements of B
        # could probably call self.append_grid() here, but while code is in flux keep
        # separate to avoid confusion
        nB_offset = self.Npoints()
        jB_offset = self.Nedges()
        iB_offset = self.Ncells()
        gridB.edges['nodes'] += nB_offset
        gridB.edges['cells'] += iB_offset
        gridB.cells['nodes'] += nB_offset

        self.nodes = concatenate_safe_dtypes( (self.nodes,gridB.nodes) )
        self.edges = concatenate_safe_dtypes( (self.edges,gridB.edges) )
        self.cells = concatenate_safe_dtypes( (self.cells,gridB.cells) )
        ## end append_grid()
        
        # Now add back in those edges
        for join_i in range(len(to_join)):
            nA,nB = to_join[join_i]
            Bedgedata = leaf_edgesBdata[join_i]
            nB_nbr = setdiff1d(Bedgedata['nodes'],[nB])[0]

            j = self.add_edge(nA,nB_nbr+nB_offset)

            # copy back the elevation info
            for field in ['marker','ridge_z','region_z','relief']:
                self.edges[j][field] = Bedgedata[field]
        self.refresh_metadata()

    def append_grid(self,gridB):
        """ Add all elements of gridB to the current grid
        compared to stitch_in_grid, this should be much faster, but no attempt
        is made to check for duplicate node locations, edges, etc.,

        also caller is responsible for self.refresh_metadata()

        modifies gridB, shifting the indices, and saves the index shift
        in gridB.{node,edge,cell}_offset
        """
        # merge grid arrays - have to update index changes for elements of B
        gridB.node_offset = self.Nnodes()
        gridB.edge_offset = self.Nedges()
        gridB.cell_offset = self.Ncells()
        gridB.edges['nodes'] += gridB.node_offset
        gridB.edges['cells'] += gridB.cell_offset
        gridB.cells['nodes'] += gridB.node_offset

        if 0:
            self.nodes = concatenate_safe_dtypes( (self.nodes,gridB.nodes) )
            self.edges = concatenate_safe_dtypes( (self.edges,gridB.edges) )
            self.cells = concatenate_safe_dtypes( (self.cells,gridB.cells) )
        else:
            self.nodes = array_concatenate( (self.nodes,gridB.nodes) )
            self.edges = array_concatenate( (self.edges,gridB.edges) )
            self.cells = array_concatenate( (self.cells,gridB.cells) )

        # clear out any bad state:
        self._node_to_cells = None
        self._node_to_edges = None
        self.index = None
        self.edge_index = None
        self._calc_vcenters = False
        self._components = None

    @staticmethod
    def stitch_grids(grids,use_envelope=True,envelope_tol=0.01,join_tolerance=0.25):
        """ grids: an iterable of TriGrid2 instances
        combines all grids together, removing duplicate points, joining coincident vertices

        use_envelope: use the rectangular bounding box of each grid to determine joinable
          nodes.   
        envelope_tol: if using the grid bounds to determine joinable leaf nodes, the tolerance
          for determining that a node does lie on the boundary.
        join_tolerance: leaf nodes from adjacent grids within this distance range will be
          considered coincident, and joined.
          
        """
        # for each grid, an array of node indices which will be considered
        all_leaves = []

        accum_grid = None

        for i,gridB in enumerate(grids):
            if i % 100 == 0:
                print "%d / %d"%(i,len(grids)-1)

            if gridB.Npoints() == 0:
                print "empty"
                continue

            gridB.verbose = 0
            gridB.renumber()

            Bleaves = array(gridB.leaf_nodes(use_envelope=use_envelope,
                                             tolerance=envelope_tol,
                                             use_degree=False),
                            int32)

            if i == 0:
                accum_grid = gridB
                if len(Bleaves):
                    all_leaves.append( Bleaves )
            else:
                accum_grid.append_grid(gridB)
                if len(Bleaves):
                    all_leaves.append( Bleaves + gridB.node_offset)

        all_leaves = concatenate(all_leaves)

        # build an index of the leaf nodes to speed up joining
        lf = field.XYZField(X=accum_grid.nodes['x'][all_leaves],
                            F=all_leaves)
        lf.build_index()

        to_join=[] # [ (i,j), ...] , with i<j, and i,j indexes into accum_grid.nodes
        for i,l in enumerate(all_leaves):
            nbrs = lf.nearest(lf.X[i],count=4)

            for nbr in nbrs:
                if nbr <= i:
                    continue
                dist = norm(lf.X[i] - lf.X[nbr])
                if dist < join_tolerance:
                    print "Joining with distance ",dist
                    to_join.append( (all_leaves[i], all_leaves[nbr]) )

        # okay - so need to allow for multiple joins with a single node.
        # done - but is the joining code going to handle that okay?

        _remapped = {} # for joined nodes, track who they became
        def canonicalize(n): # recursively resolve remapped nodes
            while _remapped.has_key(n):
                n = _remapped[n]
            return n

        for a,b in to_join:
            a = canonicalize(a)
            b = canonicalize(b)
            if a==b:
                continue

            a_nbrs = accum_grid.node_neighbors(a)
            accum_grid.delete_node(a)
            for a_nbr in a_nbrs:
                # with edges along a boundary, it's possible that
                # the new edge already exists
                try:
                    accum_grid.nodes_to_edge( [a_nbr,b] )
                except NoSuchEdgeError:
                    accum_grid.add_edge(a_nbr,b)
        return accum_grid

    def clip_to_geom(self,geom):
        """ remove all nodes which do not intersect geom, and adjacent edges, too.
        returns the number of nodes deleted
        """
        nodes_to_delete = []
        for n in self.valid_node_iter():
            p = geometry.Point( self.nodes['x'][n] )
            if not geom.contains(p):
                nodes_to_delete.append( n )
        print "clipping %d nodes"%len(nodes_to_delete)
        for n in nodes_to_delete:
            self.delete_node(n)
        return len(nodes_to_delete)
            
    def crop(self,xxyy):
        """ Truncate a grid to the given bounding box.  Currently the bounds are
        inclusive, but edges will not be subdivided, so if an edge straddles the boundary
        it will be lost entirely, but if a node is on the boundary, it will remain
        """
        x = self.nodes['x'][:,0]
        y = self.nodes['x'][:,1]
        valid = self.nodes['stat'] == STAT_OK
        
        bad_nodes = valid & ( (x < xxyy[0]) | (x > xxyy[1]) | (y < xxyy[2]) | (y > xxyy[3]) )
        for i in nonzero(bad_nodes)[0]:
            self.delete_node(i)
            
    # Undo-history management - very generic. well, hopefully.
    op_stack_serial = 10
    op_stack = None
    def checkpoint(self):
        if self.op_stack is None:
            self.op_stack_serial += 1
            self.op_stack = []
        return self.op_stack_serial,len(self.op_stack)

    def revert(self,cp):
        serial,frame = cp
        if serial != self.op_stack_serial:
            raise ValueError,"The current op stack has serial %d, but your checkpoint is %s"%(self.op_stack_serial,
                                                                                              serial)
        while len(self.op_stack) > frame:
            self.pop_op()

    def commit(self):
        self.op_stack = None
        self.op_stack_serial += 1
    
    def push_op(self,meth,*data,**kwdata):
        if self.op_stack is not None:
            self.op_stack.append( (meth,data,kwdata) )

    def pop_op(self):
        f = self.op_stack.pop()
        if self.verbose > 3:
            print "popping: ",f
        meth = f[0]
        args = f[1]
        kwargs = f[2]
        
        meth(*args,**kwargs)
    ###
        
    def unmove_node(self,i,orig_val):
        # update point index:
        if self.index is not None:
            curr_coords = self.nodes['x'][i,xxyy]
            orig_coords = orig_val[xxyy]
            
            self.index.delete(i, curr_coords )
            self.index.insert(i, orig_coords )

        self.nodes['x'][i] = orig_val

    def move_node(self,i,new_pnt):
        self.push_op(self.unmove_node,i,self.nodes['x'][i].copy())
                     
        # update point index:
        if self.index is not None:
            old_coords = self.nodes['x'][i,xxyy]
            new_coords = new_pnt[xxyy]
            self.index.delete(i, old_coords )
            self.index.insert(i, new_coords )

        self.nodes['x'][i] = new_pnt

        self.updated_node(i)
        
        for e in self.node_to_edges(i):
            if self.edge_index is not None:
                self.edge_index.delete(e,self.edges[j]['center'][xxyy])

            self.update_edge_center(e)
            
            if self.edge_index is not None:
                self.edge_index.insert(e,self.edges[j]['center'][xxyy])
            
            self.updated_edge(e)

    def updated_node(self,i):
        for cb in self._update_node_listeners.values():
            cb(i)

    def updated_edge(self,e):
        for cb in self._update_edge_listeners.values():
            cb(e)

    def updated_cell(self,c):
        for cb in self._update_cell_listeners.values():
            cb(c)

    def created_node(self,i):
        for cb in self._create_node_listeners.values():
            cb(i)

    def created_edge(self,e):
        # fix up the edge index
        ec = self.nodes['x'][self.edges['nodes'][e]].mean(axis=0)
        
        self.edges['center'][e] = ec

        if self.edge_index is not None:
            self.edge_index.insert( e, ec[xxyy] )
            
        for cb in self._create_edge_listeners.values():
            cb(e)

    def created_cell(self,c):
        for cb in self._create_cell_listeners.values():
            cb(c)

    def deleted_cell(self,c):
        for cb in self._delete_cell_listeners.values():
            cb(c)

    def deleted_node(self,i):
        for cb in self._delete_node_listeners.values():
            cb(i)
            
    def deleted_edge(self,e):
        for cb in self._delete_edge_listeners.values():
            cb(e)
        

    # subscriber interface for updates:
    listener_count = 0
    def init_listeners(self):
        self._update_node_listeners = {}
        self._update_edge_listeners = {}
        self._update_cell_listeners = {}
        self._create_node_listeners = {}
        self._create_edge_listeners = {}
        self._create_cell_listeners = {}
        self._delete_node_listeners = {}
        self._delete_edge_listeners = {}
        self._delete_cell_listeners = {}
    
    def listen(self,event,cb):
        cb_id = self.listener_count
        if event == 'update_node':
            self._update_node_listeners[cb_id] = cb
        elif event == 'update_edge':
            self._update_edge_listeners[cb_id] = cb
        elif event == 'update_cell':
            self._update_cell_listeners[cb_id] = cb
        elif event == 'create_node':
            self._create_node_listeners[cb_id] = cb
        elif event == 'create_edge':
            self._create_edge_listeners[cb_id] = cb
        elif event == 'delete_node':
            self._delete_node_listeners[cb_id] = cb
        elif event == 'delete_edge':
            self._delete_edge_listeners[cb_id] = cb
        elif event == 'delete_cell':
            self._delete_cell_listeners[cb_id] = cb
        else:
            raise Exception,"unknown event %s"%event
            
        self.listener_count += 1
        return cb_id

    def unlisten(self,cb_id):
        for l in [ self._update_node_listeners,
                   self._update_edge_listeners,
                   self._update_cell_listeners,
                   self._create_node_listeners,
                   self._create_edge_listeners,
                   self._create_cell_listeners,
                   self._delete_node_listeners,
                   self._delete_edge_listeners,
                   self._delete_cell_listeners]:
            if l.has_key(cb_id):
                del l[cb_id]
                return
        print "Failed to remove cb_id %d"%cb_id

    ## Native I/O - pickling 
    def save(self,fn):
        """ Save the complete (hopefully) state of the Paving to a file.
        The hope is that this will be useful for checkpointing a paving
        and restarting it later.
        The first attempt here is to use the Pickle protocol, and take care of
        all the gritty details in __getstate__ and __setstate__
        """
        fp = file(fn,'wb')
        pickle.dump(self,fp,-1)
        fp.close()

    @staticmethod
    def load(fn):
        fp = file(fn,'rb')
        obj = pickle.load(fp)
        fp.close()

        # seems that there have been some issues with _node_to_edges
        obj._node_to_edges = None
        
        return obj

    def __getstate__(self):
        d = self.__dict__.copy()

        # Clear out things that we don't need
        d['_create_node_listeners'] = {}
        d['_create_cell_listeners'] = {}
        d['_create_edge_listeners'] = {}
        d['_update_node_listeners'] = {}
        d['_update_edge_listeners'] = {}
        d['_delete_edge_listeners'] = {}
        d['_delete_cell_listeners'] = {}
        d['_delete_node_listeners'] = {}
        d['op_stack_serial'] = 10
        d['op_stack'] = None

        #PAVER# d['last_fill_iter'] = None

        # probably ought to be in live_dt
        #LIVEDT# d['DT'] = 'rebuild'
        #LIVEDT# d['vh'] = 'rebuild'

        #PAVER# d['poly'] = 'rebuild'
        
        #PAVER# d['plot_history'] = []
        #PAVER# d['click_handler_id'] = None
        #PAVER# d['showing_history'] = None
        d['edge_collection'] = None
        
        d['edge_inedx'] = None
        d['index'] = None
        
        return d
    
    def __setstate__(self,d):
        self.__dict__.update(d)

        #PAVER# # probably ought to be in live_dt
        #PAVER# try:
        #PAVER#     if len(self.shapely_rings) > 0:
        #PAVER#         self.poly = geo.Polygon( self.shapely_rings[0],self.shapely_rings[1:] )
        #PAVER#     else:
        #PAVER#         self.poly = None
        #PAVER# except AttributeError:
        #PAVER#     self.shapely_rings = []
        #PAVER#     self.poly = None
            
        self.refresh_metadata()
        
