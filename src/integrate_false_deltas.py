""" Code to compute a distribution of cell depths in a false delta, attempting to retain
the hypsometry of the overall region.


Inputs
1  grid
2  bathymetry on the grid
3  polygons defining the false deltas
4  polygons defining the regions the false deltas represent
5  DEM which covers the region to be integrated

In general, 3 and 4 are the same set of polygons.

Output
  bathymetry with cells in the false delta region updated to
  reflect the integrated bathy
"""


from osgeo import ogr

import field
import trigrid
from shapely import wkb
from shapely import geometry

## Try prepared geometries:
from shapely.geometry import Point
from shapely.prepared import prep

from numpy import *

class FalseDelta(object):
    """ Represents a single false delta region
    """
    def __init__(self,name,geom):
        self.name = name
        self.geom = geom
        
    def integrate(self,dem,max_elevation=inf):
        # the total dem - 3.4M valid cells
        xmin,ymin,xmax,ymax = self.geom.bounds
        cropped_dem = dem.crop( [xmin,xmax,ymin,ymax] )

        # in the cropped DEM, down to 400k cells.
        x,y = cropped_dem.xy()

        valid = ~isnan(cropped_dem.F)
        
        if isfinite(max_elevation):
            valid = valid & (cropped_dem.F < max_elevation)
            
        i,j = nonzero( valid )
        
        prepared_polygon = prep(self.geom)
        # Seems that on MacOS the prepared geometry doesn't always work.
        # prepared_polygon = self.geom

        hits = []
        for k in range(len(i)):
            if k % 10000 == 0:
                print "%d/%d"%(k,len(i))

            point = geometry.Point( x[j[k]],y[i[k]] )
            # If this line hangs, you might have a bum version of libgeos.
            # Try upgrading.
            if prepared_polygon.contains(point):
                depth = cropped_dem.F[i[k],j[k]]
                hits.append( depth )

        hits = array(hits)
        Acell = dem.dx * dem.dy
        total_area = Acell * len(hits)
        
        self.hits = hits
        self.Acell = Acell
        self.total_area = total_area

    def fill_in_depths_on_grid(self,g,new_depths,offset=0.0):
        """ For cells of the given that fall within the integration polygon for this
        false delta, overwrite their bathymetry in new_depths [Nc] by ordering the
        DEM hits deepest to shallowest, and assigning those depths nearest to farthest
        based on straight-line distance from where the grid intersects the integration
        polygon
        """
        # which cells fall in the integration polygon for this delta?
        vc = g.vcenters()

        prepared_polygon = prep(self.geom)

        cells = []
        cells_in_poly = {}
        for i in range(g.Ncells()):
            if i % 10000 == 0:
                print "%d/%d [offset=%f]"%(i,g.Ncells(),offset)

            point = geometry.Point( vc[i,0], vc[i,1] )
            if prepared_polygon.contains(point):
                cells.append(i)
                cells_in_poly[i] = 1

        cells = array(cells)

        # find where the false channel exists the integration polygon -
        outside_neighbor = None
        for i in cells:
            # look for a cell with a neighbor that isn't in the polygon
            for nbr in g.cell_neighbors(i,adjacent_only=1):
                if not cells_in_poly.has_key(nbr):
                    outside_neighbor = nbr
                    break
            if outside_neighbor is not None:
                break

        # order all the hits according to distance from that point -
        p0 = vc[outside_neighbor]
        d = sum( (p0 - vc[cells])**2, axis=1)
        ordering = argsort(d)
        ordered_cells = cells[ordering]

        # make an array of depths, deep to shallow, and respective areas
        self.hits.sort()

        ####
        bin_i = 0 # the next depth bin with some area to contribute
        depth_bins = array( [self.hits, self.Acell*ones_like(self.hits)]).T

        big_area = 1e10
        depth_bins = concatenate( (depth_bins,[[depth_bins[-1,0],big_area]]),axis=0 )

        # the depth of each of the false delta cells
        areas = g.areas()
        cell_area_sum = areas[ordered_cells].sum()

        for c in ordered_cells:
            Ac = areas[c]
            Ac_remaining = Ac
            new_depths[c] = 0.0
            while Ac_remaining>0:
                if Ac_remaining < depth_bins[bin_i,1]:
                    area_to_add = Ac_remaining
                    Ac_remaining=0.0
                    depth_bins[bin_i,1] -= area_to_add
                else:
                    area_to_add = depth_bins[bin_i,1]
                    depth_bins[bin_i,1] = 0.0
                    Ac_remaining-=area_to_add

                # negate here, because depth_bins has elevations, while new_depths
                # has soundings (typically...)
                new_depths[c] += area_to_add/Ac * (-depth_bins[bin_i,0])

                if depth_bins[bin_i,1] == 0.0:
                    bin_i += 1
                    
                if isnan(new_depths[c]):
                    raise Exception,"Somehow cell %d has depth ",new_depths[c]
            # and apply the offset -
            new_depths[c] -= offset

        # How much area was taken from the sentinel depth bin?    
        print "Area taken from fake last bin: %f [%f%% of grid area] at depth %f"%( big_area - depth_bins[-1,1],
                                                                100* (big_area - depth_bins[-1,1]) / cell_area_sum,
                                                                depth_bins[-1,0])
        print "Area remaining in integration of DEM: %f [%f%% of grid area] with min elevation %f"%( depth_bins[:-1,1].sum(),
                                                                                 100 * depth_bins[:-1,1].sum() / cell_area_sum,
                                                                                 depth_bins[bin_i,0] )
        # save the cells that were updated so we can later fix up edge depths.
        self.cells = cells
        
        
class FalseDeltas(object):
    """ Extract stats/hypsometry for False Deltas, given DEM
    and polygon shapefile
    """
    def __init__(self,integration_poly_shp,bathy_file):
        self.integration_poly_shp = integration_poly_shp
        self.bathy_file = bathy_file
        
        self.deltas = self.read_shps(self.integration_poly_shp)

    def read_shps(self,integration_poly_shp):
        ods = ogr.Open(integration_poly_shp)
        layer = ods.GetLayer(0)

        deltas = []
        while 1:
            feat = layer.GetNextFeature()
            if feat is None:
                break
            name = feat.GetField('name')
            geom = wkb.loads( feat.GetGeometryRef().ExportToWkb() )
            deltas.append(FalseDelta(name=name,geom=geom))
        return deltas

    def integrate(self):
        print "Loading DEM"
        dem = field.GdalGrid(self.bathy_file)

        for d in self.deltas:
            print "integrating one delta..."
            d.integrate(dem)

    def fill_in_depths_on_grid(self,g,new_depths,offset):
        for fd in self.deltas:
            fd.fill_in_depths_on_grid(g,new_depths,offset)
            
    def report(self):
        for d in self.deltas:
            print "%s: total area: %f  length at 500m wide: %f"%(d.name,
                                                                 d.total_area,
                                                                 d.total_area / 500 )





class FalseDeltaGenerator(object):
    def __init__(self,
                 delta_poly_shp, # shapefile with polygon defining the false deltas
                 dem_fn,         # a Gdal file giving the bathymetry in the region
                 depths,         # the current [Nc,3] cell depths - will be modified!
                 grid,           # trigrid instance
                 offset=0.0,     # constant to be added to all bathymetry (probably want -5)
                 # If specified, these will be updated where the false delta affects them
                 edge_depths=None,
                 depths_from_edges=None):
                 
        self.delta_shp = delta_poly_shp
        self.dem_fn = dem_fn
        self.depths = depths
        self.grid = grid
        self.offset = offset
        self.edge_depths = edge_depths
        self.depths_from_edges = depths_from_edges
        
    def process(self):
        print "Integrating volume of the 'real' Delta"
        self.fds = FalseDeltas(self.delta_shp,self.dem_fn)
        self.fds.integrate()

        self.fds.fill_in_depths_on_grid(self.grid,self.depths[:,2],offset = self.offset)

        # Remaining: update edge depths for edges within the false delta, then recreate the cell centered
        # depth from edges
        if self.edge_depths is not None:
            for fd in self.fds.deltas:
                cells = fd.cells
                edges_to_update = zeros( self.grid.Nedges(), bool )

                for c in cells:
                    e = self.grid.cell2edges(c)
                    edges_to_update[e] = True
                nc = self.grid.edges[edges_to_update,3:5]    
                self.edge_depths[edges_to_update,2] = self.depths[nc,2].min(axis=1)

            if self.depths_from_edges is not None:
                # And recalculate the cell depths from edges:
                celldepths = self.depths[:,2].copy()
                for i in range(self.grid.Ncells()):
                    jlist = self.grid.cell2edges(i)
                    celldepths[i] = self.edge_depths[jlist,2].max()
                self.depths_from_edges[:,2] = celldepths


## Useful methods for splicing grids together:
# Get those ordered, so the open area is on the left of the node chain
def order_along_boundary(g,bnodes,internal_to="left"):
    """ g: paver.Paving or Trigrid (probably okay)
    bnodes: array of contiguous nodes on the boundary.

    returns: reordered bnodes in order as they appear on the boundary, with the
    internal/paved side of the polyline to the given direction
    """

    print "Top of order_along_boundary: bnodes=",bnodes
    
    # Find the end points - 
    # presumably all of these nodes are connected to each other, and if we
    # consider only edges which are unpaved on one side, then there should
    # be only one polyline connecting them all.
    end_nodes = []

    for n in bnodes:
        edges_this_node = 0
        for j in g.pnt2edges(n):
            if all(g.edges[j,3:]>=0):
                continue
            if g.edges[j,0] in bnodes and g.edges[j,1] in bnodes:
                edges_this_node += 1
        if edges_this_node == 1:
            end_nodes.append(n)
        elif edges_this_node == 0:
            print "ERROR: node %d didn't find any edges"%n
        elif edges_this_node == 2:
            pass
        else:
            print "ERROR: node %d has too many connection at the edge"%n


    trav = end_nodes[0]
    chained = [trav]
    while 1:
        print "Checking from node %d"%trav
        for j in g.pnt2edges(trav):
            # restrict ourselves to edges joining two of the boundary
            # nodes, which are open on one side.
            print "Checking node %d-%d"%(g.edges[j,0],g.edges[j,1])
            if all(g.edges[j,3:]>=0): # internal edge
                print "Internal edge - ",g.edges[j,3:]
                other = None
                continue
            if g.edges[j,0] == trav:
                other = g.edges[j,1]
            else:
                other = g.edges[j,0]
            print "other node is ",other
            if other in bnodes and other not in chained:
                break
            other = None
        if other is not None:
            trav = other
            chained.append(other)
        else:
            break
    chained = array(chained)
    
    # check the internal side based on the first edge:
    j = g.find_edge( (chained[0],chained[1]) )
    if g.edges[j,0] == chained[0]:
        j_dir = 1 # it's "forward" w.r.t. the chain
    else:
        j_dir = -1 # backwards
    if internal_to=='left':
        req = 1
    else:
        req = -1
    if g.edges[j,3] >= 0:
        open_side = 1
    else:
        open_side = -1

    if j_dir * req * open_side < 0:
        chained = chained[::-1]
    
    return chained
                



# ################################################################################
# integration_poly_shp = "false_delta_integration_polys-v001.shp"
# bathy = "/home/rusty/classes/research/spatialdata/us/ca/suntans/bathymetry/usgs/northbay/delta_navd88.tif"
# depths = loadtxt("depth.dat")
# grid_path = "."
# g = trigrid.TriGrid(suntans_path=grid_path)        


