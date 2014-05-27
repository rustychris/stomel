# trying to make the code for generating depth.dat, edgedepth.dat, depth_from_edges.dat
# more independent

import sys,os

import sunreader, domain, trigrid
from numpy import *


from pylab import *
from gmtColormap import load_gradient
from matplotlib import colors


def main(args,depth_field_thunk):
    grid_dir = sys.argv[1]
    depths_onto_grid(grid_dir=grid_dir,
                     depth_field_thunk = depth_field_thunk,
                     output_dir=grid_dir,
                     do_plot=True)

def depths_onto_grid(grid_dir,depth_field_thunk,output_dir=None,do_plot=False):
    gridpath = grid_dir
    output_dir = output_dir or grid_dir
    
    g = trigrid.TriGrid(suntans_path=gridpath)

    depth_field = depth_field_thunk()
    
    bathy_offset = sunreader.read_bathymetry_offset()

    vc = g.vcenters()

    elevations = depth_field(vc)

    depths = bathy_offset - elevations

    # combine to one array:
    xyz = concatenate( (vc,depths[:,newaxis]), axis=1)
    savetxt(os.path.join(output_dir,'depth.dat'),xyz)

    # try putting depths on edges, then back-calculate cell depths:
    edgedepths = zeros( (g.Nedges(),3), float64)
    edgedepths[:,:2] = g.edge_centers()

    for j in range(g.Nedges()):
        if j % 1000==0:
            print "%d/%d"%(j,g.Nedges())
        edgedepths[j,2] = bathy_offset - depth_field.value_on_edge(g.points[g.edges[j,:2]])
    savetxt(os.path.join(output_dir,'edgedepth.dat'),edgedepths)    

    # and write a new cell depths file which takes its depth from the edges.
    celldepths = xyz.copy()
    for i in range(g.Ncells()):
        jlist = g.cell2edges(i)
        celldepths[i,2] = edgedepths[jlist,2].max()
    savetxt(os.path.join(output_dir,'depth-fromedges.dat'),celldepths)

    delta = celldepths[:,2] - xyz[:,2]

    if do_plot:
        if 1:
            depths = celldepths
            lognorm = colors.LogNorm(vmin=max(1.0,nanmin(depths[:,2])),
                                     vmax=nanmax( depths[:,2] ) )
            # cmap = gmtcm("BkBlAqGrYeOrReViWh200")
            cmap = load_gradient("oc-sst.cpt")

            figure(figsize=(100,100))
            ax = axes([0,0,1,1])
            ax.set_axis_bgcolor('gray')

            coll = g.plot_scalar(depths[:,2])
            coll.norm = lognorm
            coll.set_cmap(cmap)

            ax.xaxis.set_visible(0)
            ax.yaxis.set_visible(0)

            cax = axes([0.87,0.2,0.03,0.6])
            cbar = colorbar(coll,cax=cax)
        if 1:
            ticks = [2,5,10,20,50,100,200,500]
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(ticks)
            ax.axis('equal')
            savefig(os.path.join(output_dir,'grid_and_bathy.pdf'))


