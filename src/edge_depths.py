# Some utilities for dealing with specifying edge depths:
import getopt
import sys,os
sys.path.append(os.path.join(os.environ['HOME'], 'python'))

from numpy import *

import sunreader
#reload(sunreader)

import pdb

class EdgeDepthWriter(object):
    """  For a grid that has been partitioned, write out new edge-based bathymetry
    and update the cell bathymetry accordingly.

    If a depth field is specified to run, it will be queried (via value_on_edge) for
    the edge depths.  

    If no depth field is given, but the sunreader instance has a global edgedepths.dat,
    it will be mapped to the per-processor edges.
    
    Otherwise the previous cell depths will be interpolated to get
    the new edge depths.

    Either way, new cell depths are taken as the minimum depth of adjacent edges.
    """
    
    def __init__(self,delete_friction_strips=False):
        """ delete_friction_strips: if true, a final step will lower edges to the elevation
        of the higher of the two cell neighbors. Note that this is done *after* the
        cells were set to the lowest neighboring edge, so the overall effect is to
        increase flux areas.

        by default, the edges will have already been set to be the elevation of the
        shallower neighboring cell, based on cell depths.  That tends to smear out
        edge depths, and on average make the flux areas too small.
        """
        self.delete_friction_strips = delete_friction_strips

    def run(self,argv=[],sun=None,depth_field=None):
        """ Specify one of sun, or command-line style argv
        at the moment only a single command-line arguments is supported, specifying the
        datadir to read from

        """
        if sun is not None:
            self.sun = sun
        else:
            ## Handle command line
            datadir = '.'
            opts,rest = getopt.getopt(argv,'')

            for opt,val in opts:
                pass

            if len(rest):
                datadir = rest[0]

            ## Prep 
            self.sun = sunreader.SunReader(datadir)

        try:
            self.sun.file_path('edgedepths',0)
        except:
            print "Maybe you need to set edgedepths in suntans.dat?"
            raise

        ## Do it.
        self.process(depth_field=depth_field)

    def process(self,depth_field=None):
        np = self.sun.num_processors()

        global_edge_depths = None
        global_grid = None
        
        ## Establish the source for edge depths
        if depth_field is None:
            fn = self.sun.file_path('edgedepths')
            if os.path.exists( fn ):
                print "Loading global edge depths and mapping global->local"
                global_edge_depths = loadtxt(fn)
                global_grid = self.sun.grid()
            else:
                print "Looking for edge depth file %s"%fn
                print "Will resort to linear interpolation of cell depths for edge depth"
        else:
            print "Will use explicitly given depth field for edge elevations"
            

        # will have to hold all data at one time.. 
        proc_edgedepths = [None] * np
        proc_edgedata   = [None] * np
        proc_celldata   = [None] * np
        
        # 1 - get edge elevations and Nke.
        for proc in range(self.sun.num_processors()):
            edgedata = self.sun.edgedata(proc)
            celldata = self.sun.celldata(proc)

            g = self.sun.grid(proc)
            
            edgedepths = zeros( len(edgedata), float64 )
            # Save for inter-proc use:
            proc_edgedata[proc] = edgedata
            proc_celldata[proc] = celldata
            proc_edgedepths[proc] = edgedepths

            # if depth_field exists, query directly for edge depths
            # if global_edge_depths exists, take pre-calculated edge depths from there
            # otherwise, interpolate between neighboring cells.
            for j in range(len(edgedata)):
                if depth_field is not None:
                    de = depth_field.value_on_edge( g.points[g.edges[j,:2]] )
                elif global_edge_depths is not None:
                    global_j = global_grid.find_edge( g.edges[j,:2] )
                    de = global_edge_depths[global_j,2]
                else:                
                    nc1,nc2 = edgedata[j,8:10].astype(int32)
                    face1,face2 = edgedata[j,10:12].astype(int32)
                    if nc2 < 0:
                        nc2 = nc1
                        face2 = face1
                    elif nc1 < 0:
                        nc1 = nc2
                        face1 = face2

                    df1 = celldata[nc1,14+face1]
                    df2 = celldata[nc2,14+face2]

                    # linearly interpolate
                    de = (df1*celldata[nc2,3] + df2*celldata[nc1,3]) / (df1+df2)

                edgedepths[j] = de
            ## Update Nke
            # de are as soundings.
            # h_to_ctop expects elevations, and by default uses dzmin to get the
            # surface behavior.
            # add 1 because Nke is the *number* of levels, whereas h_to_ctop gives
            # the *index* of the level where this elevation belongs.
            # for example, if the z-levels are every 0.5m, and the edgedepth is exactly 4.5,
            # this should give us Nke =
            # Add the epsilon to make sure that any roundoff from the interpolation above doesn't
            # screw up an exact comparison

            offenders = nonzero(edgedepths > self.sun.z_levels()[-1])[0]
            # allow a bit of leeway in case of ascii roundoff - fix it regardless, but only
            # report the issue if it's significant.
            bad_offenders = nonzero(edgedepths > 1e-5 + self.sun.z_levels()[-1])[0]
            if len(offenders) > 0:
                if len(bad_offenders) > 0:
                    print "Bottom of lowest z-level is %f"%self.sun.z_levels()[-1]
                    print "There were %d edges given a depth below this"%len(offenders)
                    print "And %d of those are significant "%len(bad_offenders)
                    print "too deep by these distances: "
                    print edgedepths[bad_offenders] - self.sun.z_levels()[-1]
                    # raise Exception,"Whoa there - edges are too deep!"
                    print "WARNING: these edges will have their depth truncated. "
                    print "         to avoid this, specify a vertspace.dat.in that goes"
                    print "         deep enough"
                edgedepths[offenders] = self.sun.z_levels()[-1]
            edgedata[:,6] = searchsorted(self.sun.z_levels()+1e-8,edgedepths) + 1

            # double check
            nkmax = self.sun.conf_int('Nkmax')
            if nkmax>1 and any(edgedata[:,6]>nkmax):
                raise Exception,"How did a deep edge get through?"
            
        # Interprocessor - only exchange Nke from edgedata:
        self.sun.sendrecv_edges([a[:,6] for a in proc_edgedata])
        self.sun.sendrecv_edges(proc_edgedepths)

        ## Set cell depth, Nk from deepest edge
        for proc in range(self.sun.num_processors()):
            edgedata = proc_edgedata[proc] 
            celldata = proc_celldata[proc] 
            edgedepths = proc_edgedepths[proc]
            
            # Update cell depths and Nk[] from the edges
            for i in range(len(celldata)):
                js = celldata[i,5:8].astype(int32)

                # depth = max(depth)
                celldata[i,3] = edgedepths[js].max()
                # Nk = max(Nke)
                celldata[i,4] = edgedata[js,6].max()
            
            # And update the edge's Nkc
            for j in range(len(edgedata)):
                nc1,nc2 = edgedata[j,8:10].astype(int32)
                if nc2 < 0:
                    nc2 = nc1
                elif nc1 < 0:
                    nc1 = nc2

                edgedata[j,7] = max( celldata[nc1,4],
                                     celldata[nc2,4] )

        # Now Nkc on marker 6 edges is corrupt - fix it
        # those edges could only see one cell locally, but the neighbor
        # proc has both cells, so get Nkc from the neighbor
        self.sun.sendrecv_edges([a[:,7] for a in proc_edgedata])

        ## Optionally lower high edges to their higher cell neighbor.
        if self.delete_friction_strips:
            print "Removing friction strips by lowering edge depths"
            for proc in range(self.sun.num_processors()):
                edgedata = proc_edgedata[proc] 
                celldata = proc_celldata[proc] 
                edgedepths = proc_edgedepths[proc]
                
                for j in range(len(edgedata)):
                    nc1,nc2 = edgedata[j,8:10].astype(int32)
                    face1,face2 = edgedata[j,10:12].astype(int32)
                    if nc2 < 0:
                        nc2 = nc1
                        face2 = face1
                    elif nc1 < 0:
                        nc1 = nc2
                        face1 = face2

                    df1 = celldata[nc1,14+face1]
                    df2 = celldata[nc2,14+face2]

                    # take the shallower cell
                    edgedepths[j] = min(celldata[nc1,3],celldata[nc2,3])
                    edgedata[j,6] = min(celldata[nc1,4],celldata[nc2,4])
            # And we still have to fix those interprocessor edges:
            self.sun.sendrecv_edges([a[:,6] for a in proc_edgedata]) # Nke
            self.sun.sendrecv_edges(proc_edgedepths) # de
            
        ## Write it out
        for proc in range(self.sun.num_processors()):
            edgedata = proc_edgedata[proc] 
            celldata = proc_celldata[proc] 
            edgedepths = proc_edgedepths[proc]

            # WRITE IT OUT
            # maybe a little funky - sunreader keeps a reference to edgedata,
            # and we have been passing this reference around the whole time -
            # so no need to pass back in the modified edgedata.
            # remove them first to avoid overwriting symlinked copies:
            os.unlink(self.sun.file_path('celldata',proc))
            os.unlink(self.sun.file_path('edgedata',proc))
            self.sun.write_edgedata(proc)
            self.sun.write_celldata(proc)

            fp = open(self.sun.file_path('edgedepths',proc),'wt')
            for j in range(len(edgedepths)):
                fp.write("%.6f %.6f %.6f\n"%(edgedata[j,4], edgedata[j,5], edgedepths[j]))
            fp.close()


if __name__ == '__main__':        
    edw = EdgeDepthWriter()
    edw.run()

