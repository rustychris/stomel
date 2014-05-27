# create bathymetry for the false delta that is more specific to the real delta's
# geometry

from numpy import *
from pylab import *
import field

    
# load the real delta bathymetry
delta = field.GdalGrid('/home/rusty/classes/research/spatialdata/us/ca/suntans/bathymetry/usgs/northbay/delta_v2.asc')
# the file is in tenths of feet, with -9999 meaning no data, and 100 meaning inside leveed pond
delta.F[ delta.F == 100 ] = -9999

# cheat a little - instead of constructing the full polygon
# for the partitions, we get the same result by just cropping
# at easting = 609137.

crop = delta.bounds()
crop[0] = 609137

delta_cropped = delta.crop(crop)

xyz = delta_cropped.to_xyz()

# now convert tenths of a foot to m (and float64)
xyz.F = xyz.F * 0.03048

import navd88
xyz.F = navd88.convert_ngvd29_to_navd88(xyz.X[:,0],
                                        xyz.X[:,1],
                                        xyz.F)


min_depth = floor(xyz.F.min())
max_depth = ceil( xyz.F.max() )

dA = 100 # 10m x 10m
depths = arange(min_depth,max_depth,0.05)
areas = zeros_like(depths)
vols  = zeros_like(depths)


for i in range(len(depths)):
    wet = xyz.F < depths[i]

    areas[i] = dA * sum(wet)
    vols[i] = (dA * (depths[i] - xyz.F[wet])).sum()

    print "%g: %g %g"%(depths[i],areas[i],vols[i])


if 0:
    clf()
    subplot(2,1,1)
    plot(depths,areas)
    subplot(2,1,2)
    plot(depths,vols)

# load the false deltas:
import trigrid

delta1grid = trigrid.TriGrid(sms_fname='/home/rusty/classes/research/meshing/delta_1.grd')
delta2grid = trigrid.TriGrid(sms_fname='/home/rusty/classes/research/meshing/delta_2.grd')


dpoints = []
delevs  = []

total_false_area = abs(delta1grid.areas()).sum() + abs(delta2grid.areas()).sum()

depth_factor = areas[-1] / total_false_area
print "Multiplying depths by: ",depth_factor

for dgrid in [delta1grid,delta2grid]:
    # they're tilings, but go ahead and take the mean area:
    dA = abs(dgrid.areas()).mean() # the area of one triangle

    totalA = abs(dgrid.areas()).sum()

    # normalized areas, such 
    norm_areas = arange(1,dgrid.Ncells+1)*dA/totalA * areas[-1]

    # Algorithm:
    #   the tiling has to quantize area into one-triangle increments
    #   so we add in a triangle, which increases the total area by dA
    #  so fraction_areas[i] is the fraction of total delta area that
    #   is covered by the first i triangles.  use

    grid_depths = zeros_like(norm_areas)

    for i in range(len(norm_areas)):
        # linearly interpolate to get the estimated depth at
        # which the area becomes this large:
        grid_depths[i] = depth_factor * interp(norm_areas[i],areas,depths,right=depths[-1])

        # bug in interp:
        if isnan(grid_depths[i]):
            grid_depths[i] = depth_factor * interp(0.0001+norm_areas[i],areas,depths,
                                                   right=depths[-1],left=depths[0])


    # come up with an ordering of the cells, based on vcenter locations:
    vcenters = dgrid.vcenters()
    midY = mean(vcenters[:,1])
    leftX = vcenters[:,0].min()

    # it's 85km long, 1km wide, 
    skew_dist = sqrt(  (vcenters[:,0] - leftX)**2 + 80*80*(vcenters[:,1] - midY)**2 )
    ordering = argsort(skew_dist)


    delta_points = vcenters
    delta_elevs  = zeros( delta_points.shape[0], float64 )
    delta_elevs[ordering] = grid_depths

    dpoints.append(delta_points)
    delevs.append(delta_elevs)


# combine the new depths into a XYZ field
dpoints = concatenate(dpoints)
delevs  = concatenate(delevs)

false_field = field.XYZField(dpoints,delevs)

false_field.plot()
false_field.write('/home/rusty/classes/research/spatialdata/us/ca/suntans/bathymetry/compiled2/false_delta.bin')

    



