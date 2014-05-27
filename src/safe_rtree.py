import sys

Rtree = None
RtreeNative = None
RtreeQgis = None

# newer Qgis doesn't interfere with Rtree -
if Rtree is None:
    try:
        from rtree.index import Rtree as RtreeNative
        Rtree = RtreeNative
    except ImportError:
        pass

if Rtree is None:
    try:
        # okay - we might be forced to use qgis:
        from qgis_spatialindex import RtreeQgis
        Rtree = RtreeQgis
    except ImportError:
        pass

if Rtree is None:
    try:
        from kdtree_spatialindex import RtreeKDTree
        Rtree = RtreeKDTree
    except ImportError:
        print "Failed to find any spatial index implementations"
        raise
    
# if we are in Qgis - try to check for whether Rtree is corrupted
if sys.modules.has_key('qgis') and RtreeNative is not None:
    from qgis_spatialindex import RtreeQgis
    # Now test whether RtreeNative is functional -
    import rtree
    p = rtree.index.Property()
    p.set_storage(0)
    if p.get_storage() != 0:
        print "Running under Qgis, and Rtree is corrupt - will use Qgis SpatialIndex"
        Rtree = RtreeQgis
    else:
        print "Running under Qgis, but Rtree /appears/ safe..."
