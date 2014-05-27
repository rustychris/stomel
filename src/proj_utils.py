"""
Bundle up common projection operations
"""
from osgeo import osr
from numpy import *

def to_srs(s):
    if not isinstance(s,osr.SpatialReference):
        srs = osr.SpatialReference() ; srs.SetFromUserInput(s)
        return srs
    else:
        return s
    
def xform(src,dest):
    src = to_srs(src)
    dest = to_srs(dest)
    return osr.CoordinateTransformation(src,dest)

def mapper(src,dest):
    trans = xform(src,dest)
    def process(X):
        X = asarray(X,float64)
        out = zeros(X.shape,float64)
        inx = X[...,0]
        iny = X[...,1]

        # careful - these have to remain as views
        outx = out[...,0]
        outy = out[...,1]

        # for idx in range(len(inx)):
        for idx in ndindex(inx.shape):
            outx[idx],outy[idx],dummy = trans.TransformPoint(inx[idx],iny[idx],0)
        return out
    return process

#m = mapper("WGS84","EPSG:26910")
#print m(array([-122,35]))
