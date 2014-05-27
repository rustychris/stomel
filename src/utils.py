import functools
import numpy as np
from collections import OrderedDict

def add_to(instance):
    def decorator(f):
        import types
        f = types.MethodType(f, instance, instance.__class__)
        setattr(instance, f.func_name, f)
        return f
    return decorator

class Bucket(object):
    def __init__(self,**kwargs):
        self.__dict__.update(kwargs)

def records_to_array(records):
    # Convert that to an array
    rectype = []
    if len(records) == 0:
        recarray = np.zeros(0)
    else:        
        for k in records[0].keys():
            if k=='date':
                t=object
            elif k in ['inst','line']:
                t=np.int32
            else:
                t=np.float64
            rectype.append((k,t))

        recarray = np.zeros(len(records),dtype=rectype)
        for i,rec in enumerate(records):
            for k,v in rec.iteritems():
                recarray[i][k] = v
    return recarray

def mag(vec):
    vec = np.asarray(vec)
    return np.sqrt( (vec**2).sum(axis=-1))

def center_to_interval(c):
    d=np.ones_like(c)
    d[1:-1] = abs(0.5*(c[2:] - c[:-2]))
    d[0]=d[1] ; d[-1]=d[-2]
    return d

def center_to_edge(c):
    d=np.ones(len(c)+1)
    d[1:-1] = 0.5*(c[1:] - c[:-1])
    d[0]=2*d[1] - d[2] # d[1] - (d[2]-d[1])
    d[-1]=2*d[-2]-d[-3]
    return d
    
class BruteKDE(object):
    def __init__(self,values,weights,bw):
        self.values=values
        self.weights=weights
        self.bw=bw
        self.norm_factor=np.sum(self.weights)*np.sqrt(np.pi)*bw
    def __call__(self,x):
        res=np.zeros_like(x)
        for idx in np.ndindex(res.shape):
            res[idx]=np.sum( self.weights*np.exp( -((self.values-x[idx])/self.bw)**2 ) )
        return res/self.norm_factor

def quantize(a,stride,axis=0,reducer=np.mean):
    # first truncate to an even multiple of stride
    N=(a.shape[axis]//stride)*stride
    slices=[ slice(None) ]*a.ndim
    slices[axis]=slice(N)
    a=a[slices]
    dims=list(a.shape)
    dims[axis:axis+1] = [N//stride,stride]
    return reducer( a.reshape(dims), axis=axis+1)
