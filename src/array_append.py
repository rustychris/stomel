""" Utility functions for manipulating arrays
"""
from numpy import zeros,ndarray
import numpy

def array_append( A, b ):
    """ append b to A, where b.shape == A.shape[1:]
    Attempts to make this fast by dynamically resizing the base array of
    A, and returning the appropriate slice.
    """

    # a bit more complicated because A may have a different column ordering
    # than A.base (due to a previous transpose, probably)
    # can compare strides to see if our orderings are the same.

    # possible that it's not a view, or
    # the base array isn't big enough, or
    # the layout is different and it would just get confusing, or
    # A is a slice on other dimensions, too, which gets too confusing.
    
    if A.base is None or type(A.base) == str \
           or A.base.size == A.size or A.base.strides != A.strides \
           or A.shape[1:] != A.base.shape[1:]:
        new_shape = list(A.shape)

        # make it twice as long as A, and in case the old shape was 0, add 10
        # in for good measure.
        new_shape[0] = new_shape[0]*2 + 10  
        
        base = zeros( new_shape, dtype=A.dtype)
        base[:len(A)] = A

        # print "resized based array to %d elements"%(len(base))
    else:
        base = A.base

    A = base[:len(A)+1]
    if A.dtype.isbuiltin:
        A[-1] = b
    else:
        if type(b) == numpy.void:
            A[-1] = b.tolist()
        else:
            A[-1] = list(b) # was .tolist(), but I think this is more permissive
    return A


def array_concatenate( AB ):
    """
    similiar to array_append, but B.shape[1:] == A.shape[1:]

    while the calling convention is similar to concatenate, it currently only supports
    2 arrays
    """
    A,B = AB
    
    if A.base is None or type(A.base) == str \
           or A.base.size == A.size or A.base.strides != A.strides \
           or len(A) + len(B) > len(A.base) \
           or A.shape[1:] != A.base.shape[1:]:
        new_shape = list(A.shape)

        # make it twice as long as A, and in case the old shape was 0, add 10
        # in for good measure.
        new_shape[0] = max(new_shape[0]*2 + 10,
                           new_shape[0] + len(B))
        
        base = zeros( new_shape, dtype=A.dtype)
        base[:len(A)] = A
    else:
        base = A.base

    lenA = len(A)
    A = base[:lenA+len(B)]
    A[lenA:] = B
        
    return A
    
    
def concatenate_safe_dtypes( ab ):
    """ Concatenate two arrays, but allow for the dtypes to be different.  The
    fields are taken from the first array - matching fields in subsequent arrays
    are copied, others discarded.
    """
    a,b = ab # for now, force just two arrays
    
    result = zeros( len(a)+len(b), a.dtype)
    result[:len(a)] = a

    for name in b.dtype.names:
        if name in a.dtype.names:
            result[ len(a):len(a)+len(b)][name] = b[name]
    return result

def recarray_add_fields(A,new_fields):
    """ A: a record array
    new_fields: tuples of ('name',data)
    where data must be the same length as A.  So far, no support for
    non-scalar values
    """
    new_dtype=A.dtype.descr
    for name,val in new_fields:
        # handle non-scalar fields
        # assume that the first dimension is the "record" dimension
        new_dtype.append( (name,val.dtype,val.shape[1:] ) )
    new_names=[name for name,val in new_fields]
    new_values=[val for name,val in new_fields]
    new_A=numpy.zeros( len(A), dtype=new_dtype)
    
    for name in new_A.dtype.names:
        try:
            new_A[name]=new_values[new_names.index(name)]
        except ValueError:
            new_A[name]=A[name]

    return new_A
