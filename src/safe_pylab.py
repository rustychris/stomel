import os

if not os.environ.has_key('DISPLAY'):
    import matplotlib
    matplotlib.use('Agg')

try:
    from pylab import *
except RuntimeError:
    pass
except ImportError:
    pass
