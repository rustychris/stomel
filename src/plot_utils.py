from safe_pylab import *
import time
from matplotlib.collections import LineCollection
from matplotlib import pyplot as pl

# convenience function for getting coordinates from the plot:
def pick_points(n):
    count = [0]
    pick_points.results = zeros( (n,2), float64)

    fig = gcf()
    cid = None
    
    
    def click_handler(event):
        print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
            event.button, event.x, event.y, event.xdata, event.ydata)
        if event.xdata:
            pick_points.results[count[0]] = [event.xdata, event.ydata]
            count[0] += 1
            if count[0] >= n:
                fig.canvas.mpl_disconnect(cid)
            

    cid = fig.canvas.mpl_connect('button_press_event', click_handler)

# A rehash of pick_points:
from array_append import array_append
def ax_picker(ax):
    fig = ax.figure
    
    if hasattr(ax,'pick_cids'):
        for cid in ax.pick_cids:
            fig.canvas.mpl_disconnect(cid)

    def init_picked():
        ax.picked = zeros( (0,4), float64)
        ax.pick_start = None
    init_picked()
    
    def on_press(event):
        if fig.canvas.toolbar.mode != '':
            return
        if event.button==1 and event.xdata:
            ax.pick_start = [event.xdata,event.ydata]
        elif event.button==3:
            print ax.picked
            init_picked()
    def on_release(event):
        if fig.canvas.toolbar.mode != '':
            return
        if event.xdata and ax.pick_start is not None:
            new_pnt = array([ax.pick_start[0],event.xdata,ax.pick_start[1],event.ydata])
            ax.picked=array_append( ax.picked,new_pnt )

    cid_p = fig.canvas.mpl_connect('button_press_event', on_press)
    cid_r = fig.canvas.mpl_connect('button_release_event', on_release)
    ax.pick_cids = [cid_p,cid_r]

    

def plotyy( x1, y1, x2, y2, color1='b', color2='g', fun=None, **kwargs ):
    """
    A work-alike of the Matlab (TM) function of the same name.  This
    places two curves on the same axes using the same x-axis, but
    different y-axes.

    Call signature::

    ax, h1, h2 = plotyy( x1, y2, x2, y2, color1='b', color2='g',
                         fun=None, **kwargs )

    color1 and color2 are the colors to make respective curves and y-axes.

    fun is the function object to use for plotting.  Must accept calls
    of the form fun(x,y,color='color',**kwargs).  Typically, something
    like plot, semilogy, semilogx or loglog.  If *None*, defaults to
    pyplot.plot.

    **kwargs is any list of keyword arguments accepted by fun.

    ax is a 2 element list with the handles for the first and second
    axes.  h1 is the handle to the first curve, h2 to the second
    curve.

    NOTE that this function won't scale two curves so that y-ticks are
    in the same location as the Matlab (TM) version does.
    """
    if fun == None: fun = plot

    ax1 = gca()
    ax1.clear()

    # Get axes location
    try:
        rect = ax1.get_position().bounds
    except AttributeError:
        rect = array( ax1.get_position() )
        rect[2:] += rect[:2]

    # Add first curve
    h1 = fun( x1, y1, color=color1, **kwargs )

    # Add second axes on top of first with joined x-axis
    ax2 = twinx(ax1)

    # Plot second curve initially
    h2 = fun( x2, y2, color=color2, **kwargs )

    # Set axis properties
    setp( ax2.get_xticklabels(), visible=False)

    # Change colors appropriately
    def recolor( obj, col ):
        try: obj.set_color( col )
        except: pass
        try: obj.set_facecolor( col )
        except: pass
        try: obj.set_edgecolor( col )
        except: pass
        try:
            ch = obj.get_children()
            for c in ch:
                recolor( c, col )
        except: pass

    recolor( ax1.yaxis, color1 )
    recolor( ax2.yaxis, color2 )

    draw_if_interactive()

    return ( [ax1,ax2], h1, h2 ) 


if 0:
    def dump_figure(fname,fig=None):
        """ Save as much of the info about a figure as possible, such that
        it can be replayed later, without dependence on the original data
        but with zoom-able interactions.
        """
        if fig is None:
            fig = gcf()


    # testing for dump_figure:
    figure(1)


        



# remove parts of the plot that extend beyond the x limits of the
# axis - assumes that the x-data for each line is non-decreasing
def trim_xaxis(ax=None):
    if ax is None:
        ax = gca()
        
    xmin,xmax,ymin,ymax = ax.axis()
    for line in ax.lines:
        xdata = line.get_xdata()
        ydata = line.get_ydata()

        i_start = searchsorted(xdata,xmin) - 1
        if i_start < 0:
            i_start = 0
        i_end = searchsorted(xdata,xmax) + 1

        xdata = xdata[i_start:i_end]
        ydata = ydata[i_start:i_end]

        line.set_xdata(xdata)
        line.set_ydata(ydata)
        

def plot_tri(tri,**kwargs):
    # DEPRECATED: matplotlib now has triplot and friends
    # compile list of edges, then create the collection, and plot
    
    ex = tri.x[tri.edge_db]
    ey = tri.y[tri.edge_db]
    edges = concatenate( (ex[:,:,newaxis], ey[:,:,newaxis]), axis=2)
    colors = ones( (len(edges),4), float32 )
    colors[:,:3] = 0
    colors[:,3] = 1.0
    
    coll = LineCollection(edges,colors=colors)

    ax = gca()
    ax.add_collection(coll)


def scalebar(xy,L=None,aspect=0.05,unit_factor=1,fmt="%.0f",label_txt=None,fractions=[0,0.5,1.0],
             ax=None):
    """ Draw a simple scale bar with labels - bottom left
    is given by xy.  
    """
    ax = ax or gca()
    
    if L is None:
        xmin,xmax,ymin,ymax = ax.axis()
        L = 0.2 * (xmax - xmin)
    xmin,ymin = xy
    dx = L
    dy = aspect*L
    # xmax = xmin + L
    ymax = ymin + dy

    objs = []
    txts = []
    objs.append( ax.fill([xmin,xmin+dx,xmin+dx,xmin],
                         [ymin,ymin,ymax,ymax],
                         'k', edgecolor='k') )
    objs.append( ax.fill([xmin,xmin+0.5*dx,xmin+0.5*dx,xmin],
                         [ymin,ymin,ymax,ymax],
                         'w', edgecolor='k') )
    for frac in fractions:
        txts.append( ax.annotate(fmt%(unit_factor* frac*L),
                                 [xmin+frac*dx,ymax + 0.25*dy],
                                 ha='center' )
                     )
    # annotate(fmt%(unit_factor*L), [xmin+dx,ymax+0.25*dy], ha='center')
    
    if label_txt:
        txts.append( ax.annotate(label_txt,[xmin+1.08*dx,ymin],ha='left') )
    return objs,txts
        

def show_slopes(ax=None,slopes=[-5./3,-1],xfac=5,yfac=3):
    ax = ax or pl.gca()
    x = np.median( [l.get_xdata()[-1] for l in ax.lines] )
    y = np.max( [l.get_ydata()[-1] for l in ax.lines] )

    y *= yfac # set the legend above the plotted lines

    xs = np.array([x/xfac,x])

    for s in slopes:
        ys = np.array([y/xfac**s,y])
        ax.loglog(xs,ys,c='0.5')
        pl.annotate("%g"%s,[xs[0],ys[0]])


        
# interactive log-log slope widget:
class Sloper(object):
    def __init__(self,ax=None,slope=-5./3,xfac=5,yfac=3,xlog=True,ylog=True,x=None,y=None):
        self.slope = slope
        self.ax = ax or pl.gca()
        if x is None:
            x = np.median( [l.get_xdata()[-1] for l in self.ax.lines] )
        if y is None:
            y = np.max( [l.get_ydata()[-1] for l in self.ax.lines] )

        y *= yfac # set the legend above the plotted lines

        self.xlog = xlog
        self.ylog = ylog

        xs = np.array([x/xfac,x])
        ys = np.array([y/xfac**slope,y])

        if self.xlog and self.ylog:
            self.line = self.ax.loglog(xs,ys,c='0.5',picker=5)[0]
        elif not self.xlog and not self.ylog:
            self.line = self.ax.plot(xs,ys,c='0.5',picker=5)[0]

        self.text = self.ax.text(xs[0],1.5*ys[0],"%g"%self.slope,transform=self.ax.transData)
        
        self.ax.figure.canvas.mpl_connect('pick_event',self.onpick)

        self.drag = dict(cid=None,x=None,y=None)

        
    def onpick(self,event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        print 'onpick points:', zip(xdata[ind], ydata[ind])
        print ' mouse point: ', event.mouseevent.xdata,event.mouseevent.ydata

        cid = self.ax.figure.canvas.mpl_connect('button_release_event',self.drag_end)

        if self.drag['cid'] is not None:
            self.ax.figure.canvas.mpl_disconnect(self.drag['cid'])
            
        self.drag = dict(cid=cid,x=event.mouseevent.xdata,y=event.mouseevent.ydata)

    yoff = 1.5
    def update_text_pos(self):
        x = self.line.get_xdata()[0]
        y = self.line.get_ydata()[0]
        self.text.set_x(x)
        if self.ylog:
            self.text.set_y(self.yoff*y)
        else:
            self.text.set_y(self.yoff+y)
        
    def drag_end(self,event):
        print "drag end"
        self.ax.figure.canvas.mpl_disconnect(self.drag['cid'])
        xdata = self.line.get_xdata()
        ydata = self.line.get_ydata()
        
        if self.xlog:
            xdata *= event.xdata / self.drag['x']
        else:
            xdata += (event.xdata - self.drag['x'])
        if self.ylog:
            ydata *= event.ydata / self.drag['y']
        else:
            ydata += event.ydata - self.drag['y']
        self.line.set_xdata(xdata)
        self.line.set_ydata(ydata)
        self.update_text_pos()
        event.canvas.draw()


class LogLogSlopeGrid(object):
    """ draw evenly spaced lines, for now in log-log space, at a given slope.
    y=mx+b
    """
    def __init__(self,ax=None,slopes=[-5/3.],intervals=[10]):
        """ Note that intervals is linear!  
        """
        self.ax = ax or pl.gca()
        self.slopes = slopes
        self.intervals = intervals
        self.colls = []
        self.xlog = self.ylog = True
        self.draw()
    def draw(self):
        for c in self.colls:
            self.ax.collections.remove(c)
        self.colls = []
        xmin,xmax,ymin,ymax = self.ax.axis()
        if self.xlog:
            xmin = np.log(xmin) ; xmax = np.log(xmax)
        if self.ylog:
            ymin = np.log(ymin) ; ymax = np.log(ymax)
            
        for s,interval in zip(self.slopes,self.intervals):
            corners = np.array( [[xmin,ymin],
                                 [xmax,ymin],
                                 [xmax,ymax],
                                 [xmin,ymax]] )
            corner_b = corners[:,1] - s*corners[:,0]

            if self.ylog:
                interval = np.log(interval)
            all_b = np.arange(corner_b.min(),corner_b.max(),interval)
            segs = np.zeros( (len(all_b),2,2), np.float64)
            segs[:,0,0] = xmin
            segs[:,1,0] = xmax
            segs[:,0,1] = s*xmin+all_b
            segs[:,1,1] = s*xmax+all_b
            if self.xlog:
                segs[...,0] = np.exp(segs[...,0])
            if self.ylog:
                segs[...,1] = np.exp(segs[...,1])
                
            coll = LineCollection(segs,color='0.75',zorder=-10)
            self.ax.add_collection(coll)
            self.colls.append(coll)

            
def enable_picker(coll,ax=None,cb=None):
    """ minimal wrapper for selecting indexes from a collection, like a
    scatter plot.  cb gets the first index chosen, and it returns an
    object which when called always returns the most recent index chosen
    """
    ax = ax or pl.gca()
    coll.set_picker(5) # should be 5 points

    class dummy(object):
        idx = None
        def __call__(self):
            return self.idx
    my_dummy = dummy()
    def onpick(event):
        if event.artist == coll:
            idx = event.ind[0]
            my_dummy.idx = idx
            if cb:
                cb(idx)
            else:
                print "Picked: I=%d"%idx

    my_dummy.cid = ax.figure.canvas.mpl_connect('pick_event',onpick)
    return my_dummy

def function_contours(f=lambda x,y: x-y,ax=None,Nx=20,Ny=20,V=10,
                      fmt=None):
    """ Cheap way to draw contours of a function and label them.
    Just evaluates the function on a grid and calls contour
    """
    ax = ax or pl.gca()
    xxyy = ax.axis()
    x = np.linspace(xxyy[0],xxyy[1],Nx)
    y = np.linspace(xxyy[2],xxyy[3],Ny)
    X,Y = np.meshgrid(x,y)
    Z = f(X,Y)
    ctr = pl.contour(X,Y,Z,V,colors='k')
    if fmt:
        ax.clabel(ctr,fmt=fmt)
        
    return ctr

