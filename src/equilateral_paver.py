# Make a grid with all equilateral triangles
# Currently only supports a rectangular domain, constant density,
# and either vertical or horizontal orientation

import trigrid
import numpy as np

class EquilateralPaver(trigrid.TriGrid):
    def __init__(self,L,W,dens,orientation='horizontal',**kwargs):
        super(EquilateralPaver,self).__init__(**kwargs)
        self.L = L # x dimension
        self.W = W # y dimension
        self.dens = dens
        self.orientation = orientation

        if self.orientation == 'vertical':
            self.L,self.W = self.W,self.L
        self.create_grid()
        if self.orientation == 'vertical':
            self.L,self.W = self.W,self.L
            self.points = self.points[:,::-1]
            self.cells = self.cells[:,::-1]

        self.renumber()

    def create_grid(self):
        # first, how many rows - here we assume orientation is horizontal,
        # so the left and right sides are ragged.
        cos30 = np.cos(30*np.pi/180.)
        
        n_rows = self.W / (cos30 * self.dens)
        # to make sure that the first and last points line up, we need an
        # even number of rows of cells:
        n_rows = 2 * int( (n_rows+1.0)/ 2 )
        
        self.n_rows = n_rows
        
        # Let the length L be fudge-able - as in we prefer perfectly equilateral triangles
        # over a perfectly L-length grid.  the width W can still be exact.
        dens = self.W / (n_rows * cos30)

        print "That will make n_rows=%d and adjusted edge length %f"%(n_rows,dens)

        # this is the number of cells...
        n_cols = int(self.L / dens)

        self.n_cols = n_cols

        # Stack them up
        for r in range(n_rows+1):
            y = self.W * float(r)/n_rows
            odd = r%2
            x_off = odd * 0.5*dens
                
            for c in range(n_cols+1):
                x = x_off + dens*float(c)
                
                n = self.add_node( np.array([x,y]) )

                if c > 0:
                    if r==0:
                        self.add_edge(n-1,n,cright=-1,marker=1)
                    elif r==n_rows:
                        self.add_edge(n-1,n,cleft=-1,marker=1)
                    else:
                        self.add_edge(n,n-1)
                    
                # HERE: need to finish adding in the markers and closed boundary code.
                if r>0:
                    cright=-2
                    cleft=-2
                    marker = 0
                    
                    if odd:
                        if c==0:
                            cleft=-1
                            marker=1
                        elif c==n_cols:
                            cright=-1
                            marker=1
                            
                        self.add_edge(n-(n_cols+1),n,marker=marker,cleft=cleft,cright=cright)
                            
                        if c<n_cols:
                            self.add_edge(n,n-n_cols)
                    else:
                        if c==0:
                            cleft=-1
                            marker=1
                        elif c==n_cols:
                            cright=-1
                            marker=1
                        self.add_edge(n-(n_cols+1),n,cleft=cleft,cright=cright,marker=marker)
                        
                        if c>0:
                            self.add_edge(n,n-(n_cols+1)-1)
                

class RotatedEquilateralPaver(EquilateralPaver):
    """ Create a ragged-edged grid where the triangles are rotated the given
    angle, in radians, CCW from parallel to the x-axis.  
    """
    def __init__(self,L,W,dens,angle=0,**kwargs):
        self.final_L = L
        self.final_W = W

        # find the L and W needed to still be big enough after we've rotated -
        # adding a bit of extra to avoid funny edge effects:
        Lprime = L*np.cos(angle) + W*np.sin(angle) + 4*dens
        Wprime = W*np.cos(angle) + L*np.sin(angle) + 4*dens
        
        super(RotatedEquilateralPaver,self).__init__(L=Lprime, W=Wprime, dens=dens, **kwargs)

        self.rotate_grid(angle)
        self.trim_grid()
        self.renumber()

    def rotate_grid(self,angle):
        """ rotates the oversized grid and translates to get the origin in the right place.
        """
        # translate to get centered on the extra bit we asked for:
        self.points[:] -= 2*self.dens
        # rotate
        self.points[:] = trigrid.rot(angle,self.points)
        # and get our origin to a nice place
        self.points[:,0] += self.final_L * np.sin(angle)**2
        self.points[:,1] -= self.final_L * np.sin(angle)*np.cos(angle)

    def trim_grid(self):
        """ with the oversized grid created, and the origin correctly placed, remove points
        and associated edges/cells that fall outside the actual footprint
        """
        to_delete = (self.points[:,0] < 0) | (self.points[:,0]>self.final_L) | \
                    (self.points[:,1] < 0) | (self.points[:,1]>self.final_W)

        for n in np.nonzero(to_delete)[0]:
            self.delete_node(n,remove_edges=True)

            

if __name__ == '__main__':
    #ep = EquilateralPaver(10000.,5000.,500.,orientation='horizontal')
    #ep.plot()
    ep = RotatedEquilateralPaver(10000.,5000.,510.,angle=15*pi/180.)
    cla()
    ep.plot()
