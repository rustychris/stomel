# Node data fields:
STAT = 0
ORIG_RING = 1 # only set for nodes on the original_ring
ALPHA = 2
BETA = 3

# Node statuses:
HINT = 555  # node is on boundary, but just a hint of where the boundary is
SLIDE = 666  #  attached, but allowed to move along boundary
FREE =  777  # internal node, can move in 2D
RIGID = 888 #  may not be moved
DELETED = 999

# Beta rules:
BETA_NEVER=0
BETA_RESCUE=1
BETA_ALWAYS=2

# Nonlocal methods:
PROACTIVE_DELAUNAY=1
SHEEPISH=2
