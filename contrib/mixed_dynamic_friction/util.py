from getfem import *
from numpy import *

with_graphics=True
try:
    import getfem_tvtk
except:
    print "\n** Could NOT import getfem_tvtk -- graphical output disabled **\n"
    import time
    time.sleep(2)
    with_graphics=False

print 'Some tests with python',


m=Mesh('load', '../../tests/meshes/disc_P2_h4.mesh')


if with_graphics:
  fig = getfem_tvtk.Figure()
  fig.show(m)
  print "Press Q to continue.."
  fig.set_colormap('tripod')
  fig.loop()
