import getfem
from numpy import *

mfu=getfem.MeshFem('load','tank_3D.mfu')
m=mfu.linked_mesh()
mfp=getfem.MeshFem('load','tank_3D.mfp',m)
U = fromfile('tank_3D.U', 'd')
P = fromfile('tank_3D.P', 'd')

sl=getfem.Slice(('boundary',('intersection',('planar',+1,[0,0,0],[0,1,0]),('planar',+1,[0,0,0],[1,0,0]))),m,3);

import locale # (bug #13014)
print "importing tvtk.."
import getfem_tvtk
print "import done"
locale.setlocale(locale.LC_NUMERIC,('en_US','UTF8')) # (bug #13014)

fig = getfem_tvtk.Figure(gui='tvtk')

fig.show(sl, data=(mfp, P), vdata=(mfu,U), edges=False)

fig.show(sl, data=(mfp, P), edges=False)


old=fig.scalar_range()

sl=getfem.Slice(('boundary',('intersection',('planar',+1,[0,0,6],[0,0,-1]),
                             ('planar',+1,[0,0,0],[0,1,0]))),m,3);
fig.show(sl, data=(mfp, P), scalar_bar=True, edges=False)
fig.scalar_range((-40,40));

#print fig.scalar_range()

m.set_region(42, m.outer_faces());
m.region_substract(42, 3);

sl2=getfem.Slice(('boundary',('planar',+1,[0,0,0],[0,1,0])),m,3, m.region(42))
fig.show(sl2, faces=False, edges_color=(1,1,1))


hh=array([[1, 5, 9, 12.5, 16, 19.5]]);
H=concatenate((zeros((2,hh.size)), hh));

tsl=getfem.Slice('streamlines', mfu, U, H);

fig.show(tsl, tube_color=(1,1,1));

fig.loop()

