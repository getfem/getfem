""".
This demo use levelset to impose (weakly) a Dirichlet condition on an
implicit boundary defined by the zero of the levelset
"""
import getfem as gf
import numpy as np
from scipy import rand,setdiff1d

NX=40
ls_degree = 2

m = gf.Mesh('cartesian', np.arange(-.5,.5+1.0/NX,1.0/NX),np.arange(-.5,.5+1.0/NX,1.0/NX))
ls = gf.LevelSet(m, ls_degree)
ls2 = gf.LevelSet(m, ls_degree, 'with_secondary')

mf_ls = ls.mf()
mf_ls2 = ls2.mf()

P = mf_ls.basic_dof_nodes()
x = P[0,:]
y = P[1,:]

ULS = 1000*np.ones((1,x.size))


if 0:
  for ix in xrange(5):
    for iy in xrange(5):
      xc = (ix/4) * 0.8 - 0.4
      yc = (iy/4) * 0.8 - 0.4
      if iy%2==0:
        xc = xc + 0.05
      else:
        xc = xc - 0.05
      R = 0.03 + 0.005*(iy-1)
      ULS = np.minimum(ULS, ((x - xc)**2 + (y - yc)**2) - R**2)
else:
  for i in xrange(8):
    xc = rand() - 0.5
    yc = rand() - 0.5
    R = rand() * 0.09 + 0.02
    ULS = np.minimum(ULS, ((x - xc)**2 + (y - yc)**2) - R**2);

ls.set_values(ULS)

ULS2 = 1000*np.ones((1,x.size));
ULS2s = 1000*np.ones((1,x.size));

for i in xrange(1):
  xc = 0 # rand() - 0.5
  yc = 0 # rand() - 0.5
  theta = np.pi/3 #np.pi*rand()
  n = np.array([-np.sin(theta), np.cos(theta)])

  R = 0.19 #rand() * 0.09 + 0.02
  ULS2 = np.minimum(ULS2, ((x-xc)*n[0] + (y-yc)*n[1]))
  #ULS2s = np.minimum(ULS2s, ((x - xc).^2 + (y - yc).^2) - R^2)
  ULS2s = np.minimum(ULS2s, (abs(y - yc)+abs(x-xc) - R))

ls2.set_values(ULS2,ULS2s) # '-y-x+.2') # '(y-.2)**2 - 0.04')

mls = gf.MeshLevelSet(m)
mls.add(ls)
mls.add(ls2)
mls.adapt()

mim_bound = gf.MeshIm('levelset',mls,'boundary(a+b)', gf.Integ('IM_TRIANGLE(6)')) #, gf.Integ('IM_QUAD(5)'))
mim = gf.MeshIm('levelset',mls,'all(a+b)', gf.Integ('IM_TRIANGLE(6)'))
mim.set_integ(4)

mfu0 = gf.MeshFem(m,2)
mfu0.set_fem(gf.Fem('FEM_QK(2,3)'))

mfdu = gf.MeshFem(m,1)
mfdu.set_fem(gf.Fem('FEM_QK_DISCONTINUOUS(2,2)'))

mf_mult = gf.MeshFem(m,2)
mf_mult.set_fem(gf.Fem('FEM_QK(2,1)'))

A = gf.asm_volumic('V()+=comp()',mim_bound)

mls.cut_mesh().export_to_pos("mls.pos")
mf_ls.export_to_pos("ULS.pos",ULS,'uls')

dof_out = mfu0.dof_from_im(mim)
cv_out = mim.convex_index()
cv_in = setdiff1d(m.cvid(),cv_out)

#mfu = gf.MeshFem('partial', mfu0, dof_out, cv_in)

mfu0.export_to_pos('mesh.pos')

md = gf.Model('real')
md.add_fem_variable('u',mfu0)
md.add_initialized_data('lambda', [1])
md.add_initialized_data('mu', [1])
md.add_isotropic_linearized_elasticity_brick(mim, 'u', 'lambda', 'mu')

md.solve()

""".
U = md.variable('u')

VM = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u', 'lambda', 'mu', mfdu)

mfdu.export_to_pos('sol.pos',VM,'Von Mises or Tresca',U,'deformation')

mf_ls.export_to_pos('LS.pos',ls.values(0),'ls values')
mf_ls2.export_to_pos('LS2.pos',ls2.values(0),'ls2 values')
"""
