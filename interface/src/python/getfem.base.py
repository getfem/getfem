#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
#
#
# Getfem-Python interface
#
# Date : March, 2004.
# Author : Julien Pommier, pommier@gmm.insa-tlse.fr
#                                                       
# This file is a part of GETFEM++                                         
#                                                                         
# Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 2.1 of the License,  or
# (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program;  if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
#                                                                         
"""Getfem-interface classes.

Provides access to the pseudo-objects exported by the getfem-python interface.
"""


__version__ = "$Revision$"
# $Source: /var/lib/cvs/getfem_matlab/src/python/getfem.base.py,v $

import sys
#import Numeric
import numarray

sys.path.insert(0,'build/lib.linux-i686-2.3/');
from _getfem import *
obj_count = {}
getfem('workspace','clear all')

def generic_constructor(self,clname,*args):
    """Internal function -- acts as a constructor for all getfem objects."""
#    print 'generic_constructor.'+clname+'('+str(args)+')'
    if (len(args)==1 and type(args[0]) is GetfemObject):
        self.id = args[0]
    else:
        self.id = getfem_from_constructor(clname,*args)
    obj_count[self.id] = obj_count.get(self.id,0)+1

def generic_destructor(self,destructible=True):
    """Internal function -- acts as a destructor for all getfem objects."""
    if (not hasattr(self,'id')):
        return
#    print "Mesh.__del__       ",self.id,'count=',obj_count[self.id]
    if (obj_count.has_key(self.id)):
        obj_count[self.id] = obj_count[self.id]-1
        if (destructible and obj_count[self.id] == 0):
#            print "effective deletion"
            getfem('delete',self.id)

# stub classes for getfem-interface objects

class Mesh:
    """Class for getfem mesh objects."""
    def __init__(self, *args):
        """General constructor for Mesh objects.

        @INIT MESH:INIT ('empty')
        @INIT MESH:INIT ('cartesian')
        @INIT MESH:INIT ('regular simplices')
        @INIT MESH:INIT ('triangles grid')
        @INIT MESH:INIT ('curved')
        @INIT MESH:INIT ('prismatic')
        @INIT MESH:INIT ('pt2D')
        @INIT MESH:INIT ('ptND')
        @INIT MESH:INIT ('load')
        @INIT MESH:INIT ('from string')
        @INIT MESH:INIT ('import')
        @INIT MESH:INIT ('clone')
        """
        generic_constructor(self,'mesh',*args)
    def __del__(self):
        generic_destructor(self,destructible=True)
    def get(self, *args):
        return getfem('mesh_get',self.id, *args)
    def set(self, *args):
        return getfem('mesh_set',self.id, *args)
    def __str__(self):
        return self.char()
    def __repr__(self):
        return '<getfem.Mesh %dD, %d points, %d convexes, %d bytes>' % \
               (self.dim(),self.nbpts(),self.nbcvs(),self.memsize())
    #@RDATTR MESH:GET('dim')
    #@GET    MESH:GET('pts')
    #@RDATTR MESH:GET('nbpts')
    #@RDATTR MESH:GET('nbcvs')
    #@GET    MESH:GET('pid')
    #@GET    MESH:GET('cvid')
    #@GET    MESH:GET('max pid')
    #@GET    MESH:GET('max cvid')
    #@GET    MESH:GET('pid from cvid')
    #@GET    MESH:GET('pid from coords')
    #@GET    MESH:GET('orphaned pid')
    #@GET    MESH:GET('cvid from pid')
    #@GET    MESH:GET('faces from pid')
    #@GET    MESH:GET('faces from cvid')
    #@GET    MESH:GET('outer faces')
    #@GET    MESH:GET('edges')
    #@GET    MESH:GET('curved edges')
    #@GET    MESH:GET('triangulated surface')
    #@GET    MESH:GET('normal of face')
    #@GET    MESH:GET('normal of faces')
    #@GET    MESH:GET('quality')
    #@GET    MESH:GET('cvstruct')
    #@GET    MESH:GET('geotrans')
    #@GET    MESH:GET('regions')
    #@GET    MESH:GET('region')
    #@GET    MESH:GET('save')
    #@GET    MESH:GET('char')
    #@GET    MESH:GET('export to vtk')
    #@GET    MESH:GET('export to dx')
    #@GET    MESH:GET('memsize')

    #@SET MESH:SET('pts')
    #@SET MESH:SET('add point')
    #@SET MESH:SET('del point')
    #@SET MESH:SET('add convex')
    #@SET MESH:SET('del convex')
    #@SET MESH:SET('del convex of dim')
    #@SET MESH:SET('translate')
    #@SET MESH:SET('transform')
    #@SET MESH:SET('merge')
    #@SET MESH:SET('optimize structure')
    #@SET MESH:SET('refine')
    #@SET MESH:SET('region')
    #@SET MESH:SET('region_intersect')
    #@SET MESH:SET('region_merge')
    #@SET MESH:SET('region_substract')
    #@SET MESH:SET('delete region')

class MeshFem:
    def __init__(self, *args):
        """General constructor for MeshFem objects.

* MeshFem(mesh M [, int Qdim=1])
Build a new MeshFem object. The Qdim parameter is optional.
        @INIT MESHFEM:INIT('load')
        @INIT MESHFEM:INIT('from string')
        @INIT MESHFEM:INIT('clone')
        """
        generic_constructor(self,'mesh_fem',*args)
    def __del__(self):
        generic_destructor(self,destructible=True)
    def get(self, *args):
        return getfem('mesh_fem_get',self.id, *args)
    def set(self, *args):
        return getfem('mesh_fem_set',self.id, *args)
    def __str__(self):
        return self.char()
    def __repr__(self):
        return '<getfem.MeshFem Q=%d, %d dof, %d (+%d) bytes>' % \
               (self.qdim(),self.nbdof(),self.memsize(), \
                self.linked_mesh().memsize())
    #@RDATTR MESHFEM:GET('nbdof')
    #@GET    MESHFEM:GET('dof from cv')
    #@GET    MESHFEM:GET('dof from cvid')
    #@GET    MESHFEM:GET('non conformal dof')
    #@RDATTR MESHFEM:GET('qdim')
    #@GET    MESHFEM:GET('fem')
    #@GET    MESHFEM:GET('is_lagrangian')
    #@GET    MESHFEM:GET('is_equivalent')
    #@GET    MESHFEM:GET('is_polynomial')
    #@GET    MESHFEM:GET('dof on region')
    #@GET    MESHFEM:GET('dof nodes')
    #@GET    MESHFEM:GET('dof partition')
    #@GET    MESHFEM:GET('interpolate_convex_data')
    #@GET    MESHFEM:GET('save')
    #@GET    MESHFEM:GET('char')
    #@GET    MESHFEM:GET('linked mesh')
    #@GET    MESHFEM:GET('export to vtk')
    #@GET    MESHFEM:GET('export to dx')
    #@GET    MESHFEM:GET('memsize')
    #@SET MESHFEM:SET('fem')
    #@SET MESHFEM:SET('classical fem')
    #@SET MESHFEM:SET('classical discontinuous fem')
    #@SET MESHFEM:SET('qdim')
    #@SET MESHFEM:SET('dof partition')
    def eval(self, expression):
        """interpolate an expression on the (lagrangian) MeshFem.

        Examples:
          mf.eval('x[0]*x[1]') interpolates the function 'x*y'
          mf.eval('[x[0],x[1]]') interpolates the vector field '[x,y]'
        """
        P=self.dof_nodes()
        nbd = P.shape[1];

        if not self.is_lagrangian:
            raise RuntimeError('cannot eval on a non-Lagragian MeshFem')
        if self.qdim() != 1:
            raise RuntimeError('only works (for now) with qdim == 1')
        x=P[:,0]; r=numarray.array(eval(expression))
        Z=numarray.zeros(r.shape + (nbd,),'d')
        for i in range(0,nbd):
            x=P[:,i]
            Z[...,i]=eval(expression)
        return Z


class MeshIm:
    def __init__(self, *args):
        """General constructor for MeshIm objects.

* MeshIm(mesh M, [{Integ Im|int IM_DEGREE}])
Build a new MeshIm object. For convenience, optional arguments (IM or IM_DEGREE) can be provided, in that case a call to MESHIM:SET('integ') is issued with these arguments.

        @INIT MESHIM:INIT('load')
        @INIT MESHIM:INIT('from string')
        @INIT MESHIM:INIT('clone')
        """
        generic_constructor(self,'mesh_im',*args)
    def __del__(self):
        generic_destructor(self,destructible=True)
    def get(self, *args):
        return getfem('mesh_im_get',self.id, *args)
    def set(self, *args):
        return getfem('mesh_im_set',self.id, *args)
    def __str__(self):
        return self.char()
    def __repr__(self):
        return '<getfem.MeshIm %d (+%d) bytes>' % \
               (self.memsize(), \
                self.linked_mesh().memsize())
    #@GET    MESHIM:GET('integ')
    #@GET    MESHIM:GET('eltm')
    #@GET    MESHIM:GET('save')
    #@GET    MESHIM:GET('char')
    #@GET    MESHIM:GET('linked mesh')
    #@GET    MESHIM:GET('memsize')
    #@SET MESHIM:SET('integ')


class MdBrick:
    def __init__(self, *args):
        """General constructor for MdBrick objects.

        @INIT MDBRICK:INIT ('constraint')
        @INIT MDBRICK:INIT ('dirichlet')
        @INIT MDBRICK:INIT ('dirichlet on normal component')
        @INIT MDBRICK:INIT ('dirichlet on normal derivative')
        @INIT MDBRICK:INIT ('generalized dirichlet')
        @INIT MDBRICK:INIT ('source term')
        @INIT MDBRICK:INIT ('normal source term')
        @INIT MDBRICK:INIT ('normal derivative source term')
        @INIT MDBRICK:INIT ('neumann KirchhoffLove source term')
        @INIT MDBRICK:INIT ('qu term')
        @INIT MDBRICK:INIT ('mass matrix')
        @INIT MDBRICK:INIT ('generic elliptic')
        @INIT MDBRICK:INIT ('helmholtz')
        @INIT MDBRICK:INIT ('isotropic linearized elasticity')
        @INIT MDBRICK:INIT ('linear incompressibility term')
        @INIT MDBRICK:INIT ('nonlinear elasticity')
        @INIT MDBRICK:INIT ('nonlinear elasticity incompressibility term')
        @INIT MDBRICK:INIT ('small deformations plasticity')
        @INIT MDBRICK:INIT ('dynamic')
        @INIT MDBRICK:INIT ('navier stokes')
        @INIT MDBRICK:INIT ('bilaplacian')
        @INIT MDBRICK:INIT ('isotropic_linearized_plate')
        @INIT MDBRICK:INIT ('mixed_isotropic_linearized_plate')
        @INIT MDBRICK:INIT ('plate_source_term')
        @INIT MDBRICK:INIT ('plate_simple_support')
        @INIT MDBRICK:INIT ('plate_clamped_support')
        @INIT MDBRICK:INIT ('plate_closing')
        """
        generic_constructor(self,'mdbrick',*args)
    def __del__(self):
        generic_destructor(self,destructible=True)
    def get(self, *args):
        return getfem('mdbrick_get',self.id, *args)
    def set(self, *args):
        return getfem('mdbrick_set',self.id, *args)
    def __str__(self):
        return self.char()
    def __repr__(self):
        return '<getfem.MdBrick %d bytes>' % \
               (self.memsize(),)
    #@RDATTR MDBRICK:GET('nbdof')
    #@RDATTR MDBRICK:GET('dim')
    #@RDATTR MDBRICK:GET('is_linear')
    #@RDATTR MDBRICK:GET('is_symmetric')
    #@RDATTR MDBRICK:GET('is_coercive')
    #@RDATTR MDBRICK:GET('is_complex')
    #@GET MDBRICK:GET('mixed_variables') 
    #@RDATTR MDBRICK:GET('subclass')
    #@GET MDBRICK:GET('param_list')
    #@GET MDBRICK:GET('param')
    #@GET MDBRICK:GET('solve')
    #@GET MDBRICK:GET('von mises')
    #@GET MDBRICK:GET('memsize')
    #@GET MDBRICK:GET('tresca')
    #@SET MDBRICK:SET('param')
    #@SET MDBRICK:SET('constraints');
    #@SET MDBRICK:SET('constraints_rhs');
    #@SET MDBRICK:SET('penalization_epsilon');

class MdState:
    def __init__(self, *args):
        """General constructor for MdState objects.

        These objects hold the global model data of a chain of
        MdBricks, such as the right hand side, the tangent matrix and
        the constraints.

        * MDS=MdState(mdbrick B)
        Build a modelstate for the brick B (selects the real or
        complex state from the complexity of B).

        @INIT MDSTATE:INIT ('real')

        @INIT MDSTATE:INIT ('complex')
        """
        generic_constructor(self,'mdstate',*args)
    def __del__(self):
        generic_destructor(self,destructible=True)
    def get(self, *args):
        return getfem('mdstate_get',self.id, *args)
    def set(self, *args):
        return getfem('mdstate_set',self.id, *args)
    def __str__(self):
        return self.char()
    def __repr__(self):
        return '<getfem.MdState %d bytes>' % \
               (self.memsize(),)
    #@RDATTR MDSTATE:GET('is_complex')
    #@GET MDSTATE:GET('tangent_matrix')
    #@GET MDSTATE:GET('constraints_matrix')
    #@GET MDSTATE:GET('reduced_tangent_matrix')
    #@GET MDSTATE:GET('constraints_nullspace')
    #@GET MDSTATE:GET('state')
    #@GET MDSTATE:GET('residual')
    #@GET MDSTATE:GET('reduced_residual')
    #@GET MDSTATE:GET('unreduce')
    #@GET MDSTATE:GET('memsize')
    #@SET MDSTATE:SET('compute_reduced_system')
    #@SET MDSTATE:SET('compute_reduced_residual')
    #@SET MDSTATE:SET('compute_residual')
    #@SET MDSTATE:SET('compute_tangent_matrix')
    #@SET MDSTATE:SET('clear')
    
class GeoTrans:
    def __init__(self, *args):
        """
        @TEXT GEOTRANS:INIT('GEOTRANS_init')
        """
        generic_constructor(self,'geotrans',*args)
    def __del__(self):
        generic_destructor(self,destructible=False)
    def get(self, *args):
        return getfem('geotrans_get',self.id, *args)
    def __str__(self):
        return self.get('char')
    def __repr__(self):
        return '<getfem.Geotrans '+str(self)+'>'
    #@RDATTR GEOTRANS:GET('dim')
    #@RDATTR GEOTRANS:GET('is_linear')
    #@RDATTR GEOTRANS:GET('nbpts')
    #@GET GEOTRANS:GET('pts')
    #@GET GEOTRANS:GET('normals')
    #@GET GEOTRANS:GET('transform')
    #@GET GEOTRANS:GET('char')

class Fem:
    """FEM (Finite Element Method) objects."""
    def __init__(self, fem_name):
        """Build a FEM object from a string description.

        @TEXT FEM:INIT('FEM_list')

        @INIT FEM:INIT('interpolated_fem')
        """
        generic_constructor(self,'fem',fem_name)
    def __del__(self):
        generic_destructor(self,destructible=False)
    def get(self, *args):
        return getfem('fem_get',self.id, *args)
    def __str__(self):
        return self.get('char')
    def __repr__(self):
        return '<getfem.Fem '+str(self)+'>'
    #@RDATTR FEM:GET('nbdof')
    #@RDATTR FEM:GET('dim')
    #@RDATTR FEM:GET('target_dim')
    #@GET FEM:GET('pts')
    #@RDATTR FEM:GET('is_equivalent')
    #@RDATTR FEM:GET('is_lagrange')
    #@RDATTR FEM:GET('is_polynomial')
    #@RDATTR FEM:GET('estimated_degree')
    #@GET FEM:GET('base_value')
    #@GET FEM:GET('grad_base_value')
    #@GET FEM:GET('hess_base_value')
    #@GET FEM:GET('poly_str')
    #@GET FEM:GET('char')
    
class Integ:
    """Integration Method Objects."""
    def __init__(self, *args):
        """
        @TEXT INTEG:INIT('INTEG_init')
        """
        generic_constructor(self,'integ',*args)
    def __del__(self):
        generic_destructor(self,destructible=False)
    def get(self, *args):
        return getfem('integ_get',self.id, *args)
    def __str__(self):
        return self.get('char')
    def __repr__(self):
        return '<getfem.Integ '+str(self)+'>'
    #@RDATTR INTEG:GET('is_exact')
    #@RDATTR INTEG:GET('dim')
    #@RDATTR INTEG:GET('nbpts')
    #@GET INTEG:GET('pts')
    #@GET INTEG:GET('coeffs')
    #@GET INTEG:GET('face_pts')
    #@GET INTEG:GET('face_coeffs')
    #@GET INTEG:GET('char')

class Eltm:
    """Descriptor for an elementary matrix type."""
    def __init__(self, *args):
        """
        @TEXT ELTM:INIT('ELTM_init')
        """
        generic_constructor(self,'eltm',*args)
    def __del__(self):
        generic_destructor(self,destructible=False)

    
class CvStruct:
    """Descriptor for a convex structure."""
    def __init__(self, *args):
        generic_constructor(self,'cvstruct',*args)
    def __del__(self):
        generic_destructor(self,destructible=False)
    def get(self, *args):
        return getfem('cvstruct_get',self.id, *args)
    def __repr__(self):
        return '<getfem.CvStruct %dD, %d pts>' % (self.dim(),self.nbpts())
    #@RDATTR CVSTRUCT:GET('nbpts')
    #@RDATTR CVSTRUCT:GET('dim')
    #@RDATTR CVSTRUCT:GET('basic structure')
    #@RDATTR CVSTRUCT:GET('face')
    #@GET CVSTRUCT:GET('facepts')

class Poly:
    pass

class Slice:
    """Mesh slices.

    The slices may be considered as a (non-conformal) mesh of
    simplexes which provides fast interpolation on a P1-discontinuous
    MeshFem.
    
    It is used mainly for post-processing purposes.
    """
    def __init__(self, *args):
        """General constructor for Slice objects.

@TEXT SLICE:INIT('constructor description')
"""
        generic_constructor(self,'slice',*args)
    def __del__(self):
        generic_destructor(self,destructible=True)
    def get(self, *args):
        return getfem('slice_get',self.id, *args)
    def set(self, *args):
        return getfem('slice_set',self.id, *args)
    def mesh(self):
        return self.get('linked_mesh')
    #@RDATTR SLICE:GET('dim')
    #@GET SLICE:GET('area')
    #@GET SLICE:GET('cvs')
    #@RDATTR SLICE:GET('nbpts')
    #@GET SLICE:GET('pts')
    #@RDATTR SLICE:GET('nbsplxs')
    #@GET SLICE:GET('splxs')
    #@GET SLICE:GET('edges')
    #@GET SLICE:GET('interpolate_convex_data')
    #@GET SLICE:GET('linked mesh')
    #@GET SLICE:GET('export to vtk')
    #@GET SLICE:GET('export to pov')
    #@GET SLICE:GET('export to dx')
    #@GET SLICE:GET('memsize')
    #@SET SLICE:SET('pts')

class Spmat:
    """Getfem sparse matrix."""
    def __init__(self, *args):
        """
        @TEXT SPMAT:INIT('SPMAT_init')
        """
        generic_constructor(self,'spmat',*args)
    def __del__(self):
        generic_destructor(self,destructible=True)
    def __str__(self):
        return self.get('info')
    def __repr__(self):
        return '<getfem.Spmat '+str(self)+'>'
    def __getitem__(self, key):
        return getfem('spmat_get',self.id, 'full',*key)
    def __setitem__(self, key, keyval):
        getfem('spmat_set', self.id, 'assign', key[0], key[1], keyval)
    def __neg__(self):
        m=Spmat('copy',self)
        m.scale(-1)
        return m
    def __add__(self, other):
        return Spmat('add',self,other)
    def __sub__(self, other):
        return Spmat('add',self,other.__neg__())
    def __mul__(self, other):
        """Multiplication of a Spmat with another Spmat or a vector or a scalar.

        The result is another Spmat object.
        """
        if (isinstance(other,(int,float,complex))):
            m=Spmat('copy',self)
            m.set('scale',other)
        elif (isinstance(other,list) or isinstance(other, numarray.NDArray)):
            m=self.mult(other)
        else:
            m=Spmat('mult',self,other)
        return m;
    def __rmul__(self, other):
        if (isinstance(other,(int,float,complex))):
            m=Spmat('copy',self)
            m.set('scale',other)
        elif (isinstance(other,list) or isinstance(other, numarray.NDArray)):
            m=self.tmult(other)
        else:
            m=Spmat('mult',other,self)
        return m;
    def get(self, *args):
        return getfem('spmat_get',self.id, *args)
    def set(self, *args):
        return getfem('spmat_set',self.id, *args)
    #@GET SPMAT:GET('size')
    #@GET SPMAT:GET('nnz')
    #@GET SPMAT:GET('is_complex')
    #@GET SPMAT:GET('storage')
    #@GET SPMAT:GET('full')
    #@GET SPMAT:GET('mult')
    #@GET SPMAT:GET('tmult')
    #@GET SPMAT:GET('diag')
    #@GET SPMAT:GET('csc_ind')
    #@GET SPMAT:GET('csc_val')
    #@GET SPMAT:GET('dirichlet nullspace')
    #@GET SPMAT:GET('info')
    #@GET SPMAT:GET('save')
    #@SET SPMAT:SET('clear')
    #@SET SPMAT:SET('scale')
    #@SET SPMAT:SET('transpose')
    #@SET SPMAT:SET('conjugate')
    #@SET SPMAT:SET('transconj')
    #@SET SPMAT:SET('to_csc')
    #@SET SPMAT:SET('to_wsc')
    #@SET SPMAT:SET('to_complex')
    #@SET SPMAT:SET('diag')
    #@SET SPMAT:SET('assign')
    #@SET SPMAT:SET('add')

    
class Precond:
    """Getfem preconditioner."""
    def __init__(self, *args):
        """
        @TEXT PRECOND:INIT('PRECOND_init')
        """
        generic_constructor(self,'precond',*args)
    def __del__(self):
        generic_destructor(self,destructible=True)
    def get(self, *args):
        return getfem('precond_get',self.id, *args)
    def set(self, *args):
        return getfem('precond_set',self.id, *args)
    #@GET PRECOND:GET('mult')
    #@GET PRECOND:GET('tmult')
    #@GET PRECOND:GET('type')
    #@GET PRECOND:GET('size')
    #@GET PRECOND:GET('is_complex')
    #@GET PRECOND:GET('info')



#@FUNC ::LINSOLVE('gmres')
#@FUNC ::LINSOLVE('cg')
#@FUNC ::LINSOLVE('bicgstab')
#@FUNC ::LINSOLVE('lu')
#@FUNC ::LINSOLVE('superlu')

#@FUNC ::COMPUTE('L2 norm')
#@FUNC ::COMPUTE('H1 semi norm')
#@FUNC ::COMPUTE('H1 norm')
#@FUNC ::COMPUTE('H2 semi norm')
#@FUNC ::COMPUTE('H2 norm')
#@FUNC ::COMPUTE('gradient')
#@FUNC ::COMPUTE('hessian')
#@FUNC ::COMPUTE('interpolate on')
#@FUNC ::COMPUTE('extrapolate on')
#@FUNC ::COMPUTE('error estimate')

#@FUNC ::ASM('volumic source')
#@FUNC ::ASM('boundary source')
#@FUNC ::ASM('mass matrix')
#@FUNC ::ASM('laplacian')
#@FUNC ::ASM('linear elasticity')
#@FUNC ::ASM('stokes')
#@FUNC ::ASM('helmholtz')
#@FUNC ::ASM('bilaplacian')
#@FUNC ::ASM('dirichlet')
#@FUNC ::ASM('boundary qu term')
#@FUNC ::ASM('volumic')
#@FUNC ::ASM('boundary')
#@FUNC ::ASM('interpolation matrix')
#@FUNC ::ASM('extrapolation matrix')

#@FUNC ::UTIL('trace level')
#@FUNC ::UTIL('warning level')

class LevelSet:
    """Getfem Level-set."""
    def __init__(self, *args):
        """
        @TEXT LEVELSET:INIT('LEVELSET_init')
        """
        generic_constructor(self,'levelset',*args)
    def __del__(self):
        generic_destructor(self,destructible=True)
    def get(self, *args):
        return getfem('levelset_get',self.id, *args)
    def set(self, *args):
        return getfem('levelset_set',self.id, *args)


def memstats():
    print "*** Getfem view of the workspace:"
    getfem('workspace','stats')
    print "*** Python view of the workspace:"
    for id,c in obj_count.iteritems():
        if (c):
            name=str(factory(id).__class__)
            print "%s class %d, id %d : instances=%d" % (name,id.classid,id.objid,c)


def linsolve(what, *args):
    return getfem('linsolve', what, *args)
def compute(mf, U, what, *args):
    return getfem('compute', mf, U, what, *args)
def asm(what, *args):
    return getfem('asm', what, *args)
def util(what, *args):
    return getfem('util', what, *args)


###
#register_types(Mesh,MeshFem,GeoTrans,Fem,Integ,Eltm,CvStruct,Poly,Slice)

def factory(id):
    # must be in the same order than enum getfemint_class_id in gfi_array.h
    t = (Mesh,MeshFem,MeshIm,MdBrick,MdState,GeoTrans,Fem,Integ,Eltm,CvStruct,Poly,Slice,Spmat,Precond,LevelSet)[id.classid]
    return t(id)

register_python_factory(factory)
