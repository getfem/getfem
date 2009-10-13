#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
#
#
# Python GetFEM++ interface
#
# Copyright (C) 2004-2009 Yves Renard, Julien Pommier.
#
# This file is a part of GetFEM++
#
# GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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
  Provides access to the pseudo-objects exported by the python-getfem interface.

  $Id$
"""

#__version__ = "$Revision$"
# $Source: getfem++/interface/src/python/getfem.base.py,v $

import sys
import numpy
import numbers

from numpy import *

from _getfem import *
obj_count = {}
getfem('workspace','clear all')

def generic_constructor(self,clname,*args):
    """Internal function -- acts as a constructor for all getfem objects."""
    #print 'generic_constructor.'+clname+'('+str(args)+')'
    if (len(args)==1 and type(args[0]) is GetfemObject):
      if hasattr(self,'id'):
        print "warning: hasattr(self,'id')!"
        print "self.id: ",self.id
        print "args[0]: ",args[0]
      else:
        self.id = args[0]
        if obj_count.get(self.id,0)==0:
          #print "Reviviendo objeto..."
          #print "self: ",self
          #print "self.id: ",self.id
          #if hasattr(self.id,'classid'):
          #  print "self.id.classid: ",self.id.classid
          #else:
          #  print "self.id.classid not found!"
          getfem("undelete",self.id)
          #print "self.id: ",self.id
          #pass
    else:
      self.id = getfem_from_constructor(clname,*args)
    obj_count[self.id] = obj_count.get(self.id,0)+1

def generic_destructor(self,destructible=True):
    """Internal function -- acts as a destructor for all getfem objects."""
    if (not hasattr(self,'id')):
      return
    #print "Mesh.__del__       ",self.id,'count=',obj_count[self.id]
    if (obj_count.has_key(self.id)):
      obj_count[self.id] = obj_count[self.id]-1
      if (destructible and obj_count[self.id] == 0):
        getfem('delete',self.id)
        #print "effective deletion"

# stub classes for getfem-interface objects

class Mesh:
    """Getfem Mesh Object.

Thos object is able to store any element in any dimension even
if you mix elements with different dimensions.
    """
    def __init__(self, *args):
      """General constructor for Mesh objects.

      @INIT MESH:INIT ('empty')
      @INIT MESH:INIT ('cartesian')
      @INIT MESH:INIT ('triangles grid')
      @INIT MESH:INIT ('regular simplices')
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
    #@RDATTR MESH:GET('nbpts')
    #@RDATTR MESH:GET('nbcvs')
    #@GET    MESH:GET('pts')
    #@GET    MESH:GET('pid')
    #@GET    MESH:GET('pid in faces')
    #@GET    MESH:GET('pid in cvids')
    #@GET    MESH:GET('pid in regions')
    #@GET    MESH:GET('pid from coords')
    #@GET    MESH:GET('pid from cvid')
    #@GET    MESH:GET('pts from cvid')
    #@GET    MESH:GET('cvid')
    #@GET    MESH:GET('max pid')
    #@GET    MESH:GET('max cvid')
    #@GET    MESH:GET('edges')
    #@GET    MESH:GET('curved edges')
    #@GET    MESH:GET('orphaned pid')
    #@GET    MESH:GET('cvid from pid')
    #@GET    MESH:GET('faces from pid')
    #@GET    MESH:GET('outer faces')
    #@GET    MESH:GET('faces from cvid')
    #@GET    MESH:GET('triangulated surface')
    #@GET    MESH:GET('normal of face')
    #@GET    MESH:GET('normal of faces')
    #@GET    MESH:GET('quality')
    #@GET    MESH:GET('convex area')
    #@GET    MESH:GET('cvstruct')
    #@GET    MESH:GET('geotrans')
    #@GET    MESH:GET('regions')
    #@GET    MESH:GET('region')
    #@GET    MESH:GET('save')
    #@GET    MESH:GET('char')
    #@GET    MESH:GET('export to vtk')
    #@GET    MESH:GET('export to dx')
    #@GET    MESH:GET('export to pos')
    #@GET    MESH:GET('memsize')

    #@SET    MESH:SET('pts')
    #@SET    MESH:SET('add point')
    #@SET    MESH:SET('del point')
    #@SET    MESH:SET('add convex')
    #@SET    MESH:SET('del convex')
    #@SET    MESH:SET('del convex of dim')
    #@SET    MESH:SET('translate')
    #@SET    MESH:SET('transform')
    #@SET    MESH:SET('region')
    #@SET    MESH:SET('region intersect')
    #@SET    MESH:SET('region merge')
    #@SET    MESH:SET('region substract')
    #@SET    MESH:SET('delete region')
    #@SET    MESH:SET('merge')
    #@SET    MESH:SET('optimize structure')
    #@SET    MESH:SET('refine')


class MeshFem:
    def __init__(self, *args):
      """General constructor for MeshFem objects.

      @INIT MESHFEM:INIT('load')
      @INIT MESHFEM:INIT('from string')
      @INIT MESHFEM:INIT('clone')
      @INIT MESHFEM:INIT('sum')
      @INIT MESHFEM:INIT('levelset')
      @INIT MESHFEM:INIT('global function')
      @INIT MESHFEM:INIT('partial')
      @INIT MESHFEM:INIT('.mesh')
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
    #@RDATTR MESHFEM:GET('nb basic dof')
    #@GET    MESHFEM:GET('basic dof from cv')
    #@GET    MESHFEM:GET('dof from cv')
    #@GET    MESHFEM:GET('basic dof from cvid')
    #@GET    MESHFEM:GET('dof from cvid')
    #@GET    MESHFEM:GET('non conformal basic dof')
    #@GET    MESHFEM:GET('non conformal dof')
    #@RDATTR MESHFEM:GET('qdim')
    #@GET    MESHFEM:GET('fem')
    #@GET    MESHFEM:GET('convex_index')
    #@GET    MESHFEM:GET('is_lagrangian')
    #@GET    MESHFEM:GET('is_equivalent')
    #@GET    MESHFEM:GET('is_polynomial')
    #@RDATTR MESHFEM:GET('is_reduced')
    #@GET    MESHFEM:GET('reduction matrix')
    #@GET    MESHFEM:GET('extension matrix')
    #@GET    MESHFEM:GET('dof on region')
    #@GET    MESHFEM:GET('basic dof on region')
    #@GET    MESHFEM:GET('basic dof nodes')
    #@GET    MESHFEM:GET('dof nodes')
    #@GET    MESHFEM:GET('dof partition')
    #@GET    MESHFEM:GET('save')
    #@GET    MESHFEM:GET('char')
    #@GET    MESHFEM:GET('linked mesh')
    #@GET    MESHFEM:GET('export to vtk')
    #@GET    MESHFEM:GET('export to dx')
    #@GET    MESHFEM:GET('export to pos')
    #@GET    MESHFEM:GET('dof_from_im')
    #@GET    MESHFEM:GET('interpolate_convex_data')
    #@GET    MESHFEM:GET('memsize')
    #@GET    MESHFEM:GET('has_linked_mesh_levelset')
    #@GET    MESHFEM:GET('linked_mesh_levelset')

    #@SET    MESHFEM:SET('fem')
    #@SET    MESHFEM:SET('classical fem')
    #@SET    MESHFEM:SET('classical discontinuous fem')
    #@SET    MESHFEM:SET('qdim')
    #@SET    MESHFEM:SET('reduction')
    #@SET    MESHFEM:SET('reduction matrices')
    #@SET    MESHFEM:SET('dof partition')
    def eval(self, expression, gl={}, lo={}):
      """interpolate an expression on the (lagrangian) MeshFem.

Examples:

>>> mf.eval('x[0]*x[1]') # interpolates the function 'x*y'
>>> mf.eval('[x[0],x[1]]') # interpolates the vector field '[x,y]'

>>> import numpy as np
>>> mf.eval('np.sin(x[0])',globals(),locals()) # interpolates the function sin(x)
      """
      P = self.basic_dof_nodes()
      nbd = P.shape[1]

      if not self.is_lagrangian:
        raise RuntimeError('cannot eval on a non-Lagragian MeshFem')
      if self.qdim() != 1:
        Ind = numpy.arange(0,nbd,self.qdim()) # = sdof
        P   = P[:,Ind]
        nbd = P.shape[1] # = nb_sdof
      x = P[:,0]
      gl['x'] = P[:,0]
      lo['x'] = P[:,0]
      r = numpy.array(eval(expression,gl,lo))
      Z = numpy.zeros(r.shape + (nbd,), r.dtype)
      for i in xrange(0,nbd):
        gl['x'] = P[:,i]
        lo['x'] = P[:,i]
        Z[...,i] = eval(expression,gl,lo)
      return Z


class MeshIm:
    def __init__(self, *args):
      """General constructor for MeshIm objects.

      @INIT MESHIM:INIT('load')
      @INIT MESHIM:INIT('from string')
      @INIT MESHIM:INIT('clone')
      @INIT MESHIM:INIT('levelset')
      @INIT MESHIM:INIT('.mesh')
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
    #@GET MESHIM:GET('integ')
    #@GET MESHFEM:GET('convex_index')
    #@GET MESHIM:GET('eltm')
    #@GET MESHIM:GET('im_nodes')
    #@GET MESHIM:GET('save')
    #@GET MESHIM:GET('char')
    #@GET MESHIM:GET('linked mesh')
    #@GET MESHIM:GET('memsize')

    #@SET MESHIM:SET('integ')


class MdBrick:
    """Getfem MdBrick Object.

A model brick is basically an object which modifies a global tangent
matrix and its associated right hand side. Typical modifications are
insertion of the stiffness matrix for the problem considered (linear
elasticity, laplacian, ...), handling of a set of contraints, Dirichlet
condition, addition of a source term to the right hand side, etc. The
global tangent matrix and its right hand side are stored in a MdState
object.
    """
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
      @INIT MDBRICK:INIT ('bilaplacian')
      @INIT MDBRICK:INIT ('navier stokes')
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
    #@RDATTR MDBRICK:GET('nb_constraints')
    #@RDATTR MDBRICK:GET('is_linear')
    #@RDATTR MDBRICK:GET('is_symmetric')
    #@RDATTR MDBRICK:GET('is_coercive')
    #@RDATTR MDBRICK:GET('is_complex')
    #@GET    MDBRICK:GET('mixed_variables')
    #@RDATTR MDBRICK:GET('subclass')
    #@GET    MDBRICK:GET('param_list')
    #@GET    MDBRICK:GET('param')
    #@GET    MDBRICK:GET('solve')
    #@GET    MDBRICK:GET('von mises')
    #@GET    MDBRICK:GET('tresca')
    #@GET    MDBRICK:GET('memsize')

    #@SET    MDBRICK:SET('param')
    #@SET    MDBRICK:SET('penalization_epsilon');
    #@SET    MDBRICK:SET('constraints');
    #@SET    MDBRICK:SET('constraints_rhs');

class MdState:
    """Getfem MdState Object.

A model state is an object which store the state data for a chain of
model bricks. This includes the global tangent matrix, the right hand
side and the constraints.
    """
    def __init__(self, *args):
      """General constructor for MdState objects.
There are two sorts of model states, the 'real' and the 'complex'
model states.

      @INIT MDSTATE:INIT('.mdbrick')
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
    #@GET    MDSTATE:GET('tangent_matrix')
    #@GET    MDSTATE:GET('constraints_matrix')
    #@GET    MDSTATE:GET('reduced_tangent_matrix')
    #@GET    MDSTATE:GET('constraints_nullspace')
    #@GET    MDSTATE:GET('state')
    #@GET    MDSTATE:GET('residual')
    #@GET    MDSTATE:GET('reduced_residual')
    #@GET    MDSTATE:GET('unreduce')
    #@GET    MDSTATE:GET('memsize')

    #@SET    MDSTATE:SET('compute_reduced_system')
    #@SET    MDSTATE:SET('compute_reduced_residual')
    #@SET    MDSTATE:SET('compute_residual')
    #@SET    MDSTATE:SET('compute_tangent_matrix')
    #@SET    MDSTATE:SET('state')
    #@SET    MDSTATE:SET('clear')


class Model:
    """Getfem Model Object.

A model is an object which store all the state variable and the data of a
model and a list of bricks. A brick is a component of the model, i.e. a
term wich link some state variables. The model object  includes the global
tangent matrix, the right hand side and the constraints.
    """
    def __init__(self, *args):
      """General constructor for Model objects.
There are two sorts of model states, the 'real' and the 'complex'
model states.

      @INIT MODEL:INIT ('real')
      @INIT MODEL:INIT ('complex')
      """
      generic_constructor(self,'model',*args)
    def __del__(self):
      generic_destructor(self,destructible=True)
    def get(self, *args):
      return getfem('model_get',self.id, *args)
    def set(self, *args):
      return getfem('model_set',self.id, *args)
    def __str__(self):
      return self.char()
    def __repr__(self):
      return '<getfem.Model %d bytes>' % \
             (self.memsize(),)
    #@GET    MODEL:GET('is_complex')
    #@GET    MODEL:GET('tangent_matrix')
    #@GET    MODEL:GET('rhs')
    #@GET    MODEL:GET('memsize')
    #@GET    MODEL:GET('listvar')
    #@GET    MODEL:GET('listbricks')
    #@GET    MODEL:GET('variable')
    #@GET    MODEL:GET('mult varname Dirichlet')
    #@GET    MODEL:GET('from variables')
    #@GET    MODEL:GET('assembly')
    #@GET    MODEL:GET('solve')
    #@GET    MODEL:GET('compute isotropic linearized Von Mises or Tresca')

    #@SET    MODEL:SET('clear')
    #@SET    MODEL:SET('add fem variable')
    #@SET    MODEL:SET('add variable')
    #@SET    MODEL:SET('add multiplier')
    #@SET    MODEL:SET('add fem data')
    #@SET    MODEL:SET('add initialized fem data')
    #@SET    MODEL:SET('add data')
    #@SET    MODEL:SET('add initialized data')
    #@SET    MODEL:SET('variable')
    #@SET    MODEL:SET('to variables')
    #@SET    MODEL:SET('add Laplacian brick')
    #@SET    MODEL:SET('add generic elliptic brick')
    #@SET    MODEL:SET('add source term brick')
    #@SET    MODEL:SET('add normal source term brick')
    #@SET    MODEL:SET('add Dirichlet condition with multipliers')
    #@SET    MODEL:SET('add Dirichlet condition with penalization')
    #@SET    MODEL:SET('change penalization coeff')
    #@SET    MODEL:SET('add Helmholtz brick')
    #@SET    MODEL:SET('add Fourier Robin brick')
    #@SET    MODEL:SET('add constraint with multipliers')
    #@SET    MODEL:SET('add constraint with penalization')
    #@SET    MODEL:SET('add explicit matrix')
    #@SET    MODEL:SET('add explicit rhs')
    #@SET    MODEL:SET('set private matrix')
    #@SET    MODEL:SET('set private rhs')
    #@SET    MODEL:SET('add isotropic linearized elasticity brick')
    #@SET    MODEL:SET('add linear incompressibility brick')
    #@SET    MODEL:SET('add mass brick')
    #@SET    MODEL:SET('add basic d on dt brick')
    #@SET    MODEL:SET('add basic d2 on dt2 brick')
    #@SET    MODEL:SET('add theta method dispatcher')
    #@SET    MODEL:SET('add midpoint dispatcher')
    #@SET    MODEL:SET('velocity update for order two theta method')
    #@SET    MODEL:SET('velocity update for Newmark scheme')
    #@SET    MODEL:SET('disable bricks')
    #@SET    MODEL:SET('unable bricks')
    #@SET    MODEL:SET('first iter')
    #@SET    MODEL:SET('next iter')


class GeoTrans:
    """General function for building descriptors to geometric transformations.

The geometric transformation must be used when you are building
a custom mesh convex by convex (see the add_convex() function of
getfem.Mesh): it also defines the kind of convex (triangle,
hexahedron, prism, etc..)
    """
    def __init__(self, *args):
      """Build a GeoTrans object from a string description.

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
    #@GET    GEOTRANS:GET('pts')
    #@GET    GEOTRANS:GET('normals')
    #@GET    GEOTRANS:GET('transform')
    #@GET    GEOTRANS:GET('char')


class Fem:
    """FEM (Finite Element Method) objects."""
    def __init__(self, fem_name):
      """Build a FEM object from a string description.

      @TEXT FEM:INIT('FEM_list')

**SPECIAL FEM:**

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
    #@GET    FEM:GET('pts')
    #@RDATTR FEM:GET('is_equivalent')
    #@RDATTR FEM:GET('is_lagrange')
    #@RDATTR FEM:GET('is_polynomial')
    #@RDATTR FEM:GET('estimated_degree')
    #@GET    FEM:GET('base_value')
    #@GET    FEM:GET('grad_base_value')
    #@GET    FEM:GET('hess_base_value')
    #@GET    FEM:GET('poly_str')
    #@GET    FEM:GET('char')


class Integ:
    """Integration Method Objects.

General object for obtaining handles to various integrations
methods on convexes (used when the elementary matrices are built).
    """
    def __init__(self, *args):
      """Return a FEM Integration Method from a string description.

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
    #@GET    INTEG:GET('pts')
    #@GET    INTEG:GET('face_pts')
    #@GET    INTEG:GET('coeffs')
    #@GET    INTEG:GET('face_coeffs')
    #@GET    INTEG:GET('char')


class GlobalFunction:
    """Getfem Global Function Object.

    @TEXT GLOBALFUNCTION:INIT('GLOBALFUNCTION_init')
    """
    def __init__(self, *args):
      """General constructor for GlobalFunction object:

      @INIT GLOBALFUNCTION:INIT('cutoff')
      @INIT GLOBALFUNCTION:INIT('crack')
      @INIT GLOBALFUNCTION:INIT('parser')
      """
      if isinstance(args[0],str):
        if args[0]=='parser' and not(getfem_var('muParser')=="1"):
          raise RuntimeError("Option \'parser\' need the package muParser.")
      generic_constructor(self,'global_function',*args)
    def __del__(self):
      generic_destructor(self,destructible=False)
    def get(self, *args):
      return getfem('global_function_get',self.id, *args)
    def set(self, *args):
      return getfem('global_function_set',self.id, *args)
    def __mul__(self,other):
      if isinstance(other,numbers.Number):
        return GlobalFunction('product',self,GlobalFunction('parser',"%e"%(other)))
      return GlobalFunction('product',self,other)
    def __add__(self,other):
      if isinstance(other,numbers.Number):
        return GlobalFunction('add',self,GlobalFunction('parser',"%e"%(other)))
      return GlobalFunction('add',self,other)
    def __call__(self,Pts):
      return getfem('global_function_get',self.id, 'val', Pts)
    #@GET    GLOBALFUNCTION:GET('val')
    #@GET    GLOBALFUNCTION:GET('grad')
    #@GET    GLOBALFUNCTION:GET('hess')


class Eltm:
    """Descriptor for an elementary matrix type.

If you have very particular assembling needs, or if you just want to
check the content of an elementary matrix, this function might be
useful. But the generic assembly abilities of getfem.asm_* should
suit most needs.
    """
    def __init__(self, *args):
      """Generates a descriptor for an elementary matrix type.

      @TEXT ELTM:INIT('ELTM_init')
      """
      generic_constructor(self,'eltm',*args)
    def __del__(self):
      generic_destructor(self,destructible=False)


class CvStruct:
    """Descriptor for a convex structure.

The convex structures are internal structures of getfem++. They do
not contain points positions. These structures are recursive, since
the faces of a convex structures are convex structures.
    """
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
    #@GET    CVSTRUCT:GET('facepts')

class Poly:
    pass

class Slice:
    """Mesh slices.

The slices may be considered as a (non-conformal) mesh of simplexes
which provides fast interpolation on a P1-discontinuous MeshFem.

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
    #@GET    SLICE:GET('area')
    #@GET    SLICE:GET('cvs')
    #@RDATTR SLICE:GET('nbpts')
    #@RDATTR SLICE:GET('nbsplxs')
    #@GET    SLICE:GET('pts')
    #@GET    SLICE:GET('splxs')
    #@GET    SLICE:GET('edges')
    #@GET    SLICE:GET('interpolate_convex_data')
    #@GET    SLICE:GET('linked mesh')
    #@GET    SLICE:GET('memsize')
    #@GET    SLICE:GET('export to vtk')
    #@GET    SLICE:GET('export to pov')
    #@GET    SLICE:GET('export to dx')
    #@GET    SLICE:GET('export to pos')

    #@SET    SLICE:SET('pts')

class Spmat:
    """Getfem sparse matrix."""
    def __init__(self, *args):
      """General constructor for getfem sparse matrices.

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
      if isinstance(other,numbers.Number):
        m = Spmat('copy',self)
        m.set('scale',other)
      elif (isinstance(other,list) or isinstance(other, numpy.ndarray)):
        m = self.mult(other)
      else:
        m = Spmat('mult',self,other)
      return m
    def __rmul__(self, other):
      if isinstance(other,numbers.Number):
        m=Spmat('copy',self)
        m.set('scale',other)
      elif (isinstance(other,list) or isinstance(other, numpy.ndarray)):
        m=self.tmult(other)
      else:
        m=Spmat('mult',other,self)
      return m;
    def get(self, *args):
      return getfem('spmat_get',self.id, *args)
    def set(self, *args):
      return getfem('spmat_set',self.id, *args)
    #@GET SPMAT:GET('nnz')
    #@GET SPMAT:GET('full')
    #@GET SPMAT:GET('mult')
    #@GET SPMAT:GET('tmult')
    #@GET SPMAT:GET('diag')
    #@GET SPMAT:GET('storage')
    #@GET SPMAT:GET('size')
    #@GET SPMAT:GET('is_complex')
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
      """General constructor for getfem preconditioners.

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


class LevelSet:
    """Getfem Level-Set Object.

    @TEXT LEVELSET:INIT('LEVELSET_init')
    """
    def __init__(self, *args):
      """General constructor for LevelSet objects.

      @INIT LEVELSET:INIT('.mesh')
      """
      generic_constructor(self,'levelset',*args)
    def __del__(self):
      generic_destructor(self,destructible=True)
    def get(self, *args):
      return getfem('levelset_get',self.id, *args)
    def set(self, *args):
      return getfem('levelset_set',self.id, *args)

    #@GET    LEVELSET:GET('values')
    #@RDATTR LEVELSET:GET('degree')
    #@GET    LEVELSET:GET('mf')
    #@RDATTR LEVELSET:GET('memsize')

    #@SET    LEVELSET:SET('values')
    #@SET    LEVELSET:SET('simplify')

class MeshLevelSet:
    """Getfem Mesh-Level-Set Object.
    """
    def __init__(self, *args):
      """General constructor for MeshLevelSet objects.
      
      @INIT MESHLEVELSET:INIT('.mesh')
      """
      generic_constructor(self,'mesh_levelset',*args)
    def __del__(self):
      generic_destructor(self,destructible=True)
    def get(self, *args):
      return getfem('mesh_levelset_get',self.id, *args)
    def set(self, *args):
      return getfem('mesh_levelset_set',self.id, *args)

    #@GET    MESHLEVELSET:GET('cut_mesh')
    #@GET    MESHLEVELSET:GET('linked_mesh')
    #@GET    MESHLEVELSET:GET('nb_ls')
    #@GET    MESHLEVELSET:GET('levelsets')
    #@GET    MESHLEVELSET:GET('crack_tip_convexes')
    #@GET    MESHLEVELSET:GET('memsize')

    #@SET    MESHLEVELSET:SET('add')
    #@SET    MESHLEVELSET:SET('sup')
    #@SET    MESHLEVELSET:SET('adapt')


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

#@FUNC ::ASM('mass matrix')
#@FUNC ::ASM('laplacian')
#@FUNC ::ASM('linear elasticity')
#@FUNC ::ASM('nonlinear elasticity')
#@FUNC ::ASM('stokes')
#@FUNC ::ASM('helmholtz')
#@FUNC ::ASM('bilaplacian')
#@FUNC ::ASM('volumic source')
#@FUNC ::ASM('boundary source')
#@FUNC ::ASM('dirichlet')
#@FUNC ::ASM('boundary qu term')
#@FUNC ::ASM('volumic')
#@FUNC ::ASM('boundary')
#@FUNC ::ASM('interpolation matrix')
#@FUNC ::ASM('extrapolation matrix')

#@FUNC ::UTIL('trace level')
#@FUNC ::UTIL('warning level')

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

def factory(id):
    # must be in the same order than enum getfemint_class_id in gfi_array.h
    t = ( Mesh,
          MeshFem,
          MeshIm,
          MdBrick,
          MdState,
          Model,
          GeoTrans,
          Fem,
          Integ,
          Eltm,
          CvStruct,
          Poly,
          Slice,
          Spmat,
          Precond,
          LevelSet,
          MeshLevelSet,
          GlobalFunction)[id.classid]
    return t(id)

register_python_factory(factory)
