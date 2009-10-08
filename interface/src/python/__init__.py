"""python-getfem Copyright (C) 2004-2009 Julien Pommier, insa, toulouse

The main package for getfem++ support for Python.  Normally used by importing
the package, and perhaps a particular module inside it.

This package is licensed under the GNU LGPL 2.1.

See the COPYING file in the getfem installation directory.
"""
from getfem import *

__revision__ = getfem_var('revision')
__version__ = getfem_var('version')
__author__ = getfem_var('author')
__url__ = getfem_var('url')
__license__ = getfem_var('license')
__docformat__ = 'restructuredtext'

gf_class = ['Mesh',
            'MeshFem',
            'MeshIm',
            'MdBrick',
            'MdState',
            'Model',
            'GeoTrans',
            'Fem',
            'Integ',
            'GlobalFunction',
            'Eltm',
            'CvStruct',
            'Poly',
            'Slice',
            'Spmat',
            'Precond',
            'LevelSet',
            'MeshLevelSet']

gf_linsolve = ['linsolve_gmres',
               'linsolve_cg',
               'linsolve_bicgstab',
               'linsolve_lu',
               'linsolve_superlu']

gf_compute = ['compute_L2_norm',
              'compute_H1_semi_norm',
              'compute_H1_norm',
              'compute_H2_semi_norm',
              'compute_H2_norm',
              'compute_gradient',
              'compute_hessian',
              'compute_interpolate_on',
              'compute_extrapolate_on',
              'compute_error_estimate']

gf_asm = ['asm_mass_matrix',
          'asm_laplacian',
          'asm_linear_elasticity',
          'asm_nonlinear_elasticity',
          'asm_stokes',
          'asm_helmholtz',
          'asm_bilaplacian',
          'asm_volumic_source',
          'asm_boundary_source',
          'asm_dirichlet',
          'asm_boundary_qu_term',
          'asm_volumic',
          'asm_boundary',
          'asm_interpolation_matrix',
          'asm_extrapolation_matrix']

gf_util = ['util_trace_level',
           'util_warning_level']

__all__ = ['getfem_var']+gf_class+gf_linsolve+gf_compute+gf_asm+gf_util
