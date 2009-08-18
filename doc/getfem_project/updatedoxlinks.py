#!/usr/bin/python
# Copyright (C) 2001-2009 Yves Renard
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

import re
import glob

def doxrename(f):
    latexmacro = re.sub('[._:]','',f)
    latexmacro = re.sub('0','zero',latexmacro)
    latexmacro = re.sub('1','one',latexmacro)
    latexmacro = re.sub('2','two',latexmacro)

    doxname = re.sub('_','__',f)
    doxname = re.sub('\.','_8',doxname)
    doxname = re.sub(':','_1',doxname)

    escapedname = re.sub('_', '\\_', f)
    
    return (latexmacro,doxname,escapedname)

flist=glob.glob1('../../src/', '*.h') + glob.glob1('../../src/', '*.cc')

out=file('doxygenlinks.tex','wt')

for f in flist:
    n = doxrename(f)
    print "doing file %s" % (n,)
    out.write('\\newcommand{\\%s}{\\doxfilename{%s}{%s}}\\xspace\n' % (n[0],n[2], n[1]))

classes="""
dal::bit_vector
dal::bv_visitor
bgeot::convex_structure/bgeot::pconvex_structure
bgeot::convex_ref/bgeot::pconvex_ref
bgeot::geometric_trans/bgeot::pgeometric_trans
getfem::virtual_fem/getfem::pfem
getfem::mesh
getfem::mesh_region
getfem::mr_visitor
bgeot::mesh_structure
getfem::generic_assembly
getfem::mesh_im
getfem::mesh_fem
getfem::stored_mesh_slice
getfem::slicer_action
getfem::mesh_slice_cv_dof_data_base
getfem::slicer_none
getfem::slicer_boundary
getfem::slicer_apply_deformation
getfem::slicer_half_space
getfem::slicer_sphere
getfem::slicer_cylinder
getfem::slicer_isovalues
getfem::slicer_mesh_with_mesh
getfem::slicer_union
getfem::slicer_intersect
getfem::slicer_complementary
getfem::slicer_build_mesh
getfem::slicer_build_edges_mesh
getfem::slicer_build_stored_mesh_slice
getfem::slicer_explode
getfem::mesh_slicer
getfem::dx_export
getfem::vtk_export
getfem::level_set
getfem::mesh_level_set
getfem::mesh_im_level_set
getfem::mesh_fem_level_set
getfem::model_state
getfem::mdbrick_abstract_common_base
getfem::mdbrick_abstract
getfem::mdbrick_parameter
getfem::mdbrick_abstract_linear_pde
getfem::mdbrick_generic_elliptic
getfem::mdbrick_source_term
getfem::mdbrick_constraint
getfem::mdbrick_Dirichlet
getfem::mdbrick_isotropic_linearized_elasticity
getfem::mdbrick_QU_term
getfem::mdbrick_linear_incomp
getfem::mdbrick_plasticity
getfem::mdbrick_isotropic_linearized_plate
getfem::mdbrick_mixed_isotropic_linearized_plate
getfem::mdbrick_plate_source_term
getfem::mdbrick_plate_simple_support
getfem::mdbrick_plate_clamped_support
getfem::mdbrick_plate_closing
getfem::mdbrick_nonlinear_elasticity
getfem::mdbrick_nonlinear_incomp
struct getfem::abstract_hyperelastic_law
struct getfem::SaintVenant_Kirchhoff_hyperelastic_law
struct getfem::Ciarlet_Geymonat_hyperelastic_law
struct getfem::Mooney_Rivlin_hyperelastic_law 
gmm::iteration
"""

for c in classes.split('\n'):
    if (len(c) == 0):
        continue
    ftype = "class";
    ll=c.split(' ')
    if (len(ll)>1):
        ftype = ll[0]
        ll=ll[1]
    else:
        ll=ll[0]

    ll=ll.split('/')
    n = doxrename(ll[0])

    print "doing class %s" % (n,)


    out.write('\\newcommand{\\%s}{\\doxref{%s}{%s%s}}\\xspace\n' % (n[0],n[2], ftype, n[1]))

    if (len(ll)>1):
        m = doxrename(ll[1])
        print "doing alias %s" % (ll[1],)
        out.write('\\newcommand{\\%s}{\\doxref{%s}{%s%s}}\\xspace\n' % (m[0],m[2], ftype, n[1]))
