// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2008 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

#include <getfemint_misc.h>
#include <getfemint_pfem.h>
#include <getfemint_mesh_fem.h>
#include <getfemint_mesh_im.h>
#include <getfem/getfem_interpolated_fem.h>
#include <getfemint_workspace.h>

using namespace getfemint;
/*MLABCOM
  FUNCTION F=gf_fem(string FEM_NAME)
    Returns a handle F to one of the various Finite Elements Method defined
    in Getfem.

    example of FEM names are:
    @TEXT FEM:INIT('FEM_list')


    @INIT FEM:INIT('interpolated_fem')
MLABCOM*/

/*@TEXT FEM:INIT('FEM_list')
    * FEM_PK(N,K)<par>
      classical Lagrange element PK on a simplex of dimension N<par>
    * FEM_PK_DISCONTINUOUS(N,K[,alpha])} <par>
       discontinuous Lagrange element PK on a simplex of dim N<par>
    * FEM_QK(N,K)}      <par>
      classical Lagrange element QK on quadrangles, hexahedrons etc<par>
    * FEM_QK_DISCONTINUOUS(N,K[,alpha])}      <par>
      discontinuous Lagrange element QK on quadrangles, hexahedrons etc<par>
    * FEM_Q2_INCOMPLETE<par>
      incomplete 2D Q2 element with 8 dof (serendipity Quad 8 element) 
    * FEM_PK_PRISM(N,K) <par>
      classical Lagrange element PK on a prism<par>
    * FEM_PK_PRISM_DISCONTINUOUS(N,K[,alpha])<par>
      classical discontinuous Lagrange element PK on a prism.<par>

    * FEM_PK_WITH_CUBIC_BUBBLE(N,K)   <par>
      classical Lagrange element PK on a simplex with an additional volumic bubble function.<par>
    * FEM_P1_NONCONFORMING <par>
      non-conforming P1 method on a triangle.<par>
    * FEM_P1_BUBBLE_FACE(N)   <par>
      P1 method on a simplex with an additional bubble function on face 0.<par>
    * FEM_P1_BUBBLE_FACE_LAG  <par>
      P1 method on a simplex with an additional lagrange dof on face 0.<par>

    * FEM_PK_HIERARCHICAL(N,K)  <par>
      PK element with a hierarchical basis<par>
    * FEM_QK_HIERARCHICAL(N,K)  <par>
      QK element with a hierarchical basis<par>
    * FEM_PK_PRISM_HIERARCHICAL(N,K)   <par>
      PK element on a prism with a hierarchical basis<par>
    * FEM_STRUCTURED_COMPOSITE(FEM, K) <par>
      Composite fem on a grid with K divisions<par>
    * FEM_PK_HIERARCHICAL_COMPOSITE(N,K,S)<par>
      PK composite element on a grid with S subdivisions and with a hierarchical basis<par>
    * FEM_PK_FULL_HIERARCHICAL_COMPOSITE(N,K,S) <par>
      PK composite element with S subdivisions and a hierarchical basis on both degree and subdivision<par>
    * FEM_PRODUCT(FEM1,FEM2)<par>
      tensorial product of two polynomial elements<par>
    * FEM_HERMITE(N)<par>
      Hermite element P3 on a simplex of dimension N=1,2,3<par>
    * FEM_ARGYRIS<par>
      Argyris element P5 on the triangle.<par>
    * FEM_HCT_TRIANGLE<par> 
      Hsieh-Clough-Tocher element on the triangle (composite P3
      element which is C^1), should be used with IM_HCT_COMPOSITE()
      integration method.<par>
    * FEM_QUADC1_COMPOSITE<par>
      Quadrilateral element, composite P3 element and C^1 (16 dof).<par>
    * FEM_REDUCED_QUADC1_COMPOSITE<par>
      Quadrilateral element, composite P3 element and C^1 (12 dof).<par>
    * FEM(:_(BRT0(N)<par>
      Raviart-Thomas element of order 0 on a simplex of dimension N.<par>
    * FEM(:_(BNEDELEC(N)<par>
      Nedelec edge element of order 0 on a simplex of dimension N.<par>
      @*/

void gf_fem(getfemint::mexargs_in& in, getfemint::mexargs_out& out) {
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  std::string cmd = in.pop().to_string();
  getfemint_pfem *gfi_pf = 0;
  if (check_cmd(cmd, "interpolated fem", in, out, 2, 3, 0, 1)) {
    /*@INIT f=FEM:INIT('interpolated_fem', @mf MF1, @mim MIM2, [@ivec blocked_dof])
      Build a special fem which is interpolated from another mesh_fem.

      Using this special finite element, it is possible to interpolate
      a given mesh_fem MF1 on another mesh, given the integration
      method that will be used on this mesh.

      Note that this finite element may be quite slow, and eats much memory.
    @*/
    getfemint_mesh_fem *gfi_mf  = in.pop().to_getfemint_mesh_fem();
    getfemint_mesh_im  *gfi_mim = in.pop().to_getfemint_mesh_im();
    dal::bit_vector blocked_dof;
    if (in.remaining())
      blocked_dof = in.pop().to_bit_vector();
    getfem::pfem pf = 
      getfem::new_interpolated_fem(gfi_mf->mesh_fem(),
				   gfi_mim->mesh_im(), 
				   0, blocked_dof);
    gfi_pf = getfemint_pfem::get_from(pf);
    gfi_pf->nbdof_need_convex_number() = true;
    workspace().set_dependance(gfi_pf, gfi_mim);
    workspace().set_dependance(gfi_pf, gfi_mf);
  } else {
    getfem::pfem pf = getfem::fem_descriptor(cmd);
    gfi_pf = getfemint_pfem::get_from(pf, STATIC_OBJ | CONST_OBJ);
  }
  out.pop().from_object_id(gfi_pf->get_id(), FEM_CLASS_ID);
}
