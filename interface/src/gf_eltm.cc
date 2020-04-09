/*===========================================================================

 Copyright (C) 2006-2020 Yves Renard, Julien Pommier.

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/

#include <getfemint_misc.h>
#include <getfem/getfem_fem.h>
#include <getfem/getfem_mat_elem_type.h>

using namespace getfemint;

/*@GFDOC

  This object represents a type of elementary matrix. In order to obtain a
  numerical value of these matrices, see MESH_IM:GET('eltm').

  If you have very particular assembling needs, or if you just want to check
  the content of an elementary matrix, this function might be useful. But
  the generic assembly abilities of ::ASM(...) should suit most needs.
@*/


void gf_eltm(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  std::string cmd                  = in.pop().to_string();
  getfem::pmat_elem_type pme = 0;
  if (check_cmd(cmd, "base", in, out, 1, 1, 0, 1)) {
    /*@INIT E = ('base', @tfem FEM)
      return a descriptor for the integration of shape functions on
      elements, using the @tfem `FEM`. @*/
    pme = getfem::mat_elem_base(to_fem_object(in.pop()));
  } else if (check_cmd(cmd, "grad", in, out, 1, 1, 0, 1)) {
    /*@INIT E = ('grad', @tfem FEM)
      return a descriptor for the integration of the gradient of shape
      functions on elements, using the @tfem `FEM`.@*/
    pme = getfem::mat_elem_grad(to_fem_object(in.pop()));
  } else if (check_cmd(cmd, "hessian", in, out, 1, 1, 0, 1)) {
    /*@INIT E = ('hessian', @tfem FEM)
      return a descriptor for the integration of the hessian of shape
      functions on elements, using the @tfem `FEM`.@*/
    pme = getfem::mat_elem_hessian(to_fem_object(in.pop()));
  } else if (check_cmd(cmd, "normal", in, out, 0, 0, 0, 1)) {
    /*@INIT E = ('normal')
      return a descriptor for the unit normal of convex faces.@*/
    pme = getfem::mat_elem_unit_normal();
  } else if (check_cmd(cmd, "grad_geotrans", in, out, 0, 0, 0, 1)) {
    /*@INIT E = ('grad_geotrans')
      return a descriptor to the gradient matrix of the geometric
      transformation.@*/
    pme = getfem::mat_elem_grad_geotrans(false);
  } else if (check_cmd(cmd, "grad_geotrans_inv", in, out, 0, 0, 0, 1)) {
    /*@INIT E = ('grad_geotrans_inv')
      return a descriptor to the inverse of the gradient matrix of the
      geometric transformation (this is rarely used).@*/
    pme = getfem::mat_elem_grad_geotrans(true);
  } else if (check_cmd(cmd, "product", in, out, 2, 2, 0, 1)) {
    /*@INIT E = ('product', @teltm A, @teltm B)
      return a descriptor for the integration of the tensorial product of
      elementary matrices `A` and `B`.@*/
    getfem::pmat_elem_type  m1 = to_eltm_object(in.pop());
    getfem::pmat_elem_type  m2 = to_eltm_object(in.pop());
    pme = getfem::mat_elem_product(m1,m2);
  } else bad_cmd(cmd);
  out.pop().from_object_id(store_eltm_object(pme), ELTM_CLASS_ID);
}
