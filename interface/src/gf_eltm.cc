// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2006-2006 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

#include <getfemint_misc.h>
#include <getfemint_matelemtype.h>
#include <getfemint_pfem.h>

using namespace getfemint;

/*MLABCOM

  FUNCTION ELTM = gf_eltm(elt_matrix_type, args..)
  @TEXT ELTM:INIT('ELTM_init')
MLABCOM*/

/*@TEXT ELTM:INIT('ELTM_init')
  Generates a descriptor for an elementary matrix type.<Par>

  * ELTM:INIT('base', @tfem FEM)<par>
  Integration of shape functions on elements, using the fem FEM.<Par>

  * ELTM:INIT('grad', @tfem FEM)<par>
  Integration of gradient of shape functions on elements, using the fem FEM.<Par>

  * ELTM:INIT('hessian', @tfem FEM)<par>
  Integration of hessian of shape functions on elements, using the fem FEM.<Par>

  * ELTM:INIT('normal')<par>
  The unit normal to the current convex face<Par>

  * ELTM:INIT('grad_geotrans')<par>
  The gradient of the geometric transformation.

  * ELTM:INIT('grad_geotrans_inv')<par>
  The inverse of the gradient of the geometric transformation.

  * ELTM:INIT('product', @eltm A, @eltm B)<par>
  Integration of the tensorial product of elementary matrices A and B.<Par>
  
  In order to obtain a numerical value of theses matrices, see MESHIM:GET('eltm').
@*/


void gf_eltm(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  std::string cmd                  = in.pop().to_string();
  getfem::pmat_elem_type pme = 0;
  if (check_cmd(cmd, "base", in, out, 1, 1, 0, 1)) {  
    pme = getfem::mat_elem_base(in.pop().to_fem());
  } else if (check_cmd(cmd, "grad", in, out, 1, 1, 0, 1)) {  
    pme = getfem::mat_elem_grad(in.pop().to_fem());
  } else if (check_cmd(cmd, "hessian", in, out, 1, 1, 0, 1)) {  
    pme = getfem::mat_elem_hessian(in.pop().to_fem());
  } else if (check_cmd(cmd, "normal", in, out, 0, 0, 0, 1)) {  
    pme = getfem::mat_elem_unit_normal();
  } else if (check_cmd(cmd, "grad_geotrans", in, out, 0, 0, 0, 1)) {
    pme = getfem::mat_elem_grad_geotrans(false); 
  } else if (check_cmd(cmd, "grad_geotrans_inv", in, out, 0, 0, 0, 1)) {
    pme = getfem::mat_elem_grad_geotrans(true); 
  } else if (check_cmd(cmd, "product", in, out, 2, 2, 0, 1)) {  
    getfem::pmat_elem_type  m1 = in.pop().to_mat_elem_type();
    getfem::pmat_elem_type  m2 = in.pop().to_mat_elem_type();
    pme = getfem::mat_elem_product(m1,m2);
  } else bad_cmd(cmd);
  out.pop().from_object_id(getfemint::ind_matelemtype(pme), ELTM_CLASS_ID);
}

