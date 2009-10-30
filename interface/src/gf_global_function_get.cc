// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2009 Luis Saavedra.
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
// $Id$
#include <getfemint.h>
#include <getfemint_global_function.h>

using namespace getfemint;

/*MLABCOM
  FUNCTION I = gf_global_function_get(GF, ...)
    General function for querying information about global_function objects.

  $Id$
MLABCOM*/

void gf_global_function_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) THROW_BADARG("Wrong number of input arguments");

  getfem::abstract_xy_function *paf = in.pop().to_global_function();

  std::string cmd = in.pop().to_string();
  if (check_cmd(cmd, "val", in, out, 0, 1, 0, 1)) {
    /*@GET VALs = GLOBALFUNCTION:GET('val',@mat PTs)
    Return `val` function evaluation in `PTs` (column points).@*/
    darray P = in.pop().to_darray(2,-1);// warning: 2 = paf->nd;
    darray V = out.pop().create_darray_h(P.getn());

    for (unsigned i=0; i < unsigned(P.getn()); i++) {
      V[i] = paf->val(P(0,i),P(1,i));
    }
  } else if (check_cmd(cmd, "grad", in, out, 0, 1, 0, 1)) {
    /*@GET GRADs = GLOBALFUNCTION:GET('grad',@mat PTs)
    Return `grad` function evaluation in `PTs` (column points).

    On return, each column of `GRADs` is of the
    form [Gx,Gy].@*/

    darray P = in.pop().to_darray(2,-1);// warning: 2 = paf->nd;
    darray G = out.pop().create_darray(2,P.getn());// warning: 2 = paf->nd;

    for (unsigned i=0; i < unsigned(P.getn()); i++) {
      getfem::base_small_vector g = paf->grad(P(0,i),P(1,i));
      G(0,i) = g[0];
      G(1,i) = g[1];
    }
  } else if (check_cmd(cmd, "hess", in, out, 0, 1, 0, 1)) {
    /*@GET HESSs = GLOBALFUNCTION:GET('hess',@mat PTs)
    Return `hess` function evaluation in `PTs` (column points).

    On return, each column of `HESSs` is of the
    form [Hxx,Hxy,Hyx,Hyy].@*/
    darray P = in.pop().to_darray(2,-1);// warning: 2 = paf->nd;
    darray H = out.pop().create_darray(4,P.getn());// warning: 4 = (paf->nd)*(paf->nd);

    for (unsigned i=0; i < unsigned(P.getn()); i++) {
      getfem::base_matrix h = paf->hess(P(0,i),P(1,i));
      H(0,i) = h(0,0);
      H(1,i) = h(0,1);
      H(2,i) = h(1,0);
      H(3,i) = h(1,1);
    }
  } else bad_cmd(cmd);
}
