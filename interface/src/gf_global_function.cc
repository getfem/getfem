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

#include <getfemint.h>
#include <getfemint_workspace.h>
#include <getfemint_global_function.h>
#include <getfem/getfem_mesh_fem_global_function.h>
#include <getfem/getfem_arch_config.h>

using namespace getfemint;

/*@TEXT GLOBALFUNCTION:INIT('GLOBALFUNCTION_init')
Global function object is represented by three functions:<Par>

* The global function `val`.<par>
* The global function gradient `grad`.<par>
* The global function Hessian `hess`.<Par>

this type of function is used as local and global enrichment<par>
function. The global function Hessian is an optional parameter<par>
(only for fourth order derivative problems).@*/

/*MLABCOM

  FUNCTION GF = gf_global_function(...)

  General constructor for global function object. Returns a getfem
  handle to the newly created global function object. Note that for
  recent (> 6.0) versions of matlab, you should replace the calls to
  'gf_global_function' with 'gfGlobalFunction' (this will instruct
  Matlab to consider the getfem global function as a regular matlab
  object that can be manipulated with get() and set() methods).

  @INIT GLOBALFUNCTION:INIT('cutoff')
  @INIT GLOBALFUNCTION:INIT('crack')
  @INIT GLOBALFUNCTION:INIT('parser')
  @INIT GLOBALFUNCTION:INIT('product')

  $Id: gf_global_function.cc 2963 2009-03-30 19:36:09Z lsaavedr $
MLABCOM*/

void gf_global_function(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) THROW_BADARG( "Wrong number of input arguments");
  getfemint_global_function *ggf = NULL;

  std::string cmd = in.pop().to_string();
  if (check_cmd(cmd, "cutoff", in, out, 4, 4, 0, 1)) {
    /*@INIT GF = GLOBALFUNCTION:INIT('cutoff', @int fn, @scalar r, @scalar r1, @scalar r0)
    Create a cutoff global function.@*/
    size_type   fn = in.pop().to_integer(-1,2);
    scalar_type r  = in.pop().to_scalar();
    scalar_type r1 = in.pop().to_scalar();
    scalar_type r0 = in.pop().to_scalar();

    getfem::abstract_xy_function *cutoff = new getfem::cutoff_xy_function(fn,r,r1,r0);
    ggf = getfemint_global_function::get_from(cutoff);
  } else if (check_cmd(cmd, "crack", in, out, 1, 1, 0, 1)) {
    /*@INIT GF = GLOBALFUNCTION:INIT('crack', @int fn)
    Create a near-tip asymptotic global function for modelling cracks.@*/
    size_type fn = in.pop().to_integer(0,11);

    getfem::abstract_xy_function *crack = new getfem::crack_singular_xy_function(fn);
    ggf = getfemint_global_function::get_from(crack);
#if GETFEM_HAVE_MUPARSER_MUPARSER_H
  } else if (check_cmd(cmd, "parser", in, out, 1, 3, 0, 1)) {
    /*@INIT GF = GLOBALFUNCTION:INIT('parser', @str val[, @str grad[, @str hess]])
    Create a global function from strings `val`, `grad` and `hess`.@*/
    std::string sval = in.pop().to_string();
    std::string sgrad = "0;0;";
    std::string shess = "0;0;0;0;";
    if (in.remaining() && in.front().is_string()) sgrad = in.pop().to_string();
    if (in.remaining() && in.front().is_string()) shess = in.pop().to_string();

    getfem::abstract_xy_function *parser = new getfem::parser_xy_function(sval,sgrad,shess);
    ggf = getfemint_global_function::get_from(parser);
#endif
  } else if (check_cmd(cmd, "product", in, out, 2, 2, 0, 1)) {
    /*@INIT GF = GLOBALFUNCTION:INIT('product', @tgf F, @tgf G)
    Create a product of two global functions.@*/
    getfem::abstract_xy_function *af1 = in.pop().to_global_function();
    getfem::abstract_xy_function *af2 = in.pop().to_global_function();

    getfem::abstract_xy_function *product = new getfem::product_of_xy_functions(*af1,*af2);
    ggf = getfemint_global_function::get_from(product);
  } else bad_cmd(cmd);
  out.pop().from_object_id(ggf->get_id(), GLOBAL_FUNCTION_CLASS_ID);
}
