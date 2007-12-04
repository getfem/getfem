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

#include <getfemint.h>
#include <getfemint_pfem.h>

using namespace getfemint;

static size_type get_optional_convex_number(getfemint::mexargs_in &in, 
					    getfemint_pfem *gfi_fem, 
					    const std::string cmd) {
  size_type cv = 0; 
  if (!in.remaining() && gfi_fem->nbdof_need_convex_number())
    THROW_BADARG("This FEM requires a convex number for " << cmd);
  if (in.remaining()) 
    cv = in.pop().to_integer() - config::base_index();
  return cv;
}

/*MLABCOM
  FUNCTION I = gf_fem_get(F, ...)
    General function for querying information about FEM objects.
    
    @RDATTR FEM:GET('nbdof')
    @RDATTR FEM:GET('dim')
    @RDATTR FEM:GET('target_dim')
    @GET FEM:GET('pts')
    @RDATTR FEM:GET('is_equivalent')
    @RDATTR FEM:GET('is_lagrange')
    @RDATTR FEM:GET('is_polynomial')
    @RDATTR FEM:GET('estimated_degree')
    @GET FEM:GET('base_value')
    @GET FEM:GET('grad_base_value')
    @GET FEM:GET('hess_base_value')
    @GET FEM:GET('poly_str')
    @GET FEM:GET('char')
MLABCOM*/

void gf_fem_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfemint_pfem *gfi_fem = in.pop().to_getfemint_pfem();
  getfem::pfem fem = gfi_fem->pfem();//in.pop().to_fem();
  std::string cmd = in.pop().to_string();
  if (check_cmd(cmd, "nbdof", in, out, 0, 1, 0, 1)) {
    /*@RDATTR FEM:GET('nbdof' [, @int CV])
      Return the number of DOF for the FEM.

      Some specific FEM (for example 'interpolated_fem') may require a
      convex number CV to give their result (interpolated fems for
      example). In most of the case, you can omit this convex
      number. 
    @*/
    size_type cv = get_optional_convex_number(in, gfi_fem, cmd);
    out.pop().from_scalar(fem->nb_dof(cv));
  } else if (check_cmd(cmd, "dim", in, out, 0, 0, 0, 1)) {
    /*@RDATTR FEM:GET('dim')
      Return the dimension (dimension of the reference convex) of the FEM.
      @*/
    out.pop().from_scalar(fem->dim());
  } else if (check_cmd(cmd, "target_dim", in, out, 0, 0, 0, 1)) {
    /*@RDATTR FEM:GET('target_dim')
      Return the dimension of the target space.

      The target space dimension is usually 1, except for vector FEMs
      (none of them has been implemented in getfem++ for now)@*/
    out.pop().from_scalar(fem->target_dim());
  } else if (check_cmd(cmd, "pts", in, out, 0, 1, 0, 1)) {
    /*@GET FEM:GET('pts' [, @int CV])
      Get the location of the degrees of freedom on the reference element.
      
      Some specific FEM may require a convex number CV to give their
      result (interpolated fems for example). In most of the case, you
      can omit this convex number.
      @*/
    size_type cv = get_optional_convex_number(in, gfi_fem, cmd);
    out.pop().from_vector_container(fem->node_convex(cv).points());
  } else if (check_cmd(cmd, "is_equivalent", in, out, 0, 0, 0, 1)) {
    /*@RDATTR FEM:GET('is_equivalent')
      Return 0 if the FEM is not equivalent.
      
      Equivalent FEM are evaluated on the reference convex. This is the
      case of most classical FEMs.@*/
    out.pop().from_scalar(fem->is_equivalent());
  } else if (check_cmd(cmd, "is_lagrange", in, out, 0, 0, 0, 1)) {
    /*@RDATTR FEM:GET('is_lagrange')
    Return 0 if the FEM is not of Lagrange type.@*/
    out.pop().from_scalar(fem->is_lagrange());
  } else if (check_cmd(cmd, "is_polynomial", in, out, 0, 0, 0, 1)) {
    /*@RDATTR FEM:GET('is_polynomial')
    Return 0 if the basis functions are not polynomials.@*/
    out.pop().from_scalar(fem->is_polynomial());
  } else if (check_cmd(cmd, "estimated_degree", in, out, 0, 0, 0, 1)) {
    /*@RDATTR FEM:GET('estimated_degree')
    Return an estimation of the polynomial degree of the FEM.@*/
    out.pop().from_scalar(fem->estimated_degree());
  } else if (check_cmd(cmd, "base_value", in, out, 1, 1, 0, 1)) {
    /*@GET FEM:GET('base_value',X)
    Evaluate all base functions at the given point.@*/
    getfem::base_tensor t;
    getfem::base_node x = in.pop().to_base_node(fem->dim());
    fem->base_value(x,t);
    out.pop().from_tensor(t);
  } else if (check_cmd(cmd, "grad_base_value", in, out, 1, 1, 0, 1)) {
    /*@GET FEM:GET('grad_base_value',X)
    Evaluate the gradient of all base functions at the given point.@*/
    getfem::base_tensor t;
    getfem::base_node x = in.pop().to_base_node(fem->dim());
    fem->grad_base_value(x,t);
    out.pop().from_tensor(t);
  } else if (check_cmd(cmd, "hess_base_value", in, out, 1, 1, 0, 1)) {
    /*@GET FEM:GET('hess_base_value',X)
    Evaluate the Hessian of all base functions at the given point.@*/
    getfem::base_tensor t;
    getfem::base_node x = in.pop().to_base_node(fem->dim());
    fem->hess_base_value(x,t);
    out.pop().from_tensor(t);
  } else if (check_cmd(cmd, "poly_str", in, out, 0, 0, 0, 1)) {
    /*@GET FEM:GET('poly_str')
      Return (if possible) the polynomial expressions of the base functions in the reference element.

      This function will of course fail for non-polynomial FEMs.
      @*/
    getfem::ppolyfem pf = dynamic_cast<getfem::ppolyfem>(&(*fem));
    if (pf) {
      std::vector<std::string> s(pf->base().size());
      for (size_type i=0; i < s.size(); ++i) {
	std::stringstream ss; ss << pf->base()[i];
	s[i] = ss.str();
      }
      out.pop().from_string_container(s);
    } else THROW_BADARG("Cannot return the poly_str of non-polynomial FEMs");
  } else if (check_cmd(cmd, "char", in, out, 0, 0, 0, 1)) {
    /*@GET FEM:GET('char')
    Ouput a (unique) string representation of the FEM.

    This can be used to perform comparisons between two different FEM
    objects.
    @*/
    std::string s = getfem::name_of_fem(fem);    
    out.pop().from_string(s.c_str());
  } else bad_cmd(cmd);
}
