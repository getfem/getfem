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

#include <getfemint.h>
#include <getfem/getfem_fem.h>

using namespace getfemint;

static size_type get_optional_convex_number(getfemint::mexargs_in &in,
					    const getfem::pfem &pf,
					    const std::string cmd) {
  size_type cv = 0;
  if (!in.remaining() && pf->is_on_real_element())
    THROW_BADARG("This FEM requires a convex number for " << cmd);
  if (in.remaining())
    cv = in.pop().to_integer() - config::base_index();
  return cv;
}

/*@GFDOC
  General function for querying information about FEM objects.
@*/



// Object for the declaration of a new sub-command
struct sub_gf_fem_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   const getfem::pfem &pf) = 0;
};

typedef std::shared_ptr<sub_gf_fem_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_fem_get {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       const getfem::pfem &pf)				\
      { dummy_func(in); dummy_func(out); code }				\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }                           




void gf_fem_get(getfemint::mexargs_in& m_in, getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {


    /*@RDATTR n = ('nbdof'[, @int cv])
    Return the number of dof for the @tfem.

    Some specific @tfem (for example 'interpolated_fem') may require a
    convex number `cv` to give their result. In most of the case, you
    can omit this convex number.@*/
    sub_command
      ("nbdof", 0, 1, 0, 1,
       size_type cv = get_optional_convex_number(in, pf, "nbdof");
       out.pop().from_scalar(double(pf->nb_dof(cv)));
       );

    /*@RDATTR n = ('index of global dof', cv)
    Return the index of global dof for special fems such as interpolated fem.
    @*/
    sub_command
      ("index of global dof", 2, 2, 0, 1,
       size_type cv = in.pop().to_integer() - config::base_index();
       size_type i = in.pop().to_integer() - config::base_index();
       out.pop().from_scalar(double(pf->index_of_global_dof(cv, i) + config::base_index()));
       );


    /*@RDATTR d = ('dim')
      Return the dimension (dimension of the reference convex) of the @tfem.@*/
    sub_command
      ("dim", 0, 0, 0, 1,
       out.pop().from_scalar(pf->dim());
       );


    /*@RDATTR td = ('target_dim')
    Return the dimension of the target space.

    The target space dimension is usually 1, except for vector @tfem. @*/
    sub_command
      ("target_dim", 0, 0, 0, 1,
       out.pop().from_scalar(pf->target_dim());
       );


    /*@GET P = ('pts'[, @int cv])
      Get the location of the dof on the reference element.
      
      Some specific @tfem may require a convex number `cv` to give their
      result (for example 'interpolated_fem'). In most of the case, you
      can omit this convex number. @*/
    sub_command
      ("pts", 0, 1, 0, 1,
       size_type cv = get_optional_convex_number(in, pf, "pts");
       out.pop().from_vector_container(pf->node_convex(cv).points());
       );

    /*@RDATTR b = ('is_equivalent')
      Return 0 if the @tfem is not equivalent.
      
      Equivalent @tfem are evaluated on the reference convex. This is
      the case of most classical @tfem's.@*/
    sub_command
      ("is_equivalent", 0, 0, 0, 1,
       out.pop().from_scalar(pf->is_equivalent());
       );


    /*@RDATTR b = ('is_lagrange')
      Return 0 if the @tfem is not of Lagrange type.@*/
    sub_command
      ("is_lagrange", 0, 0, 0, 1,
       out.pop().from_scalar(pf->is_lagrange());
       );


    /*@RDATTR b = ('is_polynomial')
      Return 0 if the basis functions are not polynomials.@*/
    sub_command
      ("is_polynomial", 0, 0, 0, 1,
       out.pop().from_scalar(pf->is_polynomial());
       );


    /*@RDATTR d = ('estimated_degree')
    Return an estimation of the polynomial degree of the @tfem.

    This is an estimation for fem which are not polynomials.@*/
    sub_command
      ("estimated_degree", 0, 0, 0, 1,
       out.pop().from_scalar(pf->estimated_degree());
       );


    /*@GET E = ('base_value',@mat p)
      Evaluate all basis functions of the FEM at point `p`.
      
      `p` is supposed to be in the reference convex!@*/
    sub_command
      ("base_value", 1, 1, 0, 1,
       getfem::base_tensor t;
       getfem::base_node x = in.pop().to_base_node(pf->dim());
       pf->base_value(x,t);
       out.pop().from_tensor(t);
       );


    /*@GET ED = ('grad_base_value',@mat p)
    Evaluate the gradient of all base functions of the @tfem at point `p`.

    `p` is supposed to be in the reference convex!@*/
    sub_command
      ("grad_base_value", 1, 1, 0, 1,
       getfem::base_tensor t;
       getfem::base_node x = in.pop().to_base_node(pf->dim());
       pf->grad_base_value(x,t);
       out.pop().from_tensor(t);
       );


    /*@GET EH = ('hess_base_value',@mat p)
    Evaluate the Hessian of all base functions of the @tfem at point `p`.

    `p` is supposed to be in the reference convex!@*/
    sub_command
      ("hess_base_value", 1, 1, 0, 1,
       getfem::base_tensor t;
       getfem::base_node x = in.pop().to_base_node(pf->dim());
       pf->hess_base_value(x,t);
       out.pop().from_tensor(t);
       );


    /*@GET ('poly_str')
      Return the polynomial expressions of its basis functions in
      the reference convex.

      The result is expressed as a @MATLAB{cell array}@SCILAB{cell array}@PYTHON{tuple} of
      strings. Of course this will fail on non-polynomial @tfem's. @*/
    sub_command
      ("poly_str", 0, 0, 0, 1,
       getfem::ppolyfem ppf = dynamic_cast<getfem::ppolyfem>(pf.get());
       if (ppf) {
	 std::vector<std::string> s(ppf->base().size());
	 for (size_type i=0; i < s.size(); ++i) {
	   std::stringstream ss; ss << ppf->base()[i];
	   s[i] = ss.str();
	 }
	 out.pop().from_string_container(s);
       }
       else THROW_BADARG("Cannot return the poly_str of non-polynomial FEMs");
       );


    /*@GET @str = ('char')
    Ouput a (unique) string representation of the @tfem.

    This can be used to perform comparisons between two different @tfem
    objects.@*/
    sub_command
      ("char", 0, 0, 0, 1,
       std::string s = getfem::name_of_fem(pf);
       out.pop().from_string(s.c_str());
       );


    /*@GET ('display')
    displays a short summary for a @tfem object.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       infomsg() << "gfFem object " << getfem::name_of_fem(pf)
       << " in dimension " << int(pf->dim())
       << ", with target dim " << int(pf->target_dim()) << " dof number "
       << pf->nb_dof(0);
       if (pf->is_equivalent()) infomsg() << " EQUIV ";
       else infomsg() << " NOTEQUIV ";
       if (pf->is_polynomial()) infomsg() << " POLY ";
       else infomsg() << " NOTPOLY ";
       if (pf->is_lagrange()) infomsg() << " LAGRANGE ";
       else infomsg() << " NOTLAGRANGE ";
       infomsg() << endl;
       );

  }



  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");


  getfem::pfem pf = to_fem_object(m_in.pop());
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out, pf);
  }
  else bad_cmd(init_cmd);



}
