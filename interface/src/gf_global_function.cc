/*===========================================================================

 Copyright (C) 2009-2020 Luis Saavedra.

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
// $Id$
#include <getfemint.h>
#include <getfemint_workspace.h>
#include <getfem/getfem_global_function.h>
#include <getfem/getfem_arch_config.h>

using namespace getfemint;

/*@GFDOC
  Global function object is represented by three functions:

   * The function `val`.
   * The function gradient `grad`.
   * The function Hessian `hess`.

  this type of function is used as local and global enrichment function. The
  global function Hessian is an optional parameter (only for fourth order
  derivative problems). @*/


// Object for the declaration of a new sub-command.

struct sub_gf_globfunc : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   getfem::pxy_function &ggf) = 0;
};

typedef std::shared_ptr<sub_gf_globfunc> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_globfunc {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       getfem::pxy_function &ggf)			\
      { dummy_func(in); dummy_func(out); code }				\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }


void gf_global_function(getfemint::mexargs_in& m_in,
			getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {


    /*@INIT GF = ('cutoff', @int fn, @scalar r, @scalar r1, @scalar r0)
      Create a cutoff global function.@*/
    sub_command
      ("cutoff", 4, 4, 0, 1,
       size_type   fn = in.pop().to_integer(-1,2);
       scalar_type r  = in.pop().to_scalar();
       scalar_type r1 = in.pop().to_scalar();
       scalar_type r0 = in.pop().to_scalar();

       ggf = std::make_shared<getfem::cutoff_xy_function>(int(fn),r,r1,r0);
       );


    /*@INIT GF = ('crack', @int fn)
      Create a near-tip asymptotic global function for modelling cracks.@*/
    sub_command
      ("crack", 1, 1, 0, 1,
       size_type fn = in.pop().to_integer(0,11);
       ggf = std::make_shared<getfem::crack_singular_xy_function>(unsigned(fn));
       );

    /*@INIT GF = ('parser', @str val[, @str grad[, @str hess]])
      Create a global function from strings `val`, `grad` and `hess`.
      This function could be improved by using the derivation of the generic
      assembly language ... to be done. @*/
    sub_command
      ("parser", 1, 3, 0, 1,
       std::string sval = in.pop().to_string();
       std::string sgrad = "[0;0]";
       std::string shess = "[0,0;0,0]";
       if (in.remaining() && in.front().is_string())
	 sgrad = in.pop().to_string();
       if (in.remaining() && in.front().is_string())
	 shess = in.pop().to_string();
       ggf = std::make_shared<getfem::parser_xy_function>(sval,sgrad,shess);
       );

    /*@INIT GF = ('product', @tgf F, @tgf G)
      Create a product of two global functions.@*/
    sub_command
      ("product", 2, 2, 0, 1,
       getfem::pxy_function af1 = to_global_function_object(in.pop());
       getfem::pxy_function af2 = to_global_function_object(in.pop());
       ggf = std::make_shared<getfem::product_of_xy_functions>(af1,af2);
       );


    /*@INIT GF = ('add', @tgf gf1, @tgf gf2)
      Create a add of two global functions.@*/
    sub_command
      ("add", 2, 2, 0, 1,
       getfem::pxy_function af1 = to_global_function_object(in.pop());
       getfem::pxy_function af2 = to_global_function_object(in.pop());
       ggf = std::make_shared<getfem::add_of_xy_functions>(af1,af2);
       );

  }



  if (m_in.narg() < 1)  THROW_BADARG( "Wrong number of input arguments");

  getfem::pxy_function ggf;
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);


  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out, ggf);
  }
  else bad_cmd(init_cmd);

  m_out.pop().from_object_id(store_global_function_object(ggf),
			     GLOBAL_FUNCTION_CLASS_ID);
}


/*@PYTHONEXT
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
@*/
