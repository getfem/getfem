/*===========================================================================

 Copyright (C) 2013-2018 Yves Renard

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/** @file   getfem_generic_assembly_tree.h
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @date   November 18, 2013.
    @brief  Definition of functions and basic operators of the assembly language
 */


#ifndef GETFEM_GENERIC_ASSEMBLY_FUNC_OP_H__
#define GETFEM_GENERIC_ASSEMBLY_FUNC_OP_H__

#include "getfem/getfem_generic_assembly_tree.h"

namespace getfem {

  struct ga_instruction_set;
  using instruction_set = omp_distribute<ga_instruction_set>;
  
  class ga_predef_function {
    size_type ftype_; // 0 : C++ function with C++ derivative(s)
                      // 1 : function defined by an string expression.

    size_type dtype_; // 0 : no derivative(s)
                      // 1 : derivative(s) given by C++ functions
                      // 2 : derivatives(s) given by string expression(s)
                      // 3 : derivatives(s) to be symbolically computed.
    size_type nbargs_;         // One or two arguments
    pscalar_func_onearg f1_;   // Function pointer for a one argument function
    pscalar_func_twoargs f2_;  // Function pointer for a two arguments function
    std::string expr_;
    std::string derivative1_, derivative2_;
    mutable omp_distribute<base_vector> t, u;
    mutable omp_distribute<ga_workspace> workspace;
    copyable_ptr<instruction_set> gis;

    friend void ga_define_function(const std::string &name, size_type nbargs,
                                   const std::string &expr,
                                   const std::string &der1,
                                   const std::string &der2);
    friend void ga_define_function(const std::string &name,
                                   pscalar_func_onearg f,
                                   const std::string &der);
    friend void ga_define_function(const std::string &name,
                                   pscalar_func_twoargs f,
                                   const std::string &der1,
                                   const std::string &der2);
  public:
    scalar_type operator()(scalar_type t_, scalar_type u_ = 0.) const;

    bool is_affine(const std::string &varname) const;

    size_type ftype() const { return ftype_;}
    size_type dtype() const { return dtype_;}
    size_type nbargs() const { return nbargs_;}
    const std::string &derivative1() const { return derivative1_;}
    const std::string &derivative2() const { return derivative2_;}
    const std::string &expr() const { return expr_;}
    pscalar_func_onearg f1() const { return f1_;}
    pscalar_func_twoargs f2() const { return f2_;}

    ga_predef_function();
    ga_predef_function(pscalar_func_onearg f, size_type dtype__ = 0,
                       const std::string &der = "");
    ga_predef_function(pscalar_func_twoargs f, size_type dtype__ = 0,
                       const std::string &der1 = "",
                       const std::string &der2 = "");
    ga_predef_function(const std::string &expr__);
  };

  struct ga_predef_function_tab
    : public std::map<std::string, ga_predef_function> {

    ga_predef_function_tab(); 
  };

  struct ga_spec_function_tab : public std::set<std::string> {
    ga_spec_function_tab();
  };

} /* end of namespace */


#endif /* GETFEM_GENERIC_ASSEMBLY_FUNC_OP_H__  */
