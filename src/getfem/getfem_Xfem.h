// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2003-2008 Yves Renard
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
// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.
//
//===========================================================================

/**@file getfem_Xfem.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date April 8, 2003.
   @brief eXtended Finite Element Method.

   what is an Xfem ?

   It is a "base" fem (for example PK(2,2)), with additional base
   functions.
   These additionnal base functions are the product of:
   - a global function (the virtual_Xfem_func)
   - base functions of another fem (for example PK(2,1))

   The Xfem is built using the add_func member function, which takes
   as parameters, a global function and a fem.

   example of use: enrichment of the finite elements space with
   particular functions, which may represent discontinuities of the
   field in some elements, or singularities in the field (crack tip
   functions..)
*/

#include "getfem_fem.h"
#include "getfem_mesh_fem.h"
/* Works only for tau-equivalent elements                                  */

namespace getfem
{
  struct Xfem_func_context {
    base_node xreal, xref;
    pfem pf;
    bgeot::pgeometric_trans pgt;
    size_type base_num; /* number of the current base function of pf */
    const base_matrix& G;
    Xfem_func_context(const fem_interpolation_context &c) :
      xreal(c.xreal()), xref(c.xref()), pgt(c.pgt()), G(c.G()) {}
  };

  // Object representing global functions. To be derived.

  struct virtual_Xfem_func {
    /*
    */
    virtual scalar_type val(const Xfem_func_context&)
    { GMM_ASSERT1(false, "this Xfem_func has no value"); }
    virtual base_small_vector grad(const Xfem_func_context&)
    { GMM_ASSERT1(false, "this Xfem_func has no gradient"); }  
    virtual base_matrix hess(const Xfem_func_context&)
    { GMM_ASSERT1(false, "this Xfem_func has no hessian"); }
    virtual ~virtual_Xfem_func() {}
  };
  typedef virtual_Xfem_func *pXfem_func;
  
  // Xfem definition
  
  class Xfem : public virtual_fem {
  
  protected:
    pfem pfb; // base fem
    std::vector<pfem> uniq_pfe; 
    std::vector<size_type> func_pf; // nb_func fems which are enriched (indexes in the array uniq_pfe)
    bool is_valid;
    size_type nb_func;
    std::vector<pXfem_func> funcs; // List of functions to be added
    std::vector<size_type> func_indices;

    void get_fem_interpolation_context_tab(
		 const fem_interpolation_context& c0,
		 std::vector<fem_interpolation_context>& vc) const;
    pfem pfe(size_type k) const { return uniq_pfe[func_pf[k]]; }
  public:

    void valid(void);

    virtual size_type nb_dof(size_type) const;

    /* ind should be > 0 */
    void add_func(pfem pf, pXfem_func pXf,
		  size_type ind = size_type(-1));
    
    void base_value(const base_node &x, base_tensor &t) const;
    void grad_base_value(const base_node &x, base_tensor &t) const;
    void hess_base_value(const base_node &x, base_tensor &t) const;

    void real_base_value(const fem_interpolation_context& c, 
			 base_tensor &t, bool = true) const;    
    void real_grad_base_value(const fem_interpolation_context& c, 
			      base_tensor &t, bool = true) const;
    void real_hess_base_value(const fem_interpolation_context& c, 
			      base_tensor &t, bool = true) const;
    void init(pfem pf);
    Xfem(pfem pf);
  };


}  /* end of namespace getfem.                                            */
