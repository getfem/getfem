// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2000-2008 Yves Renard
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

/**@file getfem_mat_elem.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date December 21, 2000.
   @brief elementary computations (used by the generic assembly).

   This is the kernel of getfem..
*/

#ifndef GETFEM_MAT_ELEM_H__
#define GETFEM_MAT_ELEM_H__

#include "getfem_mat_elem_type.h"
#include "getfem_fem.h"

namespace getfem {
  /** @internal (optional) callback to be called for each point of the
      integration (i.e. only with approximate integrations). It is
      used by getfem_assembling_tensors to perform reductions before
      integration.
  */
  struct mat_elem_integration_callback {
    /* this vector is filled by mat_elem, it contains a pointer to each tensor
       build from a term in the mat_elem_type (base functions, gradients,
       nonlinear terms etc)
    */
    std::vector<const bgeot::base_tensor*> eltm;
    /** executes the callback

        @param t the destination tensor

        @param first indicates if this is the first integration point
        (in that case, t should be set to the correct size and filled
        with zeros)

        @param c the current coefficient (contains the norm of
        Jacobian and the integration weight)
    */
    virtual void exec(bgeot::base_tensor &t, bool first, scalar_type c) = 0;
    virtual ~mat_elem_integration_callback() {}
  };

  /** @internal
      this class (whose instances are returned by the mat_elem
      function, see below) holds all computations of elementary
      integrals over convexes or faces of convexes. 
  */
  class mat_elem_computation : virtual public dal::static_stored_object {
  protected : 
    
    bgeot::pgeometric_trans pgt;
    pmat_elem_type pme;
    mutable base_matrix pa;
    
  public :
    
    virtual void compute(base_tensor &t, const base_matrix &a,
                         size_type elt, 
			 mat_elem_integration_callback *icb = 0) const = 0;
    virtual void compute_on_face(base_tensor &t, const base_matrix &a,
                                 short_type f, size_type elt, 
				 mat_elem_integration_callback *icb = 0)
      const = 0;
    /**
       perform the integration on the volume of a convex.
       @param t     the destination tensor
       @param a     the list of vertices of the convex
       @param elt   the convex number
       @param icb   an optional callback which will be called for each
                    integration point
    */
    template <class CONT> void 
    gen_compute(base_tensor &t, const CONT &a,  size_type elt, 
		mat_elem_integration_callback *icb = 0) const { 
      bgeot::vectors_to_base_matrix(pa, a); 
      compute(t, pa, elt, icb); 
    }
    /** 
        perform the integration on a face of the convex
    */
    template <class CONT> void 
    gen_compute_on_face(base_tensor &t,
                        const CONT &a, short_type f, size_type elt, 
			mat_elem_integration_callback *icb = 0) const {
      bgeot::vectors_to_base_matrix(pa, a); 
      compute_on_face(t, pa, f, elt, icb); 
    }
    
    virtual ~mat_elem_computation() {}
    virtual size_type memsize() const = 0;
  };

  typedef boost::intrusive_ptr<const mat_elem_computation>
  pmat_elem_computation;

  /** 
      allocate a structure for computation (integration over elements
      or faces of elements) of elementary tensors. Internally this
      structure is linked to a "cache" which stores some pre-computed
      data.
  */ 
  pmat_elem_computation mat_elem(pmat_elem_type pm, 
				 pintegration_method pi,
				 bgeot::pgeometric_trans pg,
                                 bool prefer_comp_on_real_element = false);


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MAT_ELEM_H__                                              */
