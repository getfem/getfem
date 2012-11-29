/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2011-2012 Yves Renard, Konstantinos Poulios.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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

/** @file getfem_contact_and_friction_common.h
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @author Konstantinos Poulios <logari81@googlemail.com>
    @date November, 2011.
    @brief Comomon tools for unilateral contact and Coulomb friction bricks.
 */
#ifndef GETFEM_CONTACT_AND_FRICTION_COMMON_H__
#define GETFEM_CONTACT_AND_FRICTION_COMMON_H__

#include "getfem_models.h"
#include "getfem_assembling_tensors.h"

namespace getfem {

  //=========================================================================
  //
  //  Projection on a ball and gradient of the projection.
  //
  //=========================================================================

  template<typename VEC> void ball_projection(const VEC &x,
					      scalar_type radius) {
    scalar_type a = gmm::vect_norm2(x);
    if (radius <= 0) gmm::clear(const_cast<VEC&>(x));
    else if (a > radius) gmm::scale(const_cast<VEC&>(x), radius/a);
  }
  
  template<typename VEC, typename VECR>
  void ball_projection_grad_r(const VEC &x, scalar_type radius,
                                     VECR &g) {
    scalar_type a = gmm::vect_norm2(x);
    if (radius > 0 && a >= radius) {
      gmm::copy(x, g); gmm::scale(g, scalar_type(1)/a);
    }
    else gmm::clear(g);
  }

  template <typename VEC, typename MAT>
  void ball_projection_grad(const VEC &x, double radius, MAT &g) {
    if (radius <= scalar_type(0)) { gmm::clear(g); return; }
    gmm::copy(gmm::identity_matrix(), g);
    scalar_type a = gmm::vect_norm2(x);
    if (a >= radius) {
      gmm::scale(g, radius/a);
      // gmm::rank_one_update(g, gmm::scaled(x, -radius/(a*a*a)), x);
      for (size_type i = 0; i < x.size(); ++i)
        for (size_type j = 0; j < x.size(); ++j)
          g(i,j) -= radius*x[i]*x[j] / (a*a*a);
    }
  }



  template <typename VEC, typename VECR>
  void coupled_projection(const VEC &x, const VEC &n,
			  scalar_type f, VECR &g) {
    scalar_type xn = gmm::vect_sp(x, n);
    scalar_type xnm = gmm::neg(xn);
    scalar_type th = f * xnm;
    scalar_type xtn = gmm::sqrt(gmm::vect_norm2_sqr(x) - xn*xn);
    
    gmm::copy(gmm::scaled(n, -xnm), g);
    if (th > scalar_type(0)) {
      if (xtn <= th) {
	gmm::add(x, g);
	gmm::add(gmm::scaled(n, -xn), g);
      } else {
	gmm::add(gmm::scaled(x, f*xnm/xtn), g);
	gmm::add(gmm::scaled(n, -f*xnm*xn/xtn), g);
      }
    }
  }


  template <typename VEC, typename MAT>
  void coupled_projection_grad(const VEC &x, const VEC &n,
			       scalar_type f, MAT &g) {
    scalar_type xn = gmm::vect_sp(x, n);
    scalar_type xnm = gmm::neg(xn);
    scalar_type th = f * xnm;
    scalar_type xtn = gmm::sqrt(gmm::vect_norm2_sqr(x) - xn*xn);
    size_type N = gmm::vect_size(x);
    gmm::clear(g);
    
    if (th > scalar_type(0)) {
      if (xtn <= th) {
	gmm::copy(gmm::identity_matrix(), g);
	gmm::rank_one_update(g, gmm::scaled(n, -scalar_type(1)), n);
      } else if (xn < scalar_type(0)) {
	static base_small_vector t; gmm::resize(t, N);
	gmm::add(x, gmm::scaled(n, -xn), t);
        gmm::scale(t, scalar_type(1)/xtn);
	if (N > 2) {
	  gmm::copy(gmm::identity_matrix(), g);
	  gmm::rank_one_update(g, gmm::scaled(t, -scalar_type(1)), t);
	  gmm::rank_one_update(g, gmm::scaled(n, -scalar_type(1)), n);
	  gmm::scale(g, -xn*th/xtn);
        }
        gmm::rank_one_update(g, gmm::scaled(t, -f), n);
      }
    }

    if (xn < scalar_type(0)) gmm::rank_one_update(g, n, n);
  }

  //=========================================================================
  //
  //  De Saxce projection and its gradients.
  //
  //=========================================================================


  template<typename VEC>
  void De_Saxce_projection(const VEC &x, const VEC &n_, scalar_type f) {
    static base_small_vector n; // For more robustness, n_ is not supposed unitary
    size_type N = gmm::vect_size(x);
    gmm::resize(n, N);
    gmm::copy(gmm::scaled(n_, scalar_type(1)/gmm::vect_norm2(n_)), n);
    scalar_type xn = gmm::vect_sp(x, n);
    scalar_type nxt = sqrt(gmm::abs(gmm::vect_norm2_sqr(x) - xn*xn));
    if (xn >= scalar_type(0) && f * nxt <= xn) {
      gmm::clear(const_cast<VEC&>(x));
    } else if (xn > scalar_type(0) || nxt > -f*xn) {
      gmm::add(gmm::scaled(n, -xn), const_cast<VEC&>(x));
      gmm::scale(const_cast<VEC&>(x), -f / nxt);
      gmm::add(n, const_cast<VEC&>(x));
      gmm::scale(const_cast<VEC&>(x), (xn - f * nxt) / (f*f+scalar_type(1)));
    }
  }

  template<typename VEC, typename MAT>
  void De_Saxce_projection_grad(const VEC &x, const VEC &n_,
				scalar_type f, MAT &g) {
    static base_small_vector n;
    size_type N = gmm::vect_size(x);
    gmm::resize(n, N);
    gmm::copy(gmm::scaled(n_, scalar_type(1)/gmm::vect_norm2(n_)), n);
    scalar_type xn = gmm::vect_sp(x, n);
    scalar_type nxt = sqrt(gmm::abs(gmm::vect_norm2_sqr(x) - xn*xn));


    if (xn > scalar_type(0) && f * nxt <= xn) {
      gmm::clear(g);
    } else if (xn > scalar_type(0) || nxt > -f*xn) {
      static base_small_vector xt;
      gmm::resize(xt, N);
      gmm::add(x, gmm::scaled(n, -xn), xt);
      gmm::scale(xt, scalar_type(1)/nxt);

      if (N > 2) {
        gmm::copy(gmm::identity_matrix(), g);
        gmm::rank_one_update(g, gmm::scaled(n, -scalar_type(1)), n);
        gmm::rank_one_update(g, gmm::scaled(xt, -scalar_type(1)), xt);
        gmm::scale(g, f*(f - xn/nxt));
      } else {
        gmm::clear(g);
      }

      gmm::scale(xt, -f); gmm::add(n, xt);
      gmm::rank_one_update(g, xt, xt);
      gmm::scale(g, scalar_type(1) / (f*f+scalar_type(1)));
    } else {
      gmm::copy(gmm::identity_matrix(), g);
    }
  }


  template<typename VEC, typename MAT>
  static void De_Saxce_projection_gradn(const VEC &x, const VEC &n_,
					scalar_type f, MAT &g) {
    static base_small_vector n;
    size_type N = gmm::vect_size(x);
    scalar_type nn = gmm::vect_norm2(n_);
    gmm::resize(n, N);
    gmm::copy(gmm::scaled(n_, scalar_type(1)/nn), n);
    scalar_type xn = gmm::vect_sp(x, n);
    scalar_type nxt = sqrt(gmm::abs(gmm::vect_norm2_sqr(x) - xn*xn));
    gmm::clear(g);

    if (!(xn > scalar_type(0) && f * nxt <= xn)
	&& (xn > scalar_type(0) || nxt > -f*xn)) {
      static base_small_vector xt, aux;
      gmm::resize(xt, N); gmm::resize(aux, N);
      gmm::add(x, gmm::scaled(n, -xn), xt);
      gmm::scale(xt, scalar_type(1)/nxt);

      scalar_type c = (scalar_type(1) + f*xn/nxt)/nn;
      for (size_type i = 0; i < N; ++i) g(i,i) = c;
      gmm::rank_one_update(g, gmm::scaled(n, -c), n);
      gmm::rank_one_update(g, gmm::scaled(n, f/nn), xt);
      gmm::rank_one_update(g, gmm::scaled(xt, -f*xn/(nn*nxt)), xt);
      gmm::scale(g, xn - f*nxt);

      gmm::add(gmm::scaled(xt, -f), n, aux);
      gmm::rank_one_update(g, aux, gmm::scaled(xt, (nxt+f*xn)/nn));
     
      gmm::scale(g, scalar_type(1) / (f*f+scalar_type(1)));
    }
  }


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONTACT_AND_FRICTION_COMMON_H__ */
