/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_integration.h : exact and approximative integration.   */
/*                                                                         */
/* Date : December 17, 2000.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; version 2 of the License.                 */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software Foundation, */
/* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.         */
/*                                                                         */
/* *********************************************************************** */


#ifndef __BGEOT_INTEGRATION_H
#define __BGEOT_INTEGRATION_H

#include <bgeot_poly.h>
#include <bgeot_convex_ref.h>
#include <dal_tree_sorted.h>

namespace bgeot
{
  /** Description of an exact integration of polynomials.
   *  This class is not to be manipulate by itself. Use ppoly_integration and
   *  the functions written to produce the basic descriptions.
   */
  class poly_integration
  {
    protected :

      pconvex_structure cvs;

      std::vector<scalar_type> int_monomials;
      std::vector< std::vector<scalar_type> > int_face_monomials;

    public :

      /// Dimension of convex of reference.
      dim_type dim(void) const { return cvs->dim(); }
      /// {Structure of convex of reference.  
      pconvex_structure structure(void) const { return cvs; }

      virtual scalar_type int_monomial(const power_index &power) const = 0;
      virtual scalar_type int_monomial_on_face(const power_index &power, 
					       short_type f) const = 0;

      /// Evaluate the integral of the polynomial P on the reference element.
      scalar_type int_poly(const polynomial<scalar_type> &P) const;
      /** Evaluate the integral of the polynomial P on the face f of the
       *    reference element.
       */
      scalar_type int_poly_on_face(const polynomial<scalar_type> &P,
				   short_type f) const;

      virtual ~poly_integration() {}

  };

  typedef const poly_integration *ppoly_integration;

   /** Description of an approximate integration of polynomials of
   *  several variables on a reference element.
   *  This class is not to be manipulate by itself. Use papprox\_integration
   *  and the functions written to produce the basic descriptions.
   */
  class approx_integration
  {
    protected :

      pconvex_ref cvr;
      pstored_point_tab pint_points;
      std::vector<scalar_type> int_coeffs;
      std::vector<size_type> repartition;

    public :

      /// Dimension of reference convex.
      dim_type dim(void) const { return cvr->structure()->dim(); }
      size_type nb_points(void) const { return int_coeffs.size(); }
      /// Number of integration nodes on the reference element.
      size_type nb_points_on_convex(void) const { return repartition[0]; }
      /// Number of integration nodes on the face f of the reference element.
      size_type nb_points_on_face(short_type f) const
      { return repartition[f+1] - repartition[f]; }
      /// Structure of the reference element.
      pconvex_structure structure(void) const
        { return cvr->structure()->basic_structure(); }
      pconvex_ref ref_convex(void) const { return cvr; }

      /// Gives an array of integration nodes.
      const stored_point_tab  &
        integration_points(void) const
      { return *(pint_points); }
      /// Gives the integration node i on the reference element.
      const base_node &point(size_type i) const
      { return (*pint_points)[i]; }
      /// Gives the integration node i of the face f.
      const base_node &
	point_on_face(short_type f, size_type i) const 
      { return (*pint_points)[repartition[f] + i]; }
      /// Gives an array of the integration coefficients.
      const std::vector<scalar_type> &integration_coefficients(void) const
      { return int_coeffs; }
      /// Gives the integration coefficient corresponding to node i.
      scalar_type coeff(size_type i) const { return int_coeffs[i]; }
      /// Gives the integration coefficient corresponding to node i of face f.
      scalar_type coeff_on_face(short_type f, size_type i) const
      { return int_coeffs[repartition[f] + i]; }
    
  };

  /** @name functions on approx integration
   */
  //@{

  typedef const approx_integration * papprox_integration;

  // methods on dimension 1
  
  /** Give the Gauss approximate integration method in dimension 1 with
   *    nbpt integration nodes. Exact on polynomial of degree nbpt * 2 - 1.
   */
  papprox_integration Gauss_approx_integration(short_type nbpt);
  papprox_integration parallelepiped_Gauss_approx_integration(dim_type N,
							      short_type nbpt);

  // methods on dimension 2

  /// Integration on a triangle of order 1 with 1 point
  papprox_integration triangle1_approx_integration(void);
  /// Integration on a triangle of order 2 with 3 points
  papprox_integration triangle2_approx_integration(void);
  /// Integration on a triangle of order 2 with 3 points
  papprox_integration triangle2bis_approx_integration(void);
  /// Integration on a triangle of order 3 with 4 points
  papprox_integration triangle3_approx_integration(void);
  /// Integration on a triangle of order 4 with 6 points
  papprox_integration triangle4_approx_integration(void);
  /// Integration on a triangle of order 5 with 7 points
  papprox_integration triangle5_approx_integration(void);
  /// Integration on a triangle of order 6 with 12 points
  papprox_integration triangle6_approx_integration(void);
  /// Integration on a triangle of order 7 with 13 points
  papprox_integration triangle7_approx_integration(void);
  /// Integration on quadrilaterals of order 2 with 3 points   
  papprox_integration quad2_approx_integration(void);
  /// Integration on quadrilaterals of order 3 with 4 points   
  papprox_integration quad3_approx_integration(void);    
  /// Integration on quadrilaterals of order 5 with 7 points   
  papprox_integration quad5_approx_integration(void);


  // methods on dimension 3

  /// Integration on a tetrahedron of order 1 with 1 point
  papprox_integration tetrahedron1_approx_integration(void);
  /// Integration on a tetrahedron of order 2 with 4 points
  papprox_integration tetrahedron2_approx_integration(void);
  /// Integration on a tetrahedron of order 3 with 5 points
  papprox_integration tetrahedron3_approx_integration(void);
  /// Integration on a tetrahedron of order 5 with 15 points
  papprox_integration tetrahedron5_approx_integration(void);
 
  /** Give the approximate integration method corresponding to the
   *    tensorial product of *pai1 and *pai2.
   */
  papprox_integration convex_product_approx_integration(papprox_integration,
							papprox_integration);

  // methods on simplexes

  /** Give the Newton Cotes approximate integration method in simplex of
   *    dimension n and of degree k (corresponds to lagrange interpolation).
   */
  papprox_integration Newton_Cotes_approx_integration(dim_type n,short_type k);
  papprox_integration parallelepiped_Newton_Cotes_approx_integration
  (dim_type N, short_type k);
  inline papprox_integration prism_Newton_Cotes_approx_integration(dim_type N,
							       short_type k) {
    return convex_product_approx_integration
      (Newton_Cotes_approx_integration(N - 1, k),
       Newton_Cotes_approx_integration(1, k));
  }

  //@}


  struct integration_method
  {
    union
    {
      ppoly_integration ppi;
      papprox_integration pai;
    } method;
    bool is_ppi;

    const stored_point_tab &integration_points(void) const { 
      if (is_ppi)
	return *(org_stored_point_tab(method.ppi->structure()->dim()));
      else 
	return method.pai->integration_points();
    }

    pconvex_structure structure(void) const { 
      if (is_ppi) return method.ppi->structure();
      else return method.pai->structure();
    }

    integration_method(ppoly_integration p)
    { method.ppi = p; is_ppi = true; }

    integration_method(papprox_integration p)
    { method.pai = p; is_ppi = false; }

    bool operator >(const integration_method& p) const
    { return method.ppi > p.method.ppi; } 
    bool operator <(const integration_method& p) const
    { return method.ppi < p.method.ppi; } 
    bool operator !=(const integration_method& p) const
    { return method.ppi != p.method.ppi; } 
    bool operator ==(const integration_method& p) const
    { return method.ppi == p.method.ppi; } 

    integration_method(void) { method.pai = 0; }
    ~integration_method(void) {
      if (method.pai != 0)
	if (is_ppi) delete method.ppi; else delete method.pai;
    }

  };

  typedef const integration_method *pintegration_method;


   pintegration_method int_method_descriptor(std::string name);
   /* List :
    * IM_EXACT_SIMPLEX(n)        : exact integration on simplexes.
    * IM_PRODUCT(a, b)           : product of two integration methods
    * IM_EXACT_PARALLELEPIPED(n) : exact integration on parallelepipeds
    * IM_EXACT_PRISM(n)          : exact integration on prisms
    * IM_GAUSS1D(K)       : Gauss method on the segment, order K
    * IM_NC(N,K)          : Newton-Cotes approximative integration 
    *                       on simplexes, order K
    * IM_NC_PARALLELEPIPED(N,K) :  product of Newton-Cotes integration 
    *                              on parallelepipeds
    * IM_NC_PRISM(N,K)    : product of Newton-Cotes integration on prisms
    * IM_GAUSS_PARALLELEPIPED(N,K) :  product of Gauss1D integration
    *                                 on parallelepipeds
    * IM_TRIANGLE(K)      : Gauss methods on triangles (K=1..7)
    * IM_QUAD(K)          : Gauss methods on quadrilaterons (K=2, 3 or 5)
    * IM_TETRAHEDRON(K)   : Gauss methods on tetrahedrons (K=1, 2, 3 or 5)
    */

   pintegration_method exact_simplex_im(size_type n);
  pintegration_method exact_parallelepiped_im(size_type n);
  pintegration_method exact_prism_im(size_type n);

   std::string name_of_int_method(pintegration_method p);

}  /* end of namespace bgeot.                                             */


#endif /* __BGEOT_INTEGRATION_H                                           */
