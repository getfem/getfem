/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_approx_integration.h : approximated  integration       */
/*               on convexes of reference                                  */
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


#ifndef __BGEOT_APPROX_INTEGRATION_H
#define __BGEOT_APPROX_INTEGRATION_H

#include <dal_tree_sorted.h>
#include <bgeot_convex_structure.h>
#include <bgeot_config.h>
#include <bgeot_poly_integration.h>


namespace bgeot
{

  /** Description of an approximate integration of polynomials of
   *  several variables on a reference element.
   *  This class is not to be manipulate by itself. Use papprox\_integration
   *  and the functions written to produce the basic descriptions.
   */
  class approx_integration
  {
    protected :

      pconvex_structure cvs;
      std::vector<base_node> int_points;
      std::vector<scalar_type> int_coeffs;
      std::vector<size_type> repartition;

      void optimize(void);

    public :

      /// Dimension of reference convex.
      dim_type dim(void) const { return cvs->dim(); }
      size_type nb_points(void) const { return int_coeffs.size(); }
      /// Number of integration nodes on the reference element.
      size_type nb_points_on_convex(void) const { return repartition[0]; }
      /// Number of integration nodes on the face f of the reference element.
      size_type nb_points_on_face(short_type f) const
      { return repartition[f+1] - repartition[f]; }
      /// Structure of the reference element.
      pconvex_structure structure(void) const { return cvs; }

      /// Gives an array of integration nodes.
      const std::vector<base_node>  &
        integration_points(void) const
      { return int_points; }
      /// Gives the integration node i on the reference element.
      const base_node &point(size_type i) const
      { return int_points[i]; }
      /// Gives the integration node i of the face f.
      const base_node &
	point_on_face(short_type f, size_type i) const 
      { return int_points[repartition[f] + i]; }
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
   *    nbpt integration nodes.
   */
  papprox_integration Gauss_approx_integration(short_type nbpt);

  // methods on dimension 2

  // ... a rentrer

  // methods on dimension 3

  // ... a rentrer

  // methods on simplex

  /** Give the Newton Cotes approximate integration method in simplex of
   *    dimension n and of degree k (corresponds to lagrange interpolation).
   */
  papprox_integration Newton_Cotes_approx_integration(dim_type n,short_type k);

  /** Give the approximate integration method corresponding to the
   *    tensorial product of *pai1 and *pai2.
   */
  papprox_integration convex_product_approx_integration(papprox_integration,
							papprox_integration);
    

  //@}

}  /* end of namespace bgeot.                                              */


#endif /* __BGEOT_APPROX_INTEGRATION_H                                     */
