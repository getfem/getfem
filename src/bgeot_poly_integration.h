/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_poly_integration.h : exact integration of polynomial   */
/*               on convexes of reference.                                 */
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


#ifndef __BGEOT_POLY_INTEGRATION_H
#define __BGEOT_POLY_INTEGRATION_H

#include <bgeot_poly.h>
#include <bgeot_convex_ref.h>

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

   /** @name functions on poly integration
   */
  //@{

  typedef const poly_integration * ppoly_integration;

  /// Give the exact integration method on a simplex of dimension n.
  ppoly_integration simplex_poly_integration(dim_type);
  /** Give the exact integration method corresponding to the tensorial
   *    product of *pp1 and *pp2.
   */
  ppoly_integration convex_product_poly_integration(ppoly_integration,
						    ppoly_integration);
  /// Give the exact integration method on a parallelepiped of dimension n.
  ppoly_integration parallelepiped_poly_integration(dim_type);
  /// Give the exact integration method on a prism of dimension n.
  inline ppoly_integration prism_poly_integration(dim_type nc)
  { return convex_product_poly_integration(simplex_poly_integration(nc-1),
					   simplex_poly_integration(1));
  }
    
  //@}

}  /* end of namespace bgeot.                                             */


#endif /* __BGEOT_POLY_INTEGRATION_H                                      */
