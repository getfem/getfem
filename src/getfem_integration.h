// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : getfem_integration.h : exact and approximative integration.
//           
// Date    : December 17, 2000.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2000-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================


#ifndef GETFEM_INTEGRATION_H__
#define GETFEM_INTEGRATION_H__

#include <getfem_config.h>
#include <bgeot_convex_ref.h>
#include <dal_tree_sorted.h>
#include <bgeot_geometric_trans.h>

namespace getfem
{
  /** Description of an exact integration of polynomials.
   *  This class is not to be manipulate by itself. Use ppoly_integration and
   *  the functions written to produce the basic descriptions.
   */
  class poly_integration
  {
    protected :

      bgeot::pconvex_structure cvs;

      mutable std::vector<long_scalar_type> int_monomials;
      mutable std::vector< std::vector<long_scalar_type> > int_face_monomials;

    public :

      /// Dimension of convex of reference.
      dim_type dim(void) const { return cvs->dim(); }
      /// {Structure of convex of reference.  
      bgeot::pconvex_structure structure(void) const { return cvs; }

      virtual long_scalar_type int_monomial(const bgeot::power_index &power)
	const = 0;
      virtual long_scalar_type int_monomial_on_face(const bgeot::power_index
					       &power, short_type f) const = 0;

      /// Evaluate the integral of the polynomial P on the reference element.
      long_scalar_type int_poly(const base_poly &P) const;
      /** Evaluate the integral of the polynomial P on the face f of the
       *    reference element.
       */
      long_scalar_type int_poly_on_face(const base_poly &P, short_type f) const;

      virtual ~poly_integration() {}

  };

  typedef const poly_integration *ppoly_integration;

   /** Description of an approximate integration of polynomials of
   *  several variables on a reference element.
   *  This class is not to be manipulate by itself. Use papprox\_integration
   *  and the functions written to produce the basic descriptions.
   */

  class integration_method;
  typedef const integration_method *pintegration_method;

  class approx_integration
  {
    protected :

      typedef dal::dynamic_tree_sorted<base_node,
      dal::lexicographical_less<base_node,
      dal::approx_less<scalar_type> > > PT_TAB;

      bgeot::pconvex_ref cvr;
      bgeot::pstored_point_tab pint_points;
      std::vector<scalar_type> int_coeffs;
      std::vector<size_type> repartition;

    // index 0 : points for volumic integration, index > 0 : points for faces
      std::vector<PT_TAB> pt_to_store; 
      bool valid;

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
      bgeot::pconvex_structure structure(void) const
        { return cvr->structure()->basic_structure(); }
      bgeot::pconvex_ref ref_convex(void) const { return cvr; }

      const std::vector<size_type> &repart(void) const { return repartition; }

      /// Gives an array of integration nodes.
      const bgeot::stored_point_tab  &
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

      void add_point(const base_node &pt, scalar_type w,
		     short_type f=short_type(-1));
      void add_point_norepeat(const base_node &pt, scalar_type w,
			      short_type f=short_type(-1));
      void add_point_full_symmetric(base_node pt, scalar_type w);
      void add_method_on_face(pintegration_method ppi, short_type f);
      void valid_method(void);

      approx_integration(void) : valid(false) { }
      approx_integration(bgeot::pconvex_ref cr)
	: cvr(cr), repartition(cr->structure()->nb_faces()+1),
	  pt_to_store(cr->structure()->nb_faces()+1), valid(false)
      { std::fill(repartition.begin(), repartition.end(), 0); } 
    
  };

  typedef const approx_integration * papprox_integration;

  /**
     the list of main integration method types 
  */
  typedef enum { IM_APPROX, IM_EXACT, IM_NONE } integration_method_type;

  /**
     this structure is not intended to be used directly. It is built via
     the int_method_descriptor() function
  */
  class integration_method {    
    union {
      ppoly_integration ppi; /* for exact integrations */
      papprox_integration pai; /* for approximate integrations (i.e. cubatures) */
    } method;
    integration_method_type im_type;

  public:
    integration_method_type type(void) const { return im_type; }
    papprox_integration approx_method(void) const { return method.pai; }
    ppoly_integration exact_method(void) const { return method.ppi; }

    const bgeot::stored_point_tab &integration_points(void) const { 
      if (type() == IM_EXACT)
	return *(bgeot::org_stored_point_tab(method.ppi->structure()->dim()));
      else if (type() == IM_APPROX)
	return method.pai->integration_points();
      else DAL_THROW(dal::failure_error, "IM_NONE has no points");
    }

    bgeot::pconvex_structure structure(void) const { 
      switch (type()) {
      case IM_EXACT: return method.ppi->structure();
      case IM_APPROX: return method.pai->structure();
      case IM_NONE: DAL_THROW(dal::failure_error, "IM_NONE has no structure");
	default: DAL_INTERNAL_ERROR("");
      }
    }

    integration_method(ppoly_integration p)
    { method.ppi = p; im_type = IM_EXACT; }

    integration_method(papprox_integration p)
    { method.pai = p; im_type = IM_APPROX; }

    integration_method(void) { im_type = IM_NONE; method.pai = 0; }
    ~integration_method(void) {
      switch (type()) {
      case IM_EXACT: if (method.ppi) delete method.ppi; break;
      case IM_APPROX: if (method.pai) delete method.pai; break;
      case IM_NONE: break;
      }
    }
  };


  pintegration_method int_method_descriptor(std::string name);
  /* List :
   * "IM_EXACT_SIMPLEX(n)"             : exact integration on simplexes.
   * "IM_PRODUCT(IM1, IM2)"            : product of two integration methods
   * "IM_EXACT_PARALLELEPIPED(n)"      : exact integration on parallelepipeds
   * "IM_EXACT_PRISM(n)"               : exact integration on prisms
   * "IM_GAUSS1D(K)"                   : Gauss method on the segment, order K
   * "IM_NC(N,K)"                      : Newton-Cotes approximative 
   *                                     integration on simplexes, order K
   * "IM_NC_PARALLELEPIPED(N,K)"       : product of Newton-Cotes integration 
   *                                     on parallelepipeds
   * "IM_NC_PRISM(N,K)"                : product of Newton-Cotes integration
   *                                     on prisms
   * "IM_GAUSS_PARALLELEPIPED(N,K)"    : product of Gauss1D integration
   *                                     on parallelepipeds
   * "IM_TRIANGLE(K)"                  : Gauss methods on triangles (K=1..7)
   * "IM_QUAD(K)"                      : Gauss methods on quadrilaterons
   *                                     (K=2, 3 or 5)
   * "IM_TETRAHEDRON(K)"               : Gauss methods on tetrahedrons
   *                                     (K=1, 2, 3 or 5)
   * "IM_STRUCTURED_COMPOSITE(IM1, K)" : Composite method on a grid with
   *                                     K divisions
   */
  
  pintegration_method classical_exact_im(bgeot::pgeometric_trans pgt);
  /**
     try to find an approximate integration method for the geometric
     transformation pgt which is able to integrate exactly polynomials
     of degree <= "degree". It may return a higher order integration
     method if no method match the exact degree.
  */
  pintegration_method classical_approx_im(bgeot::pgeometric_trans pgt, dim_type degree);

  pintegration_method exact_simplex_im(size_type n);
  pintegration_method exact_parallelepiped_im(size_type n);
  pintegration_method exact_prism_im(size_type n);
  pintegration_method exact_classical_im(bgeot::pgeometric_trans pgt) IS_DEPRECATED;

  std::string name_of_int_method(pintegration_method p);
  class mesh_precomposite;
  class mesh_fem;
  papprox_integration composite_approx_int_method(const mesh_precomposite &mp, 
						  const mesh_fem &mf,
						  bgeot::pconvex_ref cr);
}  /* end of namespace getfem.                                            */


#endif /* BGEOT_INTEGRATION_H__                                           */
