/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2000-2015 Yves Renard

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

/**@file getfem_integration.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date December 17, 2000.
   @brief Integration methods (exact and approximated) on convexes

   @section im_list List of integration methods for getfem::int_method_descriptor

   - "IM_EXACT_SIMPLEX(n)"             : exact integration on simplexes.
   - "IM_PRODUCT(IM1, IM2)"            : product of two integration methods
   - "IM_EXACT_PARALLELEPIPED(n)"      : exact integration on parallelepipeds
   - "IM_EXACT_PRISM(n)"               : exact integration on prisms
   - "IM_GAUSS1D(K)"                   : Gauss method on the segment, 
                                         order K=1, 3, 5, 7, ..., 99.
   - "IM_GAUSSLOBATTO1D(K)"            : Gauss-Lobatto method on the segment, 
                                         order K=1, 3, 5, 7, ..., 99.
   - "IM_NC(N,K)"                      : Newton-Cotes approximative 
                                         integration on simplexes, order K
   - "IM_NC_PARALLELEPIPED(N,K)"       : product of Newton-Cotes integration 
                                         on parallelepipeds
   - "IM_NC_PRISM(N,K)"                : product of Newton-Cotes integration
                                         on prisms
   - "IM_GAUSS_PARALLELEPIPED(N,K)"    : product of Gauss1D integration
                                         on parallelepipeds (K=1,3,..99)
   - "IM_TRIANGLE(K)"                  : Gauss methods on triangles
                                         (K=1..10,13,17,19)
   - "IM_QUAD(K)"                      : Gauss methods on quadrilaterons
                                         (K=2, 3, 5, 7, 9, 17)
   - "IM_HEXAHEDRON(K)"                : Gauss methods on hexahedrons
                                         (K=5,9,11)
   - "IM_TETRAHEDRON(K)"               : Gauss methods on tetrahedrons
                                         (K=1, 2, 3, 5, 6, or 8)
   - "IM_SIMPLEX4D(3)"                 : Gauss method on a 4-dimensional
                                         simplex.
   - "IM_CUBE4D(K)"                    : Gauss method on a 4-dimensional cube
                                         (K=5,9).
   - "IM_STRUCTURED_COMPOSITE(IM1, K)" : Composite method on a grid with
                                         K divisions
   - "IM_HCT_COMPOSITE(IM1)"           : Composite integration suited to the
                                         HCT composite finite element.
   - "IM_QUAQC1_COMPOSITE(IM1)"        : Composite integration suited to the
                                         QUADC1 composite finite element.

   - "IM_QUASI_POLAR(IM1, IP1)"        : if IM1 is an integration method on a
               square, gives an integration method on a triangle which is
               close to a polar integration with respect to vertex IP1.
               if IM1 is an integration method on a tetrahedron, gives an
               integration method on a tetrahedron which is close to a
               cylindrical integration with respect to vertex IP1 (does not work very well).
               if IM1 is an integration method on a prism. Gives an integration
               method on a tetrahedron which is close to a
               cylindrical integration with respect to vertex IP1.
   - "IM_QUASI_POLAR(IM1, IP1, IP2)"   : IM1 should be an integration method
               on a prism. Gives an integration method on a tetrahedron which
               is close to a cylindrical integration with respect to IP1-IP2
               axis.
*/
#ifndef GETFEM_INTEGRATION_H__
#define GETFEM_INTEGRATION_H__

#include "getfem_config.h"
#include "bgeot_convex_ref.h"
#include "bgeot_geometric_trans.h"
#include "bgeot_node_tab.h"
#include "bgeot_poly_composite.h"
#include "getfem/dal_naming_system.h"




namespace getfem
{
  /** Description of an exact integration of polynomials.
   *  This class is not to be manipulate by itself. Use ppoly_integration and
   *  the functions written to produce the basic descriptions.
   */
  class poly_integration {
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
      long_scalar_type int_poly_on_face(const base_poly &P,
                                        short_type f) const;

      virtual ~poly_integration() {}

  };

  typedef std::shared_ptr<const poly_integration> ppoly_integration;

   /** Description of an approximate integration of polynomials of
   *  several variables on a reference element.
   *  This class is not to be manipulate by itself. Use papprox\_integration
   *  and the functions written to produce the basic descriptions.
   */

  class integration_method;  
  typedef std::shared_ptr<const integration_method> pintegration_method;


  class approx_integration {
  protected :
    
    typedef bgeot::node_tab PT_TAB;
    bgeot::pconvex_ref cvr;
    bgeot::pstored_point_tab pint_points;
    std::vector<scalar_type> int_coeffs;
    std::vector<size_type> repartition;
    
    // index 0 : points for volumic integration, index > 0 : points for faces
    std::vector<PT_TAB> pt_to_store; 
    bool valid;
    bool built_on_the_fly;
    
  public :
    
    /// Dimension of reference convex.
    dim_type dim(void) const { return cvr->structure()->dim(); }
    size_type nb_points(void) const { return int_coeffs.size(); }
    /// Number of integration nodes on the reference element.
    size_type nb_points_on_convex(void) const { return repartition[0]; }
    /// Number of integration nodes on the face f of the reference element.
    size_type nb_points_on_face(short_type f) const
    { return repartition[f+1] - repartition[f]; }
    size_type ind_first_point_on_face(short_type f) const 
    { return repartition[f]; }
    /// Convenience method returning the number of faces of the reference convex
    short_type nb_convex_faces() const
    { return cvr->structure()->nb_faces(); /* == repartition.size() - 1;*/ }
    bool is_built_on_the_fly(void) const { return built_on_the_fly; }
    void set_built_on_the_fly(void)  { built_on_the_fly = true; }
    /// Structure of the reference element.
    bgeot::pconvex_structure structure(void) const
    { return basic_structure(cvr->structure()); }
    bgeot::pconvex_ref ref_convex(void) const { return cvr; }
    
    const std::vector<size_type> &repart(void) const { return repartition; }
    
    /// Gives an array of integration nodes.
    // const bgeot::stored_point_tab &integration_points(void) const
    // { return *(pint_points); }
    bgeot::pstored_point_tab pintegration_points(void) const
    { return pint_points; }
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
    
    approx_integration(void) : valid(false), built_on_the_fly(false) { }
    approx_integration(bgeot::pconvex_ref cr)
      : cvr(cr), repartition(cr->structure()->nb_faces()+1),
        pt_to_store(cr->structure()->nb_faces()+1), valid(false),
        built_on_the_fly(false)
    { std::fill(repartition.begin(), repartition.end(), 0); }
    virtual ~approx_integration() {}
  };

  typedef std::shared_ptr<const approx_integration> papprox_integration;

  /**
     the list of main integration method types 
  */
  typedef enum { IM_APPROX, IM_EXACT, IM_NONE } integration_method_type;

  /**
     this structure is not intended to be used directly. It is built via
     the int_method_descriptor() function
  */
  class integration_method : virtual public dal::static_stored_object {
    ppoly_integration ppi; /* for exact integrations */
    papprox_integration pai; /* for approximate integrations (i.e. cubatures) */
    integration_method_type im_type;
    void remove() { pai =  papprox_integration(); ppi = ppoly_integration(); }

  public:
    integration_method_type type(void) const { return im_type; }
    const papprox_integration &approx_method(void) const { return pai; }
    const ppoly_integration &exact_method(void) const { return ppi; }

    void set_approx_method(papprox_integration paii)
    { remove(); pai = paii; im_type = IM_APPROX; }
    void set_exact_method(ppoly_integration ppii)
    { remove(); ppi = ppii; im_type = IM_EXACT; }

    bgeot::pstored_point_tab pintegration_points(void) const { 
      if (type() == IM_EXACT) {
        size_type n = ppi->structure()->dim();
        std::vector<base_node> spt(1); spt[0] = base_node(n);
        return store_point_tab(spt);
      }
      else if (type() == IM_APPROX)
        return pai->pintegration_points();
      else GMM_ASSERT1(false, "IM_NONE has no points");
    }

    // const bgeot::stored_point_tab &integration_points(void) const
    // { return *(pintegration_points()); }

    bgeot::pconvex_structure structure(void) const { 
      switch (type()) {
      case IM_EXACT: return ppi->structure();
      case IM_APPROX: return pai->structure();
      case IM_NONE: GMM_ASSERT1(false, "IM_NONE has no structure");
      default: GMM_ASSERT3(false, "");
      }
      return 0;
    }

    integration_method(ppoly_integration p) {
      DAL_STORED_OBJECT_DEBUG_CREATED(this, "Exact integration method");
      ppi = p; im_type = IM_EXACT;
    }

    integration_method(papprox_integration p) {
      DAL_STORED_OBJECT_DEBUG_CREATED(this, "Approximate integration method");
      pai = p; im_type = IM_APPROX;
    }

    integration_method() {
      DAL_STORED_OBJECT_DEBUG_CREATED(this, "Integration method");
      im_type = IM_NONE;
    }
    
    virtual ~integration_method()
    { DAL_STORED_OBJECT_DEBUG_DESTROYED(this, "Integration method"); }
  };


  /** Get an integration method from its name .
      @see @ref im_list 
      @param name the integration method name, for example "IM_TRIANGLE(6)"
      @param throw_if_not_found choose if an exception must be thrown
      when the integration method does not exist (if no exception, a
      null pointer is returned).
  */
  pintegration_method int_method_descriptor(std::string name,
                                            bool throw_if_not_found = true);

  /** Get the string name of an integration method .
   */
  std::string name_of_int_method(pintegration_method p);
  
  /**
     return an exact integration method for convex type handled by pgt.
     If pgt is not linear, classical_exact_im will fail.
  */
  pintegration_method classical_exact_im(bgeot::pgeometric_trans pgt);
  /**
     try to find an approximate integration method for the geometric
     transformation pgt which is able to integrate exactly polynomials
     of degree <= "degree". It may return a higher order integration
     method if no method match the exact degree.
  */
  pintegration_method classical_approx_im(bgeot::pgeometric_trans pgt, dim_type degree);

  /// return IM_EXACT_SIMPLEX(n)
  pintegration_method exact_simplex_im(size_type n);
  /// return IM_EXACT_PARALLELEPIPED(n)
  pintegration_method exact_parallelepiped_im(size_type n);
  /// return IM_EXACT_PRISM(n)
  pintegration_method exact_prism_im(size_type n);
  /// use classical_exact_im instead.
  pintegration_method exact_classical_im(bgeot::pgeometric_trans pgt) IS_DEPRECATED;
  /// return IM_NONE
  pintegration_method im_none(void);


  class mesh_im;
  papprox_integration composite_approx_int_method(const bgeot::mesh_precomposite &mp, 
                                                  const mesh_im &mf,
                                                  bgeot::pconvex_ref cr);

  /* try to integrate all monomials up to order 'order' and return the 
     max. error */
  scalar_type test_integration_error(papprox_integration pim, dim_type order);

  papprox_integration get_approx_im_or_fail(pintegration_method pim);

  /* Function allowing the add of an integration method outwards
     of getfem_integration.cc */
  
  typedef dal::naming_system<integration_method>::param_list im_param_list;

  void add_integration_name(std::string name,
                            dal::naming_system<integration_method>::pfunction f);

}  /* end of namespace getfem.                                            */


#endif /* BGEOT_INTEGRATION_H__                                           */
