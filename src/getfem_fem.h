/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library : GEneric Tool for Finite Element Methods (getfem)              */
/* File    : getfem_fem.h : definition of the finite element methods       */
/*                                                                         */
/* Date : December 21, 1999.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 1999-2002  Yves Renard.                                   */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */


#ifndef GETFEM_FEM_H__
#define GETFEM_FEM_H__

#include <bgeot_geometric_trans.h>
#include <bgeot_precomp.h>
#include <getfem_integration.h>
#include <getfem_poly_composite.h>
#include <getfem_precomp.h>
#include <deque>

namespace getfem
{
  /************************************************************************/
  /*	Class for description of an interpolation dof.                    */
  /************************************************************************/
  
  struct dof_description;
  /// Pointer on a dof_description
  typedef dof_description *pdof_description;
  
  /// Description of a unique dof of lagrange type (value at the node).
  pdof_description lagrange_dof(dim_type);
  /** Description of a unique dof of derivative type. 
   *  (a derivative at the node).
   */
  pdof_description derivative_dof(dim_type, dim_type);
  /** Description of a unique dof of normal derivative type
   *  (normal derivative at the node, regarding a face).
   */  
  pdof_description norm_derivative_dof(dim_type);
  /// Description of a unique dof of mean value type.
  pdof_description mean_value_dof(dim_type);
  
  pdof_description already_numerate_dof(dim_type);
  /// Product description of the descriptions *pnd1 and *pnd2.
  pdof_description product_dof(pdof_description, pdof_description);

  pdof_description to_coord_dof(pdof_description p, dim_type ct);

  pdof_description xfem_dof(pdof_description p, size_type ind);

  /** Gives a total order on the dof description compatible with the
   *   identification.
   */
  int dof_description_compare(pdof_description a, pdof_description b);
  /// Says if the dof is linkable.
  bool dof_linkable(pdof_description);
  /// Says if the two dofs can be identified.
  bool dof_compatibility(pdof_description, pdof_description);
  /// Returns the xfem_index of dof (0 for normal dof)
  size_type dof_xfem_index(pdof_description);
  
  /* ******************************************************************** */
  /*	Classes for description of a finite element.                      */
  /* ******************************************************************** */
  
  class virtual_fem;
  typedef const virtual_fem * pfem;

  class fem_interpolation_context;
  
  class virtual_fem
  {
  protected :

    std::vector<pdof_description> dof_types_;
    bgeot::convex_structure cvs_node;
    bgeot::convex<base_node> cv_node;
    mutable bgeot::pstored_point_tab pspt;
    mutable bool pspt_valid;
    bgeot::pconvex_ref cvr; // reference element.
    dim_type ntarget_dim;
    bool is_equiv, is_lag, is_pol, is_polycomp, real_element_defined;
    short_type es_degree, hier_raff;
    
  public :
    /// Number of degrees of freedom.
    virtual size_type nb_dof(void) const { return dof_types_.size(); }
    virtual size_type nb_base(void) const { return nb_dof(); }
    /// Number of components (nb_dof() * dimension of the target space).
    size_type nb_base_components(void) const
      { return nb_base() * ntarget_dim; }
    size_type nb_components(void) const
      { return nb_dof() * ntarget_dim; }
    /// Gives the array of pointer on dof description.
    const std::vector<pdof_description> &dof_types(void) const 
      { return dof_types_; }
    short_type hierarchical_raff(void) const { return hier_raff; }
    /// dimension of the reference element.
    dim_type dim(void) const { return cvr->structure()->dim(); }
    /// dimension of the target space.
    dim_type target_dim(void) const { return ntarget_dim; }
    /// Gives the convex structure of the reference element nodes.
    bgeot::pconvex_structure structure(void) const
      { return cv_node.structure(); }
    /// Gives the convex of the reference element.
    bgeot::pconvex_structure basic_structure(void) const
      { return cvr->structure(); }
    /// Gives the convex of the reference element.
    bgeot::pconvex_ref ref_convex(void) const { return cvr; }
    bgeot::pconvex_ref &ref_convex(void) { return cvr; }
    /// Gives the convex representing the nodes on the reference element.
    const bgeot::convex<base_node> &node_convex(void) const
      { return cv_node; }
    /// Gives the node corresponding to the dof i.
    const base_node &node_of_dof(size_type i) const
      { return cv_node.points()[i];}
    bgeot::pstored_point_tab node_tab(void) const { 
      if (!pspt_valid) {
	pspt = bgeot::store_point_tab(cv_node.points());
	pspt_valid = true;
      }
      return pspt;
    }
    bool is_on_real_element(void) const { return real_element_defined; }
    bool is_equivalent(void) const { return is_equiv; }
    bool need_G(void) const
    { return !(is_equivalent()) || real_element_defined; }
    bool is_lagrange(void) const { return is_lag; }
    bool is_polynomial(void) const { return is_pol; }
    bool is_polynomialcomp(void) const { return is_polycomp; }
    bool &is_polynomialcomp(void) { return is_polycomp; }
    bool &is_equivalent(void) { return is_equiv; }
    bool &is_lagrange(void) { return is_lag; }
    bool &is_polynomial(void) { return is_pol; }
    short_type estimated_degree(void) const { return es_degree; }
    short_type &estimated_degree(void) { return es_degree; }
    
    virtual void mat_trans(base_matrix &, const base_matrix &,
			   bgeot::pgeometric_trans) const
      { DAL_THROW(internal_error, "This function should not be called."); }
    /** Function which interpolates in a arbitrary point x given on the
     *  reference element. coeff is the vector of coefficient relatively to
     *  the shape functions. G and pgt represent the geometric transformation
     *  for non-equivalent elements.
     *  This method take all cases into account (non tau-equivalent element)
     */
    template<typename CVEC, typename VVEC> 
    void interpolation(const fem_interpolation_context& c, 
		       const CVEC& coeff, VVEC &val, dim_type Qdim=1) const;

    /** Function which interpolates in the ii th point of pfp. coeff is the
     *  vector of coefficient relatively to the shape functions. G and pgt
     *  represent the geometric transformation for non-equivalent elements.
     *  Qdim takes into account a vectorisation of the element.
     *  This method take all cases into account (non tau-equivalent element)
     */
    /** Function which interpolates the gradient in a arbitrary point x
     *  given on the reference element. coeff is the vector of coefficient
     *  relatively to the shape functions. G and pgt
     *  represent the geometric transformation for non-equivalent elements.
     *  Qdim take into account a vectorisation of the element.
     *  This method take all cases into account (non tau-equivalent element
     *  and correction of the gradient via the gradient of the geometric
     *  transformation).
     *  This function is essentially used by virtual_link_fem.
     */
    template<typename CVEC, typename VMAT> 
    void interpolation_grad(const fem_interpolation_context& c, 
			    const CVEC& coeff, VMAT &val, dim_type Qdim=1) const;
      
    /** Gives the value of all components of the base functions at the
     *  point x of the reference element. Basic function used essentially
     *  by fem_precomp.
     */
    virtual void base_value(const base_node &x, base_tensor &t) const = 0;
    /** Gives the value of all gradients (on ref. element) of the components
     *  of the base functions at the point x of the reference element.
     *  Basic function used essentially by fem_precomp.
     */
    virtual void grad_base_value(const base_node &x, base_tensor &t) const = 0;
    /** Gives the value of all hessians (on ref. element) of the components
     *  of the base functions at the point x of the reference element.
     *  Basic function used essentially by fem_precomp.
     */
    virtual void hess_base_value(const base_node &x, base_tensor &t) const = 0;
    /** Gives the value of all components of the base functions at the point
     *  ii of the pgp possibly using information on pfp and G.
     *  Used by elementary computations.
     *  Matrix M for non tau-equivalent elements not taken into account.
     */
    virtual void real_base_value(const fem_interpolation_context &c, 
				 base_tensor &t) const;

    /** Gives the value of all gradients on the real element of the components
     *  of the base functions at the point ii of the pgp possibly using
     *  information on pfp, G and B. B is the matrix wich transforms the
     *  gradient on the reference element to the gradient on the real element.
     *  Used by elementary computations.
     *  Matrix M for non tau-equivalent elements not taken into account.
     */
    virtual void real_grad_base_value(const fem_interpolation_context &c, 
				      base_tensor &t) const;
/** Gives the value of all hessians on the real element of the components
     *  of the base functions at the point ii of the pgp possibly using
     *  information on pfp, G, B2 and B32. B2 and B32 are the matrices used
     *  to transform a Hessian on reference element to a Hessian on real
     *  element. Used by elementary computations.
     *  Matrix M for non tau-equivalent elements not taken into account.
     */

    virtual void real_hess_base_value(const fem_interpolation_context &c, 
				      base_tensor &t) const;
    
    virtual size_type index_of_already_numerate_dof(size_type, size_type) const
      { DAL_THROW(internal_error, "internal error."); }

    virtual_fem(void) { 
      ntarget_dim = 1; is_equiv = is_pol = is_polycomp = is_lag = false;
      pspt_valid = false; hier_raff = 0; real_element_defined = false;
      es_degree = 5;
    }

    void add_node(const pdof_description &d, const base_node &pt) ;
    void init_cvs_node(void);
    void unfreeze_cvs_node(void);

    virtual_fem &operator =(const virtual_fem &f) {
      dof_types_ = f.dof_types_;
      cvs_node = f.cvs_node;
      cv_node = f.cv_node;
      cv_node.structure() = &(cvs_node);
      pspt = 0;
      pspt_valid = false;
      cvr = f.cvr;
      ntarget_dim = f.ntarget_dim;
      is_equiv = f.is_equiv;
      is_lag = f.is_lag;
      is_pol = f.is_pol;
      is_polycomp = f.is_polycomp;
      real_element_defined = f.real_element_defined;
      es_degree = f.es_degree;
      hier_raff = f.hier_raff;
      return *this;
    }
    virtual_fem(const virtual_fem &f) { *this = f; }

    virtual ~virtual_fem() {}
  };
  
  /* the fem_interpolation_context structure is passed as the argument
     of fem interpolation functions. This structure can be partially
     filled (for example the xreal will be computed if needed as long
     as pgp+ii is known)
  */
  class fem_interpolation_context : public bgeot::geotrans_interpolation_context {
    mutable base_matrix M_;
    pfem pf_;
    pfem_precomp pfp_;
    size_type convex_num_;
  public:
    bool have_pfp() const { return pfp_ != 0; }
    bool have_pf() const { return pf_ != 0; }
    const base_matrix& M() const;
    void base_value(base_tensor& t) const;
    void grad_base_value(base_tensor& t) const;
    void hess_base_value(base_tensor& t) const;
    const pfem pf() const { return pf_; }
    size_type convex_num() const;
    pfem_precomp pfp() const { return pfp_; }
    void set_pfp(pfem_precomp newpfp);
    void set_pf(pfem newpf);
    fem_interpolation_context();
    fem_interpolation_context(bgeot::pgeotrans_precomp pgp__, 
			      pfem_precomp pfp__, size_type ii__, 
			      const base_matrix& G__, 
			      size_type convex_num__);
    fem_interpolation_context(bgeot::pgeometric_trans pgt__, 
			      pfem_precomp pfp__, size_type ii__, 
			      const base_matrix& G__, 
			      size_type convex_num__);
    fem_interpolation_context(bgeot::pgeometric_trans pgt__,
			      pfem pf__,
			      const base_node& xref__,
			      const base_matrix& G__,
			      size_type convex_num__);
  };

  template <typename CVEC, typename VVEC>
  void virtual_fem::interpolation(const fem_interpolation_context& c, 
				  const CVEC& coeff, VVEC &val,
				  dim_type Qdim) const {
    size_type Qmult = size_type(Qdim) / target_dim();
    if (val.size() != Qdim)
      DAL_THROW(dimension_error, "dimensions mismatch");
    size_type R = nb_dof(), RR = nb_base();
    
    gmm::clear(val);
    base_tensor Z; real_base_value(c,Z);
    for (size_type j = 0; j < RR; ++j) {
      for (size_type q = 0; q < Qmult; ++q) {
	scalar_type co = 0.0;
	if (is_equivalent())
	  co = coeff[j*Qmult+q];
	else
	  for (size_type i = 0; i < R; ++i)
	    co += coeff[i*Qmult+q] * c.M()(i, j);	  
	for (size_type r = 0; r < target_dim(); ++r)
	  val[r + q*target_dim()] += co * Z[j + r*R];
      } 
    }
  }

  template<typename CVEC, typename VMAT> 
  void virtual_fem::interpolation_grad(const fem_interpolation_context& c, 
			  const CVEC& coeff, VMAT &val, dim_type Qdim) const {
    size_type Qmult = size_type(Qdim) / target_dim();
    dim_type N = c.N();
    if (gmm::mat_nrows(val) != N || gmm::mat_ncols(val) != Qdim)
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    dim_type P = dim();
    base_tensor t;
    size_type R = nb_dof(), RR = nb_base();      
    
    gmm::clear(val);
    if (!is_on_real_element()) { // optimized case
      grad_base_value(c.xref(), t);
      base_matrix val2(P, Qdim);
      gmm::clear(val2);
      for (size_type q = 0; q < Qmult; ++q) {
	base_tensor::iterator it = t.begin();	
	for (size_type k = 0; k < P; ++k)
	  for (size_type r = 0; r < target_dim(); ++r)
	    for (size_type j = 0; j < RR; ++j, ++it) {
	      scalar_type co = 0.0;
	      if (is_equivalent())
		co = coeff[j*Qmult+q];
	      else
		for (size_type i = 0; i < R; ++i)
		  co += coeff[i*Qmult+q] * c.M()(i, j);	      
	      val2(k, r + q*target_dim()) += co * (*it);
	    }
      }
      gmm::mult(c.B(), val2, val);
    } else {
      real_grad_base_value(c, t);
      for (size_type q = 0; q < Qmult; ++q) {
	base_tensor::iterator it = t.begin();
	for (size_type k = 0; k < N; ++k)
	  for (size_type r = 0; r < target_dim(); ++r)
	    for (size_type j = 0; j < RR; ++j, ++it) {
	      scalar_type co = 0.0;
	      if (is_equivalent())
		co = coeff[j*Qmult+q];
	      else
		for (size_type i = 0; i < R; ++i)
		  co += coeff[i*Qmult+q] * c.M()(i, j);	      
	      val(k, r + q*target_dim()) += co * (*it);
	    }
      }
    }
  }

  
  template <class FUNC> class fem : public virtual_fem {
  protected :
    std::vector<FUNC> base_;
    
  public :
    
    /// Gives the array of basic functions (components).
    const std::vector<FUNC> &base(void) const { return base_; }
    std::vector<FUNC> &base(void) { return base_; }
    void base_value(const base_node &x, base_tensor &t) const {
      bgeot::multi_index mi(2);
      mi[1] = target_dim(); mi[0] = nb_base();
      t.adjust_sizes(mi);
      size_type R = nb_base_components();
      base_tensor::iterator it = t.begin();
      for (size_type  i = 0; i < R; ++i, ++it)
	*it = base_[i].eval(x.begin());
    }
    void grad_base_value(const base_node &x, base_tensor &t) const {
      bgeot::multi_index mi(3);
      dim_type n = dim();
      mi[2] = n; mi[1] = target_dim(); mi[0] = nb_base();
      t.adjust_sizes(mi);
      size_type R = nb_base_components();
      base_tensor::iterator it = t.begin();
      for (dim_type j = 0; j < n; ++j)
	for (size_type i = 0; i < R; ++i, ++it)
	  { FUNC f = base_[i]; f.derivative(j); *it = f.eval(x.begin()); }
    }
    void hess_base_value(const base_node &x, base_tensor &t) const {
      bgeot::multi_index mi(4);
      dim_type n = dim();
      mi[3] = n; mi[2] = n; mi[1] = target_dim(); mi[0] = nb_base();
      t.adjust_sizes(mi);
      size_type R = nb_base_components();
      base_tensor::iterator it = t.begin();
      for (dim_type k = 0; k < n; ++k)
	for (dim_type j = 0; j < n; ++j)
	  for (size_type i = 0; i < R; ++i, ++it) {
	    FUNC f = base_[i];
	    f.derivative(j); f.derivative(k);
	    *it = f.eval(x.begin());
	  }
    }
    
    
  };
  
  typedef const fem<base_poly> * ppolyfem;
  typedef const fem<polynomial_composite> * ppolycompfem;
  
  /** Gives a pointer on the structures describing the more classical fem
   *  of degree k on a geometric convex cvs (coming from the geometric trans).
   */
  pfem classical_fem(bgeot::pgeometric_trans pg, short_type k);
  /** Gives a pointer on the structures describing the more classical discontinuous fem
   *  of degree k on a geometric convex cvs (coming from the geometric trans).
   */
  pfem classical_discontinuous_fem(bgeot::pgeometric_trans pg, short_type k);

  pfem fem_descriptor(std::string name);
  /*  List of elements :
   *  "FEM_PK(N,K)"                      : classical Lagrange element PK on a
   *                                     :  simplex
   *  "FEM_PK_DISCONTINUOUS(N,K)"        : discontinuous Lagrange element PK
   *                                     :  on a simplex
   *  "FEM_QK(N,K)"                      : classical Lagrange element QK on a
   *                                     :  parellepiped
   *  "FEM_PK_PRISM(N,K)"                : classical Lagrange element PK on a
   *                                     :  prism
   *  "FEM_PK_WITH_CUBIC_BUBBLE(N,K)"    : classical Lagrange element PK on a
   *                                     :  simplex with an additional volumic
   *                                     :  bubble function.
   *  "FEM_PRODUCT(FEM1,FEM2)"           : tensorial product of two polynomial
   *                                     :  elements
   *  "FEM_P1_NONCONFORMING"             : Nonconforming P1 method on a
   *                                     :  triangle.
   *  "FEM_P1_BUBBLE_FACE(N)"            : P1 method on a simplex with an
   *                                     :  additional bubble function on
   *                                     :  face 0.
   *  "FEM_P1_BUBBLE_FACE_LAG"           : P1 method on a simplex with an
   *                                     :  additional lagrange dof on face 0.
   *  "FEM_HERMITE_SEGMENT"              : Hermite element on the segment
   *  "FEM_PK_HIERARCHICAL(N,K)"         : PK element with a hierarchical basis
   *  "FEM_QK_HIERARCHICAL(N,K)"         : QK element with a hierarchical basis
   *  "FEM_PK_PRISM_HIERARCHICAL(N,K)"   : PK element on a prism with a
   *                                     :  hierarchical basis
   *  "FEM_STRUCTURED_COMPOSITE(FEM, K)" : Composite fem on a grid with
   *                                     :  K divisions
   *  "FEM_PK_HIERARCHICAL_COMPOSITE(N,K,S)" : PK composite element on a grid
   *                                     :  with S subdivisions and with a
   *                                     : hierarchical basis
   *  "FEM_PK_FULL_HIERARCHICAL_COMPOSITE(N,K,S)" : PK composite element with
   *                                     :  S subdivisions and a hierarchical
   *                                     :  basis on both degree and
   *                                     :  subdivision
   */

  pfem PK_fem(size_type n, short_type k);
  pfem QK_fem(size_type n, short_type k);
  pfem PK_prism_fem(size_type n, short_type k);

  std::string name_of_fem(pfem p);

  
  class mesh_fem;
  pfem virtual_link_fem(mesh_fem &mf1, mesh_fem &mf2,
			pintegration_method pim);
  pfem virtual_link_fem_with_gradient(mesh_fem &mf1, mesh_fem &mf2,
				      pintegration_method pim);
  
}  /* end of namespace getfem.                                            */


#endif
