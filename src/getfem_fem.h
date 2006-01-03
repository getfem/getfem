// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_fem.h : definition of the finite element methods
//           
// Date    : December 21, 1999.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1999-2005 Yves Renard
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

/**@file getfem_fem.h 
   @brief Definition of the finite element methods.

   This file defines the getfem::virtual_fem class, which is the common
   base class of all FEM.

   @section fem_list List of FEM known by getfem::fem_descriptor :
   
   - "FEM_PK(N,K)" : classical Lagrange element PK on a simplex.

   - "FEM_PK_DISCONTINUOUS(N,K,alpha)" : discontinuous Lagrange
   element PK on a simplex.

   - "FEM_QK(N,K)" : classical Lagrange element QK on a parellepiped.

   - "FEM_QK_DISCONTINUOUS(N,K,alpha)" : discontinuous Lagrange
   element QK on a parallelepiped.

   - "FEM_Q2_INCOMPLETE" : incomplete Q2 element with 8 dof (serendipity
                           Quad 8 element) 

   - "FEM_PK_PRISM(N,K)" : classical Lagrange element PK on a prism.

   - "FEM_PK_PRISM_DISCONTINUOUS(N,K,alpha)" : classical discontinuous
   Lagrange element PK on a prism.

   - "FEM_PK_WITH_CUBIC_BUBBLE(N,K)" : classical Lagrange element PK
   on a simplex with an additional volumic bubble function.

   - "FEM_PRODUCT(FEM1,FEM2)" : tensorial product of two polynomial
   elements

   - "FEM_P1_NONCONFORMING" : Nonconforming P1 method on a
   triangle.

   - "FEM_P1_BUBBLE_FACE(N)" : P1 method on a simplex with an
   additional bubble function on face 0.

   - "FEM_P1_BUBBLE_FACE_LAG" : P1 method on a simplex with an
   additional lagrange dof on face 0.

   - "FEM_HERMITE_SEGMENT" : Hermite element on the segment

   - "FEM_HERMITE_TRIANGLE" : Hermite element on the triangle

   - "FEM_PK_HIERARCHICAL(N,K)" : PK element with a hierarchical basis

   - "FEM_QK_HIERARCHICAL(N,K)" : QK element with a hierarchical basis

   - "FEM_PK_PRISM_HIERARCHICAL(N,K)" : PK element on a prism with a
   hierarchical basis

   - "FEM_STRUCTURED_COMPOSITE(FEM, K)" : Composite fem on a grid with
   K divisions

   - "FEM_PK_HIERARCHICAL_COMPOSITE(N,K,S)" : PK composite element on
   a grid with S subdivisions and with a hierarchical basis

   - "FEM_PK_FULL_HIERARCHICAL_COMPOSITE(N,K,S)" : PK composite
   element with S subdivisions and a hierarchical basis on both degree
   and subdivision

   
*/

#ifndef GETFEM_FEM_H__
#define GETFEM_FEM_H__

#include <dal_static_stored_objects.h>
#include <bgeot_geometric_trans.h>
#include <getfem_integration.h>
#include <getfem_poly_composite.h>
#include <deque>

namespace getfem {

  /** @defgroup dofdescr Dof description
   *  This set of functions gives a pointer to dof descriptions and
   *  tests the compatibility between two dof descriptions.
   *  A dof description describes a type of dof (Lagrange type, Hermite type ...)
   *  in order to be able to say wether two dofs are to be identified or not.
   *  The construction of dof type with the tensorial product is taken into account,
   *  making the dof_description structure a little bit complex (this structure will
   *  probably evoluate in the future).
   *  @{
   */

  struct dof_description;

  /// Type representing a pointer on a dof_description
  typedef dof_description *pdof_description;
  
  /** @brief Description of a unique dof of lagrange type (value at the node).
   * @param d the dimension of the reference element (2 for triangles, 3 for tetrahedrons ...)
   */
  pdof_description lagrange_dof(dim_type d);
  
  /** Description of a unique dof of derivative type. 
   *  (a derivative at the node).
   *  @param d the dimension of the reference element.
   *  @param r corresponds to the variable number for which the derivative is taken (0 <= r < d)
   */
  pdof_description derivative_dof(dim_type d, dim_type r);
  
  /** Description of a unique dof of normal derivative type
   *  (normal derivative at the node, regarding a face).
   *  @param d the dimension of the reference element.
   */  
  pdof_description norm_derivative_dof(dim_type d);

  /** Description of a unique dof of mean value type.
   *  @param d the dimension of the reference element.
   */
  pdof_description mean_value_dof(dim_type d);
  
  /** Description of a global dof, i.e. a numbered dof having a global scope.
   *  @param d the dimension of the reference element.
   */
  pdof_description global_dof(dim_type d);

  /// Product description of the descriptions *pnd1 and *pnd2.
  pdof_description product_dof(pdof_description pnd1, pdof_description pnd2);

  pdof_description to_coord_dof(pdof_description p, dim_type ct);

  /// Description of a special dof for Xfem
  pdof_description xfem_dof(pdof_description p, size_type ind);

  /// Returns the xfem_index of dof (0 for normal dof)
  size_type dof_xfem_index(pdof_description);

  size_type reserve_xfem_index(void);
  dim_type coord_index_of_dof(pdof_description a);

  /** Gives a total order on the dof description compatible with the
   *   identification.
   */
  int dof_description_compare(pdof_description a, pdof_description b);

  /// Says if the dof is linkable.
  bool dof_linkable(pdof_description);

  /// Says if the two dofs can be identified.
  bool dof_compatibility(pdof_description, pdof_description);


  /* @} */
  
  /* ******************************************************************** */
  /*	Classes for description of a finite element.                      */
  /* ******************************************************************** */
  
  /** @defgroup pfem Finite element description
   *
   * @{
   */


  class virtual_fem;

  /** type of pointer on a fem description @see getfem::virtual_fem */ 
  typedef boost::intrusive_ptr<const getfem::virtual_fem> pfem;
  
  struct fem_precomp_; 
  typedef boost::intrusive_ptr<const getfem::fem_precomp_> pfem_precomp;

  class fem_interpolation_context;
  
  /** @brief Base class for finite element description */
  class virtual_fem : virtual public dal::static_stored_object {

  protected :

    mutable std::vector<pdof_description> dof_types_;
    /* this convex structure is "owned" by the virtual_fem 
       (a deep copy is made when virtual_fems are copied) 
       But cvs_node has to be a pointer, as bgeot::convex_structure
       inherits from dal::static_stored_object
    */
    boost::intrusive_ptr<bgeot::convex_structure> cvs_node;
    bgeot::convex<base_node> cv_node;
    mutable bgeot::pstored_point_tab pspt;
    mutable bool pspt_valid;
    bgeot::pconvex_ref cvr; /// reference element.
    mutable dim_type ntarget_dim; /// dimension of the target space
    mutable dim_type dim_; /// dimension of the reference element 
    bool is_equiv; /// true if the FEM is equivalent
    bool is_lag; /// true if the FEM is of Lagrange type
    bool is_pol; /// true if the FEM is polynomial
    bool is_polycomp; /// true if the FEM is polynomial composite
    bool real_element_defined;
    short_type es_degree; /// estimated polynomial degree of the FEM
    short_type hier_raff; /// hierarchical refinement of the FEM
    
  public :
    /** Number of degrees of freedom. 

       @param cv the convex number for this FEM. This information is
       rarely used, by is needed by some "special" FEMs, such as
       getfem::interpolated_fem.
    */
    virtual size_type nb_dof(size_type /*cv*/) const
    { return dof_types_.size(); }
    /// Number of basis functions.
    virtual size_type nb_base(size_type cv) const
    { return nb_dof(cv); }
    /// Number of components (nb_dof() * dimension of the target space).
    size_type nb_base_components(size_type cv) const
      { return nb_base(cv) * ntarget_dim; }
    size_type nb_components(size_type cv) const
      { return nb_dof(cv) * ntarget_dim; }
    /// Get the array of pointer on dof description.
    const std::vector<pdof_description> &dof_types(void) const 
      { return dof_types_; }
    short_type hierarchical_raff(void) const { return hier_raff; }
    /// dimension of the reference element.
    dim_type dim(void) const { return dim_; }
    dim_type &dim(void) { return dim_; }
    /// dimension of the target space.
    dim_type target_dim(void) const { return ntarget_dim; }
    /// Return the convex of the reference element.
    virtual bgeot::pconvex_ref ref_convex(size_type) const { return cvr; }
    /// @internal
    bgeot::pconvex_ref &mref_convex() { return cvr; }
    /// Gives the convex of the reference element.
    bgeot::pconvex_structure basic_structure(size_type cv) const
    { return ref_convex(cv)->structure(); }
    /// Gives the convex representing the nodes on the reference element.
    virtual const bgeot::convex<base_node> &node_convex(size_type) const
      { return cv_node; }
    /// Gives the convex structure of the reference element nodes.
    bgeot::pconvex_structure structure(size_type cv) const
    { return node_convex(cv).structure(); }
    virtual bgeot::pstored_point_tab node_tab(size_type) const { 
      if (!pspt_valid) {
	pspt = bgeot::store_point_tab(cv_node.points());
	pspt_valid = true;
      }
      return pspt;
    }
    /** Gives the node corresponding to the dof i.
	@param cv the convex number for this FEM. This information is
	rarely used, by is needed by some "special" FEMs, such as
	getfem::interpolated_fem.
	@param i the local dof number (<tt>i < nb_dof(cv)</tt>)
    */
    const base_node &node_of_dof(size_type cv, size_type i) const
      { return (*(node_tab(cv)))[i];}
    bool is_on_real_element(void) const { return real_element_defined; }
    bool is_equivalent(void) const { return is_equiv; }
    bool need_G(void) const
    { return !(is_equivalent()) || real_element_defined; }
    /// true if the base functions are such that @f$ \varphi_i(\textrm{node\_of\_dof(j)}) = \delta_{ij} @f$
    bool is_lagrange(void) const { return is_lag; }
    /// true if the base functions are polynomials
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
    /** Interpolate at an arbitrary point x given on the reference
	element.
	
	@param c the fem_interpolation_context, should have been
	suitably initialized for the point of evaluation.

	@param coeff is the vector of coefficient relatively to
	the base functions, its length should be @c Qdim*this->nb_dof().

	@param val contains the interpolated value, on output (its
	size should be @c Qdim*this->target_dim()).
	
	@param Qdim is the optional Q dimension, if the FEM is
	considered as a "vectorized" one.
    */
    template<typename CVEC, typename VVEC> 
    void interpolation(const fem_interpolation_context& c, 
		       const CVEC& coeff, VVEC &val, dim_type Qdim) const;

    /** Build the interpolation matrix for the interpolation at a
	fixed point x, given on the reference element.

	The matrix @c M is filled, such that for a given @c coeff
	vector, the interpolation is given by @c M*coeff.
     */
    template <typename MAT>
    void interpolation(const fem_interpolation_context& c, 
		       MAT &M, dim_type Qdim) const;

    /** Interpolation of the gradient. The output is stored in the @f$
	Q\times N@f$ matrix @c val.
     */
    template<typename CVEC, typename VMAT> 
    void interpolation_grad(const fem_interpolation_context& c, 
			    const CVEC& coeff, VMAT &val,
			    dim_type Qdim=1) const;
      
    /** Give the value of all components of the base functions at the
     *  point x of the reference element. Basic function used essentially
     *  by fem_precomp.
     */
    virtual void base_value(const base_node &x, base_tensor &t) const = 0;

    /** Give the value of all gradients (on ref. element) of the components
     *  of the base functions at the point x of the reference element.
     *  Basic function used essentially by fem_precomp.
     */
    virtual void grad_base_value(const base_node &x, base_tensor &t) const = 0;

    /** Give the value of all hessians (on ref. element) of the components
     *  of the base functions at the point x of the reference element.
     *  Basic function used essentially by fem_precomp.
     */
    virtual void hess_base_value(const base_node &x, base_tensor &t) const = 0;

    /** Give the value of all components of the base functions at the
	current point of the fem_interpolation_context.  Used by
	elementary computations.  Matrix M for non tau-equivalent
	elements not taken into account.
    */
    virtual void real_base_value(const fem_interpolation_context &c, 
				 base_tensor &t) const;

    /** Give the gradient of all components of the base functions at the
	current point of the fem_interpolation_context.  Used by
	elementary computations.  Matrix M for non tau-equivalent
	elements not taken into account.
    */
    virtual void real_grad_base_value(const fem_interpolation_context &c, 
				      base_tensor &t) const;

    /** Give the hessian of all components of the base functions at the
	current point of the fem_interpolation_context.  Used by
	elementary computations.  Matrix M for non tau-equivalent
	elements not taken into account.
    */
    virtual void real_hess_base_value(const fem_interpolation_context &c, 
				      base_tensor &t) const;
    
    virtual size_type index_of_global_dof(size_type, size_type) const
      { DAL_THROW(internal_error, "internal error."); }

    /** internal function adding a node to an element for the creation
     * of a finite element method. Important : the faces should be the faces
     * on which the corresponding base function is non zero.
     */
    void add_node(const pdof_description &d, const base_node &pt,
		  const dal::bit_vector &faces);    
    void add_node(const pdof_description &d, const base_node &pt);
    void init_cvs_node(void);
    void unfreeze_cvs_node(void);

    virtual_fem &operator =(const virtual_fem &f) { 
      copy(f); return *this;
    }
    virtual_fem(void) { 
      ntarget_dim = 1; dim_ = 1; 
      is_equiv = is_pol = is_polycomp = is_lag = false;
      pspt_valid = false; hier_raff = 0; real_element_defined = false;
      es_degree = 5;
      cvs_node = bgeot::new_convex_structure();
    }
    virtual_fem(const virtual_fem& f) : dal::static_stored_object() { copy(f); }
    virtual ~virtual_fem() {}
  private:
    void copy(const virtual_fem &f) {
      dof_types_ = f.dof_types_;

      cvs_node = bgeot::new_convex_structure();
      *cvs_node = *f.cvs_node; // deep copy
      cv_node = f.cv_node;
      cv_node.structure() = cvs_node;
      pspt = 0;
      pspt_valid = false;
      cvr = f.cvr;
      dim_ = f.dim_;
      ntarget_dim = f.ntarget_dim;
      is_equiv = f.is_equiv;
      is_lag = f.is_lag;
      is_pol = f.is_pol;
      is_polycomp = f.is_polycomp;
      real_element_defined = f.real_element_defined;
      es_degree = f.es_degree;
      hier_raff = f.hier_raff;
    }
  };

  /**
     virtual_fem implementation as a vector of generic functions. The
     class FUNC should provide "derivative" and "eval" member
     functions (this is the case for bgeot::polynomial<T>).
  */
  template <class FUNC> class fem : public virtual_fem {
  protected :
    std::vector<FUNC> base_;
    
  public :
    
    /// Gives the array of basic functions (components).
    const std::vector<FUNC> &base(void) const { return base_; }
    std::vector<FUNC> &base(void) { return base_; }
    void base_value(const base_node &x, base_tensor &t) const {
      bgeot::multi_index mi(2);
      mi[1] = target_dim(); mi[0] = nb_base(0);
      t.adjust_sizes(mi);
      size_type R = nb_base_components(0);
      base_tensor::iterator it = t.begin();
      for (size_type  i = 0; i < R; ++i, ++it)
	*it = base_[i].eval(x.begin());
    }
    void grad_base_value(const base_node &x, base_tensor &t) const {
      bgeot::multi_index mi(3);
      dim_type n = dim();
      mi[2] = n; mi[1] = target_dim(); mi[0] = nb_base(0);
      t.adjust_sizes(mi);
      size_type R = nb_base_components(0);
      base_tensor::iterator it = t.begin();
      for (dim_type j = 0; j < n; ++j)
	for (size_type i = 0; i < R; ++i, ++it)
	  { FUNC f = base_[i]; f.derivative(j); *it = f.eval(x.begin()); }
    }
    void hess_base_value(const base_node &x, base_tensor &t) const {
      bgeot::multi_index mi(4);
      dim_type n = dim();
      mi[3] = n; mi[2] = n; mi[1] = target_dim(); mi[0] = nb_base(0);
      t.adjust_sizes(mi);
      size_type R = nb_base_components(0);
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
  
  /** Classical polynomial FEM. */
  typedef const fem<base_poly> * ppolyfem;
  /** Polynomial composite FEM */
  typedef const fem<polynomial_composite> * ppolycompfem;
  
  /** Give a pointer on the structures describing the classical
      polynomial fem of degree k on a given convex type.

      @param pgt the geometric transformation (which defines the convex type).
      @param k the degree of the fem.
      @return a ppolyfem.
  */
  pfem classical_fem(bgeot::pgeometric_trans pgt, short_type k);
  
  /** Give a pointer on the structures describing the classical
      polynomial discontinuous fem of degree k on a given convex type.

      @param pgt the geometric transformation (which defines the convex type).

      @param k the degree of the fem.

      @param alpha the "inset" factor for the dof nodes: with alpha =
      0, the nodes are located as usual (i.e. with node on the convex border), 
      and for 0 < alpha < 1, they converge to the center of gravity of the convex.

      @return a ppolyfem.
  */
  pfem classical_discontinuous_fem(bgeot::pgeometric_trans pg, short_type k, scalar_type alpha=0);

  /** get a fem descriptor from its string name. */
  pfem fem_descriptor(std::string name);

  /** get the string name of a fem descriptor. */
  std::string name_of_fem(pfem p);

  pfem PK_fem(size_type n, short_type k);
  pfem QK_fem(size_type n, short_type k);
  pfem PK_prism_fem(size_type n, short_type k);


  /**
     Pre-computations on a fem (given a fixed set of points on the
     reference convex, this object computes the value/gradient/hessian
     of all base functions on this set of points and stores them.
  */
  class fem_precomp_ : virtual public dal::static_stored_object {
  protected:      
    pfem pf;
    bgeot::pstored_point_tab pspt;
    mutable std::vector<base_tensor> c; /// store values of base functions
    mutable std::vector<base_tensor> pc; /// store gradients of base functions
    mutable std::vector<base_tensor> hpc; /// store hessians of base functions
  public:    
    inline const base_tensor &val(size_type i) const
    { if (c.empty()) init_val(); return c[i]; }
    inline const base_tensor &grad(size_type i) const
    { if (pc.empty()) init_grad(); return pc[i]; }
    inline const base_tensor &hess(size_type i) const
    { if (hpc.empty()) init_hess(); return hpc[i]; }
    inline pfem get_pfem() const { return pf; }
    inline const bgeot::stored_point_tab& get_point_tab() const
    { return *pspt; }
    fem_precomp_(pfem, bgeot::pstored_point_tab);
  private:
    void init_val() const;
    void init_grad() const;
    void init_hess() const;    
  };
  

  /** @brief Handles precomputations for FEM.  statically allocates a
     fem-precomputation object, and returns a pointer to it. The
     fem_precomp_ objects are "cached", i.e. they are stored in a
     global pool and if this function is called two times with the
     same arguments, a pointer to the same object will be returned.

     @param pf a pointer to the fem object.  @param pspt a pointer to
     a list of points in the reference convex.CAUTION: this array must
     not be destroyed as long as the fem_precomp is used!!.

     Moreover pspt is supposed to identify uniquely the set of
     points. This means that you should NOT alter its content at any
     time after using this function.

     If you need a set of "temporary" getfem::fem_precomp_, create
     them via a getfem::fem_precomp_pool structure. All memory will be
     freed when this structure will be destroyed.  */
  pfem_precomp fem_precomp(pfem pf, bgeot::pstored_point_tab pspt);

  /** Request for the removal of a pfem_precomp */
  inline void delete_fem_precomp(pfem_precomp pfp)
  { dal::del_stored_object(pfp); }


  /**
     handle a pool (i.e. a set) of fem_precomp. The difference with
     the global fem_precomp function is that these fem_precomp objects
     are freed when the fem_precomp_pool is destroyed (they can eat
     much memory). An example of use can be found in the
     getfem::interpolation_solution functions of getfem_export.h
  */
  class fem_precomp_pool {
    std::set<pfem_precomp> precomps;
    
  public :
    
    /** Request a pfem_precomp. If not already in the pool, the
	pfem_precomp is computed, and added to the pool.
	
        @param pf a pointer to the fem object.
	@param pspt a pointer to a list of points in the reference convex.CAUTION:
	this array must not be destroyed as long as the fem_precomp is used!!
	
	Moreover pspt is supposed to identify uniquely the set of
	points. This means that you should NOT alter its content until
	the fem_precomp_pool is destroyed.
    */
    pfem_precomp operator()(pfem pf, bgeot::pstored_point_tab pspt) {
      pfem_precomp p = fem_precomp(pf, pspt);
      precomps.insert(p);
      return p;
    }
    void clear(void) {
      for (std::set<pfem_precomp>::iterator it = precomps.begin();
	   it != precomps.end(); ++it)
	delete_fem_precomp(*it);
    }
    ~fem_precomp_pool() { clear(); }
  };
  
    
  /** structure passed as the argument of fem interpolation
      functions. This structure can be partially filled (for example
      the xreal will be computed if needed as long as pgp+ii is known).
  */
  class fem_interpolation_context :
    public bgeot::geotrans_interpolation_context {

    mutable base_matrix M_; /// optional transformation matrix (for non tau-equivalent fems)
    pfem pf_;               /// current fem 
    pfem_precomp pfp_;      /// optional fem_precomp_ (speed up the computations)
    size_type convex_num_;  /// the convex number (info needed by some specific FEMs)
  public:
    /// true if a fem_precomp_ has been supplied.
    bool have_pfp() const { return pfp_ != 0; }
    /// true if the pfem is available.
    bool have_pf() const { return pf_ != 0; }
    /// non tau-equivalent transformation matrix.
    const base_matrix& M() const;
    /** fill the tensor with the values of the base functions (taken
	at point @c this->xref())
    */
    void base_value(base_tensor& t) const;
    /** fill the tensor with the gradient of the base functions (taken
	at point @c this->xref())
    */
    void grad_base_value(base_tensor& t) const;
    /** fill the tensor with the hessian of the base functions (taken
	at point @c this->xref())
    */
    void hess_base_value(base_tensor& t) const;
    /** get the current FEM descriptor */
    const pfem pf() const { return pf_; }
    /** get the current convex number */
    size_type convex_num() const;
    /** get the current fem_precomp_ */
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
    if (gmm::vect_size(val) != Qdim)
      DAL_THROW(dimension_error, "dimensions mismatch");
    size_type R = nb_dof(c.convex_num()), RR = nb_base(c.convex_num());
    
    gmm::clear(val);
    const base_matrix *M = 0;
    if (!is_equivalent() && !is_on_real_element()) M = &(c.M());
    base_tensor Z; real_base_value(c,Z);
    for (size_type j = 0; j < RR; ++j) {
      for (size_type q = 0; q < Qmult; ++q) {
	typename gmm::linalg_traits<CVEC>::value_type co = 0.0;
	if (is_equivalent() || is_on_real_element())
	  co = coeff[j*Qmult+q];
	else
	  for (size_type i = 0; i < R; ++i)
	    co += coeff[i*Qmult+q] * (*M)(i, j);	  
	for (size_type r = 0; r < target_dim(); ++r)
	  val[r + q*target_dim()] += co * Z[j + r*R];
      } 
    }
  }

  template <typename MAT>
  void virtual_fem::interpolation(const fem_interpolation_context& c, 
				  MAT &M, dim_type Qdim) const {
    size_type Qmult = size_type(Qdim) / target_dim();
    size_type R = nb_dof(c.convex_num()), RR = nb_base(c.convex_num());
    if (gmm::mat_nrows(M) != Qdim || gmm::mat_ncols(M) != R*Qmult)
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    const base_matrix *MM = 0;
    if (!is_equivalent() && !is_on_real_element()) MM = &(c.M());
    gmm::clear(M);
    base_tensor Z; real_base_value(c,Z);
    for (size_type j = 0; j < RR; ++j) {
      for (size_type q = 0; q < Qmult; ++q) {
	for (size_type r = 0; r < target_dim(); ++r)
	  if (is_equivalent() || is_on_real_element())
	    M(r+q*target_dim(), j*Qmult+q) = Z[j + r*R];
	  else
	    for (size_type i = 0; i < R; ++i)
	      M(r+q*target_dim(), i*Qmult+q) = (*MM)(i, j) * Z[j + r*R];
      } 
    }
  }


  template<typename CVEC, typename VMAT> 
  void virtual_fem::interpolation_grad(const fem_interpolation_context& c, 
			  const CVEC& coeff, VMAT &val, dim_type Qdim) const {
    typedef typename gmm::linalg_traits<CVEC>::value_type T;
    size_type Qmult = size_type(Qdim) / target_dim();
    dim_type N = c.N();
    if (gmm::mat_ncols(val) != N || gmm::mat_nrows(val) != Qdim)
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    dim_type P = dim();
    base_tensor t;
    size_type R = nb_dof(c.convex_num()), RR = nb_base(c.convex_num());      
    const base_matrix *M = 0;
    if (!is_equivalent() && !is_on_real_element()) M = &(c.M());

    gmm::clear(val);
    if (!is_on_real_element()) { // optimized case
      if (!c.have_pfp()) {
	grad_base_value(c.xref(), t);
      } else {
	t = c.pfp()->grad(c.ii());
      }
      gmm::dense_matrix<T> val2(P, Qdim);
      gmm::clear(val2);
      for (size_type q = 0; q < Qmult; ++q) {
	base_tensor::const_iterator it = t.begin();	
	for (size_type k = 0; k < P; ++k)
	  for (size_type r = 0; r < target_dim(); ++r)
	    for (size_type j = 0; j < RR; ++j, ++it) {
	      T co = 0.0;
	      if (is_equivalent() || is_on_real_element())
		co = coeff[j*Qmult+q];
	      else
		for (size_type i = 0; i < R; ++i)
		  co += coeff[i*Qmult+q] * (*M)(i, j);	      
	      val2(k, r + q*target_dim()) += co * (*it);
	    }
      }
      gmm::mult(c.B(), val2, gmm::transposed(val));
      // replaced by the loop because gmm does not support product
      // of a complex<> matrix with a real matrix
      /*for (size_type i=0; i < Qdim; ++i)
	for (size_type j=0; j < N; ++j) {
	  T s = 0;
	  for (size_type k=0; k < P; ++k) s += c.B()(j,k)*val2(k,i);
	  val(i,j) = s;
	}
      */
    } else {
      real_grad_base_value(c, t);
      for (size_type q = 0; q < Qmult; ++q) {
	base_tensor::const_iterator it = t.begin();
	for (size_type k = 0; k < N; ++k)
	  for (size_type r = 0; r < target_dim(); ++r)
	    for (size_type j = 0; j < RR; ++j, ++it) {
	      T co = 0.0;
	      if (is_equivalent() || is_on_real_element())
		co = coeff[j*Qmult+q];
	      else
		for (size_type i = 0; i < R; ++i)
		  co += coeff[i*Qmult+q] * (*M)(i, j);	      
	      val(r + q*target_dim(), k) += co * (*it);
	    }
      }
    }
  }

  /* @} */

}  /* end of namespace getfem.                                            */


#endif
