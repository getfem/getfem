/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 1999-2024 Yves Renard

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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

/**@file getfem_fem.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date December 21, 1999.
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

   - "FEM_Q2_INCOMPLETE(N)" : incomplete Q2 elements with 8 and 20 dof
   (serendipity Quad 8 and Hexa 20 elements)

   - "FEM_Q2_INCOMPLETE_DISCONTINUOUS(N)" : discontinuous incomplete Q2
   elements with 8 and 20 dof (serendipity Quad 8 and Hexa 20 elements)

   - "FEM_PRISM_PK(N,K)" : classical Lagrange element PK on a prism.

   - "FEM_PRISM_PK_DISCONTINUOUS(N,K,alpha)" : classical discontinuous
   Lagrange element PK on a prism.

   - "FEM_PRISM_INCOMPLETE_P2" : Incomplete Lagrange element on a
   quadratic 3D prism (serendipity, 15-node wedge element). Can be connected
   toa standard P2 Lagrange on its triangular faces and a Q2_INCOMPLETE
   Lagrange element on its quadrangular faces.

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

   - "FEM_HERMITE(N)" : Hermite element P3 on dimension N (1, 2 por 3).

   - "FEM_ARGYRIS" : Argyris element on the triangle.

   - "FEM_HCT_TRIANGLE" : Hsieh-Clough-Tocher element on the triangle
       (composite P3 element which is C^1, 12 dof).

   - "FEM_REDUCED_HCT_TRIANGLE" : Hsieh-Clough-Tocher element on the triangle
       (composite P3 element which is C^1, 9 dof).

   - "FEM_QUADC1_COMPOSITE" : quadrilateral element, composite P3 element
      and C^1 (16 dof).

   - "FEM_REDUCED_QUADC1_COMPOSITE" : quadrilateral element, composite
      P3 element and C^1 (12 dof).

   - "FEM_PK_HIERARCHICAL(N,K)" : PK element with a hierarchical basis.

   - "FEM_QK_HIERARCHICAL(N,K)" : QK element with a hierarchical basis.

   - "FEM_PRISM_PK_HIERARCHICAL(N,K)" : PK element on a prism with a
   hierarchical basis.

   - "FEM_STRUCTURED_COMPOSITE(FEM, K)" : Composite fem on a grid with
   K divisions.

   - "FEM_PK_HIERARCHICAL_COMPOSITE(N,K,S)" : PK composite element on
   a grid with S subdivisions and with a hierarchical basis.

   - "FEM_PK_FULL_HIERARCHICAL_COMPOSITE(N,K,S)" : PK composite
   element with S subdivisions and a hierarchical basis on both degree
   and subdivision.

   - "FEM_PYRAMID_QK(K)" : Lagrange element on a 3D pyramid of degree
   K=0, 1 or 2. Can be connected to a standard P1/P2 Lagrange element on its
   triangular faces and a standard Q1/Q2 Lagrange element on its quadrangular
   face.

   - "FEM_PYRAMID_QK_DISCONTINUOUS(K)" : Discontinuous Lagrange element
   on a 3D pyramid of degree K = 0, 1 or 2.

   - "FEM_PYRAMID_Q2_INCOMPLETE" : Incomplete Lagrange element on a
   quadratic 3D pyramid (serendipity, 13-node element). Can be connected to
   a standard P2 Lagrange element on its triangular faces and a Q2_INCOMPLETE
   Lagrange element on its quadrangular face.

   - "HHO(fem_interior, fem_face_1, ..., fem_face_n)" : Build a hybrid method
     with "fem_interior" on the element itself and "fem_face_1", ...,
     "fem_face_n" on each face. If only one method is given for the faces, it
     is duplicated on each face.

*/

#ifndef GETFEM_FEM_H__
#define GETFEM_FEM_H__

#include "dal_static_stored_objects.h"
#include "bgeot_geometric_trans.h"
#include "bgeot_poly_composite.h"
#include "getfem_integration.h"
#include "dal_naming_system.h"
#include <deque>

namespace getfem {

  /** @defgroup dofdescr Dof description
   *  This set of functions gives a pointer to dof descriptions and
   *  tests the compatibility between two dof descriptions.
   *  A dof description describes a type of dof (Lagrange type,
   *  Hermite type ...) in order to be able to say wether two dofs are
   *  to be identified or not. The construction of dof type with
   *  the tensorial product is taken into
   *  account,  making the dof_description structure a little bit complex
   *  (this structure will probably evoluate in the future).
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

  /** Description of a unique dof of second derivative type.
   *  @param d the dimension of the reference element.
   *  @param num_der1 corresponds to the variable number for which the first derivative is taken (0 <= r < d)
   *  @param num_der2 corresponds to the variable number for which the second derivative is taken (0 <= r < d)
   */
  pdof_description second_derivative_dof(dim_type d,
                                         dim_type num_der1,
                                         dim_type num_der2);

  /** Description of a unique dof of normal derivative type
   *  (normal derivative at the node, regarding a face).
   *  @param d the dimension of the reference element.
   */
  pdof_description normal_derivative_dof(dim_type d);

  /** Description of a unique dof of mean value type.
   *  @param d the dimension of the reference element.
   */
  pdof_description mean_value_dof(dim_type d);

  pdof_description bubble1_dof(dim_type ct);

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

  size_type reserve_xfem_index();
  dim_type coord_index_of_dof(pdof_description a);

  int dof_weak_compatibility(pdof_description a, pdof_description b);

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
  /*    Classes for description of a finite element.                      */
  /* ******************************************************************** */

  /** @defgroup pfem Finite element description
   *
   * @{
   */


  class virtual_fem;

  /** type of pointer on a fem description @see getfem::virtual_fem */
  typedef std::shared_ptr<const getfem::virtual_fem> pfem;

  class fem_precomp_;
  typedef std::shared_ptr<const getfem::fem_precomp_> pfem_precomp;

  class fem_interpolation_context;

  /** @brief Base class for finite element description */
  class virtual_fem : virtual public dal::static_stored_object,
                      public std::enable_shared_from_this<const virtual_fem> {
  public :
    enum vec_type { VECTORIAL_NOTRANSFORM_TYPE, VECTORIAL_PRIMAL_TYPE,
                    VECTORIAL_DUAL_TYPE };

  protected :

    mutable std::vector<pdof_description> dof_types_;
    /* this convex structure is "owned" by the virtual_fem
       (a deep copy is made when virtual_fems are copied)
       But cvs_node has to be a pointer, as bgeot::convex_structure
       inherits from dal::static_stored_object
    */
    std::shared_ptr<bgeot::convex_structure> cvs_node;
    std::vector<std::vector<short_type>> face_tab; // face list for each dof
    bgeot::convex<base_node> cv_node;
    mutable bgeot::pstored_point_tab pspt;
    mutable bool pspt_valid;
    bgeot::pconvex_ref cvr; // reference element.
    dim_type ntarget_dim;   // dimension of the target space
    mutable dim_type dim_;  // dimension of the reference element
    bool is_equiv;          // true if the FEM is equivalent
    bool is_lag;            // true if the FEM is of Lagrange type
    bool is_pol;            // true if the FEM is polynomial
    bool is_polycomp;       // true if the FEM is polynomial composite
    bool real_element_defined;
    bool is_standard_fem;   // is_equiv && !real_element_defined && scalar
    short_type es_degree;   // estimated polynomial degree of the FEM
    short_type hier_raff;   // hierarchical refinement of the FEM
    vec_type vtype; // for vectorial elements, type of transformation
                    // from the reference element.
    std::string debug_name_;

  public :
    /** Number of degrees of freedom.

       @param cv the convex number for this FEM. This information is
       rarely used, but is needed by some "special" FEMs, such as
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
    const std::vector<pdof_description> &dof_types() const
      { return dof_types_; }
    short_type hierarchical_raff() const { return hier_raff; }
    /// dimension of the reference element.
    dim_type dim() const { return dim_; }
    dim_type &dim() { return dim_; }
    /// dimension of the target space.
    dim_type target_dim() const { return ntarget_dim; }
    /// Type of vectorial element.
    vec_type vectorial_type() const { return vtype; }
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
    const std::string &debug_name() const { return debug_name_; }
    std::string &debug_name() { return debug_name_; }
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
    virtual const std::vector<short_type> &
    faces_of_dof(size_type /*cv*/, size_type i) const;
    bool is_on_real_element() const { return real_element_defined; }
    bool is_equivalent() const { return is_equiv; }
    bool need_G() const
    { return !(is_equivalent()) || real_element_defined; }
    /// true if the base functions are such that @f$ \varphi_i(\textrm{node\_of\_dof(j)}) = \delta_{ij} @f$
    bool is_lagrange() const { return is_lag; }
    /// true if the base functions are polynomials
    bool is_polynomial() const { return is_pol; }
    bool is_polynomialcomp() const { return is_polycomp; }
    bool is_standard() const { return is_standard_fem; }
    bool &is_polynomialcomp() { return is_polycomp; }
    bool &is_equivalent() { return is_equiv; }
    bool &is_lagrange() { return is_lag; }
    bool &is_polynomial() { return is_pol; }
    bool &is_standard() { return is_standard_fem; }
    short_type estimated_degree() const { return es_degree; }
    short_type &estimated_degree() { return es_degree; }

    virtual void mat_trans(base_matrix &, const base_matrix &,
                           bgeot::pgeometric_trans) const
    { GMM_ASSERT1(false, "This function should not be called."); }
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

    /** Interpolation of the hessian. The output is stored in the @f$
        Q\times (N^2)@f$ matrix @c val.
     */
    template<typename CVEC, typename VMAT>
    void interpolation_hess(const fem_interpolation_context& c,
                            const CVEC& coeff, VMAT &val,
                            dim_type Qdim) const;

    /** Interpolation of the divergence. The output is stored in the
        scalar @c val.
     */
    template<typename CVEC>
    void interpolation_diverg
      (const fem_interpolation_context& c, const CVEC& coeff,
       typename gmm::linalg_traits<CVEC>::value_type &val) const;

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
        elementary computations.  if withM is false the matrix M for
        non tau-equivalent elements is not taken into account.
    */
    virtual void real_base_value(const fem_interpolation_context &c,
                                 base_tensor &t, bool withM = true) const;

    /** Give the gradient of all components of the base functions at the
        current point of the fem_interpolation_context.  Used by
        elementary computations.  if withM is false the matrix M for
        non tau-equivalent elements is not taken into account.
    */
    virtual void real_grad_base_value(const fem_interpolation_context &c,
                                      base_tensor &t, bool withM = true) const;

    /** Give the hessian of all components of the base functions at the
        current point of the fem_interpolation_context.  Used by
        elementary computations. if withM is false the matrix M for
        non tau-equivalent elements is not taken into account.
    */
    virtual void real_hess_base_value(const fem_interpolation_context &c,
                                      base_tensor &t, bool withM = true) const;

    virtual size_type index_of_global_dof(size_type, size_type) const
      { GMM_ASSERT1(false, "internal error."); }

    /** internal function adding a node to an element for the creation
     * of a finite element method. Important : the faces should be the faces
     * on which the corresponding base function is non zero.
     */
    void add_node(const pdof_description &d, const base_node &pt,
                  const dal::bit_vector &faces);
    void add_node(const pdof_description &d, const base_node &pt);
    void init_cvs_node();
    void unfreeze_cvs_node();

    virtual_fem &operator =(const virtual_fem &f) {
      copy(f); return *this;
    }

    virtual_fem() {
      DAL_STORED_OBJECT_DEBUG_CREATED(this, "Fem");
      ntarget_dim = 1; dim_ = 1;
      is_equiv = is_pol = is_polycomp = is_lag = is_standard_fem = false;
      pspt_valid = false; hier_raff = 0; real_element_defined = false;
      es_degree = 5;
      vtype = VECTORIAL_NOTRANSFORM_TYPE;
      cvs_node = bgeot::new_convex_structure();
    }
    virtual_fem(const virtual_fem& f)
      : dal::static_stored_object(),
        std::enable_shared_from_this<const virtual_fem>()
    { copy(f); DAL_STORED_OBJECT_DEBUG_CREATED(this, "Fem"); }
    virtual ~virtual_fem() { DAL_STORED_OBJECT_DEBUG_DESTROYED(this, "Fem"); }
  private:
    void copy(const virtual_fem &f);
  };

  /**
     virtual_fem implementation as a vector of generic functions. The
     class FUNC should provide "derivative" and "eval" member
     functions (this is the case for bgeot::polynomial<T>).
  */
  template <class FUNC> class fem : public virtual_fem {
  protected :
    std::vector<FUNC> base_;
    mutable std::vector<std::vector<FUNC>> grad_, hess_;
    mutable bool grad_computed_ = false;
    mutable bool hess_computed_ = false;

    void compute_grad_() const {
      if (grad_computed_) return;
      GLOBAL_OMP_GUARD
      if (grad_computed_) return;
      size_type R = nb_base_components(0);
      dim_type n = dim();
      grad_.resize(R);
      for (size_type i = 0; i < R; ++i) {
        grad_[i].resize(n);
        for (dim_type j = 0; j < n; ++j) {
          grad_[i][j] = base_[i]; grad_[i][j].derivative(j);
        }
      }
      grad_computed_ = true;
    }

    void compute_hess_() const {
      if (hess_computed_) return;
      GLOBAL_OMP_GUARD
      if (hess_computed_) return;
      size_type R = nb_base_components(0);
      dim_type n = dim();
      hess_.resize(R);
      for (size_type i = 0; i < R; ++i) {
        hess_[i].resize(n*n);
        for (dim_type j = 0; j < n; ++j) {
          for (dim_type k = 0; k < n; ++k) {
            hess_[i][j+k*n] = base_[i];
            hess_[i][j+k*n].derivative(j); hess_[i][j+k*n].derivative(k);
          }
        }
      }
      hess_computed_ = true;
    }

  public :

    /// Gives the array of basic functions (components).
    const std::vector<FUNC> &base() const { return base_; }
    std::vector<FUNC> &base() { return base_; }
    /** Evaluates at point x, all base functions and returns the result in
        t(nb_base,target_dim) */
    void base_value(const base_node &x, base_tensor &t) const {
      bgeot::multi_index mi(2);
      mi[1] = target_dim(); mi[0] = short_type(nb_base(0));
      t.adjust_sizes(mi);
      size_type R = nb_base_components(0);
      base_tensor::iterator it = t.begin();
      for (size_type  i = 0; i < R; ++i, ++it)
        *it = bgeot::to_scalar(base_[i].eval(x.begin()));
    }
    /** Evaluates at point x, the gradient of all base functions w.r.t. the
        reference element directions 0,..,dim-1 and returns the result in
        t(nb_base,target_dim,dim) */
    void grad_base_value(const base_node &x, base_tensor &t) const {
      if (!grad_computed_) compute_grad_();
      bgeot::multi_index mi(3);
      dim_type n = dim();
      mi[2] = n; mi[1] = target_dim(); mi[0] = short_type(nb_base(0));
      t.adjust_sizes(mi);
      size_type R = nb_base_components(0);
      base_tensor::iterator it = t.begin();
      for (dim_type j = 0; j < n; ++j)
        for (size_type i = 0; i < R; ++i, ++it)
	  *it = bgeot::to_scalar(grad_[i][j].eval(x.begin()));
    }
    /** Evaluates at point x, the hessian of all base functions w.r.t. the
        reference element directions 0,..,dim-1 and returns the result in
        t(nb_base,target_dim,dim,dim) */
    void hess_base_value(const base_node &x, base_tensor &t) const {
      if (!hess_computed_) compute_hess_();
      bgeot::multi_index mi(4);
      dim_type n = dim();
      mi[3] = n; mi[2] = n; mi[1] = target_dim();
      mi[0] = short_type(nb_base(0));
      t.adjust_sizes(mi);
      size_type R = nb_base_components(0);
      base_tensor::iterator it = t.begin();
      for (dim_type k = 0; k < n; ++k)
        for (dim_type j = 0; j < n; ++j)
          for (size_type i = 0; i < R; ++i, ++it)
	    *it = bgeot::to_scalar(hess_[i][j+k*n].eval(x.begin()));
    }

  };

  /** Classical polynomial FEM. */
  typedef const fem<bgeot::base_poly> * ppolyfem;
  /** Polynomial composite FEM */
  typedef const fem<bgeot::polynomial_composite> * ppolycompfem;
  /** Rational fration FEM */
  typedef const fem<bgeot::base_rational_fraction> * prationalfracfem;

  /** Give a pointer on the structures describing the classical
      polynomial fem of degree k on a given convex type.

      @param pgt the geometric transformation (which defines the convex type).
      @param k the degree of the fem.
      @param complete a flag which requests complete Langrange polynomial
      elements even if the provided pgt is an incomplete one (e.g. 8-node
      quadrilateral or 20-node hexahedral).
      @return a ppolyfem.
  */
  pfem classical_fem(bgeot::pgeometric_trans pgt, short_type k,
                     bool complete=false);

  /** Give a pointer on the structures describing the classical
      polynomial discontinuous fem of degree k on a given convex type.

      @param pgt the geometric transformation (which defines the convex type).

      @param k the degree of the fem.

      @param alpha the "inset" factor for the dof nodes: with alpha =
      0, the nodes are located as usual (i.e. with node on the convex border),
      and for 0 < alpha < 1, they converge to the center of gravity of the
      convex.

      @param complete a flag which requests complete Langrange polynomial
      elements even if the provided pgt is an incomplete one (e.g. 8-node
      quadrilateral or 20-node hexahedral).

      @return a ppolyfem.
  */
  pfem classical_discontinuous_fem(bgeot::pgeometric_trans pg, short_type k,
                                   scalar_type alpha=0, bool complete=false);

  /** get a fem descriptor from its string name. */
  pfem fem_descriptor(const std::string &name);

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
    const pfem pf;
    const bgeot::pstored_point_tab pspt;
    mutable std::vector<base_tensor> c;   // stored values of base functions
    mutable std::vector<base_tensor> pc;  // stored gradients of base functions
    mutable std::vector<base_tensor> hpc; // stored hessians of base functions
  public:
    /// returns values of the base functions
    inline const base_tensor &val(size_type i) const
      { if (c.empty()) init_val(); return c[i]; }
    /// returns gradients of the base functions
    inline const base_tensor &grad(size_type i) const
      { if (pc.empty()) init_grad(); return pc[i]; }
    /// returns hessians of the base functions
    inline const base_tensor &hess(size_type i) const
      { if (hpc.empty()) init_hess(); return hpc[i]; }
    inline pfem get_pfem() const { return pf; }
    // inline const bgeot::stored_point_tab& get_point_tab() const
    //  { return *pspt; }
    inline bgeot::pstored_point_tab get_ppoint_tab() const
    { return pspt; }
    fem_precomp_(const pfem, const bgeot::pstored_point_tab);
    ~fem_precomp_() { DAL_STORED_OBJECT_DEBUG_DESTROYED(this, "Fem_precomp"); }
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
  pfem_precomp fem_precomp(pfem pf, bgeot::pstored_point_tab pspt,
                           dal::pstatic_stored_object dep);

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
        @param pspt a pointer to a list of points in the reference convex.
        
        CAUTION:
        this array must not be destroyed as long as the fem_precomp is used!!

        Moreover pspt is supposed to identify uniquely the set of
        points. This means that you should NOT alter its content until
        the fem_precomp_pool is destroyed.
    */
    pfem_precomp operator()(pfem pf, bgeot::pstored_point_tab pspt) {
      pfem_precomp p = fem_precomp(pf, pspt, 0);
      precomps.insert(p);
      return p;
    }
    void clear();
    ~fem_precomp_pool() { clear(); }
  };


  /** structure passed as the argument of fem interpolation
      functions. This structure can be partially filled (for example
      the xreal will be computed if needed as long as pgp+ii is known).
  */
  class fem_interpolation_context :
    public bgeot::geotrans_interpolation_context {

    mutable base_matrix M_; // optional transformation matrix (for non tau-equivalent fems)
    pfem pf_;               // current fem
    pfem_precomp pfp_;      // optional fem_precomp_ (speed up the computations)
    size_type convex_num_;  // The element (convex) number
    short_type face_num_;   // Face number for boundary integration
    int xfem_side_; // For the computation of a jump with fem_level_set only
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
    void base_value(base_tensor& t, bool withM = true) const;
    // Optimized function for high level generic assembly
    void pfp_base_value(base_tensor& t, const pfem_precomp &pfp__);
    /** fill the tensor with the gradient of the base functions (taken
        at point @c this->xref())
    */
    void grad_base_value(base_tensor& t, bool withM = true) const;
    // Optimized function for high level generic assembly
    void pfp_grad_base_value(base_tensor& t, const pfem_precomp &pfp__);
    /** fill the tensor with the hessian of the base functions (taken
        at point @c this->xref())
    */
    void hess_base_value(base_tensor& t, bool withM = true) const;
    /** get the current FEM descriptor */
    const pfem pf() const { return pf_; }
    /** get the current convex number */
    size_type convex_num() const;
    bool is_convex_num_valid() const;
    void invalid_convex_num() { convex_num_ = size_type(-1); }
    /** set the current face number */
    void set_face_num(short_type f) { face_num_ = f; }
    /** get the current face number */
    short_type face_num() const;
    /** On a face ? */
    bool is_on_face() const;
    /** get the current fem_precomp_ */
    pfem_precomp pfp() const { return pfp_; }
    void set_pfp(pfem_precomp newpfp);
    void set_pf(pfem newpf);
    int xfem_side() const { return xfem_side_; }
    void set_xfem_side(int side) { xfem_side_ = side; }
    void change(bgeot::pgeotrans_precomp pgp__,
		pfem_precomp pfp__, size_type ii__,
		const base_matrix& G__, size_type convex_num__,
		short_type face_num__ = short_type(-1)) {
      bgeot::geotrans_interpolation_context::change(pgp__,ii__,G__);
      convex_num_ = convex_num__; face_num_ = face_num__;  xfem_side_ = 0;
      set_pfp(pfp__);
    }
    void change(bgeot::pgeometric_trans pgt__,
		pfem_precomp pfp__, size_type ii__,
		const base_matrix& G__, size_type convex_num__,
		short_type face_num__ = short_type(-1)) {
      bgeot::geotrans_interpolation_context::change
	(pgt__, pfp__->get_ppoint_tab(), ii__, G__);
      convex_num_ = convex_num__; face_num_ = face_num__; xfem_side_ = 0;
      set_pfp(pfp__);
    }
    void change(bgeot::pgeometric_trans pgt__,
		pfem pf__, const base_node& xref__, const base_matrix& G__,
		size_type convex_num__,	short_type face_num__=short_type(-1)) {
      bgeot::geotrans_interpolation_context::change(pgt__,xref__,G__);
      pf_ = pf__; pfp_ = 0; convex_num_ = convex_num__; face_num_ = face_num__;
      xfem_side_ = 0;
    }
    fem_interpolation_context()
      : bgeot::geotrans_interpolation_context(),
	convex_num_(size_type(-1)), face_num_(short_type(-1)), xfem_side_(0) {}
    fem_interpolation_context(bgeot::pgeotrans_precomp pgp__,
                              pfem_precomp pfp__, size_type ii__,
                              const base_matrix& G__,
                              size_type convex_num__,
                              short_type face_num__ = short_type(-1))
      : bgeot::geotrans_interpolation_context(pgp__,ii__,G__),
	convex_num_(convex_num__), face_num_(face_num__), xfem_side_(0)
    { set_pfp(pfp__); }
    fem_interpolation_context(bgeot::pgeometric_trans pgt__,
                              pfem_precomp pfp__, size_type ii__,
                              const base_matrix& G__,
                              size_type convex_num__,
                              short_type face_num__ = short_type(-1))
      : bgeot::geotrans_interpolation_context(pgt__,pfp__->get_ppoint_tab(),
					      ii__, G__),
	convex_num_(convex_num__), face_num_(face_num__), xfem_side_(0)
    { set_pfp(pfp__); }
    fem_interpolation_context(bgeot::pgeometric_trans pgt__,
                              pfem pf__,
                              const base_node& xref__,
                              const base_matrix& G__,
                              size_type convex_num__,
                              short_type face_num__ = short_type(-1))
      : bgeot::geotrans_interpolation_context(pgt__,xref__,G__),
	pf_(pf__), pfp_(0), convex_num_(convex_num__), face_num_(face_num__),
	xfem_side_(0) {}
  };

  // IN : coeff(Qmult,nb_dof)
  // OUT: val(Qdim), Qdim=target_dim*Qmult
  // AUX: Z(nb_dof,target_dim)
  template <typename CVEC, typename VVEC>
  void virtual_fem::interpolation(const fem_interpolation_context& c,
                                  const CVEC& coeff, VVEC &val,
                                  dim_type Qdim) const {
    size_type Qmult = size_type(Qdim) / target_dim();
    size_type nbdof = nb_dof(c.convex_num());
    GMM_ASSERT1(gmm::vect_size(val) == Qdim, "dimensions mismatch");
    GMM_ASSERT1(gmm::vect_size(coeff) == nbdof*Qmult,
		"Wrong size for coeff vector");

    gmm::clear(val);
    base_tensor Z; real_base_value(c, Z);

    for (size_type j = 0; j < nbdof; ++j) {
      for (size_type q = 0; q < Qmult; ++q) {
        typename gmm::linalg_traits<CVEC>::value_type co = coeff[j*Qmult+q];
        for (size_type r = 0; r < target_dim(); ++r)
          val[r + q*target_dim()] += co * Z[j + r*nbdof];
      }
    }
  }

  template <typename MAT>
  void virtual_fem::interpolation(const fem_interpolation_context& c,
                                  MAT &M, dim_type Qdim) const {
    size_type Qmult = size_type(Qdim) / target_dim();
    size_type nbdof = nb_dof(c.convex_num());
    GMM_ASSERT1(gmm::mat_nrows(M) == Qdim && gmm::mat_ncols(M) == nbdof*Qmult,
                "dimensions mismatch");

    gmm::clear(M);
    base_tensor Z; real_base_value(c, Z);
    for (size_type j = 0; j < nbdof; ++j) {
      for (size_type q = 0; q < Qmult; ++q) {
        for (size_type r = 0; r < target_dim(); ++r)
          M(r+q*target_dim(), j*Qmult+q) = Z[j + r*nbdof];
      }
    }
  }


  // IN : coeff(Qmult,nb_dof)
  // OUT: val(Qdim,N), Qdim=target_dim*Qmult
  // AUX: t(nb_dof,target_dim,N)
  template<typename CVEC, typename VMAT>
  void virtual_fem::interpolation_grad(const fem_interpolation_context& c,
                                       const CVEC& coeff, VMAT &val,
                                       dim_type Qdim) const {
    size_type N = c.N();
    size_type nbdof = nb_dof(c.convex_num());
    size_type Qmult = gmm::vect_size(coeff) / nbdof;
    GMM_ASSERT1(gmm::mat_ncols(val) == N &&
                gmm::mat_nrows(val) == target_dim()*Qmult &&
                gmm::vect_size(coeff) == nbdof*Qmult,
                "dimensions mismatch");
    GMM_ASSERT1(Qdim == target_dim()*Qmult, // Qdim seems to be superfluous input, could be removed in the future
                "dimensions mismatch");
    base_tensor t;
    real_grad_base_value(c, t); // t(nbdof,target_dim,N)

    gmm::clear(val);
    for (size_type q = 0; q < Qmult; ++q) {
      base_tensor::const_iterator it = t.begin();
      for (size_type k = 0; k < N; ++k)
        for (size_type r = 0; r < target_dim(); ++r)
          for (size_type j = 0; j < nbdof; ++j, ++it)
            val(r + q*target_dim(), k) += coeff[j*Qmult+q] * (*it);
    }
  }


  template<typename CVEC, typename VMAT>
  void virtual_fem::interpolation_hess(const fem_interpolation_context& c,
                                       const CVEC& coeff, VMAT &val,
                                       dim_type Qdim) const {
    size_type Qmult = size_type(Qdim) / target_dim();
    size_type N = c.N();
    GMM_ASSERT1(gmm::mat_ncols(val) == N*N &&
                gmm::mat_nrows(val) == Qdim, "dimensions mismatch");

    base_tensor t;
    size_type nbdof = nb_dof(c.convex_num());

    gmm::clear(val);
    real_hess_base_value(c, t);
    for (size_type q = 0; q < Qmult; ++q) {
      base_tensor::const_iterator it = t.begin();
      for (size_type k = 0; k < N*N; ++k)
        for (size_type r = 0; r < target_dim(); ++r)
          for (size_type j = 0; j < nbdof; ++j, ++it)
            val(r + q*target_dim(), k) += coeff[j*Qmult+q] * (*it);
    }
  }


  // IN : coeff(Qmult,nb_dof)
  // OUT: val
  // AUX: t(nb_dof,target_dim,N), Qmult*target_dim == N
  template<typename CVEC>
  void virtual_fem::interpolation_diverg
    (const fem_interpolation_context& c, const CVEC& coeff,
     typename gmm::linalg_traits<CVEC>::value_type &val) const {
    size_type N = c.N();
    size_type nbdof = nb_dof(c.convex_num());
    size_type Qmult = gmm::vect_size(coeff) / nbdof;
    GMM_ASSERT1(gmm::vect_size(coeff) == nbdof*Qmult , "dimensions mismatch");
    GMM_ASSERT1(target_dim()*Qmult == N &&
                (Qmult == 1 || target_dim() == 1),
                "Dimensions mismatch. Divergence operator requires fem qdim equal to dim.");
    base_tensor t;
    real_grad_base_value(c, t); // t(nbdof,target_dim,N)
    // for Qmult == 1 this is sub-optimal since it evaluates all (:,i,j)
    // gradients instead of only the diagonal ones(:,i,i)

    val = scalar_type(0);
    base_tensor::const_iterator it = t.begin();
    if (Qmult == 1)
      for (size_type k = 0; k < N; ++k) {
        if (k) it += (N*nbdof + 1);
        for (size_type j = 0; j < nbdof; ++j) {
          if (j) ++it;
          val += coeff[j] * (*it);
        }
      }
    else // if (target_dim() == 1)
      for (size_type k = 0; k < N; ++k) {
        if (k) ++it;
        for (size_type j = 0; j < nbdof; ++j) {
          if (j) ++it;
          val += coeff[j*N+k] * (*it);
        }
      }
  }

  /**
     Specific function for a HHO method to obtain the method in the interior.
     If the method is not of composite type, return the argument.
  */
  pfem interior_fem_of_hho_method(pfem hho_method);


  /* Functions allowing the add of a finite element method outside
     of getfem_fem.cc */

  typedef dal::naming_system<virtual_fem>::param_list fem_param_list;

  void inline read_poly(bgeot::base_poly &p, int d, const char *s)
  { p = bgeot::read_base_poly(short_type(d), s); }

  void add_fem_name(std::string name,
                    dal::naming_system<virtual_fem>::pfunction f);

  
  /* @} */

}  /* end of namespace getfem.                                            */


#endif
