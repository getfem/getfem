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


#ifndef __GETFEM_FEM_H
#define __GETFEM_FEM_H

#include <bgeot_convex_ref.h>
#include <bgeot_geometric_trans.h>
#include <getfem_config.h>
#include <getfem_precomp.h>
#include <getfem_mesh.h>
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
  
  pdof_description change_coord_index_dof(pdof_description, dim_type);
  /** Gives a total order on the dof description compatible with the
   *   identification.
   */
  int dof_description_compare(pdof_description a, pdof_description b);
  /// Says if the dof is linkable.
  bool dof_linkable(pdof_description);
  /// Says if the two dofs can be identified.
  bool dof_compatibility(pdof_description, pdof_description);
  
  /* ******************************************************************** */
  /*	Classes for description of a finite element.                      */
  /* ******************************************************************** */
  
  /** Describe a finite element method on a reference element. At the
   *  level of the description of nodes, basic functions, ...
   *  This class is not to be manipulate by itself. Use pfem\_interpolation
   *  and the functions written to produce the descriptions from the
   *  basic descriptions.
   *   \subsubsection{General finite element method}
   *   A general finite element method is described by the set of the
   *   three following things:
   *   \begin{itemize}
   *      \item An element $T$ (geometric figure, polyhedric and convex),
   *      \item A functional vectorial space $V$, of finite dimension $n_d$,
   *      \item A set of $n_d$ linear functional \\
   *    $ l_i : V \longrightarrow {I\hspace{-0.3em}R}, \ \ \ i = 1..n_g. $ \\
   *   \end{itemize}
   *   Very often, each linear functional $l_i$ will be associated with
   *   a node $a^i$, which is a particular point of $T$. \\
   *   The shape functions of the finite element are the functions
   *   $\varphi_i \in V$ which satisfy \\
   *   $ l_i(\varphi_j) = \delta_{ij}, i,j = 1..n_g. $
   *   The finite element is said "unisolvant" when the set of shape functions
   *   are uniquely defined. \\
   *   Often $V$ is a space of polynomials, but not always. \\
   *   
   *   \subsubsection{Considered finite element methods}
   *   The finite element methods treaten in {\sc Getfem++}, are
   *   the methods satisfying the following assumptions :
   *  
   *   \begin{itemize}
   *     \item the method is equivalent throw a polynomial geometric
   *       transformation to a method defined on a reference element.
   *       This means that there exist a reference element $\overline{T}$, and
   *       a transformation \\
   *       $ \begin{array}{rcl} \tau : \overline{T} &\longrightarrow& T, \\
   *          \overline{x} &\longmapsto& x, \end{array} $ \\
   *       which is polynomial.
   *     \item An additionnal linear transformation is allowed: \\
   *     $ \tilde{\varphi}_i(x) =  \overline{\varphi}_i(\tau(\overline{x})),$\\
   *     and  \\
   *     $ \varphi_i = \sum_{j = 1}^{n_d} \tilde{M}_{ij} \tilde{\varphi}_j $\\
   *     This additionnal matrix may be dependant on the real element and allow
   *     non-equivalent element such as Hermite elements. (...)
   *   \end{itemize}
   */
  
  class virtual_fem;
  typedef const virtual_fem * pfem;
  
  class virtual_fem
  {
  protected :

    std::vector<pdof_description> _dof_types;
    bgeot::convex_structure cvs_node;
    bgeot::convex<base_node> cv_node;
    mutable bgeot::pstored_point_tab pspt;
    mutable bool pspt_valid;
    bgeot::pconvex_ref cvr; // reference element.
    dim_type ntarget_dim;
    bool is_equiv, is_lag, is_pol, do_grad;
    short_type es_degree;
    
    void add_node(const pdof_description &d, const base_node &pt) ;
    void init_cvs_node(void);
    
  public :
    /// Number of degrees of freedom.
    virtual size_type nb_dof(void) const { return _dof_types.size(); }
    /// Number of components (nb_dof() * dimension of the target space).
    virtual size_type nb_base(void) const { return nb_dof(); }
    size_type nb_base_components(void) const
      { return nb_base() * ntarget_dim; }
    size_type nb_components(void) const
      { return nb_dof() * ntarget_dim; }
    /// Gives the array of pointer on dof description.
    const std::vector<pdof_description> &dof_types(void) const 
      { return _dof_types; }
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
    bgeot::pconvex_ref ref_convex(void) const
      { return cvr; }
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
    bool is_equivalent(void) const { return is_equiv; }
    bool is_lagrange(void) const { return is_lag; }
    bool is_polynomial(void) const { return is_pol; }
    short_type estimated_degree(void) const { return es_degree; }
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
			   bgeot::pgeometric_trans pgt) const
      { DAL_THROW(internal_error, "This function should not be called."); }
    virtual void interpolation(const base_node &x, const base_matrix &G,
			       bgeot::pgeometric_trans pgt,
			       const base_vector coeff, 
			       base_node &val) const = 0;
    virtual void interpolation(pfem_precomp pfp, size_type ii,
			       const base_matrix &G,
			       bgeot::pgeometric_trans pgt, 
			       const base_vector coeff, base_node &val) const;
    virtual void interpolation_grad(const base_node &x, const base_matrix &G,
				    bgeot::pgeometric_trans pgt,
				    const base_vector coeff,
				    base_matrix &val) const = 0;
    virtual void complete_interpolation_grad(const base_node &x,
					     const base_matrix &G,
					     bgeot::pgeometric_trans pgt,
					     const base_vector coeff,
					     base_matrix &val) const;
    virtual void interpolation_grad(pfem_precomp pfp, size_type ii,
				    const base_matrix &G,
				    bgeot::pgeometric_trans pgt, 
				    const base_vector coeff,
				    base_matrix &val) const;

    /** Gives the value of all components of the base functions at the
     *  point x of the reference element.
     */
    virtual void base_value(const base_node &x, base_tensor &t) const = 0;
    /** Gives the value of all gradients of the components of the
     *  base functions at the point x of the reference element.
     */
    virtual void grad_base_value(const base_node &x, base_tensor &t) const = 0;
    /** Gives the value of all hessians of the components of the base
     *  functions at the point x of the reference element.
     */
    virtual void hess_base_value(const base_node &x, base_tensor &t) const = 0;
    
    virtual size_type index_of_already_numerate_dof(size_type, size_type) const
      { DAL_THROW(internal_error, "internal error."); return 0; }
    bool do_grad_reduction(void) const { return do_grad; }

    virtual_fem(void) { 
      ntarget_dim = 1; is_equiv = is_pol = is_lag = false;
      pspt_valid = false; do_grad = true;
    }
  };
  
  
  
  template <class FUNC> class fem : public virtual_fem
  {
  protected :
    
    std::vector<FUNC> _base;
    
  public :
    
    /// Gives the array of basic functions (components).
    const std::vector<FUNC> &base(void) const { return _base; }
    
    void interpolation(const base_node &x, const base_matrix &G, 
		       bgeot::pgeometric_trans pgt,
		       const base_vector coeff, base_node &val) const;
    void interpolation_grad(const base_node &x, const base_matrix &G,
			    bgeot::pgeometric_trans pgt,
			    const base_vector coeff, base_matrix &val) const;
    void base_value(const base_node &x, base_tensor &t) const {
      bgeot::multi_index mi(2);
      mi[1] = target_dim(); mi[0] = nb_base();
      // cout << "mi = " << mi << endl;
      t.adjust_sizes(mi);
      size_type R = nb_base_components();
      base_tensor::iterator it = t.begin();
      for (size_type  i = 0; i < R; ++i, ++it) {
	// cout << "base " << i <<  _base[i] << endl;
	*it = _base[i].eval(x.begin());
      }
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
	  { FUNC f = _base[i]; f.derivative(j); *it = f.eval(x.begin()); }
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
	    FUNC f = _base[i];
	    f.derivative(j); f.derivative(k);
	    *it = f.eval(x.begin());
	  }
    }
    
    
  };
  
  typedef const fem<base_poly> * ppolyfem;
  
  template <class FUNC>
  void fem<FUNC>::interpolation(const base_node &x, const base_matrix &G,
				bgeot::pgeometric_trans pgt, 
				const base_vector coeff,
				base_node &val) const { 
    // optimisable.   verifier et faire le vectoriel
    base_matrix M;
    if (val.size() != target_dim())
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    size_type R = nb_dof(), RR = nb_base();
    
    if (!is_equivalent()) { M.resize(RR, R); mat_trans(M, G, pgt); }
    
    val.fill(0.0);
    for (size_type j = 0; j < RR; ++j) {
      scalar_type co = 0.0;
      if (is_equivalent())
	co = coeff[j];
      else
	for (size_type i = 0; i < R; ++i)
	  co += coeff[i] * M(i, j);
      
      for (size_type r = 0; r < target_dim(); ++r)
	val[r] += co * base()[j + r*R].eval(x.begin());
    } 
  }
  
  
  template <class FUNC>
  void fem<FUNC>::interpolation_grad(const base_node &x,
				     const base_matrix &G,
				     bgeot::pgeometric_trans pgt, 
				     const base_vector coeff,
				     base_matrix &val) const { 
    // optimisable.   verifier
    base_matrix M;
    base_tensor t;
    
    if (val.nrows() != target_dim() || val.ncols() != x.size())
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    grad_base_value(x, t);
    base_tensor::iterator it = t.begin();
    
    size_type R = nb_dof(), RR = nb_base();
    
    if (!is_equivalent()) { M.resize(RR, R); mat_trans(M, G, pgt); }
    
    val.fill(0.0);
    
    for (size_type k = 0; k < x.size(); ++k)
      for (size_type r = 0; r < target_dim(); ++r)
	for (size_type j = 0; j < RR; ++j, ++it) {
	  scalar_type co = 0.0;
	  if (is_equivalent())
	    co = coeff[j];
	  else
	    for (size_type i = 0; i < R; ++i)
	      co += coeff[i] * M(i, j);
	  
	  val(r,k) += co * (*it);
	} 
  }
  
  
  /** @name functions on convex structures
   */
  //@{
  
  /** Description of the classical finite element Pk on simplex of
   *   dimension n.
   */
  ppolyfem PK_fem(dim_type n, short_type k);
  /** Description of the classical finite element Qk on parallelepiped of
   *   dimension n.
   */
  ppolyfem QK_fem(dim_type n, short_type k);
  /// Description of the tensorial product of the two elements *pi1, *pi2.
  ppolyfem product_fem(ppolyfem pi1, ppolyfem pi2);
  /** Description of the classical finite element Pk x Pk on prism of
   *   dimension n.
   */
  inline ppolyfem PK_prism_fem(dim_type nc, short_type k) {
    return product_fem(PK_fem(nc-1, k), PK_fem(1, k));
  }
  /// Description of the P1 non conforming finite element on triangles.
  ppolyfem P1_nonconforming_fem(void);
  
  /** Description of a P1 element on a simplex with the adjonction
   *  of a bubble base fonction on a face.
   */
  ppolyfem P1_with_bubble_on_a_face(dim_type n);
  ppolyfem P1_with_bubble_on_a_face_lagrange(void);
  
  /** PK element, where the dofs are not linked with the elements of 
      the neightbouring elements */
  ppolyfem PK_discontinuous_fem(dim_type n, short_type k);
  
  
  /** PK element on with a bubble base fonction has been added 
      (hence it is not a lagrange element)
  */
  ppolyfem PK_with_cubic_bubble_fem(dim_type n, short_type k);
  
  /** Hermite element on the segment
   */
  ppolyfem segment_Hermite_fem(void);
  
  /** Gives a pointer on the structures describing the more classical fem
   *  of degree k on a geometric convex cvs (coming from the geometric trans).
   */
  pfem classical_fem(bgeot::pgeometric_trans pg, short_type k);

  
  class mesh_fem;
  class pintegration_method;
  pfem virtual_link_fem(mesh_fem &mf1, mesh_fem &mf2,
			pintegration_method pim);
  pfem virtual_link_fem_with_gradient(mesh_fem &mf1, mesh_fem &mf2,
				      pintegration_method pim);
  
  //@}
  
}  /* end of namespace getfem.                                            */


#endif
