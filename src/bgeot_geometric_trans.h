// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_geometric_trans.h : geometric transformations on convex
//           
// Date    : December 20, 2000.
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



#ifndef BGEOT_GEOMETRIC_TRANSFORMATION_H__
#define BGEOT_GEOMETRIC_TRANSFORMATION_H__

#include <set>
#include <bgeot_config.h>
#include <bgeot_convex_ref.h>

namespace bgeot
{
  /**  Description of a geometric transformation between a
   * reference element and a real element.
   * Geometric  nodes  and  vector  of polynomials. This class 
   * is not  to be  manipulate  by itself.  Use pgeometric\_trans  and
   * the   functions   written   to   produce   the   basic  geometric
   * transformations.
   *   \subsubsection*{Description of the geometry}
   *     Let $T \in\ {I\hspace{-0.3em}R}^N$ be a real element and
   *     $\overline{T} \in\ {I\hspace{-0.3em}R}^P$ be a reference element,
   *      with $N >= P$. \\
   *     The geometric nodes of $\overline{T}$ are the points
   *     $\overline{g}^i \in\ {I\hspace{-0.3em}R}^P$, for $i = 0 .. n_g-1$,
   *     and the corresponding (via the geometric transformation) nodes of
   *     $T$ are the points $g^i \in\ {I\hspace{-0.3em}R}^N$. 
   *   \subsubsection*{Geometric transformation}
   *     The geometric transformation is the application \\
   *     $ \begin{array}{rl}
   *        \tau : \overline{T} & \longrightarrow \ T, \\
   *               \overline{x} & \longmapsto \ \ x,
   *     \end{array} $ \\
   *     which should be a diffeomorphism between $\overline{T}$ and $T$. It
   *     is assumed that there exists a polynomial vector \\
   *     $ \underline{\cal N}(\overline{x})
   *        = \left({\cal N}i_(\overline{x})\right)i_, \ \ i = 0 .. n_g-1, $ \\
   *     defined on $\overline{T}$ of size $n_g$, such that the transformation
   *     $\tau$ can be written \\
   *     $ \tau(\overline{x}) = \sum_{i = 0}^{n_g-1} {\cal N}i_(\overline{x})
   *     g^i. $ \\
   *     Denoting by \\
   *     $ \underline{\underline{G}} = (g^0; g^1; ...;g^{n_g-1}), $ \\
   *     The matrix in which each column is a geometric node of $T$,
   *     the transformation $\tau$ can be written as \\
   *     $ \tau(\overline{x}) = \underline{\underline{G}} \ 
   *        \underline{\cal N}(\overline{x}). $ \\
   *   \subsubsection*{Gradient of the transformation}
   *     The gradient of the transformation is \\
   *     $ \nabla \tau(\overline{x}) = 
   *     \left( \frac{\partial \tau_i}{\partial \overline{x}j_} \right)_{ij}
   *     = \left( \sum_{l = 0}^{n_g-1}g^l_i
   *     \frac{\partial {\cal N}l_(\overline{x})}{\partial \overline{x}j_}
   *     \right)_{ij} = \underline{\underline{G}}\ \nabla
   *     \underline{\cal N}(\overline{x}), $ \\
   *     Remark : $\underline{\underline{G}}$ is a $N \times n_g$ matrix,
   *       $\nabla \underline{\cal N}(\overline{x})$ is a $n_g \times P$
   *       matrix, and thus $\nabla \tau(\overline{x})$ is a $N \times P$
   *       matrix.
   *   \subsubsection*{Inverse transformation and pseudo-inverse}
   *     to do ...
   */
  class geometric_trans : public dal::static_stored_object {
  protected :
    
    bool is_lin;
    pconvex_ref cvr;
    std::vector<base_poly> trans;
  public :
    
    /// Dimension of the reference element.
    dim_type dim(void) const { return cvr->structure()->dim(); }
    /// True if the transformation is linear (affine in fact).
    bool is_linear(void) const { return is_lin; }
    /// Number of geometric nodes.
    size_type nb_points(void) const { return cvr->nb_points(); }
    /// Pointer on the convex of reference.
    pconvex_ref convex_ref(void) const { return cvr; }
    /// Structure of the reference element.
    pconvex_structure structure(void) const { return cvr->structure(); }
    /// Basic structure of the reference element.
    pconvex_structure basic_structure(void) const
    { return cvr->structure()->basic_structure(); }
    /// Gives the vector of polynomials representing the transformation.
    const std::vector<base_poly> &poly_vector(void) const { return trans; }
    /// Gives the array of geometric nodes (on reference convex)
    const stored_point_tab &geometric_nodes(void) const
    { return cvr->points(); }
    /// Gives the array of the normals to faces (on reference convex)
    const std::vector<base_small_vector> &normals(void) const
    { return cvr->normals(); }
    /** Apply the geometric transformation to point pt, 
	PTAB containsis the points of the real convex */
    template<class CONT> base_node transform(const base_node &pt,
					     const CONT &PTAB) const;
    base_node transform(const base_node &pt, const base_matrix &G) const;
    /** Compute the gradient at point x, pc is resized to [nb_points() x dim()]
	if the transformation is linear, x is not used at all */
    void gradient(const base_node& x, base_matrix& pc) const;
    virtual ~geometric_trans() {}
  };

  template<class CONT>
    base_node geometric_trans::transform(const base_node &pt,
					 const CONT &ptab) const {
    base_node P(ptab[0].size()); P.fill(0.0);
    size_type k = nb_points();
    for (size_type l = 0; l < k; ++l)
      gmm::add(gmm::scaled(ptab[l],
			   scalar_type(poly_vector()[l].eval(pt.begin()))),P);
      //P.addmul(poly_vector()[l].eval(pt.begin()),ptab[l]);
    return P;
  }

  typedef boost::intrusive_ptr<const geometric_trans> pgeometric_trans;
  class geotrans_interpolation_context;

  template<class CONT>
  void bounding_box(base_node& min, base_node& max, 
		    const CONT &ptab, pgeometric_trans pgt = 0) {
    typename CONT::const_iterator it = ptab.begin();
    min = max = *it; size_type P = min.size();
    base_node::iterator itmin = min.begin(), itmax = max.begin();
    for ( ++it; it != ptab.end(); ++it) {
      base_node pt = *it; /* need a temporary storage since cv.points()[j] may
			     not be a reference to a base_node, but a
			     temporary base_node !! (?) */
      base_node::const_iterator it2 = pt.begin();
      for (size_type i = 0; i < P; ++i) {
	itmin[i] = std::min(itmin[i], it2[i]);
	itmax[i] = std::max(itmax[i], it2[i]);
      }
    }
    /* enlarge the box for non-linear transformations .. */
    if (pgt && !pgt->is_linear()) 
      for (size_type i = 0; i < P; ++i) {
	scalar_type e = (itmax[i]-itmin[i]) * 0.2;
	itmin[i] -= e; itmax[i] += e;
      }
  }

  /** @name functions on geometric transformations
   */
  //@{

  pgeometric_trans simplex_geotrans(size_type n, short_type k);
  pgeometric_trans parallelepiped_geotrans(size_type n, short_type k);
  pgeometric_trans parallelepiped_linear_geotrans(size_type n);
  pgeometric_trans prism_geotrans(size_type n, short_type k);
  pgeometric_trans prism_linear_geotrans(size_type n);
  pgeometric_trans product_geotrans(pgeometric_trans pg1,
				    pgeometric_trans pg2);
  pgeometric_trans linear_product_geotrans(pgeometric_trans pg1,
					   pgeometric_trans pg2);

  pgeometric_trans geometric_trans_descriptor(std::string name);
  /* List :
   * GT_PK(N,K)   : Transformation on simplexes, dim N, degree K
   * GT_QK(N,K)   : Transformation on parallelepipeds, dim N, degree K
   * GT_PRISM(N,K)          : Transformation on prisms, dim N, degree K
   * GT_PRODUCT(a,b)        : tensorial product of two transformations
   * GT_LINEAR_PRODUCT(a,b) : Linear tensorial product of two transformations
   * GT_LINEAR_QK(N) : shortcut for GT_LINEAR_PRODUCT(GT_LINEAR_QK(N-1),
   *                                                  GT_PK(1,1))
   */

  std::string name_of_geometric_trans(pgeometric_trans p);

  /** norm of returned vector is the ratio between the face surface on
   *  the reel element and the face surface on the reference element 
   *  IT IS NOT UNITARY
   *
   *  pt is the position of the evaluation point on the reference element
   */
  base_small_vector compute_normal(const geotrans_interpolation_context& c,
				   size_type face);

  /** return the local basis (i.e. the norm in the first column, and the
   *  tangent vectors in the other columns 
   */
  base_matrix 
  compute_local_basis(const geotrans_interpolation_context& c,
		      size_type face);
    //@}

  /* ********************************************************************* */
  /*       Precomputation on geometric transformations.                    */
  /* ********************************************************************* */
  
  /**
   *  precomputed geometric transformation operations use this for
   *  repetitive evaluation of a geometric transformations on a set of
   *  points "pspt" in the the reference convex which do not change.
   */
  class geotrans_precomp_ : public dal::static_stored_object {
  protected:      
    pgeometric_trans pgt;
    pstored_point_tab pspt;  /* a set of points in the reference elt*/
    mutable std::vector<base_vector> c;  /* precomputed values for the     */
                                         /* transformation                 */
    mutable std::vector<base_matrix> pc; /* precomputed values for gradient*/
                                         /* of the transformation.         */
    mutable std::vector<base_matrix> hpc; /* precomputed values for hessian*/
                                          /*  of the transformation.       */
  public:
    inline const base_vector &val(size_type i) const 
    { if (c.empty()) init_val(); return c[i]; }
    inline const base_matrix &grad(size_type i) const
    { if (pc.empty()) init_grad(); return pc[i]; }
    inline const base_matrix &hessian(size_type i) const
    { if (hpc.empty()) init_hess(); return hpc[i]; }
    
    /**
     *  Apply the geometric transformation from the reference convex to
     *  the convex whose vertices are stored in G, to the set of points
     *  listed in pspt. @param G any container of vertices of the transformed
     *  convex
     *  @result pt_tab transformed points
     */
    template <typename CONT> 
    void transform(const CONT& G,
		   stored_point_tab& pt_tab) const;
    template <typename CONT, typename VEC> 
    void transform(const CONT& G, size_type ii, VEC& pt) const;

    base_node transform(size_type i, const base_matrix &G) const;
    pgeometric_trans get_trans() const { return pgt; }
    inline const stored_point_tab& get_point_tab() const { return *pspt; }
    geotrans_precomp_(pgeometric_trans pg, pstored_point_tab ps);
  private:
    void init_val() const;
    void init_grad() const;
    void init_hess() const;
  };


  template <typename CONT, typename VEC> 
  void geotrans_precomp_::transform(const CONT& G, size_type j,
				    VEC& pt) const {
    size_type k = 0;
    if (c.empty()) init_val();
    for (typename CONT::const_iterator itk = G.begin(); 
         itk != G.end(); ++itk, ++k) {
      typename CONT::value_type::const_iterator Gk = (*itk).begin();
      typename VEC::iterator ipt = pt.begin();
      for (size_type i=0; i < (*itk).size(); ++i) {
        ipt[i] += Gk[i] * c[j][k];
      }
    }
  }

  template <typename CONT> 
  void geotrans_precomp_::transform(const CONT& G,
                                    stored_point_tab& pt_tab) const {
    if (c.empty()) init_val();
    pt_tab.clear(); pt_tab.resize(c.size(), base_node(G[0].size()));
    for (size_type j = 0; j < c.size(); ++j) {
      transform(G, j, pt_tab[j]);
    }
  }
  
  typedef boost::intrusive_ptr<const geotrans_precomp_> pgeotrans_precomp;
  
  /**
   *  precomputes a geometric transformation for a fixed set of 
   *  points in the reference convex. 
   */
  pgeotrans_precomp geotrans_precomp(pgeometric_trans pg,
				     pstored_point_tab pspt);

  void delete_geotrans_precomp(pgeotrans_precomp pgp);

  /**
   *  The object geotrans_precomp_pool Allow to allocate a certain number
   *  of geotrans_precomp and automatically delete them when it is 
   *  deleted itself.
   */
  class geotrans_precomp_pool {
    std::set<pgeotrans_precomp> precomps;

  public :
    
    pgeotrans_precomp operator()(pgeometric_trans pg, pstored_point_tab pspt) {
      pgeotrans_precomp p = geotrans_precomp(pg, pspt);
      precomps.insert(p);
      return p;
    }
    ~geotrans_precomp_pool() {
      for (std::set<pgeotrans_precomp>::iterator it = precomps.begin();
	   it != precomps.end(); ++it)
	delete_geotrans_precomp(*it);
    }
  };



  /* the geotrans_interpolation_context structure is passed as the
     argument of geometric transformation interpolation
     functions. This structure can be partially filled (for example
     the xreal will be computed if needed as long as pgp+ii is known).
     See also fem_interpolation_context in getfem_fem.h.
     The name of member data, and the computations done by this structure
     are heavily described in the Getfem++ Kernel Documentation.
  */
  class geotrans_interpolation_context {
    mutable base_node xref_; /* reference point */
    mutable base_node xreal_; /* transformed point */
    const base_matrix *G_; /* pointer to the matrix of real nodes of the convex */
    mutable base_matrix B_, B3_, B32_; /* see documentation for more details */
    pgeometric_trans pgt_;
    pgeotrans_precomp pgp_;
    size_type ii_; /* index of current point in the pgp */
    mutable scalar_type J_; /* Jacobian */
    void compute_J(void) const;
  public:
    bool have_xref() const { return !xref_.empty(); }
    bool have_xreal() const { return !xreal_.empty(); }
    bool have_G() const { return G_ != 0; }
    bool have_B() const { return !B_.empty(); }
    bool have_B3() const { return !B3_.empty(); }
    bool have_B32() const { return !B32_.empty(); }
    bool have_pgt() const { return pgt_ != 0; }
    bool have_pgp() const { return pgp_ != 0; }
    const base_node& xref() const;
    const base_node& xreal() const;
    const base_matrix& B() const;
    const base_matrix& B3() const;
    const base_matrix& B32() const;
    bgeot::pgeometric_trans pgt() const { return pgt_; }
    const base_matrix& G() const { return *G_; }
    scalar_type J() const { if (J_ < scalar_type(0)) compute_J(); return J_; }
    size_type N() const { if (have_G()) return G().nrows(); 
      else if (have_xreal()) return xreal_.size(); 
      else DAL_THROW(dal::failure_error, "cannot get N"); }
    size_type ii() const { return ii_; }
    bgeot::pgeotrans_precomp pgp() const { return pgp_; }
    void set_ii(size_type ii__);
    void set_xref(const base_node& P);

    geotrans_interpolation_context();
    geotrans_interpolation_context(bgeot::pgeotrans_precomp pgp__, 
				   size_type ii__, 
				   const base_matrix& G__); 
    geotrans_interpolation_context(bgeot::pgeometric_trans pgt__,
				   const base_node& xref__,
			 	   const base_matrix& G__);
  };

}  /* end of namespace bgeot.                                             */


#endif /* BGEOT_GEOMETRIC_TRANSFORMATION_H__                              */
