/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_geometric_trans.h : geometric transformations on convex*/
/*     									   */
/*                                                                         */
/* Date : December 20, 2000.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2000-2002  Yves Renard.                                   */
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



#ifndef __BGEOT_GEOMETRIC_TRANSFORMATION_H
#define __BGEOT_GEOMETRIC_TRANSFORMATION_H

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
   *        = \left({\cal N}_i(\overline{x})\right)_i, \ \ i = 0 .. n_g-1, $ \\
   *     defined on $\overline{T}$ of size $n_g$, such that the transformation
   *     $\tau$ can be written \\
   *     $ \tau(\overline{x}) = \sum_{i = 0}^{n_g-1} {\cal N}_i(\overline{x})
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
   *     \left( \frac{\partial \tau_i}{\partial \overline{x}_j} \right)_{ij}
   *     = \left( \sum_{l = 0}^{n_g-1}g^l_i
   *     \frac{\partial {\cal N}_l(\overline{x})}{\partial \overline{x}_j}
   *     \right)_{ij} = \underline{\underline{G}}\ \nabla
   *     \underline{\cal N}(\overline{x}), $ \\
   *     Remark : $\underline{\underline{G}}$ is a $N \times n_g$ matrix,
   *       $\nabla \underline{\cal N}(\overline{x})$ is a $n_g \times P$
   *       matrix, and thus $\nabla \tau(\overline{x})$ is a $N \times P$
   *       matrix.
   *   \subsubsection*{Inverse transformation and pseudo-inverse}
   *     to do ...
   */
  class geometric_trans {
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
    /// Gives the array of geometric nodes.
    const std::vector<base_node> &geometric_nodes(void) const
    { return cvr->points(); }
    /// Gives the array of the normals to faces.
    const std::vector<base_vector> &normals(void) const
    { return cvr->normals(); }
    
    template<class CONT> base_node transform(const base_node &pt,
					     const CONT &PTAB) const;
    base_node transform(const base_node &pt, const base_matrix &G) const;
  };

  template<class CONT>
    base_node geometric_trans::transform(const base_node &pt,
					 const CONT &ptab) const {
    base_node P(ptab[0].size()); P.fill(0.0);
    size_type k = nb_points();
    for (size_type l = 0; l < k; ++l)
      P.addmul(poly_vector()[l].eval(pt.begin()),ptab[l]);
    return P;
  }

  /** @name functions on geometric transformations
   */
  //@{

  typedef const geometric_trans *pgeometric_trans;

  pgeometric_trans simplex_geotrans(size_type n, short_type k);
  pgeometric_trans parallelepiped_geotrans(size_type n, short_type k);
  pgeometric_trans parallelepiped_linear_geotrans(size_type n);
  pgeometric_trans prism_geotrans(size_type n, short_type k);
  pgeometric_trans prism_linear_geotrans(size_type n);
  pgeometric_trans product_geotrans(pgeometric_trans pg1,pgeometric_trans pg2);
  pgeometric_trans linear_product_geotrans(pgeometric_trans pg1,
					   pgeometric_trans pg2);

  pgeometric_trans geometric_trans_descriptor(std::string name);
  /* List :
   * GT_PK(N,K)   : Transformation on simplexes, dim N, degree K
   * GT_QK(N,K)   : Transformation on parallelepipeds, dim N, degree K
   * GT_PRISM(N,K)          : Transformation on prisms, dim N, degree K
   * GT_PRODUCT(a,b)        : tensorial product of two transformations
   * GT_LINEAR_PRODUCT(a,b) : Linear tensorial product of two transformations
   */

  std::string name_of_geometric_trans(pgeometric_trans p);

  /* norm of returned vector is the ratio between the face surface on
     the reel element and the face surface on the reference element 
     IT IS NOT UNITARY

     pt is the position of the evaluation point on the reference element
  */
  base_vector compute_normal(const base_matrix &G, size_type ir,
			     pgeometric_trans pgt, const base_node &pt);

   //@}

}  /* end of namespace bgeot.                                             */


#endif /* __BGEOT_GEOMETRIC_TRANSFORMATION_H                              */
