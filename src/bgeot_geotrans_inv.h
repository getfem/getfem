/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_geotrans_inv.h : Allow to inverse geometric transf-    */
/*     	      ormations and to localize a set of points.    	       	   */
/*                                                                         */
/* Date : December 20, 2000.                                               */
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


#ifndef __BGEOT_GEOTRANS_INV_H
#define __BGEOT_GEOTRANS_INV_H

#include <bgeot_geometric_trans.h>

namespace bgeot
{
  
  template<class T> inline T sfloor(T x)
  { return (x >= 0) ? floor(x) : -floor(-x); }

  /// Sort order for a rapid localization in bgeot::geotrans_inv
  struct imbricated_box_less
    : public std::binary_function<base_node, base_node, int>
  { 
    int exp_max, exp_min;
    scalar_type c_max;

    /// comparaison function
    int operator()(const base_node &x, const base_node &y) const;
    
    imbricated_box_less(int emi = -15, int ema = -2)
    { exp_max = ema; exp_min = emi;  c_max = pow(10.0, -exp_max); }
  };

  /** Localization of points with inversion of geometric transformations. \\ \\
   *  The points to be localized are sorted with "imbricated_box_less" order
   *  which is based in fact on a lexicographical order on the decimal
   *  decomposition of the coordinates. This order allows to find in a 
   *  fast manner which points are in a specified box.\\ \\
   *  To find which points are in a given element, the following algorithm is
   *  applied : \\
   *  - compute a englobing box of the element, assuming that this element,
   *    even if the geometric transfrmation in non-linear is in the englobing
   *    box of its nodes times a factor 1.2, \\
   *  - List all the points in this englobing box, \\
   *  - For each points, invert the geometric transformation (computation
   *    of a pseudo inverse in the linear case, and use of a Newton method in
   *    the non-linear case. \\ \\
   *  The inversion can be described as follows : \\ \\
   *  The geometric transformation from the reference element to the real
   *  element is described by a polynomial function \\
   *  $ \tau : {I\hspace{-0.3em}R}^P \longrightarrow {I\hspace{-0.3em}R}^N$ \\
   *  where $P$ is the dimension of the reference element and $N$ the dimension
   *  of the real element. This transformation is decomposed in \\
   *  $ \tau(\overline{x}) = \underline{\underline{a}} \ 
   *        \underline{\cal N}(\overline{x}). $ \\
   *  where $\underline{\underline{a}}$ is a $N \times n_g$ matrix representing
   *  the nodes of the real element and $\underline{\cal N}(\overline{x})$ is
   *  the vector of polynomials representing the geometric transformation. \\
   *  Then, the gradient $\nabla \tau(\overline{x})$ is the following
   *  $N \times P$ matrix: \\
   *  $ \nabla \tau(\overline{x}) = \underline{\underline{a}}\ \nabla
   *     \underline{\cal N}(\overline{x}). $\\
   *  \subsubsection*{Linear case}
   *  If $\tau$ is linear (affin in fact), then \\
   *  $y = \tau(\overline{x}) = \nabla \tau(0) \overline{x} + y_0.$ \\
   *  So that \\
   *  $\overline{x} = (\nabla \tau^T(0) \nabla \tau(0))^{-1} \nabla \tau^T(0) 
   *  (y - y_0), $ \\
   *  where \\
   *  $B(0) = (\nabla \tau^T(0) \nabla \tau(0))^{-1} \nabla \tau^T(0), $ \\
   *  is a pseudo inverse of $\nabla \tau(0)$. \\ \\
   *  if $N > P$, the residue \\
   *  $y - y_0 -  \nabla \tau(0) \overline{x}, $ \\
   *  indicates wether or not the point in on the "surface" of the convex. \\ 
   *  \\ \subsubsection*{Non-linear case}
   *  A Newton method is applied., writing \\
   *  $y = \tau(\overline{x}+\overline{h}) = \tau(\overline{x})
   *   + \nabla \tau(\overline{x})\overline{h} + o(\|\overline{h}\|^2). $ \\
   *  It gives the iterative scheme \\
   *  $ \overline{x}_{n+1} = \overline{x}_{n}
   *   + B(\overline{x}_{n})(y - \tau(\overline{x}_{n})), $  \\
   *  where \\
   *  $ B(\overline{x}_{n}) = (\nabla \tau^T(\overline{x}_{n})
   *   \nabla \tau(\overline{x}_{n}))^{-1} \nabla \tau^T(\overline{x}_{n}).$ \\
   *  The residue is \\
   *  $ y - \tau(\overline{x}_{n}).$ \\
   */
  class geotrans_inv
  {
    protected :

      typedef dal::dynamic_tree_sorted<base_node,imbricated_box_less> TAB_TYPE;
      TAB_TYPE ptab;
      scalar_type EPS;
    
    public :

      void clear(void) { ptab.clear(); }

      /// Add the points contained in c to the list of points.
      template<class CONT> void add_points(const CONT &c)
      {
	typename CONT::const_iterator it = c.begin(), ite = c.end();
	for (; it != ite; ++it) ptab.add(*it);
      }

      /// Add point p to the list of points.
      size_type add_point(base_node p) { return ptab.add(p); }
      
      /// Add point p to the list of points.
      size_type add_point_norepeat(base_node p)
      { return ptab.add_norepeat(p); }
      
      /// Find all the points present in the box between min and max.
      void points_in_box(dal::bit_vector &pt, const base_node &min,
			                      const base_node &max) const;

      /** Search all the points in the convex cv, which is the transformation
       *  of the convex cref via the geometric transformation pgt.
       *  IMPORTANT : It is assumed that the whole convex is include in the
       *     minmax box of its nodes times a factor 1.2. If the transformation
       *     is linear, the factor is 1.0.
       */
      template<class TAB, class CONT1, class CONT2>
         size_type points_in_convex(const convex<base_node, TAB> &cv,
			       pgeometric_trans pgt, CONT1 &pftab, CONT2 &itab)
      { // à optimiser
	size_type N = pgt->structure()->dim(); /* dimension of the convex.*/
	size_type P = cv.points()[0].size(); /* dimension of the image.     */
	base_node min(N), max(N);
	base_matrix a(P, pgt->nb_points());
	base_poly PO;
	base_node x, y;
	base_matrix pc(pgt->nb_points() , N);
	base_matrix grad(P, N), TMP1(N,N), B0(N, P), CS(N,N);
	size_type nbpt = 0, i;
	
	min = max = pgt->geometric_nodes()[0];
	for (size_type j = 0; j < pgt->nb_points(); ++j) // à optimiser !!
	  for (size_type i = 0; i < P; ++i)
	  { 
	    min[i] = std::min(min[i], cv.points()[j][i]);
	    max[i] = std::max(max[i], cv.points()[j][i]);
	    a(i,j) = cv.points()[j][i];
	  }
       
	dal::bit_vector pts;

	if (pgt->is_linear())
	{
	  for (i = 0; i < N; ++i) { min[i] -= EPS; max[i] += EPS; }
	  pts.clear();
	  // cout << "boxmin = " << min << " boxmax = " << max << endl;
	  
	  points_in_box(pts, min, max);
	  // pts= ptab.index();
	  
	  // On peut éviter ce calcul en faisant appel à un pre-geotrans
	  // ou en stockant le calcul qui est toujours le même.
	  // on peut aussi l'optimiser en ne faisant pas appel à derivative().
	  for (i = 0; i < pgt->nb_points(); ++i)
	    for (dim_type n = 0; n < N; ++n)
	    { PO = pgt->poly_vector()[i]; PO.derivative(n); pc(i,n) = PO[0]; }

	  // computation of the pseudo inverse
	  bgeot::mat_product(a, pc, grad);
	  if (N != P)
	  {
	    bgeot::mat_product_tn(grad, grad, CS);
	    bgeot::mat_inv_cholesky(CS, TMP1);
	    bgeot::mat_product_tt(CS, grad, B0);
	  }
	  else
	  {
	    // L'inversion peut être optimisée par le non calcul global de B0
	    // et la resolution d'un système linéaire.
	    bgeot::mat_gauss_inverse(grad, TMP1); B0 = grad;
	  }

	  // cout << "grad inverse : " << B0 << endl;
	  // cout << "points : " << pts << endl;

	  for (i << pts; i != size_type(-1); i << pts)
	  {
	    // cout << "point : " << ptab[i] << endl;
	    y = ptab[i]; y -= cv.points()[0];
	    x = B0 * y;
	    if (pgt->convex_ref()->is_in(x) < EPS)
	      if (N == P)
	      {
		// cout << "enregistré en " << nbpt << " : " << x << endl;
		pftab[nbpt] = x; itab[nbpt++] = i;
	      }
	      else
	      {
		y -= grad * x;
		if (vect_norm2(y) < EPS)
		  { pftab[nbpt] = x; itab[nbpt++] = i; }
	      }
	  }
	  // cout << "fini ... " << endl;

	}
	else
	{
	  for (i = 0; i < N; ++i)
	  { scalar_type e = (max[i]-min[i]) * 0.2;  min[i] -= e; max[i] += e; }
	  pts.clear();
	  points_in_box(pts, min, max);
     
	  base_node xn, rn;
	  scalar_type res;

	  for (i << pts; i != size_type(-1); i << pts)
	  {
	    x = pgt->geometric_nodes()[0]; y = cv.points()[0];  
	    scalar_type d = vect_dist2(y, ptab[i]);
	    for (size_type j = 1; j < pgt->nb_points(); ++j)
	    { 
	      scalar_type d2 = vect_dist2(cv.points()[j], ptab[i]);
	      if (d2 < d)
	      { d = d2; x = pgt->geometric_nodes()[j]; y = cv.points()[j]; }
	    }
	    
	    rn = ptab[i]; rn -= y; res = vect_norm2(rn);
	    while (res > EPS)
	    {
	      for (size_type k = 0; k < pgt->nb_points(); ++k)
		for (dim_type n = 0; n < N; ++n)
		{
		  PO = pgt->poly_vector()[k];
		  PO.derivative(n);
		  pc(k,n) = PO.eval(x.begin());
		}
	      
	      // computation of the pseudo inverse (it should be possible not
	      //  to compute it at each iteration).
	      bgeot::mat_product(a, pc, grad);
	      if (N != P)
	      {
		bgeot::mat_product_tn(grad, grad, CS);
		bgeot::mat_inv_cholesky(CS, TMP1);
		bgeot::mat_product_tt(CS, grad, B0);
	      }
	      else
	      {
		bgeot::mat_inv_cholesky(grad, TMP1); B0 = grad;
	      }
	      xn = x;
	      x = B0 * rn;
	      y.fill(0.0);
	      for (size_type k = 0; k < pgt->nb_points(); ++k)
		y.addmul(pgt->poly_vector()[k].eval(x.begin()),
			  cv.points()[k]);
	      rn = ptab[i]; rn -= y; res = vect_dist2(x, xn);
	    }
	    if (pgt->convex_ref()->is_in(x) < EPS
		&& (N == P || vect_norm2(rn) < EPS))
	    { pftab[nbpt] = x; itab[nbpt++] = i; }
	    // Test un peu sevère peut-être en ce qui concerne rn.
	  }
	}
	return nbpt;
      }

      
      geotrans_inv(void) { EPS = 10E-12; }
  };
  
  

}  /* end of namespace bgeot.                                             */


#endif /* __BGEOT_GEOMETRIC_TRANSFORMATION_H                              */
