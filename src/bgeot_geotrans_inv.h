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
/* Copyright (C) 2000-2001  Yves Renard.                                   */
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

#ifndef __BGEOT_GEOTRANS_INV_H
#define __BGEOT_GEOTRANS_INV_H

#include <bgeot_geometric_trans.h>

namespace bgeot
{
 
  inline scalar_type sfloor(scalar_type x)
  { return (x >= 0) ? floor(x) : -floor(-x); }

  /// Sort order for a rapid localization in bgeot::geotrans_inv
  struct imbricated_box_less
    : public std::binary_function<base_node, base_node, int>
  { 
    mutable int exp_max, exp_min;
    mutable scalar_type c_max;
    int base;

    /// comparaison function
    int operator()(const base_node &x, const base_node &y) const;
    
    imbricated_box_less(int ba = 10, int emi = -15, int ema = -2) {
      base = ba; exp_max = ema; exp_min = emi;
      c_max = pow(double(base), double(-exp_max));
    }
  };

  class geotrans_inv
  {
    protected :

      typedef dal::dynamic_tree_sorted<base_node,imbricated_box_less> TAB_TYPE;
      TAB_TYPE ptab;
      scalar_type EPS;
    
    public :

      void clear(void) { ptab.clear(); }

      /// Add the points contained in c to the list of points.
      template<class CONT> void add_points(const CONT &c) {
	typename CONT::const_iterator it = c.begin(), ite = c.end();
	for (; it != ite; ++it) ptab.add(*it);
      }

      /// Add point p to the list of points.
      size_type add_point(base_node p) { return ptab.add(p); }
      
      /// Add point p to the list of points.
      size_type add_point_norepeat(base_node p)
      { return ptab.add_norepeat(p); }
      
      /// Find all the points present in the box between min and max.
      size_type points_in_box(dal::dynamic_array<size_type> &pt,
			      const base_node &min, 
			      const base_node &max) const;

      /** Search all the points in the convex cv, which is the transformation
       *  of the convex cref via the geometric transformation pgt.
       *  IMPORTANT : It is assumed that the whole convex is include in the
       *     minmax box of its nodes times a factor 1.2. If the transformation
       *     is linear, the factor is 1.0.
       */
      template<class TAB, class CONT1, class CONT2>
      size_type points_in_convex(const convex<base_node, TAB> &cv,
				 pgeometric_trans pgt,
				 CONT1 &pftab, CONT2 &itab);

      
      geotrans_inv(int ba = 10) : ptab(imbricated_box_less(ba)) { EPS = 10E-12; }
  };



  template<class TAB, class CONT1, class CONT2>
  size_type geotrans_inv::points_in_convex(const convex<base_node, TAB> &cv,
					   pgeometric_trans pgt,
					   CONT1 &pftab, CONT2 &itab) {
    size_type N = pgt->structure()->dim(); /* dimension of the convex.*/
    size_type P = cv.points()[0].size(); /* dimension of the image.     */
    base_node min(P), max(P);
    base_matrix a(P, pgt->nb_points());
    base_poly PO;
    base_node x(N), y(P);
    base_matrix pc(pgt->nb_points() , N);
    base_matrix grad(P, N), TMP1(N,N), B0(N, P), CS(N,N);
    size_type nbpt = 0;
    dal::dynamic_array<size_type> pts;

    min = max = cv.points()[0];
    for (size_type j = 0; j < pgt->nb_points(); ++j) // à optimiser !!
      for (size_type i = 0; i < P; ++i) { 
	min[i] = std::min(min[i], cv.points()[j][i]);
	max[i] = std::max(max[i], cv.points()[j][i]);
	a(i,j) = cv.points()[j][i];
      }
    
    if (pgt->is_linear()) {
      for (size_type i = 0; i < N; ++i) { min[i] -= EPS; max[i] += EPS; }
      // cout << "boxmin = " << min << " boxmax = " << max << endl;
      
      size_type nbib = points_in_box(pts, min, max);
      // pts= ptab.index();
      
      // On peut éviter ce calcul en faisant appel à un pre-geotrans
      // ou en stockant le calcul qui est toujours le même.
      // on peut aussi l'optimiser en ne faisant pas appel à derivative().
      for (size_type i = 0; i < pgt->nb_points(); ++i)
	for (dim_type n = 0; n < N; ++n)
	  { PO = pgt->poly_vector()[i]; PO.derivative(n); pc(i,n) = PO[0]; }
      
      // computation of the pseudo inverse
      bgeot::mat_product(a, pc, grad);
      if (N != P) {
	bgeot::mat_product_tn(grad, grad, CS);
	bgeot::mat_inv_cholesky(CS, TMP1);
	bgeot::mat_product_tt(CS, grad, B0);
      }
      else {
	// L'inversion peut être optimisée par le non calcul global de B0
	// et la resolution d'un système linéaire.
	bgeot::mat_gauss_inverse(grad, TMP1); B0 = grad;
      }
      
      for (size_type l = 0; l < nbib; ++l) {
	// cout << "point : " << ptab[i] << endl;
	y = ptab[pts[l]]; y -= cv.points()[0];
	mat_vect_product(B0, y, x); // x = B0 * y;
	if (pgt->convex_ref()->is_in(x) < EPS)
	  if (N == P) {
	    // cout << "enregistré en " << nbpt << " : " << x << endl;
	    pftab[nbpt] = x; itab[nbpt++] = pts[l];
	  }
	  else {
	    y -= grad * x;
	    if (vect_norm2(y) < EPS)
	      { pftab[nbpt] = x; itab[nbpt++] = pts[l]; }
	  }
      }
      // cout << "fini ... " << endl;
      
    }
    else { // partie non testée
      for (size_type i = 0; i < N; ++i)
	{ scalar_type e = (max[i]-min[i]) * 0.2;  min[i] -= e; max[i] += e; }
      size_type nbib = points_in_box(pts, min, max);
      
      base_node xn, rn;
      scalar_type res;

      for (size_type l = 0; l < nbib; ++l) {
	size_type i = pts[l];
	x = pgt->geometric_nodes()[0]; y = cv.points()[0];  
	scalar_type d = vect_dist2(y, ptab[i]);
	for (size_type j = 1; j < pgt->nb_points(); ++j) { 
	  scalar_type d2 = vect_dist2(cv.points()[j], ptab[i]);
	  if (d2 < d)
	    { d = d2; x = pgt->geometric_nodes()[j]; y = cv.points()[j]; }
	}
	
	rn = ptab[i]; rn -= y; res = vect_norm2(rn);
	while (res > EPS) {
	  for (size_type k = 0; k < pgt->nb_points(); ++k)
	    for (dim_type n = 0; n < N; ++n) {
	      PO = pgt->poly_vector()[k];
	      PO.derivative(n);
	      pc(k,n) = PO.eval(x.begin());
	    }
	  
	  // computation of the pseudo inverse (it should be possible not
	  //  to compute it at each iteration).
	  bgeot::mat_product(a, pc, grad);
	  if (N != P) {
	    bgeot::mat_product_tn(grad, grad, CS);
	    bgeot::mat_inv_cholesky(CS, TMP1);
	    bgeot::mat_product_tt(CS, grad, B0);
	  }
	  else {
	    bgeot::mat_gauss_inverse(grad, TMP1); B0 = grad;
	  }
	  // cout << "grad = " << grad << endl;
	  xn = x;
	  mat_vect_product(B0, rn, x); // x = B0 * rn;
	  x += xn;
	  y.fill(0.0);
	  for (size_type k = 0; k < pgt->nb_points(); ++k)
	    y.addmul(pgt->poly_vector()[k].eval(x.begin()),
		     cv.points()[k]);
	  // cout << "Point : " << x << " : " << y << " ptab : " << ptab[i] << endl; getchar();
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
  

}  /* end of namespace bgeot.                                             */


#endif /* __BGEOT_GEOMETRIC_TRANSFORMATION_H                              */
