/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_geotrans_inv.C : Allow to inverse geometric transf-    */
/*     	      ormations and to localize a set of points.    	       	   */
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



#include <bgeot_geotrans_inv.h>

namespace bgeot
{ 
  /* inversion for linear geometric transformations */
  bool geotrans_inv_convex::invert_lin(const base_node& n, base_node& n_ref, scalar_type IN_EPS) {
    base_node y(n); for (size_type i=0; i < P; ++i) y[i] -= a(i,0);
    gmm::mult(B0, y, n_ref); // n_ref = B0 * y;
    if (pgt->convex_ref()->is_in(n_ref) < IN_EPS) {
      if (N == P) return true;
      else {
	gmm::mult(grad,gmm::scaled(n_ref,-1.0),y,y);
	//        y -= grad * n_ref;
        if (vect_norm2(y) < IN_EPS) return true;
      }
    }
    return false;
  }
  
  /* inversion for non-linear geometric transformations */
  bool geotrans_inv_convex::invert_nonlin(const base_node& n, base_node& x, scalar_type IN_EPS) {
    base_node xn, y;
    //cout << "invert_nonlin(" << n << "," << x << "," << IN_EPS << ")\n";
    /* find an initial guess */
    x = pgt->geometric_nodes()[0]; y = cvpts[0];  
    scalar_type d = vect_dist2_sqr(y, n);
    for (size_type j = 1; j < pgt->nb_points(); ++j) { 
      scalar_type d2 = vect_dist2_sqr(cvpts[j], n);
      if (d2 < d)
        { d = d2; x = pgt->geometric_nodes()[j]; y = cvpts[j]; }
    }
    base_node rn(n); rn -= y; 
    scalar_type res = vect_norm2(rn);
    unsigned cnt = 1000;
    while (res > EPS && --cnt) {
      if (cnt < 30) cout << "res=" << res << ", x=" << x << "\n";
      /* compute gradient */
      for (size_type k = 0; k < pgt->nb_points(); ++k) {
        for (dim_type nn = 0; nn < N; ++nn) {
          PO = pgt->poly_vector()[k];
          PO.derivative(nn);
          pc(k,nn) = PO.eval(x.begin());
        }
      }
      
      // computation of the pseudo inverse (it should be possible not
      //  to compute it at each iteration).
      gmm::mult(a, pc, grad);
      if (N != P) {
	gmm::mult(gmm::transposed(grad), grad, CS);
        gmm::lu_inverse(CS);
	gmm::mult(gmm::transposed(CS), gmm::transposed(grad), B0);
      } else {
        gmm::lu_inverse(grad); B0 = grad;
      }
      xn = x;
      gmm::mult(B0, rn, x); // x = B0 * rn;
      x += xn;
      y.fill(0.0);
      for (size_type k = 0; k < pgt->nb_points(); ++k) {
	gmm::add(gmm::scaled(cvpts[k],pgt->poly_vector()[k].eval(x.begin())),y);
      /*y.addmul(pgt->poly_vector()[k].eval(x.begin()), cvpts[k]);*/
      }
      // cout << "Point : " << x << " : " << y << " ptab : " << ptab[i] << endl; getchar();
      rn = n; rn -= y; res = vect_dist2(x, xn);
    }
    //cout << " invert_nonlin done\n";
    if (cnt == 0) 
      DAL_THROW(dal::failure_error, 
                "inversion of non-linear geometric transformation "
                "failed (too much iterations)");
    
    // Test un peu sevère peut-être en ce qui concerne rn.
    if (pgt->convex_ref()->is_in(x) < IN_EPS
        && (N == P || vect_norm2(rn) < IN_EPS)) {
      //cout << "point " << x << "in IN (" << pgt->convex_ref()->is_in(x) << ")\n";
      return true;
    } //else cout << "point IS OUT\n";
    return false;
  }

}  /* end of namespace bgeot.                                             */
