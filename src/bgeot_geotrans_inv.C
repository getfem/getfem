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
  int imbricated_box_less::operator()(const base_node &x,
				      const base_node &y) const { 
    size_type s = x.size(); 
    scalar_type c1 = c_max, c2 = c_max * scalar_type(base);
    if (y.size() != s) DAL_THROW(dimension_error, "dimension error");
    
    base_node::const_iterator itx=x.begin(), itex=x.end(), ity=y.begin();
    for (; itx != itex; ++itx, ++ity) {
      long a = long(sfloor((*itx) * c1)), b = long(sfloor((*ity) * c1));
      if ((dal::abs(a) > scalar_type(base))
	  || (dal::abs(b) > scalar_type(base))) { 
	exp_max++; c_max /= scalar_type(base);
	return (*this)(x,y);
      }
      if (a < b) return -1; else if (a > b) return 1;
    }
    
    for (int e = exp_max; e >= exp_min; --e, c1 *= scalar_type(base),
	   c2 *= scalar_type(base)) {
      itx = x.begin(), itex = x.end(), ity = y.begin();
      for (; itx != itex; ++itx, ++ity) {
	int a = int(sfloor(((*itx) * c2) - sfloor((*itx) * c1)
			   * scalar_type(base)));
	int b = int(sfloor(((*ity) * c2) - sfloor((*ity) * c1)
			   * scalar_type(base)));
	if (a < b) return -1; else if (a > b) return 1;
      }
    }
    return 0;
  }
  
  size_type geotrans_inv::points_in_box(dal::dynamic_array<size_type> &pt,
					const base_node &min,
					const base_node &max) const {
    TAB_TYPE::const_sorted_iterator it, ite;
    size_type nb = 0;
    
    it = ptab.sorted_ge(min); ite = ptab.sorted_ge(max);
    base_node::const_iterator itl, itmin, itmax, itmine = min.end();
    for(; it != ite; ++it) { 
      bool isin = true;
      itl = (*it).begin(); itmin = min.begin(); itmax = max.begin();
      for (; itmin != itmine; ++itmin, ++itmax, ++itl)
	if (*itl < *itmin || *itl > *itmax) { isin = false; break; }
      if (isin) pt[nb++] = it.index();
    }
    return nb;
  }    
    /* The following is a version with a partition, avoiding default */
    /* of the simple search, but which is slower .. in the mean.     */
    /* 	size_type s = min.size(), i; */
    /* 	base_node c(s),boxmin(s),boxmax(s),cbox(s), iboxmin(s), iboxmax(s); */
    /* 	TAB_TYPE::const_sorted_iterator it, ite; */
    
    /*  cout << "initial box : " << min << " :: " << max << endl; */
    
    /* 	for (i = 0; i < s; ++i) */
    /* 	{ */
    /* 	  c[i] = pow(10.0, */
    /* 		 rint(std::log10(std::max(EPS, max[i] - min[i])))); */
    /* 	  boxmin[i] = floor(min[i] / c[i]) * c[i]; */
    /* 	  boxmax[i] = ceil(max[i] / c[i]) * c[i]; */
    /* 	} */
    /*  cout << "max box : " << boxmin << " :: " << boxmax << endl; */
    /*  cout << "steps : " << c << endl; */
    
    /* 	cbox = boxmin; */
    /* 	while(cbox[s-1] < boxmax[s-1]-EPS) */
    /* 	{ */
    /*  intersection */
    /* 	  for (i = 0; i < s; ++i) */
    /* 	  {  */
    /* 	    iboxmin[i]=std::max(cbox[i], min[i]); */
    /* 	    iboxmax[i]=std::max(std::min(cbox[i]+c[i], max[i]), iboxmin[i]); */
    /* 	  } */
    /* 	   recherche des points entre iboxmin et iboxmax */
    /* 	  it = ptab.sorted_ge(&iboxmin); */
    /* 	  ite = ptab.sorted_ge(&iboxmax); */
    /* 	  for(; it != ite; ++it) */
    /* 	  {  */
    /* 	    bool isin = true; */
    /* 	    for (i = 0; i < s; ++i) */
    /* 	      if ((*(*it))[i] < min[i] || (*(*it))[i] > max[i]) */
    /* 		{ isin = false; break; } */
    /* 	    if (isin) pt.add(it.index()); */
    /* 	  } */
    
    /* 	  incrementation */
    /* 	  i = 0; cbox[0] += c[0]; */
    /* 	  while((cbox[i] >= boxmax[i]-EPS) && (i < s-1)) */
    /* 	  { cbox[i] = boxmin[i]; ++i; cbox[i] += c[i]; } */
    /* 	} */
    


  
  bool geotrans_inv_convex::invert_lin(const base_node& n, base_node& n_ref, scalar_type IN_EPS) {
    base_node y(n); for (size_type i=0; i < P; ++i) y[i] -= a(i,0);
    gmm::mult(B0, y, n_ref); // n_ref = B0 * y;
    if (pgt->convex_ref()->is_in(n_ref) < IN_EPS) {
      if (N == P) return true;
      else {
	gmm::mult(grad,gmm::scaled(base_vector(n_ref),-1.0),y,y);
	//        y -= grad * n_ref;
        if (vect_norm2(y) < IN_EPS) return true;
      }
    }
    return false;
  }

  bool geotrans_inv_convex::invert_nonlin(const base_node& n, base_node& x, scalar_type IN_EPS) {
    base_node xn, y;
    
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
      for (size_type k = 0; k < pgt->nb_points(); ++k)
        y.addmul(pgt->poly_vector()[k].eval(x.begin()),
                 cvpts[k]);
      // cout << "Point : " << x << " : " << y << " ptab : " << ptab[i] << endl; getchar();
      rn = n; rn -= y; res = vect_dist2(x, xn);
    }
    if (cnt == 0) 
      DAL_THROW(dal::failure_error, 
                "inversion of non-linear geometric transformation "
                "failed (too much iterations)");
    
    // Test un peu sevère peut-être en ce qui concerne rn.
    if (pgt->convex_ref()->is_in(x) < IN_EPS
        && (N == P || vect_norm2(rn) < IN_EPS))
      return true;
    return false;
  }

}  /* end of namespace bgeot.                                             */
