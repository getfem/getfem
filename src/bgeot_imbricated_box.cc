//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_imbricated_box.cc : particular point sort.
//           
// Date    : January 26, 1999.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2005-2005 Yves Renard
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

#include <bgeot_imbricated_box.h>

namespace bgeot {
#if 0
#define PBASE2 4
  static double frexpb2(double a, int& e) {
    double x = frexp(a,&e);
    if (PBASE2 != 1) {
      int n = e; e /= PBASE2;
      x *= (1<<(n - e*PBASE2));
    }
    return x;
  }
  int imbricated_box_less::operator()(const base_node &A,
				      const base_node &B) const 
  {
    base_vector x1(A.size()), x2(A.size()); 
    std::vector<int> e1(A.size()), e2(A.size());
    int emax = -1000, emax2;
    base_node::const_iterator a=A.begin(), b=B.begin();
    for (size_type i=0; i < A.size(); ++i) {
      x1[i] = frexpb2(a[i], e1[i]); x2[i] = frexpb2(b[i], e2[i]);
      emax = std::max(emax, e1[i]);
      emax = std::max(emax, e2[i]);
    }
    emax2 = emax;
    for (size_type i=0; i < A.size(); ++i) {
      if (x1[i] < 0 || x2[i] < 0) {
	x1[i] = frexpb2(a[i]+exp2(emax2), e1[i]); x2[i] = frexpb2(b[i]+exp2(emax2), e2[i]);
	//assert(x1[i]>=0); assert(x2[i]>=0);
	emax = std::max(emax, e1[i]);
	emax = std::max(emax, e2[i]);
      }
    }
    bool finished;  
    do {
      finished = true;
      for (size_type i=0; i < A.size(); ++i) {
	int d1,d2;
	if (e1[i] != emax && e2[i] != emax) continue;
	if (e1[i] == emax) { x1[i] *= PBASE2; d1 = (int)x1[i]; x1[i] -= d1; e1[i]--; } else d1 = 0;
	if (e2[i] == emax) { x2[i] *= PBASE2; d2 = (int)x2[i]; x2[i] -= d2; e2[i]--; } else d2 = 0;
	if (d1 != d2) { return (d1 < d2) ? -1 : +1; }
	finished = finished && !(e1[i] != e2[i] || x1[i] != x2[i]);
      }
      emax--;
    } while (!finished);
    return 0;
  }
#else
  int imbricated_box_less::operator()(const base_node &x,
				      const base_node &y) const 
  {
    size_type s = x.size(); 
    scalar_type c1 = c_max, c2 = c_max * scalar_type(base);
    if (y.size() != s) DAL_THROW(dimension_error, "dimension error");
    
    base_node::const_iterator itx=x.begin(), itex=x.end(), ity=y.begin();
    int ret = 0;
    for (; itx != itex; ++itx, ++ity) {
      long a = long(sfloor((*itx) * c1)), b = long(sfloor((*ity) * c1));
      if ((dal::abs(a) > scalar_type(base))
	  || (dal::abs(b) > scalar_type(base))) { 
	exp_max++; c_max /= scalar_type(base);
	return (*this)(x,y);
      }
      if (ret == 0) { if (a < b) ret = -1; else if (a > b) ret = 1; }
    }
    if (ret) return ret;
    
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
#endif

#if 0 // disable all the rest of the file

#if 1
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
#elif 0
  size_type geotrans_inv::points_in_box(dal::dynamic_array<size_type> &pt,
					const base_node &min,
					const base_node &max) const {
    /* The following is a version with a partition, avoiding default */
    /* of the simple search, but which is slower .. in the mean.     */
    size_type s = min.size(), i, nbib = 0; 
    base_node c(s),boxmin(s),boxmax(s),cbox(s), iboxmin(s), iboxmax(s); 
    TAB_TYPE::const_sorted_iterator it, ite; 
    scalar_type logbase = log(double(base()));
    cout.precision(25);
    cout << "initial box : " << min << " :: " << max << endl;
    
    for (i = 0; i < s; ++i)
      {
	c[i] = pow(base(),
		   rint(log(std::max(EPS, max[i] - min[i]))/logbase));
	boxmin[i] = floor(min[i] / c[i]) * c[i];
	boxmax[i] = ceil(max[i] / c[i]) * c[i];
      }
    cout << "max box : " << boxmin << " :: " << boxmax << endl;
    cout << "steps : " << c << endl;
    
    cbox = boxmin;
    while(cbox[s-1] < boxmax[s-1]-EPS)
      {
	/* intersection */
	for (i = 0; i < s; ++i)
    	  {
	    if (cbox[i] > max[i] || cbox[i]+c[i] < min[i]) goto aurevoir;
    	    iboxmin[i]=std::max(cbox[i], min[i])+EPS;
    	    iboxmax[i]=std::max(std::min(cbox[i]+c[i], max[i]), iboxmin[i])-EPS;
    	  }
	cout << "intersection : " << iboxmin << " : " << iboxmax << endl;
	/* recherche des points entre iboxmin et iboxmax */
	it = ptab.sorted_ge(iboxmin);
	ite = ptab.sorted_ge(iboxmax);
	cout << "ite-it=" << std::distance(it,ite) << "\n";
	cout << "pt " << ptab[0] << " < " << iboxmin << "? : " << ptab.comparator()(ptab[0], iboxmin) << "\n";
	cout << "pt " << ptab[0] << " < " << iboxmax << "? : " << ptab.comparator()(ptab[0], iboxmax) << "\n";
	for(; it != ite; ++it)
    	  {
    	    bool isin = true;
    	    for (i = 0; i < s; ++i)
    	      if ((*it)[i] < min[i] || (*it)[i] > max[i])
    		{ isin = false; break; }
    	    if (isin) pt[nbib++] = it.index();
    	  }
	cout << "nbib = " << nbib << "\n";
      aurevoir:
	/* incrementation */
	i = 0; cbox[0] += c[0];
	while((cbox[i] >= boxmax[i]-EPS) && (i < s-1))
    	  { cbox[i] = boxmin[i]; ++i; cbox[i] += c[i]; }
      }
    return nbib;
  }
#else
  size_type geotrans_inv::points_in_box(dal::dynamic_array<size_type> &pt,
					const base_node &Pmin,
					const base_node &Pmax) const {
    size_type N = Pmin.size(), nbib=0;
    int boxsz = -1000;
    base_node boxmin(N), boxmax(N);
    TAB_TYPE::const_sorted_iterator it, ite; 

    base_node::const_iterator min=Pmin.begin(), max=Pmax.begin();
    for (size_type i=0; i < N; ++i) {
      int j; frexp(max[i]-min[i], &j);
      boxsz = std::max(boxsz,j);
    }
    for (size_type i=0; i < N; ++i) {
      scalar_type e = exp2(boxsz*2);
      boxmin[i] = floor(min[i]/e)*e; 
      boxmax[i] = ceil(max[i]/e)*e - 1e-10;
    }
    cout.precision(20);
    cout << "points_in_box(" << Pmin << ", " << Pmax << "): boxsz=2^" << boxsz << "=" << exp2(boxsz) << ", boxmin=" << boxmin << ", boxmax=" << boxmax << "\n";
    it = ptab.sorted_ge(boxmin);
    ite = ptab.sorted_ge(boxmax);
    cout << "ite-it=" << std::distance(it,ite) << "\n";
    /*cout << "pt " << ptab[0] << " < " << iboxmin << "? : " << ptab.comparator()(ptab[0], iboxmin) << "\n";
      cout << "pt " << ptab[0] << " < " << iboxmax << "? : " << ptab.comparator()(ptab[0], iboxmax) << "\n";*/
    for(; it != ite; ++it) {
      bool ge_boxmin = true, le_boxmax = true;
      for (size_type i = 0; i < N; ++i) {
	if ((*it)[i]<boxmin[i]) ge_boxmin = false;
	if ((*it)[i]>boxmax[i]) le_boxmax = false;
      }
      cout << "test point " << *it << ", ge_boxmin=" << ge_boxmin << ", le_boxmax=" << le_boxmax << "\n";
      assert(ge_boxmin && le_boxmax);
      bool isin = true;
      for (size_type i = 0; i < N; ++i) {
	if ((*it)[i] < min[i] || (*it)[i] > max[i]) {
	  isin = false; break; 
	}
      }
      if (isin) pt[nbib++] = it.index();
    }
    cout << "nbib = " << nbib << "\n";
    return nbib;
  }
#endif   
#endif
}
