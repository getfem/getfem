/* *********************************************************************** */
/*                                                                         */
/* Library : GEneric Tool for Finite Element Methods (getfem)              */
/* File    : getfem_poly_composite.C : polynomials by parts                */
/*                                                                         */
/* Date : August 26, 2002.                                                 */
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


#include <getfem_poly_composite.h>

namespace getfem
{ 
  mesh_precomposite::mesh_precomposite(const getfem_mesh &m) {
    mesh = &m;
    elt.resize(m.nb_convex());
    det.resize(m.nb_convex());
    orgs.resize(m.nb_convex());
    gtrans.resize(m.nb_convex());
    dal::bit_vector nn = m.points().index();
    for (size_type i = 0; i <= nn.last_true(); ++i) {
      vertexes.add(m.points()[i]);
    }
    nn = m.convex_index();
    for (size_type cv = nn.take_first(); cv != size_type(-1); cv << nn) {
      
      bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
      size_type N = pgt->structure()->dim();
      size_type P = m.dim();
      if (!(pgt->is_linear()) || N != P) 
	DAL_THROW(internal_error, "Bad geometric transformations.");
    
      base_poly PO;
      base_matrix a(P, pgt->nb_points());
      base_matrix pc(pgt->nb_points() , N);
      base_matrix TMP1(N,N), B0(N, P);
    
      for (size_type j = 0; j < pgt->nb_points(); ++j)
	for (size_type k = 0; k < P; ++k)
	  a(k,j) = (m.points_of_convex(cv)[j])[k];
      
      for (size_type k = 0; k < pgt->nb_points(); ++k)
	for (dim_type n = 0; n < N; ++n)
	  { PO = pgt->poly_vector()[k]; PO.derivative(n); pc(k,n) = PO[0]; }
    
      bgeot::mat_product(a, pc, B0);
      bgeot::mat_gauss_inverse(B0, TMP1);
      det[cv] = 1.0 / bgeot::mat_gauss_det(B0, TMP1);
      gtrans[cv] = B0;
      orgs[cv] = m.points_of_convex(cv)[0];
    
    }
  }

  void polynomial_composite::derivative(short_type k) {
    base_poly P, Q;
    for (size_type ic = 0; ic < mp->nb_convex(); ++ic) {
      base_vector e(mp->dim());
      e.fill(0.0); e[k] = 1.0;
      e *= mp->gtrans[ic]; // vérifier si ce n'est pas la transposée
      P.clear();
      for (dim_type n = 0; n < mp->dim(); ++n)
	{ Q = polytab[ic]; Q.derivative(n); P += Q * e[n]; }
      polytab[ic] = P;
    }
  }


  struct __str_mesh_cv {
    bgeot::pconvex_structure cvs;
    short_type n;
    getfem_mesh *pm;
    mesh_precomposite *pmp;
    bool operator < (const __str_mesh_cv &ls) const {
      if (cvs < ls.cvs) return true; if (cvs > ls.cvs) return false; 
      if (n < ls.n) return true; return false;
    }
    __str_mesh_cv(void) {}
    __str_mesh_cv(bgeot::pconvex_structure c, short_type k) : cvs(c), n(k) {}
  };

  static dal::dynamic_tree_sorted<__str_mesh_cv> *__str_mesh_cv_tab = 0;


  void structured_mesh_for_convex(bgeot::pconvex_ref cvr, short_type k,
				  pgetfem_mesh &pm, pmesh_precomposite &pmp) {

    bool found = false;
    size_type n = cvr->structure()->dim();
    size_type nbp = cvr->structure()->basic_structure()->nb_points();
    
    if (__str_mesh_cv_tab == 0)
      __str_mesh_cv_tab = new dal::dynamic_tree_sorted<__str_mesh_cv>();

    __str_mesh_cv smc(cvr->structure()->basic_structure(), k);
    size_type iss = __str_mesh_cv_tab->search(smc);
    if (iss == size_type(-1)) {

      /* Identifying simplexes.                                           */

      if (nbp == n+1)
	if (cvr->structure()->basic_structure()==bgeot::simplex_structure(n)) {
	  smc.pm = new getfem_mesh();
	  switch (n) {
	  case 1 :
	    {
	      base_vector a(1), b(1);
	      for (short_type i = 0; i < k; ++i) {
		a[0] = i * (1.0 / k); b[0] = (i+1) * (1.0 / k);
		smc.pm->add_segment_by_points(a, b);
	      }
	    }
	    break;
	  case 2 :
	    for (short_type i = 0; i < k; ++i) {
	      scalar_type a = i * (1.0 / k), b = (i+1) * (1.0 / k);
	      for (short_type l = 0; l+i < k; ++l) {
		scalar_type c = l * (1.0 / k), d = (l+1) * (1.0 / k);
		smc.pm->add_triangle_by_points
		  (base_vector(a, c),
		   base_vector(a, d),
		   base_vector(b, ((l+i+k) & 1) ? c : d));
		if (l+i+1 < k)
		  smc.pm->add_triangle_by_points
		    (base_vector(b, c),
		     base_vector(b, d),
		     base_vector(a, ((l+i+k) & 1) ? d : c));
	      }
	    }
	    break;
	  default : 
	    delete smc.pm;
	    DAL_THROW(to_be_done_error, "Sorry, not implemented.");
	  }
	  // smc.pm->write_to_file(cout);

	  smc.pmp = new mesh_precomposite(*(smc.pm));
	  iss = __str_mesh_cv_tab->add(smc);
	  found = true;
	}
      
      /* Identifying parallelepipeds.                                     */

      if (!found && nbp == (size_type(1) << n))
	if (cvr->structure()->basic_structure()
	    == bgeot::parallelepiped_structure(n)) {
	  DAL_THROW(to_be_done_error, "Sorry, not implemented.");
	  found = true;
	}
      
      /* Identifying prisms.                                              */
      
      if (!found && nbp == 2 * n)
	if (cvr->structure()->basic_structure() == bgeot::prism_structure(n)) {
	  DAL_THROW(to_be_done_error, "Sorry, not implemented.");
	  found = true;
	}
    }

    if (iss == size_type(-1))
      DAL_THROW(to_be_done_error,
		"This element is not taken into account. Contact us");
    

   
    pm  = (*__str_mesh_cv_tab)[iss].pm;
    pmp = (*__str_mesh_cv_tab)[iss].pmp;

  }
  
}  /* end of namespace getfem.                                            */
