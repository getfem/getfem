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
/* Copyright (C) 2002  Yves Renard.                                        */
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

#include <bgeot_comma_init.h>
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
      base_matrix B0(N, P);
    
      for (size_type j = 0; j < pgt->nb_points(); ++j)
	for (size_type k = 0; k < P; ++k)
	  a(k,j) = (m.points_of_convex(cv)[j])[k];
      
      for (size_type k = 0; k < pgt->nb_points(); ++k)
	for (dim_type n = 0; n < N; ++n)
	  { PO = pgt->poly_vector()[k]; PO.derivative(n); pc(k,n) = PO[0]; }
    
      gmm::mult(gmm::transposed(pc), gmm::transposed(a), B0);
      det[cv] = bgeot::mat_inverse(B0);
      gtrans[cv] = B0;
      orgs[cv] = m.points_of_convex(cv)[0];
    
    }
  }

  scalar_type polynomial_composite::eval(const base_node &pt) const {
    base_node p0(mp->dim()), p1(mp->dim());
    std::fill(mp->elt.begin(), mp->elt.end(), true);
    bgeot::mesh_convex_ind_ct::const_iterator itc, itce;
    
    mesh_precomposite::PTAB::const_sorted_iterator
      it1 = mp->vertexes.sorted_ge(pt), it2 = it1;    
    size_type i1 = it1.index(), i2;

    --it2; i2 = it2.index();
    

    while (i1 != size_type(-1) || i2 != size_type(-1)) {
      if (i1 != size_type(-1)) {
	bgeot::mesh_convex_ind_ct tc = mp->linked_mesh().convex_to_point(i1);
	itc = tc.begin(); itce = tc.end();
	for (; itc != itce; ++itc) {
	  size_type ii = *itc;
	  if (mp->elt[ii]) {
	    mp->elt[ii] = false;
	    p0 = pt; p0 -= mp->orgs[ii];
	    gmm::mult(gmm::transposed(mp->gtrans[ii]), p0, p1);
	    if (mp->trans_of_convex(ii)->convex_ref()->is_in(p1) < 1E-10)
	      return  polytab[ii].eval(p1.begin());
	  }
	}
	++it1; i1 = it1.index();
      }
      if (i2 != size_type(-1)) {
	bgeot::mesh_convex_ind_ct tc = mp->linked_mesh().convex_to_point(i2);
	itc = tc.begin(); itce = tc.end();
	for (; itc != itce; ++itc) {
	  size_type ii = *itc;
	  if (mp->elt[ii]) {
	    mp->elt[ii] = false;
	    p0 = pt; p0 -= mp->orgs[ii];
	    gmm::mult(gmm::transposed(mp->gtrans[ii]), p0, p1);
	    if (mp->trans_of_convex(ii)->convex_ref()->is_in(p1) < 1E-10)
	      return  polytab[ii].eval(p1.begin());
	  }
	}
	--it2; i2 = it2.index();
      }
    }
    DAL_THROW(internal_error, "Element not found in composite polynomial: " << pt);
  }


  polynomial_composite::polynomial_composite(const mesh_precomposite &m)
      : mp(&m), polytab(m.nb_convex()) {
    std::fill(polytab.begin(), polytab.end(), base_poly(m.dim(), 0));
  }

  void polynomial_composite::derivative(short_type k) {
    dim_type N = mp->dim();
    base_poly P(N, 0), Q;
    base_vector e(N), f(N);
    for (size_type ic = 0; ic < mp->nb_convex(); ++ic) {
      e.fill(0.0); e[k] = 1.0;
      gmm::mult(gmm::transposed(mp->gtrans[ic]), e, f);
      P.clear();
      for (dim_type n = 0; n < N; ++n)
	{ Q = polytab[ic];
	Q.derivative(n);
	P += Q * f[n];  }
      polytab[ic] = P;
    }
  }



  static void _structured_mesh_for_simplex(bgeot::pconvex_structure cvs, 
					   bgeot::pgeometric_trans opt_gt, const std::vector<base_node> *opt_gt_pts,
					   short_type k, pgetfem_mesh &pm) {
    scalar_type h = 1./k;
    switch (cvs->dim()) {
    case 1 :
      {
	base_node a(1), b(1);
	for (short_type i = 0; i < k; ++i) {
	  a[0] = i * h; b[0] = (i+1) * h;
	  if (opt_gt) a = opt_gt->transform(a, *opt_gt_pts);
	  if (opt_gt) b = opt_gt->transform(b, *opt_gt_pts);
	  pm->add_segment_by_points(a, b);
	}
      }
      break;
    case 2 :
      {
	base_node A,B,C,D;
	for (short_type i = 0; i < k; ++i) {
	  scalar_type a = i * h, b = (i+1) * h;
	  for (short_type l = 0; l+i < k; ++l) {
	    scalar_type c = l * h, d = (l+1) * h;
	    bgeot::sc(A) = a,c;
	    bgeot::sc(B) = b,c;
	    bgeot::sc(C) = a,d;
	    bgeot::sc(D) = b,d;
	    if (opt_gt) { 
	      A = opt_gt->transform(A, *opt_gt_pts); 
	      B = opt_gt->transform(B, *opt_gt_pts); 
	      C = opt_gt->transform(C, *opt_gt_pts); 
	      D = opt_gt->transform(D, *opt_gt_pts); 
	    }
	    pm->add_triangle_by_points(A,B,C);
	    if (l+i+1 < k)
	      pm->add_triangle_by_points(C,B,D);
	  }
	}
      }
      break;
    case 3 : 
      {
	/* based on decompositions of small cubes 
	   the number of tetrahedrons is k^3
	*/
	base_node A,B,C,D,E,F,G,H;
	for (short_type ci = 0; ci < k; ci++) {
	  scalar_type x = ci*h;
	  for (short_type cj = 0; cj < k-ci; cj++) {
	    scalar_type y=cj*h;
	    for (short_type ck = 0; ck < k-ci-cj; ck++) {
	      scalar_type z=ck*h;
	      
	      bgeot::sc(A) = x,y,z;
	      bgeot::sc(B) = x+h,y,z;
	      bgeot::sc(C) = x,y+h,z;
	      bgeot::sc(D) = x+h,y+h,z;
	      bgeot::sc(E) = x,y,z+h;
	      bgeot::sc(F) = x+h,y,z+h;
	      bgeot::sc(G) = x,y+h,z+h;
	      bgeot::sc(H) = x+h,y+h,z+h;
	      if (opt_gt) { 
		A = opt_gt->transform(A, *opt_gt_pts); 
		B = opt_gt->transform(B, *opt_gt_pts); 
		C = opt_gt->transform(C, *opt_gt_pts); 
		D = opt_gt->transform(D, *opt_gt_pts); 
		E = opt_gt->transform(E, *opt_gt_pts); 
		F = opt_gt->transform(F, *opt_gt_pts); 
		G = opt_gt->transform(G, *opt_gt_pts); 
		H = opt_gt->transform(H, *opt_gt_pts); 
	      }
	      size_type t[8];
	      t[0] = pm->add_point(A);
	      t[1] = pm->add_point(B);
	      t[2] = pm->add_point(C);
	      t[4] = pm->add_point(E);
	      if (k > 1 && ci+cj+ck < k-1) {
		t[3] = pm->add_point(D);
		t[5] = pm->add_point(F);
		t[6] = pm->add_point(G);
	      }
	      if (k > 2 && ci+cj+ck < k-2) {
		t[7] = pm->add_point(H);
	      }
	      /**
		 Note that the orientation of each tetrahedron is the same
	      */
	      pm->add_tetrahedron(t[0], t[1], t[2], t[4]);
	      if (k > 1 && ci+cj+ck < k-1) {
		pm->add_tetrahedron(t[1], t[2], t[4], t[5]);
		pm->add_tetrahedron(t[6], t[4], t[2], t[5]);
		pm->add_tetrahedron(t[2], t[3], t[5], t[1]);
		pm->add_tetrahedron(t[2], t[5], t[3], t[6]);
	      }
	      if (k > 2 && ci+cj+ck < k-2) {
		pm->add_tetrahedron(t[3], t[5], t[7], t[6]);
	      }
	    }
	  }
	}	
      }
      break;
    default : 
      //delete pm; pm = NULL;
      DAL_THROW(to_be_done_error, "Sorry, not implemented for simplices of dimension " << int(cvs->dim()));
    }
  }

  static void _structured_mesh_for_parallelepiped(bgeot::pconvex_structure cvs, 
						  bgeot::pgeometric_trans opt_gt, const std::vector<base_node> *opt_gt_pts,
						  short_type k, pgetfem_mesh &pm) {
    scalar_type h = 1./k;
    size_type n = cvs->dim();
    size_type pow2n = (size_type(1) << n);
    std::vector<size_type> strides(n);
    size_type nbpts = 1; for (size_type i=0; i < n; ++i) { strides[i] = nbpts; nbpts *= (k+1); }
    std::vector<short_type>  kcnt(n,0);
    std::vector<size_type> pids; pids.reserve(nbpts);
    base_node pt(n);

    /* insert nodes and fill pids with their numbers */
    while (true) {
      short_type kk=0;
      while (kk < n) { kcnt[kk]++; if (kcnt[kk] == k+1) { kcnt[kk] = 0; kk++; } else break; }
      for (size_type z = 0; z < n; ++z) {
	pt[z] = h*kcnt[kk];
	if ((k & (size_type(1) << z))) pt[z] += h;
      }
      if (opt_gt) pt = opt_gt->transform(pt, *opt_gt_pts);	  
      pids.push_back(pm->add_point(pt));
    }

    /* 
       insert convexes using node ids stored in 'pids'
       kcnt is again filled with 0,no need to do it 
    */
    std::vector<size_type> ppts(pow2n);
    while (true) {
      short_type kk=0;
      while (kk < n) { kcnt[kk]++; if (kcnt[kk] == k+1) { kcnt[kk] = 0; kk++; } else break; }
      
      for (k = 0; k < pow2n; ++k) {
	ppts[k] = 0;
	for (size_type z = 0; z < n; ++z) if ((k & (size_type(1) << z))) ppts[k] += strides[z];
      }
      pm->add_parallelepiped(n,ppts.begin());
    }
  }

  static void _structured_mesh_for_convex(bgeot::pconvex_structure cvs, 
					  bgeot::pgeometric_trans opt_gt, const std::vector<base_node> *opt_gt_pts,
					  short_type k, pgetfem_mesh &pm) {
    size_type nbp = cvs->basic_structure()->nb_points();
    size_type n = cvs->dim();
    /* Identifying simplexes.                                           */    
    if (nbp == n+1 && 
	cvs->basic_structure()==bgeot::simplex_structure(n)) {
      // smc.pm->write_to_file(cout);
      _structured_mesh_for_simplex(cvs,opt_gt,opt_gt_pts,k,pm);
    /* Identifying parallelepipeds.                                     */
    } else if (nbp == (size_type(1) << n) && 
	       cvs->basic_structure() == bgeot::parallelepiped_structure(n)) {
      _structured_mesh_for_parallelepiped(cvs,opt_gt,opt_gt_pts,k,pm);
    } else if (nbp == 2 * n && 
	       cvs->basic_structure() == bgeot::prism_structure(n)) {
      DAL_THROW(to_be_done_error, "Sorry, structured_mesh not implemented for prisms.");
    } else {
      DAL_THROW(to_be_done_error,
		"This element is not taken into account. Contact us");
    }
  }

  /* extract the mesh_structure on faces */
  static void _structured_mesh_of_faces(bgeot::pconvex_ref cvr, dim_type f, const getfem_mesh &m, bgeot::mesh_structure &facem)
  {
    //cerr << "structured_mesh_of_faces: face " << int(f) << ", le maillage init a " << m.nb_points() << " pts, " << m.nb_convex() << " convexes" << endl;
    facem.clear();
    dal::bit_vector on_face;
    bgeot::mesh_point_st_ct::const_iterator b = m.point_structures().begin(), e = m.point_structures().end();
    for (size_type i = 0; b != e; ++b, ++i) {
      if ((*b).is_valid()) {
        if (cvr->is_in_face(f, m.points()[i]) < 1e-12)
          on_face.add(i);
      }
    }
    //cerr << "on_face=" << on_face << endl;
    dal::bit_vector bv = m.convex_index();
    for (size_type cv = bv.take_first(); cv != size_type(-1); cv << bv) {
      for (size_type ff = 0; ff < m.structure_of_convex(cv)->nb_faces(); ++ff) {
        bgeot::ind_ref_mesh_point_ind_ct ipts = m.ind_points_of_face_of_convex(cv,ff);
        bool allin = true;
        for (size_type i=0; i < ipts.size(); ++i) if (!on_face[ipts[i]]) { allin = false; break; }
        if (allin) {
          /*cerr << "ajout de la face " << ff << " du convexe " << cv << ":";
          for (size_type i=0; i < ipts.size(); ++i) cerr << on_face[ipts[i]] << "/" << ipts[i] << " ";
          cerr << endl;*/
          facem.add_convex(m.structure_of_convex(cv)->faces_structure()[ff], ipts.begin());
        }
      }
    }
  }


  struct __str_mesh_cv {
    bgeot::pconvex_structure cvs;
    short_type n;
    bool simplex_mesh; /* true if the convex has been splited into simplexes, which were refined */
    getfem_mesh *pm;
    std::vector<bgeot::mesh_structure *> pfacem; /* array of mesh_structures for faces */
    dal::bit_vector nodes_on_edges;
    mesh_precomposite *pmp;
    bool operator < (const __str_mesh_cv &ls) const {
      if (cvs < ls.cvs) return true; if (cvs > ls.cvs) return false; 
      if (n < ls.n) return true; return false;
    }
    __str_mesh_cv(void) {}
    __str_mesh_cv(bgeot::pconvex_structure c, short_type k, bool _smesh) : 
      cvs(c), n(k), simplex_mesh(_smesh) {}
  };

  static dal::dynamic_tree_sorted<__str_mesh_cv> *__str_mesh_cv_tab = 0;

  /**
   * This function returns a mesh in pm which contains a refinement of the convex cvr
   * if force_simplexification is false, refined convexes have the same basic_structure than cvr,
   * if it is set to true, the cvr is decomposed into simplexes which are then refined.
   * TODO: move it into another file and separate the pmesh_precomposite part ?
   **/
  void structured_mesh_for_convex(bgeot::pconvex_ref cvr, short_type k,
				  pgetfem_mesh &pm, pmesh_precomposite &pmp, 
                                  bool force_simplexification) {
    size_type n = cvr->structure()->dim();
    size_type nbp = cvr->structure()->basic_structure()->nb_points();
    
    if (__str_mesh_cv_tab == 0)
      __str_mesh_cv_tab = new dal::dynamic_tree_sorted<__str_mesh_cv>();
    
    __str_mesh_cv smc(cvr->structure()->basic_structure(), k, (force_simplexification || nbp == n+1));
    
    size_type iss = __str_mesh_cv_tab->search(smc);
    if (iss == size_type(-1)) {
      smc.pm = new getfem_mesh();
      
      if (force_simplexification) {
	const bgeot::mesh_structure* splx_mesh = cvr->basic_convex_ref()->simplexified_convex();
	for (size_type ic=0; ic < splx_mesh->nb_convex(); ++ic) {
	  std::vector<base_node> cvpts(splx_mesh->nb_points_of_convex(ic));
	  bgeot::pgeometric_trans sgt = bgeot::simplex_geotrans(cvr->structure()->dim(), 1);
	  for (size_type j=0; j < cvpts.size(); ++j) {
	    cvpts[j] = cvr->basic_convex_ref()->points()[splx_mesh->ind_points_of_convex(ic)[j]];
	    //cerr << "cvpts[" << j << "]=" << cvpts[j] << endl;
	  }
	  _structured_mesh_for_convex(splx_mesh->structure_of_convex(ic), sgt, &cvpts, k, 
				      smc.pm);
	}
      } else {
	_structured_mesh_for_convex(cvr->structure(), 0, 0, k, smc.pm);
      }
      smc.pfacem.resize(cvr->structure()->nb_faces());
      for (dim_type f=0; f < cvr->structure()->nb_faces(); ++f) {
        smc.pfacem[f] = new bgeot::mesh_structure();
        _structured_mesh_of_faces(cvr, f, *smc.pm, *smc.pfacem[f]);
      }

      smc.pmp = new mesh_precomposite(*(smc.pm));
      iss = __str_mesh_cv_tab->add(smc);
    }
    pm  = (*__str_mesh_cv_tab)[iss].pm;
    pmp = (*__str_mesh_cv_tab)[iss].pmp;
  }

  const getfem_mesh *
  refined_simplex_mesh_for_convex(bgeot::pconvex_ref cvr, 
                                  short_type k) {
    pgetfem_mesh pm; pmesh_precomposite pmp; 
    structured_mesh_for_convex(cvr,k,pm,pmp,true);
    return pm;
  }

  const std::vector<bgeot::mesh_structure*>& 
  refined_simplex_mesh_for_convex_faces(bgeot::pconvex_ref cvr, 
                                  short_type k) {
    __str_mesh_cv smc(cvr->structure()->basic_structure(), k, true);    
    size_type iss = __str_mesh_cv_tab->search(smc);
    if (iss == size_type(-1)) DAL_THROW(dal::internal_error, "call refined_simplex_mesh_for_convex first (or fix me)");
    return (*__str_mesh_cv_tab)[iss].pfacem;
  }

}  /* end of namespace getfem.                                            */
