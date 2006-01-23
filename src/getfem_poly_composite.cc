// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_poly_composite.cc : polynomials by parts
//           
// Date    : August 26, 2002.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2002-2005 Yves Renard
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

#include <dal_singleton.h>
#include <bgeot_comma_init.h>
#include <getfem_poly_composite.h>

namespace getfem
{ 
  mesh_precomposite::mesh_precomposite(const mesh &m) {
    msh = &m;
    elt.resize(m.nb_convex());
    det.resize(m.nb_convex());
    orgs.resize(m.nb_convex());
    gtrans.resize(m.nb_convex());
    for (size_type i = 0; i <= m.points().index().last_true(); ++i) {
      vertexes.add(m.points()[i]);
    }
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv) {
      
      bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
      size_type N = pgt->structure()->dim();
      size_type P = m.dim();
      if (!(pgt->is_linear()) || N != P) 
	DAL_THROW(internal_error, "Bad geometric transformations.");
    
      base_matrix G(P, pgt->nb_points());
      base_matrix pc(pgt->nb_points() , N);
      base_matrix B0(N, P);
    
      vectors_to_base_matrix(G, m.points_of_convex(cv));
      
      base_node x(N);x.fill(0.);
      pgt->gradient(x, pc);
    
      gmm::mult(gmm::transposed(pc), gmm::transposed(G), B0);
      cout << "B0 = " << B0 << endl;
      det[cv] = gmm::lu_inverse(B0);
      gtrans[cv] = B0;
      orgs[cv] = m.points_of_convex(cv)[0];
    }
  }

  scalar_type polynomial_composite::eval(const base_node &pt) const {
    base_node p0(mp->dim()), p1(mp->dim());
    std::fill(mp->elt.begin(), mp->elt.end(), true);
    bgeot::mesh_structure::ind_cv_ct::const_iterator itc, itce;
    
    mesh_precomposite::PTAB::const_sorted_iterator
      it1 = mp->vertexes.sorted_ge(pt), it2 = it1;    
    size_type i1 = it1.index(), i2;

    --it2; i2 = it2.index();
    

    while (i1 != size_type(-1) || i2 != size_type(-1)) {
      if (i1 != size_type(-1)) {
	const bgeot::mesh_structure::ind_cv_ct &tc
	  = mp->linked_mesh().convex_to_point(i1);
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
	const bgeot::mesh_structure::ind_cv_ct &tc
	  = mp->linked_mesh().convex_to_point(i2);
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


  /* build a regularly refined mesh for a simplex of dimension <= 3.
     All simplexes have the same orientation as the original simplex.
   */
  static void structured_mesh_for_simplex_(bgeot::pconvex_structure cvs, 
				      bgeot::pgeometric_trans opt_gt,
				      const std::vector<base_node> *opt_gt_pts,
				      short_type k, pmesh &pm) {
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
	base_node A(2),B(2),C(2),D(2);
	for (short_type i = 0; i < k; ++i) {
	  scalar_type a = i * h, b = (i+1) * h;
	  for (short_type l = 0; l+i < k; ++l) {
	    scalar_type c = l * h, d = (l+1) * h;
	    A[0] = a; A[1] = c;
	    B[0] = b; B[1] = c;
	    C[0] = a; C[1] = d;
	    D[0] = b; D[1] = d;
	    if (opt_gt) { 
	      A = opt_gt->transform(A, *opt_gt_pts); 
	      B = opt_gt->transform(B, *opt_gt_pts); 
	      C = opt_gt->transform(C, *opt_gt_pts); 
	      D = opt_gt->transform(D, *opt_gt_pts); 
	    }
	    /* add triangle with correct orientation (det [B-A;C-A] > 0) */
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
	      t[3] = t[5] = t[6] = t[7] = size_type(-1);
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
		 (det > 0)
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

  static void structured_mesh_for_parallelepiped_(bgeot::pconvex_structure cvs, 
						  bgeot::pgeometric_trans opt_gt, const std::vector<base_node> *opt_gt_pts,
						  short_type k, pmesh &pm) {
    scalar_type h = 1./k;
    size_type n = cvs->dim();
    size_type pow2n = (size_type(1) << n);
    size_type powkn = 1; for (size_type i=0; i < n; ++i) powkn *= k;
    std::vector<size_type> strides(n);
    size_type nbpts = 1; for (size_type i=0; i < n; ++i) { strides[i] = nbpts; nbpts *= (k+1); }
    std::vector<short_type> kcnt; kcnt.resize(n+1,0);
    std::vector<size_type> pids; pids.reserve(nbpts);
    base_node pt(n);
    size_type kk;
    /* insert nodes and fill pids with their numbers */
    while (kcnt[n] == 0) {
      for (size_type z = 0; z < n; ++z)
	pt[z] = h*kcnt[z];
      if (opt_gt) pt = opt_gt->transform(pt, *opt_gt_pts);	  
      pids.push_back(pm->add_point(pt));
      kk=0; while (kk <= n) { kcnt[kk]++; if (kcnt[kk] == k+1) { kcnt[kk] = 0; kk++; } else break; }
    }

    /* 
       insert convexes using node ids stored in 'pids'
    */
    std::vector<size_type> ppts(pow2n);
    kcnt[n] = 0;
    while (kcnt[n] == 0) {
      for (kk = 0; kk < pow2n; ++kk) {
	size_type pos = 0;
	for (size_type z = 0; z < n; ++z) {
          pos += kcnt[z]*strides[z];
          if ((kk & (size_type(1) << z))) pos += strides[z];
        }
        ppts[kk] = pids.at(pos);
      }
      pm->add_convex(bgeot::parallelepiped_linear_geotrans(n), ppts.begin());
      kcnt[(kk = 0)]++;
      while (kk < n && kcnt[kk] == k) { kcnt[kk] = 0; kcnt[++kk]++; }
    }
  }

  static void structured_mesh_for_convex_(bgeot::pconvex_structure cvs, 
					  bgeot::pgeometric_trans opt_gt, const std::vector<base_node> *opt_gt_pts,
					  short_type k, pmesh &pm) {
    size_type nbp = cvs->basic_structure()->nb_points();
    size_type n = cvs->dim();
    /* Identifying simplexes.                                           */    
    if (nbp == n+1 && 
	cvs->basic_structure()==bgeot::simplex_structure(n)) {
      // smc.pm->write_to_file(cout);
      structured_mesh_for_simplex_(cvs,opt_gt,opt_gt_pts,k,pm);
    /* Identifying parallelepipeds.                                     */
    } else if (nbp == (size_type(1) << n) && 
	       cvs->basic_structure() == bgeot::parallelepiped_structure(n)) {
      structured_mesh_for_parallelepiped_(cvs,opt_gt,opt_gt_pts,k,pm);
    } else if (nbp == 2 * n && 
	       cvs->basic_structure() == bgeot::prism_structure(n)) {
      DAL_THROW(to_be_done_error, "Sorry, structured_mesh not implemented for prisms.");
    } else {
      DAL_THROW(to_be_done_error,
		"This element is not taken into account. Contact us");
    }
  }

  /* extract the mesh_structure on faces */
  static void structured_mesh_of_faces_(bgeot::pconvex_ref cvr, dim_type f,
					const mesh &m,
					bgeot::mesh_structure &facem) {
    facem.clear();
    dal::bit_vector on_face;
    for (size_type i = 0; i < m.nb_max_points(); ++i) {
      if (m.is_point_valid(i)) {
        if (cvr->is_in_face(f, m.points()[i]) < 1e-12) on_face.add(i);
      }
    }
    //cerr << "on_face=" << on_face << endl;
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv) {
      for (size_type ff = 0; ff < m.structure_of_convex(cv)->nb_faces(); ++ff) {
        bgeot::mesh_structure::ind_pt_face_ct ipts = m.ind_points_of_face_of_convex(cv,ff);
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

  DAL_TRIPLE_KEY(str_mesh_key, bgeot::pconvex_structure, short_type, bool);

  struct str_mesh_cv__  : virtual public dal::static_stored_object {
    bgeot::pconvex_structure cvs;
    short_type n;
    bool simplex_mesh; /* true if the convex has been splited into simplexes, which were refined */
    mesh *pm;
    std::vector<bgeot::mesh_structure *> pfacem; /* array of mesh_structures for faces */
    dal::bit_vector nodes_on_edges;
    mesh_precomposite *pmp;
    str_mesh_cv__(void) : pm(0), pmp(0) {}
    str_mesh_cv__(bgeot::pconvex_structure c, short_type k, bool smesh_) : 
      cvs(c), n(k), simplex_mesh(smesh_) {}
    ~str_mesh_cv__() { 
      if (pm) delete pm; if (pmp) delete pmp; pm = 0; pmp = 0;
      for (size_type i=0; i < pfacem.size(); ++i) delete pfacem[i];
    }
  };

  typedef boost::intrusive_ptr<const str_mesh_cv__> pstr_mesh_cv__;

  /**
   * This function returns a mesh in pm which contains a refinement of the convex cvr
   * if force_simplexification is false, refined convexes have the same basic_structure than cvr,
   * if it is set to true, the cvr is decomposed into simplexes which are then refined.
   * TODO: move it into another file and separate the pmesh_precomposite part ?
   **/
  void structured_mesh_for_convex(bgeot::pconvex_ref cvr, short_type k,
				  pmesh &pm, pmesh_precomposite &pmp, 
                                  bool force_simplexification) {
    size_type n = cvr->structure()->dim();
    size_type nbp = cvr->structure()->basic_structure()->nb_points();

    force_simplexification = force_simplexification || (nbp == n+1);

    dal::pstatic_stored_object o
      = dal::search_stored_object(str_mesh_key(cvr->structure()
					       ->basic_structure(), k,
					       force_simplexification));
    pstr_mesh_cv__ psmc;
    if (o) {
      psmc = dal::stored_cast<str_mesh_cv__>(o);
    } 
    else {
      str_mesh_cv__ &smc(*(new str_mesh_cv__(cvr->structure()
					     ->basic_structure(), k,
					     force_simplexification)));
      psmc = &smc;
      
      smc.pm = new mesh();
      
      if (force_simplexification) {
	// cout << "cvr = " << int(cvr->structure()->dim()) << " : "
	//      << cvr->structure()->nb_points() << endl;
	const bgeot::mesh_structure* splx_mesh
	  = cvr->basic_convex_ref()->simplexified_convex();
	/* splx_mesh is correctly oriented */
	for (size_type ic=0; ic < splx_mesh->nb_convex(); ++ic) {
	  std::vector<base_node> cvpts(splx_mesh->nb_points_of_convex(ic));
	  bgeot::pgeometric_trans sgt
	    = bgeot::simplex_geotrans(cvr->structure()->dim(), 1);
	  for (size_type j=0; j < cvpts.size(); ++j) {
	    cvpts[j] = cvr->basic_convex_ref()->points()
	      [splx_mesh->ind_points_of_convex(ic)[j]];
	    //cerr << "cvpts[" << j << "]=" << cvpts[j] << endl;
	  }
	  structured_mesh_for_convex_(splx_mesh->structure_of_convex(ic),
				      sgt, &cvpts, k, smc.pm);
	}
      } else {
	structured_mesh_for_convex_(cvr->structure(), 0, 0, k, smc.pm);
      }
      smc.pfacem.resize(cvr->structure()->nb_faces());
      for (dim_type f=0; f < cvr->structure()->nb_faces(); ++f) {
        smc.pfacem[f] = new bgeot::mesh_structure();
        structured_mesh_of_faces_(cvr, f, *smc.pm, *smc.pfacem[f]);
      }

      smc.pmp = new mesh_precomposite(*(smc.pm));
      dal::add_stored_object(new str_mesh_key(cvr->structure()
					      ->basic_structure(), k,
					      force_simplexification), psmc,
			     cvr->structure()->basic_structure());
    }
    pm  = psmc->pm;
    pmp = psmc->pmp;
  }

  const mesh *
  refined_simplex_mesh_for_convex(bgeot::pconvex_ref cvr, short_type k) {
    pmesh pm; pmesh_precomposite pmp; 
    structured_mesh_for_convex(cvr,k,pm,pmp,true);
    return pm;
  }

  const std::vector<bgeot::mesh_structure*>& 
  refined_simplex_mesh_for_convex_faces(bgeot::pconvex_ref cvr, 
                                  short_type k) {
    dal::pstatic_stored_object o
      = dal::search_stored_object(str_mesh_key(cvr->structure()
					       ->basic_structure(), k, true));
    if (o) {
      pstr_mesh_cv__ psmc = dal::stored_cast<str_mesh_cv__>(o);
      return psmc->pfacem;
    } 
    else DAL_THROW(dal::internal_error,
		   "call refined_simplex_mesh_for_convex first (or fix me)");
  }

}  /* end of namespace getfem.                                            */
