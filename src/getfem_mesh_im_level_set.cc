// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2005-2008 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesher.h"
#include "getfem/bgeot_kdtree.h"

namespace getfem {

  void mesh_im_level_set::update_from_context(void) const
  { is_adapted = false; }

  void mesh_im_level_set::receipt(const MESH_CLEAR &)
  { clear(); is_adapted = false; }

  void mesh_im_level_set::receipt(const MESH_DELETE &)
  { clear(); is_adapted = false; }

  void mesh_im_level_set::clear_build_methods() {
    for (size_type i = 0; i < build_methods.size(); ++i)
      del_stored_object(build_methods[i]);
    build_methods.clear();
    cut_im.clear();
  }

  void mesh_im_level_set::clear(void) {
    mesh_im::clear();
    clear_build_methods();
    is_adapted = false;
  }  

  mesh_im_level_set::mesh_im_level_set(mesh_level_set &me,
				       int integrate_where_, 
				       pintegration_method reg,
				       pintegration_method sing) 
    : mesh_im(me.linked_mesh()), mls(me), cut_im(me.linked_mesh()),
      integrate_where(integrate_where_) {
    set_simplex_im(reg, sing);
    this->add_dependency(mls);
    is_adapted = false;
  }

  pintegration_method 
  mesh_im_level_set::int_method_of_element(size_type cv) const {
    if (!is_adapted) const_cast<mesh_im_level_set *>(this)->adapt();
    if (cut_im.convex_index().is_in(cv)) 
      return cut_im.int_method_of_element(cv); 
    else {
      if (ignored_im.is_in(cv)) //integrate_where == INTEGRATE_BOUNDARY)
	return getfem::im_none();

      return mesh_im::int_method_of_element(cv);
    }
  }

  DAL_SIMPLE_KEY(special_imls_key, papprox_integration);

    /* only for INTEGRATE_INSIDE or INTEGRATE_OUTSIDE */
  mesh_im_level_set::bool2 mesh_im_level_set::is_point_in_selected_area2
  (const std::vector<mesher_level_set> &mesherls0,
   const std::vector<mesher_level_set> &mesherls1,
				  const base_node& P) {
    bool isin = true;
    int isbin = 0;
    for (unsigned i = 0; i < mls.nb_level_sets(); ++i) {
      isin = isin && ((mesherls0[i])(P) < 0);
      if (gmm::abs((mesherls0[i])(P)) < 1e-7)
	isbin = i+1;
      if (mls.get_level_set(i)->has_secondary())
	isin = isin && ((mesherls1[i])(P) < 0);
    }
    bool2 b; 
    b.in = ((integrate_where & INTEGRATE_OUTSIDE)) ? !isin : isin;
    b.bin = isbin;
    return b;
  }
  

  /* very rustic set operations evaluator */
  struct is_in_eval {
    dal::bit_vector in;  // levelsets for which the point is "inside"
    dal::bit_vector bin; // levelsets for which the point is on the boundary
    typedef mesh_im_level_set::bool2 bool2;
    bool2 do_expr(const char *&s) {
      bool2 r;
      if (*s == '(') {
	r = do_expr(++s);
	GMM_ASSERT1(*s++ == ')', 
		    "expecting ')' in csg expression at '" << s-1 << "'");
      } else if (*s == '!') { // complementary
	r = do_expr(++s); r.in = !r.in;
      } else if (*s >= 'a' && *s <= 'z') {
	unsigned idx = (*s) - 'a';
	r.in  = in.is_in(idx);
	r.bin = bin.is_in(idx) ? idx+1 : 0; 
	++s;
      } else 
	GMM_ASSERT1(false, "parse error in csg expression at '" << s << "'");
      if (*s == '+') { // Union
	//cerr << "s = " << s << ", r = " << r << "\n";
	bool2 a = r, b = do_expr(++s);
	//cerr << "->b = " << b << "\n";
	r.in = b.in || a.in;
	if      (b.bin && !a.in) r.bin = b.bin;
	else if (a.bin && !b.in) r.bin = a.bin;
	else r.bin = 0;
      } else if (*s == '-') { // Set difference
	bool2 a = r, b = do_expr(++s);
	r.in  = a.in && !b.in;
	if      (a.bin && !b.in) r.bin = a.bin;
	else if (a.in  && b.bin) r.bin = b.bin;
	else r.bin = 0;
      } else if (*s == '*') { // Intersection
	bool2 a = r, b = do_expr(++s);
	r.in  = a.in && b.in;
	if      (a.bin && b.in)  r.bin = a.bin;
	else if (a.in  && b.bin) r.bin = b.bin;
	else r.bin = 0;
      }
      return r;
    }
    bool2 is_in(const char*s) { 
      bool2 b = do_expr(s); 
      GMM_ASSERT1(!(*s), "parse error in CSG expression at " << s);
      return b;
    }
    void check() {
      const char *s[] = { "a*b*c*d",
			  "a+b+c+d",
			  "(a+b+c+d)",
			  "d*(a+b+c)",
			  "(a+b)-(c+d)",
			  "((a+b)-(c+d))",
			  "!a",
			  0 };
      for (const char **p = s; *p; ++p) 
	cerr << *p << "\n";
      for (unsigned c=0; c < 16; ++c) {
	in[0] = (c&1); bin[0] = 1;
	in[1] = (c&2); bin[1] = 1;
	in[2] = (c&4); bin[2] = 1;
	in[3] = (c&8); bin[3] = 1;
	cerr << in[0] << in[1] << in[2] << in[3] << ": ";
	for (const char **p = s; *p; ++p) {
	  bool2 b = is_in(*p);
	  cerr << b.in << "/" << b.bin << " ";
	}
	cerr << "\n";
      }
    }
  };
  
  mesh_im_level_set::bool2 
  mesh_im_level_set::is_point_in_selected_area
           (const std::vector<mesher_level_set> &mesherls0,
	    const std::vector<mesher_level_set> &mesherls1,
	    const base_node& P) {
    is_in_eval ev;
    for (unsigned i = 0; i < mls.nb_level_sets(); ++i) {
      bool sec = mls.get_level_set(i)->has_secondary();
      scalar_type d1 = (mesherls0[i])(P);
      scalar_type d2 = (sec ? (mesherls1[i])(P) : -1);
      if (d1 < 0 && d2 < 0) ev.in.add(i);
      if ((integrate_where & INTEGRATE_OUTSIDE) /*&& !sec*/)
	ev.in[i].flip();

      if (gmm::abs(d1) < 1e-7 && d2 < 1e-7) 
	ev.bin.add(i);
      /*if (integrate_where == INTEGRATE_BOUNDARY)
	cerr << "is_point_in_selected_area: i=" << i << ", in = " << ev.in.is_in(i) << ", bin=" << ev.bin.is_in(i) << " d1 = " << d1 << ", d2=" << d2 << " where=" << integrate_where << " csg=" << ls_csg_description << " ls=" << (void*)mls.get_level_set(i) << ", sec=" << sec << "\n";*/
    }
    

    bool2 r;
    if (ls_csg_description.size()) 
      r = ev.is_in(ls_csg_description.c_str());
    else {
      r.in  = (ev.in.card() == mls.nb_level_sets());
      r.bin = (ev.bin.card() >= 1 && ev.in.card() >= mls.nb_level_sets()-1);
    }
    /*bool2 r2 = is_point_in_selected_area2(mesherls0,mesherls1,P);
    if (r2.in != r.in || r2.bin != r.bin) {
      cerr << "ev.in = " << ev.in << ", bin=" << ev.bin<<"\n";
      cerr << "is_point_in_selected_area2("<<P <<"): r="<<r.in<<"/"<<r.bin
	   << ", r2=" << r2.in<<"/"<<r2.bin <<"\n";
      assert(0);
      }*/

    return r;
  }
  
  void mesh_im_level_set::build_method_of_convex(size_type cv) {
    const mesh &msh(mls.mesh_of_convex(cv));
    GMM_ASSERT3(msh.convex_index().card() != 0, "");
    base_matrix G;
    base_node B;

    std::vector<mesher_level_set> mesherls0(mls.nb_level_sets());
    std::vector<mesher_level_set> mesherls1(mls.nb_level_sets());
    dal::bit_vector convexes_arein;

    //std::fstream totof("totof", std::ios::out | std::ios::app);
    for (unsigned i = 0; i < mls.nb_level_sets(); ++i) {
      mesherls0[i] =  mls.get_level_set(i)->mls_of_convex(cv, 0, false);
      if (mls.get_level_set(i)->has_secondary())
	mesherls1[i] =  mls.get_level_set(i)->mls_of_convex(cv, 1, false);
    }

    if (integrate_where != (INTEGRATE_ALL)) {
      for (dal::bv_visitor scv(msh.convex_index()); !scv.finished(); ++scv) {
	B = gmm::mean_value(msh.points_of_convex(scv));
	convexes_arein[scv] = 
	  is_point_in_selected_area(mesherls0, mesherls1, B).in;
      }
    }
    
    bgeot::pgeometric_trans pgt = linked_mesh().trans_of_convex(cv);
    bgeot::pgeometric_trans pgt2
      = msh.trans_of_convex(msh.convex_index().first_true());
    dim_type n = pgt->dim();

    if (base_singular_pim) GMM_ASSERT1
      ((n != 2 ||
	base_singular_pim->structure()== bgeot::parallelepiped_structure(2))
       && (n != 3
	   || base_singular_pim->structure() == bgeot::prism_structure(3))
       && (n >= 2) && (n <= 3),
       "Base integration method for quasi polar integration not convenient");

    approx_integration *new_approx = new approx_integration(pgt->convex_ref());
    base_matrix KK(n,n), CS(n,n);
    base_matrix pc(pgt2->nb_points(), n);
    std::vector<size_type> ptsing;

    //cerr << "testing convex " << cv << ", " << msh.convex_index().card() << " subconvexes\n";

    for (dal::bv_visitor i(msh.convex_index()); !i.finished(); ++i) {
      papprox_integration pai = regular_simplex_pim->approx_method();
      
      if ((integrate_where != INTEGRATE_ALL) &&
	  !convexes_arein[i]) continue;
      
      if (base_singular_pim && mls.crack_tip_convexes().is_in(cv)) {
	ptsing.resize(0);
	unsigned sing_ls = unsigned(-1);

	for (unsigned ils = 0; ils < mls.nb_level_sets(); ++ils)
	  if (mls.get_level_set(ils)->has_secondary()) {
	    for (unsigned ipt = 0; ipt <= n; ++ipt) {
	      if (gmm::abs((mesherls0[ils])(msh.points_of_convex(i)[ipt]))
		  < 1E-10
		  && gmm::abs((mesherls1[ils])(msh.points_of_convex(i)[ipt]))
		  < 1E-10) {
		if (sing_ls == unsigned(-1)) sing_ls = ils;
		GMM_ASSERT1(sing_ls == ils, "Two singular point in one "
			    "sub element. To be done.");
		ptsing.push_back(ipt);
	      }
	    }
	  }
	assert(ptsing.size() < n);

	if (ptsing.size() > 0) {
	  std::stringstream sts;
	  sts << "IM_QUASI_POLAR(" << name_of_int_method(base_singular_pim)
	      << ", " << ptsing[0];
	  if (ptsing.size() > 1) sts << ", " <<  ptsing[1];
	  sts << ")";
	  //cout << "Singular int method : " << sts.str() << endl;
	  pai = int_method_descriptor(sts.str())->approx_method();
	}
      }

      base_matrix G2;
      vectors_to_base_matrix(G2, linked_mesh().points_of_convex(cv));
      bgeot::geotrans_interpolation_context
	cc(linked_mesh().trans_of_convex(cv), pai->point(0), G2);

      if (integrate_where & (INTEGRATE_INSIDE | INTEGRATE_OUTSIDE)) {
		
	vectors_to_base_matrix(G, msh.points_of_convex(i));
	bgeot::geotrans_interpolation_context c(msh.trans_of_convex(i),
						pai->point(0), G);

	for (size_type j = 0; j < pai->nb_points_on_convex(); ++j) {
	  c.set_xref(pai->point(j));
	  pgt2->poly_vector_grad(pai->point(j), pc);
	  gmm::mult(G,pc,KK);
	  scalar_type J = gmm::lu_det(KK);
	  new_approx->add_point(c.xreal(), pai->coeff(j) * gmm::abs(J));

	  /*if (integrate_where == INTEGRATE_INSIDE) {
	    cc.set_xref(c.xreal());
	    totof << cc.xreal()[0] << "\t" << cc.xreal()[1] << "\t"
	    << pai->coeff(j) * gmm::abs(J) << "\n";
	    }*/
	}
      }

      // pgt2 = msh.trans_of_convex(i);

      for (unsigned f = 0; f < pgt2->structure()->nb_faces(); ++f) {
	short_type ff = short_type(-1);
	unsigned isin = unsigned(-1);

	if (integrate_where == INTEGRATE_BOUNDARY) {
	  bool lisin = true;
	  for (unsigned ipt = 0; ipt < 
		 pgt2->structure()->nb_points_of_face(f); ++ipt) {
	    const base_node P = msh.points_of_face_of_convex(i, f)[ipt];
	    isin = is_point_in_selected_area(mesherls0, mesherls1, P).bin;
	    //cerr << P << ":" << isin << " ";
	    if (!isin) { lisin = false; break; }
	  }
	  /*cerr << "testing face " << f << " sub " << i << " of cv " << cv << " "
	    << " -> lisin=" << lisin << "\n";*/
	  if (!lisin) continue;
	  else isin--;
	} else {
	  B = gmm::mean_value(msh.points_of_face_of_convex(i, f));
	  if (pgt->convex_ref()->is_in(B) < -1E-7) continue;
	  for (short_type fi = 0; fi < pgt->structure()->nb_faces(); ++fi) {
	    if (gmm::abs(pgt->convex_ref()->is_in_face(fi, B)) < 1E-6) ff = fi;
	  }
	  GMM_ASSERT3(ff != short_type(-1), "");
	}
	  
	vectors_to_base_matrix(G, msh.points_of_convex(i));
	bgeot::geotrans_interpolation_context c(msh.trans_of_convex(i),
						pai->point(0), G);
	
	
	for (size_type j = 0; j < pai->nb_points_on_face(f); ++j) {
	  c.set_xref(pai->point_on_face(f, j));
	  base_small_vector un = pgt2->normals()[f], up(msh.dim());
	  gmm::mult(c.B(), un, up);
	  scalar_type nup = gmm::vect_norm2(up);

	  scalar_type nnup(1);
	  if (integrate_where == INTEGRATE_BOUNDARY) {
	    cc.set_xref(c.xreal());
	    mesherls0[isin].grad(c.xreal(), un);
	    un /= gmm::vect_norm2(un);
	    gmm::mult(cc.B(), un, up);
	    nnup = gmm::vect_norm2(up);
	  }
	  new_approx->add_point(c.xreal(), pai->coeff_on_face(f, j)
				* gmm::abs(c.J()) * nup * nnup, ff);

	  /*if (integrate_where == INTEGRATE_BOUNDARY) {
	    static double ssum = 0.0;
	    ssum += pai->coeff_on_face(f, j) * gmm::abs(c.J()) * nup * nnup;
	    cout << "add crack point " << c.xreal() << " : "
		 << pai->coeff_on_face(f, j) * gmm::abs(c.J()) * nup * nnup
		 << " sum = " << ssum << endl;
	  }*/
	  /*if (integrate_where == INTEGRATE_BOUNDARY) {
	    cc.set_xref(c.xreal());
	    totof << cc.xreal()[0] << "\t" << cc.xreal()[1] << "\n";
	  }
	  */
	} 
      }
    }

    new_approx->valid_method();

    if (new_approx->nb_points()) {
      pintegration_method pim = new integration_method(new_approx);
      dal::add_stored_object(new special_imls_key(new_approx), pim,
			     new_approx->ref_convex(),
			     &(new_approx->integration_points()));
      build_methods.push_back(pim);
      cut_im.set_integration_method(cv, pim);
    }
    else {
      delete new_approx;
    }
  }

  void mesh_im_level_set::adapt(void) {
    context_check();
    clear_build_methods();
    ignored_im.clear();
    for (dal::bv_visitor cv(linked_mesh().convex_index()); 
	 !cv.finished(); ++cv) {
      if (mls.is_convex_cut(cv)) {
	build_method_of_convex(cv);
      }

      if (!cut_im.convex_index().is_in(cv)) {
	/* not exclusive with mls.is_convex_cut ... sometimes, cut cv
	   contains no integration points.. */

	if (integrate_where == INTEGRATE_BOUNDARY) {
	  ignored_im.add(cv);
	} else if (integrate_where != (INTEGRATE_OUTSIDE|INTEGRATE_INSIDE)) {
	  /* remove convexes that are not in the integration area */
	  std::vector<mesher_level_set> mesherls0(mls.nb_level_sets());
	  std::vector<mesher_level_set> mesherls1(mls.nb_level_sets());
	  for (unsigned i = 0; i < mls.nb_level_sets(); ++i) {
	    mesherls0[i] = mls.get_level_set(i)->mls_of_convex(cv, 0, false);
	    if (mls.get_level_set(i)->has_secondary())
	      mesherls1[i] = mls.get_level_set(i)->mls_of_convex(cv, 1, false);
	  }

	  base_node B(gmm::mean_value(linked_mesh().trans_of_convex(cv)
				      ->convex_ref()->points()));
	  if (!is_point_in_selected_area(mesherls0, mesherls1, B).in)
	    ignored_im.add(cv);
	}
      }
    }
    // cerr << "mesh_im_level_set: integrate = " << integrate_where
    //      << ", ignored = " << ignored_im << "\n";
    is_adapted = true; touch();
  }

  
}  /* end of namespace getfem.                                             */



