// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2007-2007 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

#include <getfem_trace_mesh_fem_level_set.h>

namespace getfem {

  //======================================================================
  // Fem part (sub space fem)
  //======================================================================

  void sub_space_fem::init() {
    cvr = org_fem->ref_convex(cv);
    dim_ = cvr->structure()->dim();
    is_equiv = real_element_defined = true;
    is_polycomp = is_pol = is_lag = false;
    es_degree = 5; /* humm ... */
    ntarget_dim = 1;

    std::stringstream nm;
    nm << "FEM_SUB_SPACE(" << org_fem->debug_name() << ", " << B << ")";
    debug_name_ = nm.str();
    
    init_cvs_node();
    if (org_fem->target_dim() != 1)
      DAL_THROW(to_be_done_error, "Vectorial fems not (yet) supported");

    base_node P(dim()); gmm::fill(P, 1./4.);
    for (size_type k = 0; k < ind.size(); ++k)
	add_node(global_dof(dim()), P);
  }

  size_type sub_space_fem::index_of_global_dof(size_type, size_type j) const
  { return ind[j]; }

  void sub_space_fem::base_value(const base_node &, 
				 base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void sub_space_fem::grad_base_value(const base_node &, 
				      base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void sub_space_fem::hess_base_value(const base_node &, 
			     base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }

  void sub_space_fem::real_base_value(const fem_interpolation_context &c,
				    base_tensor &t, bool) const {
    fem_interpolation_context c0 = c;
    base_tensor val_e;
    
    if (c0.have_pfp()) {
      c0.set_pfp(fem_precomp(org_fem, &c0.pfp()->get_point_tab()));
    } else { c0.set_pf(org_fem); }
    c0.base_value(val_e);
    
    t.mat_transp_reduction(val_e, B, 0);
  }

  void sub_space_fem::real_grad_base_value(const fem_interpolation_context &c,
					   base_tensor &t, bool) const {
    fem_interpolation_context c0 = c;
    base_tensor grad_e;
    
    if (c0.have_pfp()) {
      c0.set_pfp(fem_precomp(org_fem, &c0.pfp()->get_point_tab()));
    } else { c0.set_pf(org_fem); }
    c0.grad_base_value(grad_e);

    t.mat_transp_reduction(grad_e, B, 0);
  }
  
  void sub_space_fem::real_hess_base_value(const fem_interpolation_context &c,
				  base_tensor &t, bool) const {
    fem_interpolation_context c0 = c;
    base_tensor hess_e;
   
    if (c0.have_pfp()) {
      c0.set_pfp(fem_precomp(org_fem, &c0.pfp()->get_point_tab()));
    } else { c0.set_pf(org_fem); }
    c0.hess_base_value(hess_e);
    
    t.mat_transp_reduction(hess_e, B, 0);
  }

  //======================================================================
  // Mesh fem part
  //======================================================================

  
  void trace_mesh_fem_level_set::receipt(const MESH_CLEAR &)
  { clear(); is_adapted = false; }
  void trace_mesh_fem_level_set::receipt(const MESH_DELETE &)
  { clear(); is_adapted = false; }
  void trace_mesh_fem_level_set::clear_build_methods() {
    for (size_type i = 0; i < build_methods.size(); ++i)
      del_stored_object(build_methods[i]);
    build_methods.clear();
  }
  void trace_mesh_fem_level_set::clear(void) {
    mesh_fem::clear();
    clear_build_methods();
    is_adapted = false;
  }

  static void clear_pairs(gmm::row_matrix<gmm::rsvector<bool> > &pairs,
			  size_type i) {
    gmm::linalg_traits<gmm::rsvector<bool> >::const_iterator
      it = gmm::vect_const_begin(pairs[i]),
      ite = gmm::vect_const_end(pairs[i]);
    for (; it != ite; ++it) if (it.index() != i) pairs[i].sup(it.index());
    gmm::clear(pairs[i]);
  }
  
  trace_mesh_fem_level_set::trace_mesh_fem_level_set(const mesh_level_set &me,
						     const mesh_fem &mef,
						     unsigned degree_)
    : mesh_fem(mef.linked_mesh()), mls(me), mf(mef), degree(degree_) {
    if (mf.get_qdim() != 1)
      DAL_THROW(to_be_done_error, "base mesh_fem for mesh_fem_level_set has "
		"to be of qdim one for the moment ...");
    this->add_dependency(mls);
    is_adapted = false;
  }

  DAL_SIMPLE_KEY(special_tracemf_key, pfem);

  void trace_mesh_fem_level_set::adapt(void) {
    context_check();
    clear();
    
    pfem pf
      = classical_fem(bgeot::simplex_geotrans(linked_mesh().dim(), 1), degree);
    base_matrix G;
    std::vector<base_node> pts;
    gmm::row_matrix<gmm::rsvector<bool> > pairs(mf.nb_dof(), mf.nb_dof());
    for (dal::bv_visitor cv(linked_mesh().convex_index());
	 !cv.finished(); ++cv) {
      
      if (mls.is_convex_cut(cv)) {
      pts.resize(0);
	
	// Building a set of points of the intersection of the element with
	// the level-set.
	const mesh &msh(mls.mesh_of_convex(cv));
	bgeot::pgeometric_trans pgt = linked_mesh().trans_of_convex(cv);
	bgeot::pgeometric_trans pgt2
	  = msh.trans_of_convex(msh.convex_index().first_true());
	std::vector<mesher_level_set> mesherls0(mls.nb_level_sets());
	std::vector<mesher_level_set> mesherls1(mls.nb_level_sets());
	
	for (unsigned i = 0; i < mls.nb_level_sets(); ++i) {
	  mesherls0[i] =  mls.get_level_set(i)->mls_of_convex(cv, 0, false);
	  if (mls.get_level_set(i)->has_secondary())
	    mesherls1[i] =  mls.get_level_set(i)->mls_of_convex(cv, 1, false);
	}

	for (dal::bv_visitor i(msh.convex_index()); !i.finished(); ++i) {
	  unsigned f; bool lisin = false;
	  for (f = 0; f < pgt2->structure()->nb_faces(); ++f) {
	    
	    // for each face of the sub-elements, testing if this face is on
	    // the level-set.
	    for (unsigned ils = 0; ils < mls.nb_level_sets(); ++ils) {
	      lisin = true;
	      for (unsigned ipt = 0;
		   ipt < pgt2->structure()->nb_points_of_face(f); ++ipt) {
		lisin = lisin && (gmm::abs((mesherls0[ils])
			(msh.points_of_face_of_convex(i, f)[ipt])) < 1E-7);
		if (mls.get_level_set(ils)->has_secondary())
		  lisin = lisin && ((mesherls1[i])
			(msh.points_of_face_of_convex(i, f)[ipt]) < -1E-7);
	      }
	      if (lisin) break;
	    }
	    if (!lisin) continue;
	    // A face is detected to be on the level-set,
	    // building a set of points on this face.
	    
	    vectors_to_base_matrix(G, msh.points_of_convex(i));
	    bgeot::geotrans_interpolation_context c(msh.trans_of_convex(i),
						    pf->node_of_dof(0,0), G);
	    
	    for (size_type j=0; j<pf->structure(0)->nb_points_of_face(f);++j) {
	      c.set_xref(pf->node_convex(0).points_of_face(f)[j]);
	      cout << "Adding to point stock " << c.xreal() << endl;
	      pts.push_back(c.xreal());
	    }
	  }
	}
	// Selecting the right number of points
	std::vector<size_type> selection;
	cout << "begining selection\n";
	for (size_type k = 0; k < pf->structure(0)->nb_points_of_face(0);++k) {
	  size_type ipt = size_type(-1);
	  scalar_type maxdmin = scalar_type(0);
	  
	  for (size_type i = 0; i < pts.size(); ++i) {
	    scalar_type dmin = (selection.size() == 0)
	      ? gmm::vect_norm2(pts[i])
	      : gmm::vect_dist2(pts[i], pts[selection[0]]);
	    for (size_type j = 1; j < selection.size(); ++j) {
	      dmin = std::min(dmin,gmm::vect_dist2(pts[i], pts[selection[j]]));
	      if (dmin <= maxdmin) break;
	    }
	    if (dmin > maxdmin) { maxdmin = dmin; ipt = i; }
	  }
	  
	  if (ipt ==  size_type(-1))
	    DAL_THROW(failure_error, "Trace fem: not enought points");
	  cout << "adding " << pts[ipt] << endl;
	  selection.push_back(ipt);
	}
	

	// Now that a set of points is selected, evaluating the base function
	// of the fem on it and retaining the two base functions having the
	// largest value on each point. Adding them to the list of candidate
	// pairs of dof.

	pfem pfcv = mf.fem_of_element(cv);
	vectors_to_base_matrix(G, linked_mesh().points_of_convex(cv));

	fem_interpolation_context c(pgt, pfcv, pts[selection[0]], G, cv);
	base_tensor t;

	for (size_type i = 0; i < selection.size(); ++i) { 
	  c.set_xref(pts[selection[i]]);
	  c.base_value(t);
	  
	  size_type n1 = size_type(-1), n2(0);
	  scalar_type y1(0), y2(0);
	  
	  for (size_type j = 0; j < pfcv->nb_dof(cv); ++j) {
	    if (t[j] > y1) { n1 = j; y1 = t[j]; }
	    else if (t[j] > y2) { n2 = j; y2 = t[j]; }
	  }
	  if (n1 == size_type(-1)) DAL_INTERNAL_ERROR("");
	  size_type nd1 = mf.ind_dof_of_element(cv)[n1];
	  size_type nd2 = mf.ind_dof_of_element(cv)[n2];
	  
	  if (n2 < scalar_type(1e-4)) pairs(nd1, nd1) = true;
	  else pairs(nd1, nd2) = pairs(nd2, nd1) = true;
	}
	

      }
      
    }

    // At this stage, a certain number of pairs of dof has been produced.
    // Now, a selection is made based on the following criteria:
    //   - Singletons are selected, other pairs containing
    //     this dof are eliminated.
    //   - Between pairs having a dof exclusively for their own and attached
    //     to the same other dof, one pair is arbitrary selected and the others
    //     are eliminated.
    //   - If all the pairs are linked with both the two dofs, one pair is
    //     arbitrary selected, the pairs having a common dof with this pair are
    //     eliminated.

    std::vector<size_type> indpairing(mf.nb_dof(), size_type(-1));

    bool ttouched;
    size_type nbdof(0);
    do {
      int nb_nnz(0), lasti(0), lastj(0);
      ttouched = false;
      for (size_type i = 0; i < mf.nb_dof(); ++i) { // to be optimized ... 
	if (gmm::nnz(pairs[i]) > 0) {
	  lasti = i;
	  lastj = gmm::vect_const_begin(pairs[i]).index();
	  ++nb_nnz;
	  if (pairs(i,i) == true) {
	    indpairing[i] = nbdof++;
	    ttouched = true;
	    clear_pairs(pairs, i);
	  }
	  else {
	    if (gmm::nnz(pairs[i]) == 1) {
	      indpairing[lasti] = indpairing[lastj] = nbdof++;
	      clear_pairs(pairs, lasti); clear_pairs(pairs, lastj);
	      ttouched = true;
	    }
	  }
	  
	}
      }
      if (nb_nnz > 0 && !ttouched) {
	indpairing[lasti] = indpairing[lastj] = nbdof++;
	clear_pairs(pairs, lasti); clear_pairs(pairs, lastj);
	ttouched = true;
      }
    } while (ttouched);

    // Now that the convenient pairs of dofs are selected, the special
    // fem are built.

    std::vector<size_type> glob_dof, indlist1, indlist2;
    for (dal::bv_visitor cv(linked_mesh().convex_index());
	 !cv.finished(); ++cv) {
      if (mls.is_convex_cut(cv)) {
	glob_dof.resize(0); indlist1.resize(0); indlist2.resize(0); 

	for (size_type i = 0; i < mf.nb_dof_of_element(cv); ++i) {
	  size_type id = indpairing[mf.ind_dof_of_element(cv)[i]];
	  if (id != size_type(-1)) {
	    std::vector<size_type>::iterator it = 
	      std::find(glob_dof.begin(), glob_dof.end(), id);
	    if (it  == glob_dof.end()) {
	      glob_dof.push_back(id);
	      indlist1.push_back(i);
	      indlist2.push_back(i);
	    }
	    else
	      indlist2[it - glob_dof.begin()] = i;
	  }
	}

	base_matrix B(glob_dof.size(), mf.nb_dof_of_element(cv));
	for (size_type i = 0; i < glob_dof.size(); ++i) {
	  B(i, indlist1[i]) = scalar_type(1);
	  B(i, indlist2[i]) = scalar_type(1);
	}
	
	pfem pfnew = new sub_space_fem(mf.fem_of_element(cv), glob_dof, B, cv);
	build_methods.push_back(pfnew);
	set_finite_element(cv, pfnew);
      }
    }
    is_adapted = true; touch();
  }


}  /* end of namespace getfem.                                            */

