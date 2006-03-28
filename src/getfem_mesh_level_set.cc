// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2005-2006 Julien Pommier
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

#include <getfem_mesh_level_set.h>


namespace getfem {
  

#ifdef DEBUG_LS
#include <asm/msr.h>
struct Chrono {
  unsigned long long tprev;
  unsigned long long acc;
  unsigned long long tmax;
  bool running; unsigned cnt;
  Chrono() : acc(0), tmax(0), running(false), cnt(0) {}
  void tic() { rdtscll(tprev); running = true; ++cnt; }
  double toc() {
    if (!running) return 0.; running = false;
    unsigned long long t; rdtscll(t);
    t -= tprev;
    acc += t; tmax = std::max(tmax, t);
    return double(t)/2.8e9;
  }
  double total() { return double(acc) / 2.8e9; }
  double max() { return double(tmax) / 2.8e9; }
  double mean() { return cnt ? total() / cnt : 0.; }
  unsigned count() { return cnt; }
};
  
  Chrono interpolate_face_chrono;
#endif

  static bool noisy = false;
  void getfem_mesh_level_set_noisy(void) { noisy = true; }


  void mesh_level_set::receipt(const MESH_CLEAR &)
  { clear(); is_adapted_ = false; }
  void mesh_level_set::receipt(const MESH_DELETE &) {
    clear(); is_adapted_ = false;
    sup_sender(linked_mesh_->lmsg_sender());
  }
  void mesh_level_set::receipt(const MESH_SUP_CONVEX &m) { 
    cut_cv.erase(m.icv);
  }
  void mesh_level_set::receipt(const MESH_SWAP_CONVEX &) { 
    is_adapted_ = false;
  }

  void mesh_level_set::clear(void) {
    cut_cv.clear();
    is_adapted_ = false; touch();
  }


  const dal::bit_vector &mesh_level_set::crack_tip_convexes() const {
    return crack_tip_convexes_;
  }


  mesh_level_set::mesh_level_set(mesh &me) {
    linked_mesh_ = &me; is_adapted_ = false; 
    this->add_dependency(me);
    add_sender(me.lmsg_sender(), *this,
	       lmsg::mask(MESH_CLEAR::ID) | lmsg::mask(MESH_SUP_CONVEX::ID) |
	       lmsg::mask(MESH_SWAP_CONVEX::ID) | lmsg::mask(MESH_DELETE::ID));
  }

  mesh_level_set::~mesh_level_set() {}

  static void interpolate_face(mesh &m, dal::bit_vector& ptdone, 
			       const std::vector<size_type>& ipts,
			       bgeot::pconvex_structure cvs, 
			       size_type nb_vertices,
			       const std::vector<dal::bit_vector> &constraints,
			       const std::vector<const
			       mesher_signed_distance *> &list_constraints) {
    if (cvs->dim() == 0) return;
    else if (cvs->dim() > 1) {
      std::vector<size_type> fpts;
      for (size_type f=0; f < cvs->nb_faces(); ++f) {
	fpts.resize(cvs->nb_points_of_face(f));
	for (size_type k=0; k < fpts.size(); ++k)
	  fpts[k] = ipts[cvs->ind_points_of_face(f)[k]];
	interpolate_face(m,ptdone,fpts,cvs->faces_structure()[f], nb_vertices,
			 constraints, list_constraints);
      }
    }
    dal::bit_vector cts; size_type cnt = 0;
    for (size_type i=0; i < ipts.size(); ++i) {
      if (ipts[i] < nb_vertices) { 
	if (cnt == 0) cts = constraints[ipts[i]];
	else cts &= constraints[ipts[i]];
	++cnt;
      }
    }
    if (cts.card()) {
      // dal::bit_vector new_cts;
      for (size_type i=0; i < ipts.size(); ++i) {
	if (ipts[i] >= nb_vertices && !ptdone[ipts[i]]) { 
	  base_node &P = m.points()[ipts[i]];
	  // if (cts.card() > 1)
	  //   cout << "WARNING, projection sur " << cts << endl;
	  pure_multi_constraint_projection(list_constraints, P, cts);
	  // dist(P, new_cts);
	}
      }
    }
  }


  struct point_stock {
    
    dal::dynamic_tree_sorted
    <base_node, dal::lexicographical_less<base_node,
		dal::approx_less<scalar_type> > > points;
    std::vector<dal::bit_vector> constraints_;
    std::vector<scalar_type> radius_;
    const std::vector<const mesher_signed_distance*> &list_constraints;
    
    void clear(void) { points.clear(); constraints_.clear();radius_.clear(); }
    scalar_type radius(size_type i) const { return radius_[i]; }
    const dal::bit_vector &constraints(size_type i) const
    { return constraints_[i]; }
    size_type size(void) const { return points.size(); }
    const base_node &operator[](size_type i) const { return points[i]; }

    size_type add(const base_node &pt, const dal::bit_vector &bv,
		    scalar_type r) {
      size_type j = points.search(pt);
      if (j == size_type(-1)) {
	j = points.add(pt);
	constraints_.push_back(bv);
	radius_.push_back(r);
      }
      return j;
    }
    size_type add(const base_node &pt, scalar_type r) {
      size_type j = points.search(pt);
      if (j == size_type(-1)) {
	dal::bit_vector bv;
	for (size_type i = 0; i < list_constraints.size(); ++i)
	  if (gmm::abs((*(list_constraints[i]))(pt)) < 1E-8) bv.add(i);
	j = points.add(pt);
	constraints_.push_back(bv);
	radius_.push_back(r);
      }
      return j;
    }
    size_type add(const base_node &pt) {
      size_type j = points.search(pt);
      if (j == size_type(-1)) {
	dal::bit_vector bv;
	for (size_type i = 0; i < list_constraints.size(); ++i)
	  if (gmm::abs((*(list_constraints[i]))(pt)) < 1E-8) bv.add(i);
	j = points.add(pt);
	constraints_.push_back(bv);
	scalar_type r = min_curvature_radius_estimate(list_constraints,pt,bv);
	radius_.push_back(r);
      }
      return j;
    }

    dal::bit_vector decimate(const mesher_signed_distance &ref_element,
			     scalar_type dmin) const {
      dal::bit_vector retained_points;
      if (size() != 0) {
	size_type n = points[0].size();

	bgeot::kdtree points_tree;
	points_tree.reserve(size());
	for (size_type i = 0; i < size(); ++i)
	  points_tree.add_point_with_id(points[i], i);

	for (size_type nb_co = n; nb_co != size_type(-1); nb_co --) {
	  for (size_type i = 0; i < size(); ++i) {
	    if (!(retained_points.is_in(i)) &&
		(constraints(i).card() == nb_co ||
		 (nb_co == n && constraints(i).card() > n)) &&
		ref_element(points[i]) < 1E-8) {
	      bool kept = true;
	      scalar_type d = (nb_co == 0) ? (dmin * 1.5)
		: std::min(radius(i)*0.25, dmin);
	      base_node min = points[i], max = min;
	      for (size_type m = 0; m < n; ++m) { min[m]-=d; max[m]+=d; }
	      bgeot::kdtree_tab_type inpts;
	      points_tree.points_in_box(inpts, min, max);
	      for (size_type j = 0; j < inpts.size(); ++j)
		if (retained_points.is_in(inpts[j].i) &&
		    constraints(inpts[j].i).contains(constraints(i))
		    && gmm::vect_dist2(points[i], points[inpts[j].i]) < d)
		  { kept = false; break; }
	      if (kept) {
		if (noisy) cout << "kept point : " << points[i] << " co = "
				<< constraints(i) << endl;
		retained_points.add(i);
	      }
	    }
	  }
	}
      }
      return retained_points;
    }

    point_stock(const std::vector<const mesher_signed_distance*> &ls)
      : points(dal::lexicographical_less<base_node,
	       dal::approx_less<scalar_type> >(1E-7)), list_constraints(ls) {}
  };


  static mesher_signed_distance *new_ref_element(bgeot::pgeometric_trans pgt) {
    size_type n = pgt->structure()->dim();
    size_type nbp = pgt->basic_structure()->nb_points();
    /* Identifying simplexes.                                          */
    if (nbp == n+1 && pgt->basic_structure() == bgeot::simplex_structure(n)) {
	return new mesher_simplex_ref(n);
    }
    
    /* Identifying parallelepiped.                                     */
    if (nbp == (size_type(1) << n) &&
	pgt->basic_structure() == bgeot::parallelepiped_structure(n)) {
      base_node rmin(n), rmax(n);
      std::fill(rmax.begin(), rmax.end(), scalar_type(1));
      return new mesher_rectangle(rmin, rmax);
    }
    
    /* Identifying prisms.                                             */
    if (nbp == 2 * n && pgt->basic_structure() == bgeot::prism_structure(n)) {
      return new mesher_prism_ref(n);
    }
    
    DAL_THROW(to_be_done_error,
	      "This element is not taken into account. Contact us");
  }


  struct global_mesh_for_mesh_level_set : public mesh {};
  static mesh& global_mesh() {
    return dal::singleton<global_mesh_for_mesh_level_set>::instance();
  }

  void mesh_level_set::run_delaunay(std::vector<base_node> &fixed_points,
				    gmm::dense_matrix<size_type> &simplexes,
				    std::vector<dal::bit_vector> &
				    /* fixed_points_constraints */) {
    double t0=dal::uclock_sec();
    if (noisy) cout << "running delaunay with " << fixed_points.size()
		    << " points.." << std::flush;
    delaunay(fixed_points, simplexes);
    if (noisy) cout << " -> " << gmm::mat_ncols(simplexes)
		    << " simplexes [" << dal::uclock_sec()-t0 << "sec]\n";
  }

  static bool intersects(const mesh_level_set::zone &z1, 
			 const mesh_level_set::zone &z2) {
    for (std::set<const std::string *>::const_iterator it = z1.begin(); it != z1.end();
	 ++it)
      if (z2.find(*it) != z2.end()) return true;
    return false;
  }

  void mesh_level_set::merge_zoneset(zoneset &zones1,
				     const zoneset &zones2) const {
    for (std::set<const zone *>::const_iterator it2 = zones2.begin(); 
	 it2 != zones2.end(); ++it2) {
      zone z = *(*it2);
      for (std::set<const zone *>::iterator it1 = zones1.begin(); 
	   it1 != zones1.end(); ) {
	if (intersects(z, *(*it1))) {
	  z.insert((*it1)->begin(), (*it1)->end());
	  zones1.erase(it1++);
	}
	else ++it1;
      }
      zones1.insert(&(*(allzones.insert(z).first)));
    }
  }

  /* recursively replace '0' by '+' and '-', and the add the new zones */
  static void add_sub_zones_no_zero(std::string &s, 
				    mesh_level_set::zone &z, 
				    std::set<std::string> &allsubzones) {
    size_t i = s.find('0');
    if (i != size_t(-1)) {
      s[i] = '+'; add_sub_zones_no_zero(s, z, allsubzones);
      s[i] = '-'; add_sub_zones_no_zero(s, z, allsubzones);
    } else {
      z.insert(&(*(allsubzones.insert(s).first)));
    }
  }

  void mesh_level_set::merge_zoneset(zoneset &zones1,
				     const std::string &subz) const {
    // very sub-optimal
    zone z; std::string s(subz);
    add_sub_zones_no_zero(s, z, allsubzones);
    zoneset zs;
    zs.insert(&(*(allzones.insert(z).first)));
    merge_zoneset(zones1, zs);
  }

  /* prezone was filled for the whole convex by find_crossing_level_set. 
     This information is now refined for each sub-convex.
  */
  void mesh_level_set::find_zones_of_element(size_type cv,
					     std::string &prezone) {
    convex_info &cvi = cut_cv[cv];
    cvi.zones.clear();
    for (dal::bv_visitor i(cvi.pmsh->convex_index()); !i.finished();++i) {
      std::string subz = prezone;
      //cout << "prezone for convex " << cv << " : " << subz << endl;
      for (size_type j = 0; j < level_sets.size(); ++j) {
	if (subz[j] == '*' || subz[j] == '0') {
	  int s = sub_simplex_is_not_crossed_by(cv, level_sets[j], i);
	  subz[j] = (s < 0) ? '-' : ((s > 0) ? '+' : '0');
	}
      }
      merge_zoneset(cvi.zones, subz);
    }
    if (noisy) cout << "Number of zones for convex " << cv << " : "
		    << cvi.zones.size() << endl;
  }


  void mesh_level_set::cut_element(size_type cv,
				   const dal::bit_vector &primary,
				   const dal::bit_vector &secondary) {
    
    cut_cv[cv] = convex_info();
    cut_cv[cv].pmsh = pmesh(new mesh);
    if (noisy) cout << "cutting element " << cv << endl;
    bgeot::pgeometric_trans pgt = linked_mesh().trans_of_convex(cv);
    std::auto_ptr<mesher_signed_distance> ref_element(new_ref_element(pgt));
    std::vector<mesher_level_set> mesher_level_sets;
    
    size_type n = pgt->structure()->dim();
    size_type nbtotls = primary.card() + secondary.card();
    pintegration_method exactint = classical_exact_im(pgt);
    if (nbtotls > 16)
      DAL_THROW(failure_error,
		"An element is cut by more than 16 level_set, aborting");

		/* 
     * Step 1 : build the signed distances, estimate the curvature radius
     *          and the degree K.
     */
    dim_type K = 0; // Max. degree of the level sets.
    scalar_type r0 = 1E+10; // min curvature radius
    std::vector<const mesher_signed_distance*> list_constraints;
    point_stock mesh_points(list_constraints);

    ref_element->register_constraints(list_constraints);
    size_type nbeltconstraints = list_constraints.size();
    mesher_level_sets.reserve(nbtotls);
    for (size_type ll = 0; ll < level_sets.size(); ++ll) {
      if (primary[ll]) {
	base_node X(n);
	K = std::max(K, (level_sets[ll])->degree());
	mesher_level_sets.push_back(level_sets[ll]->mls_of_convex(cv, 0));
	mesher_level_set &mls(mesher_level_sets.back());
	list_constraints.push_back(&mesher_level_sets.back());
	r0 = std::min(r0, curvature_radius_estimate(mls, X, true));
	if (secondary[ll]) {
	  mesher_level_sets.push_back(level_sets[ll]->mls_of_convex(cv, 1));
	  mesher_level_set &mls2(mesher_level_sets.back());
	  list_constraints.push_back(&mesher_level_sets.back());
	  r0 = std::min(r0, curvature_radius_estimate(mls2, X, true));
	}
      }
    }
    
    /* 
     * Step 2 : projection of points of a grid on each signed distance.
     */
    scalar_type h0 = std::min(1./(K+1), 0.2 * r0), dmin = 0.55;
    bool h0_is_ok = true;

    do {
      h0_is_ok = true;
      mesh_points.clear();
      scalar_type geps = .001*h0; 
      size_type nbpt = 1;
      std::vector<size_type> gridnx(n);
      for (size_type i=0; i < n; ++i)
	{ gridnx[i]=1+(size_type)(1./h0); nbpt *= gridnx[i]; }
      base_node P(n), Q(n), V(n);
      dal::bit_vector co;
      std::vector<size_type> co_v;

      for (size_type i=0; i < nbpt; ++i) {
	for (size_type k=0, r = i; k < n; ++k) { // building grid point
	  unsigned p =  r % gridnx[k];
	  P[k] = p * 1. / (gridnx[k]-1);
	  r /= gridnx[k];
	}
	co.clear(); co_v.resize(0);
	if ((*ref_element)(P) < geps) {
	  bool kept = false;
	  for (size_type k=list_constraints.size()-1; k!=size_type(-1); --k) {
	    // Rerversed to project on level set before on convex boundaries
	    gmm::copy(P, Q);
	    scalar_type d = list_constraints[k]->grad(P, V);
	    if (gmm::vect_norm2(V)*h0*7 > gmm::abs(d))
	      if (try_projection(*(list_constraints[k]), Q, true)) {
		if (k >= nbeltconstraints
		    && gmm::vect_dist2(P, Q) < h0 * 3.5) kept = true;
		if (gmm::vect_dist2(P, Q) < h0 / 1.5) {
		  co.add(k); co_v.push_back(k); 
		  if ((*ref_element)(Q) < 1E-8) {
		    scalar_type r =
		      curvature_radius_estimate(*(list_constraints[k]), Q);
		    r0 = std::min(r0, r);
		    if (r0 < 4.*h0) { h0 = 0.2 * r0; h0_is_ok = false; break; }
		    if (k >= nbeltconstraints || kept) mesh_points.add(Q, r);
		  }
		}
	      }
	  }
	  if (kept) mesh_points.add(P, 1.E10);
	}
	size_type nb_co = co.card();
	dal::bit_vector nn;
	if (h0_is_ok)
	  for (size_type count = 0; count < (size_type(1) << nb_co); ++count) {
	    nn.clear();
	    for (size_type j = 0; j < nb_co; ++j)
	      if (count & (size_type(1) << j)) nn.add(co_v[j]);
	    if (nn.card() > 1 && nn.card() <= n) {
	      if (noisy) cout << "trying set of constraints" << nn << endl;
	      gmm::copy(P, Q);
	      bool ok=pure_multi_constraint_projection(list_constraints,Q,nn);
	      if (ok && (*ref_element)(Q) < 1E-8) {
		if (noisy) cout << "Intersection point found " << Q << " with "
		     << nn << endl;
		mesh_points.add(Q);
	      }
	    }
	  }
      }
      
      if (!h0_is_ok) continue;
      
    /* 
     * Step 3 : Delaunay, test if a simplex cut a level set and adapt
     *          the mesh to the curved level sets.
     */

      dmin = 2.*h0;
      if (noisy) cout << "dmin = " << dmin << " h0 = " << h0 << endl;
      if (noisy) cout << "cetait le convexe " << cv << endl;
      if (noisy) getchar(); 
      
      dal::bit_vector retained_points
	= mesh_points.decimate(*ref_element, dmin);
      
    delaunay_again :
      std::vector<base_node> fixed_points;
      std::vector<dal::bit_vector> fixed_points_constraints;
      fixed_points.reserve(retained_points.card());
      fixed_points_constraints.reserve(retained_points.card());
      for (dal::bv_visitor i(retained_points); !i.finished(); ++i) {
	fixed_points.push_back(mesh_points[i]);
	fixed_points_constraints.push_back(mesh_points.constraints(i));
      }

      gmm::dense_matrix<size_type> simplexes;
      run_delaunay(fixed_points, simplexes, fixed_points_constraints);

      mesh &msh(*(cut_cv[cv].pmsh));

      mesh_region &ls_border_faces(cut_cv[cv].ls_border_faces);

      msh.clear(); 

      for (size_type i = 0; i <  fixed_points.size(); ++i) {
	size_type j = msh.add_point(fixed_points[i], false);
	assert(j == i);
      }

      std::vector<base_node> cvpts;
      for (size_type i = 0; i < gmm::mat_ncols(simplexes); ++i) {
	size_type j = msh.add_convex(bgeot::simplex_geotrans(n,1),
				      gmm::vect_const_begin(gmm::mat_col(simplexes, i)));
	/* remove flat convexes => beware the final mesh may not be conformal ! */
 	if (msh.convex_quality_estimate(j) < 1E-12) msh.sup_convex(j);
 	else {
	  std::vector<scalar_type> signs(list_constraints.size());
	  std::vector<size_type> prev_point(list_constraints.size());
	  for (size_type ii = 0; ii <= n; ++ii) {
	    for (size_type jj = 0; jj < list_constraints.size(); ++jj) {
	      scalar_type dd =
		(*(list_constraints[jj]))(msh.points_of_convex(j)[ii]);
	      if (gmm::abs(dd) > 1E-7) {
		if (dd * signs[jj] < 0.0) {
		  if (noisy) cout << "Intersection trouvee ... \n";
		  // calcul d'intersection
		  base_node X = msh.points_of_convex(j)[ii], G;
		  base_node VV = msh.points_of_convex(j)[prev_point[jj]] - X;
		  if (dd > 0.) gmm::scale(VV, -1.);
		  dd = (*(list_constraints[jj])).grad(X, G);
		  size_type nbit = 0;
		  while (gmm::abs(dd) > 1e-15) {
		    if (++nbit > 1000) {
		      if (noisy) cout << "Intersection not found";
		      assert(false);
		    }
		    scalar_type nG = std::max(1E-8, gmm::vect_sp(G, VV));
		    gmm::add(gmm::scaled(VV, -dd / nG), X);
		    dd = (*(list_constraints[jj])).grad(X, G);
		  }
		  size_type kk = mesh_points.add(X);
		  if (!(retained_points[kk])) {
		    retained_points.add(kk);
		    goto delaunay_again;
		  }
		}
		if (signs[jj] == 0.0) { signs[jj] = dd; prev_point[jj] = ii; }
	      }
	    }
	  }
	}
      }

      if (K>1) {
	for (dal::bv_visitor_c j(msh.convex_index()); !j.finished(); ++j) {
	  bgeot::pgeometric_trans pgt2 = bgeot::simplex_geotrans(n, K);
	  cvpts.resize(pgt2->nb_points());
	  for (size_type k=0; k < pgt2->nb_points(); ++k) {
	    cvpts[k] = bgeot::simplex_geotrans(n,1)->transform
	      (pgt2->convex_ref()->points()[k], msh.points_of_convex(j));
	  }
	  msh.sup_convex(j);
	  /* some new point will be added to the mesh, but the initial ones (those on the
	     convex vertices) are kepts */
	  size_type jj = msh.add_convex_by_points(pgt2, cvpts.begin());
	  assert(jj == j);
	}
      }

      if (noisy) {
	getfem::stored_mesh_slice sl;
	sl.build(msh, getfem::slicer_none(), 1);
	char s[512]; sprintf(s, "totobefore%d.dx", int(cv));
	getfem::dx_export exp(s);
	exp.exporting(sl);
	exp.exporting_mesh_edges();
	exp.write_mesh();
      }

      /* detect the faces lying on level_set boundaries
	 (not exact, some other faces maybe be found, use this only for vizualistion) */
      for (dal::bv_visitor i(msh.convex_index()); !i.finished(); ++i) {
	for (size_type f = 0; f <= n; ++f) {
	  const mesh::ind_pt_face_ct &fpts
	    = msh.ind_points_of_face_of_convex(i, f);
	  
	  dal::bit_vector cts; bool first = true;
	  for (unsigned k=0; k < fpts.size(); ++k) {
	    if (fpts[k] >= fixed_points_constraints.size()) {
	      assert(K>1);
	      continue; /* ignore additional point inserted when K>1 */
	    }
	    if (first) { cts = fixed_points_constraints[fpts[k]]; first = false; }
	    else cts &= fixed_points_constraints[fpts[k]];
	  }
	  for (dal::bv_visitor ii(cts); !ii.finished(); ++ii) {
	    if (ii >= nbeltconstraints)
	      ls_border_faces.add(i, f);
	  }
	}
      }

      
      if (K > 1) { // à ne faire que sur les convexes concernés ..
	dal::bit_vector ptdone;
	std::vector<size_type> ipts;
	for (dal::bv_visitor i(msh.convex_index()); !i.finished(); ++i) {
	  for (size_type f = 0; f <= n; ++f) {
	    const mesh::ind_pt_face_ct &fpts
	      = msh.ind_points_of_face_of_convex(i, f);
	    ipts.assign(fpts.begin(), fpts.end());
#ifdef DEBUG_LS
	    interpolate_face_chrono.tic();
#endif
	    
	    interpolate_face(msh, ptdone, ipts,
			     msh.trans_of_convex(i)->structure()
			     ->faces_structure()[f], fixed_points.size(),
			     fixed_points_constraints, list_constraints);
#ifdef DEBUG_LS
	    interpolate_face_chrono.toc();
#endif
	  }
	}
      }

      if (noisy) {
	getfem::stored_mesh_slice sl;
	sl.build(msh, getfem::slicer_none(), 12);
	char s[512]; sprintf(s, "toto%d.dx", int(cv));
	getfem::dx_export exp(s);
	exp.exporting(sl);
	exp.exporting_mesh_edges();
	exp.write_mesh();
      }
    
      /* 
       * Step 4 : Building the new integration method.
       */
      base_matrix G;
      bgeot::pgeometric_trans pgt2 = bgeot::simplex_geotrans(n, K);
      papprox_integration pai = classical_approx_im(pgt2,2*K)->approx_method();
      approx_integration new_approx(pgt->convex_ref());
      base_matrix KK(n,n), CS(n,n);
      base_matrix pc(pgt2->nb_points(), n); 
      for (dal::bv_visitor i(msh.convex_index()); !i.finished(); ++i) {
	vectors_to_base_matrix(G, msh.points_of_convex(i));
	bgeot::geotrans_interpolation_context c(msh.trans_of_convex(i),
						pai->point(0), G);
	scalar_type sign = 0.0;
	for (size_type j = 0; j < pai->nb_points_on_convex(); ++j) {
	  c.set_xref(pai->point(j));
	  pgt2->poly_vector_grad(pai->point(j), pc);
	  gmm::mult(G,pc,KK);
	  scalar_type J = gmm::lu_det(KK);
	  if (noisy && J * sign < 0)
	    cout << "ATTENTION retrournement de situation en convexe " << i
		 << "sign = " << sign << " J = " << J << " p1 = "
		 << msh.points_of_convex(i)[0] << " p2 = "
		 << msh.points_of_convex(i)[1] << " p3 = "
		 << msh.points_of_convex(i)[2] << " p4 = "
		 << msh.points_of_convex(i)[3] << " p5 = "
		 << msh.points_of_convex(i)[4] << " p6 = "
		 << msh.points_of_convex(i)[5] << " K = " << int(K) << endl;
	  if (sign == 0 && gmm::abs(J) > 1E-14) sign = J;
	  new_approx.add_point(c.xreal(), pai->coeff(j) * gmm::abs(J));
	}
      }

      getfem::mesh_region border_faces;
      getfem::outer_faces_of_mesh(msh, border_faces);
      for (getfem::mr_visitor it(border_faces); !it.finished(); ++it) {
	vectors_to_base_matrix(G, msh.points_of_convex(it.cv()));
	bgeot::geotrans_interpolation_context c(msh.trans_of_convex(it.cv()),
						pai->point(0), G);
	for (size_type j = 0; j < pai->nb_points_on_face(it.f()); ++j) {
	  c.set_xref(pai->point_on_face(it.f(), j));
	  new_approx.add_point(c.xreal(), pai->coeff_on_face(it.f(), j)
			       * gmm::abs(c.J()), it.f() /* faux */);
	}
      }
      new_approx.valid_method();

      /* 
       * Step 5 : Test the validity of produced integration method.
       */
      
      scalar_type error = test_integration_error(&new_approx, 1);
      if (noisy) cout << " max monomial integration err: " << error << "\n";
      if (error > 1e-5) {
	if (noisy) cout << "PAS BON NON PLUS\n"; if (noisy) getchar();
	if (dmin > 3*h0) { dmin /= 2.; }
	else { h0 /= 2.0; dmin = 2.*h0; }
	h0_is_ok = false;
      }
 

      if (h0_is_ok && noisy) { // ajout dans global mesh pour visu
	vectors_to_base_matrix(G, linked_mesh().points_of_convex(cv));
	std::vector<size_type> pts(msh.nb_points());
	for (size_type i = 0; i < msh.nb_points(); ++i)
	  pts[i] = global_mesh().add_point(pgt->transform(msh.points()[i], G));
	
	for (dal::bv_visitor i(msh.convex_index()); !i.finished(); ++i)
	  global_mesh().add_convex(msh.trans_of_convex(i), 
				 dal::index_ref_iterator(pts.begin(),
				     msh.ind_points_of_convex(i).begin()));
      }

    } while (!h0_is_ok);

#ifdef DEBUG_LS
    cout << "Interpolate face: " << interpolate_face_chrono.total()
	 << " moyenne: " << interpolate_face_chrono.mean() << "\n";
#endif
  }


  void mesh_level_set::global_cut_mesh(mesh &m) const {
    m.clear();
    base_matrix G;
    for (dal::bv_visitor cv(linked_mesh().convex_index()); 
	 !cv.finished(); ++cv) {
      if (is_convex_cut(cv)) {
	const convex_info &ci = (cut_cv.find(cv))->second;
	mesh &msh = *(ci.pmsh);
	bgeot::pgeometric_trans pgt = linked_mesh().trans_of_convex(cv);
	vectors_to_base_matrix(G, linked_mesh().points_of_convex(cv));
	std::vector<size_type> pts(msh.nb_points());
	for (size_type i = 0; i < msh.nb_points(); ++i)
	  pts[i] = m.add_point(pgt->transform(msh.points()[i], G));
	
	std::vector<size_type> ic2(msh.convex_index().last_true()+1);

	for (dal::bv_visitor i(msh.convex_index()); !i.finished(); ++i) {
	  ic2[i] = m.add_convex(msh.trans_of_convex(i), 
				dal::index_ref_iterator(pts.begin(),
			        msh.ind_points_of_convex(i).begin()));
	  
	}
	for (mr_visitor i(ci.ls_border_faces); !i.finished(); ++i) {
	  m.region(0).add(ic2[i.cv()], i.f());
	}

      } else {
	m.add_convex_by_points(linked_mesh().trans_of_convex(cv),
			       linked_mesh().points_of_convex(cv).begin());
      }
    }


  }

  void mesh_level_set::update_crack_tip_convexes() {
    crack_tip_convexes_.clear();
     
    for (std::map<size_type, convex_info>::const_iterator it = cut_cv.begin(); 
	 it != cut_cv.end(); ++it) {
      size_type cv = it->first;
      mesh &msh = *(it->second.pmsh);      
      for (unsigned ils = 0; ils < nb_level_sets(); ++ils) {
	if (get_level_set(ils)->has_secondary()) {
	  mesher_level_set mesherls0 =  get_level_set(ils)->mls_of_convex(cv, 0, false);
	  mesher_level_set mesherls1 =  get_level_set(ils)->mls_of_convex(cv, 1, false);
	  for (dal::bv_visitor ii(msh.convex_index()); !ii.finished(); ++ii) {
	    for (unsigned ipt = 0; ipt < msh.nb_points_of_convex(ii); ++ipt) {
	      if (gmm::abs(mesherls0(msh.points_of_convex(ii)[ipt])) < 1E-10
		  && gmm::abs(mesherls1(msh.points_of_convex(ii)[ipt])) < 1E-10) {
		crack_tip_convexes_.add(cv);
		goto next_convex;
	      }
	    }
	  }
	}
      }
    next_convex:
      ;
    }    
  }

  void mesh_level_set::adapt(void) {

    // compute the elements touched by each level set
    // for each element touched, compute the sub mesh
    //   then compute the adapted integration method
    cut_cv.clear();
    allsubzones.clear();
    zones_of_convexes.clear();
    allzones.clear();

    std::string z;
    for (dal::bv_visitor cv(linked_mesh().convex_index()); 
	 !cv.finished(); ++cv) {
      dal::bit_vector prim, sec;
      find_crossing_level_set(cv, prim, sec, z);
      zones_of_convexes[cv] = &(*(allsubzones.insert(z).first));
      if (noisy) cout << "element " << cv << " cut level sets : "
		      << prim << " zone : " << z << endl;
      if (prim.card()) {
	cut_element(cv, prim, sec);
	find_zones_of_element(cv, z);
      }
    }
    if (noisy) {
      getfem::stored_mesh_slice sl;
      sl.build(global_mesh(), getfem::slicer_none(), 6);
      getfem::dx_export exp("totoglob.dx");
      exp.exporting(sl);
      exp.exporting_mesh_edges();
      exp.write_mesh();
    }

    update_crack_tip_convexes();
    is_adapted_ = true;
  }
  
  int mesh_level_set::sub_simplex_is_not_crossed_by(size_type cv,
						    plevel_set ls,
						    size_type sub_cv) {
    scalar_type EPS = 1e-8;
    bgeot::pgeometric_trans pgt = linked_mesh().trans_of_convex(cv);
    convex_info &cvi = cut_cv[cv];
    bgeot::pgeometric_trans pgt2 = cvi.pmsh->trans_of_convex(sub_cv);

    mesher_level_set mls0 = ls->mls_of_convex(cv, 0), mls1(mls0);
    if (ls->has_secondary()) mls1 = ls->mls_of_convex(cv, 1);
    int p = 0;
    bool cutted = false;
    scalar_type d2 = 0, d1 = 1, d0 = 0, d0min = 0;
    for (size_type i = 0; i < pgt2->nb_points(); ++i) {
      d0 = mls0(cvi.pmsh->points_of_convex(sub_cv)[i]);
      if (i == 0) d0min = gmm::abs(d0);
      else d0min = std::min(d0min, gmm::abs(d0));
      if (ls->has_secondary())
	d1 = std::min(d1, mls1(cvi.pmsh->points_of_convex(sub_cv)[i]));
     
      int p2 = ( (d0 < -EPS) ? -1 : ((d0 > EPS) ? +1 : 0));
      if (p == 0) p = p2;
      if (gmm::abs(d0) > gmm::abs(d2)) d2 = d0;
      if (!p2 || p*p2 < 0) cutted = true;
    }
    if (cutted && d1 > +EPS)
      { cout << "ooops, je retourne 0: d1 = " << d1 << endl; return 0; }
    if (d0min < EPS &&  d1 > -EPS) return 0;
    return (d2 < 0.) ? -1 : 1;
  }


  int mesh_level_set::is_not_crossed_by(size_type cv, plevel_set ls,
					unsigned lsnum) {
    const mesh_fem &mf = ls->get_mesh_fem();
    bgeot::pgeometric_trans pgt = linked_mesh().trans_of_convex(cv);

    const mesh_fem::ind_dof_ct &dofs = mf.ind_dof_of_element(cv);
    pfem pf = mf.fem_of_element(cv);
    int p = -2;
    mesher_level_set mls1 = ls->mls_of_convex(cv, lsnum, false);
    scalar_type EPS = 1e-8;

    /* easy cases: 
     - different sign on the dof nodes => intersection for sure
     - min value of the levelset greater than the radius of the convex
     => no intersection
    */ 
    for (mesh_fem::ind_dof_ct::const_iterator it=dofs.begin();
	 it != dofs.end(); ++it) {
      scalar_type v = ls->values(lsnum)[*it];
      int p2 = ( (v < -EPS) ? -1 : ((v > EPS) ? +1 : 0));
      if (p == -2) p=p2;
      if (!p2 || p*p2 < 0) return 0;
    }

    base_node X(pf->dim()), G(pf->dim());
    gmm::fill_random(X); X *= 1E-2;
    scalar_type d = mls1.grad(X, G);
    if (gmm::vect_norm2(G)*2. < gmm::abs(d)) return p;

    std::auto_ptr<mesher_signed_distance> ref_element(new_ref_element(pgt));
    
    gmm::fill_random(X); X *= 1E-2;
    mesher_intersection mi1(*ref_element, mls1);
    if (!try_projection(mi1, X)) return p;
    if ((*ref_element)(X) > 1E-8) return p;
    
    gmm::fill_random(X); X *= 1E-2;
    mesher_level_set mls2 = ls->mls_of_convex(cv, lsnum, true);
    mesher_intersection mi2(*ref_element, mls2);
    if (!try_projection(mi2, X)) return p;
    if ((*ref_element)(X) > 1E-8) return p;
   
    return 0;
  }

  void mesh_level_set::find_crossing_level_set(size_type cv,
					       dal::bit_vector &prim,
					       dal::bit_vector &sec,
					       std::string &z) {
    prim.clear(); sec.clear();
    z = std::string(level_sets.size(), '*');
    unsigned lsnum = 0;
    for (size_type k = 0; k < level_sets.size(); ++k, ++lsnum) {
      if (noisy) cout << "testing cv " << cv << " with level set "
		      << k << endl;
      int s = is_not_crossed_by(cv, level_sets[k], 0);
      if (!s) {
	if (noisy) cout << "is cut \n";
	if (level_sets[k]->has_secondary()) {
	  s = is_not_crossed_by(cv, level_sets[k], 1);
	  if (!s) { sec.add(lsnum); prim.add(lsnum); }
	  else if (s < 0) prim.add(lsnum); else z[k] = '0';
	}
	else prim.add(lsnum);
      }
      else z[k] = (s < 0) ? '-' : '+';
    }
  }


}  /* end of namespace getfem.                                             */



