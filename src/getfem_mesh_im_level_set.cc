// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Elements Methods (getfem)
// File    : getfem_mesh_im_level_set.cc : adaptable integration methods
//           on convex meshes.
// Date    : February 02, 2005.
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

#include <getfem_mesh_im_level_set.h>
#include <getfem_mesher.h>
#include <bgeot_kdtree.h>

namespace getfem {
  
  void mesh_im_level_set::receipt(const MESH_CLEAR &) { clear(); }
  void mesh_im_level_set::receipt(const MESH_DELETE &) { clear(); }
  void mesh_im_level_set::clear(void) { mesh_im::clear(); level_sets.clear(); }

  mesh_im_level_set::mesh_im_level_set(getfem_mesh &me,
				       pintegration_method reg,
				       pintegration_method sing) 
    : mesh_im(me), cut_im(me) {
    regular_simplex_pim = reg;
    singular_simplex_pim = (sing == 0) ? reg : sing;
  }

  pintegration_method 
  mesh_im_level_set::int_method_of_element(size_type cv) const {
    if (cut_im.convex_index().is_in(cv)) 
      return cut_im.int_method_of_element(cv); 
    else return mesh_im::int_method_of_element(cv);
  }

  
  static void interpolate_face(getfem_mesh &m, dal::bit_vector& ptdone, 
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
	      scalar_type d = (nb_co == 0) ? (dmin * 3.)
		: std::min(radius(i)*0.25, dmin);
	      // if (nb_co == n) d = dmin / 10.; // utile ?
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
		cout << "kept point : " << points[i] << " co = "
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


  static getfem_mesh global_mesh; // to visualize the result with Open dx

  void mesh_im_level_set::cut_element(size_type cv,
				      const dal::bit_vector &primary,
				      const dal::bit_vector &secondary) {
    cout << "cutting element " << cv << endl;
    std::auto_ptr<mesher_signed_distance> ref_element;
    std::vector<mesher_level_set> mesher_level_sets;
    
    bgeot::pgeometric_trans pgt = linked_mesh().trans_of_convex(cv);
    size_type n = pgt->structure()->dim();
    size_type nbp = pgt->basic_structure()->nb_points();
    size_type nbtotls = primary.card() + secondary.card();
    pintegration_method exactint;
    if (nbtotls > 16)
      DAL_THROW(failure_error,
		"An element is cutted by more than 16 level_set, aborting");
      
    /*
     * Step 1 : build the signed distance for the reference element.
     */
    unsigned found = 0;

    /* Identifying simplexes.                                          */
    if (nbp == n+1 && pgt->basic_structure() == bgeot::simplex_structure(n)) {
	ref_element.reset(new mesher_simplex_ref(n)); found = 1;
	exactint = exact_simplex_im(n);
    }
    
    /* Identifying parallelepiped.                                     */
    if (!found && nbp == (size_type(1) << n) &&
	pgt->basic_structure() == bgeot::parallelepiped_structure(n)) {
      base_node rmin(n), rmax(n);
      std::fill(rmax.begin(), rmax.end(), scalar_type(1));
      ref_element.reset(new mesher_rectangle(rmin, rmax)); found = 2;
      exactint = exact_parallelepiped_im(n);
    }

    /* Identifying prisms.                                             */
    if (!found && nbp == 2 * n &&
	pgt->basic_structure() == bgeot::prism_structure(n)) {
      ref_element.reset(new mesher_prism_ref(n)); found = 3;
      exactint = exact_prism_im(n);
    }
    
    if (!found) 
      DAL_THROW(to_be_done_error,
		"This element is not taken into account. Contact us");

    /* 
     * Step 2 : build the signed distances, estimate the curvature radius
     *          and the degree K and find points on intersections of level sets.
     */
    dim_type K = 0; // Max. degree of the level sets.
    scalar_type r0 = 1E+10; // min curvature radius
    std::vector<const mesher_signed_distance*> list_constraints, signed_dist;
    point_stock mesh_points(list_constraints);
    cout.precision(16);

    for (size_type count = 0; count < (size_type(1) << nbtotls); ++count) {
      signed_dist.clear();  signed_dist.reserve(1+nbtotls); 
      mesher_level_sets.clear(); mesher_level_sets.reserve(1+nbtotls); 
      signed_dist.push_back(ref_element.get());
      unsigned k = 0, l = 0; K = 0;
      pfem max_deg_fem(0);
      for (std::set<plevel_set>::const_iterator it = level_sets.begin();
	   it != level_sets.end(); ++it, ++l) {
	if (primary[l]) {
	  if ((*it)->degree() > K) {
	    K = (*it)->degree();
	    max_deg_fem = (*it)->get_mesh_fem().fem_of_element(cv);
	  }
	  mesher_level_sets.push_back((*it)->mls_of_convex(cv, 0,
					     count & (size_type(1) << k)));
	  signed_dist.push_back(&mesher_level_sets.back());
	  ++k;
	  if (secondary[l]) {
	    mesher_level_sets.push_back((*it)->mls_of_convex(cv, 1,
						 count & (size_type(1)<< k)));
	    signed_dist.push_back(&mesher_level_sets.back());
	    ++k;
	  }
	}
      }
      assert(k == nbtotls);

      mesher_intersection final_dist(signed_dist);
      list_constraints.clear();
      final_dist.register_constraints(list_constraints);
      base_node local_box_min(n), local_box_max(n);
      bool local_box_init = false;
      for (size_type i=0; i < max_deg_fem->node_convex(cv).nb_points();++i) {
	base_node X = max_deg_fem->node_convex(cv).points()[i];
	if (try_projection(final_dist, X)) {
	  // computation of the curvature radius
	  dal::bit_vector bv;
	  final_dist(X, bv);
	  cout << "constraints : " << bv << endl;
	  scalar_type r
	    = min_curvature_radius_estimate(list_constraints, X, bv);
	  cout << "rayon de courbure de : " << r << endl;
	  r0 = std::min(r, r0);
	  if (!local_box_init)
	    { local_box_min = local_box_max = X; local_box_init = true; }
	  else for (size_type j = 0; j < n; ++j) {
	    local_box_min[j] = std::min(local_box_min[j], X[j]);
	    local_box_max[j] = std::max(local_box_max[j], X[j]);
	  }
	  mesh_points.add(X, r);
	}
      }
      cout << "Local_box_min = " << local_box_min << " local_box_max = "
	   << local_box_max << endl;
    }
    
    /* 
     * Step 3 : projection of points of a grid on each signed distance.
     */
    scalar_type h0 = std::min(1./(K+1), 0.2 * r0), dmin = 0.55;
    bool h0_is_ok = true;

    do {
      h0_is_ok = true;
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
	  mesh_points.add(P, 1.E10);
	  for (size_type k=0; k < list_constraints.size(); ++k) {
	    gmm::copy(P, Q);
	    scalar_type d = list_constraints[k]->grad(P, V);
	    if (gmm::vect_norm2(V)*h0 > gmm::abs(d))
	      if (try_projection(*(list_constraints[k]), Q, true))
		if (gmm::vect_dist2(P, Q) < h0 / 1.5) {
		  co.add(k); co_v.push_back(k); 
		  if ((*ref_element)(Q) < 1E-8) {
		    scalar_type r =
		      curvature_radius_estimate(*(list_constraints[k]), Q);
		    r0 = std::min(r0, r);
		    if (r0 < 4.*h0) { h0 = 0.2 * r0; h0_is_ok = false; break; }
		    mesh_points.add(Q, r);
		  }
		}
	  }
	}
	size_type nb_co = co.card();
	dal::bit_vector nn;
	if (h0_is_ok)
	  for (size_type count = 0; count < (size_type(1) << nb_co); ++count) {
	    nn.clear();
	    for (size_type j = 0; j < nb_co; ++j)
	      if (count & (size_type(1) << j)) nn.add(co_v[j]);
	    if (nn.card() > 1 && nn.card() <= n) {
	      cout << "trying set of constraints" << nn << endl;
	      gmm::copy(P, Q);
	      bool ok=pure_multi_constraint_projection(list_constraints,Q,nn);
	      if (ok && (*ref_element)(Q) < 1E-8) {
		cout << "Intersection point found " << Q << " with "
		     << nn << endl;
		mesh_points.add(Q);
	      }
	    }
	  }
      }
      
      if (!h0_is_ok) continue;
      
    /* 
     * Step 4 : Delaunay, test if a simplex cut a level set and adapt
     *          the mesh to the curved level sets.
     */

      dmin = 2.*h0;
      cout << "dmin = " << dmin << " h0 = " << h0 << endl;
      cout << "cetait le convexe " << cv << endl;
      getchar(); 
      
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
      delaunay(fixed_points, simplexes);
      cout << "Nb simplexes = " << gmm::mat_ncols(simplexes) << endl;
      getfem_mesh mesh;
      for (size_type i = 0; i <  fixed_points.size(); ++i) {
	size_type j = mesh.add_point(fixed_points[i], false);
	assert(j == i);
      }
      std::vector<base_node> cvpts;
      for (size_type i = 0; i < gmm::mat_ncols(simplexes); ++i) {
	size_type j = mesh.add_convex(bgeot::simplex_geotrans(n,1),
				      gmm::vect_const_begin(gmm::mat_col(simplexes, i)));
	if (mesh.convex_quality_estimate(j) < 1E-18) mesh.sup_convex(j);
	else {
	  std::vector<scalar_type> signs(list_constraints.size());
	  std::vector<size_type> prev_point(list_constraints.size());
	  for (size_type ii = 0; ii <= n; ++ii) {
	    for (size_type jj = 0; jj < list_constraints.size(); ++jj) {
	      scalar_type dd =
		(*(list_constraints[jj]))(mesh.points_of_convex(j)[ii]);
	      if (gmm::abs(dd) > 1E-7) {
		if (dd * signs[jj] < 0.0) {
		  cout << "Intersection trouvee ... \n";
		  // calcul d'intersection
		  base_node X = mesh.points_of_convex(j)[ii], G;
		  base_node VV = mesh.points_of_convex(j)[prev_point[jj]] - X;
		  if (dd > 0.) gmm::scale(VV, -1.);
		  dd = (*(list_constraints[jj])).grad(X, G);
		  size_type nbit = 0;
		  while (gmm::abs(dd) > 1e-15) {
		    if (++nbit > 1000)
		      { cout << "Intersection not found"; assert(false);}
		    scalar_type nG = std::max(1E-8, gmm::vect_sp(G, VV));
		    gmm::add(gmm::scaled(VV, -dd / nG), X);
		    dd = (*(list_constraints[jj])).grad(X, G);
		  }
		  size_type ii = mesh_points.add(X);
		  if (!(retained_points[ii])) {
		    retained_points.add(ii);
		    goto delaunay_again;
		  }
		}
		if (signs[jj] == 0.0) { signs[jj] = dd; prev_point[jj] = ii; }
	      }
	    }
	  }
	  if (K > 1) {
	    bgeot::pgeometric_trans pgt2 = bgeot::simplex_geotrans(n, K);
	    cvpts.resize(pgt2->nb_points());
	    for (size_type k=0; k < pgt2->nb_points(); ++k) {
	      cvpts[k] = bgeot::simplex_geotrans(n,1)->transform
		(pgt2->convex_ref()->points()[k], mesh.points_of_convex(j));
	    }
	    mesh.sup_convex(j);
	    mesh.add_convex_by_points(pgt2, cvpts.begin());
	  }
	}
      }
      
      cout << "storing mesh" << endl;
      
      {
	getfem::stored_mesh_slice sl;
	sl.build(mesh, getfem::slicer_none(), 1); //getfem::slicer_noneexplode(0.8), 8);
	char s[512]; sprintf(s, "totobefore%d.dx", cv);
	getfem::dx_export exp(s);
	exp.exporting(sl);
	exp.exporting_mesh_edges();
	exp.write_mesh();
      }
      
      cout << "end of storing mesh" << endl;
      
      if (K > 1) { // à ne faire que sur les convexes concernés ..
	dal::bit_vector ptdone;
	std::vector<size_type> ipts;
	for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i) {
	  for (size_type f = 0; f <= n; ++f) {
	    bgeot::ind_ref_mesh_point_ind_ct fpts
	      = mesh.ind_points_of_face_of_convex(i, f);
	    ipts.assign(fpts.begin(), fpts.end());
	    interpolate_face(mesh, ptdone, ipts,
			     mesh.trans_of_convex(i)->structure()->faces_structure()[f],
			     fixed_points.size(), fixed_points_constraints, list_constraints);
	  }
	}
      }
      
      {
	getfem::stored_mesh_slice sl;
	sl.build(mesh, getfem::slicer_none(), 12); //getfem::slicer_noneexplode(0.8), 8);
	char s[512]; sprintf(s, "toto%d.dx", cv);
	getfem::dx_export exp(s);
	exp.exporting(sl);
	exp.exporting_mesh_edges();
	exp.write_mesh();
      }
      
      std::vector<base_node> gauss_points;
      std::vector<scalar_type> gauss_weights;
      for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i) {
	base_matrix G;
	vectors_to_base_matrix(G, mesh.points_of_convex(i));
	
	papprox_integration pai = regular_simplex_pim->approx_method();
	for (size_type j = 0; j < pai->nb_points_on_convex(); ++j) {
	  bgeot::geotrans_interpolation_context c(mesh.trans_of_convex(i), pai->point(j), G); // sortir de la boucle et faire des chgt de points
	  // base_matrix B = c.B(); // for J to be computed;
	  // cout << "adding Gauss point " << c.xreal() << " with weight " << pai->coeff(j) * c.J() << "(J = " << c.J() << ")\n";
	  gauss_points.push_back(c.xreal());
	  gauss_weights.push_back(pai->coeff(j) * gmm::abs(c.J()));
	}
      }
      
      // + test ...
      scalar_type wtot(0);
      cout << "Number of gauss points : " << gauss_weights.size();
      for (size_type k = 0; k < gauss_weights.size(); ++k) wtot += gauss_weights[k];
      base_poly poly = bgeot::one_poly(n);
      scalar_type exactvalue = exactint->exact_method()->int_poly(poly);
      cout.precision(16);
      cout << "The Result : " << wtot << " compared to " << exactvalue << endl; 
      
      if (gmm::abs(wtot-exactvalue) > 1E-7){
	cout << "PAS BON NON PLUS\n"; getchar();
	if (dmin > 3*h0) { dmin /= 2.; }
	else { h0 /= 2.0; dmin = 2.*h0; }
	h0_is_ok = false;
      }
 
      if (h0_is_ok) { // ajout dans global mesh pour visu
	base_matrix G;
	vectors_to_base_matrix(G, linked_mesh().points_of_convex(cv));
	std::vector<size_type> pts(mesh.nb_points());
	for (size_type i = 0; i < mesh.nb_points(); ++i)
	  pts[i] = global_mesh.add_point(pgt->transform(mesh.points()[i], G));
	
	for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i)
	  global_mesh.add_convex(mesh.trans_of_convex(i), 
				 dal::index_ref_iterator(pts.begin(),
				     mesh.ind_points_of_convex(i).begin()));
      }

    } while (!h0_is_ok);
  }


  void mesh_im_level_set::adapt(void) {

    // compute the elements touched by each level set
    // for each element touched, compute the sub mesh
    //   then compute the adapted integration method
    for (dal::bv_visitor cv(linked_mesh().convex_index()); 
	 !cv.finished(); ++cv) {
      dal::bit_vector prim, sec;
      find_crossing_level_set(cv, prim, sec);
      cout << "element " << cv << " cutted level sets : " << prim << endl;
      if (prim.card()) cut_element(cv, prim, sec);
    }
    getfem::stored_mesh_slice sl;
    sl.build(global_mesh, getfem::slicer_none(), 6); //getfem::slicer_noneexplode(0.8), 8);
    getfem::dx_export exp("totoglob.dx");
    exp.exporting(sl);
    exp.exporting_mesh_edges();
    exp.write_mesh();

  }

  bool mesh_im_level_set::is_crossed_by(size_type cv, plevel_set ls, unsigned lsnum) {
    const mesh_fem &mf = ls->get_mesh_fem();
    ref_mesh_dof_ind_ct dofs = mf.ind_dof_of_element(cv);
    base_vector v; v.reserve(dofs.size());
    int p = -2;
    base_node G = dal::mean_value(linked_mesh().points_of_convex(cv));
    scalar_type radius2 = 0, dmin = 1e30;

    scalar_type EPS = 1e-8;

    /* easy cases: 
     - different sign on the dof nodes => intersection for sure
     - min value of the levelset greater than the radius of the convex
     => no intersection (assuming the level set is locally a signed
     distance, and the geotrans of the convex is not too strong..)
    */ 
    for (ref_mesh_dof_ind_ct::const_iterator it = dofs.begin();
	 it != dofs.end(); ++it) {
      radius2 = std::max(radius2, 
			 gmm::vect_dist2_sqr(G, mf.point_of_dof(*it)));
      v.push_back(ls->values()[*it]);
      int p2 = ( (v.back() < -EPS) ? -1 :
		 ((v.back() > EPS) ? +1 : 0));
      if (p == -2) p=p2; else
	if (p*p2 < 0) return true;
      // il faudrait obtenir le gradient ... dist.grad(X, grad);
      // dmin = std::min(dmin, gmm::abs(v.back() / gmm::vect_norm2(grad)));
      dmin = std::min(dmin, gmm::abs(v.back()));
    }
    cout << "v = " << v << endl;
    cout << "G = " << G << endl;
    cout << "dmin = " << dmin << endl;
    cout << "radius2 = " << radius2 << endl;
    
    
    if (dmin*dmin > radius2*2.) return false;

    pfem pf = mf.fem_of_element(cv);
    size_type nbdof = mf.nb_dof_of_element(cv);
    assert(pf->target_dim() == 1);
    mesher_level_set mls = ls->mls_of_convex(cv, lsnum);
    for (size_type i=0; i < nbdof; ++i) {
      base_node Pt = pf->node_convex(cv).points()[i];
     
      cout << "trying node " << i << " pt = " << Pt << " ptreal = " << mf.point_of_dof(dofs[i]) << " val = " << ls->values()[dofs[i]] << endl;
      try_projection(mls, Pt, true);
      cout << "  -> proj=" << Pt << " isin=" << pf->ref_convex(cv)->is_in(Pt) << "\n";
      if (pf->ref_convex(cv)->is_in(Pt) < EPS) return true;
    }
    return false;
  }

  void mesh_im_level_set::find_crossing_level_set(size_type cv, dal::bit_vector &prim, dal::bit_vector &sec) {
    prim.clear(); sec.clear();
    unsigned lsnum = 0;
    for (std::set<plevel_set>::const_iterator it = level_sets.begin(); 
	 it != level_sets.end(); ++it, ++lsnum) {
      if (is_crossed_by(cv, *it, 0)) {
	prim.add(lsnum);
	if ((*it)->has_secondary() && is_crossed_by(cv, *it, 1)) sec.add(lsnum);
      }
    }
  }


  
}  /* end of namespace getfem.                                             */



