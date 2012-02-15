// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2004-2011 Julien Pommier, Yves Renard
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

#include "getfem/getfem_mesher.h"

namespace getfem {

  void mesher_level_set::init_grad(void) const {
    gradient.resize(base.dim());
    for (dim_type d=0; d < base.dim(); ++d) {
      gradient[d] = base; gradient[d].derivative(d);
    }
    initialized = 1; 
  }

  void mesher_level_set::init_hess(void) const {
    if (initialized < 1) init_grad();
    hessian.resize(base.dim()*base.dim());
    for (dim_type d=0; d < base.dim(); ++d) {
      for (dim_type e=0; e < base.dim(); ++e) {
	hessian[d*base.dim()+e] = gradient[d];
	hessian[d*base.dim()+e].derivative(e);
      }
    }
    initialized = 2;
  }

  scalar_type mesher_level_set::grad(const base_node &P,
				     base_small_vector &G) const {
    if (initialized < 1) init_grad();
    gmm::resize(G, P.size());
    for (size_type i = 0; i < P.size(); ++i)
      G[i] = gradient[i].eval(P.begin());
    return (*this)(P);
  }

  void mesher_level_set::hess(const base_node &P, base_matrix &H) const {
    if (initialized < 2) init_hess();
    gmm::resize(H, P.size(), P.size()); 
    for (size_type i = 0; i < base.dim(); ++i)
      for (size_type j = 0; j < base.dim(); ++j) {
	H(i,j) = hessian[i*P.size()+j].eval(P.begin());
      }
  }


  //
  // Exported functions
  //

  bool try_projection(const mesher_signed_distance& dist, base_node &X,
		      bool on_surface) {
    base_small_vector G; base_node Y = X;
    scalar_type d = dist.grad(X, G), dmin = gmm::abs(d);
    size_type iter(0), count_falt(0);
    if (on_surface || d > 0.0)
      while (iter == 0 || dmin > 1e-15 || gmm::vect_dist2(X, Y) > 1e-15 ) {
	gmm::copy(X, Y);
	if (++iter > 1000) {
	  GMM_WARNING4("Try projection failed, 1000 iterations\n\n");
	  return false;
	  // return (gmm::abs(d) < 1E-10); // is there a possibility to detect
	} 	// the impossibility without making 1000 iterations ?
	gmm::scale(G, -d / std::max(1E-8, gmm::vect_norm2_sqr(G)));
	gmm::add(G, X);
//	scalar_type alpha(1);
//	scalar_type dn = dist(X);
// 	while (0 && gmm::abs(dn) > 2. * gmm::abs(d) && alpha > 0.1) {
// 	  alpha /= 2;
// 	  gmm::add(gmm::scaled(G, -alpha), X);
// 	  dn = dist(X);
// 	}
//	if (gmm::abs(alpha) < 1E-15) return (gmm::abs(d) < 1E-10);
	d = dist.grad(X, G);
	if (gmm::abs(d) >= dmin*0.95) ++count_falt;
	else { count_falt = 0; dmin = gmm::abs(d); }
	// if (count_falt > 20) return (gmm::abs(d) < 1E-10);
	if (count_falt > 20) return false;
      }
    return true;
  }
  
  // Try to find an intersection of a set of signed distance d_i.
  // Newton method on v solution to d_i(X+Gv) = 0, where the column of
  // G are the gradients of d_i. 
  bool pure_multi_constraint_projection
  (const std::vector<const mesher_signed_distance*> &list_constraints,
   base_node &X, const dal::bit_vector &cts) {
    size_type nbco = cts.card(), i(0), info, N = X.size();
    if (!nbco) return true;
    // cout << "nbco = " << nbco << endl;
    std::vector<const mesher_signed_distance*> ls(nbco);
    std::vector<scalar_type> d(nbco), v(nbco);
    std::vector<int> ipvt(nbco);
    gmm::col_matrix<base_node> G(N, nbco);
    gmm::dense_matrix<scalar_type> H(nbco, nbco);
    base_small_vector dd(N);
    for (dal::bv_visitor ic(cts); !ic.finished(); ++ic, ++i)
      { ls[i] = list_constraints[ic]; d[i] = -(ls[i]->grad(X, G[i])); }
    base_node oldX, aux(N);
    size_type iter = 0;
    scalar_type residual(0), alpha(0), det(0);

    if (nbco == 1) {
      // cout << "one constraint" << endl;
      try_projection(*(ls[0]), X, true);
      d[0] = -(ls[0]->grad(X, G[0]));
    } else {
      do {
	oldX = X;
	det = scalar_type(0); info = 1;
	alpha = -scalar_type(1);
	
	if (nbco <= N) {
	  if (nbco < N)
	    gmm::mult(gmm::transposed(G), G, H);
	  else
	    gmm::copy(gmm::transposed(G), H);
	  // cout << "H = " << H << endl;
	  info = gmm::lu_factor(H, ipvt);
	  det = scalar_type(1);
	  for (i = 0; i < nbco; ++i) det *= H(i,i);
	}
	// cout << "G = " << G << endl;
	
	if (info || gmm::abs(det) < 1E-40) {
	  dal::bit_vector cts_red = cts;
	  int eliminated = 0;
	  i = 0;
	  for (dal::bv_visitor ic(cts); !ic.finished(); ++ic, ++i) {
	    for (size_type j = 0; j < i; ++j)
	      gmm::add(gmm::scaled(G[j], -gmm::vect_sp(G[j], G[i])), G[i]);
	    scalar_type norm_gi = gmm::vect_norm2(G[i]);
	    // cout << "norm_gi = " << norm_gi << endl;
	    if (norm_gi < 1E-10)
	      { cts_red[ic] = false; eliminated++; }
	    else
	      gmm::scale(G[i], scalar_type(1)/norm_gi);
	  }
	  // cout << "G after = " << G << endl;
	  if (eliminated >= 1) {
	    // cout << "rec call with " << eliminated << " eliminated constraints" << endl;
	    pure_multi_constraint_projection(list_constraints, X, cts_red); 
	    for (i = 0; i < nbco; ++i) d[i] = -(ls[i]->grad(X, G[i]));
	    
	    det = scalar_type(0); info = 1;
	    if (nbco <= N) {
	      if (nbco < N)
		gmm::mult(gmm::transposed(G), G, H);
	      else
		gmm::copy(G, H);
	      info = gmm::lu_factor(H, ipvt);
	      det = scalar_type(1);
	      for (i = 0; i < nbco; ++i) det *= H(i,i);
	    }
	    
	  }
	}
	
	if (gmm::vect_norm2(d) > 1e-14) {
	  if (info || gmm::abs(det) < 1E-40) {
	    for (i = 0; i < nbco; ++i)
	      try_projection(*(ls[i]), X, true);
	    for (i = 0; i < nbco; ++i) d[i] = -(ls[i]->grad(X, G[i]));
	  }
	  else {
	    gmm::lu_solve(H, ipvt, v, d);
	    if (nbco < N)
	      gmm::mult(G, v, dd);
	    else
	      gmm::copy(v, dd);
	    GMM_ASSERT1(nbco <= N, "internal error");
	    gmm::add(dd, X);
	    for (i = 0; i < nbco; ++i) d[i] = -(ls[i]->grad(X, G[i]));
	    alpha = scalar_type(1);
	    if (iter > 0)
	      while (gmm::vect_norm2(d) > residual && alpha > 1E-10) {
		alpha /= scalar_type(2);
		gmm::add(gmm::scaled(dd, -alpha), X);
		for (i = 0; i < nbco; ++i) d[i] = -(ls[i]->grad(X, G[i]));
	      }
	    if (alpha < 1E-15) break;
	  }
	}
	
	++iter;
	residual = gmm::vect_norm2(d);
	// cout << "residual = " << residual << " alpha = " << alpha;
	// cout << " gmm::vect_dist2(oldX,X) = " << gmm::vect_dist2(oldX,X) << endl;
      } while (residual > 1e-14 && (residual > 1e-11 || iter < 15)
	       && iter < 200);
      // cout << "final residual = " << residual << endl;
    }
    for (i = 0; i < nbco; ++i) if (gmm::abs(d[i]) > SEPS) {
      //cout << "PURE MULTI HAS FAILED for " << cts << " nb iter = " << iter << endl;
      return false;
    }
    return true;
  }




  // return the eigenvalue of maximal abolute value for a
  // symmetric base_matrix
  static scalar_type max_vp(const base_matrix& M) {
    size_type m = mat_nrows(M);
    GMM_ASSERT1(is_hermitian(M), "Matrix is not symmetric");
    std::vector<scalar_type> eig(m);
    gmm::symmetric_qr_algorithm(M, eig);
    scalar_type emax(0);
    for (size_type i = 0; i < m; ++i) emax = std::max(emax, gmm::abs(eig[i]));
    return emax;
  }

  scalar_type curvature_radius_estimate(const mesher_signed_distance &dist,
					base_node X, bool proj) {  
    if (proj) try_projection(dist, X, true);
    base_small_vector V;
    base_matrix H;
    dist.grad(X, V);
    dist.hess(X, H);
    return gmm::vect_norm2(V) / std::max(1E-10, max_vp(H));
  }

  scalar_type min_curvature_radius_estimate
  (const std::vector<const mesher_signed_distance*> &list_constraints,
   const base_node &X, const dal::bit_vector &cts, size_type hide_first) {
    scalar_type r0 = 1E+10;
    for (dal::bv_visitor j(cts); !j.finished(); ++j) 
      if (j >= hide_first) {
      scalar_type r
	= curvature_radius_estimate(*(list_constraints[j]), X);
      r0 = std::min(r, r0);
    }
    return r0;
  }


  //
  // local functions
  //


  template <typename MAT, typename MAT2> void
  Frobenius_condition_number_sqr_gradient(const MAT& M, MAT2& G) { 
    typedef typename gmm::linalg_traits<MAT>::value_type T;
    typedef typename gmm::number_traits<T>::magnitude_type R;
    
    size_type n = mat_ncols(M);
    gmm::dense_matrix<T> B(n,n), C(n,n);
    gmm::mult(gmm::transposed(M), M, B);
    R trB = gmm::mat_trace(B);
    gmm::lu_inverse(B);
    R trBinv = gmm::mat_trace(B);
    gmm::mult(B,B,C);
    gmm::mult(gmm::scaled(M, T(-2)*trB), C, G);
    gmm::add(gmm::scaled(M, T(2)*trBinv), G);
  }


  struct pt_attribute {
    bool fixed;
    dal::bit_vector constraints;
    bool operator<(const pt_attribute &other) const {
      if (fixed && !other.fixed) return true;
      else if (!fixed && other.fixed) return false;
      else {
	if (constraints.last_true() > other.constraints.last_true())
	  return false;
	else if (constraints.last_true() < other.constraints.last_true())
	  return true;
	else if (constraints.card() > other.constraints.card()) return true;
	else if (constraints.card() < other.constraints.card()) return false;
	else for (dal::bv_visitor i1(constraints), i2(other.constraints);
		  !i1.finished(); ++i1, ++i2) {
	  if (i1 < i2) return true;
	  else if (i2 > i1) return false;
	}
      }
      return false;
    }
  };

  struct mesher {
    const mesher_signed_distance& dist;
    const mesher_virtual_function& edge_len;
    scalar_type h0, dist_point_hull, boundary_threshold_flatness;
    size_type N, K, iter_max, iter_wtcc;
    int prefind, noisy;
    base_node bounding_box_min, bounding_box_max;
    base_vector L, L0;

    std::vector<base_node> pts, pts_prev;
    std::vector<const pt_attribute*> pts_attr;
    std::set<pt_attribute> attributes_set;
    gmm::dense_matrix<size_type> t;    
    
    scalar_type ptol, ttol, L0mult, deltat, geps, deps;

    std::vector<const mesher_signed_distance*> constraints;

    gmm::dense_matrix<scalar_type> W;
    bgeot::mesh_structure edges_mesh;

    std::vector<size_type> attracted_points;
    std::vector<base_node> attractor_points;

    mesher(size_type K_,
	   const mesher_signed_distance& dist_, 
	   const mesher_virtual_function& edge_len_, 
	   scalar_type h0_,
	   mesh &m, const std::vector<base_node> &fixed_points,
	   int noise, size_type itm, int pref,
	   scalar_type dph, scalar_type btf)
      : dist(dist_), edge_len(edge_len_), dist_point_hull(dph),
	boundary_threshold_flatness(btf), iter_max(itm), prefind(pref),
	noisy(noise) {
      if (noise == -1) noisy = gmm::traces_level::level() - 2;
      K=K_; h0=h0_;
      ptol = 0.0025;
      ttol = .1;
      dist.bounding_box(bounding_box_min,bounding_box_max);
      N = bounding_box_min.size();
      if (N == 2) { 
	L0mult = 1.2; deltat = .2; geps = .001*h0; 
      } else {
	L0mult=1+0.4/pow(scalar_type(2),scalar_type(N-1));
	deltat=.1; geps=1e-1*h0;
      }
      deps=sqrt(1e-8)*h0;
      dist.register_constraints(this->constraints);

      bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(N,1);
      gmm::resize(W,N,N);
      base_matrix G(N,N+1); 
      vectors_to_base_matrix(G,
	    bgeot::equilateral_simplex_of_reference(dim_type(N))->points());
      gmm::mult(G, bgeot::geotrans_precomp(pgt,&pgt->convex_ref()
					   ->points(), 0)->grad(0), W);
      gmm::lu_inverse(W);
      do_build_mesh(m,fixed_points);
    }

    template<class TAB> scalar_type simplex_quality(const TAB &ppts) {
      base_matrix G(N,N), GW(N, N);
      for (size_type i=0; i < N; ++i) {
	base_node P = ppts[i+1] - ppts[0];
	std::copy(P.const_begin(), P.const_begin()+N, G.begin()+i*N);
      }
      gmm::mult(G, W, GW);
      return gmm::abs(1./gmm::condition_number(GW));
    }

    scalar_type worst_element, best_element;
    
    scalar_type fbcond_cost_function(const base_vector &c) {
      unsigned nbt = unsigned(gmm::mat_ncols(t));
      scalar_type cost = 0;
      base_matrix S(N,N), SW(N,N);
      worst_element = 1.; best_element = 1E40;
      for (unsigned i=0; i < nbt; ++i) {
	for (size_type j=0; j < N; ++j) {
	  for (size_type k=0; k < N; ++k) {
	    S(k,j) = c[t(j+1,i)*N+k] - c[t(0,i)*N+k];
	  }
	}
	gmm::mult(S,W,SW);
	if (gmm::lu_det(SW) < 1E-16) cost += 1e30;
	else {
	  scalar_type qual = gmm::Frobenius_condition_number_sqr(SW);
	  cost += qual;
	  worst_element = std::max(worst_element, qual / scalar_type(N*N));
	  best_element = std::min(best_element, qual / scalar_type(N*N));
	}
      }
      return cost / scalar_type(N * N);
    }

    void fbcond_cost_function_derivative(const base_vector& c,
					 base_vector &grad) {
      gmm::clear(grad);
      base_matrix Dcond(N,N), G(N,N), S(N,N), SW(N,N);
      unsigned nbt = unsigned(gmm::mat_ncols(t));
      
      for (unsigned i=0; i < nbt; ++i) {
	for (size_type j=0; j < N; ++j) {
	  for (size_type k=0; k < N; ++k) {
	    S(k,j) = c[t(j+1,i)*N+k] - c[t(0,i)*N+k];
	  }
	}
	gmm::mult(S,W,SW);
	Frobenius_condition_number_sqr_gradient(SW,Dcond);
	gmm::mult(Dcond, gmm::transposed(W), G);
	// gmm::mult(W, gmm::transposed(Dcond), G);
	for (size_type j=0; j < N; ++j) {
	  for (size_type k=0; k < N; ++k) {
	    grad[t(j+1,i)*N+k] += G(k,j);
	    grad[t(0,i)*N+k] -= G(k,j);
	  }
	}
      }
       for (unsigned i=0; i < pts.size(); ++i) {
   	if (pts_attr[i]->constraints.card() || pts_attr[i]->fixed) 
    	  for (size_type k=0; k < N; ++k) {
    	    grad[i*N+k] = 0;
    	  }
       }
      gmm::scale(grad, scalar_type(1) / scalar_type(N * N));
      
    }

    struct fbcond_cost_function_object {
      mesher &m;
      fbcond_cost_function_object(mesher &m_) : m(m_) {}
      scalar_type operator()(const base_vector& c)
      { return m.fbcond_cost_function(c); }
    };

    struct fbcond_cost_function_derivative_object {
      mesher &m;
      fbcond_cost_function_derivative_object(mesher &m_) : m(m_) {}
      void operator()(const base_vector& c, base_vector &grad)
      { m.fbcond_cost_function_derivative(c, grad); }
    };

    void optimize_quality() {

      size_type nbt = gmm::mat_ncols(t);

      if (noisy > 0) cout << "Quality post-optimization\n";
      base_vector X(pts.size() * N);
      for (unsigned i=0; i < pts.size(); ++i)
	gmm::copy_n(pts[i].const_begin(), N, X.begin() + i*N);

      base_matrix S(N,N), SW(N,N);
      for (unsigned i=0; i < nbt; ++i) {
	for (size_type j=0; j < N; ++j) {
	  for (size_type k=0; k < N; ++k) {
	    S(k,j) = X[t(j+1,i)*N+k] - X[t(0,i)*N+k];
	  }
	}
	if (gmm::lu_det(S) < 0) {
	  std::swap(t(0,i), t(1,i));
	  for (size_type j=0; j < N; ++j) {
	    for (size_type k=0; k < N; ++k) {
	      S(k,j) = X[t(j+1,i)*N+k] - X[t(0,i)*N+k];
	    }
	  }
	}
	if (noisy > 0 && gmm::abs(lu_det(S)) < 1e-10)
	  cout << "Element " << i << " is very bad, det = "
	       << gmm::abs(lu_det(S)) << "\n";
	gmm::mult(S,W,SW);
      }
      
      if (noisy > 0) {     
	cout << "Initial quality: "
	     << fbcond_cost_function(X)/scalar_type(nbt);
	cout << ", best element: " << sqrt(best_element)
	     << " worst element: " << sqrt(worst_element) << endl;      
      }  
      gmm::iteration iter; iter.set_noisy(noisy-1); iter.set_maxiter(1000);
      iter.set_resmax(5E-2*sqrt(scalar_type(nbt)));
      gmm::bfgs(fbcond_cost_function_object(*this), 
		fbcond_cost_function_derivative_object(*this),
		X, 10, iter, 0, 0.001, float(gmm::mat_ncols(t)));

      if (noisy > 0) cout << "Final quality: "
			  << fbcond_cost_function(X)/scalar_type(nbt)
			  << ", best element: " << sqrt(best_element)
			  << " worst element: " << sqrt(worst_element) << endl;

      for (unsigned i=0; i < pts.size(); ++i)
	gmm::copy_n(X.begin() + i*N, N, pts[i].begin());      
    }

    void delete_element(size_type ic) {
      if (ic != gmm::mat_ncols(t)-1) {
	for (size_type k=0; k < N+1; ++k)
	  std::swap(t(k,ic), t(k,gmm::mat_ncols(t)-1));
      }
      t.resize(N+1,gmm::mat_ncols(t)-1);
    }

    scalar_type quality_of_element(size_type i) {
      return simplex_quality(gmm::index_ref_iterator
			     (pts.begin(), gmm::mat_col(t,i).begin()));
    }

    void control_mesh_surface(void) {
	mesh m;
	adapt_mesh(m, 1);
	dal::bit_vector ii = m.convex_index(), points_to_project;
	size_type ic, ipt;	
	for (ic << ii; ic != size_type(-1); ic << ii) {
	  for (short_type f = 0; f <= N; ++f) {
	    if (!m.is_convex_having_neighbour(ic,f)) {
	      for (unsigned i = 0; i < N; ++i) {
		ipt = m.ind_points_of_face_of_convex(ic, f)[i];
		if (pts_attr[ipt]->constraints.card() == 0)
		  points_to_project.add(ipt);
		else if (dist(pts[ipt]) < -1e-2) 
		  cout << "WARNING, point " << ipt 
		       << " incoherent !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
	      }
	    }
	  }
	}
	if (points_to_project.card()) {
	  iter_wtcc = 0;
	  if (noisy > 1)
	    cout << "points to project : " << points_to_project << endl;
	  ii = points_to_project;
	  for (ipt << ii; ipt != size_type(-1); ipt << ii)
	    surface_projection_and_update_constraints(ipt);
	}
    }

    void suppress_flat_boundary_elements(void) {
      size_type nb_deleted = 0;
      do {
	mesh m;
	adapt_mesh(m, 1);
	std::vector<base_small_vector> normals;
	dal::bit_vector ii = m.convex_index();
	dal::bit_vector flat_bound_elemnts;
	size_type ic;
	
	for (ic << ii; ic != size_type(-1); ic << ii) {
	  scalar_type max_flatness = -2.0;
	  normals.resize(0);
	  for (short_type f = 0; f <= N; ++f) {
	    if (!m.is_convex_having_neighbour(ic,f)) {
	      if (quality_of_element(ic) < 1E-8) max_flatness = 1E-8;
	      else {
		base_small_vector n = m.normal_of_face_of_convex(ic, f);
		normals.push_back(n / gmm::vect_norm2(n));
	      }
	    }
	  }
	  
	  if (noisy > 1 && max_flatness != -2.0)
	    cout << "flatness of element " << ic << " : " << max_flatness;
	  if (normals.size() >= 2) {
	    if (noisy > 1) cout << "flatness of element " << ic << " : ";
	    
	    for (unsigned i = 1; i < normals.size(); ++i)
	      for (unsigned j = 0; j < i; ++j) {
		scalar_type flatness=1.0-gmm::vect_sp(normals[i], normals[j]);
		max_flatness = std::max(max_flatness, flatness);
		if (noisy > 1) cout << flatness << " ";
	      }
	  }
	  if (max_flatness < boundary_threshold_flatness
	      && max_flatness!=-2.0) {
	    flat_bound_elemnts.add(ic); 
	    if (noisy > 1) cout << " -> deleting";
	  }
	  if ((normals.size() >= 2 || max_flatness != -2.0) && noisy > 1)
	    cout << endl;
	  
	}
	ii = flat_bound_elemnts; nb_deleted = flat_bound_elemnts.card();
	for (ic = ii.take_last(); ic != size_type(-1); ic = ii.take_last())
	  delete_element(ic);
      } while (nb_deleted > 0);
    }
    

    bool try_projection(base_node &X)
    { return getfem::try_projection(dist, X); }

    void projection(base_node &X) {
      base_small_vector G(X.size());
      scalar_type d = dist.grad(X, G);
      size_type it(0);
      if (d > 0.0)
	while (gmm::abs(d) > 1e-10) {
	  ++it;
	  GMM_ASSERT1(it <= 10000, "Object empty, or bad signed distance");
// 	  cout << "iter " << it << " X = " << X << " dist = " << d << 
// 	    " grad = " << G << endl;
	  gmm::add(gmm::scaled(G, -d / gmm::vect_norm2_sqr(G)), X);
	  d = dist.grad(X, G);
	}
    }
    
    void surface_projection(base_node &X) { 
      base_small_vector G;
      scalar_type d = dist.grad(X, G);
      size_type it(0);
      while (gmm::abs(d) > 1e-10) {
	++it;
	GMM_ASSERT1(it <= 10000, "Object empty, or bad signed distance X="
		    << X << ", G=" << G << " d = " << d);
//  	if (it > 9980) 
// 	  cout << "iter " << it << " X = " << X << " dist = " << d << 
//  	    " grad = " << G << endl;
	gmm::add(gmm::scaled(G, -d / gmm::vect_norm2_sqr(G)), X);
	d = dist.grad(X, G);
      }
    }

    void projection(base_node &X, dal::bit_vector& bv)
    { projection(X); bv.clear(); dist(X,bv); }

    void constraint_projection(base_node &X, size_type cnum) {
      base_small_vector G;
      scalar_type d = constraints[cnum]->grad(X, G);
      while (gmm::abs(d) > 1e-10) {
	gmm::add(gmm::scaled(G, -d / gmm::vect_norm2_sqr(G)), X);
	d=constraints[cnum]->grad(X, G);
      }
    }

    bool pure_multi_constraint_projection(base_node &X,
					  const dal::bit_vector &cts) {
      getfem::pure_multi_constraint_projection(constraints, X, cts);

//       base_node oldX;
//       size_type cnt = 0;
//       do {     
// 	oldX = X;
// 	for (dal::bv_visitor ic(cts); !ic.finished(); ++ic)
// 	  constraint_projection(X,ic);
// 	++cnt;
//       } while (cts.card() && gmm::vect_dist2(oldX,X) > 1e-14 && cnt < 1000);
//       if (cnt >= 1000) return false;
//       if (cts.card()) {
      dal::bit_vector ct2;
      dist(X, ct2);
      return ct2.contains(cts);
//       }
//       return true;
    }

    bool multi_constraint_projection(base_node &X,
				     const dal::bit_vector &cts) {
      if (!(cts.card())) { projection(X); return true; }
      else {
	base_node oldX;
	size_type cnt = 0;
	do {     
	  oldX = X;
	  for (dal::bv_visitor ic(cts); !ic.finished(); ++ic)
	    constraint_projection(X,ic);
	  projection(X);
	  ++cnt;
	} while (gmm::vect_dist2(oldX,X) > 1e-14 && cnt < 1000);
	if (cnt >= 1000) return false;
	dal::bit_vector ct2;
	dist(X, ct2);
	return ct2.contains(cts);
      }
    }

    void tangential_displacement(base_small_vector &X,
				 const dal::bit_vector &cts) {
      base_small_vector G;
      base_matrix normals(N, cts.card());
      size_type cnt = 0;
      for (dal::bv_visitor ic(cts); !ic.finished(); ++ic) {
	constraints[ic]->grad(X, G);
	gmm::copy(G, gmm::mat_col(normals, cnt));
	for (size_type k=0; k < cnt; ++k) {
	  gmm::add(gmm::scaled(gmm::mat_col(normals, k),
			       -gmm::vect_sp(gmm::mat_col(normals,k),
					     gmm::mat_col(normals,cnt))),
		   gmm::mat_col(normals,cnt));
	  scalar_type n = gmm::vect_norm2(gmm::mat_col(normals, cnt));
	  if (n < 1e-8) continue;
	  gmm::scale(gmm::mat_col(normals, cnt), 1./n);
	  gmm::add(gmm::scaled(gmm::mat_col(normals, cnt),
			       -gmm::vect_sp(X,gmm::mat_col(normals, cnt))),X);
	  ++cnt;
	}
      }
    }

    const pt_attribute* get_attr(bool fixed, const dal::bit_vector &bv) {
      pt_attribute a; 
      a.fixed = fixed;
      a.constraints = bv; 
      return &(*attributes_set.insert(a).first);
    }

    void surface_projection_and_update_constraints(size_type ip) { 
      surface_projection(pts[ip]);
      dal::bit_vector new_cts;
      dist(pts[ip], new_cts);
      pts_attr[ip] = get_attr(pts_attr[ip]->fixed, new_cts);
    }

    void project_and_update_constraints(size_type ip) {
      const dal::bit_vector& cts = pts_attr[ip]->constraints;
      dal::bit_vector new_cts;
      multi_constraint_projection(pts[ip], cts);
      dist(pts[ip], new_cts);
      if (noisy > 1 && !new_cts.contains(cts)) {
	cout << "Point #" << ip << " has been downgraded from "
	     << cts << " to " << new_cts << endl;
      }
      else if (noisy > 1 && new_cts.card() > cts.card()) {
	cout << "Point #" << ip << " has been upgraded from " 
	     << cts << " to " << new_cts << endl;
      }
      if (new_cts != cts) {
	pts_attr[ip] = get_attr(pts_attr[ip]->fixed, new_cts);
	iter_wtcc = 0;
      }
      
    }
    
    template <class VECT> void move_carefully(size_type ip, const VECT &VV) {
      base_node V(N); gmm::copy(VV, V);
//       if (pts_attr[ip]->constraints.card() != 0) {
// 	base_small_vector grad;
//      dist.grad(pts[ip], grad);
// 	scalar_type r = gmm::vect_norm2_sqr(grad);
// 	gmm::add(gmm::scaled(grad, -gmm::vect_sp(V, grad) / r), V);
//       }
      scalar_type norm = gmm::vect_norm2(V);
      if (norm <= h0/scalar_type(4))
	gmm::add(V, pts[ip]);
      else
	gmm::add(gmm::scaled(V, h0 / (scalar_type(4) * norm)), pts[ip]);
      project_and_update_constraints(ip);
    }

     template <class VECT> void move_carefully(const VECT &V) {
       scalar_type norm_max(0), lambda(1);
       size_type npt = gmm::vect_size(V) / N;
       for (size_type i = 0; i < npt; ++i)
	 norm_max = std::max(norm_max, gmm::vect_norm2
			     (gmm::sub_vector(V, gmm::sub_interval(i*N, N))));
       if (norm_max > h0/scalar_type(3.7))
       	 lambda = h0 / (scalar_type(3.7) * norm_max);
       
       for (size_type i = 0; i < npt; ++i)
       	 move_carefully(i, gmm::scaled(gmm::sub_vector
				       (V, gmm::sub_interval(i*N, N)),lambda));
     }

    void distribute_points_regularly(const std::vector<base_node>
				     &fixed_points) {
      size_type nbpt = 1;
      std::vector<size_type> gridnx(N);
      mesh m;
      base_node eff_boxmin(N), eff_boxmax(N);
      bool eff_box_init = false;

      for (size_type i=0; i < N; ++i) 
	h0 = std::min(h0, bounding_box_max[i] - bounding_box_min[i]);
      GMM_ASSERT1(h0 >= 1E-10, "h0 = " << h0 << " too small, aborting.");

      for (size_type i=0; i < N; ++i) {
	scalar_type h = h0;
	if (N == 2 && i == 1) h = sqrt(3.)/2. * h0;
	gridnx[i]=1+(size_type)((bounding_box_max[i]-bounding_box_min[i])/h);
	nbpt *= gridnx[i];
      }

      /* build the regular grid and filter points outside */
      for (size_type i=0; i < fixed_points.size(); ++i) {
	if (dist(fixed_points[i]) < geps &&
	    m.search_point(fixed_points[i]) == size_type(-1)) {
	  m.add_point(fixed_points[i]); 
	  pts.push_back(fixed_points[i]); 
	  pts_attr.push_back(get_attr(true,dal::bit_vector()));
	}
	else if (noisy > 0)
	  cout << "Removed duplicate fixed point: "<<fixed_points[i]<<"\n";
      }
      base_node P(N), Q(N);
      for (size_type i=0; i < nbpt; ++i) {
	for (size_type k=0, r = i; k < N; ++k) {
	  unsigned p =  unsigned(r % gridnx[k]);
	  P[k] = p * (bounding_box_max[k] - bounding_box_min[k]) / 
	    scalar_type((gridnx[k]-1)) + bounding_box_min[k];
	  if (N==2 && k==0 && ((r/gridnx[0])&1)==1) P[k] += h0/2;
	  r /= gridnx[k];
	}

	dal::bit_vector co;
	if ((prefind == 1 && dist(P) < 0) || prefind == 2) {
	  for (size_type k = 0; k < constraints.size() && co.card() < N; ++k) {
	    gmm::copy(P, Q);
	    if (gmm::abs((*(constraints[k]))(Q)) < h0) {
	      constraint_projection(Q, k);
	      if (dist(Q) > -geps && 
		  (gmm::vect_dist2(P, Q) < h0 / scalar_type(2))) co.add(k);
	    }
	  }
	}
	gmm::copy(P, Q);
	if (co.card() > 0) { 
	  bool ok = pure_multi_constraint_projection(Q, co);
	  if (!ok || gmm::abs(dist(Q)) > geps) { gmm::copy(P, Q); co.clear(); }
	}

	if (prefind == 3) {
	  try_projection(Q);
	}

	if (dist(Q) < geps) {
	  if (m.search_point(Q) == size_type(-1)) {
	    //cout << "adding point : " << Q << endl;
	    if (!eff_box_init)
	      { eff_boxmin = eff_boxmax = Q; eff_box_init = true; }
	    else for (size_type k = 0; k < N; ++k) {
	      eff_boxmin[k] = std::min(eff_boxmin[k], Q[k]);
	      eff_boxmax[k] = std::max(eff_boxmax[k], Q[k]);
	    }
	    m.add_point(Q); pts.push_back(Q);
	    pts_attr.push_back(get_attr(false, co));
	  }
	}
      }
      if (noisy > 0) cout << "effective bounding box : " << eff_boxmin << " : " << eff_boxmax << endl;
      if (prefind == 3) {
	h0 = std::min(h0, gmm::vect_dist2(eff_boxmin, eff_boxmax)
		      / scalar_type(2));
      }
    }

    void add_point_hull(void) { 
      if (dist_point_hull > 0) {
	size_type nbpt = pts.size(), nbadd(0);
	base_node P, Q;
	base_small_vector V;
	for (unsigned i=0; i < nbpt; ++i) {
	  if (pts_attr[i]->constraints.card()) {
	    P = pts[i];
	    dist.grad(P, V);
	    scalar_type d = gmm::vect_norm2(V);
	    if (d > 0) {
	      P += V * (dist_point_hull*h0/d);
	      if (dist(P)*sqrt(scalar_type(N)) > dist_point_hull*h0) {
		Q = P;
		projection(Q);
		if (gmm::vect_dist2(P, Q) > dist_point_hull*h0/scalar_type(2))
		  { pts.push_back(P); ++nbadd; }
	      }
	    }
	  }
	}
	if (noisy > 1) cout << "point hull: " << nbadd << " points added\n";
      }
    }


    scalar_type pts_dist_max(const std::vector<base_node> &A, 
		      const std::vector<base_node> &B) {
      scalar_type dist_max = 0;
      for (size_type i=0; i < pts.size(); ++i) 
	dist_max = std::max(dist_max, gmm::vect_dist2_sqr(A[i],B[i]));
      return sqrt(dist_max);
    }

    struct cleanup_points_compare {
      const std::vector<base_node> &pts;
      const std::vector<const pt_attribute*> &attr;
      cleanup_points_compare(const std::vector<base_node> &pts_, 
			     const std::vector<const pt_attribute*> &attr_)
	: pts(pts_), attr(attr_) {}
      bool operator()(size_type a, size_type b) {
	if (attr[a] != attr[b]) return attr[a] < attr[b];
	return pts[a] < pts[b];
      }
    };
    void cleanup_points() {
      std::vector<size_type> idx(pts.size());
      for (size_type i=0; i < idx.size(); ++i) idx[i] = i;
   
      std::sort(idx.begin(), idx.end(), cleanup_points_compare(pts,pts_attr));
      bgeot::kdtree tree;
      bgeot::kdtree_tab_type neighbours;
      dal::bit_vector keep_pts; keep_pts.add(0,idx.size());
      for (size_type i=0, i0=0; i < idx.size(); ++i) {
	const base_node &P = pts[idx[i]];
	const pt_attribute *a = pts_attr[idx[i]];
	tree.add_point_with_id(P,i);
	if (i == idx.size()-1 || a != pts_attr[idx[i+1]]) {
	  for (size_type j=i0; j < i+1; ++j) {
	    base_node bmin = P, bmax = P;
	    scalar_type h = h0*edge_len(P);
	    for (size_type k = 0; k < N; ++k)
	      { bmin[k] -= h/20.; bmax[k] += h/20.; }
	    
	    tree.points_in_box(neighbours, bmin, bmax);
	    for (size_type k=0; k < neighbours.size(); ++k) {
	      if (neighbours[k].i != i && keep_pts.is_in(neighbours[k].i)
		  && keep_pts.is_in(i)) {
		if (noisy > 0)
		  cout << "point #" << i << " " << P
		       << " is too near from point #"  << neighbours[k].i
		       << pts[idx[neighbours[k].i]] << " : will be removed\n";
		keep_pts.sup(i);
	      }
	    }
	  }
	  tree.clear();
	  i0 = i+1;
	}
      }
      pts_prev.resize(keep_pts.card());
      size_type cnt = 0;
      std::vector<const pt_attribute*> pts_attr2(keep_pts.card());
      for (dal::bv_visitor i(keep_pts); !i.finished(); ++i, ++cnt) {
	pts_prev[cnt].swap(pts[idx[i]]);
	pts_attr2[cnt] = pts_attr[idx[i]];
      }
      pts_attr.swap(pts_attr2);
      pts.resize(pts_prev.size()); 
      std::copy(pts_prev.begin(), pts_prev.end(), pts.begin());
    }

    scalar_type worst_q;
    base_node worst_q_P;

    void select_elements(int version) {
      size_type nbpt = pts.size();

      worst_q = 1.;
      base_node weights(N+1);
      for (size_type i=0; i < gmm::mat_ncols(t); )  {
	bool ext_simplex = false;
	bool boundary_simplex = true;
	bool on_boundary_simplex = false;
	bool is_bridge_simplex = false;
	scalar_type q(0), dG(0);
	base_node G;
	
	for (size_type k=0; k <= N; ++k)
	  if (t(k, i) >= nbpt) ext_simplex = true;

	if (!ext_simplex) {
	  G = pts[t(0,i)];
	  for (size_type k=1; k <= N; ++k) G += pts[t(k,i)];
	  gmm::scale(G, scalar_type(1)/scalar_type(N+1));
	  dG = dist(G);
	  gmm::clear(weights);
	  
	  q = quality_of_element(i);
	  
	  for (size_type k=0; k <= N; ++k) {
	    if (pts_attr[t(k,i)]->constraints.card() == 0)
	      boundary_simplex = false;
	    else
	      on_boundary_simplex = true;
	  }
	  
	  if (version == 1 && on_boundary_simplex) 
	    for (size_type k=1; k < N+1; ++k) 
	      for (size_type l=0; l < k; ++l) {
		dal::bit_vector all_cts = pts_attr[t(k,i)]->constraints
		  | pts_attr[t(l,i)]->constraints;
		if (/* gmm::vect_dist2(pts[t(k,i)], pts[t(l,i)]) > h0 && */
		    !(pts_attr[t(k,i)]->constraints.contains(all_cts))
		    && !(pts_attr[t(l,i)]->constraints.contains(all_cts))
		    && dist(0.5*(pts[t(k,i)] + pts[t(l,i)])) > 0.)
		  is_bridge_simplex = true;
	      }
	}
	if (ext_simplex || dG > 0 || is_bridge_simplex || q < 1e-14) {
	  delete_element(i);
	} else {
	  ++i;
	  if (q < worst_q) {
	    worst_q = q;
	    worst_q_P = G*(scalar_type(1)/scalar_type(N+1));
	  }
	}
      }
      
    }

    void adapt_mesh(mesh &m, size_type degree) {
      std::vector<base_node> cvpts(N+1), cvpts2;
      size_type cvnum;
      m.clear();
      for (size_type ip=0; ip < pts.size(); ++ip) {
	size_type z;
	bgeot::small_vector<scalar_type> P = pts[ip];
	while ((z = m.add_point(P)) != ip) {
	  if (noisy > 0) cout << "WARNING : points are too near ...\n";
	  bgeot::small_vector<scalar_type> Z(N); gmm::fill_random(Z);
	  gmm::add(gmm::scaled(Z, h0/1000.0), P);
	}
      }
      for (size_type i=0; i < t.size()/(N+1); ++i) {
	for (size_type k=0; k < N+1; ++k) cvpts[k] = pts[t[i*(N+1)+k]];
	if (degree == 1) {
	  cvnum = m.add_convex(bgeot::simplex_geotrans(N,1), &t[i*(N+1)]);
	  assert(cvnum == i);
	} else {
	  bgeot::pgeometric_trans pgt =
	    bgeot::simplex_geotrans(N,short_type(degree));
	  cvpts2.resize(pgt->nb_points());
	  for (size_type k=0; k < pgt->nb_points(); ++k) {
	    cvpts2[k] = bgeot::simplex_geotrans(N,1)->transform
	      (pgt->convex_ref()->points()[k], cvpts);
	  }
	  cvnum = m.add_convex_by_points(pgt, cvpts2.begin());
	  assert(cvnum == i);
	}
      }
      if (degree>1) {
	//m.optimize_structure();
	getfem::mesh_region border_faces;
	getfem::outer_faces_of_mesh(m, border_faces);
	dal::bit_vector ptdone; ptdone.sup(0,m.points_index().last_true());
	for (getfem::mr_visitor it(border_faces); !it.finished(); ++it) {
	  mesh::ind_pt_face_ct fpts_
	    = m.ind_points_of_face_of_convex(it.cv(), it.f());
	  std::vector<size_type> fpts(fpts_.size());
	  std::copy(fpts_.begin(), fpts_.end(), fpts.begin());
	  interpolate_face(m, ptdone, fpts, 
			   m.trans_of_convex(it.cv())->structure()
			   ->faces_structure()[it.f()]);
	}
      }
    }



    void interpolate_face(mesh &m, dal::bit_vector& ptdone, 
			  const std::vector<size_type>& ipts,
			  bgeot::pconvex_structure cvs) {
      if (cvs->dim() == 0) return;
      else if (cvs->dim() > 1) {
	std::vector<size_type> fpts;
	for (short_type f=0; f < cvs->nb_faces(); ++f) {
	  fpts.resize(cvs->nb_points_of_face(f));
	  for (size_type k=0; k < fpts.size(); ++k)
	    fpts[k] = ipts[cvs->ind_points_of_face(f)[k]];
	  interpolate_face(m,ptdone,fpts,cvs->faces_structure()[f]);
	}
      }
      dal::bit_vector cts; size_type cnt = 0;
      for (size_type i=0; i < ipts.size(); ++i) {
	if (ipts[i] < pts.size()) { 
	  if (cnt == 0) cts = pts_attr[ipts[i]]->constraints;
	  else cts &= pts_attr[ipts[i]]->constraints;
	  ++cnt;
	}
      }
      if (cts.card()) {
	// dal::bit_vector new_cts;
	for (size_type i=0; i < ipts.size(); ++i) {
	  if (ipts[i] >= pts.size() && !ptdone[ipts[i]]) { 
	    base_node &P = m.points()[ipts[i]];
	    multi_constraint_projection(P, cts);
	    // dist(P, new_cts);
	  }
	}
      }
    }

    void special_constraints_management(void) {

      bgeot::kdtree tree;
      bgeot::kdtree_tab_type neighbours;
      bool tree_empty = true;
      mesh::ind_set iAneighbours, iBneighbours, common_pts;
	  
      attractor_points.resize(0); attracted_points.resize(0);
      
      for (dal::bv_visitor ie(edges_mesh.convex_index());
	   !ie.finished(); ++ie) {
	size_type iA = edges_mesh.ind_points_of_convex(ie)[0];
	size_type iB = edges_mesh.ind_points_of_convex(ie)[1];
	// if (L[ie] > L0[ie]) continue;
	if (pts_attr[iA] == pts_attr[iB] ||
	    pts_attr[iA]->constraints.card() == 0 ||
	    pts_attr[iB]->constraints.card() == 0) continue;
	if (pts_attr[iA]->constraints == pts_attr[iB]->constraints)
	  continue;
	dal::bit_vector bv1(pts_attr[iA]->constraints);
	bv1.setminus(pts_attr[iB]->constraints);
	dal::bit_vector bv2(pts_attr[iB]->constraints);
	bv2.setminus(pts_attr[iA]->constraints);
	if (bv1.card() && bv2.card()) {
	  bv1 |= bv2;
	  edges_mesh.ind_points_to_point(iA, iAneighbours);
	  edges_mesh.ind_points_to_point(iB, iBneighbours);
	  common_pts.resize(0);
	  for (size_type i = 0; i < iAneighbours.size(); ++i)
	    if (std::find(iBneighbours.begin(), iBneighbours.end(),
			  iAneighbours[i]) != iBneighbours.end())
	      common_pts.push_back(iAneighbours[i]);
	  bool do_projection = true;
	  if (dist(.5*(pts[iA]+pts[iB])) < 0) {
	    for (mesh::ind_set::iterator it = common_pts.begin();
		 it != common_pts.end(); ++it) {
	      if (pts_attr[*it]->constraints.contains(bv1)) {
		do_projection = false;
		break;
	      }
	    }
	  }
	  if (do_projection) {
	    
	    if (pts_attr[iA]->constraints.card()
		< pts_attr[iB]->constraints.card()) std::swap(iA,iB);
	    
	    base_node PA = pts[iA];
	    bool okA = multi_constraint_projection(PA, bv1);
	    base_node PB = pts[iB];
	    bool okB = multi_constraint_projection(PB, bv1);
	    
	    if (okB && !okA)
	      { std::swap(PA,PB); std::swap(iA,iB); std::swap(okB, okA); }
	    
	    if (okB && (gmm::vect_dist2(PA,pts[iA])
			> 1.1*gmm::vect_dist2(PB,pts[iB]))) {
	      // 1.5 au lieu de 1.1 ?
	      std::swap(iA,iB); std::swap(PA,PB);
	    }
	    
	    
	    if (okA && gmm::vect_dist2(PA, pts[iA]) < h0*0.75) {
	      
	      if (tree_empty) {
		for (size_type i=0; i < pts.size(); ++i)
		  tree.add_point_with_id(pts[i],i);
		tree_empty = false;
	      }
	      
	      base_node bmin = PA, bmax = PA;
	      for (size_type k = 0; k < N; ++k)
		{ bmin[k] -= h0/1.8; bmax[k] += h0/1.8; }
	      tree.points_in_box(neighbours, bmin, bmax);
	      for (size_type k=0; k < neighbours.size(); ++k) {
		if (neighbours[k].i != iA && neighbours[k].i != iB)
		  do_projection = false;
	      }
	      
	      if (do_projection) {
		attractor_points.push_back(PA);
		attracted_points.push_back(iA);
	      }
	    }
	  }
	}
      }
    }


    void running_delaunay(bool mct) {
      if (noisy > 0)
	cout << "NEW DELAUNAY, running on " << pts.size() << " points\n";
      size_type nbpt = pts.size();
      add_point_hull();
      delaunay(pts, t);
      pts.resize(nbpt);
      if (noisy > 1) cout << "number of elements before selection = "
			  << gmm::mat_ncols(t) << "\n";
      if (mct) {
	select_elements(0);
	edges_mesh.clear();
	for (size_type i=0; i < gmm::mat_ncols(t); ++i)
	  for (size_type j=0; j < N+1; ++j)
	    for (size_type k=j+1; k < N+1; ++k)
	      edges_mesh.add_segment(t(j,i), t(k,i));
	special_constraints_management();
      }
      select_elements(1);
      if (noisy > 0) cout << "number of elements after selection = "
			  << gmm::mat_ncols(t) << "\n";
      edges_mesh.clear();
      for (size_type i=0; i < gmm::mat_ncols(t); ++i)
	for (size_type j=0; j < N+1; ++j)
	  for (size_type k=j+1; k < N+1; ++k)
	    edges_mesh.add_segment(t(j,i), t(k,i));
    }

    void standard_move_strategy(base_vector &X) {
      for (dal::bv_visitor ie(edges_mesh.convex_index());
	   !ie.finished(); ++ie) {
	size_type iA = edges_mesh.ind_points_of_convex(ie)[0];
	size_type iB = edges_mesh.ind_points_of_convex(ie)[1];
	base_node bar = pts[iB]-pts[iA];
	scalar_type F = std::max(L0[ie]-L[ie], 0.);
	
	if (F) {
	  base_node Fbar = (bar)*(F/L[ie]);
	  
	  if (!pts_attr[iA]->fixed) // pts[iA] -= deltat*Fbar; 
	    gmm::add(gmm::scaled(Fbar, -deltat),
		     gmm::sub_vector(X, gmm::sub_interval(iA*N, N)));
	  if (!pts_attr[iB]->fixed) // pts[iB] += deltat*Fbar;
	    gmm::add(gmm::scaled(Fbar, deltat),
		     gmm::sub_vector(X, gmm::sub_interval(iB*N, N)));
	}
      }
    }

    void do_build_mesh(mesh &m,
		       const std::vector<base_node> &fixed_points) {

      distribute_points_regularly(fixed_points);
      std::vector<base_node> pts2(pts.size(),base_node(N));
      size_type count = 0, count_id = 0, count_ct = 0;
      bool pt_changed = false;
      iter_wtcc = 0;
      attracted_points.resize(0);
      attractor_points.resize(0);

      do {
	if (noisy > 1) {
	  cout << "Iter " << count << " / " << count + iter_max - iter_wtcc;
	  if (count && pts_prev.size() == pts.size())
	    cout << ", dist_max since last delaunay ="
		 << pts_dist_max(pts, pts_prev) << ", tol=" << ttol*h0;
	  cout << endl;
	}
	if (count==0 || pts_prev.size() != pts.size()
	    || (pts_dist_max(pts, pts_prev) > ttol*h0 && count_id >= 5)) {
	  size_type nbpt = pts.size();
	  cleanup_points(); /* and copy pts to pts_prev */
	  if (noisy == 1) cout << "Iter " << count << " ";
	  bool mct = false;
	  if (count_ct >= 20 || pt_changed) {
	    mct = true;
	    count_ct = 0;
	  }
	  running_delaunay(mct);
	  pt_changed = nbpt != pts.size();
	  count_id = 0;
	}
	++count_id; ++count_ct;
	pts2 = pts;

	// computation of L and L0.
	size_type nbcv = edges_mesh.convex_index().card();
	GMM_ASSERT1(nbcv != 0, "no more edges!");
	L.resize(nbcv); L0.resize(nbcv);
	scalar_type sL = 0, sL0 = 0;
	for (dal::bv_visitor ie(edges_mesh.convex_index());
	     !ie.finished(); ++ie) {
	  const base_node &A = pts[edges_mesh.ind_points_of_convex(ie)[0]];
	  const base_node &B = pts[edges_mesh.ind_points_of_convex(ie)[1]];
	  base_node C(A); C+=B; C /= scalar_type(2);
	  L[ie] = gmm::vect_dist2(A, B);
	  L0[ie] = edge_len(C);
	  sL += pow(L[ie],scalar_type(N));
	  sL0 += pow(L0[ie],scalar_type(N));
	}
	gmm::scale(L0, L0mult * pow(sL/sL0, scalar_type(1)/scalar_type(N)));

	// Moving the points with standard strategy
	base_vector X(pts.size() * N);
	standard_move_strategy(X);
	for (size_type i = 0; i < attracted_points.size(); ++i) {
	  size_type npt = attracted_points[i];
	  gmm::copy(attractor_points[i] - pts[npt],
		    gmm::sub_vector(X, gmm::sub_interval(npt*N, N)));
	}

	move_carefully(X);

	scalar_type maxdp = pts_dist_max(pts,pts2);
	if (noisy > 1) 
	  cout << ", maxdp = " << maxdp << ", ptol = "
	       << ptol << " CV=" << sqrt(maxdp)*deltat/h0 << "\n";
	++count; ++iter_wtcc;

	if (iter_wtcc == 100) control_mesh_surface();

	// m.clear();
	// for (size_type i=0; i < t.size()/(N+1); ++i)
	//  m.add_convex_by_points(bgeot::simplex_geotrans(N,1),
	//		 dal::index_ref_iterator(pts.begin(), &t[i*(N+1)]));
	// char s[50]; sprintf(s, "toto%02d.mesh", count);
	// m.write_to_file(s);

	if ( (count > 40 && sqrt(maxdp)*deltat < ptol * h0)
	     || iter_wtcc>iter_max || count > 10000) {

// 	  {
// 	    m.clear();
// 	    adapt_mesh(m,K);
// 	    m.optimize_structure();
// 	    getfem::vtk_export exp("toto1.vtk");
// 	    exp.exporting(m);
// 	    exp.write_mesh_quality(m);
// 	  }

	  control_mesh_surface();
	  size_type nbpt = pts.size();
	  add_point_hull();	  
	  delaunay(pts, t);
	  pts.resize(nbpt);
	  select_elements((prefind == 3) ? 0 : 1);
	  suppress_flat_boundary_elements();

// 	  {
// 	    m.clear();
// 	    adapt_mesh(m,K);
// 	    m.optimize_structure();
// 	    getfem::vtk_export exp("toto2.vtk");
// 	    exp.exporting(m);
// 	    exp.write_mesh_quality(m);
// 	  }

	  if (prefind != 3) optimize_quality();
	  
	  // ajout d'un point au barycentre des elements trop plats : 
	  
	  if (prefind != 3)
 	    for (unsigned cv = 0; cv < gmm::mat_ncols(t); ++cv) {
	      
	      if (quality_of_element(cv) < 0.05) {
		base_node G = pts[t(0,cv)];
		for (size_type k=1; k <= N; ++k) G += pts[t(k,cv)];
		gmm::scale(G, scalar_type(1)/scalar_type(N+1));
		pts.push_back(G);
		pts_attr.push_back(get_attr(false, dal::bit_vector()));
	      }
	    }
	  
	  if (pts.size() != nbpt) {
	    control_mesh_surface();
	    nbpt = pts.size();
	    add_point_hull();
	    delaunay(pts, t);
	    pts.resize(nbpt);
	    select_elements((prefind == 3) ? 0 : 1);
	    suppress_flat_boundary_elements();

// 	    {
// 	      m.clear();
// 	      adapt_mesh(m,K);
// 	      m.optimize_structure();
// 	      getfem::vtk_export exp("toto3.vtk");
// 	      exp.exporting(m);
// 	      exp.write_mesh_quality(m);
// 	    }
	    
	    if (prefind != 3) optimize_quality();
	  }
	  break;
	}


      } while (true);

      m.clear();
      adapt_mesh(m,K);
      // m.write_to_file("toto.mesh");
      
      m.optimize_structure();


//       getfem::vtk_export exp("toto4.vtk");
//       exp.exporting(m);
//       exp.write_mesh_quality(m);

//       getfem::stored_mesh_slice sl;
//       sl.build(m, getfem::slicer_explode(0.8), 8);
//       getfem::vtk_export exp2("totoq.vtk");
//       exp2.exporting(sl);
//       exp2.write_mesh();
//       exp2.write_mesh_quality(m);
//       getfem::dx_export exp3("totoq.dx");
//       exp3.exporting(sl);
//       exp3.write_mesh();




//       getfem::stored_mesh_slice slb; slb.build(m, getfem::slicer_boundary(m), 4);
//       getfem::stored_mesh_slice sl2;
//       getfem::mesh_slicer ms(m); 


//       getfem::slicer_build_stored_mesh_slice bb(sl2);
//       ms.push_back_action(bb);
//       getfem::convex_face_ct cvlst;
//       for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv) {
// 	scalar_type q = m.convex_quality_estimate(cv);
// 	if (q< 0.2) 
// 	  cvlst.push_back(getfem::convex_face(cv));
// 	//cout << "cv " << cv << ": q=" << q << "\n";
//       }
//       //ms.exec(3, cvlst);
//       sl2.merge(slb);
      

//       getfem::vtk_export exp3("totoq2.vtk");
//       exp3.exporting(sl2);
//       exp3.write_mesh();
//       exp3.write_mesh_quality(m);
    }
    

  };


  const mesher_half_space void_signed_distance(base_node(0.0, 0.0),
					       base_small_vector(0.0, 1.0));



  void build_mesh(mesh &m, const mesher_signed_distance& dist_,
		  scalar_type h0, const std::vector<base_node> &fixed_points,
		  size_type K, int noise, size_type iter_max, int prefind,
		  scalar_type dist_point_hull,
		  scalar_type boundary_threshold_flatness) {
    mesher mg(K, dist_, getfem::mvf_constant(1), h0, m, fixed_points, noise,
	      iter_max, prefind, dist_point_hull, boundary_threshold_flatness);
  }
  

  // ******************************************************************
  //    Interface with qhull
  // ******************************************************************

# ifndef GETFEM_HAVE_QHULL_QHULL_H
  void delaunay(const std::vector<base_node> &,
		gmm::dense_matrix<size_type>&) {
    GMM_ASSERT1(false, "Qhull header files not installed. "
		"Install qhull library and reinstall Getfem++ library.");
  }
# else

  extern "C" {
#ifdef _MSC_VER
# include <libqhull/qhull_a.h>
#else
# include <qhull/qhull.h>
//# include <qhull/mem.h>
# include <qhull/qset.h>
//# include <qhull/geom.h>
//# include <qhull/merge.h>
//# include <qhull/poly.h>
//# include <qhull/io.h>
//# include <qhull/stat.h>
#endif
}

  void delaunay(const std::vector<base_node> &pts,
		gmm::dense_matrix<size_type>& simplexes) {
    // cout << "running delaunay with " << pts.size() << " points\n";
    size_type dim = pts[0].size();   /* points dimension.           */
    if (pts.size() <= dim) { gmm::resize(simplexes, dim+1, 0); return; }
    if (pts.size() == dim+1) {
      gmm::resize(simplexes, dim+1, 1);
      for (size_type i=0; i <= dim; ++i) simplexes(i, 0) = i;
      return;
    }
    std::vector<coordT> Pts(dim * pts.size());
    for (size_type i=0; i < pts.size(); ++i)
      gmm::copy(pts[i], gmm::sub_vector(Pts, gmm::sub_interval(i*dim, dim)));
    boolT ismalloc=0;  /* True if qhull should free points in
			* qh_freeqhull() or reallocation      */
    /* Be Aware: option QJ could destabilizate all, it can break everything. */
    /* option Qbb -> QbB (????) */
    /* option flags for qhull, see qh_opt.htm */
    char flags[]= "qhull QJ d Qbb Pp T0"; //QJ s i TO";//"qhull Tv";
    FILE *outfile= 0;    /* output from qh_produce_output()
                          *  use NULL to skip qh_produce_output() */ 
    FILE *errfile= stderr;    /* error messages from qhull code */ 
    int exitcode;             /* 0 if no error from qhull */
    facetT *facet;	          /* set by FORALLfacets */
    int curlong, totlong;	  /* memory remaining after qh_memfreeshort */
    vertexT *vertex, **vertexp;
    exitcode = qh_new_qhull (int(dim), int(pts.size()), &Pts[0], ismalloc,
                             flags, outfile, errfile);
    if (!exitcode) { /* if no error */ 
      size_type nbf=0;
      FORALLfacets { if (!facet->upperdelaunay) nbf++; }
      gmm::resize(simplexes, dim+1, nbf);
	/* 'qh facet_list' contains the convex hull */
      nbf=0;
      FORALLfacets {
        if (!facet->upperdelaunay) {
	  size_type s=0;
          FOREACHvertex_(facet->vertices) {
	    assert(s < (unsigned)(dim+1));
            simplexes(s++,nbf) = qh_pointid(vertex->point);
	  }
	  nbf++;
        }
      }
    }
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort (&curlong, &totlong);
    if (curlong || totlong)
      cerr << "qhull internal warning (main): did not free " << totlong << 
        " bytes of long memory (" << curlong << " pieces)\n";

  }

#endif

}

