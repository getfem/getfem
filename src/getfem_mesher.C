/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_mesher.C : experimental mesher.                       */
/*     									   */
/*                                                                         */
/* Date : May 1, 2004.                                                     */
/* Author : Julien Pommier, Julien.Pommier@insa-toulouse.fr                */
/*          Yves Renard, Yves.Renard@insa-toulouse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2004    Julien Pommier, Yves Renard.                      */
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

#include <getfem_mesher.h>

namespace getfem {

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
	if (constraints.last_true() > other.constraints.last_true()) return false;
	else if (constraints.last_true() < other.constraints.last_true()) return true;
	else if (constraints.card() > other.constraints.card()) return true;
	else if (constraints.card() < other.constraints.card()) return false;
	else for (dal::bv_visitor i1(constraints), i2(other.constraints); !i1.finished(); ++i1, ++i2) {
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
    scalar_type h0;
    size_type N,K;
    base_node bounding_box_min, bounding_box_max;

    std::vector<base_node> pts, pts_prev;
    std::vector<const pt_attribute*> pts_attr;
    std::set<pt_attribute> attributes_set;
    gmm::dense_matrix<size_type> t;    
    
    scalar_type ptol, ttol, L0mult, deltat, geps, deps;

    std::vector<const mesher_signed_distance*> constraints;

    gmm::dense_matrix<scalar_type> W;

    mesher(size_type K_,
	   const mesher_signed_distance& dist_, 
	   const mesher_virtual_function& edge_len_, 
	   scalar_type h0_,
	   getfem_mesh &m, const std::vector<base_node> &fixed_points)
      : dist(dist_), edge_len(edge_len_) {
      K=K_; h0=h0_;
      ptol = 0.001;
      ttol = .1;
      dist.bounding_box(bounding_box_min,bounding_box_max);
      N = bounding_box_min.size();
      if (N == 2) { 
	L0mult = 1.2; deltat = .2; geps = .001*h0; 
      } else {
	L0mult=1+.4/pow(2,(N-1)); deltat=.1; geps=1e-1*h0;
      }
      deps=sqrt(1e-8)*h0;
      dist.register_constraints(this->constraints);

      bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(N,1);
      gmm::resize(W,N,N);
      base_matrix G(N,N+1); 
      vectors_to_base_matrix(G, bgeot::equilateral_simplex_of_reference(N)->points());
      gmm::mult(G, bgeot::geotrans_precomp(pgt,&pgt->convex_ref()->points())->grad(0), W);
      gmm::lu_inverse(W);
      run(m,fixed_points);
    }

    template<class TAB> scalar_type simplex_quality(const TAB &ppts) {
      base_matrix G(N,N), GW(N, N);
      base_matrix::iterator it = G.begin();
      for (size_type i=0; i < N; ++i) {
	base_node P = ppts[i+1] - ppts[0];
	std::copy(P.const_begin(), P.const_begin()+N, G.begin()+i*N);
      }
      gmm::mult(G, W, GW);
      return dal::abs(1./gmm::condition_number(GW));
    }

    scalar_type worst_element, best_element;
    
    scalar_type fbcond_cost_function(const base_vector &c) {
      unsigned nbt = gmm::mat_ncols(t);
      scalar_type cost = 0;
      base_matrix S(N,N), SW(N,N);
      worst_element = 1; best_element = 1E40;
      for (unsigned i=0; i < nbt; ++i) {
	for (size_type j=0; j < N; ++j) {
	  for (size_type k=0; k < N; ++k) {
	    S(k,j) = c[t(j+1,i)*N+k] - c[t(0,i)*N+k];
	  }
	}
	gmm::mult(S,W,SW);
	if (gmm::lu_det(SW) < 0) cost += 1e30;
	else {
	  scalar_type qual = gmm::Frobenius_condition_number_sqr(SW);
	  cost += qual;
// 	  if (qual > 20000.0) { 
// 	    cout << "element " << i << "is horrible" << endl;
// 	    cout << "real precond : " << gmm::condition_number(SW) << endl;
// 	    for (uint k = 0; k < N+1; ++k)
// 	      cout << "point " << k << " : " << pts[t(k, i)] << endl;
// 	  }
	  bool att = true;
	  for (uint k = 0; k < N+1; ++k)
	    if (!(pts_attr[t(k, i)]->constraints.card() || pts_attr[t(k, i)]->fixed))
	      att = false;
//	  if (att) cout << "Element " << i << "is comp. attached\n";

	  worst_element = std::max(worst_element, qual / (N*N));
	  best_element = std::min(best_element, qual / (N*N));
	}
      }
      return cost / scalar_type(N * N);
    }

    void fbcond_cost_function_derivative(const base_vector& c,
					 base_vector &grad) {
      gmm::clear(grad);
      base_matrix Dcond(N,N), G(N,N), S(N,N), SW(N,N);
      unsigned nbt = gmm::mat_ncols(t);
      
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
//	    cout << " " << grad[i*N+k];
    	  }
       }
      gmm::scale(grad, scalar_type(1) / scalar_type(N * N));
      
    }

    struct fbcond_cost_function_object {
      mesher &m;
      fbcond_cost_function_object(mesher &m_) : m(m_) {}
      scalar_type operator()(const base_vector& c) { return m.fbcond_cost_function(c); }
    };

    struct fbcond_cost_function_derivative_object {
      mesher &m;
      fbcond_cost_function_derivative_object(mesher &m_) : m(m_) {}
      void operator()(const base_vector& c, base_vector &grad) { m.fbcond_cost_function_derivative(c, grad); }
    };

    void optimize_quality() {

      size_type nbt = gmm::mat_ncols(t);

      cout << "Opt qual ...\n";
      base_vector X(pts.size() * N);
      for (unsigned i=0; i < pts.size(); ++i)
	dal::copy_n(pts[i].const_begin(), N, X.begin() + i*N);

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
	if (dal::abs(lu_det(S)) < 1e-10)
	  cout << "oulala " << i << ", " << dal::abs(lu_det(S)) << "\n";
	gmm::mult(S,W,SW);
	
	assert(gmm::lu_det(S) >= 0);
	assert(gmm::lu_det(SW) >= 0);
      }
      
      cout << "Initial quality: " << fbcond_cost_function(X)/nbt << "\n";
      cout << "best element : " << sqrt(best_element) << " worst element : "
	   << sqrt(worst_element) << endl;
      gmm::iteration iter; iter.set_noisy(1); iter.set_maxiter(1000);
      iter.set_resmax(1E-2);
      gmm::bfgs(fbcond_cost_function_object(*this), 
		fbcond_cost_function_derivative_object(*this),
		X, 20, iter, 0, 0.001, float(gmm::mat_ncols(t)));

      cout << "Final quality: " << fbcond_cost_function(X)/nbt << "\n";
      cout << "best element : " << sqrt(best_element) << " worst element : "
	   << sqrt(worst_element) << endl;

      for (unsigned i=0; i < pts.size(); ++i)
	dal::copy_n(X.begin() + i*N, N, pts[i].begin());      
    }

    void projection(base_node &X) { 
      scalar_type d = dist(X);
      if (d > 0) { 
	while (dal::abs(d) > 1e-10) {
	  X -= d*dist.grad(X); 
	  d=dist(X);
	}
      }
    }
    
    void projection(base_node &X, dal::bit_vector& bv) {
      projection(X);
      bv.clear();
      dist(X,bv);
    }


    void constraint_projection(base_node &X, size_type cnum) {
      scalar_type d = (*constraints[cnum])(X);
      while (dal::abs(d) > 1e-10) {
	X -= d*(*constraints[cnum]).grad(X); 
	d=(*constraints[cnum])(X);
      }
    }

    bool pure_multi_constraint_projection(base_node &X,
					  const dal::bit_vector &cts) {
      base_node oldX;
      size_type cnt = 0;
      do {     
	oldX = X;
	for (dal::bv_visitor ic(cts); !ic.finished(); ++ic)
	  constraint_projection(X,ic);
	++cnt;
      } while (cts.card() && gmm::vect_dist2(oldX,X) > 1e-16 && cnt < 10000);
      if (cnt >= 10000) {
	cout << "c'est dur : X=" << X << ", cts=" << cts << "dist = "
	     << gmm::vect_dist2(oldX,X) << "\n"; getchar();
	return false;
      }
      if (cts.card()) {
	dal::bit_vector ct2;
	dist(X, ct2);
	return ct2.contains(cts);
      }
      return true;
    }

    bool multi_constraint_projection(base_node &X,
				     const dal::bit_vector &cts) {
      base_node oldX;
      size_type cnt = 0;
      do {     
	oldX = X;
	for (dal::bv_visitor ic(cts); !ic.finished(); ++ic)
	  constraint_projection(X,ic);
	projection(X);
	++cnt;
      } while (cts.card() && gmm::vect_dist2(oldX,X) > 1e-16 && cnt < 10000);
      if (cnt > 10000) { 
	cout << "c'est dur : X=" << X << ", cts=" << cts << "dist = "
	     << gmm::vect_dist2(oldX,X) << "\n"; getchar();
	return false;
      }
      if (cts.card()) {
	dal::bit_vector ct2;
	dist(X, ct2);
	return ct2.contains(cts);
      }
      return true;
    }

    void tangential_displacement(base_small_vector &X, const dal::bit_vector &cts) {
      base_matrix normals(N, cts.card());
      size_type cnt = 0;
      for (dal::bv_visitor ic(cts); !ic.finished(); ++ic) {
	gmm::copy(constraints[ic]->grad(X), gmm::mat_col(normals, cnt));
	for (size_type k=0; k < cnt; ++k) {
	  gmm::add(gmm::scaled(gmm::mat_col(normals, k), -gmm::vect_sp(gmm::mat_col(normals,k),gmm::mat_col(normals,cnt))),
		   gmm::mat_col(normals,cnt));
	  scalar_type n = gmm::vect_norm2(gmm::mat_col(normals, cnt));
	  if (n < 1e-8) continue;
	  gmm::scale(gmm::mat_col(normals, cnt), 1./n);
	  gmm::add(gmm::scaled(gmm::mat_col(normals, cnt), -gmm::vect_sp(X, gmm::mat_col(normals, cnt))), X);
	  ++cnt;
	}
      }
    }

    const pt_attribute*
    get_attr(bool fixed, const dal::bit_vector &bv) {
      pt_attribute a; 
      a.fixed = fixed;
      a.constraints = bv; 
      return &(*attributes_set.insert(a).first);
    }

    void project_and_update_constraints(size_type ip) {
      const dal::bit_vector& cts = pts_attr[ip]->constraints;
      dal::bit_vector new_cts;
      base_node oldX = pts[ip];
      multi_constraint_projection(pts[ip], cts);
      dist(pts[ip], new_cts);
      if (new_cts.card() > cts.card()) {
	cout << "Point #" << ip << " " << oldX << "->" << pts[ip]
	     << " has been upgraded from " 
	     << cts << " to " << new_cts
	     << ", congratulations to the winner\n";
	pts_attr[ip] = get_attr(pts_attr[ip]->fixed, new_cts);
      }
      if (!new_cts.contains(cts)) {
	cout << "\n\n\nPoint #" << ip << " " << oldX << "->" << pts[ip]
	     << " is doing n'importe quoi: "  << cts << " to " << new_cts
	     << ", congratulations to the looser\n\n\n"; // getchar();
	//assert(0);
      }
    }
    
    template <class VECT> void move_carefully(size_type ip, const VECT &V) {
      scalar_type norm = gmm::vect_norm2(V);
      if (norm <= h0/scalar_type(6))
	gmm::add(V, pts[ip]);
      else
	gmm::add(gmm::scaled(V, h0 / (scalar_type(6) * norm)), pts[ip]);
      project_and_update_constraints(ip);
    }

     template <class VECT> void move_carefully(const VECT &V) {
       size_type npt = gmm::vect_size(V) / N;
       for (size_type i = 0; i < npt; ++i)
	 move_carefully(i, gmm::sub_vector(V, gmm::sub_interval(i*N, N)));
     }

    void distribute_points_regularly(getfem_mesh &m, 
			      const std::vector<base_node> &fixed_points) {
      size_type nbpt = 1;
      std::vector<size_type> gridnx(N);

      for (size_type i=0; i < N; ++i) 
	h0 = std::min(h0, bounding_box_max[i] - bounding_box_min[i]);
      if (h0 < 1E-20) { 
	cout << "h0 = " << h0 << " too small, aborting.\n";
	return;
      }
      for (size_type i=0; i < N; ++i) {
	scalar_type h = h0;
	if (N == 2 && i == 1) h = sqrt(3.)/2. * h0;
	gridnx[i]=1+(size_type)((bounding_box_max[i]-bounding_box_min[i])/h);
	nbpt *= gridnx[i];
      }
      m.clear();

      /* build the regular grid and filter points outside */
      for (size_type i=0; i < fixed_points.size(); ++i) {
	if (dist(fixed_points[i]) < geps &&
	    m.search_point(fixed_points[i]) == size_type(-1)) {
	  m.add_point(fixed_points[i]); 
	  pts.push_back(fixed_points[i]); 
	  pts_attr.push_back(get_attr(true,dal::bit_vector()));
	} else cout << "removed duplicate fixed point : " << fixed_points[i] << "\n";
      }
      cout << "gridnx=" << gridnx << "\n";
      base_node P(N), Q(N);
      for (size_type i=0; i < nbpt; ++i) {
	// cout << "dealing with point " << i << " on " << nbpt << " : ";
	for (size_type k=0, r = i; k < N; ++k) {
	  unsigned p =  r % gridnx[k];
	  P[k] = p * (bounding_box_max[k] - bounding_box_min[k]) / 
	    (gridnx[k]-1) + bounding_box_min[k];
	  if (N==2 && k==0 && ((r/gridnx[0])&1)==1) P[k] += h0/2;
	  r /= gridnx[k];
	}

	dal::bit_vector co;
	// if (dist(P) < 0) {  
	if (false) {
	  for (size_type k = 0; k < constraints.size() && co.card() < N; ++k) {
	    gmm::copy(P, Q);
	    if (gmm::abs((*(constraints[k]))(Q)) < h0) {
	      constraint_projection(Q, k);
	      if (dist(Q) > -geps && 
		  gmm::vect_dist2(P, Q) < h0 / scalar_type(2)) co.add(k);
	    }
	  }
	}
	gmm::copy(P, Q);
	if (co.card() > 0) { 
	  bool ok = pure_multi_constraint_projection(Q, co);
	  if (!ok || gmm::abs(dist(Q)) > geps) { gmm::copy(P, Q); co.clear(); }
	}

	if (dist(Q) < geps) {
	  if (m.search_point(Q) == size_type(-1)) {
	    size_type q = m.add_point(Q); pts.push_back(Q);
	    pts_attr.push_back(get_attr(false, co));
	    if (co.card()) cout << "adding points with constraints : "
				<< co << " : " << q << endl;
	  }
	}
      }
    }

    void add_point_hull(void) { 
      scalar_type dist_hull(3);
      size_type nbpt = pts.size(), nbadd(0);
      base_node P, Q;
      base_small_vector V;
      for (unsigned i=0; i < nbpt; ++i) {
	if (pts_attr[i]->constraints.card()) {
	  P = pts[i];
	  V = dist.grad(P);
	  scalar_type d = gmm::vect_norm2(V);
	  if (d > 0) {
	    P += V * (dist_hull*h0/d);
	    if (dist(P)*sqrt(scalar_type(N)) > dist_hull*h0) {
	    Q = P;
	    projection(Q);
	    if (gmm::vect_dist2(P, Q) > dist_hull*h0/scalar_type(2))
	      { pts.push_back(P); ++nbadd; }
	    }
	  }
	}
      }
      cout << "point hull : " << nbadd << " points added\n";
    }


    scalar_type pts_dist_max(const std::vector<base_node> &A, 
		      const std::vector<base_node> &B) {
      scalar_type dist_max = 0;
      for (size_type i=0; i < pts.size(); ++i) 
	dist_max = std::max(dist_max, bgeot::vect_dist2_sqr(A[i],B[i]));
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
      unsigned nb_bb_pts = 0;
      for (size_type i=0; i < idx.size(); ++i) idx[i] = i;
   
      std::sort(idx.begin(), idx.end(), cleanup_points_compare(pts,pts_attr));
      bgeot::kdtree tree;
      bgeot::kdtree_tab_type neighbours;
      dal::bit_vector keep_pts; keep_pts.add(0,idx.size());
      cout << "cleanup points : in the beginning there were " << pts.size()
	   << " points ( " << nb_bb_pts << " on boundary)\n";
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
		cout << "point #" << i << " " << P << " is too near from point #" 
		     << neighbours[k].i << pts[idx[neighbours[k].i]] << " : will be removed\n"; getchar();
               keep_pts.sup(i);
	      }
	    }
	  }
	  tree.clear();
	  i0 = i+1;
	}
      }
      cout << "cleanup points : at the end, only " << keep_pts.card() << " remain\n";      
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
      /*pts_prev.resize(pts.size());
	std::copy(pts.begin(), pts.end(), pts_prev.begin());*/
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
	bool is_bridge_simplex = false;
	scalar_type q(0), dWG(0), dG(0);
	base_node G;
	
	for (size_type k=0; k <= N; ++k)
	  if (t(k, i) >= nbpt) ext_simplex = true;

	if (!ext_simplex) {
	  G = pts[t(0,i)];
	  for (size_type k=1; k <= N; ++k) G += pts[t(k,i)];
	  gmm::scale(G, 1./(N+1));
	  dG = dist(G);
	  gmm::clear(weights);
//  	  scalar_type tot(0);
//  	  for (size_type k = 0; k <= N; ++k) {
//  	    for (size_type l = 0; l <= N; ++l) {
//  	      scalar_type d = gmm::vect_dist2(pts[t(k,i)], pts[t(l,i)]);
//  	      weights[k] += d; tot +=d;
//  	    }
//  	  }
//  	  gmm::clear(G);
//  	  for (size_type k = 0; k <= N; ++k)
//  	    gmm::add(gmm::scaled(pts[t(k,i)], weights[k] / tot), G);
//  	  dWG = dist(G);

	  dWG = -1.0; // desactivated
	  
	  q = simplex_quality(dal::index_ref_iterator(pts.begin(),
						      gmm::mat_col(t,i).begin()));
	  
	  for (size_type k=0; k <= N; ++k) {
	    if (pts_attr[t(k,i)]->constraints.card() == 0) {
	      boundary_simplex = false;
	    }
	  }
	  if (boundary_simplex) {
	    dal::bit_vector all_cts;
	    for (size_type k=0; k < N+1; ++k) 
	      all_cts |= pts_attr[t(k,i)]->constraints;
	    for (size_type k=0; k < N+1; ++k)
	      if (pts_attr[t(k,i)]->constraints.contains(all_cts))
		is_bridge_simplex = true;
	  }
	  //qf << q << "\n";
	}
	if ( (version == 0) && (ext_simplex || q < 1e-2 || dG > 0 || dWG > 0
				// || (is_bridge_simplex && dG>-geps)
				)
	     || ((version == 1) && (ext_simplex || dG > 0 || dWG > 0 || q < 1e-30))) {
	  if (i != gmm::mat_ncols(t)-1) {
	    for (size_type k=0; k < N+1; ++k)
	      std::swap(t(k,i), t(k,gmm::mat_ncols(t)-1));
	  }
	  t.resize(N+1,gmm::mat_ncols(t)-1);
	} else {
	  ++i;
	  if (q < worst_q) { worst_q = q; worst_q_P = G*(1./(N+1)); }
	}
      }
      
    }


    void run(getfem_mesh &m,
	     const std::vector<base_node> &fixed_points) {

      distribute_points_regularly(m, fixed_points);
      
      std::vector<base_node> pts2(pts.size(),base_node(N));    
      bool first = true;
      bgeot::mesh_structure edges_mesh;
      base_vector L, L0;
      scalar_type maxdp;
      size_type count = 0, nbpt;
      do {
	scalar_type dist_max = first ? 0 : pts_dist_max(pts, pts_prev);
	cout << "**** Iter " << count << ", dist_max since last delaunay ="
	     << dist_max << ", tol=" << ttol*h0 << "\n";
	if ((dist_max > ttol*h0) || first) {
	  cout << "-----> NEW DELAUNAY\n";
	  first = false;

	  //	  for (unsigned repeat=0; repeat < 2; ++repeat) { // pourquoi 2 ?
	    cleanup_points(); /* and copy pts to pts_prev */
	    cout << "running delaunay on " << pts.size() << " points\n";
	    nbpt = pts.size();
	    add_point_hull();
	    delaunay(pts, t);
	    pts.resize(nbpt);
	    cout << "nb splx avant suppr = " << gmm::mat_ncols(t) << "\n";
	    select_elements(0);
	    // optimize_quality();
	    cout << " worst q = " << worst_q << "\n";
	    //	  }
	  edges_mesh.clear();
	  cout << "nb splx = " << gmm::mat_ncols(t) << "\n";
	  for (size_type i=0; i < gmm::mat_ncols(t); ++i)
	    for (size_type j=0; j < N+1; ++j)
	      for (size_type k=j+1; k < N+1; ++k)
		edges_mesh.add_segment(t(j,i), t(k,i));
	  //edges_mesh.optimize_structure()	
	}

	size_type nbcv = edges_mesh.convex_index().card();
	if (nbcv == 0) DAL_THROW(dal::failure_error, "no more edges!");
	L.resize(nbcv); L0.resize(nbcv);
	scalar_type sL = 0, sL0 = 0;
	size_type nbBoundary = 0; scalar_type sLboundary = 0.;
	dal::bit_vector boundary_edges;
	boundary_edges.sup(0,edges_mesh.convex_index().last_true());
	for (dal::bv_visitor ie(edges_mesh.convex_index());
	     !ie.finished(); ++ie) {
	  const base_node &A = pts[edges_mesh.ind_points_of_convex(ie)[0]];
	  const base_node &B = pts[edges_mesh.ind_points_of_convex(ie)[1]];
	  L[ie] = gmm::vect_norm2(B-A);
	  L0[ie] = edge_len(.5*(A+B));
	  sL += pow(L[ie],N);
	  sL0 += pow(L0[ie],N);
                
	  if (dist(A) >= -geps && dist(B) >= -geps)
	    { nbBoundary++; sLboundary += L[ie]; boundary_edges.add(ie); }
	}
	if (nbBoundary) sLboundary /= nbBoundary; 
	cout << "L diff = " << pow(sL/sL0, 1./N)
	     << " Lboundary=" << sLboundary/h0 << "\n";

	gmm::scale(L0, L0mult * pow(sL/sL0, 1./N));
	//cout << "L0 = " << L0 << "\n";
	pts2.resize(pts.size());
	std::copy(pts.begin(),pts.end(),pts2.begin());

	// Move the points

	base_vector X(pts.size() * N);
	
	for (dal::bv_visitor ie(edges_mesh.convex_index());
	     !ie.finished(); ++ie) {
	  size_type iA = edges_mesh.ind_points_of_convex(ie)[0];
	  size_type iB = edges_mesh.ind_points_of_convex(ie)[1];
	  base_node bar = pts2[iB]-pts2[iA];
	  scalar_type F = std::max(L0[ie]-L[ie], 0.);
	  // scalar_type F = L0[ie]-L[ie]; if (F < 0.) F /= 10.0;
	  
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
	move_carefully(X);






	// à ne faire que de temps en temps.

	if (count % 10 == 0) {

	  bgeot::kdtree tree;
	  bgeot::kdtree_tab_type neighbours;
	  bool tree_empty = true;
	  
	  
	  for (dal::bv_visitor ie(edges_mesh.convex_index());
	     !ie.finished(); ++ie) {
	    size_type iA = edges_mesh.ind_points_of_convex(ie)[0];
	    size_type iB = edges_mesh.ind_points_of_convex(ie)[1];
	    if (L[ie] > L0[ie]) continue;
	    if (pts_attr[iA] == pts_attr[iB] ||
		pts_attr[iA]->constraints.card() == 0 ||
		pts_attr[iB]->constraints.card() == 0) continue;
	    if (pts_attr[iA]->constraints == pts_attr[iB]->constraints) continue;
	    dal::bit_vector bv1(pts_attr[iA]->constraints);
	    bv1.setminus(pts_attr[iB]->constraints);
	    dal::bit_vector bv2(pts_attr[iB]->constraints);
	    bv2.setminus(pts_attr[iA]->constraints);
	    if (bv1.card() && bv2.card()) {
	      bv1 |= bv2;
	      bgeot::mesh_point_search_ind_ct
		iAneighbours = edges_mesh.ind_points_to_point(iA);
	      bgeot::mesh_point_search_ind_ct
		iBneighbours = edges_mesh.ind_points_to_point(iB);
	      std::vector<size_type>
		common_pts(iAneighbours.size()+iBneighbours.size());
	      std::sort(iAneighbours.begin(),iAneighbours.end());
	      std::sort(iBneighbours.begin(),iBneighbours.end());
	      std::vector<size_type>::iterator ite = 
		std::set_intersection(iAneighbours.begin(), iAneighbours.end(),
				      iBneighbours.begin(), iBneighbours.end(),
				      common_pts.begin());
	      common_pts.resize(ite-common_pts.begin());
	      bool do_projection = true;
	      if (dist(.5*(pts[iA]+pts[iB])) < 0) {
		for (std::vector<size_type>::iterator it = common_pts.begin();
		     it != ite; ++it) {
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
		    { bmin[k] -= h0/3.; bmax[k] += h0/3.; }
		  tree.points_in_box(neighbours, bmin, bmax);
		  for (size_type k=0; k < neighbours.size(); ++k) {
		    if (neighbours[k].i != iA && neighbours[k].i != iB)
		      do_projection = false;
		    cout << "neighbours found : " << neighbours[k].i
			 << " iA = " << iA << " iB = " << iB <<endl; 
		  }
		  
		  if (do_projection) {
		    cout << "cts_A = " << pts_attr[iA]->constraints
			 << ", cts_B = " << pts_attr[iB]->constraints
			 << ", L[ie]=" << L[ie] << "\n";
		    cout << bv1 << bv2 << " Promotion !!!!!! pt#" << iA << " "
			 << pts[iA] << " : " 
			 << pts_attr[iA]->constraints << "->" << bv1 << "\n";
		    
		    pts_attr[iA] = get_attr(pts_attr[iA]->fixed, bv1);
		    project_and_update_constraints(iA);
		    cout << " ---> new value: "<<pts[iA]<< " PA = "<<PA<<"\n";
		  }
		}
	    }
	    }
	  }
	}

	/* 7. Bring outside points back to boundary && eval term criterion */
	maxdp = pts_dist_max(pts,pts2);
	cout << "count = " << count << ", maxdp = " << maxdp << ", ptol = "
	     << ptol << " CV=" << sqrt(maxdp)*deltat/(ptol * h0) << "\n";
	++count;
	// m.clear();
	// for (size_type i=0; i < t.size()/(N+1); ++i)
	//  m.add_convex_by_points(bgeot::simplex_geotrans(N,1),
	//		 dal::index_ref_iterator(pts.begin(), &t[i*(N+1)]));
	// char s[50]; sprintf(s, "toto%02d.mesh", count);
	// m.write_to_file(s);
      } while ((count < 40 || sqrt(maxdp)*deltat > ptol * h0) && count<1000);

      { m.clear();
	build_simplex_mesh(m,K);
	m.write_to_file("toto.mesh");
      }

      m.optimize_structure();
      {
        getfem::vtk_export exp("toto2.vtk");
	exp.exporting(m);
        exp.write_mesh_quality(m);
      }
      nbpt = pts.size();
      add_point_hull();
      delaunay(pts, t);
      pts.resize(nbpt);
      select_elements(0);
      optimize_quality();
       { m.clear();
 	build_simplex_mesh(m,K);
 	m.write_to_file("toto.mesh");
       }

       m.optimize_structure();
       {
         getfem::vtk_export exp("toto3.vtk");
 	exp.exporting(m);
         exp.write_mesh_quality(m);
       }

      cout << "mesh done! (" << count << " iter )\n";
      base_vector simplex_q; simplex_q.reserve(m.convex_index().card());
      for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv) {
	simplex_q.push_back(simplex_quality(dal::index_ref_iterator(pts.begin(), t.begin()+cv*(N+1))));
      }
      scalar_type mean_q = dal::mean_value(simplex_q.begin(), simplex_q.end());
      scalar_type std_dev_q = 0, max_q = 0., min_q = 1.;
      for (size_type i=0; i < simplex_q.size(); ++i) {
	std_dev_q += dal::sqr(simplex_q[i]-mean_q); 
	max_q = std::max(max_q,simplex_q[i]); 
	min_q = std::min(min_q, simplex_q[i]);
      }
      std_dev_q = sqrt(std_dev_q/simplex_q.size());
      cout << "number of simplexes : " << simplex_q.size() << "\n";
      cout << "best quality : " << max_q << "\n";
      cout << "worst quality : " << min_q << "\n";
      cout << "mean quality : " << mean_q << ", std_dev = " << std_dev_q << "\n";
      m.optimize_structure();
    }

    void build_simplex_mesh(getfem_mesh &m, size_type degree) {
      std::vector<base_node> cvpts(N+1), cvpts2;
      size_type cvnum;
      m.clear();
      for (size_type ip=0; ip < pts.size(); ++ip) {
	size_type z;
	bgeot::small_vector<scalar_type> P = pts[ip];
	while ((z = m.add_point(P)) != ip) {
	  cout << "WARNING : points are too near ...\n";
	  bgeot::small_vector<scalar_type> Z(N); gmm::fill_random(Z);
	  gmm::add(gmm::scaled(Z, h0/1000.0), P);
	}
      }
      for (size_type i=0; i < t.size()/(N+1); ++i) {
	for (size_type k=0; k < N+1; ++k) cvpts[k] = pts[t[i*(N+1)+k]];
	if (degree == 1) {
	  //cvnum = m.add_convex_by_points(bgeot::simplex_geotrans(N,1), cvpts.begin());
	  cvnum = m.add_convex(bgeot::simplex_geotrans(N,1), &t[i*(N+1)]);
	  assert(cvnum == i);
	} else {
	  bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(N,degree);
	  cvpts2.resize(pgt->nb_points());
	  for (size_type k=0; k < pgt->nb_points(); ++k) {
	    cvpts2[k] = bgeot::simplex_geotrans(N,1)->transform(pgt->convex_ref()->points()[k], 
								cvpts);
	  }
	  cvnum = m.add_convex_by_points(pgt, cvpts2.begin());
	  assert(cvnum == i);
	}
      }
//       while (1) {
// 	getfem::convex_face_ct border_faces;
// 	getfem::outer_faces_of_mesh(m, border_faces);
// 	size_type nbrm = 0;
// 	for (getfem::convex_face_ct::const_iterator it = border_faces.begin();
// 	     it != border_faces.end(); ++it) {
// 	  if (!m.convex_index()[it->cv]) continue;
// 	  scalar_type q = simplex_quality(dal::index_ref_iterator(pts.begin(), &t[it->cv*(N+1)]));
// 	  scalar_type dG = dist(dal::mean_value(m.points_of_convex(it->cv).begin(), m.points_of_convex(it->cv).end()));
// 	  if (q < .01) {
// 	    //cout << "removing flat border convex " << it->cv << " (q=" << q << ")\n";
// 	    //for (size_type i=0; i < N+1; ++i) cout << " " << m.points_of_convex(it->cv)[i]; cout << "\n";
// 	    m.sup_convex(it->cv); nbrm++;
// 	  } else if (dG > -0.0) {
// 	    //cout << "removing convex because of its gravity center is too near" << it->cv << " (dG=" << dG << ")\n";
// 	    m.sup_convex(it->cv); nbrm++;  
// 	  }
// 	}
// 	if (nbrm == 0) break;
// 	//else cout << "\n\n !!!!! ON CONTINUE ENCORE UN COUP..\n\n";
//       }
      if (degree>1) {
	//m.optimize_structure();
	getfem::convex_face_ct border_faces;
	getfem::outer_faces_of_mesh(m, border_faces);
	dal::bit_vector ptdone; ptdone.sup(0,m.points_index().last_true());
	for (getfem::convex_face_ct::const_iterator it = border_faces.begin();
	     it != border_faces.end(); ++it) {
	  bgeot::ind_ref_mesh_point_ind_ct fpts_ = m.ind_points_of_face_of_convex(it->cv, it->f);
	  std::vector<size_type> fpts(fpts_.size());
	  std::copy(fpts_.begin(), fpts_.end(), fpts.begin());
	  interpolate_face(m, ptdone, fpts, 
			   m.trans_of_convex(it->cv)->structure()->faces_structure()[it->f]);
	  /*
	  std::copy(fpts_.begin(), fpts_.end(), fpts.begin());
	  cout << "\nBOUNDARY FACE (cv#" << it->cv << ", face#" << it->f << ") : pts = ";
	  for (size_type i=0; i < m.points_of_face_of_convex(it->cv,it->f).size(); ++i) 
	    cout << " [pt#" << fpts[i] << "]:" << m.points_of_face_of_convex(it->cv,it->f)[i];
	  cout << "\n";
	  
	  

	  for (size_type iface=0; iface < fpts.size(); ++iface) {
	    size_type ip = fpts[iface];
	    if (ptdone[ip]) continue;
	    scalar_type d = dist(m.points()[ip]);
	    base_node old = m.points()[ip]; scalar_type oldd = d;
	    while (dal::abs(d) > 1e-3) {
	      m.points()[ip] -= d*dist.grad(m.points()[ip]);
	      d=dist(m.points()[ip]);
	    }
	    if (bgeot::vect_dist2(m.points()[ip],old) > 1e-12) {
	      cout << "projected point " << old << ", d=" << oldd << " (cv#" << it->cv 
		   << ", face#" << it->f << ", face pt #" << iface << ") :";
	      cout << " -> " << m.points()[ip] << ", d=" << d << " -> distance = "<< bgeot::vect_dist2(m.points()[ip],old) << "\n";
	    }
	    ptdone.add(ip);
	  }
	  */
	}
      }
    }



    void interpolate_face(getfem_mesh &m, dal::bit_vector& ptdone, 
			  const std::vector<size_type>& ipts, bgeot::pconvex_structure cvs) {
      if (cvs->dim() == 0) return;
      else if (cvs->dim() > 1) {
	std::vector<size_type> fpts;
	for (size_type f=0; f < cvs->nb_faces(); ++f) {
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
      cout << "interpolate face! cts = " << cts << "\n";
      if (cts.card()) {
	dal::bit_vector new_cts;
	for (size_type i=0; i < ipts.size(); ++i) {
	  if (ipts[i] >= pts.size() && !ptdone[ipts[i]]) { 
	    base_node &P = m.points()[ipts[i]];
	    cout << "interpolation of #"<< ipts[i] << "=" << P << " on " << cts;
	    multi_constraint_projection(P, cts);
	    dist(P, new_cts);
	    cout << " -> new_cts=" << new_cts << ", P=" << P << "=" << m.points()[ipts[i]] << "\n";
	  }
	}
      }
    }
  };


   void build_mesh(getfem_mesh &m, const mesher_signed_distance& dist_, scalar_type h0,
		   const std::vector<base_node> &fixed_points, size_type K) {
     mesher mg(K, dist_, getfem::mvf_constant(1), h0, m, fixed_points);
   }


  // ******************************************************************
  //    Interface with qhull
  // ******************************************************************

extern "C" {
#include <qhull/qhull.h>
#include <qhull/mem.h>
#include <qhull/qset.h>
#include <qhull/geom.h>
#include <qhull/merge.h>
#include <qhull/poly.h>
#include <qhull/io.h>
#include <qhull/stat.h>
}

  void delaunay(const std::vector<base_node> &pts,
		gmm::dense_matrix<size_type>& simplexes) {
    
    if (pts.size() == 0) return;
    int dim = pts[0].size();   /* points dimension.           */
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
    exitcode = qh_new_qhull (dim, pts.size(), &Pts[0], ismalloc,
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

}



#if 0

  // ******************************************************************
  //    Archives ....
  // ******************************************************************

#define REVC 0.
    static scalar_type force(scalar_type l, scalar_type l0) { 
      //return 1./l - 1./l0;
      return l0 - l;
      /*if (l0 - l > 0) return l0-l; else return REVC*(l0-l);
	return std::max(l0 - l, 0.); */
    }
    static scalar_type dforce(scalar_type, scalar_type) {
      return -1;
      //return (l < l0) ? -1 : -REVC;
    }
    static scalar_type hforce(scalar_type, scalar_type) {
      return 0;
    }
    
    void point_newton(const base_vector &L0, size_type ip, 
		      const bgeot::mesh_structure &edges_mesh) {
      const bgeot::mesh_convex_ind_ct edges = edges_mesh.convex_to_point(ip);
      if (edges.empty()) return;
      size_type NN = pts[0].size();
      base_matrix H(NN,NN);
      base_node B(NN); gmm::clear(B);
      scalar_type J = 0;
      scalar_type mL0 = 0;
      size_type nedges = 0;
      dal::bit_vector bv_edges;
      scalar_type maxdl = 100;
      for (bgeot::mesh_convex_ind_ct::const_iterator ie = edges.begin();
	   ie != edges.end(); ++ie, ++nedges) {
	size_type i1 = edges_mesh.ind_points_of_convex(*ie)[0];
	size_type i2 = edges_mesh.ind_points_of_convex(*ie)[1]; 
	if (i1 != ip) std::swap(i1,i2);
	assert(i1 == ip);
	scalar_type l0 = L0[*ie]; maxdl = std::min(maxdl,l0/20);
	bv_edges.add(nedges);
	base_node X(pts[i1]-pts[i2]);
	scalar_type n = gmm::vect_norm2(X), n3 = n*dal::sqr(n);
	scalar_type F = force(n,l0), DF=dforce(n,l0), HF=hforce(n,l0);
	mL0 += l0;
	J += .5 * dal::sqr(F);
	maxdl = std::min(maxdl, dal::abs(n-l0));
	//cout << "ie = " << *ie << ", n = " << n << ", L0=" << l0 << "\n";
	if (n > 1e-10) {
	  /*for (size_type i=0; i < NN; ++i) {
	    H(i,i) += -F/n;
	    for (size_type j=0; j < NN; ++j) {
	    H(i,j) += (l0/n3)*(X[i]*X[j]);
	    }
	    }*/
	  for (size_type i=0; i < NN; ++i) {
	    H(i,i) += F*DF/n;
	    for (size_type j=0; j < NN; ++j) {
	      H(i,j) += (((DF*DF+F*HF)*n - F*DF)/n3)*(X[i]*X[j]);
	    }
	  }
	  B += (-F*DF*1./n)*X;
	}	    
      }
      if (bv_edges.card() == 0) return;
      mL0 /= nedges;
      base_node V(NN), X(NN);
      //gmm::copy(B,V);
      tangential_displacement(B, pts_attr[ip]->constraints);
      gmm::scale(B,100*deltat);

      if (gmm::vect_norm2(B) > h0/10)
	gmm::scale(B, h0/(10));
      pts[ip] += B; return;

      gmm::lu_solve(H,V,B);

      /*cout << "|V|=" << gmm::vect_norm2(V) << ", .2*|B|=" << gmm::vect_norm2(gmm::scaled(B,.2)) 
	<< ", maxdl=" << maxdl << "\n";*/

      if (gmm::vect_norm2(V) > maxdl) { gmm::scale(V, maxdl/gmm::vect_norm2(V)); }

      //cout << "H=" << H << ", B=" << B << ", -> X=" << V << "\n";

      scalar_type step = 1., Jf;
      do {
	X =  pts[ip] + step * V;
	//cout << "J initial : " << J << "\n";
	Jf = 0;
	size_type ecnt = 0;
	for (bgeot::mesh_convex_ind_ct::const_iterator ie = edges.begin();
	     ie != edges.end(); ++ie, ++ecnt) {
	  if (!bv_edges[ecnt]) { /*cout << "SKIP " << ecnt << "\n";*/ continue; }
	  size_type i1 = edges_mesh.ind_points_of_convex(*ie)[0];
	  size_type i2 = edges_mesh.ind_points_of_convex(*ie)[1];
	  if (i1 != ip) std::swap(i1,i2);
	  scalar_type n = gmm::vect_norm2(X-pts[i2]);
	  scalar_type l0 = L0[*ie];
	  scalar_type F = force(n,l0);
	  Jf += .5 * dal::sqr(F);
	}
	step /= 2.0;
      } while (Jf > J && step > 0.0001);
      pts[ip] = X;
    }

  };













// 	else if (ORIGINAL == 2) {
//           scalar_type desth = bgeot::equilateral_simplex_of_reference(N)->points()[N][N-1] * L0mult * pow(sL/sL0, 1./N);
//           for (size_type it=0; it < gmm::mat_ncols(t); ++it) {
//             base_node sG(N); sG.fill(0);
//             for (size_type ip=0; ip < N+1; ++ip) { 
//               sG += pts[t(ip,it)];
//               cout << "t#" << it << "/pt#" << t(ip,it) << " = " << pts[t(ip,it)] << "\n";
//             }
            
//             cout << "t#" << it << " G=" << (1./(N+1))*sG << "\n";
//             for (size_type face=0; face < N+1; ++face) {
//               cout << "t#" << it << "/f#" << face << "\n";
//               std::vector<size_type> fpts(N+1);              
//               dal::bit_vector cts;
//               for (size_type ip=(face+1)%(N+1), cnt = 0; cnt < N+1; ++cnt, ip = (ip+1)%(N+1)) {
//                 fpts[cnt] = t(ip,it);
//                 cout << "fpts[" << cnt << "] = " << fpts[cnt] << " = " << pts[fpts[cnt]] << "\n";
//                 if (cnt == 0) cts = pts_attr[fpts[cnt]]->constraints;
//                 else cts &= pts_attr[fpts[cnt]]->constraints;
//               }
//               if (cts.card()) continue;
//               std::vector<base_small_vector> base(N);
//               for (size_type i=0; i < N; ++i) { 
//                 base_node v = pts[fpts[i+1]] - pts[fpts[0]];
//                 for (size_type j=0; j < i; ++j) 
//                   v -= gmm::vect_sp(v,base[j])*base[j];
//                 scalar_type n = gmm::vect_norm2(v);
//                 if (n > 1e-8) v *= 1./n;
//                 base[i] = v;
//                 cout << "base[" << i << "]=" << base[i] << "\n";
//               }
              
//               base_node faceG = (1./N)*(sG - pts[fpts[0]]);
//               base_node target = faceG + desth * edge_len(faceG) * base[N-1];
//               base_node Fbar = target - pts[fpts[N]];
//               cout << "faceG=" << faceG << ", target = " << target << ", Fbar = " << Fbar << "\n";
//               scalar_type d = gmm::vect_norm2(Fbar);
//               if (d > 1e-8) 
//                 pts[fpts[N]] += 0.1*deltat * Fbar;
//             }
//           }
//           for (size_type ip=0; ip < pts.size(); ++ip) {
// 	    project_and_update_constraints(ip);
// 	  }
// 	} 
// 	else {
// 	  //for (size_type ip=0; ip < pts.size(); ++ip) {
// 	  for (size_type ip=pts.size()-1; ip != size_type(-1); --ip) {
// 	    if (pts_attr[ip]->fixed) continue;
// 	    point_newton(L0,ip, edges_mesh);
// 	    project_and_update_constraints(ip);
// 	  }
// 	}

#endif
