/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_mesher : experimental mesher.                         */
/*     									   */
/*                                                                         */
/* Date : May 1, 2004.                                                     */
/* Author : Julien Pommier, Julien.Pommier@insa-toulouse.fr                */
/*          Yves Renard, Yves.Renard@insa-toulouse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2004  Julien Pommier, Yves Renard.                        */
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


#include <getfem_mesh.h>
#include <gmm_condition_number.h>
#include <bgeot_comma_init.h>
#include <gmm_solver_bfgs.h>
#include <getfem_export.h>
#include <bgeot_kdtree.h>
#include <typeinfo>

namespace getfem {

  class mesher_virtual_function {
  public:
    virtual scalar_type operator()(const base_node &P) const = 0;
    virtual ~mesher_virtual_function() {}
  };

  class mvf_constant : public mesher_virtual_function {
    scalar_type c;
  public: 
    mvf_constant(scalar_type c_) : c(c_) {}
    scalar_type operator()(const base_node &) const { return c; }
  };

  // Signed distance definition. Should not be a real signed distance but
  // should satisfy that dist(P) is less than the euclidean distance from
  // P to the boundary and more than 1/sqrt(N) time this distance.

#define SEPS 1e-5

  class mesher_signed_distance : public mesher_virtual_function {
  protected:
    mutable size_type id;
  public:
    mesher_signed_distance() : id(size_type(-1)) {}
    virtual bool bounding_box(base_node &bmin, base_node &bmax) const = 0;
    virtual scalar_type operator()(const base_node &P,
				   dal::bit_vector &bv) const = 0;
    virtual base_small_vector grad(const base_node &P) const = 0;
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const=0;
    virtual scalar_type operator()(const base_node &P) const  = 0;
  };

  class mesher_half_space : public mesher_signed_distance {
    base_node x0; base_node n; scalar_type xon;
  public:
    mesher_half_space(base_node x0_, base_node n_) : x0(x0_), n(n_)
    { n /= gmm::vect_norm2(n); xon = bgeot::vect_sp(x0, n); }
    bool bounding_box(base_node &, base_node &) const
    { return false; }
    virtual scalar_type operator()(const base_node &P) const
    { return xon - bgeot::vect_sp(P,n); }
    virtual scalar_type operator()(const base_node &P,
				   dal::bit_vector &bv) const {
      scalar_type d = xon - bgeot::vect_sp(P,n);
      bv[id] = (dal::abs(d) < SEPS);
      return d;
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      id = list.size(); list.push_back(this);
    }
    virtual base_small_vector grad(const base_node &) const {
      return -1.*n;
    }
  };

  class mesher_tube : public mesher_signed_distance {
    base_node x0; base_node n; scalar_type R;
  public:
      mesher_tube(base_node x0_, base_node n_, scalar_type R_)
	: x0(x0_), n(n_), R(R_)
    { n /= gmm::vect_norm2(n); }
    bool bounding_box(base_node &, base_node &) const
    { return false; }
    virtual scalar_type operator()(const base_node &P) const {
      base_node v(P); v -= x0;
      gmm::add(gmm::scaled(n, -gmm::vect_sp(v, n)), v);
      return gmm::vect_norm2(v) - R;
    }
    virtual scalar_type operator()(const base_node &P,
				   dal::bit_vector &bv) const {
      scalar_type d = (*this)(P);
      bv[id] = (dal::abs(d) < SEPS);
      return d;
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      id = list.size(); list.push_back(this);
    }
    virtual base_small_vector grad(const base_node &P) const {
      base_node v(P); v -= x0;
      scalar_type norm(0);
      for ( ; norm == 0.; gmm::fill_random(v)) {
	gmm::add(gmm::scaled(n, -gmm::vect_sp(v, n)), v);
	norm = gmm::vect_norm2(v);
      }
      v /= norm;
      return v;
    }
  };


  class mesher_ball : public mesher_signed_distance {
    base_node x0; scalar_type R;
  public:
    mesher_ball(base_node x0_, scalar_type R_) : x0(x0_), R(R_) {}
    bool bounding_box(base_node &bmin, base_node &bmax) const { 
      bmin = bmax = x0; 
      for (size_type i=0; i < x0.size(); ++i) { bmin[i] -= R; bmax[i] += R; }
      return true;
    }
    virtual scalar_type operator()(const base_node &P,
				   dal::bit_vector &bv) const {
      scalar_type d = bgeot::vect_dist2(P,x0)-R;
      bv[id] = (dal::abs(d) < SEPS);
      return d;
    }
    virtual scalar_type operator()(const base_node &P) const
    { return bgeot::vect_dist2(P,x0)-R; }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      id = list.size(); list.push_back(this);
    }
    virtual base_small_vector grad(const base_node &P) const {
      base_small_vector g(P - x0);
      scalar_type d= gmm::vect_norm2(g);
      if (d != scalar_type(0)) { g /= d; return g; }
      else { gmm::fill_random(g); g /= gmm::vect_norm2(g); return g; }
    }
  };

  class mesher_rectangle : public mesher_signed_distance {
    // ajouter une rotation rigide et translation ..
    base_node rmin, rmax;
    std::vector<mesher_half_space> hfs;
  public:
    mesher_rectangle(base_node rmin_, base_node rmax_)
      : rmin(rmin_), rmax(rmax_) {
      base_node n(rmin_.size());
      for (unsigned k = 0; k < rmin.size(); ++k) {
	n[k] = 1.0;
	hfs.push_back(mesher_half_space(rmin, n));
	n[k] = -1.0;
	hfs.push_back(mesher_half_space(rmax, n));
	n[k] = 0.0;
      }
    }
    bool bounding_box(base_node &bmin, base_node &bmax) const {
      bmin = rmin; bmax = rmax;
      return true;
    }
    virtual scalar_type operator()(const base_node &P) const {
      size_type N = rmin.size();
      scalar_type d = rmin[0] - P[0];
      for (size_type i=0; i < N; ++i) {
	d = std::max(d, rmin[i] - P[i]);
	d = std::max(d, P[i] - rmax[i]);
      }
      return d;
    }

    virtual scalar_type operator()(const base_node &P, dal::bit_vector &bv)
      const {
      scalar_type d = this->operator()(P);
      if (gmm::abs(d) < SEPS)
	for (int k = 0; k < 2*rmin.size(); ++k) hfs[k](P, bv);
      return d;
    }
    virtual base_small_vector grad(const base_node &P) const {
      unsigned i = 0; scalar_type di = hfs[i](P);
      for (int k = 1; k < 2*rmin.size(); ++k) {
	scalar_type dk = hfs[k](P);
	if (dk > di) { i = k; di = dk; }
      }
      return hfs[i].grad(P);
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      for (int k = 0; k < 2*rmin.size(); ++k)
	hfs[k].register_constraints(list);
    }
  };

  class mesher_union : public mesher_signed_distance {
    const mesher_signed_distance &a, &b;
  public:
    mesher_union(const mesher_signed_distance& a_,
		 const mesher_signed_distance &b_) : a(a_), b(b_) {}
    bool bounding_box(base_node &bmin, base_node &bmax) const {
      base_node bmin2(bmin.size()), bmax2(bmin.size());
      bool ba = a.bounding_box(bmin, bmax);
      bool bb = b.bounding_box(bmin2, bmax2);
      if (!ba || !bb) return false;
      for (unsigned i=0; i < bmin.size(); ++i) { 
        bmin[i] = std::min(bmin[i],bmin2[i]);
	bmax[i] = std::max(bmax[i],bmax2[i]);
      }
      return true;
    }
    scalar_type operator()(const base_node &P, dal::bit_vector &bv) const {
      scalar_type da = a(P), db = b(P);
      if (da > -SEPS && db > -SEPS) {
	if (da < SEPS) a(P, bv);
	if (db < SEPS) b(P, bv);
      }
      return std::min(da,db);
    }
    virtual scalar_type operator()(const base_node &P) const
    { return std::min(a(P),b(P)); }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      a.register_constraints(list); b.register_constraints(list);
    }
    virtual base_small_vector grad(const base_node &P) const {
      if (a(P) < b(P)) return a.grad(P);
      else return b.grad(P);
    }
  };

  class mesher_intersection : public mesher_signed_distance {
    const mesher_signed_distance &a, &b;
  public:
    mesher_intersection(const mesher_signed_distance &a_,
			const mesher_signed_distance &b_) : a(a_), b(b_) {}
    bool bounding_box(base_node &bmin, base_node &bmax) const {
      bool ba = a.bounding_box(bmin, bmax);
      base_node bbmin(bmin.size()), bbmax(bmin.size());
      bool bb = b.bounding_box(bbmin, bbmax);
      if (!ba) { bmin = bbmin; bmax = bbmax; }
      if (ba && bb)
	for (unsigned k = 0; k < bmin.size(); ++k) {
	  bmin[k] = std::max(bmin[k], bbmin[k]);
	  bmax[k] = std::max(std::min(bmax[k], bbmax[k]), bmin[k]);
	}
      return (ba || bb);
    }
    scalar_type operator()(const base_node &P, dal::bit_vector &bv) const {
      scalar_type da = a(P), db = b(P);
      if (da < SEPS && db < SEPS) {
	if (da > -SEPS) a(P, bv);
	if (db > -SEPS) b(P, bv);
      }
      return std::max(da,db);
    }
    scalar_type operator()(const base_node &P) const
    { return std::max(a(P),b(P)); }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      a.register_constraints(list); b.register_constraints(list);
    }
    virtual base_small_vector grad(const base_node &P) const {
      scalar_type da = a(P), db = b(P);
      if (da > db) return a.grad(P);
      else return b.grad(P);
    }
  };

  class mesher_setminus : public mesher_signed_distance {
    const mesher_signed_distance &a, &b;
  public:
    mesher_setminus(const mesher_signed_distance& a_,
		    const mesher_signed_distance &b_) : a(a_), b(b_) {}
    bool bounding_box(base_node &bmin, base_node &bmax) const
    { return a.bounding_box(bmin,bmax); }
    scalar_type operator()(const base_node &P, dal::bit_vector &bv) const {
      scalar_type da = a(P), db = -b(P);
      if (da < SEPS && db < SEPS) {
	if (da > -SEPS) a(P, bv);
	if (db > -SEPS) b(P, bv);
      }
      return std::max(da, db);
    }
    scalar_type operator()(const base_node &P) const
    { return std::max(a(P),-b(P)); }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      a.register_constraints(list); b.register_constraints(list);
    }
    virtual base_small_vector grad(const base_node &P) const {
      scalar_type da = a(P), db = -b(P);
      if (da > db) return a.grad(P);
      else { base_small_vector G(b.grad(P)); G *= -1.; return G; }
    }
  };
  
  class mesher_cylinder : public mesher_signed_distance {
    base_node x0; base_small_vector n;
    scalar_type L, R;
    mesher_tube t;
    mesher_half_space p1, p2;
    mesher_intersection i1, i2;
  public:
    mesher_cylinder(const base_node &center, const base_small_vector &no,
		    scalar_type L_, scalar_type R_)
      : x0(center), n(no/gmm::vect_norm2(no)), L(L_), R(R_), t(x0, n, R_),
	p1(x0, n), p2(x0+n*L, -1.0 * n), i1(p1, p2), i2(i1, t) {}
    bool bounding_box(base_node &bmin, base_node &bmax) const {
      base_node x1(x0+n*L);
      bmin = bmax = base_node(3);
      for (unsigned i = 0; i < 3; ++i) {
	bmin[i] = std::min(x0[i], x1[i]) - R;
	bmax[i] = std::max(x0[i], x1[i]) + R;
      }
      return true;
    }
    virtual scalar_type operator()(const base_node &P) const { return i2(P); }
    virtual scalar_type operator()(const base_node &P,
				   dal::bit_vector& bv) const
    { return i2(P, bv); }
    virtual base_small_vector grad(const base_node &P) const
    { return i2.grad(P); }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const
    { i2.register_constraints(list); }
  };


class mesher_ellipse : public mesher_signed_distance { // TODO
  base_node x0; base_small_vector n, t;
  scalar_type r, R, a;
  public:
    mesher_ellipse(const base_node &center, const base_small_vector &no,
		    scalar_type r_, scalar_type R_)
      : x0(center), n(no/gmm::vect_norm2(no)), r(r_), R(R_) {
      t[0] = -n[1]; t[1] = n[0];
      if (R < r) { std::swap(r, R); std::swap(n, t); }
      a = sqrt(R*R - r*r);
    }
    bool bounding_box(base_node &bmin, base_node &bmax) const {
      bmin = bmax = x0;
      for (unsigned i = 0; i < 2; ++i) { bmin[i] -= R; bmax[i] += R; }
      return true;
    }
    virtual scalar_type operator()(const base_node &P) const { 
      base_small_vector v(P); v -= x0;
      scalar_type vt = gmm::vect_sp(v, t);
      vt = std::max(-a, std::min(a, vt));
      base_node x1 = x0 + t*vt;
      base_small_vector v1(P); v1 -= x1;
      scalar_type v1n = gmm::vect_sp(v1, n), v1t = gmm::vect_sp(v1, t);
      scalar_type x1n = gmm::vect_sp(x1, n), x1t = gmm::vect_sp(x1, t);
      scalar_type ea = v1n*v1n / (r*r) + v1t * v1t / (R*R);
      scalar_type eb = 2. * (x1n*v1n / (r*r) + x1t*v1t / (R*R));
      scalar_type ec = x1n*x1n / (r*r) + x1t * x1t / (R*R);

      scalar_type delta = eb*eb - 4 * ea * ec;
      assert(delta >= 0);
      scalar_type lambda = (-eb + sqrt(delta)) / (2. * ea);
      base_node x2 = lambda*P + (1-lambda)*x1;
      return (1.-lambda)*gmm::vect_norm2(v1);
    }
    virtual scalar_type operator()(const base_node &P,
				   dal::bit_vector& bv) const {
      scalar_type d = this->operator()(P);
      bv[id] = (dal::abs(d) < SEPS);
      return d;
    }
    virtual base_small_vector grad(const base_node &P) const
  { assert(0); }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const
    { id = list.size(); list.push_back(this); }
  };


  class mesher_torus : public mesher_signed_distance { // to be done
    // ajouter une rotation rigide et translation ..
    scalar_type R, r;
  public:
    mesher_torus(scalar_type RR, scalar_type rr) : R(RR), r(rr) {}
    bool bounding_box(base_node &bmin, base_node &bmax) const {
      bmin = base_node(3); bmax = base_node(3);
      bmin[0] = bmin[1] = -R-r; bmin[2] = -r;
      bmax[0] = bmax[1] = +R+r; bmax[2] = +r;
      return true;
    }
    virtual scalar_type operator()(const base_node &P) const {
      scalar_type x = P[0], y = P[1], z = P[2], c = sqrt(x*x + y*y);
      return (c == 0.) ? R - r : sqrt(gmm::sqr(c-R) + z*z) - r;
    }
    virtual scalar_type operator()(const base_node &P, dal::bit_vector&bv)
      const {
      scalar_type d = this->operator()(P);
      bv[id] = (dal::abs(d) < SEPS);
      return d;
    }
    virtual base_small_vector grad(const base_node &P) const {
      base_node G(3);
      scalar_type x = P[0], y = P[1], z = P[2], c = sqrt(x*x + y*y);
      if (c == 0.) { 
	gmm::fill_random(G); G[2] = 0.0; G /= gmm::vect_norm2(G);
      }
      else {
	scalar_type w = 1. - R / c, e = sqrt(gmm::sqr(c-R) + z*z);
	if (e == 0.) {
	  gmm::fill_random(G); G[0] = x; G[1] = y; G /= gmm::vect_norm2(G);
	}
	else {
	  G[0] = x * w / e; G[1] = y * w / e; G[2] = z / e;
	}
      }
      return G;
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>&list) const
    { id = list.size(); list.push_back(this); }
  };


  

  // interface with qhull
  void delaunay(const std::vector<base_node> &pts,
		gmm::dense_matrix<size_type>& simplexes);

  
  // mesher
  void build_mesh(getfem_mesh &m, const mesher_signed_distance& dist_,
		  scalar_type h0, const std::vector<base_node> &fixed_points
		  = std::vector<base_node>(), size_type K = 1, int noise = 1,
		  size_type iter_max = 1000, scalar_type dist_point_hull = 5,
		  scalar_type boundary_threshold_flatness = 0.11);


}

