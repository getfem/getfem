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


#ifdef HAVE_QHULL_QHULL_H
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
    virtual scalar_type operator()(const base_node &P, dal::bit_vector &bv) const = 0;
    virtual base_small_vector grad(const base_node &P) const = 0;
    virtual void register_constraints(std::vector<const mesher_signed_distance*>& list) const = 0;
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
    virtual void register_constraints(std::vector<const mesher_signed_distance*>& list) const {
      id = list.size(); list.push_back(this);
    }
    virtual base_small_vector grad(const base_node &) const {
      return -1.*n;
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
    virtual scalar_type operator()(const base_node &P, dal::bit_vector &bv) const {
      scalar_type d = bgeot::vect_dist2(P,x0)-R;
      bv[id] = (dal::abs(d) < SEPS);
      return d;
    }
    virtual scalar_type operator()(const base_node &P) const
    { return bgeot::vect_dist2(P,x0)-R; }
    virtual void register_constraints(std::vector<const mesher_signed_distance*>& list) const {
      id = list.size(); list.push_back(this);
    }
    virtual base_small_vector grad(const base_node &P) const {
      base_small_vector g(P - x0);
      scalar_type d= gmm::vect_norm2(g);
      if (d != scalar_type(0)) return g / d; else return P;
      return g/gmm::vect_norm2(g);
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
    mesher_union(const mesher_signed_distance& a_, const mesher_signed_distance &b_) : a(a_), b(b_) {}
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
    mesher_intersection(const mesher_signed_distance& a_,
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
  
  class mesher_cylinder : public mesher_signed_distance { // to be done
  public:
    mesher_cylinder() {}
    bool bounding_box(base_node &bmin, base_node &bmax) const {
      bmin = base_node(3); bmax = base_node(3); 
      for (size_type i=0; i < 3; ++i) { bmin[i] = -1.; bmax[i] = 1.; }
      return true;
    }
    virtual scalar_type operator()(const base_node &P) const {
      scalar_type d1 = sqrt(dal::sqr(P[0])+dal::sqr(P[1])) - 1;
      scalar_type d2 = -1 - P[2];
      scalar_type d3 = P[2] - 1;
      scalar_type d = std::max(std::max(d1,d2),d3);
      if (d1 > 0 && d2 > 0) d = sqrt(dal::sqr(d1)+dal::sqr(d2));
      else if (d1 > 0 && d3 > 0) d = sqrt(dal::sqr(d1)+dal::sqr(d3));
      return d;
    }
    virtual scalar_type operator()(const base_node &P, dal::bit_vector& )
      const {
      scalar_type d = this->operator()(P);
      // to be done
      return d;
    }
    virtual base_small_vector grad(const base_node &) const { assert(0); }
    virtual void register_constraints(std::vector<const mesher_signed_distance*>&) const { assert(0); }
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
      scalar_type a = P[0];
      scalar_type b = P[1];
      scalar_type c = sqrt(a*a + b*b);
      if (c == 0.) return R - r;
      scalar_type d = 1. - R / c;
      return sqrt(gmm::sqr(d*P[0]) + gmm::sqr(d*P[1]) + gmm::sqr(P[2])) - r;
    }
    virtual scalar_type operator()(const base_node &P, dal::bit_vector&bv)
      const {
      scalar_type d = this->operator()(P);
      bv[id] = (dal::abs(d) < SEPS);
      return d;
    }
    virtual base_small_vector grad(const base_node &P) const {
      base_node G(3);
      scalar_type a = P[0];
      scalar_type b = P[1];
      scalar_type c = sqrt(a*a + b*b);
      if (c == 0.) return G;
      scalar_type d = 1. - R / c;
      scalar_type e = sqrt(gmm::sqr(d*P[0])+gmm::sqr(d*P[1]) + gmm::sqr(P[2]));
      if (e == 0) return G;

      G[0] = (gmm::sqr(d) + d * R*gmm::sqr(P[0]) / (c*c*c)) * P[0]/ e;
      G[1] = (gmm::sqr(d) + d * R*gmm::sqr(P[1]) / (c*c*c)) * P[1]/ e;
      G[2] = P[2] / e;
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
		  = std::vector<base_node>(), size_type K = 1);


}

#endif
