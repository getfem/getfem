/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_.                       */
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

#define SEPS 1e-5

  class mesher_signed_distance : public mesher_virtual_function {
  protected:
    mutable size_type id;
  public:
    mesher_signed_distance() : id(size_type(-1)) {}
    virtual void bounding_box(base_node &bmin, base_node &bmax) const = 0;
    virtual scalar_type operator()(const base_node &P, dal::bit_vector &bv) const = 0;
    virtual base_small_vector grad(const base_node &P) const = 0;
    virtual void register_constraints(std::vector<const mesher_signed_distance*>& list) const = 0;
    scalar_type operator()(const base_node &P) const {
      dal::bit_vector bv;
      return (*this)(P,bv);
    }
  };

  class mesher_half_space : public mesher_signed_distance {
    base_node x0; base_node n;
  public:
    mesher_half_space(base_node x0_, base_node n_) : x0(x0_), n(n_) { n /= gmm::vect_norm2(n); }
    void bounding_box(base_node &, base_node &) const { 
      /* to be done ... */
    }
    virtual scalar_type operator()(const base_node &P, dal::bit_vector &bv) const {
      scalar_type d = bgeot::vect_sp(x0-P,n);
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


  class mesher_sphere : public mesher_signed_distance {
    base_node x0; scalar_type R;
  public:
    mesher_sphere(base_node x0_, scalar_type R_) : x0(x0_), R(R_) {}
    void bounding_box(base_node &bmin, base_node &bmax) const { 
      bmin = bmax = x0; 
      for (size_type i=0; i < x0.size(); ++i) { bmin[i] -= R; bmax[i] += R; }
    }
    virtual scalar_type operator()(const base_node &P, dal::bit_vector &bv) const {
      scalar_type d = bgeot::vect_dist2(P,x0)-R;
      bv[id] = (dal::abs(d) < SEPS);
      return d;
    }
    virtual void register_constraints(std::vector<const mesher_signed_distance*>& list) const {
      id = list.size(); list.push_back(this);
    }
    virtual base_small_vector grad(const base_node &P) const {
      base_small_vector g(P - x0);
      return g/gmm::vect_norm2(g);
    }
  };

  class mesher_rectangle : public mesher_signed_distance {
    base_node rmin, rmax;
  public:
    mesher_rectangle(base_node rmin_, base_node rmax_) : rmin(rmin_), rmax(rmax_) {}
    void bounding_box(base_node &bmin, base_node &bmax) const {
      bmin = rmin; bmax = rmax; 
    }
    virtual scalar_type operator()(const base_node &P, dal::bit_vector &) const {
      size_type N = rmin.size();
      scalar_type d = -1000000000000.;
      for (size_type i=0; i < N; ++i) {
	d = std::max(d, rmin[i] - P[i]);
	d = std::max(d, P[i] - rmax[i]);
      }
      return d;
      if (d > 0) { /* handle distance to the vertices */
	for (size_type i=1; i < pow(3,N); ++i) {
	  base_node Q(rmin);
	  for (size_type m=i,k=0; m; ++k,m/=3) { 
	    if (m % 3 == 0) Q[k] = P[k];
	    else if (m % 3 == 1) Q[k] = rmax[k]; 
	  }
	  d = std::min(d, bgeot::vect_dist2(P,Q));
	}
      }
      return d;
    }
    virtual base_small_vector grad(const base_node &) const { assert(0); }
    virtual void register_constraints(std::vector<const mesher_signed_distance*>&) const { cout << "to be done\n"; assert(0); }
  };

  class mesher_union : public mesher_signed_distance {
    const mesher_signed_distance &a, &b;
  public:
    mesher_union(const mesher_signed_distance& a_, const mesher_signed_distance &b_) : a(a_), b(b_) {}
    void bounding_box(base_node &bmin, base_node &bmax) const {
      base_node bmin2, bmax2;
      a.bounding_box(bmin,bmax);b.bounding_box(bmin2,bmax2);
      for (size_type i=0; i < bmin.size(); ++i) { 
        bmin[i] = std::min(bmin[i],bmin2[i]);bmax[i] = std::max(bmax[i],bmax2[i]);
      }
    }
    scalar_type operator()(const base_node &P, dal::bit_vector &bv) const {
      return std::min(a(P,bv),b(P,bv));
    }
    virtual void register_constraints(std::vector<const mesher_signed_distance*>& list) const {
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
    mesher_intersection(const mesher_signed_distance& a_, const mesher_signed_distance &b_) : a(a_), b(b_) {}
    void bounding_box(base_node &bmin, base_node &bmax) const {
      a.bounding_box(bmin,bmax);b.bounding_box(bmin,bmax);
      /* TODO */
    }
    scalar_type operator()(const base_node &P, dal::bit_vector &bv) const {
      return std::max(a(P,bv),b(P,bv));
    }
    virtual void register_constraints(std::vector<const mesher_signed_distance*>& list) const {
      a.register_constraints(list); b.register_constraints(list);
    }
    virtual base_small_vector grad(const base_node &P) const {
      scalar_type da = a(P), db = b(P);
      if (da > db) return a.grad(P);
      else return b.grad(P);
    }
  };
  
  class mesher_cylinder : public mesher_signed_distance {
  public:
    mesher_cylinder() {}
    void bounding_box(base_node &bmin, base_node &bmax) const {
      bmin = base_node(3); bmax = base_node(3); 
      for (size_type i=0; i < 3; ++i) { bmin[i] = -1.; bmax[i] = 1.; }
    }
    virtual scalar_type operator()(const base_node &P, dal::bit_vector& ) const {
      scalar_type d1 = sqrt(dal::sqr(P[0])+dal::sqr(P[1])) - 1;
      scalar_type d2 = -1 - P[2];
      scalar_type d3 = P[2] - 1;
      scalar_type d = std::max(std::max(d1,d2),d3);
      if (d1 > 0 && d2 > 0) d = sqrt(dal::sqr(d1)+dal::sqr(d2));
      else if (d1 > 0 && d3 > 0) d = sqrt(dal::sqr(d1)+dal::sqr(d3));
      return d;
    }
    virtual base_small_vector grad(const base_node &) const { assert(0); }
    virtual void register_constraints(std::vector<const mesher_signed_distance*>&) const { assert(0); }
  };


}
