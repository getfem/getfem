/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2004-2015 Julien Pommier, Yves Renard

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/**@file getfem_mesher.h
   @author  Julien Pommier <Julien.Pommier@insa-toulouse.fr>, Yves Renard <Yves.Renard@insa-lyon.fr>
   @date May 1, 2004.
   @brief An experimental mesher.
*/

#ifndef GETFEM_MESHER_H__
#define GETFEM_MESHER_H__


#include "getfem_mesh.h"
#include "gmm/gmm_condition_number.h"
#include "bgeot_comma_init.h"
#include "gmm/gmm_solver_bfgs.h"
#include "getfem_export.h"
#include "bgeot_kdtree.h"
#include <typeinfo>

namespace getfem {

  class mesher_virtual_function : public dal::static_stored_object {
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

#define SEPS 1e-8

  class mesher_signed_distance : public mesher_virtual_function {
  protected:
    mutable size_type id;
  public:
    mesher_signed_distance() : id(size_type(-1)) {}
    virtual ~mesher_signed_distance() {}
    virtual bool bounding_box(base_node &bmin, base_node &bmax) const = 0;
    virtual scalar_type operator()(const base_node &P,
				   dal::bit_vector &bv) const = 0;
    virtual scalar_type grad(const base_node &P,
			     base_small_vector &G) const = 0;
    virtual void hess(const base_node &P, base_matrix &H) const = 0;
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const=0;
    virtual scalar_type operator()(const base_node &P) const  = 0;
  };

  typedef std::shared_ptr<const mesher_signed_distance> pmesher_signed_distance;

  class mesher_half_space : public mesher_signed_distance {
    base_node x0; base_small_vector n; scalar_type xon;
  public:
    mesher_half_space(const base_node &x0_, const base_small_vector &n_)
      : x0(x0_), n(n_)
    { n /= gmm::vect_norm2(n); xon = gmm::vect_sp(x0, n); }
    bool bounding_box(base_node &, base_node &) const
    { return false; }
    virtual scalar_type operator()(const base_node &P) const
    { return xon - gmm::vect_sp(P,n); }
    virtual scalar_type operator()(const base_node &P,
				   dal::bit_vector &bv) const {
      scalar_type d = xon - gmm::vect_sp(P,n);
      bv[id] = (gmm::abs(d) < SEPS);
      return d;
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      id = list.size(); list.push_back(this);
    }
    scalar_type grad(const base_node &P, base_small_vector &G) const {
      G = n; G *= scalar_type(-1); 
      return xon - gmm::vect_sp(P,n);
    }
    void hess(const base_node &P, base_matrix &H) const {
      gmm::resize(H, P.size(), P.size()); gmm::clear(H);
    }

  };

  inline pmesher_signed_distance  new_mesher_half_space
  (const base_node &x0, const base_small_vector &n)
  { return std::make_shared<mesher_half_space>(x0, n); }

  class mesher_level_set : public mesher_signed_distance {
    bgeot::base_poly base;
    mutable std::vector<base_poly> gradient;
    mutable std::vector<base_poly> hessian;
    const fem<base_poly> *pf;
    mutable int initialized;
    scalar_type shift_ls;     // for the computation of a gap on a level_set.
  public:
    bool is_initialized(void) const { return initialized; }
    mesher_level_set() : initialized(0) {}
    template <typename VECT>
    mesher_level_set(pfem pf_, const VECT &coeff_,
		     scalar_type shift_ls_ = scalar_type(0)) {
      init_base(pf_, coeff_);
      set_shift(shift_ls_);
    }
    void set_shift(scalar_type shift_ls_) {
      shift_ls = shift_ls_; 
      if (shift_ls != scalar_type(0)) {
	base_node P(pf->dim()); base_small_vector G(pf->dim());
	grad(P, G);
	shift_ls *= gmm::vect_norm2(G);
      }
    }
    template <typename VECT> void init_base(pfem pf_, const VECT &coeff_);
    void init_grad(void) const;
    void init_hess(void) const;

    bool bounding_box(base_node &, base_node &) const
    { return false; }
    virtual scalar_type operator()(const base_node &P) const
    {  return bgeot::to_scalar(base.eval(P.begin())) + shift_ls; }
    virtual scalar_type operator()(const base_node &P,
				   dal::bit_vector &bv) const
    { scalar_type d = (*this)(P); bv[id] = (gmm::abs(d) < SEPS); return d; }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      id = list.size(); list.push_back(this);
    }
    scalar_type grad(const base_node &P, base_small_vector &G) const;
    void hess(const base_node &P, base_matrix &H) const;
  };

  template <typename VECT> 
  void mesher_level_set::init_base(pfem pf_, const VECT &coeff_) {
    std::vector<scalar_type> coeff(coeff_.begin(), coeff_.end());
    GMM_ASSERT1(gmm::vect_norm2(coeff) != 0, "level is zero!");
    pf = dynamic_cast<const fem<base_poly>* > (pf_.get());
    GMM_ASSERT1(pf, "PK fem are required for level set (got "
		<< typeid(pf_).name() << ")");
    base = base_poly(pf->base()[0].dim(), pf->base()[0].degree());
    for (unsigned i=0; i < pf->nb_base(0); ++i) {
      base += pf->base()[i] * coeff[i];
    }
    initialized = 0; 
  }

  template <typename VECT> 
  inline pmesher_signed_distance new_mesher_level_set
  (pfem pf, const VECT &coeff, scalar_type shift_ls = scalar_type(0))
  { return std::make_shared<mesher_level_set>(pf, coeff, shift_ls); }


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
      scalar_type d = gmm::vect_dist2(P,x0)-R;
      bv[id] = (gmm::abs(d) < SEPS);
      return d;
    }
    virtual scalar_type operator()(const base_node &P) const
    { return gmm::vect_dist2(P,x0)-R; }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      id = list.size(); list.push_back(this);
    }
    scalar_type grad(const base_node &P, base_small_vector &G) const {
      G = P; G -= x0;
      scalar_type e= gmm::vect_norm2(G), d = e - R;
      while (e == scalar_type(0))
	{ gmm::fill_random(G); e = gmm::vect_norm2(G); }
      G /= e;
      return d;
    }
    void hess(const base_node &, base_matrix &) const {
      GMM_ASSERT1(false, "Sorry, to be done");
    }
  };

  inline pmesher_signed_distance new_mesher_ball(base_node x0, scalar_type R)
  { return std::make_shared<mesher_ball>(x0, R); }

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
    scalar_type grad(const base_node &P, base_small_vector &G) const {
      unsigned i = 0; scalar_type di = hfs[i](P);
      for (int k = 1; k < 2*rmin.size(); ++k) {
	scalar_type dk = hfs[k](P);
	if (dk > di) { i = k; di = dk; }
      }
      return hfs[i].grad(P, G);
    }
    void hess(const base_node &P, base_matrix &H) const {
      gmm::resize(H, P.size(), P.size()); gmm::clear(H);
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      for (int k = 0; k < 2*rmin.size(); ++k)
	hfs[k].register_constraints(list);
    }
  };

  inline pmesher_signed_distance new_mesher_rectangle(base_node rmin,
						      base_node rmax)
  { return std::make_shared<mesher_rectangle>(rmin, rmax); }

  class mesher_simplex_ref : public mesher_signed_distance {
    // ajouter une rotation rigide, homothetie,  et translation ..
    std::vector<mesher_half_space> hfs;
    unsigned N;
    base_node org;
  public:
    mesher_simplex_ref(unsigned N_) : N(N_) {
      base_node no(N);
      org = no;
      for (unsigned k = 0; k < N; ++k) {
	no[k] = scalar_type(1);
	hfs.push_back(mesher_half_space(org, no));
	no[k] = scalar_type(0);
      }
      std::fill(org.begin(), org.end(), scalar_type(1)/scalar_type(N));
      no = -org;
      hfs.push_back(mesher_half_space(org, no));
    }
    bool bounding_box(base_node &bmin, base_node &bmax) const {
      bmin.resize(N); bmax.resize(N);
      std::fill(bmin.begin(), bmin.end(), scalar_type(0));
      std::fill(bmax.begin(), bmax.end(), scalar_type(1));
      return true;
    }
    virtual scalar_type operator()(const base_node &P) const {
      scalar_type d = - P[0];
      for (size_type i=1; i < N; ++i) d = std::max(d, - P[i]);
      d = std::max(d, gmm::vect_sp(P - org, org) / gmm::vect_norm2(org));
      return d;
    }

    virtual scalar_type operator()(const base_node &P, dal::bit_vector &bv)
      const {
      scalar_type d = this->operator()(P);
      if (gmm::abs(d) < SEPS) for (unsigned k = 0; k < N+1; ++k) hfs[k](P, bv);
      return d;
    }
    scalar_type grad(const base_node &P, base_small_vector &G) const {
      unsigned i = 0; scalar_type di = hfs[i](P);
      for (unsigned k = 1; k < N+1; ++k) {
	scalar_type dk = hfs[k](P);
	if (dk > di) { i = k; di = dk; }
      }
      return hfs[i].grad(P, G);
    }
    void hess(const base_node &P, base_matrix &H) const {
      gmm::resize(H, P.size(), P.size()); gmm::clear(H);
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const
    { for (unsigned k = 0; k < N+1; ++k) hfs[k].register_constraints(list); }
  };

  inline pmesher_signed_distance new_mesher_simplex_ref(unsigned N)
  { return std::make_shared<mesher_simplex_ref>(N); }


  class mesher_prism_ref : public mesher_signed_distance {
    // ajouter une rotation rigide, homothetie,  et translation ..
    std::vector<mesher_half_space> hfs;
    unsigned N;
    base_node org;
  public:
    mesher_prism_ref(unsigned N_) : N(N_) {
      base_node no(N);
      org = no;
      for (unsigned k = 0; k < N; ++k) {
	no[k] = scalar_type(1);
	hfs.push_back(mesher_half_space(org, no));
	no[k] = scalar_type(0);
      }
      no[N-1] = -scalar_type(1);
      org[N-1] = scalar_type(1);
      hfs.push_back(mesher_half_space(org, no));
      std::fill(org.begin(), org.end(), scalar_type(1)/scalar_type(N));
      org[N-1] = scalar_type(0);
      no = -org;
      hfs.push_back(mesher_half_space(org, no));
      
    }
    bool bounding_box(base_node &bmin, base_node &bmax) const {
      bmin.resize(N); bmax.resize(N);
      std::fill(bmin.begin(), bmin.end(), scalar_type(0));
      std::fill(bmax.begin(), bmax.end(), scalar_type(1));
      return true;
    }
    virtual scalar_type operator()(const base_node &P) const {
      scalar_type d = - P[0];
      for (size_type i=1; i < N; ++i) d = std::max(d, - P[i]);
      d = std::max(d, P[N-1] - scalar_type(1));
      d = std::max(d, gmm::vect_sp(P - org, org) / gmm::vect_norm2(org));
      return d;
    }

    virtual scalar_type operator()(const base_node &P, dal::bit_vector &bv)
      const {
      scalar_type d = this->operator()(P);
      if (gmm::abs(d) < SEPS) for (unsigned k = 0; k < N+2; ++k) hfs[k](P, bv);
      return d;
    }
    scalar_type grad(const base_node &P, base_small_vector &G) const {
      unsigned i = 0; scalar_type di = hfs[i](P);
      for (unsigned k = 1; k < N+2; ++k) {
	scalar_type dk = hfs[k](P);
	if (dk > di) { i = k; di = dk; }
      }
      return hfs[i].grad(P, G);
    }
    void hess(const base_node &P, base_matrix &H) const {
      gmm::resize(H, P.size(), P.size()); gmm::clear(H);
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const
    { for (unsigned k = 0; k < N+2; ++k) hfs[k].register_constraints(list); }
  };

  inline pmesher_signed_distance new_mesher_prism_ref(unsigned N)
  { return std::make_shared<mesher_prism_ref>(N); }



  // Rem : It would be better for the convergence of Newton's methods to
  //       take somthing more smooth than min. for instance tthe square root
  //       of positive parts or negative parts depending on the position of
  //       point (in or out the domain).

  class mesher_union : public mesher_signed_distance {
    std::vector<pmesher_signed_distance> dists;
    mutable std::vector<scalar_type> vd;
    mutable bool isin;
    bool with_min;
  public:
    mesher_union(const std::vector<pmesher_signed_distance>
			&dists_) : dists(dists_) 
    { vd.resize(dists.size()); with_min = true; }

    mesher_union
    (const pmesher_signed_distance &a,
     const pmesher_signed_distance &b,
     const pmesher_signed_distance &c = pmesher_signed_distance(),
     const pmesher_signed_distance &d = pmesher_signed_distance(),
     const pmesher_signed_distance &e = pmesher_signed_distance(),
     const pmesher_signed_distance &f = pmesher_signed_distance(),
     const pmesher_signed_distance &g = pmesher_signed_distance(),
     const pmesher_signed_distance &h = pmesher_signed_distance(),
     const pmesher_signed_distance &i = pmesher_signed_distance(),
     const pmesher_signed_distance &j = pmesher_signed_distance(),
     const pmesher_signed_distance &k = pmesher_signed_distance(),
     const pmesher_signed_distance &l = pmesher_signed_distance(),
     const pmesher_signed_distance &m = pmesher_signed_distance(),
     const pmesher_signed_distance &n = pmesher_signed_distance(),
     const pmesher_signed_distance &o = pmesher_signed_distance(),
     const pmesher_signed_distance &p = pmesher_signed_distance(),
     const pmesher_signed_distance &q = pmesher_signed_distance(),
     const pmesher_signed_distance &r = pmesher_signed_distance(),
     const pmesher_signed_distance &s = pmesher_signed_distance(),
     const pmesher_signed_distance &t = pmesher_signed_distance()) {
      dists.push_back(a); dists.push_back(b);
      with_min = true;
      if (c) dists.push_back(c); if (d) dists.push_back(d);
      if (e) dists.push_back(e); if (f) dists.push_back(f);
      if (g) dists.push_back(g); if (h) dists.push_back(h);
      if (i) dists.push_back(i); if (j) dists.push_back(j);
      if (k) dists.push_back(k); if (l) dists.push_back(l);
      if (m) dists.push_back(m); if (n) dists.push_back(n);
      if (o) dists.push_back(o); if (p) dists.push_back(p);
      if (q) dists.push_back(q); if (r) dists.push_back(r);
      if (s) dists.push_back(s); if (t) dists.push_back(t);
      vd.resize(dists.size());
    }
    
    bool bounding_box(base_node &bmin, base_node &bmax) const {
      base_node bmin2, bmax2;
      bool b = dists[0]->bounding_box(bmin, bmax);
      if (!b) return false;
      for (size_type k = 1; k < dists.size(); ++k) {
	b = dists[k]->bounding_box(bmin2, bmax2);
	if (!b) return false;
	for (unsigned i=0; i < bmin.size(); ++i) { 
	  bmin[i] = std::min(bmin[i],bmin2[i]);
	  bmax[i] = std::max(bmax[i],bmax2[i]);
	}
      }
      return true;
    }
    virtual scalar_type operator()(const base_node &P) const {
      scalar_type d, f(0), g(1);
      if (with_min) {
	d = (*(dists[0]))(P);
	for (size_type k = 1; k < dists.size(); ++k)
	  d = std::min(d, (*(dists[k]))(P));
      }
      else { // essai raté ...
	isin = false;
	for (size_type k = 0; k < dists.size(); ++k) {
	  vd[k] = (*(dists[k]))(P);
	  if (vd[k] <= scalar_type(0)) isin = true;
	  f += gmm::sqr(gmm::neg(vd[k]));
	  g *= gmm::pos(vd[k]);
	}
	d = isin ? -gmm::sqrt(f)
	  : pow(g, scalar_type(1) / scalar_type(dists.size()));
      }
      return d;
    }
    scalar_type operator()(const base_node &P, dal::bit_vector &bv) const {
      if (with_min) {
	scalar_type d = vd[0] = (*(dists[0]))(P);
	bool ok = (d > -SEPS);
	for (size_type k = 1; k < dists.size(); ++k) {
	  vd[k] = (*(dists[k]))(P); if (vd[k] <= -SEPS) ok = false;
	  d = std::min(d,vd[k]);
	}
	for (size_type k = 0; ok && k < dists.size(); ++k) {
	  if (vd[k] < SEPS) (*(dists[k]))(P, bv);
	}
	return d;
      }
      else { // essai raté ...
	vd[0] = (*(dists[0]))(P);
	bool ok = (vd[0] > -SEPS);
	for (size_type k = 1; k < dists.size(); ++k) {
	  vd[k] = (*(dists[k]))(P); if (vd[k] <= -SEPS) ok = false;
	}
	for (size_type k = 0; ok && k < dists.size(); ++k) {
	  if (vd[k] < SEPS) (*(dists[k]))(P, bv);
	}
	return operator()(P);
      }
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      for (size_type k = 0; k < dists.size(); ++k)
	dists[k]->register_constraints(list); 
    }
    scalar_type grad(const base_node &P, base_small_vector &G) const {
      scalar_type d;
      if (with_min /* || gmm::abs(d) < SEPS*/) {
	d = (*(dists[0]))(P);
	size_type i = 0;
	for (size_type k = 1; k < dists.size(); ++k) {
	  scalar_type d2 = (*(dists[k]))(P);
	  if (d2 < d) { d = d2; i = k; }
	}
	return dists[i]->grad(P, G);
      }
      else { // essai raté ...
	d = operator()(P);
	base_small_vector Gloc;
	for (size_type k = 0; k < dists.size(); ++k) {
	  dists[k]->grad(P, Gloc);
	  if (isin) 
	    Gloc *= -gmm::neg(vd[k]);
	  else 
	    Gloc *= pow(d, scalar_type(dists.size())) / vd[k];
	  if (!k) G = Gloc; else G += Gloc;
	}
	if (isin)
	  G *= scalar_type(1) / d;
	else
	  G /= pow(d, scalar_type(dists.size()-1)) * scalar_type(dists.size());
	return d;
      }
    }
    void hess(const base_node &P, base_matrix &H) const {
      scalar_type d = (*(dists[0]))(P);
      if (with_min || gmm::abs(d) < SEPS) {
	size_type i = 0;
	for (size_type k = 1; k < dists.size(); ++k) {
	  scalar_type d2 = (*(dists[k]))(P);
	  if (d2 < d) { d = d2; i = k; }
	}
	dists[i]->hess(P, H);
      }
      else {
	GMM_ASSERT1(false, "Sorry, to be done");
      }
    }
  };

  inline pmesher_signed_distance new_mesher_union
  (const std::vector<pmesher_signed_distance> &dists)
  { return std::make_shared<mesher_union>(dists); }

  inline pmesher_signed_distance new_mesher_union
  (const pmesher_signed_distance &a,
   const pmesher_signed_distance &b,
   const pmesher_signed_distance &c = pmesher_signed_distance(),
   const pmesher_signed_distance &d = pmesher_signed_distance(),
   const pmesher_signed_distance &e = pmesher_signed_distance(),
   const pmesher_signed_distance &f = pmesher_signed_distance(),
   const pmesher_signed_distance &g = pmesher_signed_distance(),
   const pmesher_signed_distance &h = pmesher_signed_distance(),
   const pmesher_signed_distance &i = pmesher_signed_distance(),
   const pmesher_signed_distance &j = pmesher_signed_distance(),
   const pmesher_signed_distance &k = pmesher_signed_distance(),
   const pmesher_signed_distance &l = pmesher_signed_distance(),
   const pmesher_signed_distance &m = pmesher_signed_distance(),
   const pmesher_signed_distance &n = pmesher_signed_distance(),
   const pmesher_signed_distance &o = pmesher_signed_distance(),
   const pmesher_signed_distance &p = pmesher_signed_distance(),
   const pmesher_signed_distance &q = pmesher_signed_distance(),
   const pmesher_signed_distance &r = pmesher_signed_distance(),
   const pmesher_signed_distance &s = pmesher_signed_distance(),
   const pmesher_signed_distance &t = pmesher_signed_distance()) {
    return std::make_shared<mesher_union>
      (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t);
  }
 


  class mesher_intersection : public mesher_signed_distance {
    std::vector<pmesher_signed_distance> dists;
    mutable std::vector<scalar_type> vd;

    // const mesher_signed_distance &a, &b;
  public:
    
    mesher_intersection(const std::vector<pmesher_signed_distance>
			&dists_) : dists(dists_) 
    { vd.resize(dists.size()); }
    
    mesher_intersection
    (const pmesher_signed_distance &a,
     const pmesher_signed_distance &b,
     const pmesher_signed_distance &c = pmesher_signed_distance(),
     const pmesher_signed_distance &d = pmesher_signed_distance(),
     const pmesher_signed_distance &e = pmesher_signed_distance(),
     const pmesher_signed_distance &f = pmesher_signed_distance(),
     const pmesher_signed_distance &g = pmesher_signed_distance(),
     const pmesher_signed_distance &h = pmesher_signed_distance(),
     const pmesher_signed_distance &i = pmesher_signed_distance(),
     const pmesher_signed_distance &j = pmesher_signed_distance(),
     const pmesher_signed_distance &k = pmesher_signed_distance(),
     const pmesher_signed_distance &l = pmesher_signed_distance(),
     const pmesher_signed_distance &m = pmesher_signed_distance(),
     const pmesher_signed_distance &n = pmesher_signed_distance(),
     const pmesher_signed_distance &o = pmesher_signed_distance(),
     const pmesher_signed_distance &p = pmesher_signed_distance(),
     const pmesher_signed_distance &q = pmesher_signed_distance(),
     const pmesher_signed_distance &r = pmesher_signed_distance(),
     const pmesher_signed_distance &s = pmesher_signed_distance(),
     const pmesher_signed_distance &t = pmesher_signed_distance()) {
      dists.push_back(a); dists.push_back(b);
      if (c) dists.push_back(c); if (d) dists.push_back(d);
      if (e) dists.push_back(e); if (f) dists.push_back(f);
      if (g) dists.push_back(g); if (h) dists.push_back(h);
      if (i) dists.push_back(i); if (j) dists.push_back(j);
      if (k) dists.push_back(k); if (l) dists.push_back(l);
      if (m) dists.push_back(m); if (n) dists.push_back(n);
      if (o) dists.push_back(o); if (p) dists.push_back(p);
      if (q) dists.push_back(q); if (r) dists.push_back(r);
      if (s) dists.push_back(s); if (t) dists.push_back(t);
      vd.resize(dists.size());
    }
    bool bounding_box(base_node &bmin, base_node &bmax) const {
      base_node bmin2, bmax2;
      bool first;
      bool b = dists[0]->bounding_box(bmin, bmax); first = !b;
      for (size_type k = 1; k < dists.size(); ++k) {
	bool bb = dists[k]->bounding_box(bmin2, bmax2);
	for (unsigned i=0; i < bmin.size() && bb && !first; ++i) { 
	  bmin[i] = std::max(bmin[i],bmin2[i]);
	  bmax[i] = std::max(std::min(bmax[i],bmax2[i]), bmin[i]);
	}
	if (first && bb) { bmin = bmin2; bmax = bmax2; first = false; }
	b = b || bb;
      }
      return b;
    }
    virtual scalar_type operator()(const base_node &P) const {
      scalar_type d = (*(dists[0]))(P);
      for (size_type k = 1; k < dists.size(); ++k)
	d = std::max(d, (*(dists[k]))(P));
      return d;

    }
    scalar_type operator()(const base_node &P, dal::bit_vector &bv) const {
      scalar_type d = vd[0] = (*(dists[0]))(P);
      bool ok = (d < SEPS);
      for (size_type k = 1; k < dists.size(); ++k) {
	vd[k] = (*(dists[k]))(P); if (vd[k] >= SEPS) ok = false;
	d = std::min(d,vd[k]);
      }
      for (size_type k = 0; ok && k < dists.size(); ++k) {
	if (vd[k] > -SEPS) (*(dists[k]))(P, bv);
      }
      return d;
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      for (size_type k = 0; k < dists.size(); ++k) {
	dists[k]->register_constraints(list); 
      }
    }
    scalar_type grad(const base_node &P, base_small_vector &G) const {
      scalar_type d = (*(dists[0]))(P);
      size_type i = 0;
      for (size_type k = 1; k < dists.size(); ++k) {
	scalar_type d2 = (*(dists[k]))(P);
	if (d2 > d) { d = d2; i = k; }
      }
      return dists[i]->grad(P, G);
    }
    void hess(const base_node &P, base_matrix &H) const {
      scalar_type d = (*(dists[0]))(P);
      size_type i = 0;
      for (size_type k = 1; k < dists.size(); ++k) {
	scalar_type d2 = (*(dists[k]))(P);
	if (d2 > d) { d = d2; i = k; }
      }
      dists[i]->hess(P, H);
    }
  };

  inline pmesher_signed_distance new_mesher_intersection
  (const std::vector<pmesher_signed_distance> &dists)
  { return std::make_shared<mesher_intersection>(dists); }

  inline pmesher_signed_distance new_mesher_intersection
  (const pmesher_signed_distance &a,
   const pmesher_signed_distance &b,
   const pmesher_signed_distance &c = pmesher_signed_distance(),
   const pmesher_signed_distance &d = pmesher_signed_distance(),
   const pmesher_signed_distance &e = pmesher_signed_distance(),
   const pmesher_signed_distance &f = pmesher_signed_distance(),
   const pmesher_signed_distance &g = pmesher_signed_distance(),
   const pmesher_signed_distance &h = pmesher_signed_distance(),
   const pmesher_signed_distance &i = pmesher_signed_distance(),
   const pmesher_signed_distance &j = pmesher_signed_distance(),
   const pmesher_signed_distance &k = pmesher_signed_distance(),
   const pmesher_signed_distance &l = pmesher_signed_distance(),
   const pmesher_signed_distance &m = pmesher_signed_distance(),
   const pmesher_signed_distance &n = pmesher_signed_distance(),
   const pmesher_signed_distance &o = pmesher_signed_distance(),
   const pmesher_signed_distance &p = pmesher_signed_distance(),
   const pmesher_signed_distance &q = pmesher_signed_distance(),
   const pmesher_signed_distance &r = pmesher_signed_distance(),
   const pmesher_signed_distance &s = pmesher_signed_distance(),
   const pmesher_signed_distance &t = pmesher_signed_distance()) {
    return std::make_shared<mesher_intersection>
      (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t);
  }
 


  class mesher_setminus : public mesher_signed_distance {
    pmesher_signed_distance a, b;
  public:
    mesher_setminus(const pmesher_signed_distance& a_,
		    const pmesher_signed_distance &b_) : a(a_), b(b_) {}
    bool bounding_box(base_node &bmin, base_node &bmax) const
    { return a->bounding_box(bmin,bmax); }
    scalar_type operator()(const base_node &P, dal::bit_vector &bv) const {
      scalar_type da = (*a)(P), db = -(*b)(P);
      if (da < SEPS && db < SEPS) {
	if (da > -SEPS) (*a)(P, bv);
	if (db > -SEPS) (*b)(P, bv);
      }
      return std::max(da, db);
    }
    scalar_type operator()(const base_node &P) const
    { return std::max((*a)(P),-(*b)(P)); }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      a->register_constraints(list); b->register_constraints(list);
    }
    scalar_type grad(const base_node &P, base_small_vector &G) const {
      scalar_type da = (*a)(P), db = -(*b)(P);
      if (da > db) return a->grad(P, G);
      else { b->grad(P, G); G *= scalar_type(-1); return db; }
    }
    void hess(const base_node &P, base_matrix &H) const {
      scalar_type da = (*a)(P), db = -(*b)(P);
      if (da > db) a->hess(P, H);
      else { b->hess(P, H); gmm::scale(H, scalar_type(-1)); }
    }

  };
  
  inline pmesher_signed_distance new_mesher_setminus
  (const pmesher_signed_distance &a, const pmesher_signed_distance &b)
  { return std::make_shared<mesher_setminus>(a, b); }


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
      bv[id] = (gmm::abs(d) < SEPS);
      return d;
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      id = list.size(); list.push_back(this);
    }
    scalar_type grad(const base_node &P, base_small_vector &G) const {
      G = P; G -= x0;
      gmm::add(gmm::scaled(n, -gmm::vect_sp(G, n)), G);
      scalar_type e = gmm::vect_norm2(G), d = e - R;
      while (e == scalar_type(0)) {
	gmm::fill_random(G);
	gmm::add(gmm::scaled(n, -gmm::vect_sp(G, n)), G);
	e = gmm::vect_norm2(G);
      }
      G /= e;
      return d;
    }
    void hess(const base_node &, base_matrix &) const {
      GMM_ASSERT1(false, "Sorry, to be done");
    }
  };

  inline pmesher_signed_distance new_mesher_tube(base_node x0, base_node n,
						 scalar_type R)
  { return std::make_shared<mesher_tube>(x0,n,R); }


  class mesher_cylinder : public mesher_signed_distance {
    base_node x0; base_small_vector n;
    scalar_type L, R;
    pmesher_signed_distance t, p1, p2, i1;
  public:
    mesher_cylinder(const base_node &c, const base_small_vector &no,
		    scalar_type L_, scalar_type R_)
      : x0(c), n(no/gmm::vect_norm2(no)), L(L_), R(R_),
	t(new_mesher_tube(x0, n, R)),
	p1(new_mesher_half_space(x0, n)),
	p2(new_mesher_half_space(x0+n*L, -1.0 * n)),
	i1(new_mesher_intersection(p1, p2, t)) {}
    bool bounding_box(base_node &bmin, base_node &bmax) const {
      base_node x1(x0+n*L);
      bmin = bmax = x0;
      for (unsigned i = 0; i < gmm::vect_size(x0); ++i) {
	bmin[i] = std::min(x0[i], x1[i]) - R;
	bmax[i] = std::max(x0[i], x1[i]) + R;
      }
      return true;
    }
    virtual scalar_type operator()(const base_node &P) const {return (*i1)(P); }
    virtual scalar_type operator()(const base_node &P,
				   dal::bit_vector& bv) const
    { return (*i1)(P, bv); }
    scalar_type grad(const base_node &P, base_small_vector &G) const
      { return i1->grad(P, G); }
    void hess(const base_node &, base_matrix &) const {
      GMM_ASSERT1(false, "Sorry, to be done");
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const
    { i1->register_constraints(list); }
  };

  inline pmesher_signed_distance new_mesher_cylinder
  (const base_node &c, const base_small_vector &no,
   scalar_type L, scalar_type R)
  { return std::make_shared<mesher_cylinder>(c,no,L,R); }


  class mesher_infinite_cone : public mesher_signed_distance {
    // uses the true distance to a cone.
    base_node x0; base_node n; scalar_type alpha;
  public:
      mesher_infinite_cone(base_node x0_, base_node n_, scalar_type alpha_)
	: x0(x0_), n(n_), alpha(alpha_)
    { n /= gmm::vect_norm2(n); }
    bool bounding_box(base_node &, base_node &) const
    { return false; }
    virtual scalar_type operator()(const base_node &P) const {
      base_node v(P); v -= x0;
      scalar_type v_n = gmm::vect_sp(v, n);
      gmm::add(gmm::scaled(n, -v_n), v);
      return gmm::vect_norm2(v) * cos(alpha) - gmm::abs(v_n) * sin(alpha);
    }
    virtual scalar_type operator()(const base_node &P,
				   dal::bit_vector &bv) const {
      scalar_type d = (*this)(P);
      bv[id] = (gmm::abs(d) < SEPS);
      return d;
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const {
      id = list.size(); list.push_back(this);
    }
    scalar_type grad(const base_node &P, base_small_vector &v) const {
      v = P; v -= x0;
      scalar_type v_n = gmm::vect_sp(v, n);
      gmm::add(gmm::scaled(n, -v_n), v);
      scalar_type no = gmm::vect_norm2(v);
      scalar_type d = no * cos(alpha) - gmm::abs(v_n) * sin(alpha);
      while (no ==  scalar_type(0)) {
	gmm::fill_random(v);
	gmm::add(gmm::scaled(n, -gmm::vect_sp(v, n)), v);
	no = gmm::vect_norm2(v);
      }
      v *= cos(alpha) / no;
      v -= (sin(alpha) * gmm::sgn(v_n)) * n;
      return d;
    }
    void hess(const base_node &, base_matrix &) const {
      GMM_ASSERT1(false, "Sorry, to be done");
    }
  };

  inline pmesher_signed_distance new_mesher_infinite_cone
  (base_node x0, base_node n, scalar_type alpha)
  { return std::make_shared<mesher_infinite_cone>(x0,n,alpha); }


  class mesher_cone : public mesher_signed_distance {
    base_node x0; base_small_vector n;
    scalar_type L, alpha;
    pmesher_signed_distance t, p1, p2, i1;
  public:
    mesher_cone(const base_node &c, const base_small_vector &no,
		    scalar_type L_, scalar_type alpha_)
      : x0(c), n(no/gmm::vect_norm2(no)), L(L_), alpha(alpha_),
	t(new_mesher_infinite_cone(x0, n, alpha)),
	p1(new_mesher_half_space(x0, n)),
	p2(new_mesher_half_space(x0+n*L, -1.0 * n)),
	i1(new_mesher_intersection(p1, p2, t)) {}
    bool bounding_box(base_node &bmin, base_node &bmax) const {
      base_node x1(x0+n*L);
      scalar_type R = L * sin(alpha);
      bmin = bmax = x0;
      for (unsigned i = 0; i < gmm::vect_size(x0); ++i) {
	bmin[i] = std::min(x0[i], x1[i]) - R;
	bmax[i] = std::max(x0[i], x1[i]) + R;
      }
      return true;
    }
    virtual scalar_type operator()(const base_node &P) const {return (*i1)(P); }
    virtual scalar_type operator()(const base_node &P,
				   dal::bit_vector& bv) const
    { return (*i1)(P, bv); }
    scalar_type grad(const base_node &P, base_small_vector &G) const
      { return i1->grad(P, G); }
    void hess(const base_node &, base_matrix &) const {
      GMM_ASSERT1(false, "Sorry, to be done");
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const
    { i1->register_constraints(list); }
  };

  inline pmesher_signed_distance new_mesher_cone
  (const base_node &c, const base_small_vector &no,
   scalar_type L, scalar_type alpha)
  { return std::make_shared<mesher_cone>(c,no,L,alpha); }


  class mesher_ellipse : public mesher_signed_distance {
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
      return (1.-lambda)*gmm::vect_norm2(v1);
    }
    virtual scalar_type operator()(const base_node &P,
				   dal::bit_vector& bv) const {
      scalar_type d = this->operator()(P);
      bv[id] = (gmm::abs(d) < SEPS);
      return d;
    }
    scalar_type grad(const base_node &, base_small_vector &) const
    { GMM_ASSERT1(false, "Sorry, to be done"); return 0.; }
    void hess(const base_node &, base_matrix &) const {
      GMM_ASSERT1(false, "Sorry, to be done");
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>& list) const
    { id = list.size(); list.push_back(this); }
  };

  inline pmesher_signed_distance new_mesher_ellipse
  (const base_node &center, const base_small_vector &no,
   scalar_type r, scalar_type R)
  { return std::make_shared<mesher_ellipse>(center,no,r,R); }



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
      bv[id] = (gmm::abs(d) < SEPS);
      return d;
    }
    scalar_type grad(const base_node &P, base_small_vector &G) const {
      G.resize(3);
      scalar_type x = P[0], y = P[1], z = P[2], c = sqrt(x*x + y*y), d(0);
      if (c == 0.) {
	d = R - r;
	gmm::fill_random(G); G[2] = 0.0; G /= gmm::vect_norm2(G);
      }
      else {
	scalar_type w = 1. - R / c, e = sqrt(gmm::sqr(c-R) + z*z);
	d = e - r;
	if (e == 0.) {
	  gmm::fill_random(G); G[0] = x; G[1] = y; G /= gmm::vect_norm2(G);
	}
	else {
	  G[0] = x * w / e; G[1] = y * w / e; G[2] = z / e;
	}
      }
      return d;
    }
    void hess(const base_node &, base_matrix &) const {
      GMM_ASSERT1(false, "Sorry, to be done");
    }
    virtual void register_constraints(std::vector<const
				      mesher_signed_distance*>&list) const
    { id = list.size(); list.push_back(this); }
  };


  inline pmesher_signed_distance new_mesher_torus(scalar_type R, scalar_type r)
  { return std::make_shared<mesher_torus>(R,r); }


  // interface with qhull
  void delaunay(const std::vector<base_node> &pts,
		gmm::dense_matrix<size_type>& simplexes);

  
  // mesher
  void build_mesh(mesh &m, const pmesher_signed_distance& dist_,
		  scalar_type h0, const std::vector<base_node> &fixed_points
		  = std::vector<base_node>(), size_type K = 1, int noise = -1,
		  size_type iter_max = 500, int prefind = 1,
		  scalar_type dist_point_hull = 4,
		  scalar_type boundary_threshold_flatness = 0.11);

  // exported functions
  bool try_projection(const mesher_signed_distance& dist, base_node &X,
		      bool on_surface = false);
  bool pure_multi_constraint_projection
  (const std::vector<const mesher_signed_distance*> &list_constraints,
   base_node &X, const dal::bit_vector &cts);
  scalar_type curvature_radius_estimate(const mesher_signed_distance &dist,
					base_node X, bool proj = false);
  scalar_type min_curvature_radius_estimate
  (const std::vector<const mesher_signed_distance*> &list_constraints,
   const base_node &X, const dal::bit_vector &cts, size_type hide_first = 0);

}

#endif /* GETFEM_MESHER_H__ */
