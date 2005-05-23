// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_fem.cc : definition of the finite element methods
//           
// Date    : December 21, 1999.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1999-2005 Yves Renard
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


#include <dal_singleton.h>
#include <dal_tree_sorted.h>
#include <dal_algobase.h>
#include <dal_naming_system.h>
#include <getfem_fem.h>

/* do not read this file ! */
#include <getfem_gauss_lobatto_fem_coef.h>
#include <getfem_integration.h> /* for gauss-lobatto points */
namespace getfem
{
  typedef dal::naming_system<virtual_fem>::param_list fem_param_list;

  const base_matrix& fem_interpolation_context::M() const {
    if (!have_pgt() || !have_G() || !have_pf())
      DAL_THROW(dal::failure_error, "cannot compute M");
    M_.resize(pf_->nb_base(convex_num()), pf_->nb_dof(convex_num()));
    pf_->mat_trans(M_,G(),pgt());
    return M_;
  }
  
  size_type fem_interpolation_context::convex_num() const { 
    if (convex_num_ == size_type(-1)) 
      DAL_INTERNAL_ERROR(""); 
    return convex_num_; 
  }

  void fem_interpolation_context::base_value(base_tensor& t) const {
    if (pf()->is_on_real_element())
      pf()->real_base_value(*this, t);
    else {
      if (have_pfp()) t=pfp_->val(ii());
      else pf()->base_value(xref(), t);
    }
  }

  void fem_interpolation_context::grad_base_value(base_tensor& t) const {
    if (pf()->is_on_real_element())
      pf()->real_grad_base_value(*this, t);
    else {
      if (have_pfp()) {
	t.mat_transp_reduction(pfp_->grad(ii()), B(), 2);
      } else {
	base_tensor u;
	pf()->grad_base_value(xref(), u);
	if (u.size()) /* only if the FEM can provide grad_base_value */
	  t.mat_transp_reduction(u, B(), 2);
      }
    }
  }

  void fem_interpolation_context::hess_base_value(base_tensor& t) const {
    if (pf()->is_on_real_element())
      pf()->real_hess_base_value(*this, t);
    else {
      base_tensor tt;
      if (have_pfp()) {
	tt = pfp()->hess(ii());
      } else {
	pf()->hess_base_value(xref(), tt);
      }
      if (tt.size()) { /* only if the FEM can provide grad_base_value */
	bgeot::multi_index mim(3);
	mim[2] = gmm::sqr(tt.sizes()[2]); mim[1] = tt.sizes()[1];
	mim[0] = tt.sizes()[0];
	tt.adjust_sizes(mim);
	t.mat_transp_reduction(tt, B3(), 2);
	if (have_pfp()) {
	  tt.mat_transp_reduction(pfp()->grad(ii()), B32(), 2);
	} else {
	  base_tensor u;
	  pf()->grad_base_value(xref(), u);
	  tt.mat_transp_reduction(u, B32(), 2);
	}
	t -= tt;
      }
    }
  }

  void fem_interpolation_context::set_pfp(pfem_precomp newpfp) {
    pfp_ = newpfp;
    if (pfp_) { pf_ = pfp()->get_pfem(); }
    else pf_ = 0;
    M_.resize(0,0);
  }

  void fem_interpolation_context::set_pf(pfem newpf) {
    set_pfp(0);
    pf_ = newpf;
  }

  fem_interpolation_context::fem_interpolation_context() :
    bgeot::geotrans_interpolation_context(), pf_(0), pfp_(0), 
    convex_num_(size_type(-1)) {}
  fem_interpolation_context::fem_interpolation_context
  (bgeot::pgeotrans_precomp pgp__, pfem_precomp pfp__, size_type ii__, 
   const base_matrix& G__, size_type convex_num__) : 
    bgeot::geotrans_interpolation_context(pgp__,ii__,G__), 
    convex_num_(convex_num__) { set_pfp(pfp__); }
  fem_interpolation_context::fem_interpolation_context
  (bgeot::pgeometric_trans pgt__, pfem_precomp pfp__, size_type ii__, 
   const base_matrix& G__, size_type convex_num__) :
    bgeot::geotrans_interpolation_context(pgt__,&pfp__->get_point_tab(), ii__, G__),
    convex_num_(convex_num__) {
      set_pfp(pfp__); 
    }
  fem_interpolation_context::fem_interpolation_context(
   bgeot::pgeometric_trans pgt__, pfem pf__,
   const base_node& xref__,const base_matrix& G__, size_type convex_num__) :
    bgeot::geotrans_interpolation_context(pgt__,xref__,G__),
    pf_(pf__), pfp_(0), convex_num_(convex_num__) {}
 
 
  void virtual_fem::real_base_value(const fem_interpolation_context &c, 
				    base_tensor &t) const {
    c.base_value(t);
  }
  void virtual_fem::real_grad_base_value(const fem_interpolation_context &c, 
				    base_tensor &t) const {
    c.grad_base_value(t);
  }
  void virtual_fem::real_hess_base_value(const fem_interpolation_context &c, 
					 base_tensor &t) const {
    c.hess_base_value(t);
  }

  /* ******************************************************************** */
  /*	Class for description of an interpolation dof.                    */
  /* ******************************************************************** */

  enum ddl_type { LAGRANGE, NORM_DERIVATIVE, DERIVATIVE, MEAN_VALUE, BUBBLE1, 
		  LAGRANGE_NONCONFORMING, GLOBAL_DOF};

  enum coord_type { NORMAL = -2, TANGENTIAL = -1, FIRST = 0, SECOND, THIRD };

  struct ddl_elem {
    ddl_type t;
    dal::int16_type hier_degree;
    short_type hier_raff;
    bool operator < (const ddl_elem &l) const {
      if (t < l.t) return true; if (t > l.t) return false; 
      if (hier_degree < l.hier_degree) return true; 
      if (hier_degree > l.hier_degree) return false;
      if (hier_raff < l.hier_raff) return true; return false;
    }
    ddl_elem(ddl_type s = LAGRANGE, dal::int16_type k = -1, short_type l = 0)
      : t(s), hier_degree(k), hier_raff(l) {}
  };

  struct dof_description {
    std::vector<ddl_elem> ddl_desc;
    bool linkable;
    coord_type coord_index;
    size_type xfem_index;
    bool all_faces;

    dof_description(void)
    { linkable = true; all_faces = false; coord_index = FIRST; xfem_index = 0; }
  };

  struct dof_description_comp__ {
    int operator()(const dof_description &m, const dof_description &n) const;
  };

  // ATTENTION : en cas de modif, changer aussi dof_description_compare,
  //             product_dof, et dof_hierarchical_compatibility.
  int dof_description_comp__::operator()(const dof_description &m,
					 const dof_description &n) const { 
    int nn = dal::lexicographical_less<std::vector<ddl_elem> >()
      (m.ddl_desc, n.ddl_desc);
    if (nn < 0) return -1; if (nn > 0) return 1;
    nn = int(m.linkable) - int(n.linkable);
    if (nn < 0) return -1; if (nn > 0) return 1;
    nn = int(m.coord_index) - int(n.coord_index);
    if (nn < 0) return -1; if (nn > 0) return 1;
    nn = int(m.xfem_index) - int(n.xfem_index);
    if (nn < 0) return -1; if (nn > 0) return 1;
    nn = int(m.all_faces) - int(n.all_faces);
    if (nn < 0) return -1; if (nn > 0) return 1;
    return 0;
  }

  typedef dal::dynamic_tree_sorted<dof_description, dof_description_comp__> dof_d_tab;

  pdof_description lagrange_dof(dim_type n) {
    static dim_type n_old = dim_type(-2);
    static pdof_description p_old = 0;
    if (n != n_old) {
      dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
      dof_description l;
      l.ddl_desc.resize(n);
      std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(LAGRANGE));
      p_old = &(tab[tab.add_norepeat(l)]);
      n_old = n;
    }
    return p_old;
  } 

  pdof_description lagrange_0_dof(dim_type n) {
    static dim_type n_old = dim_type(-2);
    static pdof_description p_old = 0;
    if (n != n_old) {
      dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
      dof_description l;
      l.all_faces = true;
      l.ddl_desc.resize(n);
      std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(LAGRANGE));
      p_old = &(tab[tab.add_norepeat(l)]);
      n_old = n;
    }
    return p_old;
  }

  pdof_description deg_hierarchical_dof(pdof_description p, int deg) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l = *p;
    for (size_type i = 0; i < l.ddl_desc.size(); ++i)
      l.ddl_desc[i].hier_degree = deg;
    return &(tab[tab.add_norepeat(l)]);
  }

  size_type reserve_xfem_index(void) {
    static size_type ind = 100;
    return ind += 1000;
  }

  pdof_description xfem_dof(pdof_description p, size_type ind) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l = *p; l.xfem_index = ind;
    return &(tab[tab.add_norepeat(l)]);
  }


  pdof_description to_coord_dof(pdof_description p, dim_type ct) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l = *p;
    l.coord_index = coord_type(ct);
    return &(tab[tab.add_norepeat(l)]);
  }
  
  pdof_description raff_hierarchical_dof(pdof_description p, short_type deg) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l = *p;
    for (size_type i = 0; i < l.ddl_desc.size(); ++i)
      l.ddl_desc[i].hier_raff = deg;
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description lagrange_nonconforming_dof(dim_type n) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l; l.linkable = false;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(LAGRANGE));
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description bubble1_dof(dim_type n) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(BUBBLE1));
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description derivative_dof(dim_type n, dim_type num_der) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(LAGRANGE));
    l.ddl_desc[num_der] = ddl_elem(DERIVATIVE);
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description norm_derivative_dof(dim_type n) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(NORM_DERIVATIVE));
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description mean_value_dof(dim_type n) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(MEAN_VALUE));
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description global_dof(dim_type n) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l;
    l.all_faces = true;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(),ddl_elem(GLOBAL_DOF));
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description product_dof(pdof_description a, pdof_description b) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    size_type nb1 = a->ddl_desc.size(), nb2 = b->ddl_desc.size();
    dof_description l;
    l.linkable = a->linkable && b->linkable;
    l.coord_index = std::max(a->coord_index, b->coord_index); // logique ?
    l.xfem_index = a->xfem_index;
    l.all_faces = a->all_faces || b->all_faces;
    if (a->xfem_index != b->xfem_index)
      DAL_THROW(failure_error, "Invalid product of dof");
    l.ddl_desc.resize(nb1+nb2);
    std::copy(a->ddl_desc.begin(), a->ddl_desc.end(), l.ddl_desc.begin());
    std::copy(b->ddl_desc.begin(), b->ddl_desc.end(), l.ddl_desc.begin()+nb1);
    
    {
      dal::int16_type deg = -1;
      for (size_type i = 0; i < l.ddl_desc.size(); ++i)
	deg = std::max(deg, l.ddl_desc[i].hier_degree);
      for (size_type i = 0; i < l.ddl_desc.size(); ++i)
	l.ddl_desc[i].hier_degree = deg;
    }
    { 
      short_type deg = 0;
      for (size_type i = 0; i < l.ddl_desc.size(); ++i)
	deg = std::max(deg, l.ddl_desc[i].hier_raff);
      for (size_type i = 0; i < l.ddl_desc.size(); ++i)
	l.ddl_desc[i].hier_raff = deg;
    }
    return &(tab[tab.add_norepeat(l)]);
  }

  // ATTENTION : en cas de modif, changer aussi
  //             dof_description_comp__::operator,
  //             product_dof, et dof_hierarchical_compatibility.
  int dof_description_compare(pdof_description a, pdof_description b) {
    int nn;
    if ((nn = int(a->coord_index) - int(b->coord_index)) != 0) return nn;
    if ((nn = int(a->linkable) - int(b->linkable)) != 0) return nn;
    if ((nn = int(a->xfem_index) - int(b->xfem_index)) != 0) return nn;
    std::vector<ddl_elem>::const_iterator
      ita = a->ddl_desc.begin(), itae = a->ddl_desc.end(),
      itb = b->ddl_desc.begin(), itbe = b->ddl_desc.end();
    for (; ita != itae && itb != itbe; ++ita, ++itb)
    {
      if ((nn = int(ita->t) - int (itb->t)) != 0) return nn;
      if ((nn = int(ita->hier_degree) - int (itb->hier_degree)) != 0)
	return nn;
    }
    for (; ita != itae; ++ita) if (ita->t != LAGRANGE) return 1;
    for (; itb != itbe; ++itb) if (itb->t != LAGRANGE) return -1;
    return 0;
  }

  bool dof_linkable(pdof_description a)
  { return a->linkable; }

  bool dof_compatibility(pdof_description a, pdof_description b)
  { return (dof_linkable(a) && dof_description_compare(a, b) == 0); }
  
  size_type dof_xfem_index(pdof_description a) 
  { return a->xfem_index; }

  bool dof_hierarchical_compatibility(pdof_description a, pdof_description b)
  { 
    if (a->coord_index != b->coord_index) return false;
    if (a->linkable != b->linkable) return false;
    if (a->xfem_index != b->xfem_index) return false;
    std::vector<ddl_elem>::const_iterator
      ita = a->ddl_desc.begin(), itae = a->ddl_desc.end(),
      itb = b->ddl_desc.begin(), itbe = b->ddl_desc.end();
    for (; ita != itae && itb != itbe; ++ita, ++itb)
    { if (ita->t != itb->t) return false; }
    for (; ita != itae; ++ita) if (ita->t != LAGRANGE) return false;
    for (; itb != itbe; ++itb) if (itb->t != LAGRANGE) return false;
    return true;
  }

  void virtual_fem::add_node(const pdof_description &d, const base_node &pt) {
    short_type nb = cv_node.nb_points();
    cv_node.points().resize(nb+1);
    cv_node.points()[nb] = pt;
    dof_types_.resize(nb+1);
    dof_types_[nb] = d;
    cvs_node.add_point_adaptative(nb, short_type(-1));
    for (short_type f = 0; f < cvs_node.nb_faces(); ++f)
      if (d->all_faces || gmm::abs(cvr->is_in_face(f, pt)) < 1.0E-7)
	cvs_node.add_point_adaptative(nb, f);
    pspt_valid = false;
  }

  void virtual_fem::init_cvs_node(void) {
    cvs_node.init_for_adaptative(cvr->structure());
    cv_node = bgeot::convex<base_node>(&cvs_node);
    pspt_valid = false;
  }

  void virtual_fem::unfreeze_cvs_node(void) {
    cv_node.structure() = &cvs_node;
    pspt_valid = false;
  }

  /* ******************************************************************** */
  /*	PK class.                                                         */
  /* ******************************************************************** */

  class PK_fem_ : public fem<base_poly> {
  public :
    void calc_base_func(base_poly &p, size_type i, short_type K) const;
    PK_fem_(dim_type nc, short_type k);
  };
  
  void PK_fem_::calc_base_func(base_poly &p, size_type i, short_type K) const {
    dim_type N = dim();
    base_poly l0(N, 0), l1(N, 0);
    bgeot::power_index w(N+1);
    l0.one(); l1.one(); p = l0;
    
    if (K != 0) {
      for (int nn = 0; nn < N; ++nn) l0 -= base_poly(N, 1, nn);
      
      w[0] = K;
      for (int nn = 1; nn <= N; ++nn) { 
	w[nn]=int(floor(0.5+((cv_node.points()[i])[nn-1]*opt_long_scalar_type(K))));
	w[0]-=w[nn];
      }
      
      for (int nn = 0; nn <= N; ++nn)
	for (int j = 0; j < w[nn]; ++j) {
	  if (nn == 0)
	    p *= (l0 * (opt_long_scalar_type(K) / opt_long_scalar_type(j+1))) 
	      - (l1 * (opt_long_scalar_type(j) / opt_long_scalar_type(j+1)));
	  else
	    p *= (base_poly(N,1,nn-1) * (opt_long_scalar_type(K)/opt_long_scalar_type(j+1))) 
	      - (l1 * (opt_long_scalar_type(j) / opt_long_scalar_type(j+1)));
	}
    }
  }
  
  PK_fem_::PK_fem_(dim_type nc, short_type k) {
    cvr = bgeot::simplex_of_reference(nc);
    dim_ = cvr->structure()->dim();
    is_equiv = is_pol = is_lag = true;
    es_degree = k;
    
    init_cvs_node();
    bgeot::pconvex_ref cvn = bgeot::simplex_of_reference(nc, k);
    size_type R = cvn->nb_points();
    for (size_type i = 0; i < R; ++i)
      add_node(k==0 ? lagrange_0_dof(nc) : lagrange_dof(nc), cvn->points()[i]);
    
    base_.resize(R);
    for (size_type r = 0; r < R; r++) calc_base_func(base_[r], r, k);
  }

  static pfem PK_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    if (n <= 0 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");
    virtual_fem *p = new PK_fem_(n, k);
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }


  /* ******************************************************************** */
  /*	Tensorial product of fem (for polynomial fem).                    */
  /* ******************************************************************** */

  struct tproduct_femi : public fem<base_poly> { 
    tproduct_femi(ppolyfem fi1, ppolyfem fi2);
  };

  tproduct_femi::tproduct_femi(ppolyfem fi1, ppolyfem fi2) {
    if (fi2->target_dim() != 1) std::swap(fi1, fi2);
    if (fi2->target_dim() != 1) 
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    is_pol = true;
    is_equiv = fi1->is_equivalent() && fi2->is_equivalent();
    if (!is_equiv) 
      DAL_THROW(to_be_done_error, 
		"Product of non equivalent elements not available, sorry.");
    is_lag = fi1->is_lagrange() && fi2->is_lagrange();;
    es_degree = fi1->estimated_degree() + fi2->estimated_degree();
    
    bgeot::convex<base_node> cv 
      = bgeot::convex_direct_product(fi1->node_convex(0), fi2->node_convex(0));
    cvr = bgeot::convex_ref_product(fi1->ref_convex(0), fi2->ref_convex(0));
    dim_ = cvr->structure()->dim();
    init_cvs_node();
    
    ntarget_dim = fi2->target_dim();
    base_.resize(cv.nb_points() * ntarget_dim);
    size_type i, j, r;
    for (j = 0, r = 0; j < fi2->nb_dof(0); ++j)
      for (i = 0; i < fi1->nb_dof(0); ++i, ++r)
	add_node(product_dof(fi1->dof_types()[i], fi2->dof_types()[j]),
		 cv.points()[r]);
    
    for (j = 0, r = 0; j < fi2->nb_base_components(0); j++)
      for (i = 0; i < fi1->nb_base_components(0); i++, ++r) {
	base_[r] = fi1->base()[i];
	base_[r].direct_product(fi2->base()[j]); 
      }
  }

  static pfem product_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	  "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 1 || params[1].type() != 1)
      DAL_THROW(failure_error, "Bad type of parameters");
    pfem pf1 = params[0].method();
    pfem pf2 = params[1].method();
    if (!(pf1->is_polynomial() && pf2->is_polynomial()))
      DAL_THROW(failure_error, "Bad parameters");
    virtual_fem *p = new tproduct_femi(ppolyfem(pf1.get()),
				       ppolyfem(pf2.get()));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*    Generic Hierarchical fem (for polynomial fem). To be interfaced.  */
  /* ******************************************************************** */

  struct thierach_femi : public fem<base_poly> { 
    thierach_femi(ppolyfem fi1, ppolyfem fi2);
  };


  thierach_femi::thierach_femi(ppolyfem fi1, ppolyfem fi2)
    : fem<base_poly>(*fi1) {
    if (fi2->target_dim() != fi1->target_dim())
      DAL_THROW(dimension_error, "dimensions mismatch.");
    if (fi2->basic_structure(0) != fi1->basic_structure(0))
      DAL_THROW(failure_error, "Incompatible elements.");
    if (!(fi1->is_equivalent() &&  fi2->is_equivalent()))
      DAL_THROW(to_be_done_error,
	   "Sorry, no hierachical construction for non tau-equivalent fem.");
    es_degree = fi2->estimated_degree();
    is_lag = false;
    unfreeze_cvs_node();
    for (size_type i = 0; i < fi2->nb_dof(0); ++i) {
      bool found = false;
      for (size_type j = 0; j < fi1->nb_dof(0); ++j) {
	if ( dal::lexicographical_less<base_node,
	     dal::approx_less<scalar_type> >()
	     (fi2->node_of_dof(0,i), fi1->node_of_dof(0,j)) == 0
	     && dof_hierarchical_compatibility(fi2->dof_types()[i],
					       fi1->dof_types()[j]))
	    { found = true; break; }
      }
      if (!found) {
	add_node(deg_hierarchical_dof(fi2->dof_types()[i], 
				      fi1->estimated_degree()),
		 fi2->node_of_dof(0,i));
	base_.resize(nb_dof(0));
	base_[nb_dof(0)-1] = (fi2->base())[i];
      }
    }
  }

  struct thierach_femi_comp : public fem<polynomial_composite>
  { 
    thierach_femi_comp(ppolycompfem fi1, ppolycompfem fi2);  
  };

  thierach_femi_comp::thierach_femi_comp(ppolycompfem fi1, ppolycompfem fi2)
    : fem<polynomial_composite>(*fi1) {
    if (fi2->target_dim() != fi1->target_dim())
      DAL_THROW(dimension_error, "dimensions mismatch.");
    if (fi2->basic_structure(0) != fi1->basic_structure(0))
      DAL_THROW(failure_error, "Incompatible elements.");
    if (!(fi1->is_equivalent() && fi2->is_equivalent()))
      DAL_THROW(to_be_done_error,
	   "Sorry, no hierachical construction for non tau-equivalent fem.");
    es_degree = std::max(fi2->estimated_degree(), fi1->estimated_degree());
    
    is_lag = false;
    hier_raff = fi1->hierarchical_raff() + 1;
    unfreeze_cvs_node();
    for (size_type i = 0; i < fi2->nb_dof(0); ++i) {
      bool found = false;
      for (size_type j = 0; j < fi1->nb_dof(0); ++j) {
	if ( dal::lexicographical_less<base_node,
	     dal::approx_less<scalar_type> >()
	     (fi2->node_of_dof(0,i), fi1->node_of_dof(0,j)) == 0
	     && dof_hierarchical_compatibility(fi2->dof_types()[i],
					       fi1->dof_types()[j]))
	  { found = true; break; }
      }
      if (!found) {
	add_node(raff_hierarchical_dof(fi2->dof_types()[i], hier_raff),
		 fi2->node_of_dof(0,i));
	base_.resize(nb_dof(0));
	base_[nb_dof(0)-1] = (fi2->base())[i];
      }
    }
  }

  static pfem gen_hierarchical_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	  "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 1 || params[1].type() != 1)
      DAL_THROW(failure_error, "Bad type of parameters");
    pfem pf1 = params[0].method();
    pfem pf2 = params[1].method();
    if (pf1->is_polynomial() && pf2->is_polynomial())
      return new thierach_femi(ppolyfem(pf1.get()),ppolyfem(pf2.get()));
    if (pf1->is_polynomialcomp() && pf2->is_polynomialcomp()) {
      virtual_fem *p = new thierach_femi_comp(ppolycompfem(pf1.get()),
					      ppolycompfem(pf2.get()));
      dependencies.push_back(p->ref_convex(0));
      dependencies.push_back(p->node_tab(0));
      return p;
    }
    DAL_THROW(failure_error, "Bad parameters");
  }

  /* ******************************************************************** */
  /* PK hierarchical fem.                                                 */
  /* ******************************************************************** */

  static pfem PK_hierarch_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01)), s;
    if (n <= 0 || n >= 100 || k <= 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");
    std::stringstream name;
    if (k == 1) 
      name << "FEM_PK(" << n << ",1)";
    else {
      for (s = 2; s <= k; ++s) if ((k % s) == 0) break;
      name << "FEM_GEN_HIERARCHICAL(FEM_PK_HIERARCHICAL(" << n << ","
	   << k/s << "), FEM_PK(" << n << "," << k << "))";
    }
    return fem_descriptor(name.str());
  }

  static pfem QK_hierarch_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    if (n <= 0 || n >= 100 || k <= 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");
    std::stringstream name;
    if (n == 1)
      name << "FEM_PK_HIERARCHICAL(1," << k << ")";
    else
      name << "FEM_PRODUCT(FEM_PK_HIERARCHICAL(" << n-1 << "," << k
	   << "),FEM_PK(1_HIERARCHICAL," << k << "))";
    return fem_descriptor(name.str());
  }

  static pfem PK_prism_hierarch_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    if (n <= 1 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");
    std::stringstream name;
    if (n == 2)
      name << "FEM_QK_HIERARCHICAL(1," << k << ")";
    else 
      name << "FEM_PRODUCT(FEM_PK_HIERARCHICAL(" << n-1 << "," << k
	   << "),FEM_PK_HIERARCHICAL(1," << k << "))";
    return fem_descriptor(name.str());
  }


  /* ******************************************************************** */
  /* parallelepiped fems.                                                 */
  /* ******************************************************************** */

  static pfem QK_fem_(fem_param_list &params, bool discontinuous) {
    const char *fempk = discontinuous ? "FEM_PK_DISCONTINUOUS" : "FEM_PK";
    const char *femqk = discontinuous ? "FEM_QK_DISCONTINUOUS" : "FEM_QK";
    if (!(params.size() == 2 || (discontinuous && params.size() == 3)))
      DAL_THROW(failure_error, 
		"Bad number of parameters : " << params.size() << " should be 2.");
    if ((params[0].type() != 0 || params[1].type() != 0) ||
	params.size() == 3 && params[2].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");

    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    char alpha[128]; alpha[0] = 0;
    if (discontinuous && params.size() == 3) {
      scalar_type v = params[2].num();
      if (v < 0 || v > 1) DAL_THROW(failure_error, "Bad value for alpha: " << v);
      sprintf(alpha, ",%g", v);
    }
    if (n <= 0 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");
    std::stringstream name;
    if (n == 1)
      name << fempk << "(1," << k << alpha << ")";
    else 
      name << "FEM_PRODUCT(" << femqk << "(" << n-1 << ","
	   << k << alpha << ")," << fempk << "(1," << k << alpha << "))";
    return fem_descriptor(name.str());
  }
  static pfem QK_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &) {
    return QK_fem_(params, false);
  }
  static pfem QK_discontinuous_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &) {
    return QK_fem_(params, true);
  }

  
  /* ******************************************************************** */
  /* prims fems.                                                          */
  /* ******************************************************************** */

  static pfem PK_prism_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    if (n <= 1 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");
    std::stringstream name;
    if (n == 2)
      name << "FEM_QK(1," << k << ")";
    else 
      name << "FEM_PRODUCT(FEM_PK(" << n-1 << "," << k << "),FEM_PK(1,"
	   << k << "))";
    return fem_descriptor(name.str());
  }

  static pfem PK_prism_discontinuous_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &) {
    if (params.size() != 2 && params.size() != 3)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0 ||
	(params.size() == 3 && params[2].type() != 0))
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    char alpha[128]; alpha[0] = 0;
    if (params.size() == 3) {
      scalar_type v = params[2].num();
      if (v < 0 || v > 1) DAL_THROW(failure_error, "Bad value for alpha: " << v);
      sprintf(alpha, ",%g", v);
    }
    if (n <= 1 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");
    std::stringstream name;
    if (n == 2)
      name << "FEM_QK_DISCONTINUOUS(1," << k << alpha << ")";
    else 
      name << "FEM_PRODUCT(FEM_PK_DISCONTINUOUS(" << n-1 << "," << k << alpha
	   << "),FEM_PK_DISCONTINUOUS(1,"
	   << k << alpha << "))";
    return fem_descriptor(name.str());
  }

  /* ******************************************************************** */
  /*	P1 NON CONFORMING (dim 2)                                         */
  /* ******************************************************************** */

   static pfem P1_nonconforming_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {
    if (params.size() != 0)
      DAL_THROW(failure_error, "Bad number of parameters");
    fem<base_poly> *p = new fem<base_poly>;
    p->mref_convex() = bgeot::simplex_of_reference(2);
    p->dim() = 2;
    p->is_equivalent() = p->is_polynomial() = p->is_lagrange() = true;
    p->estimated_degree() = 1;
    p->init_cvs_node();
    base_poly one(2, 0), x(2, 1, 0), y(2, 1, 1); one.one();
    p->base().resize(3);

    p->add_node(lagrange_dof(2), base_small_vector(0.5, 0.5));
    p->base()[0] = x * 2.0 + y * 2.0 - one;
    p->add_node(lagrange_dof(2), base_small_vector(0.0, 0.5));
    p->base()[1] = one - x * 2.0;
    p->add_node(lagrange_dof(2), base_small_vector(0.5, 0.0));
    p->base()[2] = one - y * 2.0;
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));

    return p;
  }

   /* ******************************************************************** */
   /*	P1 element with a bubble base fonction on a face                   */
   /* ******************************************************************** */

   struct P1_wabbfoaf_ : public PK_fem_ {
     P1_wabbfoaf_(dim_type nc);
   };

  P1_wabbfoaf_::P1_wabbfoaf_(dim_type nc) : PK_fem_(nc, 1) {
    is_lag = false; es_degree = 2;
    base_node pt(nc); pt.fill(0.5);
    unfreeze_cvs_node();
    add_node(bubble1_dof(nc), pt);
    base_.resize(nb_dof(0));
    base_[nc+1] = base_[1]; base_[nc+1] *= scalar_type(1 << nc);
    for (int i = 2; i <= nc; ++i) base_[nc+1] *= base_[i];
    // Le raccord assure la continuite
    // des possibilités de raccord avec du P2 existent mais il faudrait
    // modifier qlq chose (transformer les fct de base P1) 
  }

  static pfem P1_with_bubble_on_a_face(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n <= 1 || n >= 100 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");
    virtual_fem *p = new P1_wabbfoaf_(n);
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*	P1 element with a bubble base fonction on a face : type lagrange  */
  /* ******************************************************************** */

  struct P1_wabbfoafla_ : public PK_fem_
  { // idem elt prec mais avec raccord lagrange. A faire en dim. quelconque ..
    P1_wabbfoafla_(void);
  };

  P1_wabbfoafla_::P1_wabbfoafla_(void) : PK_fem_(2, 1) {
    unfreeze_cvs_node();
    es_degree = 2;
    base_node pt(2); pt.fill(0.5);
    add_node(lagrange_dof(2), pt);
    base_.resize(nb_dof(0));
    
    base_poly one(2, 0); one.one();
    base_poly x(2, 1, 0), y(2, 1, 1);
    
    base_[0] = one - y - x;
    base_[1] = x * (one - y * 2.0);
    base_[2] = y * (one - x * 2.0);
    base_[3] = y * x * 4.0;
  }
  
  static pfem P1_with_bubble_on_a_face_lagrange(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {
    if (params.size() != 0)
      DAL_THROW(failure_error, "Bad number of parameters");
    virtual_fem *p = new P1_wabbfoafla_;
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }


  /* ******************************************************************** */
  /*	PK Gauss-Lobatto element on the segment                           */
  /* ******************************************************************** */

  class PK_GL_fem_ : public fem<base_poly> {
  public :
    PK_GL_fem_(unsigned k);
  };
  
  PK_GL_fem_::PK_GL_fem_(unsigned k) {
    cvr = bgeot::simplex_of_reference(1);
    dim_ = cvr->structure()->dim();
    is_equiv = is_pol = is_lag = true;
    es_degree = k;
    
    if (k >= fem_coeff_gausslob_max_k || !fem_coeff_gausslob[k]) 
      DAL_THROW(dal::failure_error, "try another degree");
    
    init_cvs_node();
    std::stringstream sstr; sstr << "IM_GAUSSLOBATTO1D(" << k*2-1 << ")";
    pintegration_method gl_im = int_method_descriptor(sstr.str());
    for (size_type i = 0; i < k+1; ++i) {
      add_node(lagrange_dof(1), gl_im->approx_method()->point(i));
    }
    base_.resize(k+1);
    const double *coefs = fem_coeff_gausslob[k];
    for (size_type r = 0; r < k+1; r++) {
      base_[r] = base_poly(1,k);
      std::copy(coefs, coefs+k+1, base_[r].begin());
      coefs += k+1;
      //cerr << "base_[" << r << "]=" << base_[r] << " @ " << gl_im->approx_method()->point(r) << "\n";
    }
  }

  static pfem PK_GL_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int k = int(::floor(params[0].num() + 0.01));
    virtual_fem *p = new PK_GL_fem_(k);
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*	Hermite element on the segment                                    */
  /* ******************************************************************** */

  struct hermite_segment__ : public fem<base_poly> {
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
			   bgeot::pgeometric_trans pgt) const;
    hermite_segment__(void);
  };

  void hermite_segment__::mat_trans(base_matrix &M,
					    const base_matrix &G,
					 bgeot::pgeometric_trans pgt) const { 
    dim_type P = 1, N = G.nrows();
    base_matrix K(N, P); // optimisable : eviter l'allocation.
    M.fill(1.0);
    bgeot::pgeotrans_precomp pgp = bgeot::geotrans_precomp(pgt, node_tab(0));
    if (N != 1)
      DAL_THROW(failure_error, "This element cannot be used for Q > 1");
    // gradient au pt 0
    gmm::mult(G, pgp->grad(0), K);
    M(2,2) = K(0,0);
    cout << "K(0,0) = " << K(0,0) << endl;
    // gradient au pt 1
    if (!(pgt->is_linear())) gmm::mult(G, pgp->grad(1), K);
    M(3,3) = K(0,0);
  }

  hermite_segment__::hermite_segment__(void) { 
    base_node pt(1);
    base_poly one(1, 0), x(1, 1, 0); one.one();
    cvr = bgeot::simplex_of_reference(1);
    dim_ = cvr->structure()->dim();
    init_cvs_node();
    es_degree = 3;
    is_pol = true;
    is_lag = false;
    // is_equiv = false;
    is_equiv = true;
    base_.resize(4);
    
    pt[0] = 0.0; add_node(lagrange_dof(1), pt);
    base_[0] = one + x * x * (-one * 3.0  + x * 2.0);
    // base_[0] = one - x;
    // base_[0] = (one - x) * (one - x * 2.0);
    pt[0] = 1.0; add_node(lagrange_dof(1), pt);
    base_[1] = x * x * (one * 3.0 - x * 2.0);
    // base_[1] = x * (x * 2.0 - one);
    pt[0] = 0.0; add_node(derivative_dof(1, 0), pt);
    // pt[0] = 0.1; add_node(lagrange_nonconforming_dof(1), pt);
    base_[2] = x * (one + x * (-one * 2.0 + x));
    pt[0] = 1.0; add_node(derivative_dof(1, 0), pt);
    // pt[0] = 0.9; add_node(lagrange_nonconforming_dof(1), pt);
    base_[3] = x * x * (-one + x);
    
    cout << "base(0) = " << base_[0] << endl;
    cout << "base(1) = " << base_[1] << endl;
    cout << "base(2) = " << base_[2] << endl;
    cout << "base(3) = " << base_[3] << endl;
  }

  static pfem segment_Hermite_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {
    if (params.size() != 0)
      DAL_THROW(failure_error, "Bad number of parameters");
    virtual_fem *p = new hermite_segment__;
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*	DISCONTINUOUS PK                                                  */
  /* ******************************************************************** */

  struct PK_discont_ : public PK_fem_ {
  public :
    
    PK_discont_(dim_type nc, short_type k, scalar_type alpha=0.)
      : PK_fem_(nc, k) {
      std::fill(dof_types_.begin(), dof_types_.end(),
		lagrange_nonconforming_dof(nc));

      if (alpha != 0.) {
	base_node G = 
	  dal::mean_value(cv_node.points().begin(), cv_node.points().end());
      for (size_type i=0; i < cv_node.nb_points(); ++i) 
	cv_node.points()[i] = (1-alpha)*cv_node.points()[i] + alpha*G;
      // cout << "creation of PK_discont(" << int(nc) << "," << k
      //      << "," << alpha << ")\n";
	for (size_type d = 0; d < nc; ++d) {
	  base_poly S(1,2); 
	  S[0] = -alpha * G[d] / (1-alpha);
	  S[1] = 1. / (1-alpha);
	  for (size_type j=0; j < nb_base(0); ++j) {
	    base_[j] = bgeot::poly_substitute_var(base_[j],S,d);
	  }
	}
      }
    }
  };
  
  static pfem PK_discontinuous_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {
    if (params.size() != 2 && params.size() != 3)
      DAL_THROW(failure_error, "Bad number of parameters : " << params.size()
		<< " should be 2 or 3.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    scalar_type alpha = 0.;
    if (params.size() == 3) alpha = params[2].num();
    if (n <= 0 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num() || alpha < 0 || alpha >= 1)
      DAL_THROW(failure_error, "Bad parameters");
    virtual_fem *p = new PK_discont_(n, k, alpha);
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }


  /* ******************************************************************** */
  /*	PK element with a bubble base fonction                            */
  /* ******************************************************************** */
  
  struct PK_with_cubic_bubble_ : public PK_fem_ {
    PK_with_cubic_bubble_(dim_type nc, short_type k);
  };
  
  PK_with_cubic_bubble_::PK_with_cubic_bubble_(dim_type nc, short_type k)
    : PK_fem_(nc, k) {
    unfreeze_cvs_node();
    is_lag = false; es_degree = nc+1;
    base_node pt(nc); 
    size_type j;
    PK_fem_ P1(nc, 1);
    
    pt.fill(1./(nc+1)); /* barycenter of the convex */
    
    add_node(bubble1_dof(nc), pt);
    base_.resize(nb_dof(0));
    
    j = nb_dof(0) - 1;
    base_[j] = base_poly(nc, 0);
    base_[j].one();
    for (size_type i = 0; i < P1.nb_dof(0); i++) base_[j] *= P1.base()[i];
    // cout << "buble = " << base_[j] << endl;
  }

  static pfem PK_with_cubic_bubble(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    if (k >= n+1) DAL_THROW(dimension_error, "dimensions mismatch");
    if (n <= 0 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");
    virtual_fem *p = new PK_with_cubic_bubble_(n, k);
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*	classical fem                                                     */
  /* ******************************************************************** */

  static pfem classical_fem_(const char *suffix, const char *arg, 
			     bgeot::pgeometric_trans pgt,
			     short_type k) {
    static bgeot::pgeometric_trans pgt_last = 0;
    static short_type k_last = short_type(-1);
    static pfem fm_last = 0;
    bool found = false;

    if (pgt_last == pgt && k_last == k)
      return fm_last;


    size_type n = pgt->structure()->dim();
    size_type nbp = pgt->basic_structure()->nb_points();
    std::stringstream name;

    /* Identifying P1-simplexes.                                          */

    if (nbp == n+1)
      if (pgt->basic_structure() == bgeot::simplex_structure(n))
    	{ name << "FEM_PK" << suffix << "("; found = true; }
    
    /* Identifying Q1-parallelepiped.                                     */

    if (!found && nbp == (size_type(1) << n))
      if (pgt->basic_structure() == bgeot::parallelepiped_structure(n))
    	{ name << "FEM_QK" << suffix << "("; found = true; }

    /* Identifying Q1-prisms.                                             */
 
    if (!found && nbp == 2 * n)
      if (pgt->basic_structure() == bgeot::prism_structure(n))
     	{ name << "FEM_PK_PRISM" << suffix << "("; found = true; }
     
    // To be completed

    if (found) {
      name << int(n) << ',' << int(k) << arg << ')';
      fm_last = fem_descriptor(name.str());
      pgt_last = pgt;
      k_last = k;
      return fm_last;
    }
 

    DAL_THROW(to_be_done_error,
	      "This element is not taken into account. Contact us");
  }

  pfem classical_fem(bgeot::pgeometric_trans pgt, short_type k) {
    return classical_fem_("", "", pgt, k);
  }
  
  pfem classical_discontinuous_fem(bgeot::pgeometric_trans pgt, short_type k, scalar_type alpha) {
    char arg[128]; arg[0] = 0;
    if (alpha) sprintf(arg, ",%g", alpha); 
    return classical_fem_("_DISCONTINUOUS", arg, pgt, k);
  }
  
  /* ******************************************************************** */
  /*    Naming system                                                     */
  /* ******************************************************************** */

  pfem structured_composite_fem_method(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies);
  pfem PK_composite_hierarch_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies);
  pfem PK_composite_full_hierarch_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies);

  struct fem_naming_system : public dal::naming_system<virtual_fem> {
    fem_naming_system() : dal::naming_system<virtual_fem>("FEM") {
      add_suffix("HERMITE_SEGMENT", segment_Hermite_fem);
      add_suffix("PK", PK_fem);
      add_suffix("QK", QK_fem);
      add_suffix("QK_DISCONTINUOUS", QK_discontinuous_fem);
      add_suffix("PK_PRISM", PK_prism_fem);
      add_suffix("PK_DISCONTINUOUS", PK_discontinuous_fem);
      add_suffix("PK_PRISM_DISCONTINUOUS", PK_prism_discontinuous_fem);
      add_suffix("PK_WITH_CUBIC_BUBBLE", PK_with_cubic_bubble);
      add_suffix("PRODUCT", product_fem);
      add_suffix("P1_NONCONFORMING", P1_nonconforming_fem);
      add_suffix("P1_BUBBLE_FACE", P1_with_bubble_on_a_face);
      add_suffix("P1_BUBBLE_FACE_LAG", P1_with_bubble_on_a_face_lagrange);
      add_suffix("GEN_HIERARCHICAL", gen_hierarchical_fem);
      add_suffix("PK_HIERARCHICAL", PK_hierarch_fem);
      add_suffix("QK_HIERARCHICAL", QK_hierarch_fem);
      add_suffix("PK_PRISM_HIERARCHICAL", PK_prism_hierarch_fem);
      add_suffix("STRUCTURED_COMPOSITE", structured_composite_fem_method);
      add_suffix("PK_HIERARCHICAL_COMPOSITE", PK_composite_hierarch_fem);
      add_suffix("PK_FULL_HIERARCHICAL_COMPOSITE",
		 PK_composite_full_hierarch_fem);
      add_suffix("PK_GAUSSLOBATTO1D", PK_GL_fem);
    }
  };
  
  pfem fem_descriptor(std::string name) {
    size_type i = 0;
    return dal::singleton<fem_naming_system>::instance().method(name, i);
  }

  std::string name_of_fem(pfem p) {
    try {
      return dal::singleton<fem_naming_system>::instance().shorter_name_of_method(p);
    } 
    catch (dal::failure_error&) {
      return std::string("FEM_UNKNOWN()");
    }
  }

  /* ******************************************************************** */
  /*    Aliases functions                                                 */
  /* ******************************************************************** */

  pfem PK_fem(size_type n, short_type k) {
    static pfem pf = 0;
    static size_type d = size_type(-2);
    static short_type r = short_type(-2);
    if (d != n || r != k) {
      std::stringstream name;
      name << "FEM_PK(" << n << "," << k << ")";
      pf = fem_descriptor(name.str());
      d = n; r = k;
    }
    return pf;
  }

  pfem QK_fem(size_type n, short_type k) {
    static pfem pf = 0;
    static size_type d = size_type(-2);
    static short_type r = short_type(-2);
    if (d != n || r != k) {
      std::stringstream name;
      name << "FEM_QK(" << n << "," << k << ")";
      pf = fem_descriptor(name.str());
      d = n; r = k;
    }
    return pf;
  }

  pfem PK_prism_fem(size_type n, short_type k) {
    static pfem pf = 0;
    static size_type d = size_type(-2);
    static short_type r = short_type(-2);
    if (d != n || r != k) {
      std::stringstream name;
      name << "FEM_PK_PRISM(" << n << "," << k << ")";
      pf = fem_descriptor(name.str());
      d = n; r = k;
    }
    return pf;
  }

  /* ********************************************************************* */
  /*       Precomputation on fem.                                          */
  /* ********************************************************************* */

  DAL_DOUBLE_KEY(pre_fem_key_, pfem, bgeot::pstored_point_tab);

  fem_precomp_::fem_precomp_(pfem pff, bgeot::pstored_point_tab ps) :
    pf(pff), pspt(ps) {
      for (size_type i = 0; i < pspt->size(); ++i)
	if ((*pspt)[i].size() != pf->dim())
	  DAL_THROW(dimension_error, "dimensions mismatch");
    }

  //  fem_precomp_::fem_precomp_() : pf(0), pspt(0) {}
  
  void fem_precomp_::init_val() const {
    c.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i) 
      pf->base_value((*pspt)[i], c[i]);
  }

  void fem_precomp_::init_grad() const {
    pc.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i)
      pf->grad_base_value((*pspt)[i], pc[i]);
  }

  void fem_precomp_::init_hess() const {
    hpc.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i)
      pf->hess_base_value((*pspt)[i], hpc[i]);
  }

  pfem_precomp fem_precomp(pfem pf, bgeot::pstored_point_tab pspt) {
    dal::pstatic_stored_object o
      = dal::search_stored_object(pre_fem_key_(pf, pspt));
    if (o) return dal::stored_cast<fem_precomp_>(o);
    pfem_precomp p = new fem_precomp_(pf, pspt);
    dal::add_stored_object(new pre_fem_key_(pf, pspt), p, pspt,
			   dal::AUTODELETE_STATIC_OBJECT);
    if (dal::exists_stored_object(pf)) dal::add_dependency(p, pf);
    return p;
    
  }

}  /* end of namespace getfem.                                            */
