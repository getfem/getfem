/* *********************************************************************** */
/*                                                                         */
/* Library : GEneric Tool for Finite Element Methods (getfem)              */
/* File    : getfem_fem.C : definition of the finite element methods       */
/*                                                                         */
/* Date : December 21, 1999.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 1999-2002  Yves Renard.                                   */
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


#include <getfem_fem.h>
#include <dal_tree_sorted.h>
#include <dal_algobase.h>
#include <ftool_naming.h>

namespace getfem
{
  typedef ftool::naming_system<virtual_fem>::param_list fem_param_list;

  void virtual_fem::interpolation(pfem_precomp pfp, size_type ii,
				  const base_matrix &G,
				  bgeot::pgeometric_trans pgt, 
				  const base_vector &coeff, base_node &val,
				  dim_type Qdim) const {
    // optimisable.   verifier et faire le vectoriel
    base_matrix M;
    size_type Qmult = size_type(Qdim) / target_dim();
    if (val.size() != Qdim)
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    size_type R = nb_dof(), RR = nb_base();

    if (!is_equivalent()) { M.resize(RR, R); mat_trans(M, G, pgt); }

    val.fill(0.0);
    for (size_type j = 0; j < RR; ++j) {
      for (size_type q = 0; q < Qmult; ++q) {
	scalar_type co = 0.0;
	if (is_equivalent())
	  co = coeff[j*Qmult+q];
	else
	  for (size_type i = 0; i < R; ++i)
	    co += coeff[i*Qmult+q] * M(i, j);
      
	for (size_type r = 0; r < target_dim(); ++r)
	  val[r*Qmult+q] += co * pfp->val(ii)[j + r*R];
      } 
    }
  }

  void virtual_fem::interpolation_grad(const base_node &x,
				       const base_matrix &G,
				       bgeot::pgeometric_trans pgt,
				       const base_vector &coeff,
				       base_matrix &val) const {
    dim_type N = G.nrows();
    dim_type P = dim();
    size_type npt = G.ncols();
    base_matrix pc(npt , P);
    base_matrix grad(N, P), B0(P,N), CS(P,P), M;
    base_matrix val2(target_dim(), P);
    base_poly PP;
    base_tensor t;
    
    for (size_type i = 0; i < npt; ++i)
      for (dim_type n = 0; n < P; ++n) {
	PP = pgt->poly_vector()[i];
	PP.derivative(n);
	pc(i, n) = PP.eval(x.begin());
      }
    
    gmm::mult(G, pc, grad);
    if (P != N) {
      gmm::mult(gmm::transposed(grad), grad, CS);
      gmm::lu_inverse(CS);
      gmm::mult(gmm::transposed(CS), gmm::transposed(grad), B0);
    }
    else {
      gmm::lu_inverse(grad);
      B0 = grad;
    }
     
    grad_base_value(x, t);
    base_tensor::iterator it = t.begin();
    
    size_type R = nb_dof(), RR = nb_base();
    
    if (!is_equivalent()) { M.resize(RR, R); mat_trans(M, G, pgt); }
    
    val2.fill(0.0);
    
    for (size_type k = 0; k < x.size(); ++k)
      for (size_type r = 0; r < target_dim(); ++r)
 	for (size_type j = 0; j < RR; ++j, ++it) {
 	  scalar_type co = 0.0;
 	  if (is_equivalent())
 	    co = coeff[j];
	  else
 	    for (size_type i = 0; i < R; ++i)
 	      co += coeff[i] * M(i, j);
	  
 	  val2(r,k) += co * (*it);
 	} 
    
    gmm::mult(val2, B0, val);
  }


  void virtual_fem::real_base_value(pgeotrans_precomp, pfem_precomp pfp,
				    size_type ip, const base_matrix &,
				    base_tensor &t, size_type) const
  { t = pfp->val(ip); }
  void virtual_fem::real_grad_base_value(pgeotrans_precomp,
					 pfem_precomp pfp,
					 size_type ip, const base_matrix &,
					 const base_matrix &B,
					 base_tensor &t, size_type) const
  { t.mat_transp_reduction(pfp->grad(ip), B, 2); }
  void virtual_fem::real_hess_base_value(pgeotrans_precomp,
					 pfem_precomp pfp,
					 size_type ip, const base_matrix &,
					 const base_matrix &B3,
					 const base_matrix &B32,
					 base_tensor &t, size_type) const {
    base_tensor tt = pfp->hess(ip);
    bgeot::multi_index mim(3);
    mim[2] = dal::sqr(tt.sizes()[2]); mim[1] = tt.sizes()[1];
    mim[0] = tt.sizes()[0];
    tt.adjust_sizes(mim);
    t.mat_transp_reduction(tt, B3, 2);
    tt.mat_transp_reduction(pfp->grad(ip), B32, 2);
    t -= tt;
  }


  /* ******************************************************************** */
  /*	Class for description of an interpolation dof.                    */
  /* ******************************************************************** */

  enum ddl_type { LAGRANGE, NORM_DERIVATIVE, DERIVATIVE, MEAN_VALUE, BUBBLE1, 
		  LAGRANGE_NONCONFORMING, ALREADY_NUMERATE};

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

    dof_description(void)
    { linkable = true; coord_index = FIRST; xfem_index = 0; }
  };

  struct __dof_description_comp {
    int operator()(const dof_description &m, const dof_description &n) const;
  };

  // ATTENTION : en cas de modif, changer aussi dof_description_compare,
  //             product_dof, et dof_hierarchical_compatibility.
  int __dof_description_comp::operator()(const dof_description &m,
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
    return 0;
  }

  static dal::dynamic_tree_sorted<dof_description, __dof_description_comp>
                                                            *_dof_d_tab = 0;

  static void init_tab(void) // because of problem with initialization
  {                          // in dynamic libraries.
    
    if (!_dof_d_tab)
      _dof_d_tab = new dal::dynamic_tree_sorted<dof_description,
	                                        __dof_description_comp>();
  }

  pdof_description lagrange_dof(dim_type n) {
    dim_type n_old = dim_type(-2);
    pdof_description p_old = 0;
    if (n != n_old) {
      init_tab();
      dof_description l;
      l.ddl_desc.resize(n);
      std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(LAGRANGE));
      size_type i = _dof_d_tab->add_norepeat(l);
      p_old = &((*_dof_d_tab)[i]);
      n_old = n;
    }
    return p_old;
  }

  pdof_description deg_hierarchical_dof(pdof_description p, int deg) {
    init_tab();
    dof_description l = *p;
    for (size_type i = 0; i < l.ddl_desc.size(); ++i)
      l.ddl_desc[i].hier_degree = deg;
    size_type i = _dof_d_tab->add_norepeat(l);
    return &((*_dof_d_tab)[i]);
  }

  pdof_description xfem_dof(pdof_description p, size_type ind) {
    init_tab(); dof_description l = *p; l.xfem_index = ind;
    size_type i = _dof_d_tab->add_norepeat(l);
    return &((*_dof_d_tab)[i]);
  }


  pdof_description to_coord_dof(pdof_description p, dim_type ct) {
    init_tab();
    dof_description l = *p;
    l.coord_index = coord_type(ct);
    size_type i = _dof_d_tab->add_norepeat(l);
    return &((*_dof_d_tab)[i]);
  }
  
  pdof_description raff_hierarchical_dof(pdof_description p, short_type deg) {
    init_tab();
    dof_description l = *p;
    for (size_type i = 0; i < l.ddl_desc.size(); ++i)
      l.ddl_desc[i].hier_raff = deg;
    size_type i = _dof_d_tab->add_norepeat(l);
    return &((*_dof_d_tab)[i]);
  }

  pdof_description lagrange_nonconforming_dof(dim_type n) {
    init_tab();
    dof_description l; l.linkable = false;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(LAGRANGE));
    size_type i = _dof_d_tab->add_norepeat(l);
    return &((*_dof_d_tab)[i]);
  }

  pdof_description bubble1_dof(dim_type n) {
    init_tab();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(BUBBLE1));
    size_type i = (*_dof_d_tab).add_norepeat(l);
    return &((*_dof_d_tab)[i]);
  }

  pdof_description derivative_dof(dim_type n, dim_type num_der) {
    init_tab();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(LAGRANGE));
    l.ddl_desc[num_der] = ddl_elem(DERIVATIVE);
    size_type ii = (*_dof_d_tab).add_norepeat(l);
    return &((*_dof_d_tab)[ii]);
  }

  pdof_description norm_derivative_dof(dim_type n) {
    init_tab();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(NORM_DERIVATIVE));
    size_type ii = (*_dof_d_tab).add_norepeat(l);
    return &((*_dof_d_tab)[ii]);
  }

  pdof_description mean_value_dof(dim_type n) {
    init_tab();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(MEAN_VALUE));
    size_type ii = (*_dof_d_tab).add_norepeat(l);
    return &((*_dof_d_tab)[ii]);
  }

  pdof_description already_numerate_dof(dim_type n) {
    init_tab();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(),ddl_elem(ALREADY_NUMERATE));
    size_type i = _dof_d_tab->add_norepeat(l);
    return &((*_dof_d_tab)[i]);
  }

  pdof_description product_dof(pdof_description a, pdof_description b) {
    init_tab();
    size_type nb1 = a->ddl_desc.size(), nb2 = b->ddl_desc.size();
    dof_description l;
    l.linkable = a->linkable && b->linkable;
    l.coord_index = std::max(a->coord_index, b->coord_index); // logique ?
    l.xfem_index = a->xfem_index;
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
    size_type ii = (*_dof_d_tab).add_norepeat(l);
    return &((*_dof_d_tab)[ii]);
  }

  // ATTENTION : en cas de modif, changer aussi
  //             __dof_description_comp::operator,
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
  { return (dof_description_compare(a, b) == 0 && dof_linkable(a)); }

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
    _dof_types.resize(nb+1);
    _dof_types[nb] = d;
    cvs_node.add_point_adaptative(nb, short_type(-1));
    for (short_type f = 0; f < cvs_node.nb_faces(); ++f)
      if (dal::abs(cvr->is_in_face(f, pt)) < 1.0E-7)
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

  class _PK_fem : public fem<base_poly> {
  public :
    void calc_base_func(base_poly &p, size_type i, short_type K) const;
    _PK_fem(dim_type nc, short_type k);
  };
  
  void _PK_fem::calc_base_func(base_poly &p, size_type i, short_type K) const {
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
  
  _PK_fem::_PK_fem(dim_type nc, short_type k) {
    cvr = bgeot::simplex_of_reference(nc);
    is_equiv = is_pol = is_lag = true;
    es_degree = k;
    
    init_cvs_node();
    bgeot::pconvex_ref cvn = bgeot::simplex_of_reference(nc, k);
    size_type R = cvn->nb_points();
    for (size_type i = 0; i < R; ++i)
      add_node(lagrange_dof(nc), cvn->points()[i]);
    
    _base.resize(R);
    for (size_type r = 0; r < R; r++) calc_base_func(_base[r], r, k);
  }

  static pfem PK_fem(fem_param_list &params) {
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
    return new _PK_fem(n, k);
  }


  /* ******************************************************************** */
  /*	Tensorial product of fem (for polynomial fem).                    */
  /* ******************************************************************** */

  struct tproduct_femi : public fem<base_poly>
  { 
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
      = bgeot::convex_direct_product(fi1->node_convex(), fi2->node_convex());
    cvr = bgeot::convex_ref_product(fi1->ref_convex(), fi2->ref_convex());
    init_cvs_node();
    
    ntarget_dim = fi2->target_dim();
    _base.resize(cv.nb_points() * ntarget_dim);
    size_type i, j, r;
    for (j = 0, r = 0; j < fi2->nb_dof(); ++j)
      for (i = 0; i < fi1->nb_dof(); ++i, ++r)
	add_node(product_dof(fi1->dof_types()[i], fi2->dof_types()[j]),
		 cv.points()[r]);
    
    for (j = 0, r = 0; j < fi2->nb_base_components(); j++)
      for (i = 0; i < fi1->nb_base_components(); i++, ++r) {
	_base[r] = fi1->base()[i];
	_base[r].direct_product(fi2->base()[j]); 
      }
  }

  static pfem product_fem(fem_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	  "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 1 || params[1].type() != 1)
      DAL_THROW(failure_error, "Bad type of parameters");
    pfem pf1 = params[0].method();
    pfem pf2 = params[1].method();
    if (!(pf1->is_polynomial() && pf2->is_polynomial()))
      DAL_THROW(failure_error, "Bad parameters");
    return new tproduct_femi(ppolyfem(pf1), ppolyfem(pf2));
  }

  /* ******************************************************************** */
  /*    Generic Hierarchical fem (for polynomial fem). To be interfaced.  */
  /* ******************************************************************** */

  struct thierach_femi : public fem<base_poly>
  { 
    thierach_femi(ppolyfem fi1, ppolyfem fi2);
  };


  thierach_femi::thierach_femi(ppolyfem fi1, ppolyfem fi2)
    : fem<base_poly>(*fi1) {
    if (fi2->target_dim() != fi1->target_dim())
      DAL_THROW(dimension_error, "dimensions mismatch.");
    if (fi2->basic_structure() != fi1->basic_structure())
      DAL_THROW(failure_error, "Incompatible elements.");
    if (!(fi1->is_equivalent() &&  fi2->is_equivalent()))
      DAL_THROW(to_be_done_error,
	   "Sorry, no hierachical construction for non tau-equivalent fem.");
    es_degree = fi2->estimated_degree();
    is_lag = false;
    unfreeze_cvs_node();
    for (size_type i = 0; i < fi2->nb_dof(); ++i) {
      bool found = false;
      for (size_type j = 0; j < fi1->nb_dof(); ++j) {
	if ( dal::lexicographical_less<base_node,
	     dal::approx_less<scalar_type> >()
	     (fi2->node_of_dof(i), fi1->node_of_dof(j)) == 0
	     && dof_hierarchical_compatibility(fi2->dof_types()[i],
					       fi1->dof_types()[j]))
	    { found = true; break; }
      }
      if (!found) {
	add_node(deg_hierarchical_dof(fi2->dof_types()[i], 
				      fi1->estimated_degree()),
		 fi2->node_of_dof(i));
	_base.resize(nb_dof());
	_base[nb_dof()-1] = (fi2->base())[i];
	// 	  cout << "adding base : " << _base[nb_dof()-1] << endl;
	// 	  cout << "point : " << node_of_dof(nb_dof()-1) << endl;
	// 	  cout << "Nb dof cici = " << nb_dof() << endl;
      }
    }
    //        cout << "Nb dof = " << nb_dof() << endl;
    //        for (size_type j = 0; j < nb_dof(); ++j)
    //  	cout << " base : " << j << " : " << _base[j] << endl;
  }

  struct thierach_femi_comp : public fem<polynomial_composite>
  { 
    thierach_femi_comp(ppolycompfem fi1, ppolycompfem fi2);  
  };

  thierach_femi_comp::thierach_femi_comp(ppolycompfem fi1, ppolycompfem fi2)
    : fem<polynomial_composite>(*fi1) {
    if (fi2->target_dim() != fi1->target_dim())
      DAL_THROW(dimension_error, "dimensions mismatch.");
    if (fi2->basic_structure() != fi1->basic_structure())
      DAL_THROW(failure_error, "Incompatible elements.");
    if (!(fi1->is_equivalent() && fi2->is_equivalent()))
      DAL_THROW(to_be_done_error,
	   "Sorry, no hierachical construction for non tau-equivalent fem.");
    es_degree = std::max(fi2->estimated_degree(), fi1->estimated_degree());
    
    is_lag = false;
    hier_raff = fi1->hierarchical_raff() + 1;
    unfreeze_cvs_node();
    for (size_type i = 0; i < fi2->nb_dof(); ++i) {
      bool found = false;
      for (size_type j = 0; j < fi1->nb_dof(); ++j) {
	if ( dal::lexicographical_less<base_node,
	     dal::approx_less<scalar_type> >()
	     (fi2->node_of_dof(i), fi1->node_of_dof(j)) == 0
	     && dof_hierarchical_compatibility(fi2->dof_types()[i],
					       fi1->dof_types()[j]))
	  { found = true; break; }
      }
      if (!found) {
	add_node(raff_hierarchical_dof(fi2->dof_types()[i], hier_raff),
		 fi2->node_of_dof(i));
	_base.resize(nb_dof());
	_base[nb_dof()-1] = (fi2->base())[i];
      }
    }
  }

  static pfem gen_hierarchical_fem(fem_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	  "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 1 || params[1].type() != 1)
      DAL_THROW(failure_error, "Bad type of parameters");
    pfem pf1 = params[0].method();
    pfem pf2 = params[1].method();

    if (pf1->is_polynomial() && pf2->is_polynomial())
      return new thierach_femi(ppolyfem(pf1),ppolyfem(pf2));
    if (pf1->is_polynomialcomp() && pf2->is_polynomialcomp())
     return new thierach_femi_comp(ppolycompfem(pf1),ppolycompfem(pf2));
    DAL_THROW(failure_error, "Bad parameters");
  }

  /* ******************************************************************** */
  /* PK hierarchical fem.                                                 */
  /* ******************************************************************** */

  static pfem PK_hierarch_fem(fem_param_list &params) {
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

  static pfem QK_hierarch_fem(fem_param_list &params) {
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

  static pfem PK_prism_hierarch_fem(fem_param_list &params) {
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
  /* parallelepiped structures.                                           */
  /* ******************************************************************** */

  static pfem QK_fem_(fem_param_list &params, const std::string& fempk, const std::string femqk) {
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
    std::stringstream name;
    if (n == 1)
      name << fempk << "(1," << k << ")";
    else 
      name << "FEM_PRODUCT(" << femqk << "(" << n-1 << "," << k << ")," << fempk << "(1,"
	   << k << "))";
    return fem_descriptor(name.str());
  }
  static pfem QK_fem(fem_param_list &params) {
    return QK_fem_(params, "FEM_PK", "FEM_QK");
  }
  static pfem QK_discontinuous_fem(fem_param_list &params) {
    return QK_fem_(params, "FEM_PK_DISCONTINUOUS", "FEM_QK_DISCONTINUOUS");
  }

  
  /* ******************************************************************** */
  /* prims structures.                                                    */
  /* ******************************************************************** */

  static pfem PK_prism_fem(fem_param_list &params) {
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

  /* ******************************************************************** */
  /*	P1 NON CONFORMING (dim 2)                                         */
  /* ******************************************************************** */

   static pfem P1_nonconforming_fem(fem_param_list &params) {
    if (params.size() != 0)
      DAL_THROW(failure_error, "Bad number of parameters");
    fem<base_poly> *p = new fem<base_poly>;
    p->ref_convex() = bgeot::simplex_of_reference(2);
    p->is_equivalent() = p->is_polynomial() = p->is_lagrange() = true;
    p->estimated_degree() = 1;
    p->init_cvs_node();
    base_poly one(2, 0), x(2, 1, 0), y(2, 1, 1); one.one();
    p->base().resize(3);

    p->add_node(lagrange_dof(2), base_vector(0.5, 0.5));
    p->base()[0] = x * 2.0 + y * 2.0 - one;
    p->add_node(lagrange_dof(2), base_vector(0.0, 0.5));
    p->base()[1] = one - x * 2.0;
    p->add_node(lagrange_dof(2), base_vector(0.5, 0.0));
    p->base()[2] = one - y * 2.0;

    return p;
  }

   /* ******************************************************************** */
   /*	P1 element with a bubble base fonction on a face                   */
   /* ******************************************************************** */

   struct _P1_wabbfoaf : public _PK_fem
   {
     _P1_wabbfoaf(dim_type nc);
   };

  _P1_wabbfoaf::_P1_wabbfoaf(dim_type nc) : _PK_fem(nc, 1) {
    is_lag = false; es_degree = 2;
    base_node pt(nc); pt.fill(0.5);
    unfreeze_cvs_node();
    add_node(bubble1_dof(nc), pt);
    _base.resize(nb_dof());
    _base[nc+1] = _base[1]; _base[nc+1] *= scalar_type(1 << nc);
    for (int i = 2; i <= nc; ++i) _base[nc+1] *= _base[i];
    // Le raccord assure la continuite
    // des possibilités de raccord avec du P2 existent mais il faudrait
    // modifier qlq chose (transformer les fct de base P1) 
  }

  static pfem P1_with_bubble_on_a_face(fem_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n <= 1 || n >= 100 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");
    return new _P1_wabbfoaf(n);
  }

  /* ******************************************************************** */
  /*	P1 element with a bubble base fonction on a face : type lagrange  */
  /* ******************************************************************** */

  struct _P1_wabbfoafla : public _PK_fem
  { // idem elt prec mais avec raccord lagrange. A faire en dim. quelconque ..
    _P1_wabbfoafla(void);
  };

  _P1_wabbfoafla::_P1_wabbfoafla(void) : _PK_fem(2, 1) {
    unfreeze_cvs_node();
    es_degree = 2;
    base_node pt(2); pt.fill(0.5);
    add_node(lagrange_dof(2), pt);
    _base.resize(nb_dof());
    
    base_poly one(2, 0); one.one();
    base_poly x(2, 1, 0), y(2, 1, 1);
    
    _base[0] = one - y - x;
    _base[1] = x * (one - y * 2.0);
    _base[2] = y * (one - x * 2.0);
    _base[3] = y * x * 4.0;
  }
  
  static pfem P1_with_bubble_on_a_face_lagrange(fem_param_list &params) {
    if (params.size() != 0)
      DAL_THROW(failure_error, "Bad number of parameters");
    return new _P1_wabbfoafla;
  }

  /* ******************************************************************** */
  /*	Hremite element on the segment                                    */
  /* ******************************************************************** */

  struct _hermite_segment_ : public fem<base_poly> {
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
			   bgeot::pgeometric_trans pgt) const;
    _hermite_segment_(void);
  };

  void _hermite_segment_::mat_trans(base_matrix &M,
					    const base_matrix &G,
					 bgeot::pgeometric_trans pgt) const { 
    dim_type P = 1, N = G.nrows();
    base_matrix K(N, P); // optimisable : eviter l'allocation.
    M.fill(1.0);
    pgeotrans_precomp pgp = geotrans_precomp(pgt, node_tab());
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

  _hermite_segment_::_hermite_segment_(void) { 
    base_node pt(1);
    base_poly one(1, 0), x(1, 1, 0); one.one();
    cvr = bgeot::simplex_of_reference(1);
    init_cvs_node();
    es_degree = 3;
    is_pol = true;
    is_lag = false;
    // is_equiv = false;
    is_equiv = true;
    _base.resize(4);
    
    pt[0] = 0.0; add_node(lagrange_dof(1), pt);
    _base[0] = one + x * x * (-one * 3.0  + x * 2.0);
    // _base[0] = one - x;
    // _base[0] = (one - x) * (one - x * 2.0);
    pt[0] = 1.0; add_node(lagrange_dof(1), pt);
    _base[1] = x * x * (one * 3.0 - x * 2.0);
    // _base[1] = x * (x * 2.0 - one);
    pt[0] = 0.0; add_node(derivative_dof(1, 0), pt);
    // pt[0] = 0.1; add_node(lagrange_nonconforming_dof(1), pt);
    _base[2] = x * (one + x * (-one * 2.0 + x));
    pt[0] = 1.0; add_node(derivative_dof(1, 0), pt);
    // pt[0] = 0.9; add_node(lagrange_nonconforming_dof(1), pt);
    _base[3] = x * x * (-one + x);
    
    cout << "base(0) = " << _base[0] << endl;
    cout << "base(1) = " << _base[1] << endl;
    cout << "base(2) = " << _base[2] << endl;
    cout << "base(3) = " << _base[3] << endl;
  }

  static pfem segment_Hermite_fem(fem_param_list &params) {
    if (params.size() != 0)
      DAL_THROW(failure_error, "Bad number of parameters");
    return new _hermite_segment_;
  }

  /* ******************************************************************** */
  /*	DISCONTINUOUS PK                                                  */
  /* ******************************************************************** */

  
  struct _PK_discont : public _PK_fem {
  public :
    
    _PK_discont(dim_type nc, short_type k, scalar_type alpha=0.) : _PK_fem(nc, k) {
      std::fill(_dof_types.begin(), _dof_types.end(),
		lagrange_nonconforming_dof(nc));
      base_node G = cv_node.points()[0];
      for (size_type i=0; i < cv_node.nb_points(); ++i) G += cv_node.points()[i];
      G /= scalar_type(cv_node.nb_points());
      for (size_type i=0; i < cv_node.nb_points(); ++i) 
	cv_node.points()[i] = (1-alpha)*cv_node.points()[i] + alpha*G;
      cout << "creation of PK_discont(" << nc << "," << k << "," << alpha << "\n";
      for (size_type d = 0; d < nc; ++d) {
        base_poly S(1,2); 
        S[0] = -alpha * G[d] / (1-alpha);
        S[1] = 1. / (1-alpha);
        for (size_type j=0; j < nb_base(); ++j) {
          _base[j] = bgeot::poly_substitute_var(_base[j],S,d);
        }
      }
      for (size_type j=0; j < nb_base(); ++j) cout << " base[" << j << "]=" << _base[j] << "\n";
      for (size_type i=0; i < cv_node.nb_points(); ++i) {
        cout << "node=" << cv_node.points()[i] << ":";
        for (size_type j=0; j < nb_base(); ++j) cout << _base[j].eval(cv_node.points()[i].begin()) << " ";
        cout << "\n";
      }
    }
  };
  
  static pfem PK_discontinuous_fem(fem_param_list &params) {
    if (params.size() != 2 && params.size() != 3)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2 or 3.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    scalar_type alpha = 0.;
    if (params.size() == 3) alpha = params[2].num();
    if (n <= 0 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num() || alpha < 0 || alpha >= 1)
      DAL_THROW(failure_error, "Bad parameters");
    return new _PK_discont(n, k, alpha);
  }


  /* ******************************************************************** */
  /*	PK element with a bubble base fonction                            */
  /* ******************************************************************** */
  
  struct _PK_with_cubic_bubble : public _PK_fem {
    _PK_with_cubic_bubble(dim_type nc, short_type k);
  };
  
  _PK_with_cubic_bubble::_PK_with_cubic_bubble(dim_type nc, short_type k)
    : _PK_fem(nc, k) {
    unfreeze_cvs_node();
    is_lag = false; es_degree = nc+1;
    base_node pt(nc); 
    size_type j;
    _PK_fem P1(nc, 1);
    
    pt.fill(1./(nc+1)); /* barycenter of the convex */
    
    add_node(bubble1_dof(nc), pt);
    _base.resize(nb_dof());
    
    j = nb_dof() - 1;
    _base[j] = base_poly(nc, 0);
    _base[j].one();
    for (size_type i = 0; i < P1.nb_dof(); i++) _base[j] *= P1.base()[i];
  }

  static pfem PK_with_cubic_bubble(fem_param_list &params) {
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
    return new _PK_with_cubic_bubble(n, k);
  }

  /* ******************************************************************** */
  /*	classical fem                                                     */
  /* ******************************************************************** */

  pfem classical_fem(bgeot::pgeometric_trans pgt, short_type k)
  {
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
    	{ name << "FEM_PK("; found = true; }
    
    /* Identifying Q1-parallelepiped.                                     */

    if (!found && nbp == (size_type(1) << n))
      if (pgt->basic_structure() == bgeot::parallelepiped_structure(n))
    	{ name << "FEM_QK("; found = true; }

    /* Identifying Q1-prisms.                                             */
 
    if (!found && nbp == 2 * n)
      if (pgt->basic_structure() == bgeot::prism_structure(n))
     	{ name << "FEM_PK_PRISM("; found = true; }
     
    // To be completed

    if (found) {
      name << int(n) << ',' << int(k) << ')';
      fm_last = fem_descriptor(name.str());
      pgt_last = pgt;
      k_last = k;
      return fm_last;
    }
 
    DAL_THROW(to_be_done_error,
	      "This element is not taken into account. Contact us");
  }

  
  /* ******************************************************************** */
  /*    Naming system                                                     */
  /* ******************************************************************** */

  pfem structured_composite_fem_method(fem_param_list &params);
  pfem PK_composite_hierarch_fem(fem_param_list &params);
  pfem PK_composite_full_hierarch_fem(fem_param_list &params);

  static ftool::naming_system<virtual_fem> *_fem_naming_system = 0;
  
  static void init_fem_naming_system(void) {
    _fem_naming_system = new ftool::naming_system<virtual_fem>("FEM");
    _fem_naming_system->add_suffix("HERMITE_SEGMENT", segment_Hermite_fem);
    _fem_naming_system->add_suffix("PK", PK_fem);
    _fem_naming_system->add_suffix("QK", QK_fem);
    _fem_naming_system->add_suffix("QK_DISCONTINUOUS", QK_discontinuous_fem);
    _fem_naming_system->add_suffix("PK_PRISM", PK_prism_fem);
    _fem_naming_system->add_suffix("PK_DISCONTINUOUS", PK_discontinuous_fem);
    _fem_naming_system->add_suffix("PK_WITH_CUBIC_BUBBLE",
				   PK_with_cubic_bubble);
    _fem_naming_system->add_suffix("PRODUCT", product_fem);
    _fem_naming_system->add_suffix("P1_NONCONFORMING", P1_nonconforming_fem);
    _fem_naming_system->add_suffix("P1_BUBBLE_FACE", P1_with_bubble_on_a_face);
    _fem_naming_system->add_suffix("P1_BUBBLE_FACE_LAG",
				   P1_with_bubble_on_a_face_lagrange);
    _fem_naming_system->add_suffix("GEN_HIERARCHICAL", gen_hierarchical_fem);
    _fem_naming_system->add_suffix("PK_HIERARCHICAL", PK_hierarch_fem);
    _fem_naming_system->add_suffix("QK_HIERARCHICAL", QK_hierarch_fem);
    _fem_naming_system->add_suffix("PK_PRISM_HIERARCHICAL",
				   PK_prism_hierarch_fem);
    _fem_naming_system->add_suffix("STRUCTURED_COMPOSITE",
				   structured_composite_fem_method);
    _fem_naming_system->add_suffix("PK_HIERARCHICAL_COMPOSITE",
				   PK_composite_hierarch_fem);
    _fem_naming_system->add_suffix("PK_FULL_HIERARCHICAL_COMPOSITE",
				   PK_composite_full_hierarch_fem);
  }
  
  pfem fem_descriptor(std::string name) {
    if (_fem_naming_system == 0) init_fem_naming_system();
    size_type i = 0;
    pfem res = _fem_naming_system->method(name, i);
    return res;
  }

  std::string name_of_fem(pfem p) {
    if (_fem_naming_system == 0) init_fem_naming_system();
    try {
      return _fem_naming_system->shorter_name_of_method(p);
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

}  /* end of namespace getfem.                                            */
