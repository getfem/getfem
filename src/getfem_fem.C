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
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; version 2 of the License.                 */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software Foundation, */
/* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.         */
/*                                                                         */
/* *********************************************************************** */


#include <getfem_fem.h>
#include <dal_tree_sorted.h>
#include <dal_algobase.h>

namespace getfem
{
  void virtual_fem::interpolation(pfem_precomp pfp, size_type ii,
			     const base_matrix &G,
			     const base_vector coeff, base_node &val) const {
    // optimisable.   verifier et faire le vectoriel
    base_matrix M;
    if (val.size() != target_dim())
      throw dimension_error
	("virtual_fem::interpolation : dimensions mismatch");
    
    size_type R = nb_dof();

    if (!is_equivalent()) // utilité ?
    { if (M.nrows() != R || M.ncols() != R) M.resize(R, R); mat_trans(M, G); }

    val.fill(0.0);
    
    for (size_type j = 0; j < R; ++j) {
      scalar_type co = 0.0;
      if (is_equivalent())
	co = coeff[j];
      else
	for (size_type i = 0; i < R; ++i)
	  co += coeff[i] * M(i, j);
      
      for (size_type r = 0; r < target_dim(); ++r)
	val[r] += co * pfp->val(ii)[j + r*R];
    } 
  }


  void virtual_fem::interpolation_grad(pfem_precomp pfp, size_type ii,
			     const base_matrix &G,
			     const base_vector coeff, base_matrix &val) const {
    // optimisable !!   verifier et faire le vectoriel
 
    size_type R = nb_dof();
    base_matrix M;
    size_type P = structure()->dim();
    

    if (val.nrows() != target_dim() ||
	val.ncols() != P ||
	ii >= pfp->get_point_tab()->size() ||
	coeff.size() != R)
      throw dimension_error
	("virtual_fem::interpolation_grad: dimension mismatch");

    if (pfp->get_pfem() != this) {
      throw internal_error("virtual_fem::interpolation_grad: internal error");
    }

    base_tensor::const_iterator it = pfp->grad(ii).begin();

    if (!is_equivalent())
     { if (M.nrows() != R || M.ncols() != R) M.resize(R, R); mat_trans(M, G); }

     val.fill(0.0);
     
    for (size_type k = 0; k < P; ++k) {
      for (size_type r = 0; r < target_dim(); ++r) {
	for (size_type j = 0; j < R; ++j, ++it) {
	  
	  scalar_type co = 0.0;
	  if (is_equivalent())
	    co = coeff[j];
	  else
	    for (size_type i = 0; i < R; ++i)
	      co += coeff[i] * M(i, j);
	  val(r,k) += co * (*it);
	} 
      }
    }
  }


  /* ******************************************************************** */
  /*	Class for description of an interpolation dof.                    */
  /* ******************************************************************** */

  enum ddl_type { LAGRANGE, NORM_DERIVATIVE, DERIVATIVE, MEAN_VALUE, BUBBLE1, LAGRANGE_NONCONFORMING };

  struct dof_description
  {
    std::vector<ddl_type> ddl_desc;
    bool linkable;
    dim_type coord_index;

    dof_description(void)
    { linkable = true; coord_index = 0; }
  };

  struct __dof_description_comp
  {
    int operator()(const dof_description &m, const dof_description &n) const
    { 
      int nn;
      nn = dal::lexicographical_less<std::vector<ddl_type> >()
	(m.ddl_desc, n.ddl_desc);
      if (nn < 0) return -1; if (nn > 0) return 1;
      nn = int(m.linkable) - int(n.linkable);
      if (nn < 0) return -1; if (nn > 0) return 1;
      nn = int(m.coord_index) - int(m.coord_index);
      if (nn < 0) return -1; if (nn > 0) return 1;
      return 0;
    }
  };

  static dal::dynamic_tree_sorted<dof_description, __dof_description_comp>
                                              *_dof_d_tab;

  static inline void init_tab(void) // because of problem with initialization
  {                                 // in dynamic libraries.
    static bool initialized = false;
    if (!initialized)
    { 
      initialized = true;
      _dof_d_tab = new dal::dynamic_tree_sorted<dof_description,
	                                        __dof_description_comp>();
    }
  }

  pdof_description lagrange_dof(dim_type n)
  {
    init_tab();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), LAGRANGE);
    size_type i = _dof_d_tab->add_norepeat(l);
    return &((*_dof_d_tab)[i]);
  }

  pdof_description lagrange_nonconforming_dof(dim_type n)
  {
    init_tab();
    dof_description l; l.linkable = false;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), LAGRANGE);
    size_type i = _dof_d_tab->add_norepeat(l);
    return &((*_dof_d_tab)[i]);
  }

  pdof_description bubble1_dof(dim_type n)
  {
    init_tab();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), BUBBLE1);
    size_type i = (*_dof_d_tab).add_norepeat(l);
    return &((*_dof_d_tab)[i]);
  }

  pdof_description derivative_dof(dim_type n, dim_type num_der)
  {
    init_tab();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), LAGRANGE);
    l.ddl_desc[num_der] = DERIVATIVE;
    size_type ii = (*_dof_d_tab).add_norepeat(l);
    return &((*_dof_d_tab)[ii]);
  }

  pdof_description norm_derivative_dof(dim_type n)
  {
    init_tab();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), NORM_DERIVATIVE);
    size_type ii = (*_dof_d_tab).add_norepeat(l);
    return &((*_dof_d_tab)[ii]);
  }

  pdof_description mean_value_dof(dim_type n)
  {
    init_tab();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), MEAN_VALUE);
    size_type ii = (*_dof_d_tab).add_norepeat(l);
    return &((*_dof_d_tab)[ii]);
  }

  pdof_description change_coord_index_dof(pdof_description a, dim_type n)
  {
    init_tab();
    dof_description l = *a;
    a->coord_index = n;
    size_type ii = (*_dof_d_tab).add_norepeat(l);
    return &((*_dof_d_tab)[ii]);
  }

  pdof_description product_dof(pdof_description a, pdof_description b)
  {
    init_tab();
    size_type nb1 = a->ddl_desc.size(), nb2 = b->ddl_desc.size();
    dof_description l;
    l.linkable = a->linkable && b->linkable;
    l.coord_index = std::max(a->coord_index, b->coord_index); // logique ?
    l.ddl_desc.resize(nb1+nb2);
    std::copy(a->ddl_desc.begin(), a->ddl_desc.end(), l.ddl_desc.begin());
    std::copy(b->ddl_desc.begin(), b->ddl_desc.end(), l.ddl_desc.begin()+nb1);
    size_type ii = (*_dof_d_tab).add_norepeat(l);
    return &((*_dof_d_tab)[ii]);
  }

  int dof_description_compare(pdof_description a, pdof_description b)
  {
    int nn;
    if ((nn = int(a->coord_index) - int(b->coord_index)) != 0) return nn;
    if ((nn = int(a->linkable) - int(b->linkable)) != 0) return nn;
    std::vector<ddl_type>::const_iterator
      ita = a->ddl_desc.begin(), itae = a->ddl_desc.end(),
      itb = b->ddl_desc.begin(), itbe = b->ddl_desc.end();
    for (; ita != itae && itb != itbe; ++ita, ++itb)
    { if ((nn = int(*ita) - int (*itb)) != 0) return nn; }
    for (; ita != itae; ++ita) if (*ita != LAGRANGE) return 1;
    for (; itb != itbe; ++itb) if (*itb != LAGRANGE) return -1;
    return 0;
  }

  bool dof_linkable(pdof_description a)
  { return a->linkable; }

  bool dof_compatibility(pdof_description a, pdof_description b)
  { return (dof_description_compare(a, b) == 0 && dof_linkable(a)); }


  void virtual_fem::add_node(const pdof_description &d, const base_node &pt) {
    short_type nb = cv_node.nb_points();
    cv_node.points().resize(nb+1);
    cv_node.points()[nb] = pt;
    _dof_types.resize(nb+1);
    _dof_types[nb] = d;
    cvs_node.add_point_adaptative(nb, short_type(-1));
    for (short_type f = 0; f < cvs_node.nb_faces(); ++f)
      if (dal::abs(cvr->is_in_face(f, pt)) < 1.0E-5)
	cvs_node.add_point_adaptative(nb, f);
    pspt_valid = false;
  }

  void virtual_fem::init_cvs_node(void) {
    cvs_node.init_for_adaptative(cvr->structure());
    cv_node = bgeot::convex<base_node>(&cvs_node);
    pspt_valid = false;
  }

  /* ******************************************************************** */
  /*	PK class.                                                         */
  /* ******************************************************************** */

  struct _PK_femi_light
  {
    dim_type nc; short_type K;
    bool operator < (const _PK_femi_light &l) const
    {
      if (nc < l.nc) return true; if (nc > l.nc) return false; 
      if (K < l.K) return true; return false;
    }
    _PK_femi_light(dim_type n, short_type k) { nc = n; K = k; }
    _PK_femi_light(void) { }
  };

  class _PK_fem : public fem<base_poly>
  {
    public :

      void calc_base_func(base_poly &p, size_type i, short_type K) const
      {
	dim_type N = dim();
	base_poly l0(N, 0), l1(N, 0);
	bgeot::power_index w(N+1);
	l0.one(); l1.one(); p = l0;

	if (K != 0)
	{
	  for (int nn = 0; nn < N; ++nn) l0 -= base_poly(N, 1, nn);
	  
	  w[0] = K;
	  for (int nn = 1; nn <= N; ++nn)
	  { 
	    w[nn]=int(floor(0.5+((cv_node.points()[i])[nn-1]*double(K))));
	    w[0]-=w[nn];
	  }
	
	  for (int nn = 0; nn <= N; ++nn)
	    for (int j = 0; j < w[nn]; ++j)
	      if (nn == 0)
		p *= (l0 * (scalar_type(K) / scalar_type(j+1))) 
		  - (l1 * (scalar_type(j) / scalar_type(j+1)));
	      else
		p *= (base_poly(N,1,nn-1) * (scalar_type(K)/scalar_type(j+1))) 
		  - (l1 * (scalar_type(j) / scalar_type(j+1)));
	}
      }

     _PK_fem(const _PK_femi_light &ls)
     {
       cvr = bgeot::simplex_of_reference(ls.nc);
       is_equiv = is_pol = is_lag = true;
       es_degree = ls.K;
       
       init_cvs_node();
       bgeot::pconvex_ref cvn = bgeot::simplex_of_reference(ls.nc, ls.K);
       size_type R = cvn->nb_points();
       for (size_type i = 0; i < R; ++i)
	 add_node(lagrange_dof(ls.nc), cvn->points()[i]);

       _base.resize(R);
       for (size_type r = 0; r < R; r++) calc_base_func(_base[r], r, ls.K);
     }
  };

  ppolyfem PK_fem(dim_type n, short_type k)
  {
    static dal::FONC_TABLE<_PK_femi_light, _PK_fem> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_PK_femi_light, _PK_fem>();
      isinit = true;
    }
    return tab->add(_PK_femi_light(n, k));
  }

  /* ******************************************************************** */
  /*	Tensorial product of fem (for polynomial fem).                    */
  /* ******************************************************************** */

  struct _pr_fem_int_light 
  {
    ppolyfem fi1, fi2;
    bool operator < (const _pr_fem_int_light &ls) const
    {
      if (fi1 < ls.fi1) return true; if (fi1 > ls.fi1) return false; 
      if (fi2 < ls.fi2) return true; return false;
    }
    _pr_fem_int_light(ppolyfem a, ppolyfem b)
    { fi1 = a; fi2 = b; }
    _pr_fem_int_light(void) { }
  };

  struct tproduct_femi : public fem<base_poly>
  { 
    tproduct_femi(const _pr_fem_int_light &ls)
    {
      ppolyfem fi1 = ls.fi1, fi2 = ls.fi2;
      if (fi2->target_dim() != 1) std::swap(fi1, fi2);
      if (fi2->target_dim() != 1) 
	DAL_THROW(dimension_error, "dimensions mismatch");
    
      is_pol = true;
      is_equiv = fi1->is_equivalent() && fi2->is_equivalent();
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
	  
      for (j = 0, r = 0; j < fi2->nb_components(); j++)
	for (i = 0; i < fi1->nb_components(); i++, ++r)
	  {
	    _base[r] = fi1->base()[i];
	    _base[r].direct_product(fi2->base()[j]); 
	  }
    }
  };

  ppolyfem product_fem(ppolyfem fex, ppolyfem fey) {
    static dal::FONC_TABLE<_pr_fem_int_light, tproduct_femi> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_pr_fem_int_light, tproduct_femi>();
      isinit = true;
    }
    return tab->add(_pr_fem_int_light(fex, fey));
  }

  /* ******************************************************************** */
  /* parallelepiped structures.                                           */
  /* ******************************************************************** */

  struct _QK_femi
  {
    dim_type nc; short_type K; ppolyfem pf;
    bool operator < (const _QK_femi &l) const
    {
      if (nc < l.nc) return true; if (nc > l.nc) return false; 
      if (K < l.K) return true; return false;
    }
    _QK_femi(dim_type n, short_type k) { nc = n; K = k; }
    _QK_femi(void) { }
  };

  ppolyfem QK_fem(dim_type n, short_type k)
  {
    static dal::dynamic_tree_sorted<_QK_femi> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::dynamic_tree_sorted<_QK_femi>();
      isinit = true;
    }
    if (n <= 1) return PK_fem(n, k);
    _QK_femi fi(n, k);
    size_type i = tab->search(fi);
    if (i == size_type(-1))
    {
      fi.pf = product_fem(QK_fem(n-1, k), PK_fem(1, k));
      i = tab->add(fi);
    }
    return (*tab)[i].pf;
  }



  /* ******************************************************************** */
  /*	P1 NON CONFORMING (dim 2)                                         */
  /* ******************************************************************** */

  struct _P1_nonconforming_fem : public fem<base_poly>
  {
    _P1_nonconforming_fem(void)
    {
      is_equiv = is_pol = is_lag = true;
      es_degree = 1;
      cvr = bgeot::simplex_of_reference(2);
      _base.resize(3);
      _base[0] = base_poly(2, 0);
      _base[0].one(); _base[1] = _base[2] = _base[0];
      init_cvs_node();
      base_node pt(2);
      pt[0] = pt[1] = 0.5;
      add_node(lagrange_dof(2), pt);
      _base[0] *= -1.0;
      _base[0] += (base_poly(2, 1, 0)*2.0) + (base_poly(2, 1, 1)*2.0);
      pt[0] = 0.0; pt[1] = 0.5;
      add_node(lagrange_dof(2), pt);
      _base[1] -= base_poly(2, 1, 0) * 2.0;
      pt[0] = 0.5; pt[1] = 0.0;
      add_node(lagrange_dof(2), pt);
      _base[2] -= base_poly(2, 1, 1) * 2.0;
    }
  };

  ppolyfem P1_nonconforming_fem(void) {
    static _P1_nonconforming_fem *p;
    static bool isinit = false;
    if (!isinit) {
      p = new _P1_nonconforming_fem();
      isinit = true;
    }
    return p;
  }

   /* ******************************************************************** */
   /*	P1 element with a bubble base fonction on a face                   */
   /* ******************************************************************** */

   struct _P1_wabbfoaf : public _PK_fem
   {
     public :
    
     _P1_wabbfoaf(dim_type nc) : _PK_fem(_PK_femi_light(nc, 1))
     {
       is_lag = false; es_degree = 2;
       base_node pt(nc); pt.fill(0.5);
       add_node(bubble1_dof(nc), pt);
       _base.resize(nb_dof());
       _base[nc+1] = _base[1]; _base[nc+1] *= scalar_type(1 << nc);
       for (int i = 2; i <= nc; ++i) _base[nc+1] *= _base[i];
       // Le raccord assure la continuite
       // des possibilités de raccord avec du P2 existent mais il faudrait
       // modifier qlq chose (transformer les fct de base P1) 
     }
   };  
  
  ppolyfem P1_with_bubble_on_a_face(dim_type nc)
  { 
    static dal::dynamic_array<ppolyfem, 2> *tab;
    static dal::bit_vector *exists;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::dynamic_array<ppolyfem, 2>();
      exists = new dal::bit_vector();
      isinit = true;
    }
    if (nc <= 1)
      throw dimension_error
	("P1_with_bubble_on_a_face : dimensions mismatch");
    if (!(*exists)[nc])
      { (*exists)[nc] = true; (*tab)[nc]= new _P1_wabbfoaf(nc); }
    return (*tab)[nc];
  }


  /* ******************************************************************** */
  /*	P1 element with a bubble base fonction on a face : type lagrange  */
  /* ******************************************************************** */

  struct _P1_wabbfoafla : public _PK_fem
  { // idem elt prec mais avec raccord lagrange. A faire en dim. quelconque ..
    _P1_wabbfoafla(void) : _PK_fem(_PK_femi_light(2, 1))
    {
      es_degree = 2;
      base_node pt(2); pt.fill(0.5);
      add_node(lagrange_dof(2), pt);
      _base.resize(nb_dof());
      
      base_poly one(2, 0); one.one();
      
      _base[0] = one - base_poly(2, 1, 1) - base_poly(2, 1, 0);
      _base[1] = base_poly(2, 1, 0) * (one - base_poly(2, 1, 1) * 2.0);
      _base[2] = base_poly(2, 1, 1) * (one - base_poly(2, 1, 0) * 2.0);
      _base[3] = base_poly(2, 1, 1) * base_poly(2, 1, 0) * 4.0;
    }
  };
  

  ppolyfem P1_with_bubble_on_a_face_lagrange(void)
  { 
    static bool _P1_wabbfoafla_exists = false;
    static _P1_wabbfoafla *elt;
    if (!_P1_wabbfoafla_exists)
    { _P1_wabbfoafla_exists = true; elt = new _P1_wabbfoafla; }
    return elt;
  }



  /* ******************************************************************** */
  /*	DISCONTINUOUS PK                                                  */
  /* ******************************************************************** */

  
  struct _PK_discont : public _PK_fem
   {
     public :
    
     _PK_discont(const _PK_femi_light &l) : _PK_fem(l)
     {
       std::fill(_dof_types.begin(), _dof_types.end(), lagrange_nonconforming_dof(l.nc));
     }
   };  

  ppolyfem PK_discontinuous_fem(dim_type n, short_type k)
  {
    static dal::FONC_TABLE<_PK_femi_light, _PK_discont> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_PK_femi_light, _PK_discont>();
      isinit = true;
    }
    return tab->add(_PK_femi_light(n, k));
  }


  /* ******************************************************************** */
  /*	PK element with a bubble base fonction                            */
  /* ******************************************************************** */
  
   struct _PK_with_cubic_bubble : public _PK_fem
   {
     public :
    
     _PK_with_cubic_bubble(const _PK_femi_light &l) : _PK_fem(l)
     {
       is_lag = false; es_degree = l.nc+1;
       base_node pt(l.nc); 
       int i,j;
       _PK_fem P1(_PK_femi_light(l.nc, 1));

       /* barycenter of the convex */
       pt.fill(1./(l.nc+1));
       /*
       for (i=0; i < cv_node.nb_points(); i++) {
	 pt += cv_node.points()[i] * (1./(float)(cv_node.nb_points()));
       }
       */

       add_node(bubble1_dof(l.nc), pt);
       _base.resize(nb_dof());


       j = nb_dof()-1;
       _base[j] = base_poly(l.nc, 0);
       _base[j].one();
       for (i=0; i < P1.nb_dof(); i++) {
	 _base[j] *= P1.base()[i];
       }
     }
   };

  ppolyfem PK_with_cubic_bubble_fem(dim_type n, short_type k)
  {
    static dal::FONC_TABLE<_PK_femi_light, _PK_with_cubic_bubble> *tab;
    static bool isinit = false;

    if (k >= n+1)
      throw dimension_error
	("PK_with_cubic_bubble_fem : dimensions mismatch");
    
    if (!isinit) {
      tab = new dal::FONC_TABLE<_PK_femi_light, _PK_with_cubic_bubble>();
      isinit = true;
    }
    return tab->add(_PK_femi_light(n, k));
  }



  /* ******************************************************************** */
  /*	classical fem                                                     */
  /* ******************************************************************** */

  pfem classical_fem(bgeot::pgeometric_trans pgt, short_type k)
  {
    size_type n = pgt->structure()->dim();
    size_type nbp = pgt->basic_structure()->nb_points();

    /* Identifying P1-simplexes.                                          */

    if (nbp == n+1)
      if (pgt->basic_structure() == bgeot::simplex_structure(n))
    	return PK_fem(n, k);
    
    /* Identifying Q1-parallelepiped.                                     */

    if (nbp == (1 << n))
      if (pgt->basic_structure() == bgeot::parallelepiped_structure(n))
    	return QK_fem(n, k);

    /* Identifying Q1-prisms.                                             */
 
    if (nbp == 2 * n)
      if (pgt->basic_structure() == bgeot::prism_structure(n))
     	return PK_prism_fem(n, k);
     
    // To be completed

 
    throw to_be_done_error
      ("classical_fem : This element is not taken into account. Contact us");
    return NULL;
  }

}  /* end of namespace getfem.                                            */
