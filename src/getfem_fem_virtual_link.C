/* *********************************************************************** */
/*                                                                         */
/* Library : GEneric Tool for Finite Element Methods (getfem)              */
/* File    : getfem_fem_virtual_link.C : definition of the virtual finite  */
/*           element method which interpole an other fem on another mesh.  */
/*                                                                         */
/* Date : Aout 9, 2002.                                                    */
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
#include <getfem_mesh_fem.h>
#include <bgeot_geotrans_inv.h>

namespace getfem
{

  struct mesh_fem_link_fem_light {
    mesh_fem *pmf1, *pmf2;
    bool operator < (const mesh_fem_link_fem_light &l) const {
      if (pmf1 < l.pmf1) return true; if (pmf1 > l.pmf1) return false; 
      if (pmf2 < l.pmf2) return true; return false;
    }
    mesh_fem_link_fem_light(mesh_fem *p1, mesh_fem *p2) : pmf1(p1), pmf2(p2) {}
    mesh_fem_link_fem_light(void) : pmf1(0), pmf2(0) {}
  };

  class mesh_fem_link_fem : public getfem_mesh_receiver {
 
  protected :

    struct gauss_pt_info {
      base_node localcoords;
      size_type indcv;
      gauss_pt_info(void) { indcv = size_type(-1); }
    };

    struct cv_info {
      std::vector<size_type> indgausstab;
      std::deque<size_type> doftab;
    };

    mesh_fem *pmf1, *pmf2; // pmf1 -> mef to interpolated.
                           // pmf2 -> mef to build.
    bool to_be_computed;
    std::vector<gauss_pt_info> gauss_ptab;
    std::vector<cv_info> cv_info_tab;
    size_type max_dof;

    void add_dof_to_cv(size_type cv, size_type i);
    
    void compute(void);

  public :
    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_SUP_CONVEX &m);
    void receipt(const MESH_SWAP_CONVEX &m);
    void receipt(const MESH_REFINE_CONVEX &m);
    void receipt(const MESH_UNREFINE_CONVEX &m);
    void receipt(const MESH_FEM_DELETE &m);
    void receipt(const MESH_FEM_CHANGE &m);

    size_type nb_max_dof_per_element(void);

    size_type ind_of_dof(size_type cv, size_type i) {
      if (to_be_computed) compute();
      if (i >= cv_info_tab[cv].doftab.size())
	return 0;
      else
	return cv_info_tab[cv].doftab[i];
    }

    void mat_trans(base_matrix &M, const base_matrix &G,
			   bgeot::pgeometric_trans pgt) const
    { 
      size_type npt = G.ncols(), P = G.nrows();
      base_node pt(P);
      std::vector<size_type> ind(npt);
      base_matrix::const_iterator itm = G.begin();
      for (size_type k = 0; k < npt; ++k, itm += P) {
	std::copy(itm, itm + P, pt.begin());
	ind[k] = pmf1->linked_mesh().points().search(P);
	if (ind[k] == size_type(-1))
	  DAL_THROW(internal_error, "internal error.");
      }
      // 1 - retrouver le num de convexe
      // 2 - interpoler les fct de base des dof du cv et porter les coeff dans
      //     la matrice M

      // ...

      DAL_THROW(to_be_done_error, "To be finished.");
    }

    mesh_fem_link_fem(const mesh_fem_link_fem_light &ls);
    ~mesh_fem_link_fem();

  };

  void mesh_fem_link_fem::add_dof_to_cv(size_type cv, size_type i) {
    std::deque<size_type> *p = &(cv_info_tab[cv].doftab);
    std::deque<size_type>::iterator it = p->begin(), ite = p->end();
    for (; it != ite; ++it) if (*it == i) return;
    p->push_back(i);
    max_dof = std::max(max_dof, p->size());
  }

  size_type mesh_fem_link_fem::nb_max_dof_per_element(void) {
    if (to_be_computed) compute();
    return max_dof;
  }
    
  void mesh_fem_link_fem::compute(void) {
    bgeot::geotrans_inv gti;
    dal::bit_vector nn = pmf2->convex_index();
    pintegration_method pim;
    bgeot::pgeometric_trans pgt;
    size_type cv, nbpt, maxgpt = 0;
    dal::dynamic_array<std::deque<size_type> > gauss_to_cv; // + index local
    max_dof = 0;
    
    cv_info_tab.resize(nn.last_true() + 1);
    std::fill(cv_info_tab.begin(), cv_info_tab.end(), cv_info());
    for (cv << nn; cv != size_type(-1); cv << nn) {
      pim = pmf2->int_method_of_element(cv);
      if (pim.is_ppi) 
	DAL_THROW(internal_error,
		  "Method only defined for approximate integration.");
      nbpt = pim.integration_points().size();
      pgt = pmf2->linked_mesh().trans_of_convex(cv);
      cv_info_tab[cv].indgausstab.resize(nbpt);
      for (size_type k = 0; k < nbpt; ++k) {
	size_type i = gti.add_point
	  (pgt->transform(pim.integration_points()[i],
			  pmf2->linked_mesh().points_of_convex(cv)));
	(gauss_to_cv[i]).push_back(cv);
	maxgpt = std::max(i+1, maxgpt);
	cv_info_tab[cv].indgausstab[k] = i;
      }
    }
    
    dal::dynamic_array<base_node> ptab;
    dal::dynamic_array<size_type> itab;
    nn = pmf1->convex_index();
    gauss_ptab.resize(maxgpt);
    std::fill(gauss_ptab.begin(), gauss_ptab.end(), gauss_pt_info());
    for (cv << nn; cv != size_type(-1); cv << nn) {
      nbpt = gti.points_in_convex(pmf1->linked_mesh().convex(cv),
				  pmf1->linked_mesh().trans_of_convex(cv),
				  ptab, itab);
      for (size_type k = 0; k < nbpt; ++k) {
	gauss_pt_info *p = &(gauss_ptab[itab[k]]);
	if (p->indcv == size_type(-1)) {
	  p->indcv = cv;
	  p->localcoords = ptab[k];
	  
	  for (size_type j = 0; j < pmf1->nb_dof_of_element(cv); ++j)
	    for (size_type l = 0; l < gauss_to_cv[itab[k]].size(); ++l)
	      add_dof_to_cv((gauss_to_cv[itab[k]])[l],
			    pmf1->ind_dof_of_element(cv)[j]);
	  }
      }
    }
    to_be_computed = false;
  }
 
  void mesh_fem_link_fem::receipt(const MESH_CLEAR &)
    { to_be_computed = true; }
  void mesh_fem_link_fem::receipt(const MESH_SUP_CONVEX &m)
    { to_be_computed = true; }
  void mesh_fem_link_fem::receipt(const MESH_SWAP_CONVEX &m)
    { to_be_computed = true; }
  void mesh_fem_link_fem::receipt(const MESH_REFINE_CONVEX &m) {
    // ajouter la strategie au rafinement / derafinement
    DAL_THROW(internal_error, "internal error");
  }
  void mesh_fem_link_fem::receipt(const MESH_UNREFINE_CONVEX &m){
    // ajouter la strategie au rafinement / derafinement
    DAL_THROW(internal_error, "internal error");
  }
  void mesh_fem_link_fem::receipt(const MESH_FEM_CHANGE &m) {
    if (m.ptr == (void *)(pmf1) || m.ptr == (void *)(pmf2))
      to_be_computed = true;
    if (m.ptr == (void *)(pmf1))
      pmf2->linked_mesh().lmsg_sender().send(MESH_FEM_CHANGE((void *)(pmf2)));
  }
  
  mesh_fem_link_fem::mesh_fem_link_fem(const mesh_fem_link_fem_light &ls)
    : pmf1(ls.pmf1), pmf2(ls.pmf2) { 
    to_be_computed = true;
    add_sender(pmf1->linked_mesh().lmsg_sender(), *this,
	       lmsg::mask(MESH_CLEAR()) | lmsg::mask(MESH_SUP_CONVEX()) |
	       lmsg::mask(MESH_SWAP_CONVEX()) |
	       lmsg::mask(MESH_REFINE_CONVEX()) |
	       lmsg::mask(MESH_UNREFINE_CONVEX()) |
	       lmsg::mask(MESH_FEM_CHANGE()) |
	       lmsg::mask(MESH_FEM_DELETE()));
    if (&(pmf1->linked_mesh()) != &(pmf2->linked_mesh()))
      add_sender(pmf2->linked_mesh().lmsg_sender(), *this,
		 lmsg::mask(MESH_CLEAR()) | lmsg::mask(MESH_SUP_CONVEX()) |
		 lmsg::mask(MESH_SWAP_CONVEX()) |
		 lmsg::mask(MESH_REFINE_CONVEX()) |
		 lmsg::mask(MESH_UNREFINE_CONVEX()) |
		 lmsg::mask(MESH_FEM_CHANGE()) |
		 lmsg::mask(MESH_FEM_DELETE()));
  }
  
  typedef mesh_fem_link_fem *pmesh_fem_link_fem;
  static void sup_virtual_link_fem(pmesh_fem_link_fem pmflf);

  mesh_fem_link_fem::~mesh_fem_link_fem() {
      sup_sender(pmf1->linked_mesh().lmsg_sender());
      if (&(pmf1->linked_mesh()) != &(pmf2->linked_mesh()))
	sup_sender(pmf2->linked_mesh().lmsg_sender());
      sup_virtual_link_fem(this);
  }
  
  static dal::FONC_TABLE<mesh_fem_link_fem_light, mesh_fem_link_fem> 
    *mflf_tab = 0;
  
  pmesh_fem_link_fem mf_link_fem(mesh_fem &mf1, mesh_fem &mf2) {
    if (mflf_tab == 0)
      mflf_tab = 
	new dal::FONC_TABLE<mesh_fem_link_fem_light, mesh_fem_link_fem>();
    return mflf_tab->add(mesh_fem_link_fem_light(&mf1, &mf2));
  }

  void sup_mf_link_fem(mesh_fem &mf1, mesh_fem &mf2) {
    if (mflf_tab != 0) {
      mflf_tab->sup(mesh_fem_link_fem_light(&mf1, &mf2));
    }
  }

  void mesh_fem_link_fem::receipt(const MESH_FEM_DELETE &m) {
    sup_mf_link_fem(*pmf1, *pmf2);
  }

  struct _virtual_link_fem_light {
    pmesh_fem_link_fem pmflf;
    bgeot::papprox_integration pai;
    bool operator < (const _virtual_link_fem_light &l) const {
      if (pmflf < l.pmflf) return true; if (pmflf > l.pmflf) return false; 
      if (pai < l.pai) return true; return false;
    }
    _virtual_link_fem_light(pmesh_fem_link_fem a,bgeot::papprox_integration b)
      : pmflf(a), pai(b) {}
    _virtual_link_fem_light(void) {}
  };


  class _virtual_link_fem : public virtual_fem
  {

  protected :
    pmesh_fem_link_fem pmflf;
    bgeot::papprox_integration pai;
    
    void build_dof(size_type nb) {
      init_cvs_node();
      dim_type d = pai->dim();
      base_node pt(d); pt.fill(0.0);
      for (size_type k = 0; k < nb; ++k)
	add_node(already_numerate_dof(d), pt);
    }
    
  public :

    pmesh_fem_link_fem associated_mf_link_fem(void) { return pmflf; }
    bgeot::papprox_integration associated_integration(void) { return pai; }
    
    virtual size_type nb_dof(void) const {
      size_type nb = pmflf->nb_max_dof_per_element();
      if (nb != _dof_types.size()) 
	((_virtual_link_fem *)(this))->build_dof(nb);
      return nb;
    }
    virtual size_type nb_base(void) const { return pai->nb_points(); }
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
			   bgeot::pgeometric_trans pgt) const
      { pmflf->mat_trans(M, G, pgt); }

    virtual size_type index_of_already_numerate_dof(size_type cv, size_type i)
      const { return pmflf->ind_of_dof(cv, i); }
    
    void interpolation(const base_node &x, const base_matrix &G, 
		       bgeot::pgeometric_trans pgt,
		       const base_vector coeff, base_node &val) const {
      DAL_THROW(internal_error,
	  "You cannot interpolate this element, use the original element.");
    }
    void interpolation_grad(const base_node &x, const base_matrix &G,
			    bgeot::pgeometric_trans pgt,
			    const base_vector coeff, base_matrix &val) const {
      DAL_THROW(internal_error,
	  "You cannot interpolate this element, use the original element.");
    }
    void base_value(const base_node &x, base_tensor &t) const {
      const bgeot::stored_point_tab *p = &(pai->integration_points());
      for (size_type i = 0; i < p->size(); ++i)
	if (&((*p)[i]) == &x) {
	  bgeot::multi_index mi(2);
	  mi[1] = target_dim(); mi[0] = nb_base();
	  t.adjust_sizes(mi);
	  std::fill(t.begin(), t.end(), 0.0);
	  t[i] = 1.0;
	  return;
	}
      DAL_THROW(internal_error,
	  "You cannot interpolate this element, use the original element.");
    }
    void grad_base_value(const base_node &x, base_tensor &t) const {
      DAL_THROW(to_be_done_error,
      "Sorry, for the moment, the gradient is not available on this element.");
    }
    void hess_base_value(const base_node &x, base_tensor &t) const {
      DAL_THROW(to_be_done_error,
      "Sorry, for the moment, the Hessian is not available on this element.");
    }
    
    _virtual_link_fem(const _virtual_link_fem_light &ls)
      : pmflf(ls.pmflf), pai(ls.pai) {
      is_equiv = is_pol = is_lag = false; es_degree = 5;
      cvr = pai->ref_convex();
    }

  };


  typedef dal::FONC_TABLE<_virtual_link_fem_light, _virtual_link_fem>
    virtual_link_fem_table;
  static virtual_link_fem_table *__vlf_tab = 0;

  pfem virtual_link_fem(mesh_fem &mf1, mesh_fem &mf2,
			pintegration_method pim) {
    if (pim.is_ppi) DAL_THROW(std::invalid_argument,
	     "This element is only defined on approximated integration.");
    if (__vlf_tab == 0) __vlf_tab = new virtual_link_fem_table();
    return __vlf_tab->add(_virtual_link_fem_light(mf_link_fem(mf1, mf2),
						  pim.method.pai));
  }
	
  static void sup_virtual_link_fem(pmesh_fem_link_fem pmflf) {
    if (mflf_tab != 0) {
      virtual_link_fem_table::desc_table_type::const_iterator
	it = __vlf_tab->table().begin(), ite = __vlf_tab->table().end();
      for (; it != ite; ++it)
	if (*it != 0 && (*it)->associated_mf_link_fem() == pmflf)
	  __vlf_tab->sup(_virtual_link_fem_light
			 (pmflf, (*it)->associated_integration()));
    }
  }

}  /* end of namespace getfem.                                            */
