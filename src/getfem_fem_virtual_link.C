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
  void mesh_fem_link_fem::add_dof_to_cv(size_type cv, size_type i) {
    std::deque<size_type> *p = &(cv_info_tab[cv].doftab);
    std::deque<size_type>::iterator it = p->begin(), ite = p->end();
    for (; it != ite; ++it) if (*it == i) return;
    p->push_back(i);
  }
    
  void mesh_fem_link_fem::compute(void) {
    if (!valid) DAL_THROW(internal_error, "internal error.");
    bgeot::geotrans_inv gti;
    dal::bit_vector nn = pmf2->convex_index();
    pintegration_method pim;
    bgeot::pgeometric_trans pgt;
    size_type cv, nbpt, maxgpt = 0;
    dal::dynamic_array<std::deque<size_type> > gauss_to_cv; // + index local
    
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
  void mesh_fem_link_fem::receipt(const MESH_FEM_DELETE &m) {
    if (m.ptr == (void *)(pmf1)) {
      if (pmf2 == 0 || &(pmf1->linked_mesh()) != &(pmf2->linked_mesh()))
	sup_sender(pmf1->linked_mesh().lmsg_sender());
      pmf1 = 0; valid = false;
    }
    if (m.ptr == (void *)(pmf2)) {
      if (pmf1 == 0 || &(pmf1->linked_mesh()) != &(pmf2->linked_mesh()))
	sup_sender(pmf2->linked_mesh().lmsg_sender());
      pmf2 = 0; valid = false;
    }
    
  }
  void mesh_fem_link_fem::receipt(const MESH_FEM_CHANGE &m) {
    if (m.ptr == (void *)(pmf1) || m.ptr == (void *)(pmf2))
      to_be_computed = true;
    if (m.ptr == (void *)(pmf1) && pmf2 != 0)
      pmf2->linked_mesh().lmsg_sender().send(MESH_FEM_CHANGE((void *)(pmf2)));
  }
  
  mesh_fem_link_fem::mesh_fem_link_fem(mesh_fem &mf1, mesh_fem &mf2)
    : pmf1(&mf1), pmf2(&mf2) { 
    to_be_computed = true; valid = true;
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

  mesh_fem_link_fem::~mesh_fem_link_fem() {
    if (pmf1 != 0) sup_sender(pmf1->linked_mesh().lmsg_sender());
    if ((pmf2 != 0) 
	&& (pmf1 == 0 || &(pmf1->linked_mesh()) != &(pmf2->linked_mesh())))
      sup_sender(pmf2->linked_mesh().lmsg_sender());
  }

#ifdef dksfjdsjsdlk

  class _virtual_link_fem : public virtual_fem
  {

  protected :
    mesh_fem_link_fem *pmflf;
    bgeot::papprox_integration pai;
    
    
  public :

    virtual void mat_trans(base_matrix &M, const base_matrix &G,
			   bgeot::pgeometric_trans pgt) const
    { 
      ...;
    }
    
    void interpolation(const base_node &x, const base_matrix &G, 
		       bgeot::pgeometric_trans pgt,
		       const base_vector coeff, base_node &val) const
      { DAL_THROW(internal_error, "You cannot interpolate this element."); }
    void interpolation_grad(const base_node &x, const base_matrix &G,
			    bgeot::pgeometric_trans pgt,
			    const base_vector coeff, base_matrix &val) const
      { DAL_THROW(internal_error, "You cannot interpolate this element."); }
    void base_value(const base_node &x, base_tensor &t) const {
      ...;
    }
    void grad_base_value(const base_node &x, base_tensor &t) const {
      ...;
    }
    void hess_base_value(const base_node &x, base_tensor &t) const {
      ...;
    }
    
    _virtual_link_fem(const mesh_fem_link_fem &a,
		      bgeot::papprox_integration b) : pmflf(&a), pai(b) {
      es_degree = 5;
    }

  };


  pfem virtual_link_fem(mesh_fem_link_fem, pintegration_method pim) {

    if (pim.is_ppi) DAL_THROW(std::invalid_argument,
	     "This element is only defined on approximated integration.");

    ...;
    // définir une structure qui réagit à la suppression de la structure 
    // mesh_fem_link_fem	      
    

  }




#endif



}  /* end of namespace getfem.                                            */
