/*===========================================================================

 Copyright (C) 1999-2015 Yves Renard

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

===========================================================================*/


#include "getfem/getfem_mesh_fem_product.h"

namespace getfem {
    
  void fem_product::init() {
    
    GMM_ASSERT1(pfems[0]->target_dim() == 1, "To be done");
    GMM_ASSERT1(pfems[1]->target_dim() == 1,
		"The second finite element should be scalar");

    cvr = pfems[0]->ref_convex(cv);
    dim_ = cvr->structure()->dim();
    is_equiv = real_element_defined = true;
    is_polycomp = is_pol = is_lag = false;
    es_degree = 5; /* humm ... */
    ntarget_dim = 1;
    std::stringstream nm;
    nm << "FEM_PRODUCT(" << pfems[0]->debug_name() << ","
       << pfems[1]->debug_name() << "," << cv << ")";
    debug_name_ = nm.str();
    
    init_cvs_node();
    for (dal::bv_visitor i(enriched_dof1); !i.finished(); ++i) {
      for (size_type j = 0; j < pfems[1]->nb_dof(cv); ++j)
	add_node(xfem_dof(pfems[0]->dof_types()[i], xfem_index+j),
		 pfems[0]->node_of_dof(cv,i));
    }
  }

  void fem_product::base_value(const base_node &, 
				 base_tensor &) const
  { GMM_ASSERT1(false, "No base values, real only element."); }
  void fem_product::grad_base_value(const base_node &, 
				      base_tensor &) const
  { GMM_ASSERT1(false, "No base values, real only element."); }
  void fem_product::hess_base_value(const base_node &, 
			     base_tensor &) const
  { GMM_ASSERT1(false, "No base values, real only element."); }

  void fem_product::real_base_value(const fem_interpolation_context &c,
				      base_tensor &t, bool) const {
    bgeot::multi_index mi(2);
    mi[1] = target_dim(); mi[0] = short_type(nb_dof(0));
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin(), itf;

    fem_interpolation_context c0 = c;
    std::vector<base_tensor> val_e(2);
    for (size_type k = 0; k < 2; ++k) {
      if (c0.have_pfp()) {
	c0.set_pfp(fem_precomp(pfems[k], &c0.pfp()->get_point_tab(),
			       c0.pfp()));
      } else { c0.set_pf(pfems[k]); }
      c0.base_value(val_e[k]);
    }
    
    assert(target_dim() == 1);
    for (dal::bv_visitor i(enriched_dof1); !i.finished(); ++i) {
      itf = val_e[1].begin();
      scalar_type e = val_e[0][i];
      for (size_type j = 0; j < pfems[1]->nb_dof(cv); ++j)
	*it++ = *itf++ * e;
    }
    assert(it == t.end());
  }

  void fem_product::real_grad_base_value(const fem_interpolation_context &c,
					 base_tensor &t, bool) const {
    bgeot::multi_index mi(3);
    mi[2] = short_type(c.N()); mi[1] = target_dim();
    mi[0] = short_type(nb_dof(0));
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin();
    
    fem_interpolation_context c0 = c;
    std::vector<base_tensor> grad_e(2), val_e(2);
    for (size_type k = 0; k < 2; ++k) {
      if (c0.have_pfp()) {
	c0.set_pfp(fem_precomp(pfems[k], &c0.pfp()->get_point_tab(),
			       c0.pfp()));
      } else { c0.set_pf(pfems[k]); }
      c0.grad_base_value(grad_e[k]);
      c0.base_value(val_e[k]);
    }

    assert(target_dim() == 1);
    for (dim_type k = 0; k < c.N() ; ++k) {
      for (dal::bv_visitor i(enriched_dof1); !i.finished(); ++i) {
	size_type posg0 = k * pfems[0]->nb_dof(cv);
	size_type posg1 = k * pfems[1]->nb_dof(cv);
	for (size_type j = 0; j < pfems[1]->nb_dof(cv); ++j)
	  *it++ = grad_e[0][i + posg0] * val_e[1][j]
	    + grad_e[1][j + posg1] * val_e[0][i];
      }
    }
    assert(it == t.end());
  }
  
  void fem_product::real_hess_base_value(const fem_interpolation_context &c,
				  base_tensor &t, bool) const {
    bgeot::multi_index mi(4);
    mi[3] = mi[2] = short_type(c.N()); mi[1] = target_dim();
    mi[0] = short_type(nb_dof(0));
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin();
    
    fem_interpolation_context c0 = c;
    std::vector<base_tensor> hess_e(2), grad_e(2), val_e(2);
    for (size_type k = 0; k < 2; ++k) {
      if (c0.have_pfp()) {
	c0.set_pfp(fem_precomp(pfems[k], &c0.pfp()->get_point_tab(),
			       c0.pfp()));
      } else { c0.set_pf(pfems[k]); }
      c0.hess_base_value(hess_e[k]);
      c0.grad_base_value(grad_e[k]);
      c0.base_value(val_e[k]);
    }

    assert(target_dim() == 1);
    for (dim_type k0 = 0; k0 < c.N(); ++k0) {
      for (dim_type k1 = 0; k1 < c.N() ; ++k1) {
	for (dal::bv_visitor i(enriched_dof1); !i.finished(); ++i) {
	  size_type posh0 = (k0*c.N()+k1) * pfems[0]->nb_dof(cv);
	  size_type posh1 = (k0*c.N()+k1) * pfems[1]->nb_dof(cv);
	  size_type posg00 = k0 * pfems[0]->nb_dof(cv);
	  size_type posg01 = k1 * pfems[0]->nb_dof(cv);
	  size_type posg10 = k0 * pfems[1]->nb_dof(cv);
	  size_type posg11 = k1 * pfems[1]->nb_dof(cv);
	  for (size_type j = 0; j < pfems[1]->nb_dof(cv); ++j) {
	    *it++ =  hess_e[0][i + posh0] * val_e[1][j] + 
	      hess_e[1][j + posh1] * val_e[0][i] + 
	      grad_e[0][i + posg00] * grad_e[1][j + posg11] +
	      grad_e[0][i + posg01] * grad_e[1][j + posg10];
	  }
	}
      }
    }
    assert(it == t.end());
  }

  void mesh_fem_product::clear_build_methods() {
    for (size_type i = 0; i < build_methods.size(); ++i)
      del_stored_object(build_methods[i]);
    build_methods.clear();
  }
  void mesh_fem_product::clear(void) {
    mesh_fem::clear();
    clear_build_methods();
    is_adapted = false;
  }
  
  DAL_SIMPLE_KEY(special_mflproduct_key, pfem);
  
  void mesh_fem_product::adapt(void) {
    context_check();
    clear();

    GMM_ASSERT1(!mf1.is_reduced() && !mf2.is_reduced(),
		"Sorry, mesh_fem_product not defined for reduced mesh_fems");

    for (dal::bv_visitor cv(linked_mesh().convex_index()); !cv.finished();
	 ++cv) {
      dal::bit_vector local_enriched_dof;
      for (size_type i = 0; i < mf1.nb_basic_dof_of_element(cv); ++i)
	if (enriched_dof.is_in(mf1.ind_basic_dof_of_element(cv)[i]))
	  local_enriched_dof.add(i);
      if (local_enriched_dof.card() > 0) {
	pfem pf = new fem_product(mf1.fem_of_element(cv),
				  mf2.fem_of_element(cv), cv,
				  xfem_index, local_enriched_dof);
	special_mflproduct_key *psm = new special_mflproduct_key(pf);
	dal::add_stored_object(psm, pf, pf->ref_convex(0), pf->node_tab(0));
	build_methods.push_back(pf);
	set_finite_element(cv, pf);
      }
    }
    is_adapted = true; touch();
  }


}  /* end of namespace getfem.                                             */

