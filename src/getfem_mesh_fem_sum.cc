// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 1999-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================


#include <getfem_mesh_fem_sum.h>

namespace getfem {
    
  void fem_sum::init() {
    cvr = pfems[0]->ref_convex(cv);
    dim_ = cvr->structure()->dim();
    is_equiv = !smart_global_dof_linking_;
    real_element_defined = true;
    is_polycomp = is_pol = is_lag = false;
    es_degree = 5; /* humm ... */
    ntarget_dim = 1;

    std::stringstream nm;
    nm << "FEM_SUM(" << pfems[0]->debug_name() << ",";
    for (size_type i = 1; i < pfems.size(); ++i)
      nm << pfems[i]->debug_name() << ",";
    nm << " cv:" << cv << ")";
    debug_name_ = nm.str();
    
    init_cvs_node();
    for (size_type i = 0; i < pfems.size(); ++i) {
      if (pfems[i]->target_dim() != 1)
	DAL_THROW(to_be_done_error, "Vectorial fems not supported");

      for (size_type k = 0; k < pfems[i]->nb_dof(cv); ++k) {
	add_node(pfems[i]->dof_types()[k], pfems[i]->node_of_dof(cv,k));
      }
    }
  }

  /* used only when smart_global_dof_linking_ is on */
  void fem_sum::mat_trans(base_matrix &M,
			  const base_matrix &G,
			  bgeot::pgeometric_trans pgt) const {

    //cerr << "fem_sum::mat_trans " << debug_name_ << " / smart_global_dof_linking_ = " << smart_global_dof_linking_<< "\n";
    pdof_description gdof = 0, lagdof = lagrange_dof(dim());
    std::vector<pdof_description> hermdof(dim());
    for (size_type id = 0; id < dim(); ++id)
      hermdof[id] = derivative_dof(dim(), id);
    if (pfems.size()) gdof = global_dof(dim());
    gmm::copy(gmm::identity_matrix(), M);
    base_vector val(1), val2(1);
    base_matrix grad(1, dim());
    
    // very inefficient, to be optimized ...
    for (size_type ifem1 = 0, i=0; ifem1 < pfems.size(); ++ifem1) {
      pfem pf1 = pfems[ifem1];
      /* find global dofs */
      for (size_type idof1 = 0; idof1 < pf1->nb_dof(cv); ++idof1, ++i) {
	if (pf1->dof_types()[idof1] == gdof) {
	  base_vector coeff(pfems[ifem1]->nb_dof(cv));
	  coeff[idof1] = 1.0;
	  fem_interpolation_context fic(pgt, pf1, base_node(dim()), G, cv);
	  for (size_type ifem2 = 0, j=0; ifem2 < pfems.size(); ++ifem2) {
	    pfem pf2 = pfems[ifem2];
	    fem_interpolation_context fic2(pgt, pf2, base_node(dim()), G, cv);
	    for (size_type idof2 = 0; idof2 < pf2->nb_dof(cv); ++idof2, ++j) {
	      pdof_description pdd = pf2->dof_types()[idof2];
	      bool found = false;

	      base_vector coeff2(pfems[ifem2]->nb_dof(cv));
	      coeff2[idof2] = 1.0;

	      if (pdd != gdof) {
		fic.set_xref(pf2->node_of_dof(cv, idof2));
		fic2.set_xref(pf2->node_of_dof(cv, idof2));
		if (dof_weak_compatibility(pdd, lagdof) == 0) {
		  pfems[ifem1]->interpolation(fic, coeff, val, 1);

		  pfems[ifem2]->interpolation(fic2, coeff2, val2, 1);

		  M(i, j) = -val[0]*val2[0];
		  /*if (pf2->nb_dof(cv) > 4) 
		    cout << "dof " << idof2 << " (" 
			 << pf2->node_of_dof(cv,idof2) 
			 << ") " << (void*)&(*pdd) 
			 << " compatible with lagrange\n";*/
		  found = true;
		}
		else for (size_type id = 0; id < dim(); ++id) {
		  if (dof_weak_compatibility(pdd, hermdof[id]) == 0) {
		    pfems[ifem1]->interpolation_grad(fic, coeff, grad, 1);
		    M(i, j) = -grad(0, id);
		    cout << "dof " << idof2 << "compatible with hermite " << id << "\n";
		    found = true;
		  }
		}
		if (!found)
		  DAL_THROW(failure_error,
			    "Sorry, smart_global_dof_linking not "
			    "compatible with this kind of dof");
	      }
	    }
	    //if (pf2->nb_dof(cv) > 4) cout << " M=" << M << "\n";
	  }
	}
      }
    }
    
    //static int cnt=0; if(++cnt < 40) cout << "fem = " << debug_name_ << ", M = " << M << "\n";
  }

  size_type fem_sum::index_of_global_dof(size_type , size_type j) const {
    for (size_type i = 0; i < pfems.size(); ++i) {
      size_type nb = pfems[i]->nb_dof(cv);
      if (j >= nb) j -= pfems[i]->nb_dof(cv);
      else return pfems[i]->index_of_global_dof(cv, j);
    }
    DAL_THROW(failure_error, "incoherent global dof.");
  }

  void fem_sum::base_value(const base_node &, 
				 base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void fem_sum::grad_base_value(const base_node &, 
				      base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void fem_sum::hess_base_value(const base_node &, 
			     base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }

  void fem_sum::real_base_value(const fem_interpolation_context &c,
				base_tensor &t, 
				bool withM) const {
    bgeot::multi_index mi(2);
    mi[1] = target_dim(); mi[0] = nb_dof(0);
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin(), itf;

    fem_interpolation_context c0 = c;
    std::vector<base_tensor> val_e(pfems.size());
    for (size_type k = 0; k < pfems.size(); ++k) {
      if (c0.have_pfp()) {
	c0.set_pfp(fem_precomp(pfems[k], &c0.pfp()->get_point_tab()));
      } else { c0.set_pf(pfems[k]); }
      c0.base_value(val_e[k]);
    }
    
    for (dim_type q = 0; q < target_dim(); ++q) {
      for (size_type k = 0; k < pfems.size(); ++k) {
	itf = val_e[k].begin() + q * pfems[k]->nb_dof(cv);
	for (size_type i = 0; i <  pfems[k]->nb_dof(cv); ++i)
	  *it++ = *itf++;
      }
    }
    assert(it == t.end());
    if (!is_equivalent() && withM) { 
      base_tensor tt(t); 
      t.mat_transp_reduction(tt, c.M(), 0); 
    }
    //cerr << "fem_sum::real_base_value(" << c.xreal() << ")\n";
  }

  void fem_sum::real_grad_base_value(const fem_interpolation_context &c,
				     base_tensor &t, bool withM) const {
    bgeot::multi_index mi(3);
    mi[2] = c.N(); mi[1] = target_dim(); mi[0] = nb_dof(0);
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin(), itf;
    
    fem_interpolation_context c0 = c;
    std::vector<base_tensor> grad_e(pfems.size());
    for (size_type k = 0; k < pfems.size(); ++k) {
      if (c0.have_pfp()) {
	c0.set_pfp(fem_precomp(pfems[k], &c0.pfp()->get_point_tab()));
      } else { c0.set_pf(pfems[k]); }
      c0.grad_base_value(grad_e[k]);
    }

    for (dim_type k = 0; k < c.N() ; ++k) {
      for (dim_type q = 0; q < target_dim(); ++q) {
	for (size_type f = 0; f < pfems.size(); ++f) {
	  itf = grad_e[f].begin()
	    + (k * target_dim() + q) * pfems[f]->nb_dof(cv); 
	  for (size_type i = 0; i < pfems[f]->nb_dof(cv); ++i)
	    *it++ = *itf++;
	}
      }
    }
    assert(it == t.end());
    if (!is_equivalent() && withM) { 
      base_tensor tt(t); 
      t.mat_transp_reduction(tt, c.M(), 0); 
    }
  }
  
  void fem_sum::real_hess_base_value(const fem_interpolation_context &c,
				     base_tensor &t, bool withM) const {
    bgeot::multi_index mi(4);
    mi[3] = mi[2] = c.N(); mi[1] = target_dim(); mi[0] = nb_dof(0);
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin(), itf;
    
    fem_interpolation_context c0 = c;
    std::vector<base_tensor> hess_e(pfems.size());
    for (size_type k = 0; k < pfems.size(); ++k) {
      if (c0.have_pfp()) {
	c0.set_pfp(fem_precomp(pfems[k], &c0.pfp()->get_point_tab()));
      } else { c0.set_pf(pfems[k]); }
      c0.hess_base_value(hess_e[k]);
    }

    for (dim_type j = 0; j < c.N() ; ++j) {
      for (dim_type k = 0; k < c.N() ; ++k) {
	for (dim_type q = 0; q < target_dim(); ++q) {
	  for (size_type f = 0; f < pfems.size(); ++f) {
	    itf = hess_e[f].begin()
	      + ((j * c.N() + k) * target_dim() + q) * pfems[f]->nb_dof(cv); 
	    for (size_type i = 0; i < pfems[f]->nb_dof(cv); ++i)
	      *it++ = *itf++;
	  }
	}
      }
    }
    assert(it == t.end());
    if (!is_equivalent() && withM) { 
      base_tensor tt(t); 
      t.mat_transp_reduction(tt, c.M(), 0); 
    }
  }


  void mesh_fem_sum::receipt(const MESH_CLEAR &)
  { clear(); is_adapted = false; }
  void mesh_fem_sum::receipt(const MESH_DELETE &)
  { clear(); is_adapted = false; }
  void mesh_fem_sum::clear_build_methods() {
    for (size_type i = 0; i < build_methods.size(); ++i)
      del_stored_object(build_methods[i]);
    build_methods.clear();
  }
  void mesh_fem_sum::clear(void) {
    mesh_fem::clear();
    clear_build_methods();
    situations.clear();
    is_adapted = false;
  }
  
  DAL_SIMPLE_KEY(special_mflsum_key, pfem);
  
  void mesh_fem_sum::adapt(void) {
    context_check();
    clear();
    for (dal::bv_visitor i(linked_mesh().convex_index()); !i.finished(); ++i) {
      std::vector<pfem> pfems;
      bool is_cv_dep = false;
      for (size_type j = 0; j < mfs.size(); ++j) {
	if (mfs[j]->convex_index().is_in(i)) {
	  pfem pf = mfs[j]->fem_of_element(i);
	  if (pf->nb_dof(i)) {
	    pfems.push_back(pf);
	    if (pf->is_on_real_element()) is_cv_dep = true;
	  }
	}
      }
      if (pfems.size() == 1) {
	set_finite_element(i, pfems[0]);
      }
      else if (pfems.size() > 0) {
	if (situations.find(pfems) == situations.end() || is_cv_dep) {
	  pfem pf = new fem_sum(pfems, i, smart_global_dof_linking_);
	  dal::add_stored_object(new special_mflsum_key(pf), pf,
				 pf->ref_convex(0),
				 pf->node_tab(0));
	  build_methods.push_back(pf);
	  situations[pfems] = pf;
	}
	set_finite_element(i, situations[pfems]);
      }
    }
    is_adapted = true; touch();
  }


}  /* end of namespace getfem.                                             */

