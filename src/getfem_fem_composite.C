/* *********************************************************************** */
/*                                                                         */
/* Library : GEneric Tool for Finite Element Methods (getfem)              */
/* File    : getfem_fem_composite.C : composite fem                        */
/*                                                                         */
/* Date : August 26, 2002.                                                 */
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


#include <getfem_poly_composite.h>
#include <getfem_integration.h>
#include <getfem_mesh_fem.h>

namespace getfem
{ 
  typedef const fem<polynomial_composite> * ppolycompfem;

  ppolycompfem composite_fe_method(const mesh_precomposite &mp, 
				   const mesh_fem &mf, bgeot::pconvex_ref cr) {
    
    if (&(mf.linked_mesh()) != &(mp.linked_mesh()))
      DAL_THROW(failure_error, "Meshes are different.");
    dal::bit_vector nn = mf.convex_index();
    fem<polynomial_composite> *p = new fem<polynomial_composite>;

    p->ref_convex() = cr;
    p->is_equivalent() = true;
    p->is_polynomial() = p->is_lagrange() = false;
    p->estimated_degree() = 0;
    p->init_cvs_node();

    std::vector<polynomial_composite> base(mf.nb_dof());
    std::fill(base.begin(), base.end(), polynomial_composite(mp));
    std::vector<pdof_description> dofd(mf.nb_dof());
    
    for (size_type cv = nn.take_first(); cv != size_type(-1); cv << nn) {
      pfem pf1 = mf.fem_of_element(cv);
      if (!(pf1->is_equivalent() && pf1->is_polynomial())) {
	delete p;
	DAL_THROW(failure_error, "Only for polynomial and equivalent fem.");
      }
      ppolyfem pf = ppolyfem(pf1);
      p->estimated_degree() = std::max(p->estimated_degree(),
				       pf->estimated_degree());
      for (size_type k = 0; k < pf->nb_dof(); ++k) {
	size_type igl = mf.ind_dof_of_element(cv)[k];
	base_poly fu = pf->base()[k];
	base[igl].poly_of_subelt(cv) = fu;
	dofd[igl] = pf->dof_types()[k];
      }
    }
    p->base().resize(mf.nb_dof());
    for (size_type k = 0; k < mf.nb_dof(); ++k) {  
      p->add_node(dofd[k], mf.point_of_dof(k));
      p->base()[k] = base[k];
    }
    return p;
  }
  
}  /* end of namespace getfem.                                            */
