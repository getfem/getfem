/* *********************************************************************** */
/*                                                                         */
/* Library : GEneric Tool for Finite Element Methods (getfem)              */
/* File    : getfem_integration_composite.C : composite int methods        */
/*                                                                         */
/* Date : August 26, 2002.                                                 */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
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


#include <getfem_poly_composite.h>
#include <getfem_integration.h>
#include <getfem_mesh_fem.h>
#include <ftool_naming.h>

namespace getfem
{ 
  papprox_integration
  composite_approx_int_method(const mesh_precomposite &mp, const mesh_fem &mf,
			      bgeot::pconvex_ref cr) {
    approx_integration *p = new approx_integration(cr);
    base_vector w;
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
      pintegration_method pim = mf.int_method_of_element(cv);
      bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
      if (pim->is_exact() || !(pgt->is_linear())) {
	delete p;
	DAL_THROW(failure_error,
		  "Approx integration and linear transformation only.");
      }
      papprox_integration pai = pim->approx_method();
      
      for (size_type j = 0; j < pai->nb_points_on_convex(); ++j) {
	base_node pt = pgt->transform(pim->integration_points()[j],
				      mf.linked_mesh().points_of_convex(cv));
	p->add_point(pt, pai->coeff(j) * dal::abs(mp.det[cv]));
      }
      for (short_type f = 0; f < pgt->structure()->nb_faces(); ++f) {

	base_node barycentre = dal::mean_value(mf.linked_mesh().points_of_face_of_convex(cv, f).begin(),
                                               mf.linked_mesh().points_of_face_of_convex(cv, f).end());
	short_type f2 = short_type(-1);
	for (short_type f3 = 0; f3 < cr->structure()->nb_faces(); ++f3) {
	  if (dal::abs(cr->is_in_face(f3, barycentre)) < 1.0E-7)
	    { f2 = f3; break;}
	}
	if (f2 != short_type(-1)) {
	  w.resize(gmm::mat_nrows(mp.gtrans[cv]));
	  gmm::mult(mp.gtrans[cv], pgt->normals()[f], w);
	  scalar_type coeff_mul = dal::abs(gmm::vect_norm2(w) * mp.det[cv]);
	  for (size_type j = 0; j < pai->nb_points_on_face(f); ++j) {
	    base_node pt = pgt->transform
	      (pai->point_on_face(f, j), 
	       mf.linked_mesh().points_of_convex(cv));
	    p->add_point(pt, pai->coeff_on_face(f, j) * coeff_mul, f2);
	  }
	}
      }
    }
    p->valid_method();	

    return p;
  }

  typedef ftool::naming_system<integration_method>::param_list im_param_list;

  pintegration_method structured_composite_int_method(im_param_list &params) {

    if (params.size() != 2)
      DAL_THROW(failure_error, 
	  "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 1 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    pintegration_method pim = params[0].method();
    int k = int(::floor(params[1].num() + 0.01));
    if (pim->is_exact() || k <= 0 || k > 150 || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");

    pgetfem_mesh pm;
    pmesh_precomposite pmp;

    structured_mesh_for_convex(pim->approx_method()->ref_convex(), k, pm, pmp);
    mesh_fem mf(*pm);
    mf.set_finite_element(pm->convex_index(),
			  classical_fem(pm->trans_of_convex(0), 0), pim);

    return new integration_method
      (composite_approx_int_method(*pmp, mf,
				   pim->approx_method()->ref_convex()));
  }
  
}  /* end of namespace getfem.                                            */
