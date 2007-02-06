// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2002-2007 Yves Renard
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


#include "getfem/bgeot_poly_composite.h"
#include "getfem/getfem_integration.h"
#include "getfem/getfem_mesh_im.h"
#include "getfem/dal_naming_system.h"

namespace getfem {
 
  papprox_integration
  composite_approx_int_method(const bgeot::mesh_precomposite &mp,
			      const mesh_im &mi,
			      bgeot::pconvex_ref cr) {
    approx_integration *p = new approx_integration(cr);
    base_vector w;
    for (dal::bv_visitor cv(mi.convex_index()); !cv.finished(); ++cv) {
      pintegration_method pim = mi.int_method_of_element(cv);
      bgeot::pgeometric_trans pgt = mi.linked_mesh().trans_of_convex(cv);
      if (pim->type() != IM_APPROX || !(pgt->is_linear())) {
	delete p;
	GMM_ASSERT1(false, "Approx integration and linear transformation "
		    "are required.");
      }
      papprox_integration pai = pim->approx_method();
      
      for (size_type j = 0; j < pai->nb_points_on_convex(); ++j) {
	base_node pt = pgt->transform(pim->integration_points()[j],
				      mi.linked_mesh().points_of_convex(cv));
	p->add_point(pt, pai->coeff(j) * gmm::abs(mp.det[cv]));
      }
      for (short_type f = 0; f < pgt->structure()->nb_faces(); ++f) {

	base_node barycentre = gmm::mean_value(mi.linked_mesh().points_of_face_of_convex(cv, f).begin(),
                                               mi.linked_mesh().points_of_face_of_convex(cv, f).end());
	short_type f2 = short_type(-1);
	for (short_type f3 = 0; f3 < cr->structure()->nb_faces(); ++f3) {
	  if (gmm::abs(cr->is_in_face(f3, barycentre)) < 1.0E-7)
	    { f2 = f3; break;}
	}
	if (f2 != short_type(-1)) {
	  w.resize(gmm::mat_nrows(mp.gtrans[cv]));
	  gmm::mult(mp.gtrans[cv], pgt->normals()[f], w);
	  scalar_type coeff_mul = gmm::abs(gmm::vect_norm2(w) * mp.det[cv]);
	  for (size_type j = 0; j < pai->nb_points_on_face(f); ++j) {
	    base_node pt = pgt->transform
	      (pai->point_on_face(f, j), 
	       mi.linked_mesh().points_of_convex(cv));
	    p->add_point(pt, pai->coeff_on_face(f, j) * coeff_mul, f2);
	  }
	}
      }
    }
    p->valid_method();	

    return p;
  }

  typedef dal::naming_system<integration_method>::param_list im_param_list;

  pintegration_method structured_composite_int_method(im_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
		<< params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 1 && params[1].type() == 0,
		"Bad type of parameters");
    pintegration_method pim = params[0].method();
    int k = int(::floor(params[1].num() + 0.01));
    GMM_ASSERT1(pim->type() == IM_APPROX && k > 0 && k <= 150 &&
		double(k) == params[1].num(), "Bad parameters");

    bgeot::pbasic_mesh pm;
    bgeot::pmesh_precomposite pmp;

    structured_mesh_for_convex(pim->approx_method()->ref_convex(), k, pm, pmp);
    mesh m(*pm);
    mesh_im mi(m);
    mi.set_integration_method(pm->convex_index(), pim);

    integration_method *p
      = new integration_method
      (composite_approx_int_method(*pmp, mi,
				   pim->approx_method()->ref_convex()));
    dependencies.push_back(p->approx_method()->ref_convex());
    dependencies.push_back(&(p->approx_method()->integration_points()));
    return p;
  }



  struct just_for_singleton_HCT__ { mesh m; bgeot::mesh_precomposite mp; };
  
  pintegration_method HCT_composite_int_method(im_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {

    just_for_singleton_HCT__ &jfs
      = dal::singleton<just_for_singleton_HCT__>::instance();

    GMM_ASSERT1(params.size() == 1, "Bad number of parameters : "
		<< params.size() << " should be 1.");
    GMM_ASSERT1(params[0].type() == 1, "Bad type of parameters");
    pintegration_method pim = params[0].method();
    GMM_ASSERT1(pim->type() == IM_APPROX, "Bad parameters");

    
    jfs.m.clear();
    size_type i0 = jfs.m.add_point(base_node(1.0/3.0, 1.0/3.0));
    size_type i1 = jfs.m.add_point(base_node(0.0, 0.0));
    size_type i2 = jfs.m.add_point(base_node(1.0, 0.0));
    size_type i3 = jfs.m.add_point(base_node(0.0, 1.0));
    jfs.m.add_triangle(i0, i2, i3);
    jfs.m.add_triangle(i0, i3, i1);
    jfs.m.add_triangle(i0, i1, i2);
    jfs.mp = bgeot::mesh_precomposite(jfs.m);

    mesh_im mi(jfs.m);
    mi.set_integration_method(jfs.m.convex_index(), pim);

    integration_method *p
      = new integration_method
      (composite_approx_int_method(jfs.mp, mi,
				   pim->approx_method()->ref_convex()));
    dependencies.push_back(p->approx_method()->ref_convex());
    dependencies.push_back(&(p->approx_method()->integration_points()));
    return p;
  }


  struct just_for_singleton_QUADC1__ { mesh m; bgeot::mesh_precomposite mp; };
  
  pintegration_method QUADC1_composite_int_method(im_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {

    just_for_singleton_QUADC1__ &jfs
      = dal::singleton<just_for_singleton_QUADC1__>::instance();

    GMM_ASSERT1(params.size() == 1, "Bad number of parameters : "
		<< params.size() << " should be 1.");
    GMM_ASSERT1(params[0].type() == 1, "Bad type of parameters");
    pintegration_method pim = params[0].method();
    GMM_ASSERT1(pim->type() == IM_APPROX, "Bad parameters");
    
    jfs.m.clear();
    size_type i0 = jfs.m.add_point(base_node(0.0, 0.0));
    size_type i1 = jfs.m.add_point(base_node(1.0, 0.0));
    size_type i2 = jfs.m.add_point(base_node(0.0, 1.0));
    size_type i3 = jfs.m.add_point(base_node(1.0, 1.0));
    size_type i4 = jfs.m.add_point(base_node(0.5, 0.5));
    jfs.m.add_triangle(i1, i3, i4);
    jfs.m.add_triangle(i2, i0, i4);
    jfs.m.add_triangle(i3, i2, i4);
    jfs.m.add_triangle(i0, i1, i4);
    jfs.mp = bgeot::mesh_precomposite(jfs.m);

    mesh_im mi(jfs.m);
    mi.set_integration_method(jfs.m.convex_index(), pim);

    integration_method *p = new integration_method
      (composite_approx_int_method(jfs.mp, mi,
				   bgeot::parallelepiped_of_reference(2)));
    dependencies.push_back(p->approx_method()->ref_convex());
    dependencies.push_back(&(p->approx_method()->integration_points()));
    return p;
  }




  
}  /* end of namespace getfem.                                            */
