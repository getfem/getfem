// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Elements Methods (getfem)
// File    : getfem_mesh_im_level_set.cc : adaptable integration methods
//           on convex meshes.
// Date    : February 02, 2005.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2005-2005 Yves Renard
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

#include <getfem_mesh_im_level_set.h>
#include <getfem_mesher.h>

#ifndef GETFEM_HAVE_QHULL_QHULL_H

    DAL_THROW(failure_error, "Qhull header files not installed. "
	      "This part of getfem++ needs Qhull."
	      "Install qhull library and reinstall Getfem++ library.");
    
#else
    
    


namespace getfem {
  
  void mesh_im_level_set::receipt(const MESH_CLEAR &) { clear(); }
  void mesh_im_level_set::receipt(const MESH_DELETE &) { clear(); }
  void mesh_im_level_set::clear(void) { mesh_im::clear(); level_sets.clear(); }

  mesh_im_level_set::mesh_im_level_set(getfem_mesh &me,
				       pintegration_method reg,
				       pintegration_method sing) 
    : mesh_im(me), cut_im(me) {
    regular_simplex_pim = reg;
    singular_simplex_pim = (sing == 0) ? reg : sing;
  }

  pintegration_method 
  mesh_im_level_set::int_method_of_element(size_type cv) const {
    if (cut_im.convex_index().is_in(cv)) 
      return cut_im.int_method_of_element(cv); 
    else return mesh_im::int_method_of_element(cv);
  }

  void mesh_im_level_set::cut_element(size_type cv, const dal::bit_vector &primary,
				      const dal::bit_vector &secondary) {
    cout << "cutting element " << cv << endl;
    std::auto_ptr<mesher_signed_distance> ref_element;
    std::vector<mesher_level_set> mesher_level_sets;
    std::vector<const mesher_signed_distance *> signed_dist;
    std::vector<base_node> gauss_points;
    std::vector<scalar_type> gauss_weights;

    bgeot::pgeometric_trans pgt = linked_mesh().trans_of_convex(cv);
    unsigned found = 0;
    size_type n = pgt->structure()->dim();
    size_type nbp = pgt->basic_structure()->nb_points();

    size_type nbtotls = primary.card() + secondary.card();
    if (nbtotls > 16)
      DAL_THROW(failure_error, "An element is cutted by more than 16 level_set, aborting");

    /* Identifying simplexes.                                          */

    if (nbp == n+1 && pgt->basic_structure() == bgeot::simplex_structure(n)) {
	ref_element.reset(new mesher_simplex_ref(n)); found = 1;
    }
    
    /* Identifying parallelepiped.                                     */

    if (!found && nbp == (size_type(1) << n) &&
	pgt->basic_structure() == bgeot::parallelepiped_structure(n)) {
      base_node rmin(n), rmax(n);
      std::fill(rmax.begin(), rmax.end(), scalar_type(1));
      ref_element.reset(new mesher_rectangle(rmin, rmax)); found = 2;
    }

    /* Identifying prisms.                                             */
 
    if (!found && nbp == 2 * n &&
	pgt->basic_structure() == bgeot::prism_structure(n))
      { ref_element.reset(new mesher_prism_ref(n)); found = 3; }
    
    if (!found) 
      DAL_THROW(to_be_done_error,
		"This element is not taken into account. Contact us");

    scalar_type h0 = 0.5;
    // + Boucle est-ce que c'est bon ou il faut faire plus fin.

    gauss_points.clear();
    gauss_weights.clear();
    for (size_type count = 0; count < (size_type(1) << nbtotls); ++count) {
      // signed distances for the level sets
      size_type K = 1;
      signed_dist.clear();
      mesher_level_sets.clear();
      signed_dist.push_back(ref_element.get());
      unsigned k = 0, l = 0;

      for (std::set<plevel_set>::const_iterator it = level_sets.begin();
	   it != level_sets.end(); ++it, ++l) {
	if (primary[l]) {
	  unsigned deg = (*it)->degree() + (found == 3 ? (*it)->degree() : 0) + (found == 2 ? (*it)->degree()*(n-1) : 0);
	  K = std::max(deg, K);
	  mesher_level_sets.push_back((*it)->mls_of_convex(cv, 0, count & (1 << k)));
	  signed_dist.push_back(&mesher_level_sets[mesher_level_sets.size()-1]);
	  ++k;
	  if (secondary[l]) {
	    mesher_level_sets.push_back((*it)->mls_of_convex(cv, 1, count & (1 << k)));
	    signed_dist.push_back(&mesher_level_sets[mesher_level_sets.size()-1]);
	    ++k;
	  }
	}
      }
      
      mesher_intersection final_dist(signed_dist);
      getfem_mesh mesh; // + boucle sur 
      std::vector<base_node> fixed_points; // mettre les sommets
      build_mesh(mesh, final_dist, h0/50, fixed_points, K, 1 /*noise*/, 
		 100, 2 /*prefind*/);
      
      for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i) {
	base_matrix G;
	vectors_to_base_matrix(G, mesh.points_of_convex(cv));
	papprox_integration pai = regular_simplex_pim->approx_method();
	for (size_type j = 0; j < pai->nb_points_on_convex(); ++j) {
	  bgeot::geotrans_interpolation_context c(mesh.trans_of_convex(i), pai->point(j), G);
	  gauss_points.push_back(c.xreal());
	  gauss_weights.push_back(pai->coeff(j) * c.J());
	}
      }
    }
      
    // + test ...
    scalar_type wtot(0);
    for (size_type k = 0; k < gauss_weights.size(); ++k) wtot += gauss_weights[k];
    cout << "The Result : " << wtot << endl;

  }


  void mesh_im_level_set::adapt(void) {

    // compute the elements touched by each level set
    // for each element touched, compute the sub mesh
    //   then compute the adapted integration method
    for (dal::bv_visitor cv(linked_mesh().convex_index()); 
	 !cv.finished(); ++cv) {
      dal::bit_vector prim, sec;
      find_crossing_level_set(cv, prim, sec);
      cout << "element " << cv << " cutted level sets : " << prim << endl;
      if (prim.card()) cut_element(cv, prim, sec);
    }
  }

  bool mesh_im_level_set::is_crossed_by(size_type cv, plevel_set ls, unsigned lsnum) {
    const mesh_fem &mf = ls->get_mesh_fem();
    ref_mesh_dof_ind_ct dofs = mf.ind_dof_of_element(cv);
    base_vector v; v.reserve(dofs.size());
    int p = -2;
    base_node G = dal::mean_value(linked_mesh().points_of_convex(cv));
    scalar_type radius2 = 0, dmin = 1e30;

    scalar_type EPS = 1e-8;

    /* easy cases: 
     - different sign on the dof nodes => intersection for sure
     - min value of the levelset greater than the radius of the convex
     => no intersection (assuming the level set is locally a signed
     distance, and the geotrans of the convex is not too strong..)
    */ 
    for (ref_mesh_dof_ind_ct::const_iterator it = dofs.begin(); it != dofs.end(); ++it) {
      radius2 = std::max(radius2, 
			 gmm::vect_dist2_sqr(G, mf.point_of_dof(*it)));
      v.push_back(ls->values()[*it]);
      int p2 = ( (v.back() < -EPS) ? -1 :
		 ((v.back() > EPS) ? +1 : 0));
      if (p == -2) p=p2; else
	if (p*p2 < 0) return true;
      dmin = std::min(dmin, gmm::abs(v.back()));
    }
    cout << "v = " << v << endl;
    cout << "G = " << G << endl;
    cout << "dmin = " << dmin << endl;
    cout << "radius2 = " << radius2 << endl;
    
    
    if (dmin*dmin > radius2*2.) return false;

    pfem pf = mf.fem_of_element(cv);
    size_type nbdof = mf.nb_dof_of_element(cv);
    assert(pf->target_dim() == 1);
    mesher_level_set mls = ls->mls_of_convex(cv, lsnum);
    for (size_type i=0; i < nbdof; ++i) {
      base_node Pt = pf->node_convex(cv).points()[i];
      base_node grad;
      scalar_type d = mls.grad(Pt, grad);
      unsigned niter = 0;
      cout << "trying node " << i << " pt = " << Pt << " d = " << d << " ptreal = " << mf.point_of_dof(dofs[i]) << " val = " << ls->values()[dofs[i]] << endl;
      while (gmm::abs(d) > EPS) {
	gmm::add(gmm::scaled(grad, -d), Pt);
	d = mls.grad(Pt, grad);
	// cout << "d = " << d << " pt = " << Pt << " Grad = " << grad << endl;
	if (++niter > 100) { DAL_WARNING(0,"too much iter"); break; }
      }
      cout << "  -> proj=" << Pt << " isin=" << pf->ref_convex(cv)->is_in(Pt) << "\n";
      if (pf->ref_convex(cv)->is_in(Pt) < EPS) return true;
    }
    return false;
  }

  void mesh_im_level_set::find_crossing_level_set(size_type cv, dal::bit_vector &prim, dal::bit_vector &sec) {
    prim.clear(); sec.clear();
    unsigned lsnum = 0;
    for (std::set<plevel_set>::const_iterator it = level_sets.begin(); 
	 it != level_sets.end(); ++it, ++lsnum) {
      if (is_crossed_by(cv, *it, 0)) {
	prim.add(lsnum);
	if ((*it)->has_secondary() && is_crossed_by(cv, *it, 1)) sec.add(lsnum);
      }
    }
  }


#endif
  
}  /* end of namespace getfem.                                             */



