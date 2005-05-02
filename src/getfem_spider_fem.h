// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_spider_fem.h : definition of a finite element
//           method which interpolates a fem on a different mesh.
// Date    : October 29, 2004.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2004-2005 Yves Renard
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



#ifndef GETFEM_SPIDER_FEM_H__
#define GETFEM_SPIDER_FEM_H__

#include <getfem_interpolated_fem.h>
#include <getfem_Xfem.h>
#include <getfem_regular_meshes.h>

namespace getfem {


  struct Xfem_sqrtr : public virtual_Xfem_func {
    virtual scalar_type val(const Xfem_func_context &c)
    { return ::sqrt(c.xreal[0]); }
    virtual base_small_vector grad(const Xfem_func_context &c)
    { base_small_vector V(2); V[0] = 1. / (2.* ::sqrt(c.xreal[0])); return V; }
    virtual base_matrix hess(const Xfem_func_context &c) {
      base_matrix m(2,2); m(0,0) = -1. / (4.* ::sqrt(c.xreal[0])*c.xreal[0]);
      return m;
    }
  };

  struct interpolated_transformation : public virtual_interpolated_func{
    /* Polar transformation and its gradient. */
    base_small_vector trans;
    scalar_type theta0;
     
    virtual void val(const base_node &xreal, base_node &v) const {
      base_node w =  xreal - trans;
      v[0] = gmm::vect_norm2(w);
      v[1] = atan2(w[1], w[0]) - theta0;
    }
    virtual void grad(const base_small_vector &xreal, base_matrix &m) const {
      base_node w =  xreal - trans;
      scalar_type r = gmm::vect_norm2(w); assert(gmm::abs(r)>1e-30);
      m(0,0) = w[0] / r; m(0,1) = w[1] / r;
      m(1,0) = -w[1] / gmm::sqr(r); m(1,1) = w[0] / gmm::sqr(r);
    }
    virtual void hess(const base_node &, base_matrix &) const
    { DAL_THROW(dal::failure_error,"this interpolated_func has no hessian"); }
    
    
    virtual ~interpolated_transformation() {}
  };


  class spider_fem {
    
  protected :
    
    getfem_mesh cartesian;
    mesh_fem cartesian_fem;
    pfem Qk;
    Xfem enriched_Qk;
    scalar_type R;
    unsigned Nr, Ntheta, K;
    Xfem_sqrtr Sqrtr;
    pfem final_fem;
    interpolated_transformation itt;

  public :
    
    pfem get_pfem(void) { return final_fem; }
    
    ~spider_fem () { if (final_fem) del_interpolated_fem(final_fem); }
    
    spider_fem(scalar_type R_, mesh_im &mim, unsigned Nr_, unsigned Ntheta_,
	       unsigned K_, base_small_vector translation, scalar_type theta0)
        : cartesian_fem(cartesian), enriched_Qk(0), R(R_), Nr(Nr_),
	  Ntheta(Ntheta_), K(K_), final_fem(0) {
        
	itt.trans = translation;
	itt.theta0 = theta0;

        /* make the cartesian mesh */
        bgeot::pgeometric_trans pgt = 
	  bgeot::geometric_trans_descriptor("GT_LINEAR_QK(2)");
        std::vector<size_type> nsubdiv(2);
	nsubdiv[0] = Nr; nsubdiv[1] = Ntheta;
        getfem::regular_unit_mesh(cartesian, nsubdiv, pgt, false);
	bgeot::base_matrix M(2,2);
	M(0,0) = R;   
	M(1,1) = 2. * M_PI;
	cartesian.transformation(M);
	bgeot::base_small_vector V(2);
	V[1] = -M_PI;
	cartesian.translation(V); 

	getfem::convex_face_ct border_faces;
	getfem::outer_faces_of_mesh(cartesian, border_faces);
	for (getfem::convex_face_ct::const_iterator it = border_faces.begin();
	     it != border_faces.end(); ++it) {
	  base_node un = cartesian.normal_of_face_of_convex(it->cv, it->f);
	  un /= gmm::vect_norm2(un);
	  if (un[0] >= 0.8) cartesian.add_face_to_set(0, it->cv, it->f);
	}

	std::stringstream Qkname;
	Qkname << "FEM_QK(2," << K << ")";
	Qk = fem_descriptor(Qkname.str());
	enriched_Qk.add_func(Qk, &Sqrtr);
	enriched_Qk.valid();

	cartesian_fem.set_finite_element(cartesian.convex_index(),
					 &enriched_Qk);  
	dal::bit_vector blocked_dof = cartesian_fem.dof_on_set(0);
	//	cout << "blocked dofs = " <<  blocked_dof << endl;

	final_fem = new_interpolated_fem(cartesian_fem, mim,&itt,blocked_dof);
      }
  };



}  /* end of namespace getfem.                                            */

#endif
