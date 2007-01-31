// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2004-2007 Yves Renard
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

/**@file getfem_spider_fem.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date October 29, 2004.
   @brief work in progress...
*/

#ifndef GETFEM_SPIDER_FEM_H__
#define GETFEM_SPIDER_FEM_H__

#include "getfem_interpolated_fem.h"
#include "getfem_Xfem.h"
#include "getfem_regular_meshes.h"

namespace getfem {


  struct Xfem_sqrtr : public virtual_Xfem_func {
    virtual scalar_type val(const Xfem_func_context &c)
    { return ::sqrt(c.xreal[0]); }
    //    { return ::sqrt(c.xreal[0])*cos(log(c.xreal[0])); }
    virtual base_small_vector grad(const Xfem_func_context &c)
    { base_small_vector V(2); 
      V[0] = 1. / (2.* ::sqrt(c.xreal[0])); return V; }
      // V[0] =  (1./::sqrt(c.xreal[0]))( (cos(log(::sqrt(c.xreal[0])^2)))/2. - sin(log(::sqrt(c.xreal[0])^2)) ); return V; }
    virtual base_matrix hess(const Xfem_func_context &c) {
      base_matrix m(2,2); 
      m(0,0) = 1. / (4.* ::sqrt(c.xreal[0])*c.xreal[0]);
      //m(0,0) = ((1./::sqrt(c.xreal[0]))^(3))( (3. * cos(log(::sqrt(c.xreal[0])^2)))/4. - sin(log(::sqrt(c.xreal[0])^2)));
      return m;
    }
  };

  /*
  struct Xfem_sqrtr : public virtual_Xfem_func {
    virtual scalar_type val(const Xfem_func_context &c)
    { return 1; }
    virtual base_small_vector grad(const Xfem_func_context &c)
    { base_small_vector V(2); return V; }
    virtual base_matrix hess(const Xfem_func_context &c) {
      base_matrix m(2,2); return m;
    }
  };
  */

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
    { DAL_THROW(failure_error,"this interpolated_func has no hessian"); }
    
    
    virtual ~interpolated_transformation() {}
  };


  class spider_fem {
    
  protected :
    
    mesh cartesian;
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
	M(1,1) = 2. *  M_PI;
	cartesian.transformation(M);
	bgeot::base_small_vector V(2);
	V[1] = -M_PI;
	cartesian.translation(V); 

	getfem::mesh_region border_faces;
	getfem::outer_faces_of_mesh(cartesian, border_faces);
	for (getfem::mr_visitor it(border_faces); !it.finished(); ++it) {
	  base_node un = cartesian.normal_of_face_of_convex(it.cv(), it.f());
	  un /= gmm::vect_norm2(un);
	  if (un[0] >= 0.8) cartesian.region(0).add(it.cv(), it.f());
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

	final_fem = new_interpolated_fem(cartesian_fem, mim,&itt,blocked_dof, false);
      }
    void check() {
      const interpolated_fem &ife = dynamic_cast<const interpolated_fem&>(*final_fem);
      dal::bit_vector bv = ife.interpolated_convexes();
      cerr << "interpolated_convexes: nb=" << bv.card() << "; " << bv << "\n";
      unsigned ming, maxg;
      scalar_type meang;
      ife.gauss_pts_stats(ming,maxg,meang);
      cerr << " gauss pts in interpolated mesh_fem convexes: min=" << ming << ", max=" << maxg << ", meang=" << meang << "\n";
      if (bv.card() != cartesian.convex_index().card()) {
	cerr << cartesian.convex_index().card() - bv.card() << 
	  "convexes missed by interpolated_fem, increase the "
	  "number of integration points";
	DAL_INTERNAL_ERROR("");
      }
    }
  };



}  /* end of namespace getfem.                                            */

#endif
