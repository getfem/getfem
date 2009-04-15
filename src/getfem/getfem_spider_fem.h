// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2004-2008 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// As a special exception, you  may use  this file  as it is a part of a free
// software  library  without  restriction.  Specifically,  if   other  files
// instantiate  templates  or  use macros or inline functions from this file,
// or  you compile this  file  and  link  it  with other files  to produce an
// executable, this file  does  not  by itself cause the resulting executable
// to be covered  by the GNU Lesser General Public License.  This   exception
// does not  however  invalidate  any  other  reasons why the executable file
// might be covered by the GNU Lesser General Public License.
//
//===========================================================================

/**@file getfem_spider_fem.h
   @author Yves Renard <Yves.Renard@insa-lyon.fr>
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
    //   { return ::sqrt(c.xreal[0])*cos(log(c.xreal[0])); }
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
  
  
  struct Xfem_sqrtrcos : public virtual_Xfem_func {
    
    scalar_type eps;
    virtual scalar_type val(const Xfem_func_context &c) {
      return ::sqrt(c.xreal[0]) * cos( eps*log(c.xreal[0]) ); 
    }
    
    virtual base_small_vector grad(const Xfem_func_context &c) {
      base_small_vector V(2); 
      V[0] =  ( 1./::sqrt( c.xreal[0]) ) * ( cos( eps*log(c.xreal[0]) )/2. - eps*sin( eps*log(c.xreal[0]) ) ); 
      return V; 
    }
    
    virtual base_matrix hess(const Xfem_func_context &c) {
      base_matrix m(2,2); 
      m(0,0) = (1./::pow(sqrt(c.xreal[0]),3)) * (  ((-1./4.) - pow(eps,2)) *  cos( eps*log(c.xreal[0]) ) );
      return m;
    }
  };
  

  struct Xfem_sqrtrsin : public virtual_Xfem_func {
    
    scalar_type eps; 
    virtual scalar_type val(const Xfem_func_context &c) {
      return ::sqrt(c.xreal[0]) * sin( eps*log(c.xreal[0]) );  
    }
      
    virtual base_small_vector grad(const Xfem_func_context &c) {
      base_small_vector V(2); 
      V[0] =  ( 1./::sqrt(c.xreal[0]) ) * ( sin( eps*log(c.xreal[0]) )/2. + eps*cos( eps*log(c.xreal[0]) ) );
      return V;
    }
      
    virtual base_matrix hess(const Xfem_func_context &c) {
      base_matrix m(2,2); 
      m(0,0) = (1./::pow(sqrt(c.xreal[0]),3)) * (  ((-1./4.) - pow(eps,2)) *  sin( eps*log(c.xreal[0]) ) );
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
    { GMM_ASSERT1(false, "this interpolated_func has no hessian"); }
    
    
    virtual ~interpolated_transformation() {}
  };

  DAL_SIMPLE_KEY(special_cartesianfem_key, pfem);

  class spider_fem {
    
  protected :
    
    mesh cartesian;
    mesh_fem cartesian_fem;
    pfem Qk;
    Xfem *penriched_Qk;
    scalar_type R;
    unsigned Nr, Ntheta, K;
    Xfem_sqrtr Sqrtr;
    Xfem_sqrtrcos Sqrtrcos;
    Xfem_sqrtrsin Sqrtrsin;
    int bimat_enrichment;
    scalar_type epsilon;
    pfem final_fem;
    interpolated_transformation itt;

  public :
    
    pfem get_pfem(void) { return final_fem; }
    
    ~spider_fem () { 
      pfem pf = penriched_Qk;
      dal::del_stored_object(pf);
      if (final_fem) del_interpolated_fem(final_fem);
    }
    spider_fem(scalar_type R_, mesh_im &mim, unsigned Nr_, unsigned Ntheta_,
	       unsigned K_, base_small_vector translation, scalar_type theta0,
	       int bimat_enrichment_ = 0, scalar_type epsilon_ = scalar_type(0))
        : cartesian_fem(cartesian), R(R_), Nr(Nr_),
	  Ntheta(Ntheta_), K(K_), bimat_enrichment(bimat_enrichment_), epsilon(epsilon_), final_fem(0) {
        
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
	penriched_Qk = new Xfem(0);
	if(bimat_enrichment == 0){
	  cout << "Using SpiderFem homogenuous isotropic enrichment [sqrt(r)]..." << endl;
	  penriched_Qk->add_func(Qk, &Sqrtr);
	}
	else {
	  cout << "Using SpiderFem bimaterial enrichement..." << endl;
	  Sqrtrcos.eps = epsilon;
	  Sqrtrsin.eps = epsilon;
	  //cout << "epsilon = " << epsilon << endl;
	  penriched_Qk->add_func(Qk, &Sqrtrcos);
	  penriched_Qk->add_func(Qk, &Sqrtrsin);
	}
	penriched_Qk->valid();
	pfem pf = penriched_Qk;
	dal::add_stored_object(new special_cartesianfem_key(pf), pf,
			   pf->ref_convex(0),
			   pf->node_tab(0));

	cartesian_fem.set_finite_element(cartesian.convex_index(), pf);  
	GMM_ASSERT1(!cartesian_fem.is_reduced(), "To be adapted");
	dal::bit_vector blocked_dof = cartesian_fem.basic_dof_on_region(0);
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
	GMM_ASSERT3(false, "");
      }
    }
  };



}  /* end of namespace getfem.                                            */

#endif
