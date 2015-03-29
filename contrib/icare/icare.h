/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2002-2015 Michel Fournié, Julien Pommier,
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/
#ifndef NAVIER_STOKES_H_
#define NAVIER_STOKES_H_

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_modeling.h"
#include "gmm/gmm.h"

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

namespace getfem {


	/*
	 * Adams-Bashfort - For explict treatment
	 */
	
	template<typename VEC1, typename VEC>
	void asm_NS_uuT_Explicite(const VEC1 &V, const mesh_im &mim,
					const mesh_fem &mf, const VEC &U0,
					const mesh_region &rg = mesh_region::all_convexes()) {
		generic_assembly assem;    
		assem.set("u=data(#1); "
				  "t = comp(vBase(#1).vBase(#1).vGrad(#1));"
				  "V(#1)+=u(l).u(i).t(i,k,:,j,l,j,k)"
				   );
		assem.push_mi(mim);
		assem.push_mf(mf);
		assem.push_data(U0);
		assem.push_vec(V);
		assem.assembly(rg);
	}
	
	
	
  /*
   * div(uuT) term. 
   */
  
  template<typename MAT, typename VEC>
  void asm_NS_uuT(const MAT &M, const mesh_im &mim,
		  const mesh_fem &mf, const VEC &U0,
		  const mesh_region &rg = mesh_region::all_convexes()) {
    generic_assembly assem;    
    assem.set("u=data(#1); "
	      "t = comp(vBase(#1).vBase(#1).vGrad(#1));"
	      "M(#1, #1)+=u(i).t(i,k,:,j,:,j,k)"
	      // ";M(#1, #1)+=u(i).t(i,k,:,k,:,j,j)" // optional (or *0.5)
	      );
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(U0);
    assem.push_mat(const_cast<MAT&>(M));
    assem.assembly(rg);
  }

  /*
   * Boundary term for the pressure 
   */
  
  template<typename MAT>
  void asm_B_boundary(const MAT &M, const mesh_im &mim,
		      const mesh_fem &mfu, const mesh_fem &mfp,
		      const mesh_region &rg) {
    generic_assembly assem;    
    assem.set("M$1(#2, #1)+=-comp(Base(#2).vBase(#1).Normal())(:,:,i,i);");
    assem.push_mi(mim);
    assem.push_mf(mfu);
    assem.push_mf(mfp);
    assem.push_mat(const_cast<MAT&>(M));
    assem.assembly(rg);
  }



 //template<typename VEC, typename VECTOR>
// void traineePortance2D(VEC &Cd, VEC &Cl, const mesh_im &mim,
//			const mesh_fem &mf, const VECTOR &lambda,
//			const mesh_region &rg) {
//
//   generic_assembly assem;
//   assem.set("l=data$1(#1); V$1()+=comp(vBase(#1).Normal())(i,1,1).l(i);V$2()+=comp(vBase(#1).Normal())(i,1,2).l(i);");
//   
//   assem.push_mi(mim);
//   
//   assem.push_mf(mf);
//   assem.push_data(lambda);
//
//   assem.push_vec(Cd);
//   assem.push_vec(Cl);
//
//   assem.assembly(rg); //region = bord du cylindre
//						 
// }

	
	
	template<typename VEC, typename VECTOR>
	void ClCd2D(VEC &Cxn, VEC &Cxp,VEC &Cyn, VEC &Cyp, const mesh_im &mim,
						   const mesh_fem &mf1,const mesh_fem &mf2, const VECTOR &U0,const VECTOR &P0, const mesh_region &rg) {
		
		generic_assembly assem;
		
		assem.set("u=data$1(#1); p=data$2(#2);"
				  "t1=comp(vGrad(#1).Normal());"
				  "V$1()+=2*t1(i,1,1,1).u(i) + t1(i,1,2,2).u(i) + t1(i,2,1,2).u(i);"
				  "t2=comp(Base(#2).Normal());"
				  "V$2()+=-1*t2(i,1).p(i);"
				  "t3=comp(vGrad(#1).Normal());"
				  "V$3()+=2*t3(i,2,2,2).u(i) + t3(i,1,2,1).u(i) + t3(i,2,1,1).u(i);"
				  "t4=comp(Base(#2).Normal());"
				  "V$4()+=-1*t4(i,2).p(i);");
		
		
		//assem.set("nuDxU=data$1(#1); nuDyU=data$2(#1); nuDxV=data$3(#1); nuDyV=data$4(#1); p=data$5(#2);"
		//		  "t1=comp(vBase(#1).Normal());"
		//		  "V$1()+=2*t1(i,1,1).nuDxU(i) + t1(i,1,2).nuDyU(i) + t1(i,1,2).nuDxV(i);"
		//		  "t2=comp(Base(#2).Normal());"
		//		  "V$2()+=-1*t2(i,1).p(i);"
		//		  "t3=comp(vBase(#1).Normal());"
		//		  "V$3()+=2*t3(i,1,2).nuDyV(i) + t3(i,1,1).nuDyU(i) + t3(i,1,1).nuDxV(i);"
		//		  "t4=comp(Base(#2).Normal());"
		//		  "V$4()+=-1*t4(i,2).p(i);");
		
		
		assem.push_mi(mim);
		assem.push_mf(mf1);
		assem.push_mf(mf2);

		assem.push_data(U0);
		assem.push_data(P0);
		
		assem.push_vec(Cxn); 
		assem.push_vec(Cxp);
		assem.push_vec(Cyn); 
		assem.push_vec(Cyp);
		
		assem.assembly(rg); //region = bord du cylindre
		
	}
		

 template<typename VEC, typename VECTOR>
 void ClCd3D(VEC &Cxn, VEC &Cxp,VEC &Cyn, VEC &Cyp, const mesh_im &mim,
			const mesh_fem &mf1, const mesh_fem &mf2,const VECTOR &U0, 
			const VECTOR &P0, const mesh_region &rg) {

   generic_assembly assem;

	assem.set("u=data$1(#1); p=data$2(#2);"
			  "t1=comp(vGrad(#1).Normal());"
			  "V$1()+=2*t1(i,1,1,1).u(i) + t1(i,1,2,2).u(i) + t1(i,2,1,2).u(i) + t1(i,3,1,3).u(i) + t1(i,1,3,3).u(i);"
			  "t2=comp(Base(#2).Normal());"
			  "V$2()+=-1*t2(i,1).p(i);"
			  "t3=comp(vGrad(#1).Normal());"
			  "V$3()+=2*t3(i,2,2,2).u(i) + t3(i,1,2,1).u(i) + t3(i,2,1,1).u(i) + t3(i,3,2,3).u(i) + t3(i,2,3,3).u(i);"
			  "t4=comp(Base(#2).Normal());"
			  "V$4()+=-1*t4(i,2).p(i);");
	 
   assem.push_mi(mim);

   assem.push_mf(mf1);
   assem.push_mf(mf2);

   assem.push_data(U0);   
   assem.push_data(P0);
   
   assem.push_vec(Cxn); 
   assem.push_vec(Cxp);
   assem.push_vec(Cyn); 
   assem.push_vec(Cyp);
  
   assem.assembly(rg); //region = bord du cylindre
						 
 }





  /* ******************************************************************** */
  /*		A non-reflective condition.                               */
  /* ******************************************************************** */
  // idée: construire la condition non-réflective comme une condition de dirichlet en explicitant tous les termes: Unp1= Un+dt*Un.N*lamda




// initialisation des CL
// résolution du pb -> nouveau lambda<--| 
//mise à jour de CL de sortie        -->|  


 template<typename VECTOR, typename VECT>
  void asm_nonref_right_hand_side(VECTOR &R,
				  const mesh_im &mim,
				  const mesh_fem &mf_u,
				  const mesh_fem &mf_mult,
				  const VECT &Un,
				  const mesh_region &rg) {
    generic_assembly assem;
    // construction du terme de droite dans [M]*Unp1=F
    
    // mise en place de Un.N + Un.N*(dUn/dn).N

    std::stringstream ss;
    ss << "u=data$1(#1); "
     "V(#1)+=comp(vBase(#1).vBase(#2))(j,i,:,i).u(j)";
 
    assem.set(ss.str());
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_mult);
    assem.push_data(Un);
    assem.push_vec(R);
    assem.assembly(rg);
  }

  /* non reflective boundary conditions */


  // construction du terme de droite dans [M]*Unp1=F
  // mise en place de Un + Un.N*(dUn/dn).N
  template <typename VEC1, typename VEC2>
  void asm_basic_non_reflective_bc(VEC1 &VV, 
				   const getfem::mesh_im &mim, 
				   const getfem::mesh_fem &mf_u, 
				   const VEC2 &Un0, 
				   const getfem::mesh_fem &mf_mult,
				   scalar_type dt,
				   const getfem::mesh_region &rg) {
    getfem::generic_assembly assem;
    std::stringstream ss;
    ss << "u=data$1(#1); "
      "V(#2)+="
       << -dt << "*"
      " comp(vBase(#1).Normal().vGrad(#1).Normal().vBase(#2))(l,i,i,m,k,j,j,:,k).u(l).u(m)+"     //"(l,i,i,m,j,k,j,:,k).u(l).u(m)+"
      "comp(vBase(#1).vBase(#2))(j,i,:,i).u(j)";                                         //"comp(vBase(#1).vBase(#2).u(j))(j,i,:,i)";
    assem.set(ss.str());
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_mult);  
    assem.push_data(Un0);
    assem.push_vec(VV);
    assem.assembly(rg);
  }

  template<typename VECT1> class improved_non_reflective_bc_nonlinear_term 
    : public getfem::nonlinear_elem_term {
    const mesh_fem &mf;
    std::vector<scalar_type> U;
    scalar_type dt, nu;
    size_type N;
    base_vector coeff, valU;
    base_matrix gradU, hessU;
    bgeot::multi_index sizes_;
  public:
    improved_non_reflective_bc_nonlinear_term
    (const mesh_fem &mf_, const VECT1 &U_, scalar_type dt_, scalar_type nu_)
      : mf(mf_), U(mf_.nb_basic_dof()), dt(dt_), nu(nu_), N(mf_.get_qdim()), 
	valU(N), gradU(N, N), hessU(N,N*N) {
      sizes_.resize(1); sizes_[0] = short_type(N); /*assert(N == 2);*/
      mf.extend_vector(U_, U);
    }
    const bgeot::multi_index &sizes(size_type) const {  return sizes_; }
    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num();
      coeff.resize(mf.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf.ind_basic_dof_of_element(cv))),
		coeff);
      ctx.pf()->interpolation(ctx, coeff, valU, mf.get_qdim());
      ctx.pf()->interpolation_grad(ctx, coeff, gradU, mf.get_qdim());
      ctx.pf()->interpolation_hess(ctx, coeff, hessU, mf.get_qdim());

      if (N==2){
	t[0] = (valU[0] +  nu *dt* hessU(0,3) ) / (1 + gradU(0,0)*dt);
	t[1] =  valU[1] +   nu *dt* hessU(1,3) - dt * t[0] * gradU(1,0);
      }
      if (N==3) {
	t[0] = (valU[0] + nu *dt* hessU(0,4) ) / (1 + gradU(0,0)*dt);
	t[1] =  valU[1] + nu *dt* hessU(1,4)  - dt * t[0] * gradU(1,0);
	t[2] =  valU[2] + nu *dt* hessU(2,4)  - dt * t[0] * gradU(2,0);
      }
    }
  };


  template <typename VEC1, typename VEC2>
  void asm_improved_non_reflective_bc(VEC1 &VV, 
				      const getfem::mesh_im &mim, 
				      const getfem::mesh_fem &mf_u, 
				      const VEC2 &Un0, 
				      const getfem::mesh_fem &mf_mult,
				      scalar_type dt, scalar_type nu,
				      const getfem::mesh_region &rg) {
    getfem::generic_assembly assem;
    improved_non_reflective_bc_nonlinear_term<VEC2> nlterm(mf_u, Un0, dt, nu);
    assem.set("V(#2)+=comp(NonLin(#1).vBase(#2))(i,:,i)");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_mult);  
    assem.push_vec(VV);
    assem.push_nonlinear_term(&nlterm);
    assem.assembly(rg);
  }



}
#endif
