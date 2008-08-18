// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2007-2008 Yves Renard,  Julien Pommier
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
// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.
//
//===========================================================================

#include "crack_bilaplacian.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_fourth_order.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_superlu.h"
#include "getfem/getfem_derivatives.h"


namespace getfem {
  /* build a "ring" of convexes of given center and radius */
  dal::bit_vector 
  build_sif_ring_from_mesh(const mesh &m, 
                           base_node center, scalar_type r) {
    dal::bit_vector bv;
    scalar_type r2 = r - 1e-4;
    unsigned in = 0;
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv) {
      unsigned in1=0, out1=0;
      unsigned in2=0, out2=0;
      for (unsigned i=0; i < m.nb_points_of_convex(cv); ++i) {
        base_node P = m.points_of_convex(cv)[i];
        if (gmm::vect_dist2(P, center) < r) ++in1; else ++out1;
        if (gmm::vect_dist2(P, center) < r2) ++in2; else ++out2;
      }
      if ((in1 && out1) || (in2 && out2)) bv.add(cv);
      in += in1;
    }
    if (in < 3) GMM_WARNING1("looks like the radius is too small...");
    return bv;
  }

  /* return the crack tip in P,
     and the outgoing tangent of the crack in T,
     and the normal in N */
  void get_crack_tip_and_orientation(const level_set &/* ls */,
                                     base_node &P, 
                                     base_small_vector &T, base_small_vector &N) {
    cerr << __PRETTY_FUNCTION__ << " IS TO BE DONE\n";
    /* too lazy to do it now */
    P.resize(2); P[0] = 0.; P[1] = 0;
    T.resize(2); T[0] = 1; T[1] = 0;
    N.resize(2); N[0] = 0; N[1] = 1;
  }




  /* return a (very rough) estimate of the stress intensity factors using 
     the local displacement near the crack tip */
  template <typename VECT> 
  void estimate_crack_stress_intensity_factors(const level_set &ls,
                                               const mesh_fem &mf, const VECT &U, 
                                               scalar_type young_modulus, 
                                               scalar_type &KI, 
                                               scalar_type &KII, double EPS=1e-2) {
    base_node P(2);
    base_small_vector T(2),N(2); 
    get_crack_tip_and_orientation(ls, P, T, N);
    base_node P1 = P - EPS*T + EPS/100.*N;
    base_node P2 = P - EPS*T - EPS/100.*N;
    std::vector<double> V(4);
    getfem::mesh_trans_inv mti(mf.linked_mesh());
    mti.add_point(P1);
    mti.add_point(P2);
    cout << "P1 = " << P1 << ", P2=" << P2 << "\n";
    base_matrix M;
    getfem::interpolation(mf, mti, U, V, false);
    KI = (V[1] - V[3])/sqrt(EPS/(2*M_PI)) * young_modulus / 8.0;
    KII = (V[0] - V[2])/sqrt(EPS/(2*M_PI)) * young_modulus / 8.0;
  }

/* compute with great precision the stress intensity factors using
     integral computed on a ring around the crack tip, for BILAPLACIAN PLATE PROBLEM */
  template <typename VECT>
  void compute_crack_stress_intensity_factors_KL(const level_set &ls,
                                              const mesh_im &mim,
                                              const mesh_fem &mf, 
                                              const VECT &U,
                                              scalar_type ring_radius, 
                                              scalar_type D, scalar_type nu,
                                              scalar_type young_modulus, scalar_type epsilon,
                                              scalar_type &KI, scalar_type &KII) {
    const mesh &mring = mim.linked_mesh();
    mesh_fem_global_function mf_mode(mring, 1), mf_mode_3(mring, 1);
    mesh_fem mf_q(mring,1);

    std::vector<pglobal_function> cfun(2), cfun_3(4) ;
    for (unsigned j=4; j < 6; ++j)
        cfun[j-4] = bilaplacian_crack_singular(j, ls, nu, 0.) ;
    mf_mode.set_functions(cfun);

    cfun_3[0] = bilaplacian_crack_singular(311, ls, nu, 0.) ;
    cfun_3[1] = bilaplacian_crack_singular(312, ls, nu, 0.) ;
    cfun_3[2] = bilaplacian_crack_singular(321, ls, nu, 0.) ;
    cfun_3[3] = bilaplacian_crack_singular(322, ls, nu, 0.) ;
    mf_mode_3.set_functions(cfun_3);
    mf_mode_3.set_qdim(2) ;

    mf_mode.set_qdim(1);
    mf_q.set_classical_finite_element(1);

    base_node crack_tip;
    base_small_vector T, N;
    get_crack_tip_and_orientation(ls, crack_tip, T, N);

    dal::bit_vector cvring = build_sif_ring_from_mesh(mring, crack_tip, 
                                                      ring_radius);

    /* fill the "q" ring field with a approximately linear field, equal to 
       1 on the inner boundary, and equal to zero on the outer boundary */
    std::vector<scalar_type> q(mf_q.nb_dof());
    for (unsigned d = 0; d < mf_q.nb_dof(); ++d) {
      base_node P = mf_q.point_of_dof(d);
      q[d] = (gmm::vect_dist2(P, crack_tip) > ring_radius) ? 0 : 1;
    }

    base_vector U_mode(mf_mode.nb_dof()); assert(U_mode.size() == 2);
    base_vector U_mode_3(mf_mode_3.nb_dof()); assert(U_mode_3.size() == 8);


    // Remaillage conforme des mailles coupées par la fissure


    mesh_fem mf_du(mim.linked_mesh(), 1) ;
    mf_du.set_classical_finite_element(3);
    plain_vector GRAD(mf_du.nb_dof() * 2) ;
    getfem::compute_gradient(mf, mf_du, U, GRAD) ;
    // separation des 2 composantes 
    plain_vector D1(mf_du.nb_dof()), D2(mf_du.nb_dof()) ;
    for (unsigned i=0 ; i < mf_du.nb_dof() ; ++i){
       D1[i] = GRAD[2*i] ;
       D2[i] = GRAD[2*i+1] ;
    }
    cout << "1er derivation OK \n " ;

    mesh_fem mf_d2u(mim.linked_mesh(), 1) ;
    mf_d2u.set_classical_discontinuous_finite_element(2);
    plain_vector GRAD_U1(mf_d2u.nb_dof() * 2), GRAD_U2(mf_d2u.nb_dof() * 2) ;
    getfem::compute_gradient(mf_du, mf_d2u, D1, GRAD_U1) ;
    getfem::compute_gradient(mf_du, mf_d2u, D2, GRAD_U2) ;
    cout << "2e derivation : OK \n" ; 
    // calcul du laplacien
    plain_vector LAP(mf_d2u.nb_dof()) ;
    cout << "taille GRAD_U1 = " << GRAD_U1.size() << "\n" ;
    cout << "taille GRAD_U1 = " << GRAD_U1.size() << "\n" ;
    
    for (unsigned i=0 ; i < mf_d2u.nb_dof() ; ++i)
       LAP[i] = GRAD_U1[2*i] + GRAD_U2[2*i+1] ;

    // dérivation :
    mesh_fem mf_d3u(mim.linked_mesh(), 1) ;
    mf_d3u.set_classical_discontinuous_finite_element(1);
    plain_vector D3U(mf_d3u.nb_dof() * 2) ;
    cout << "juste avant la 3e derivation \n";
    getfem::compute_gradient(mf_d2u, mf_d3u, LAP, D3U) ;
    mf_d3u.set_qdim(2) ;
    cout << "derivations : OK ! \n" ;
    

//     mesh_fem mf_du1(mim.linked_mesh(), 1) ;
//     mf_du1.set_classical_finite_element(5);
//     //mf_du.set_qdim(2) ; // peut-etre a enlever ?...
//     plain_vector W(mf_du.nb_dof() * 2), W1(mf_du.nb_dof()), W2(W1) ;
//     getfem::compute_gradient(mf, mf_du, U, W) ;
//     for(unsigned i=0 ; i < mf_du.nb_dof() ; ++i) {
//        W1[i] = W[2*i] ;
//        W2[i] = W[2*i+1] ;}
// 
// 
// 

    /* expression for SIF computation taken from "a finite element
       method for crack growth without remeshing", moes, dolbow &
       belytschko */
//   generic_assembly
//       assem("D_1_nu=data$1(1); D_nu=data$2(1); minus=data$3(1); D=data$4(1); x1=data$5(mdim(#1)); U1=data$6(#1); U2=data$7(#2); q=data$8(#3); U3=data$9(#4); V1=data$10(#5);"
//             "T=U1(i).U2(j).q(k).comp(Hess(#1).Hess(#2).Grad(#3))(i,:,:,j,:,:,k,:);"
//             "t=T;"
//             "V()+=minus(p).D_1_nu(p).t(i,j,i,k,j).x1(k) + minus(p).D_nu(p).t(i,i,j,k,j).x1(k);"
//             "V()+=minus(p).D_1_nu(p).t(i,k,i,j,j).x1(k) + minus(p).D_nu(p).t(j,k,i,i,j).x1(k);"
//             "V()+=D_1_nu(p).t(i,j,i,j,k).x1(k) + D_nu(p).t(i,i,j,j,k).x1(k);");
    generic_assembly
      assem("D_1_nu=data$1(1); D_nu=data$2(1); minus=data$3(1); D=data$4(1); x1=data$5(mdim(#1)); U1=data$6(#1); U2=data$7(#2); q=data$8(#3); U3=data$9(#4); D3U=data$10(#5);"
            "T=U1(i).U2(j).q(k).comp(Hess(#1).Hess(#2).Grad(#3))(i,:,:,j,:,:,k,:);"
            "t=T;"
            "V()+=minus(p).D_1_nu(p).t(i,j,i,k,j).x1(k) + minus(p).D_nu(p).t(i,i,j,k,j).x1(k);"
            "V()+=minus(p).D_1_nu(p).t(i,k,i,j,j).x1(k) + minus(p).D_nu(p).t(j,k,i,i,j).x1(k);"
            "V()+=D_1_nu(p).t(i,j,i,j,k).x1(k) + D_nu(p).t(i,i,j,j,k).x1(k);"
            "T2=U3(i).U1(j).q(k).comp(vBase(#4).Grad(#1).Grad(#3))(i,:,j,:,k,:);"
            "t2=T2;"
            "V()+=D(p).T2(i,j,i).x1(j);"
            "T1=D3U(i).U2(j).q(k).comp(vBase(#5).Grad(#2).Grad(#3))(i,:,j,:,k,:);"
            "t1=T1;"
            "V()+=D(p).T1(i,j,i).x1(j);");
    assem.push_mf(mf);
    assem.push_mf(mf_mode);
    assem.push_mf(mf_q);
    assem.push_mf(mf_mode_3);
    assem.push_mf(mf_d3u);
    assem.push_mi(mim);
    base_vector vD_1_nu(1); vD_1_nu[0] =  D * (1. - nu);
    base_vector vD_nu(1); vD_nu[0] = D * nu ;
    base_vector minus(1); minus[0] = -1. ;
    base_vector vD(1); vD[0] = D ;
    assem.push_data(vD_1_nu); 
    assem.push_data(vD_nu); 
    assem.push_data(minus); 
    assem.push_data(vD);
    assem.push_data(T); // outgoing tangent of the crack
    assem.push_data(U); 
    assem.push_data(U_mode);
    assem.push_data(q);
    assem.push_data(U_mode_3);
    assem.push_data(D3U);
    base_vector V(1);
    assem.push_vec(V);

//     mesh_fem mf_du1(mim.linked_mesh(), 1) ;
//     mf_du1.set_qdim(2) ; // peut-etre a enlever ?...
//     mf_du1.set_classical_finite_element(5);
//     plain_vector DU1(mf_du1.nb_dof() * 2), D1U1(mf_du1.nb_dof()), D2U1(D1U1) ;
//     getfem::compute_gradient(mf, mf_du1, U, DU1) ;
// 
//       mesh_fem mf_du2(mim.linked_mesh(), 1) ;
//       mf_du2.set_classical_finite_element(5);
//       //mf_du.set_qdim(2) ; // peut-etre a enlever ?...
//       plain_vector DU2(mf_du2.nb_dof() * 2), D1U2(mf_du2.nb_dof()), D2U2(D1U2) ;
// 
// //     getfem::compute_gradient(mf, mf_du2, U_mode, DU2) ;
// 
//     for(unsigned i=0 ; i < mf_du1.nb_dof() ; ++i) {
//        D1U1[i] = DU1[2*i] ;
//        D2U1[i] = DU1[2*i+1] ;}
// //     for(unsigned i=0 ; i < mf_du2.nb_dof() ; ++i) {
// //        D1U2[i] = DU2[2*i] ;
// //        D2U2[i] = DU2[2*i+1] ;}
// //getfem::interpolation(mf1, mf2, U, V) ;
//     cout << "interpolations OK \n" ;
// 
//     generic_assembly
//             assem("D_1_nu=data$1(1); D_nu=data$2(1); x1=data$3(mdim(#1)); U1=data$4(#1); U2=data$5(#2); q=data$6(#3); D1U1=data$7(#4); D2U1=data$8(#4); D1U2=data$9(#5); D2U1=data$10(#5); x2=data$11(mdim(#1));"
//             "T=U1(i).U2(j).q(k).comp(Hess(#1).Hess(#2).Grad(#3))(i,:,:,j,:,:,k,:);"
//             "T1=D1U1(i).U2(j).q(k).comp(Hess(#4).Grad(#2).Grad(#3))(i,:,:,j,:,k,:);"
//             "T2=D2U1(i).U2(j).q(k).comp(Hess(#4).Grad(#2).Grad(#3))(i,:,:,j,:,k,:);"
//             "T3=D1U2(i).U1(j).q(k).comp(Hess(#5).Grad(#1).Grad(#3))(i,:,:,j,:,k,:);"
//             "T4=D2U2(i).U1(j).q(k).comp(Hess(#5).Grad(#1).Grad(#3))(i,:,:,j,:,k,:);"
//             "t=T;"
//             "V()+=D_1_nu(p).t(i,j,i,k,j).x1(k) + D_nu(p).t(i,i,j,k,j).x1(k);"
//             "V()+=D_1_nu(p).t(i,k,i,j,j).x1(k) + D_nu(p).t(j,k,i,i,j).x1(k);"
//             "V()+=D_1_nu(p).t(i,j,i,j,k).x1(k) + D_nu(p).t(i,i,j,j,k).x1(k);"
//             "t1=T1;"
//             "t2=T2;"
//             "t3=T3;"
//             "t4=T4;"
//             "V()+=D_1_nu(p).t1(k,i,i,k).x1(k) + D_nu(p).t1(i,i,k,k).x1(k);"
//             "V()+=D_1_nu(p).t2(k,i,i,k).x2(k) + D_nu(p).t2(i,i,j,k).x1(j).x2(k);");
//     assem.push_mf(mf);
//     assem.push_mf(mf_mode);
//     assem.push_mf(mf_q);
//     assem.push_mf(mf_du1);
//     assem.push_mf(mf_du2);
//     assem.push_mi(mim);
//     base_vector vD_1_nu(1); vD_1_nu[0] = - D * (nu - 1.);
//     base_vector vD_nu(1); vD_nu[0] = - D * nu ;
//     assem.push_data(vD_1_nu); 
//     assem.push_data(vD_nu); 
//     assem.push_data(T); // outgoing tangent of the crack
//     assem.push_data(U); 
//     assem.push_data(U_mode);
//     assem.push_data(q);
//     assem.push_data(D1U1);
//     assem.push_data(D2U1);
//     assem.push_data(D1U2);
//     assem.push_data(D2U2);
//     assem.push_data(N); // Normal to the crack (for x2).
//     base_vector V(1);
//     assem.push_vec(V);

    /* fill with the crack opening mode I or mode II */
    for (unsigned mode = 1; mode <= 2; ++mode) {
      base_vector::iterator it = U_mode.begin();
      base_vector::iterator it_3 = U_mode_3.begin();
      scalar_type coeff=sqrt(2.) * (1. - nu * nu) 
                    / (2. * young_modulus * epsilon * (3. + nu));
      switch(mode) {
        case 1: {
          *it++ = 1.;        /* mode I */
          *it++ = 0.;      /* mode II */
          /* "colonne" 1: ux, colonne 2: uy */
          *it_3++ = 1.;     *it_3++ = 0.; /* - cos(theta/2.) */
          *it_3++ = 0.;     *it_3++ = 1.; /* - sin(3 *theta/2.) */
          *it_3++ = 0.;     *it_3++ = 0.; /* - sin(3 *theta/2.)*/ 
          *it_3++ = 0.;     *it_3++ = 0.; /* cos(theta/2.) */
        } break;
        case 2: {
          *it++ = 0.;
          *it++ = 1.;
          *it_3++ = 0.;     *it_3++ = 0.;
          *it_3++ = 0.;     *it_3++ = 0.;
          *it_3++ = 1.;     *it_3++ = 0.;
          *it_3++ = 0.;     *it_3++ = 1.;
        } break;
      }
      gmm::scale(U_mode, coeff);
      gmm::scale(U_mode_3, coeff);
      //

/*      getfem::compute_gradient(mf, mf_du2, U_mode, DU2) ;
      for(unsigned i=0 ; i < mf_du2.nb_dof() ; ++i) {
       D1U2[i] = DU2[2*i] ;
       D2U2[i] = DU2[2*i+1] ;}*/
      //
      V[0] = 0.;
      cout << "assemblig SIFs ..." << std::flush; 
      double time = gmm::uclock_sec();
      assem.assembly(cvring);
      cout << "done (" << gmm::uclock_sec()-time << " sec)\n";
      scalar_type a = 2. * epsilon * M_PI * (1. + nu) / (3. * young_modulus * (3. + nu) ) ;
      V[0] /= 2. * a ; 
      if (mode == 1) KI = V[0]; else KII = V[0];
    }
  }

} // end of namespace ?

scalar_type young_modulus(scalar_type lambda, scalar_type mu){
  return 4*mu*(lambda + mu)/(lambda+2*mu);
}



void bilaplacian_crack_problem::compute_sif(const plain_vector &U, scalar_type ring_radius){

  cout << "Computing stress intensity factors\n";
  base_node tip; 
  base_small_vector T, N;
  get_crack_tip_and_orientation(ls, tip, T, N);
  cout << "crack tip is : " << tip << ", T=" << T << ", N=" << N << "\n";
  scalar_type E = 3. * (1. - nu * nu) * D / (2. * epsilon * epsilon * epsilon) ;
  cout << "young modulus: " << E << "\n";
  //scalar_type ring_radius = 0.2;
  
  scalar_type KI, KII;
/*  estimate_crack_stress_intensity_factors(ls, mf_u(), U, 
                                          E, 
                                          KI, KII, 1e-2);
  cout << "estimation of crack SIF: " << KI << ", " << KII << "\n";*/
  
  compute_crack_stress_intensity_factors_KL(ls, mim, mf_u(), U, ring_radius, D, nu, E, epsilon, KI, KII);
  cout << "computation of crack SIF: " << KI << ", " << KII << "\n";
}
