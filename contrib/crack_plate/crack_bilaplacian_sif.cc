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
  build_sif_ring_from_mesh_KL(const mesh &m, 
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
    if (in < 8) GMM_WARNING1("looks like the radius is too small...");
    return bv;
  }

  /* return the crack tip in P,
     and the outgoing tangent of the crack in T,
     and the normal in N */
  void get_crack_tip_and_orientation_KL(const level_set &/* ls */,
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
    get_crack_tip_and_orientation_KL(ls, P, T, N);
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
                                              const mesh_im &mim2,
                                              const mesh_im &mim3,
                                              const mesh_fem &mf_pre_u, 
                                              const mesh_fem &mf, 
                                              const VECT &U,
                                              scalar_type ring_radius, 
                                              scalar_type D, scalar_type nu,
                                              scalar_type young_modulus, scalar_type epsilon,
                                              scalar_type &KI, scalar_type &KII) {
    // build mesh and geometrical features
    const mesh &mring = mim.linked_mesh();
    base_node crack_tip;
    base_small_vector T, N;
    get_crack_tip_and_orientation_KL(ls, crack_tip, T, N);

    dal::bit_vector cvring = build_sif_ring_from_mesh_KL(mring, crack_tip, 
                                                      ring_radius);

    // set mesh_fem
    mesh_fem_global_function mf_mode(mring, 1), mf_mode_3(mring, 1);

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

    
    getfem::pfem pf ;

    /* fill the "q" ring field with a C^1 field, equal to 
       1 on the inner boundary, and equal to zero on the outer boundary */
    std::vector<scalar_type> q(mf_pre_u.nb_dof());
    for (unsigned d = 0; d < mf_pre_u.nb_dof(); ++d) {
        // getting a getfem::pdof_description about the curent dof
        const getfem::mesh::ind_cv_ct cvs = mf_pre_u.convex_to_dof(d);
        unsigned cv = cvs[0], ld = unsigned(-1);
        // getting the local index  (ld)
        for (unsigned dd = 0; dd < mf_pre_u.nb_dof_of_element(cv); dd += 1) {
           if (mf_pre_u.ind_dof_of_element(cv)[dd] == d) {
             ld = dd; 
             break;}
        }
        pf = mf_pre_u.fem_of_element(cv);
        if ( pf->dof_types().at(ld) == getfem::lagrange_dof(2)){ 
                // pf-> dof_types gives a std_vector<pdof_desciption>
//        if (ld == 0 || ld == 3 || ld == 6 || ld == 9){
           base_node P = mf_pre_u.point_of_dof(d);
           q[d] = (gmm::vect_dist2(P, crack_tip) > 0.999 * ring_radius) ? 0. : 1.; 
        }
        else
           q[d] = 0. ;
    }

      std::string datafilename_q = "q_field.vtk" ;
      cout << "export q field to " << datafilename_q << "..\n";
      getfem::vtk_export exp(datafilename_q, 1);
      exp.exporting(mf_pre_u); 
      exp.write_point_data(mf_pre_u, q, "q_field_displacement");
      cout << "export done, you can view the data file with (for example)\n"
	   << "mayavi -d " << datafilename_q
	   << " -m BandedSurfaceMap -m Outline -f WarpScalar\n";
 
    base_vector U_mode(mf_mode.nb_dof()); assert(U_mode.size() == 2);
    base_vector U_mode_3(mf_mode_3.nb_dof()); assert(U_mode_3.size() == 8);
    
  

cout << "OK avant assemblage\n" ;
  

// ESSAI : dérivation du champs solution de U3, pour prendre le laplacien
// dans l'assemblage (probleme : on perd de la régularité en route)

//     mesh_fem mf_du(mim.linked_mesh(), 1) ;
//     mf_du.set_classical_finite_element(3);
//     plain_vector GRAD(mf_du.nb_dof() * 2) ;
//     getfem::compute_gradient(mf, mf_du, U, GRAD) ;
//     // separation des 2 composantes 
//     plain_vector D1(mf_du.nb_dof()), D2(mf_du.nb_dof()) ;
//     for (unsigned i=0 ; i < mf_du.nb_dof() ; ++i){
//        D1[i] = GRAD[2*i] ;
//        D2[i] = GRAD[2*i+1] ;
//     }
//     cout << "1er derivation OK \n " ;
// 
//     mesh_fem mf_d2u(mim.linked_mesh(), 1) ;
//     mf_d2u.set_classical_discontinuous_finite_element(2);
//     plain_vector GRAD_U1(mf_d2u.nb_dof() * 2), GRAD_U2(mf_d2u.nb_dof() * 2) ;
//     getfem::compute_gradient(mf_du, mf_d2u, D1, GRAD_U1) ;
//     getfem::compute_gradient(mf_du, mf_d2u, D2, GRAD_U2) ;
//     cout << "2e derivation : OK \n" ; 
//     // calcul du laplacien
//     plain_vector LAP(mf_d2u.nb_dof()) ;
//     cout << "taille GRAD_U1 = " << GRAD_U1.size() << "\n" ;
//     cout << "taille GRAD_U1 = " << GRAD_U1.size() << "\n" ;
//     
//     for (unsigned i=0 ; i < mf_d2u.nb_dof() ; ++i)
//        LAP[i] = GRAD_U1[2*i] + GRAD_U2[2*i+1] ;
// 
//     //dérivation :
//     mesh_fem mf_d3u(mim.linked_mesh(), 1) ;
//     mf_d3u.set_classical_discontinuous_finite_element(1);
//     plain_vector D3U(mf_d3u.nb_dof() * 2) ;
//     cout << "juste avant la 3e derivation \n";
//     getfem::compute_gradient(mf_d2u, mf_d3u, LAP, D3U) ;
//     mf_d3u.set_qdim(2) ;
//     cout << "derivations : OK ! \n" ;
    

// 
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
//     generic_assembly
//       assem("D_1_nu=data$1(1); D_nu=data$2(1); minus=data$3(1); D=data$4(1); x1=data$5(mdim(#1)); U1=data$6(#1); U2=data$7(#2); q=data$8(#3); U3=data$9(#4); D3U=data$10(#5);"
//             "T=U1(i).U2(j).q(k).comp(Hess(#1).Hess(#2).Grad(#3))(i,:,:,j,:,:,k,:);"
//             "t=T;"
//             "V()+=minus(p).D_1_nu(p).t(i,j,i,k,j).x1(k) + minus(p).D_nu(p).t(i,i,j,k,j).x1(k);"
//             "V()+=minus(p).D_1_nu(p).t(i,k,i,j,j).x1(k) + minus(p).D_nu(p).t(j,k,i,i,j).x1(k);"
//             "V()+=D_1_nu(p).t(i,j,i,j,k).x1(k) + D_nu(p).t(i,i,j,j,k).x1(k);"
//             "T2=U3(i).U1(j).q(k).comp(vBase(#4).Grad(#1).Grad(#3))(i,:,j,:,k,:);"
//             "t2=T2;"
//             "V()+=D(p).T2(i,j,i).x1(j);"
//             "T1=D3U(i).U2(j).q(k).comp(vBase(#5).Grad(#2).Grad(#3))(i,:,j,:,k,:);"
//             "t1=T1;"
//             "V()+=D(p).T1(i,j,i).x1(j);");
    generic_assembly
      assem("D_1_nu=data$1(1); D_nu=data$2(1); minus=data$3(1); D=data$4(1); x1=data$5(mdim(#1)); U1=data$6(#1); U2=data$7(#2); q=data$8(#3); U3=data$9(#4);"
            "T=U1(i).U2(j).q(k).comp(Hess(#1).Hess(#2).Grad(#3))(i,:,:,j,:,:,k,:);"
            "t=T;"
            "V()+=minus(p).D_1_nu(p).t(i,j,i,k,j).x1(k) + minus(p).D_nu(p).t(i,i,j,k,j).x1(k);"
            "V()+=minus(p).D_1_nu(p).t(i,k,i,j,j).x1(k) + minus(p).D_nu(p).t(j,k,i,i,j).x1(k);"
            "V()+=D_1_nu(p).t(i,j,i,j,k).x1(k) + D_nu(p).t(i,i,j,j,k).x1(k);"
            "T2=U3(i).U1(j).q(k).comp(vBase(#4).Grad(#1).Grad(#3))(i,:,j,:,k,:);"
            "t2=T2;"
            "V()+=D(p).T2(i,j,i).x1(j);"
            "T11=U1(i).U2(j).q(k).comp(Hess(#1).Grad(#2).Hess(#3))(i,:,:,j,:,k,:,:);"
            "t11=T11;"
            "V()+=minus(p).D(p).t11(i,i,j,k,k).x1(j);"
            "V()+=minus(p).D(p).t(i,i,j,k,k).x1(j);");
    assem.push_mf(mf);
    assem.push_mf(mf_mode);
    assem.push_mf(mf_pre_u);
    assem.push_mf(mf_mode_3);
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
    base_vector V(1);
    assem.push_vec(V);

//             "TT2=q(k).comp(%2,Grad(#3))(k,:);"
// 	    "tt2=TT2;" 
// 	    "V()+=D(p).tt2(k).x2(k);"

    generic_assembly assem2("U1=data$1(#1);"
                            "U2=data$2(#2);"
			    "q=data$3(#3);"
                            "D=data$4(1);"
			    "x1=data$5(mdim(#1));"
                            "x2=data$6(mdim(#1));"
			    "minus=data$7(1);"
	    "T=U1(i).U2(j).q(k).comp(Hess(#1).Grad(#2).Grad(#3))(i,:,:,j,:,k,:);"
	    "t=T;"
	    "V()+=minus(p).D(p).t(i,i,j,k).x1(j).x2(k);") ;
    assem2.push_mf(mf);
    assem2.push_mf(mf_mode);
    assem2.push_mf(mf_pre_u);
    assem2.push_mi(mim2);
    assem2.push_data(U); 
    assem2.push_data(U_mode);
    assem2.push_data(q);
    assem2.push_data(vD); 
    assem2.push_data(T); // outgoing tangent of the crack 
    assem2.push_data(N); // Upward normal to the crack
    assem2.push_data(minus);
    base_vector W2(1);
    assem2.push_vec(W2);
    
    generic_assembly assem3("U1=data$1(#1);"
                            "U2=data$2(#2);"
			    "q=data$3(#3);"
                            "D=data$4(1);"
			    "x1=data$5(mdim(#1));"
                            "x2=data$6(mdim(#1));"
			    "minus=data$7(1);"
	    "T=U1(i).U2(j).q(k).comp(Hess(#1).Grad(#2).Grad(#3))(i,:,:,j,:,k,:);"
	    "t=T;"
	    "V()+=D(p).t(i,i,j,k).x1(j).x2(k);") ;
    assem3.push_mf(mf);
    assem3.push_mf(mf_mode);
    assem3.push_mf(mf_pre_u);
    assem3.push_mi(mim3);
    assem3.push_data(U); 
    assem3.push_data(U_mode);
    assem3.push_data(q);
    assem3.push_data(vD); 
    assem3.push_data(T); // outgoing tangent of the crack 
    assem3.push_data(N); // Upward normal to the crack
    assem3.push_data(minus);
    base_vector W3(1);
    assem3.push_vec(W3);
    
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
          *it_3++ = 1.;     *it_3++ = 0.; /* - cos(3 *theta/2.) */
          *it_3++ = 0.;     *it_3++ = 1.; /* - sin(3 *theta/2.) */
          *it_3++ = 0.;     *it_3++ = 0.; /* - sin(3 *theta/2.) */
          *it_3++ = 0.;     *it_3++ = 0.; /*   cos(3 *theta/2.) */
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
      W2[0] = 0. ;
      W3[0] = 0. ;
      assem2.assembly(cvring);
      assem3.assembly(cvring);
      cout << "done (" << gmm::uclock_sec()-time << " sec)\n";
      scalar_type a = 2. * epsilon * M_PI * (1. + nu) / (3. * young_modulus * (3. + nu) ) ;
      V[0] /= 2. * a ; 
      W2[0] /= 2. * a ;
      W3[0] /= 2. * a ;
      scalar_type W = W2[0] + W3[0] ;   
      if ( mode == 2 ) {
         cout << "mode " << mode << " ;\n  ";
	 cout << "W="  << W     << "\n" ;
         cout << "W2=" << W2[0] << "\n" ;
	 cout << "W3=" << W3[0] << "\n" ;
	 }
      if (mode == 1) KI = V[0]  ; else KII = V[0] ;
    }
  }

} // end of namespace ?

// scalar_type young_modulus(scalar_type lambda, scalar_type mu){
//   return 4*mu*(lambda + mu)/(lambda+2*mu);
// }



void bilaplacian_crack_problem::compute_sif(const plain_vector &U, scalar_type ring_radius){

  cout << "Computing stress intensity factors\n";
  base_node tip; 
  base_small_vector T, N;
  get_crack_tip_and_orientation_KL(ls, tip, T, N);
  cout << "crack tip is : " << tip << ", T=" << T << ", N=" << N << "\n";
  scalar_type E = 3. * (1. - nu * nu) * D / (2. * epsilon * epsilon * epsilon) ;
  cout << "young modulus: " << E << "\n";
  //scalar_type ring_radius = 0.2;
  
  scalar_type KI, KII ;
/*  estimate_crack_stress_intensity_factors(ls, mf_u(), U, 
                                          E, 
                                          KI, KII, 1e-2);
  cout << "estimation of crack SIF: " << KI << ", " << KII << "\n";*/
    getfem::level_set ls2(mesh, 1, true) ;
    scalar_type x, y ;
    for (size_type d = 0; d < ls2.get_mesh_fem().nb_dof(); ++d) {
      x = ls2.get_mesh_fem().point_of_dof(d)[0];
      y = ls2.get_mesh_fem().point_of_dof(d)[1];
      ls2.values(0)[d] = y + 1e-6 ;
      ls2.values(1)[d] = x ;
      }
    getfem::mesh_level_set mls2(mesh);
    mls2.add_level_set(ls2);
    mls2.adapt() ;
    
    getfem::level_set ls3(mesh, 1, true) ;
    for (size_type d = 0; d < ls3.get_mesh_fem().nb_dof(); ++d) {
      x = ls3.get_mesh_fem().point_of_dof(d)[0];
      y = ls3.get_mesh_fem().point_of_dof(d)[1];
      ls3.values(0)[d] = y - 1e-6 ;
      ls3.values(1)[d] = x ;
      }  
    getfem::mesh_level_set mls3(mesh);
    mls3.add_level_set(ls3);
    mls3.adapt() ;
    
    
    cout << "initialisations ls2, ls3, mls2, mls3 : OK\n" ;
    
    int where = PARAM.int_value("WHERE") ;
// 1	getfem::mesh_im_level_set::INTEGRATE_INSIDE (integrate over p(x)<0),
// 2	getfem::mesh_im_level_set::INTEGRATE_OUTSIDE (integrate over p(x)>0),
// 3	getfem::mesh_im_level_set::INTEGRATE_ALL,
// 4	getfem::mesh_im_level_set::INTEGRATE_BOUNDARY (integrate over p(x)=0 and s(x) <= 0)
   switch (where){
   case 3: {
   where = getfem::mesh_im_level_set::INTEGRATE_ALL ;}
   break ;
   case 4:{
   where = getfem::mesh_im_level_set::INTEGRATE_BOUNDARY ;}
   break;
   }
   getfem::mesh_im_level_set mim2(mls2, where) ;
   getfem::mesh_im_level_set mim3(mls3, where) ;
   
   
    if (!PARAM.int_value("MIXED_ELEMENTS")){ 
        // Must put 2D methods.
        std::string INTEGRATION = PARAM.string_value("INTEGRATION_LINE",
   					       "Name of integration method");
        std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION_LINE",
    					 "Name of simplex integration method");  
        //must put a method defined for quadrangles.
        //std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION_LINE");
        getfem::pintegration_method ppi = 
        getfem::int_method_descriptor(INTEGRATION);
        getfem::pintegration_method sppi = 
        getfem::int_method_descriptor(SIMPLEX_INTEGRATION);
        //getfem::pintegration_method sing_ppi = (SINGULAR_INTEGRATION.size() ?
        //  getfem::int_method_descriptor(SINGULAR_INTEGRATION) : 0);
        mim2.set_simplex_im(sppi, 0);   // No need to specify a singular sing integration,
        mim3.set_simplex_im(sppi, 0);   // since the J-integral doesn't meet the crack tip.
        mim2.set_integration_method(mesh.convex_index(), ppi);
        mim3.set_integration_method(mesh.convex_index(), ppi);
    }
    else{
        // Must put 2D methods.
        std::string TRI_INTEGRATION = PARAM.string_value("TRI_INTEGRATION_LINE",
    					       "Name of integration method");
        std::string TRI_SIMPLEX_INTEGRATION = PARAM.string_value("TRI_SIMPLEX_INTEGRATION_LINE",
    					 "Name of simplex integration method");  
        getfem::pintegration_method tri_ppi = 
        getfem::int_method_descriptor(TRI_INTEGRATION);
        getfem::pintegration_method tri_sppi = 
        getfem::int_method_descriptor(TRI_SIMPLEX_INTEGRATION);
        std::string QUAD_INTEGRATION = PARAM.string_value("QUAD_INTEGRATION_LINE",
    					       "Name of integration method");
        std::string QUAD_SIMPLEX_INTEGRATION = PARAM.string_value("QUAD_SIMPLEX_INTEGRATION_LINE",
    					 "Name of simplex integration method");  
        getfem::pintegration_method quad_ppi = 
        getfem::int_method_descriptor(QUAD_INTEGRATION);
        getfem::pintegration_method quad_sppi = 
        getfem::int_method_descriptor(QUAD_SIMPLEX_INTEGRATION);
	mim2.set_simplex_im(tri_sppi, 0);   // No need to specify a singular sing integration,
        mim3.set_simplex_im(tri_sppi, 0);   // since the J-integral doesn't meet the crack tip.
	dal::bit_vector quad_among_cvx, tri_among_cvx ;
        for (dal::bv_visitor i(mesh.convex_index()) ; !i.finished() ; ++i){
          if (mesh.points_of_convex(i).size() == 3)
             tri_among_cvx.add(i) ;
          else if (mesh.points_of_convex(i).size() == 4)
             quad_among_cvx.add(i) ;
          else cout << "WARNING : an element has nor 3 or 4 nodes ! \n" ;
        }
        mim2.set_integration_method(tri_among_cvx, tri_ppi);
        mim2.set_integration_method(quad_among_cvx, quad_ppi);
        mim3.set_integration_method(tri_among_cvx, tri_ppi);
        mim3.set_integration_method(quad_among_cvx, quad_ppi);
    }  
    mim2.adapt() ;
    mim3.adapt() ;
    cout << "mim2.set_integration_method(... -> OK\n";
    
    cout << "done.\nEntering in compute_crack_stress_intensity_factors_KL(... \n" ;
  compute_crack_stress_intensity_factors_KL(ls, mim, mim2, mim3, mf_pre_u, mf_u(), U, ring_radius, D, nu, E, epsilon, KI, KII);
  cout << "computation of crack SIF via J-Integral \nK1_J:" << KI << "\nK2_J:" << KII << "\n";
}

void bilaplacian_crack_problem::exact_sif(scalar_type &K1, scalar_type &K2){

   scalar_type E, nu2 ;
   base_small_vector k(2) ;
   nu2 = nu * nu ;
   E = 3. * (1. - nu2) * D / (2. * epsilon * epsilon * epsilon) ;

   switch (sol_ref) {
      case 0: {
      K1 = exact_sol.U[2] * (3.*(nu - 1.)) + exact_sol.U[3] * (3. * nu + 5.) ;
      K1 *= - ( sqrt(2.) * E * epsilon ) / (4. * (1. - nu2) ) ;
      K2 = (exact_sol.U[0] + 3. * exact_sol.U[1]) ;
      K2 *= - E * epsilon * sqrt(2.) * (3. + nu) / (4. * (1. + nu) * (1. + nu) ) ;
      } break;
      case 1: {
      K1 = 3. * sqrt(PARAM.real_value("CRACK_SEMI_LENGTH")) / (2. * epsilon * epsilon) ;
      K2 = 0. ;
      } break;
   }
   cout << "K1_EXACT:"  << K1  << "\n" ;
   cout << "K2_EXACT:" << K2 << "\n" ;      
}

size_type my_is_global_dof_type(getfem::pdof_description dof){
size_type global_dof = 0 ;
   for (dim_type d = 0; d < 4 ; ++d){
       if (dof == getfem::global_dof(d)) {
          global_dof = 1;
	      }
   }
return global_dof ;
}

void bilaplacian_crack_problem::sif_direct_estimation(const plain_vector &U){
   
   unsigned q = mf_u().get_qdim();
   base_small_vector tab_fic(4) ;
   unsigned cpt = 0;
   assert(enrichment_option > 2) ; // global singularities not available in XFEM 
                                   // standard or without any enrichment.

   // Looking for coefficients corresponding to global singularites
   for (size_type d=0; d < mf_u().nb_dof(); d += q) {
      size_type cv = mf_u().first_convex_of_dof(d) ;
      getfem::pfem pf = mf_u().fem_of_element(cv);
      unsigned ld = unsigned(-1);
      for (unsigned dd = 0; dd < mf_u().nb_dof_of_element(cv); dd += q) {
         if (mf_u().ind_dof_of_element(cv)[dd] == d){ 
            ld = dd/q; 
            break;
         }
      }   
      if (ld == unsigned(-1)) {
         cout << "DOF " << d << "NOT FOUND in " << cv << " BUG BUG\n";
      } 
      else {
         if ( my_is_global_dof_type(pf->dof_types().at(ld)) ){
         //cout << "coeff:" << U[d] << "\n" ;
         tab_fic[cpt] = U[d] ;
         cpt += 1 ;
         }
      }
   }

   scalar_type E, nu2 ;
   base_small_vector k(2) ;
   nu2 = nu * nu ;
   E = 3. * (1. - nu2) * D / (2. * epsilon * epsilon * epsilon) ;

   // Calculation of SIF KI and KII

   if (PARAM.int_value("SING_BASE_TYPE") == 0){
      k[0] = tab_fic[2] * (3.*(nu - 1.)) + tab_fic[3] * (3. * nu + 5.) ;
      k[0] *= - ( sqrt(2.) * E * epsilon ) / (4. * (1. - nu2) ) ;
      k[1] = (tab_fic[0] + 3. * tab_fic[1]) ;
      k[1] *= - E * epsilon * sqrt(2.) * (3. + nu) / (4. * (1. + nu) * (1. + nu) ) ;
   }
   if (PARAM.int_value("SING_BASE_TYPE") == 1){
      scalar_type A ;
      A =  sqrt(2) * epsilon * E * (nu + 3.) / (1. - nu2) ;
      k[0] = - tab_fic[0] * A ;
      k[1] =   tab_fic[1] * A ; 
   }
   cout << "K1_DIRECT:" << k[0] << "\n" ;
   cout << "K2_DIRECT:" << k[1] << "\n" ;
    
}























