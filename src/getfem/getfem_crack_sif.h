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
//===========================================================================

/**@file getfem_crack_sif.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>, J. Pommier
   @date July 2007
   @brief crack support functions for computation of SIF (stress intensity factors)
*/

#ifndef GETFEM_CRACK_SIF_H
#define GETFEM_CRACK_SIF_H

#include "getfem/getfem_mesh.h"
#include "getfem/getfem_mesh_fem_global_function.h"
#include "getfem/getfem_assembling_tensors.h"
#include "getfem/getfem_level_set.h"

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
    P.resize(2); P[0] = .5; P[1] = 0;
    T.resize(2); T[0] = 1; T[1] = 0;
    N.resize(2); N[0] = 0; N[1] = 1;
  }


  /* compute with great precision the stress intensity factors using
     integral computed on a ring around the crack tip */
  template <typename VECT>
  void compute_crack_stress_intensity_factors(const level_set &ls,
                                              const mesh_im &mim,
                                              const mesh_fem &mf, 
                                              const VECT &U,
                                              scalar_type ring_radius, 
                                              scalar_type lambda, scalar_type mu, 
                                              scalar_type young_modulus, 
                                              scalar_type &KI, scalar_type &KII) {
    const mesh &mring = mim.linked_mesh();
    mesh_fem_global_function mf_mode(mring, 1);
    mesh_fem mf_q(mring,1);
    
    std::vector<pglobal_function> cfun(4);
    for (unsigned j=0; j < 4; ++j) {
      crack_singular_xy_function *s = 
	new crack_singular_xy_function(j);
      cfun[j] = global_function_on_level_set(ls, *s);
    }
    mf_mode.set_functions(cfun);
    mf_mode.set_qdim(2);
    
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

    base_vector U_mode(mf_mode.nb_dof()); assert(U_mode.size() == 8);
      
    /* expression for SIF computation taken from "a finite element
       method for crack growth without remeshing", moes, dolbow &
       belytschko */

    generic_assembly 
      assem("lambda=data$1(1); mu=data$2(1); x1=data$3(mdim(#1)); U1=data$4(#1); U2=data$5(#2); q=data$6(#3);"
            "t=U1(i).U2(j).q(k).comp(vGrad(#1).vGrad(#2).Grad(#3))(i,:,:,j,:,:,k,:);"
            "e1=(t{1,2,:,:,:}+t{2,1,:,:,:})/2;"
            "e2=(t{:,:,3,4,:}+t{:,:,4,3,:})/2;"
            "e12=(e1{:,:,3,4,:}+e1{:,:,4,3,:})/2;"
            "V()+=2*mu(p).e1(i,j,i,k,j).x1(k) + lambda(p).e1(i,i,j,k,j).x1(k);"
            "V()+=2*mu(p).e2(i,k,i,j,j).x1(k) + lambda(p).e2(j,k,i,i,j).x1(k);"
            "V()+=-2*mu(p).e12(i,j,i,j,k).x1(k) - lambda(p).e12(i,i,j,j,k).x1(k);");
    assem.push_mf(mf);
    assem.push_mf(mf_mode);
    assem.push_mf(mf_q);
    assem.push_mi(mim);
    base_vector vlambda(1); vlambda[0] = lambda;
    base_vector vmu(1); vmu[0] = mu;
    assem.push_data(vlambda); 
    assem.push_data(vmu); 
    assem.push_data(T); // outgoing tangent of the crack
    assem.push_data(U); 
    assem.push_data(U_mode);
    assem.push_data(q);
    base_vector V(1);
    assem.push_vec(V);

    /* fill with the crack opening mode I or mode II */
    for (unsigned mode = 1; mode <= 2; ++mode) {
      base_vector::iterator it = U_mode.begin();
      scalar_type coeff=0.;
      switch(mode) {
        case 1: {
          scalar_type A=2+2*mu/(lambda+2*mu), B=-2*(lambda+mu)/(lambda+2*mu);
          /* "colonne" 1: ux, colonne 2: uy */
          *it++ = 0;       *it++ = A-B; /* sin(theta/2) */
          *it++ = A+B;     *it++ = 0;   /* cos(theta/2) */
          *it++ = -B;      *it++ = 0;   /* sin(theta/2)*sin(theta) */ 
          *it++ = 0;       *it++ = B;   /* cos(theta/2)*cos(theta) */
          coeff = 1/sqrt(2*M_PI);
        } break;
        case 2: {
          scalar_type C1 = (lambda+3*mu)/(lambda+mu); 
          *it++ = C1+2-1;   *it++ = 0;
          *it++ = 0;      *it++ = -(C1-2+1);
          *it++ = 0;      *it++ = 1;
          *it++ = 1;      *it++ = 0;
          coeff = 2*(mu+lambda)/(lambda+2*mu)/sqrt(2*M_PI);
        } break;
      }
      gmm::scale(U_mode, coeff/young_modulus);
      V[0] = 0.;
      cout << "assemblig SIFs ..." << std::flush; 
      double time = gmm::uclock_sec();
      assem.assembly(cvring);
      cout << "done (" << gmm::uclock_sec()-time << " sec)\n";
      V[0] *= young_modulus/2; /* plane stress only, should be E/(2*(1-nu)) for plane strain */
      if (mode == 1) KI = V[0]; else KII = V[0];
    }
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
}

#endif // GETFEM_CRACK_SIF_H
