// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2007-2008 Yves Renard, Julien Pommier.
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
#include <getfem/getfem_assembling.h> /* import assembly methods (and comp. of norms) */
#include <getfem/getfem_assembling.h> /* import assembly methods (and comp. of norms) */
#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_norm.h>

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small(dim < 16) vectors */
using bgeot::base_node; /* geometrical nodes (derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */

scalar_type interp_fun(base_node x) {
  scalar_type res = 0.;
  for (size_type i=0; i < x.size(); ++i) {
    res += sin(3*x[i]*scalar_type(i+1))+cos(3*x[i]/scalar_type(i+1));
  }
  return res;
}

void test_norm(bgeot::pgeometric_trans pgt, 
	       getfem::pfem pf, 
	       getfem::pintegration_method im, bool noised) {
  getfem::mesh mesh1, mesh2;
  std::vector<size_type> nsubdiv(pgt->dim());
  std::fill(nsubdiv.begin(),nsubdiv.end(),13);
  getfem::regular_unit_mesh(mesh1, nsubdiv, pgt, noised);
  std::fill(nsubdiv.begin(),nsubdiv.end(),5);
  getfem::regular_unit_mesh(mesh2, nsubdiv, pgt, noised);
  getfem::mesh_im mim1(mesh1), mim2(mesh2);
  getfem::mesh_fem mf1(mesh1), mf2(mesh2);
  mim1.set_integration_method(mesh1.convex_index(), im);
  mim2.set_integration_method(mesh2.convex_index(), im);  
  mf1.set_finite_element(mesh1.convex_index(), pf);
  mf2.set_finite_element(mesh2.convex_index(), pf);
  std::vector<scalar_type> U1(mf1.nb_dof()), U2(mf2.nb_dof());
  getfem::mesh_trans_inv gti(mf1.linked_mesh());
  for (size_type d=0; d < mf1.nb_dof(); ++d) {
    U1[d] = interp_fun(mf1.point_of_basic_dof(d));
    gti.add_point(mf1.point_of_basic_dof(d));
  }
  for (size_type d=0; d < mf2.nb_dof(); ++d)
    U2[d] = interp_fun(mf2.point_of_basic_dof(d));

  std::vector<scalar_type> V1(mf1.nb_dof());
  getfem::interpolation(mf1,gti,U1,V1);
  scalar_type iterr = gmm::vect_dist2(U1,V1);
  cout << "interpol error: " << iterr << "\n";
  assert(iterr < 1e-10);
  /*
  cout << "interpol.."; cout.flush();
  double t0;
  t0 = ftool::uclock_sec();
  getfem::interpolation(mf1,mf2,U1,U2);
  cout << ftool::uclock_sec() - t0 << " sec -- " << gmm::vect_norm2(U1) << " " << gmm::vect_norm2(U2) << "\n";
  */
  scalar_type U1_l2 = getfem::asm_L2_norm(mim1, mf1, U1);
  scalar_type U1_h1 = getfem::asm_H1_norm(mim1, mf1, U1);
  cout << "|U1|_l2=" << U1_l2 << "|U1|_h1=" << U1_h1 << "\n";
  scalar_type U2_l2 = getfem::asm_L2_norm(mim2, mf2, U2);
  scalar_type U2_h1 = getfem::asm_H1_norm(mim2, mf2, U2);
  cout << "|U2|_l2=" << U2_l2 << "|U2|_h1=" << U2_h1 << "\n";

  scalar_type d_l2, d_h1;
  getfem::solutions_distance(mf2,U2,mf1,U1,im, &d_l2, &d_h1);
  d_h1 = sqrt(gmm::sqr(d_l2) + gmm::sqr(d_h1));
  cout << "distance_l2 = " << d_l2 << ", distance_h1 = " 
       << d_h1 << "\n";
  getfem::solutions_distance(mf1,U1,mf2,U2,im, &d_l2, &d_h1);
  d_h1 = sqrt(gmm::sqr(d_l2) + gmm::sqr(d_h1));
  cout << "distance_l2 = " << d_l2 << ", distance_h1 = " 
       << d_h1 << "\n";
  std::fill(U2.begin(), U2.end(),0.);
  scalar_type d2_l2, d2_h1;
  getfem::solutions_distance(mf1,U1,mf2,U2,im, &d2_l2, &d2_h1);  
  d2_h1 = sqrt(gmm::sqr(d2_l2) + gmm::sqr(d2_h1));
  cout << "norm_l2 = " << d2_l2 << "(diff=" << d2_l2 - U1_l2 
       << "), norm_h1 = " << d2_h1 << " (diff=" << d2_h1 - U1_h1 << ")\n";
  assert(gmm::abs(d2_l2 - U1_l2) < 1e-5);
  assert(gmm::abs(d2_h1 - U1_h1) < 1e-5);
}

int main(int /*argc*/, char **/*argv*/) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {    
    test_norm(bgeot::geometric_trans_descriptor("GT_PK(2,1)"),
	      getfem::fem_descriptor("FEM_PK(2,2)"),
	      getfem::int_method_descriptor("IM_TRIANGLE(3)"),false);
    /*test_norm(bgeot::geometric_trans_descriptor("GT_PK(3,1)"),
	      getfem::fem_descriptor("FEM_PK(3,1)"),
	      getfem::int_method_descriptor("IM_TETRAHEDRON(1)"),false);*/
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
