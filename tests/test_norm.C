#include <getfem_assembling.h> /* import assembly methods (and comp. of norms) */
#include <getfem_regular_meshes.h>
#include <getfem_norm.h>

/* try to enable the SIGFPE if something evaluates to a Not-a-number of infinity
 * during computations
 */
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small(dim < 16) vectors */
using bgeot::base_node; /* geometrical nodes (derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */

scalar_type interp_fun(base_node x) {
  scalar_type res = 0.;
  for (size_type i=0; i < x.size(); ++i) {
    res += sin(3*x[i]*(i+1.))+cos(3*x[i]/scalar_type(i+1));
  }
  return res;
}

void test_norm(bgeot::pgeometric_trans pgt, 
	       getfem::pfem pf, 
	       getfem::pintegration_method im, bool noised) {
  getfem::getfem_mesh mesh1, mesh2;
  std::vector<size_type> nsubdiv(pgt->dim());
  std::fill(nsubdiv.begin(),nsubdiv.end(),43);
  getfem::regular_unit_mesh(mesh1, nsubdiv, pgt, noised);
  std::fill(nsubdiv.begin(),nsubdiv.end(),28);
  getfem::regular_unit_mesh(mesh2, nsubdiv, pgt, noised);
  getfem::mesh_fem mf1(mesh1), mf2(mesh2);
  mf1.set_finite_element(mesh1.convex_index(), pf, im);
  mf2.set_finite_element(mesh2.convex_index(), pf, im);
  std::vector<scalar_type> U1(mf1.nb_dof()), U2(mf2.nb_dof());
  getfem::mesh_trans_inv gti(mf1.linked_mesh());
  for (size_type d=0; d < mf1.nb_dof(); ++d) {
    U1[d] = interp_fun(mf1.point_of_dof(d));
    gti.add_point(mf1.point_of_dof(d));
  }
  for (size_type d=0; d < mf2.nb_dof(); ++d)
    U2[d] = interp_fun(mf2.point_of_dof(d));

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
  scalar_type U1_l2 = getfem::asm_L2_norm(mf1,U1);
  scalar_type U1_h1 = getfem::asm_H1_norm(mf1,U1);
  cout << "|U1|_l2=" << U1_l2 << "|U1|_h1=" << U1_h1 << "\n";
  scalar_type U2_l2 = getfem::asm_L2_norm(mf2,U2);
  scalar_type U2_h1 = getfem::asm_H1_norm(mf2,U2);
  cout << "|U2|_l2=" << U2_l2 << "|U2|_h1=" << U2_h1 << "\n";

  scalar_type d_l2, d_h1;
  getfem::solutions_distance(mf2,U2,mf1,U1,im->approx_method(), &d_l2, &d_h1);
  d_h1 = sqrt(dal::sqr(d_l2) + dal::sqr(d_h1));
  cout << "distance_l2 = " << d_l2 << ", distance_h1 = " 
       << d_h1 << "\n";
  getfem::solutions_distance(mf1,U1,mf2,U2,im->approx_method(), &d_l2, &d_h1);
  d_h1 = sqrt(dal::sqr(d_l2) + dal::sqr(d_h1));
  cout << "distance_l2 = " << d_l2 << ", distance_h1 = " 
       << d_h1 << "\n";
  std::fill(U2.begin(), U2.end(),0.);
  scalar_type d2_l2, d2_h1;
  getfem::solutions_distance(mf1,U1,mf2,U2,im->approx_method(), &d2_l2, &d2_h1);  
  d2_h1 = sqrt(dal::sqr(d2_l2) + dal::sqr(d2_h1));
  cout << "norm_l2 = " << d2_l2 << "(diff=" << d2_l2 - U1_l2 
       << "), norm_h1 = " << d2_h1 << " (diff=" << d2_h1 - U1_h1 << ")\n";
  assert(dal::abs(d2_l2 - U1_l2) < 1e-7);
  assert(dal::abs(d2_h1 - U1_h1) < 1e-7);
}

int main(int /*argc*/, char **/*argv*/) {
  dal::exception_callback_debug cb;
  dal::exception_callback::set_exception_callback(&cb); // In order to debug ...

#ifdef GETFEM_HAVE_FEENABLEEXCEPT /* trap SIGFPE */
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif

  try {    
    test_norm(bgeot::geometric_trans_descriptor("GT_PK(2,1)"),
	      getfem::fem_descriptor("FEM_PK(2,2)"),
	      getfem::int_method_descriptor("IM_TRIANGLE(3)"),false);
    /*test_norm(bgeot::geometric_trans_descriptor("GT_PK(3,1)"),
	      getfem::fem_descriptor("FEM_PK(3,1)"),
	      getfem::int_method_descriptor("IM_TETRAHEDRON(1)"),false);*/
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0; 
}
