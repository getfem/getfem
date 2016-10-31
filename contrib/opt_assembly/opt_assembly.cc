/*===========================================================================

 Copyright (C) 2007-2016 Yves Renard, Julien Pommier.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

===========================================================================*/
#include "getfem/getfem_assembling.h"
#include "getfem/getfem_generic_assembly.h"
#include "getfem/getfem_export.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_partial_mesh_fem.h"
#include "getfem/getfem_mat_elem.h"
#include "gmm/gmm.h"
#ifdef GETFEM_HAVE_SYS_TIMES
# include <sys/times.h>
#endif
#ifndef _MSC_VER
#include <unistd.h>
#endif
using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;
using std::flush;

using bgeot::base_vector;
using bgeot::base_matrix;
using bgeot::base_small_vector;
using bgeot::base_node;
using bgeot::scalar_type;
using bgeot::size_type;
using bgeot::short_type;
using bgeot::dim_type;

typedef gmm::wsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef std::vector<scalar_type> linalg_vector;



#ifdef HAVE_SYS_TIMES
struct chrono {
  struct ::tms t;
  ::clock_t t_elapsed;
  float cpu_, elapsed_, system_;
  float nbclocktk;
public:
  chrono() { nbclocktk = ::sysconf(_SC_CLK_TCK); }
  void init() { elapsed_=0; cpu_=0; system_ =0; }
  void tic() { t_elapsed = ::times(&t); }
  void toc() { 
    struct tms t2; ::clock_t t2_elapsed = ::times(&t2); 
    elapsed_ += (t2_elapsed - t_elapsed) / nbclocktk;
    cpu_     += (t2.tms_utime - t.tms_utime) / nbclocktk;
    system_  += (t2.tms_stime - t.tms_stime) / nbclocktk;
    memcpy(&t, &t2, sizeof(struct tms));
  }
  float cpu() const { return cpu_; }
  float elapsed() const { return elapsed_; }
  float system() const { return system_; }
};
#else
struct chrono {
  float t,cpu_;
public:
  chrono() { }
  void init() { cpu_=0; }
  void tic() { t = float(::clock())/float(CLOCKS_PER_SEC); }
  void toc() {
    float t2 = float(::clock())/float(CLOCKS_PER_SEC);
    cpu_ += t2 - t; t = t2;
  }
  float cpu() const { return cpu_; }
  float elapsed() const { return cpu_; }
  float system() const { return 0.; }
};
#endif

std::ostream& operator<<(std::ostream& o, const chrono& c) {
  o << "[elapsed=" << int(c.elapsed()*1000) << "ms, cpu="
    << int(c.cpu()*1000) << "ms, system=" << int(c.system()*1000) << "ms]";
  return o;
}





#define VEC_TEST_1(title, ndof, expr, mim_, region, I_, old_asm)        \
  cout << "\n" << title << endl;                                        \
  ch.init(); ch.tic(); workspace.clear_expressions();			\
  workspace.add_expression(expr, mim_, region);                         \
  workspace.assembly(1); ch.toc();					\
  cout << "Elapsed time for new assembly " << ch.elapsed() << endl;     \
  getfem::base_vector V(ndof), V2(ndof);                                \
  ch.init(); ch.tic(); old_asm; ch.toc();                               \
  gmm::copy(V, V2);                                                     \
  cout << "Elapsed time for old assembly " << ch.elapsed() << endl;     \
  gmm::add(gmm::scaled(gmm::sub_vector(workspace.assembled_vector(),    \
                                       I_), scalar_type(-1)), V);       \
  scalar_type norm_error = gmm::vect_norminf(V);                        \
  cout << "Error : " << norm_error << endl;

#define VEC_TEST_2(ndof, expr, mim_, region, I_)                        \
  ch.init(); ch.tic(); workspace.clear_expressions();			\
  workspace.add_expression(expr, mim_, region);                         \
  workspace.assembly(1); ch.toc();					\
  cout << "Elapsed time for new assembly, alternative expression "      \
          << ch.elapsed() << endl;                                      \
  gmm::copy(V2, V);                                                     \
  gmm::add(gmm::scaled(gmm::sub_vector(workspace.assembled_vector(),    \
                                       I_), scalar_type(-1)), V);       \
  scalar_type norm_error = gmm::vect_norminf(V);                        \
  cout << "Error : " << norm_error << endl;

#define VEC_TEST_3(title, ndof, expr, mim_, region)			\
  cout << "\n" << title << endl;                                        \
  ch.init(); ch.tic(); workspace.clear_expressions();			\
  workspace.add_expression(expr, mim_, region);                         \
  workspace.assembly(1); ch.toc();					\
  cout << "Elapsed time for new assembly " << ch.elapsed() << endl;

#define MAT_TEST_1(title, ndof1, ndof2, expr, mim_, I1_, I2_, old_asm)  \
  cout << "\n" << title << endl;					\
  getfem::model_real_sparse_matrix K(ndof1, ndof2), K2(ndof1, ndof2);   \
  getfem::model_real_sparse_matrix K3(I1_.last()+1, I2_.last()+1);	\
  ch.init(); ch.tic(); workspace.clear_expressions();			\
  workspace.set_assembled_matrix(K3);					\
  workspace.add_expression(expr, mim_);					\
  workspace.assembly(2); ch.toc();					\
  cout << "Elapsed time for new assembly " << ch.elapsed() << endl;     \
  ch.init(); ch.tic(); old_asm; ch.toc();				\
  gmm::copy(K, K2);                                                     \
  cout << "Elapsed time for old assembly " << ch.elapsed() << endl;     \
  gmm::add(gmm::scaled(gmm::sub_matrix(K3,I1_, I2_),			\
                       scalar_type(-1)), K);				\
  scalar_type norm_error = gmm::mat_norminf(K);                         \
  cout << "Error : " << norm_error << endl;


#define MAT_TEST_2(nbdof1, nbdof2, expr, mim_, I1_, I2_)                \
  gmm::clear(K3);							\
  ch.init(); ch.tic(); workspace.clear_expressions();			\
  workspace.add_expression(expr, mim_);                                 \
  workspace.assembly(2);   ch.toc();					\
  cout << "Elapsed time for new assembly, alternative expression "      \
          << ch.elapsed() << endl;                                      \
  gmm::copy(K2, K);                                                     \
  gmm::add(gmm::scaled(gmm::sub_matrix(K3, I1_, I1_),			\
		       scalar_type(-1)), K);				\
  norm_error = gmm::mat_norminf(K);                                     \
  cout << "Error : " << norm_error << endl;



static void test_new_assembly(int N, int NX, int pK) {

  
  cout << "\n\n-------------------------------------\n"
       <<     "Tests in dimension " << N << " with P" << pK << " elements"
       <<   "\n-------------------------------------"
       << endl << endl;

    

  getfem::ga_workspace workspace;
  base_vector a(1); a[0] = 3.0;
  workspace.add_fixed_size_constant("a", a);
  base_vector b(2); b[0] = 3.0; b[1] = 6.0;
  workspace.add_fixed_size_constant("b", b);
  // base_vector c(1); c[0] = 1.0;
  // workspace.add_fixed_size_variable("c", gmm::sub_interval(0, 1), c);
  
  getfem::mesh m;
  
  char Ns[5]; sprintf(Ns, "%d", N);
  char Ks[5]; sprintf(Ks, "%d", pK);
  bgeot::pgeometric_trans pgt =
    bgeot::geometric_trans_descriptor
    ((std::string("GT_PK(") + Ns + ",1)").c_str());
  std::vector<size_type> nsubdiv(N, NX);
  getfem::regular_unit_mesh(m, nsubdiv, pgt);
  
  const size_type NEUMANN_BOUNDARY_NUM = 1;
  const size_type DIRICHLET_BOUNDARY_NUM = 2;
  
  base_small_vector Dir(N); Dir[N-1] = 1.0;
  getfem::mesh_region border_faces = getfem::outer_faces_of_mesh(m);
  getfem::mesh_region Neumann_faces
    = getfem::select_faces_of_normal(m, border_faces, Dir, 0.1);
  m.region(NEUMANN_BOUNDARY_NUM) = Neumann_faces;
  m.region(DIRICHLET_BOUNDARY_NUM)
    = getfem::mesh_region::subtract(border_faces, Neumann_faces);
  
  
  getfem::mesh_fem mf_u(m);
  getfem::pfem pf_u = getfem::fem_descriptor
    ((std::string("FEM_PK(") + Ns + "," + Ks + ")").c_str());
  mf_u.set_finite_element(m.convex_index(), pf_u);
  mf_u.set_qdim(dim_type(N));
  
  getfem::mesh_fem mf_p(m);
  getfem::pfem pf_p = getfem::fem_descriptor
    ((std::string("FEM_PK(") + Ns + "," + Ks + ")").c_str());
  mf_p.set_finite_element(m.convex_index(), pf_p);
  
  getfem::mesh_im mim(m);
  mim.set_integration_method(m.convex_index(), dim_type(2*pK));
  
  getfem::mesh_im mim2(m);
  mim2.set_integration_method(m.convex_index(), dim_type(2*pK-2));
  
  std::vector<scalar_type> U(mf_u.nb_dof());
  gmm::fill_random(U);
  std::vector<scalar_type> A(mf_u.nb_dof()*N);
  gmm::fill_random(A);
  std::vector<scalar_type> P(mf_p.nb_dof());
  gmm::fill_random(P);
  size_type ndofu = mf_u.nb_dof(), ndofp = mf_p.nb_dof();
  cout << "ndofu = " << ndofu << " ndofp = " << ndofp;
  
  gmm::sub_interval Iu(0, ndofu);
  gmm::sub_interval Ip(ndofu, ndofp);
  
  workspace.add_fem_variable("u", mf_u, Iu, U);
  workspace.add_fem_constant("A", mf_u, A);
  workspace.add_fem_variable("p", mf_p, Ip, P);
  
  getfem::partial_mesh_fem mf_chi(mf_p);
  dal::bit_vector kept_dof
    = mf_p.basic_dof_on_region(DIRICHLET_BOUNDARY_NUM);
  mf_chi.adapt(kept_dof);
  
  size_type ndofchi = mf_chi.nb_dof();
  cout << " ndofchi = " << ndofchi << endl;
  std::vector<scalar_type> chi(ndofchi);
  gmm::fill_random(chi);
  gmm::sub_interval Ichi(ndofu+ndofp, ndofchi);
  workspace.add_fem_variable("chi", mf_chi, Ichi, chi);
  
  
  
  chrono ch;
  
  bool all = false;
  bool select = true;
  int only_one = 8;

  if (all || select || only_one == 1) {
    VEC_TEST_1("Test for source term", ndofu, "u.Test_u", mim, size_type(-1),
	       Iu, getfem::asm_source_term(V, mim, mf_u, mf_u, U));
    
  }
  
  if (all ||  select || only_one == 2) {
    VEC_TEST_3("Test for nonlinear residual", ndofu, "Det(Grad_u)", mim,
	       size_type(-1));
  }

  if (all || only_one == 3) {
    
    {VEC_TEST_1("Test for Neumann term", ndofu, "u.Test_u",
		mim, NEUMANN_BOUNDARY_NUM,
		Iu, getfem::asm_source_term(V, mim, mf_u, mf_u,
					    U, NEUMANN_BOUNDARY_NUM));}
    
    {VEC_TEST_1("Test for Neumann term", ndofu,
		"(((Reshape(A,meshdim,meshdim))')*Normal).Test_u",
		mim, NEUMANN_BOUNDARY_NUM,
		Iu, getfem::asm_normal_source_term(V, mim, mf_u, mf_u,
						   A, NEUMANN_BOUNDARY_NUM));}
    
    if (N == 2)
      {VEC_TEST_1("Test for Neumann term", ndofu,
                  "(A'*Normal).Test_u", mim,
                  NEUMANN_BOUNDARY_NUM,
                  Iu, getfem::asm_normal_source_term(V, mim, mf_u, mf_u,
						     A, NEUMANN_BOUNDARY_NUM));}
    if (N == 3)
      {VEC_TEST_1("Test for Neumann term", ndofu,
                  "(A'*Normal).Test_u", mim, NEUMANN_BOUNDARY_NUM,
                  Iu, getfem::asm_normal_source_term(V, mim, mf_u, mf_u,
						     A, NEUMANN_BOUNDARY_NUM));}
  }
  
  if (all || only_one == 4) {
    {VEC_TEST_1("Test for Neumann term with reduced fem", ndofchi,
		"p*Test_chi", mim, DIRICHLET_BOUNDARY_NUM,
		Ichi, getfem::asm_source_term(V, mim, mf_chi, mf_p,
					      P, DIRICHLET_BOUNDARY_NUM));}
  }
  
  
  
  
  
 
  
  if (all || select || only_one == 5) {
    MAT_TEST_1("Test for scalar Mass matrix", ndofp, ndofp, "Test_p.Test2_p",
	       mim, Ip, Ip,  getfem::asm_mass_matrix(K, mim, mf_p));
  }
 
  if (all || select || only_one == 6) {
    MAT_TEST_1("Test for vector Mass matrix", ndofu, ndofu, "Test_u.Test2_u",
	       mim, Iu, Iu,  getfem::asm_mass_matrix(K, mim, mf_u));
  }

  if (all || select || only_one == 7) {
    MAT_TEST_1("Test for Laplacian stiffness matrix", ndofp, ndofp,
	       "Grad_Test_p:Grad_Test2_p", mim2, Ip, Ip,
	       getfem::asm_stiffness_matrix_for_homogeneous_laplacian
	       (K, mim2, mf_p));
    // MAT_TEST_2(ndofp, ndofp, "(Grad_p:Grad_p)/2", mim2, Ip, Ip);
    // MAT_TEST_2(ndofp, ndofp, "sqr(Norm(Grad_p))/2", mim2, Ip, Ip);
    MAT_TEST_2(ndofp, ndofp, "Norm_sqr(Grad_p)/2", mim2, Ip, Ip);
    if (false && N == 2) {
      MAT_TEST_2(ndofp, ndofp,
		 "(sqr(Grad_p(1)) + sqr(Grad_p(2)))/2", mim2, Ip, Ip);
      MAT_TEST_2(ndofp, ndofp,
		 "(Grad_p(1)*Grad_p(1) + Grad_p(2)*Grad_p(2))/2",
		 mim2, Ip, Ip);
      MAT_TEST_2(ndofp, ndofp,
		 "([Grad_p(2); Grad_p(1)].[Grad_p(2); Grad_p(1)])/2",
		 mim2, Ip, Ip);
      MAT_TEST_2(ndofp, ndofp, "sqr(Norm([Grad_p(2); Grad_p(1)]))/2",
		 mim2, Ip, Ip);
    }
    if (false && N == 3) {
      MAT_TEST_2(ndofp, ndofp,
		 "(sqr(Grad_p(1)) + sqr(Grad_p(2)) + sqr(Grad_p(3)))/2",
		 mim2, Ip, Ip);
      MAT_TEST_2(ndofp, ndofp,
		 "(Grad_p(1)*Grad_p(1) + Grad_p(2)*Grad_p(2)"
		 "+ Grad_p(3)*Grad_p(3))/2", mim2, Ip, Ip);
      MAT_TEST_2(ndofp, ndofp,
		 "([Grad_p(1); Grad_p(3); Grad_p(2)]."
		 "[Grad_p(1); Grad_p(3); Grad_p(2)])/2",
		 mim2, Ip, Ip);
    }
  }
  
  if (all || select || only_one == 8) {
    base_vector lambda(1); lambda[0] = 3.0;
    workspace.add_fixed_size_constant("lambda", lambda);
    base_vector mu(1); mu[0] = 2.0;
    workspace.add_fixed_size_constant("mu", mu);
    
    MAT_TEST_1("Test for linear homogeneous elasticity stiffness matrix",
	       ndofu, ndofu, "(Div_Test_u*(lambda*Id(qdim(u))) "
	       "+ (2*mu)*Sym(Grad_Test_u)):Grad_Test2_u", mim2,
	       Iu, Iu,
	       getfem::asm_stiffness_matrix_for_homogeneous_linear_elasticity
	       (K, mim2, mf_u, lambda, mu));
    if (all)
      MAT_TEST_2(ndofu, ndofu, "lambda*Div_Test_u*Div_Test2_u "
		 "+ mu*(Grad_Test_u'+Grad_Test_u):Grad_Test2_u", mim2, Iu, Iu);
    
    // MAT_TEST_2(ndofu, ndofu,
    //           "lambda*((Grad_Test2_u@Grad_Test_u):Id(meshdim))"
    //           ":Id(meshdim) + mu*(Grad_Test_u'+Grad_Test_u):Grad_Test2_u",
    //           mim2, Iu, Iu);
    
    // MAT_TEST_2(ndofu, ndofu,
    //           "lambda*Id(meshdim)@Id(meshdim)*Grad_Test_u"
    //           ":Grad_Test2_u + mu*(Grad_Test_u'+Grad_Test_u):Grad_Test2_u",
    //           mim2, Iu, Iu);
    
    // MAT_TEST_2(ndofu, ndofu,
    //           "lambda*(Id(meshdim)*Id(meshdim))@Id(meshdim)"
    //           "*Grad_Test_u:Grad_Test2_u"
    //           "+ mu*(Grad_Test_u'+Grad_Test_u):Grad_Test2_u",
    //           mim2, Iu, Iu);
    
    if (N == 2 && all) {
      MAT_TEST_2(ndofu,ndofu,"lambda*Trace(Grad_Test_u)*Trace(Grad_Test2_u) "
		 "+mu*(Grad_Test_u'(:,1)"
		 "+Grad_Test_u(:,1)):Grad_Test2_u(:,1)"
		 "+mu*(Grad_Test_u'(:,2)"
		 "+Grad_Test_u(:,2)):Grad_Test2_u(:,2) ", mim2, Iu, Iu);
      
      MAT_TEST_2(ndofu,ndofu,"lambda*Trace(Grad_Test_u)*Trace(Grad_Test2_u) "
		 "+mu*(Grad_Test_u'(1,:)"
		 "+Grad_Test_u(1,:)):Grad_Test2_u(1,:)"
		 "+mu*(Grad_Test_u'(2,:)"
		 "+Grad_Test_u(2,:)):Grad_Test2_u(2,:) ", mim2, Iu, Iu);
    }
    if (N == 3 && all) {
      MAT_TEST_2(ndofu,ndofu,"lambda*Trace(Grad_Test_u)*Trace(Grad_Test2_u) "
		 "+mu*(Grad_Test_u'(:,1)"
		 "+Grad_Test_u(:,1)):Grad_Test2_u(:,1)"
		 "+mu*(Grad_Test_u'(:,2)"
		 "+Grad_Test_u(:,2)):Grad_Test2_u(:,2)"
		 "+mu*(Grad_Test_u'(:,3)"
		 "+Grad_Test_u(:,3)):Grad_Test2_u(:,3) ", mim2, Iu, Iu);
      
      MAT_TEST_2(ndofu,ndofu,"lambda*Trace(Grad_Test_u)*Trace(Grad_Test2_u) "
		 "+ mu*(Grad_Test_u'(1,:)"
		 "+Grad_Test_u(1,:)):Grad_Test2_u(1,:)"
		 "+ mu*(Grad_Test_u'(2,:)"
		 "+Grad_Test_u(2,:)):Grad_Test2_u(2,:)"
		 "+mu*(Grad_Test_u'(3,:)"
		 "+Grad_Test_u(3,:)):Grad_Test2_u(3,:) ", mim2, Iu, Iu);
    }
    
  }
  
  if (all || select || only_one == 9) {
    base_vector lambda2(ndofp, 3.0);
    workspace.add_fem_constant("lambda2", mf_p, lambda2);
    base_vector mu2(ndofp, 2.0);
    workspace.add_fem_constant("mu2", mf_p, mu2);
    
    MAT_TEST_1("Test for linear non homogeneous elasticity stiffness matrix",
	       ndofu, ndofu, "(Div_Test_u*(lambda2*Id(meshdim)) "
	       "+ mu2*Sym(Grad_Test_u)):Grad_Test2_u",
	       mim2, Iu, Iu,
	       getfem::asm_stiffness_matrix_for_linear_elasticity
	       (K, mim2, mf_u, mf_p, lambda2, mu2));
  }
  
}


int main(int /* argc */, char * /* argv */[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.
  
  bool all = true;
  int only_one = 5;

  // Mesured times for
  // - new assembly,
  // - old one,
  // - estimate of the storage in sparse matrices part for the new assembly,
  // - global assembly part (assembly instruction),
  // - ga_exec cost (instructions not executed, includes the compilation and 
  //   workspace initialization),
  // - J computation.
  // - Instructions execution except for assembly ones
  //                        new  | old  | sto  | asse | exec | Ins  |
  if (all || only_one == 1) // ndofu = 321602 ndofp = 160801 ndofchi = 1201
    test_new_assembly(2, 400, 1);
  // Vector source term   : 0.20 | 0.66 |
  // Nonlinear residual   : 0.31 |      |
  // Mass (scalar)        : 0.18 | 0.56 | 0.04 | 0.06 | 0.06 | 0.06 |
  // Mass (vector)        : 0.30 | 0.80 | 0.09 | 0.15 | 0.06 | 0.09 |
  // Laplacian            : 0.16 | 0.80 | 0.04 | 0.05 | 0.06 | 0.05 |
  // Homogeneous elas     : 0.30 | 1.88 | 0.08 | 0.13 | 0.06 | 0.09 |
  // Non-homogeneous elast: 0.35 | 2.26 | 0.09 | 0.16 | 0.06 | 0.13 |
  if (all || only_one == 2) // ndofu = 151959 ndofp =  50653 ndofchi = 6553
    test_new_assembly(3, 36, 1);
  // Vector source term   : 0.26 | 0.79 |
  // Nonlinear residual   : 0.54 |      |
  // Mass (scalar)        : 0.21 | 0.58 | 0.05 | 0.08 | 0.08 | 0.06 |
  // Mass (vector)        : 0.68 | 1.37 | 0.17 | 0.27 | 0.08 | 0.33 |
  // Laplacian            : 0.17 | 1.15 | 0.03 | 0.06 | 0.08 | 0.03 |
  // Homogeneous elas     : 0.79 | 4.25 | 0.26 | 0.33 | 0.08 | 0.38 |
  // Non-homogeneous elast: 0.83 | 6.29 | 0.26 | 0.33 | 0.08 | 0.42 |
  if (all || only_one == 3) // ndofu = 321602 ndofp = 160801 ndofchi = 1201
    test_new_assembly(2, 200, 2);
  // Vector source term   : 0.09 | 0.23 |
  // Nonlinear residual   : 0.15 |      |
  // Mass (scalar)        : 0.09 | 0.25 | 0.02 | 0.03 | 0.03 | 0.03 |
  // Mass (vector)        : 0.26 | 0.44 | 0.05 | 0.09 | 0.03 | 0.14 |
  // Laplacian            : 0.09 | 0.37 | 0.02 | 0.03 | 0.03 | 0.03 |
  // Homogeneous elas     : 0.26 | 1.28 | 0.06 | 0.09 | 0.03 | 0.14 |
  // Non-homogeneous elast: 0.28 | 2.38 | 0.06 | 0.10 | 0.03 | 0.16 |
  if (all || only_one == 4) // ndofu = 151959 ndofp =  50653 ndofchi = 6553
    test_new_assembly(3, 18, 2);
  // Vector source term   : 0.13 | 0.23 |
  // Nonlinear residual   : 0.37 |      |
  // Mass (scalar)        : 0.11 | 0.25 | 0.05 | 0.05 | 0.03 | 0.03 |
  // Mass (vector)        : 1.13 | 0.89 | 0.21 | 0.35 | 0.03 | 0.75 |
  // Laplacian            : 0.10 | 0.53 | 0.03 | 0.04 | 0.03 | 0.03 |
  // Homogeneous elas     : 1.66 | 3.35 | 0.59 | 0.73 | 0.03 | 0.90 |
  // Non-homogeneous elast: 1.70 | 9.08 | 0.59 | 0.73 | 0.03 | 0.94 |
  if (all || only_one == 5) // ndofu = 151959 ndofp =  50653 ndofchi = 6553
    test_new_assembly(3, 9, 4);
  // Vector source term   : 0.11 | 0.19 |
  // Nonlinear residual   : 0.36 |      |
  // Mass (scalar)        : 0.51 | 0.34 | 0.09 | 0.16 | 0.01 | 0.33 |
  // Mass (vector)        : 4.37 | 1.31 | 0.41 | 1.27 | 0.01 | 3.10 |
  // Laplacian            : 0.40 | 0.76 | 0.09 | 0.14 | 0.01 | 0.25 |
  // Homogeneous elas     : 8.93 | 5.23 | 0.82 | 1.41 | 0.01 | 7.49 |
  // Non-homogeneous elast: 8.94 | 47.4 | 0.82 | 1.41 | 0.01 | 7.50 |

  // Conclusions :
  // - Deactivation of debug test has no sensible effect.
  // - Compile time of assembly strings is negligible (< 0.0004)
  // - (J, K, B) computation takes half the computational time of the exec part
  // - The optimized instruction call is negligible
  // - For uniform mesh_fem, the resize operations has been suppressed and
  //   the "update pfp" has been isolated in a set  of instruction being
  //   executed only on change of integration method.
  // - The mass matrix is more expansive due to a larger number of Gauss points.
  // - Loop unrolling may have an important impact, especially for high degree 


  // Possible optimizations (focusing on a case)
  // - Optimization of J computation (especially in 2D and standard cases)
  // - Unroll loop in used instructions
  // - Assembly optimization
  // - Detection of the very simples cases where the elementary matrix do not
  //   have to be computed on each element (mass matrix, laplacian ...)
  //   on uniform mesh_fem and mesh_im ?
  // - storage optimization (matrices ...)
  // - Why such a difference between mass matrix and laplacian for 3D and P2 ?

  // Original table (r5370) :
#if 0
  //                        new  | old  | 
  test_new_assembly(2, 400, 1);
  // Vector source term   : 0.76 | 0.77 |
  // Nonlinear residual   : 1.03 |      |
  // Mass (scalar)        : 0.80 | 0.64 |
  // Mass (vector)        : 0.94 | 0.91 |
  // Laplacian            : 0.55 | 0.88 | 
  // Homogeneous elas     : 0.95 | 2.02 | 
  // Non-homogeneous elast: 1.16 | 2.41 |
  test_new_assembly(3, 36, 1);
  // Vector source term   : 0.99 | 1.26 |
  // Nonlinear residual   : 2.82 |      |
  // Mass (scalar)        : 0.92 | 0.97 |
  // Mass (vector)        : 1.70 | 1.80 |
  // Laplacian            : 0.98 | 1.54 |
  // Homogeneous elas     : 2.48 | 5.09 |
  // Non-homogeneous elast: 2.72 | 7.10 | 
  test_new_assembly(2, 200, 2);
  // Vector source term   : 0.33 | 0.26 |
  // Nonlinear residual   : 0.47 |      |
  // Mass (scalar)        : 0.35 | 0.29 |
  // Mass (vector)        : 0.57 | 0.54 |
  // Laplacian            : 0.28 | 0.42 |
  // Homogeneous elas     : 0.74 | 1.42 |
  // Non-homogeneous elast: 0.89 | 2.56 |
  test_new_assembly(3, 18, 2);
  // Vector source term   : 0.49 | 0.28 |
  // Nonlinear residual   : 1.30 |      |
  // Mass (scalar)        : 0.51 | 0.37 |
  // Mass (vector)        : 2.31 | 1.09 |
  // Laplacian            : 0.38 | 0.65 |
  // Homogeneous elas     : 3.35 | 4.13 |
  // Non-homogeneous elast: 3.48 | 10.2 |
  test_new_assembly(3, 9, 4);
  // Vector source term   : 0.40 | 0.20 |
  // Nonlinear residual   : 0.93 |      |
  // Mass (scalar)        : 0.79 | 0.47 |
  // Mass (vector)        : 7.11 | 1.59 |
  // Laplacian            : 0.96 | 0.91 |
  // Homogeneous elas     : 13.4 | 6.61 |
  // Non-homogeneous elast: 13.7 | 49.1 |
#endif

  return 0; 
}
