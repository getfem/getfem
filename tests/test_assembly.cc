/*===========================================================================

 Copyright (C) 2007-2017 Yves Renard, Julien Pommier.

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
template <typename T> std::ostream &operator <<
  (std::ostream &o, const std::vector<T>& m) { gmm::write(o,m); return o; }

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


using std::flush;
#define flushy flush

int fail_cnt = 0;

// static void classical_mesh_fem(getfem::mesh_fem& mf, getfem::short_type K) {
//   for (dal::bv_visitor cv(mf.linked_mesh().convex_index()); !cv.finished();
//        ++cv) {
//     bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
//     mf.set_finite_element(cv, getfem::classical_fem(pgt,K));
//   }
//   //mf.set_classical_finite_element(K,2*K);
// }

// old assembly functions. To compare
namespace getfem { // old assembly procedures with low level generic assembly

  template<typename VEC>
  scalar_type old_asm_L2_norm
  (const mesh_im &mim, const mesh_fem &mf, const VEC &U,
   const mesh_region &rg=mesh_region::all_convexes()) {
    return
      sqrt(old_asm_L2_norm_sqr(mim, mf, U, rg,
			   typename gmm::linalg_traits<VEC>::value_type()));
  }
  
  template<typename VEC, typename T>
  scalar_type old_asm_L2_norm_sqr(const mesh_im &mim, const mesh_fem &mf,
				  const VEC &U, const mesh_region &rg_, T) {
    mesh_region rg(rg_);
    mim.linked_mesh().intersect_with_mpi_region(rg);
    generic_assembly assem;    
    if (mf.get_qdim() == 1)
      assem.set("u=data(#1); V()+=u(i).u(j).comp(Base(#1).Base(#1))(i,j)");
    else
      assem.set("u=data(#1);"
		"V()+=u(i).u(j).comp(vBase(#1).vBase(#1))(i,k,j,k)");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(U);
    std::vector<scalar_type> v(1);
    assem.push_vec(v);
    assem.assembly(rg);
    MPI_SUM_SCALAR(v[0]);
    return v[0];
  }

  template<typename VEC, typename T>
  scalar_type old_asm_L2_norm_sqr(const mesh_im &mim, const mesh_fem &mf,
				  const VEC &U,
				  const mesh_region &rg, std::complex<T>) {
    return asm_L2_norm_sqr(mim, mf,gmm::real_part(U),rg,T()) + 
      asm_L2_norm_sqr(mim, mf,gmm::imag_part(U),rg,T());
  }


  
  template<typename VEC>
  scalar_type old_asm_H1_semi_norm
  (const mesh_im &mim, const mesh_fem &mf, const VEC &U,
   const mesh_region &rg = mesh_region::all_convexes()) {
    typedef typename gmm::linalg_traits<VEC>::value_type T;
    return sqrt(old_asm_H1_semi_norm_sqr(mim, mf, U, rg, T()));
  }

  template<typename VEC, typename T>
  scalar_type old_asm_H1_semi_norm_sqr(const mesh_im &mim, const mesh_fem &mf,
				   const VEC &U, const mesh_region &rg_, T) {

    mesh_region rg(rg_);
    mim.linked_mesh().intersect_with_mpi_region(rg);
    generic_assembly assem;    
    if (mf.get_qdim() == 1)
      assem.set("u=data(#1); V()+=u(i).u(j).comp(Grad(#1).Grad(#1))(i,d,j,d)");
    else
      assem.set("u=data(#1);"
		"V()+=u(i).u(j).comp(vGrad(#1).vGrad(#1))(i,k,d,j,k,d)");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(U);
    std::vector<scalar_type> v(1);
    assem.push_vec(v);
    assem.assembly(rg);
    return MPI_SUM_SCALAR(v[0]);
  }

  

  template<typename VEC, typename T>
  scalar_type old_asm_H1_semi_norm_sqr(const mesh_im &mim, const mesh_fem &mf,
				   const VEC &U,
				   const mesh_region &rg, std::complex<T>) {
    return old_asm_H1_semi_norm_sqr(mim, mf, gmm::real_part(U), rg, T()) + 
      old_asm_H1_semi_norm_sqr(mim, mf, gmm::imag_part(U), rg, T());
  }


  /*
    assembly of a matrix with 1 parameter (real or complex)
    (the most common here for the assembly routines below)
  */
  template <typename MAT, typename VECT>
  void old_asm_real_or_complex_1_param
  (MAT &M, const mesh_im &mim, const mesh_fem &mf_u, const mesh_fem &mf_data,
   const VECT &A, const mesh_region &rg, const char *assembly_description,
   const mesh_fem *mf_mult = 0) {
    old_asm_real_or_complex_1_param_
      (M, mim, mf_u, mf_data, A, rg, assembly_description, mf_mult,
       typename gmm::linalg_traits<VECT>::value_type());
  }

  /* real version */
  template<typename MAT, typename VECT, typename T>
  void old_asm_real_or_complex_1_param_
  (const MAT &M, const mesh_im &mim,  const mesh_fem &mf_u,
   const mesh_fem &mf_data, const VECT &A,  const mesh_region &rg,
   const char *assembly_description, const mesh_fem *mf_mult, T) {
    generic_assembly assem(assembly_description);
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_data);
    if (mf_mult) assem.push_mf(*mf_mult);
    assem.push_data(A);
    assem.push_mat_or_vec(const_cast<MAT&>(M));
    assem.assembly(rg);
  }

  /* complex version */
  template<typename MAT, typename VECT, typename T>
  void old_asm_real_or_complex_1_param_
  (MAT &M, const mesh_im &mim, const mesh_fem &mf_u, const mesh_fem &mf_data,
   const VECT &A, const mesh_region &rg,const char *assembly_description,
   const mesh_fem *mf_mult, std::complex<T>) {
    old_asm_real_or_complex_1_param_(gmm::real_part(M),mim,mf_u,mf_data,
				 gmm::real_part(A),rg,
				 assembly_description, mf_mult, T());
    old_asm_real_or_complex_1_param_(gmm::imag_part(M),mim,mf_u,mf_data,
				 gmm::imag_part(A),rg,
				 assembly_description, mf_mult, T());
  }

  
  template<typename VECT1, typename VECT2>
  void old_asm_source_term
  (const VECT1 &B, const mesh_im &mim, const mesh_fem &mf,
   const mesh_fem &mf_data, const VECT2 &F,
   const mesh_region &rg = mesh_region::all_convexes()) {
    GMM_ASSERT1(mf_data.get_qdim() == 1 ||
		mf_data.get_qdim() == mf.get_qdim(),
		"invalid data mesh fem (same Qdim or Qdim=1 required)");

    const char *st;
    if (mf.get_qdim() == 1)
      st = "F=data(#2); V(#1)+=comp(Base(#1).Base(#2))(:,j).F(j);";
    else if (mf_data.get_qdim() == 1)
      st = "F=data(qdim(#1),#2);"
	"V(#1)+=comp(vBase(#1).Base(#2))(:,i,j).F(i,j);";
    else
      st = "F=data(#2);"
	"V(#1)+=comp(vBase(#1).vBase(#2))(:,i,j,i).F(j);";
    
    old_asm_real_or_complex_1_param(const_cast<VECT1 &>(B),mim,mf,
				    mf_data,F,rg,st);
  }

  template<typename VECT1, typename VECT2>
  void old_asm_normal_source_term(VECT1 &B, const mesh_im &mim,
				  const mesh_fem &mf,
				  const mesh_fem &mf_data, const VECT2 &F,
				  const mesh_region &rg) {
    GMM_ASSERT1(mf_data.get_qdim() == 1 ||
		mf_data.get_qdim() == mf.get_qdim(),
		"invalid data mesh_fem (same Qdim or Qdim=1 required)");

    const char *st;
    if (mf.get_qdim() == 1)
      st = "F=data(mdim(#1),#2);"
	"V(#1)+=comp(Base(#1).Base(#2).Normal())(:,j,k).F(k,j);";
    else if (mf_data.get_qdim() == 1)
      st = "F=data(qdim(#1),mdim(#1),#2);"
	"V(#1)+=comp(vBase(#1).Base(#2).Normal())(:,i,j,k).F(i,k,j);";
    else
      st = "F=data(mdim(#1),#2);"
	"V(#1)+=comp(vBase(#1).vBase(#2).Normal())(:,i,j,i,k).F(k,j);";

    old_asm_real_or_complex_1_param(B, mim, mf, mf_data, F, rg, st);
  }

  template<typename MAT>
  void old_asm_mass_matrix(const MAT &M, const mesh_im &mim,
		       const mesh_fem &mf_u1,
		       const mesh_region &rg = mesh_region::all_convexes()) {
    generic_assembly assem;
    if (mf_u1.get_qdim() == 1)
      assem.set("M(#1,#1)+=sym(comp(Base(#1).Base(#1)))");
    else
      assem.set("M(#1,#1)+=sym(comp(vBase(#1).vBase(#1))(:,i,:,i));");
    assem.push_mi(mim);
    assem.push_mf(mf_u1);
    assem.push_mat(const_cast<MAT &>(M));
    assem.assembly(rg);
  }

  template<typename MAT>
  void old_asm_mass_matrix(const MAT &M, const mesh_im &mim, const mesh_fem &mf_u1,
		       const mesh_fem &mf_u2,
		       const mesh_region &rg = mesh_region::all_convexes()) {
    generic_assembly assem;
    if (mf_u1.get_qdim() == 1 && mf_u2.get_qdim() == 1)
      assem.set("M(#1,#2)+=comp(Base(#1).Base(#2))");
    else if (mf_u1.get_qdim() == 1)
      assem.set("M(#1,#2)+=comp(Base(#1).vBase(#2))(:,:,1);"); // could be i in place of 1
    else if (mf_u2.get_qdim() == 1)
      assem.set("M(#1,#2)+=comp(vBase(#1).Base(#2))(:,1,:);");
    else
      assem.set("M(#1,#2)+=comp(vBase(#1).vBase(#2))(:,i,:,i);");
    assem.push_mi(mim);
    assem.push_mf(mf_u1);
    assem.push_mf(mf_u2);
    assem.push_mat(const_cast<MAT &>(M));
    assem.assembly(rg);
  }

  /** 
      Stiffness matrix for linear elasticity, with Lamé coefficients
      @ingroup asm
  */
  template<class MAT, class VECT>
  void old_asm_stiffness_matrix_for_linear_elasticity
  (const MAT &RM_, const mesh_im &mim, const mesh_fem &mf,
   const mesh_fem &mf_data, const VECT &LAMBDA, const VECT &MU,
   const mesh_region &rg = mesh_region::all_convexes()) {
    MAT &RM = const_cast<MAT &>(RM_);
    GMM_ASSERT1(mf_data.get_qdim() == 1,
		"invalid data mesh fem (Qdim=1 required)");
    
    GMM_ASSERT1(mf.get_qdim() == mf.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");
    /* e = strain tensor,
       M = 2*mu*e(u):e(v) + lambda*tr(e(u))*tr(e(v))
    */
    generic_assembly assem("lambda=data$1(#2); mu=data$2(#2);"
			   "t=comp(vGrad(#1).vGrad(#1).Base(#2));"
			   //"e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}"
			   //"+t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/4;"
			   //"e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:})*0.5;"
			   /*"M(#1,#1)+= sym(2*e(:,i,j,:,i,j,k).mu(k)"
                             " + e(:,i,i,:,j,j,k).lambda(k))");*/
                           "M(#1,#1)+= sym(t(:,i,j,:,i,j,k).mu(k)"
			   "+ t(:,j,i,:,i,j,k).mu(k)"
			   "+ t(:,i,i,:,j,j,k).lambda(k))");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mf_data);
    assem.push_data(LAMBDA);
    assem.push_data(MU);
    assem.push_mat(RM);
    assem.assembly(rg);
  }


  /** 
      Stiffness matrix for linear elasticity, with constant Lamé coefficients
      @ingroup asm
  */
  template<class MAT, class VECT>
  void old_asm_stiffness_matrix_for_homogeneous_linear_elasticity
  (const MAT &RM_, const mesh_im &mim, const mesh_fem &mf,
   const VECT &LAMBDA, const VECT &MU,
   const mesh_region &rg = mesh_region::all_convexes()) {
    MAT &RM = const_cast<MAT &>(RM_);
    GMM_ASSERT1(mf.get_qdim() == mf.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");
    generic_assembly assem("lambda=data$1(1); mu=data$2(1);"
			   "t=comp(vGrad(#1).vGrad(#1));"
                           "M(#1,#1)+= sym(t(:,i,j,:,i,j).mu(1)"
			   "+ t(:,j,i,:,i,j).mu(1)"
			   "+ t(:,i,i,:,j,j).lambda(1))");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(LAMBDA);
    assem.push_data(MU);
    assem.push_mat(RM);
    assem.assembly(rg);
  }

  template<typename MAT>
  void old_asm_stiffness_matrix_for_homogeneous_laplacian
  (const MAT &M_, const mesh_im &mim, const mesh_fem &mf,
   const mesh_region &rg = mesh_region::all_convexes()) {
    MAT &M = const_cast<MAT &>(M_);
    generic_assembly 
      assem("M$1(#1,#1)+=sym(comp(Grad(#1).Grad(#1))(:,i,:,i))");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mat(M);
    assem.assembly(rg);
  }
  
}

typedef enum {
  DO_SCAL_VOLUMIC_SOURCE,
  DO_VEC_VOLUMIC_SOURCE,
  DO_SCAL_MASS_MATRIX,
  DO_VEC_MASS_MATRIX,
  DO_LIN_ELAST, 
  NB_TESTS} t_do_what;

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

struct g_params {
  bgeot::md_param PARAM;

  size_type NX,Ndim;
  int mesh_type;
  int K, K2, Kdata;
  bool do_new, do_old;
  int do_what;
  void init(int argc, char *argv[]);
};

g_params param;

void g_params::init(int argc, char *argv[]) {
  PARAM.add_int_param("NX", 25);
  PARAM.add_int_param("NDIM", 2);
  PARAM.add_int_param("MESH_TYPE", 0);
  PARAM.add_int_param("K", 3);
  PARAM.add_int_param("K2", 3);
  PARAM.add_int_param("KDATA", 1);
  PARAM.add_int_param("BENCH_NEW", 1);
  PARAM.add_int_param("BENCH_OLD", 1);
  PARAM.add_int_param("BENCH_WHAT", -1);

  PARAM.read_command_line(argc, argv);
  NX = PARAM.int_value("NX", "Domaine dimension");
  Ndim = PARAM.int_value("NDIM", "Number of dimensions");
  mesh_type = int(PARAM.int_value("MESH_TYPE", "Mesh type "));
  K = int(PARAM.int_value("K", "Finite element degree"));
  K2 = int(PARAM.int_value("K", "Finite element degree"));
  Kdata = int(PARAM.int_value("KDATA",
                              "Finite element degree for data meshfem"));
  do_new = PARAM.int_value("BENCH_NEW", "bench new assembly routines");
  do_old = PARAM.int_value("BENCH_OLD", "bench old assembly routines");
  do_what = int(PARAM.int_value("BENCH_WHAT",
                                "which test do you want to run?"));
}


#define SCAL_TEST_0(title, expr, mim_, val)                             \
  cout << "\n" << title << endl;                                        \
  workspace.clear_expressions();                                        \
  workspace.add_expression(expr, mim_);                                 \
  ch.init(); ch.tic(); workspace.assembly(0); ch.toc();                 \
  cout << "Elapsed time for new assembly " << ch.elapsed() << endl;     \
  { scalar_type E1 = workspace.assembled_potential();                   \
    cout << "Result=" << E1 << endl;                                    \
    scalar_type error = gmm::abs(E1-val);                               \
    cout << "Error : " << error << endl;                                \
    GMM_ASSERT1(error < 1E-8,						\
                "Error in high or low level generic assembly");         \
  }

#define SCAL_TEST_1(title, expr, mim_, old_asm)                         \
  cout << "\n" << title << endl;                                        \
  workspace.clear_expressions();                                        \
  workspace.add_expression(expr, mim_);                                 \
  ch.init(); ch.tic(); workspace.assembly(0); ch.toc();                 \
  cout << "Elapsed time for new assembly " << ch.elapsed() << endl;     \
  scalar_type E1 = workspace.assembled_potential();                     \
  ch.init(); ch.tic(); scalar_type E2 = old_asm; ch.toc();              \
  cout << "Elapsed time for old assembly " << ch.elapsed() << endl;     \
  scalar_type error = gmm::abs(E1-E2);                                  \
  cout << "Error : " << error << endl;                                  \
  GMM_ASSERT1(error < 1E-7,                                             \
              "Error in high or low level generic assembly");

#define SCAL_TEST_2(expr, mim_)                                         \
  workspace.clear_expressions();                                        \
  workspace.add_expression(expr, mim_);                                 \
  ch.init(); ch.tic(); workspace.assembly(0); ch.toc();                 \
  cout << "Elapsed time for new assembly, alternative expression "      \
       << ch.elapsed() << endl;                                         \
  E1 = workspace.assembled_potential();                                 \
  error = gmm::abs(E1-E2) / (E1+E2);                                    \
  cout << "Error : " << error << endl;                                  \
  GMM_ASSERT1(error < 1E-7,                                             \
              "Error in high or low level generic assembly");

#define VEC_TEST_1(title, ndof, expr, mim_, region, I_, old_asm)        \
  cout << "\n" << title << endl;                                        \
  workspace.clear_expressions();                                        \
  workspace.add_expression(expr, mim_, region);                         \
  ch.init(); ch.tic(); workspace.assembly(1); ch.toc();                 \
  cout << "Elapsed time for new assembly " << ch.elapsed() << endl;     \
  getfem::base_vector V(ndof), V2(ndof);                                \
  ch.init(); ch.tic(); old_asm; ch.toc();                               \
  gmm::copy(V, V2);                                                     \
  cout << "Elapsed time for old assembly " << ch.elapsed() << endl;     \
  gmm::add(gmm::scaled(gmm::sub_vector(workspace.assembled_vector(),    \
                                       I_), scalar_type(-1)), V);       \
  scalar_type norm_error = gmm::vect_norminf(V);                        \
  cout << "Error : " << norm_error << endl;                             \
  GMM_ASSERT1(norm_error < 1E-10,                                       \
              "Error in high or low level generic assembly");

#define VEC_TEST_2(ndof, expr, mim_, region, I_)                        \
  workspace.clear_expressions();                                        \
  workspace.add_expression(expr, mim_, region);                         \
  ch.init(); ch.tic(); workspace.assembly(1); ch.toc();                 \
  cout << "Elapsed time for new assembly, alternative expression "      \
          << ch.elapsed() << endl;                                      \
  gmm::copy(V2, V);                                                     \
  gmm::add(gmm::scaled(gmm::sub_vector(workspace.assembled_vector(),    \
                                       I_), scalar_type(-1)), V);       \
  scalar_type norm_error = gmm::vect_norminf(V);                        \
  cout << "Error : " << norm_error << endl;                             \
  GMM_ASSERT1(norm_error < 1E-10,                                       \
              "Error in high or low level generic assembly");


#define MAT_TEST_1(title, ndof1, ndof2, expr, mim_, I1_, I2_, old_asm)  \
  cout << "\n" << title << endl;                                        \
  workspace.clear_expressions();                                        \
  workspace.add_expression(expr, mim_);                                 \
  ch.init(); ch.tic(); workspace.assembly(2); ch.toc();                 \
  cout << "Elapsed time for new assembly " << ch.elapsed() << endl;     \
  getfem::model_real_sparse_matrix K(ndof1, ndof2), K2(ndof1, ndof2);   \
  ch.init(); ch.tic(); old_asm; ch.toc();                               \
  gmm::copy(K, K2);                                                     \
  cout << "Elapsed time for old assembly " << ch.elapsed() << endl;     \
  gmm::add(gmm::scaled(gmm::sub_matrix(workspace.assembled_matrix(),    \
                                       I1_, I2_), scalar_type(-1)), K); \
  scalar_type norm_error = gmm::mat_norminf(K);                         \
  cout << "Error : " << norm_error << endl;                             \
  GMM_ASSERT1(norm_error < 1E-10,                                       \
              "Error in high or low level generic assembly");


#define MAT_TEST_2(nbdof1, nbdof2, expr, mim_, I1_, I2_)                \
  workspace.clear_expressions();                                        \
  workspace.add_expression(expr, mim_);                                 \
  ch.init(); ch.tic(); workspace.assembly(2);   ch.toc();               \
  cout << "Elapsed time for new assembly, alternative expression "      \
          << ch.elapsed() << endl;                                      \
  gmm::copy(K2, K);                                                     \
  gmm::add(gmm::scaled(gmm::sub_matrix(workspace.assembled_matrix(),    \
                                       I1_, I1_), scalar_type(-1)), K); \
  norm_error = gmm::mat_norminf(K);                                     \
  cout << "Error : " << norm_error << endl;                             \
  GMM_ASSERT1(norm_error < 1E-10,                                       \
              "Error in high or low level generic assembly");



static void test_new_assembly(int N, int NX, int pK) {

    // std::string expr="([1,2;3,4]@[1,2;1,2])(:,2,1,1)(1)+ [1,2;3,4](1,:)(2)"; // should give 4
    // std::string expr="[1,2;3,4]@[1,2;1,2]*[2,3;2,1]/4 + [1,2;3,1]*[1;1](1)"; // should give [4, 8; 12, 13]
    // std::string expr="[1,2;3,a](2,:) + b(:)"; // should give [6, 9]
    // std::string expr="[1,1;1,2,,1,1;1,2;;1,1;1,2,,1,1;1,3](:,:,:,2)";
    // std::string expr="[sin(pi);-2] + Derivative_Norm(Grad_u) + Derivative_Norm(b) + Derivative_sin(pi)*[0;2]";
    // std::string expr = "([1,2;3,4]@[1,2;1,2]).[1;2]";
    // std::string expr = "[u.u; u(1); (u./u)(1); a*Norm(u); c]";
    // std::string expr = "(3*(1*Grad_u)).Grad_Test_u*2 + 0*[1;2].Grad_Test_u + c*Grad_Test_u(1) + [u;1](1)*Test_u";
    // std::string expr = "-(4+(2*3)+2*(1+2))/-(-3+5)"; // should give 8
    // std::string expr="[1,2;3,4]@[1,2;1,2]*(Grad_u@Grad_u)/4 + [1,2;3,1]*[1;1](1)";
    // std::string expr = "Test_u.Test2_u";

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
    m.optimize_structure();

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
    // mf_p.set_qdim(dim_type(N));

    getfem::mesh_im mim(m);
    mim.set_integration_method(m.convex_index(), 4);
    
    getfem::mesh_im mim2(m);
    mim2.set_integration_method(m.convex_index(), 2);

    std::vector<scalar_type> U(mf_u.nb_dof());
    gmm::fill_random(U);
    std::vector<scalar_type> A(mf_u.nb_dof()*N);
    gmm::fill_random(A);
    std::vector<scalar_type> P(mf_p.nb_dof());
    gmm::fill_random(P);
    size_type ndofu = mf_u.nb_dof(), ndofp = mf_p.nb_dof();
    cout << "ndofu = " << ndofu << " ndofp = " << ndofp << endl;
    
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
    cout << "ndofchi = " << ndofchi << endl;
    std::vector<scalar_type> chi(ndofchi);
    gmm::fill_random(chi);
    gmm::sub_interval Ichi(ndofu+ndofp, ndofchi);
    workspace.add_fem_variable("chi", mf_chi, Ichi, chi);
    

    
    chrono ch;

    cout << "\n\nTests in dimension " << N << endl << endl;

    bool all = true;


    if (all) {
      SCAL_TEST_0("Test on function integration 1",
                  "1", mim, 1);
      SCAL_TEST_0("Test on function integration 1",
                  "cos(pi*X(1))", mim, 0);
      cout << "N = " << N << endl;
      SCAL_TEST_0("Test on function integration 2",
                  "cos(pi*X).exp(X*0)", mim, 0);
      SCAL_TEST_0("Test on function integration 2",
                  "-Derivative_sin(pi*X).exp(X*0)", mim, 0);
      auto value = 2.0*double(N) / M_PI;
      SCAL_TEST_0("Test on function integration 2",
                  "sin(pi*X).exp(X*0)", mim, value);
      auto min_value = -value;
      SCAL_TEST_0("Test on function integration 2",
                  "Derivative_cos(pi*X).exp(X*0)", mim, min_value);
      SCAL_TEST_0("Test on function integration 3",
                  "cos(pi*X).Id(meshdim)(:,1)", mim,0);
    }

    if (all) {
      if (N == 2) {
        getfem::ga_define_function("dummyfunc", 1,
                                   "sin(pi*t/2)+2*sqr(t)-[t;t].[t;t]");
        SCAL_TEST_0("Test on user defined functions",
                    "dummyfunc(5)", mim, 1);
        getfem::ga_define_function("dummyfunc2", 1, "cos(pi*t)");
        SCAL_TEST_0("Test on user defined functions",
                    "dummyfunc2(X(1))", mim, 0);
      }
    }



    if (all) {
      SCAL_TEST_1("Test on L2 norm", "u.u", mim,
                  gmm::sqr(getfem::old_asm_L2_norm(mim, mf_u, U)));
      SCAL_TEST_2("Norm_sqr(u)", mim);

      if (N == 2) {
        SCAL_TEST_2("sqr(u(1)) + sqr(u(2))", mim);
        SCAL_TEST_2("u(1)*u(1) + u(2)*u(2)", mim);
        SCAL_TEST_2("[u(2);u(1)].[u(2);u(1)]", mim);
      }
      if (N == 3) {
        SCAL_TEST_2("u(1)*u(1) + u(2)*u(2) + u(3)*u(3)", mim);
        SCAL_TEST_2("[u(2);u(1);u(3)].[u(2);u(1);u(3)]", mim);
      }
    }

    if (all) {
      SCAL_TEST_1("Test on H1 semi-norm", "Grad_u:Grad_u", mim2,
                  gmm::sqr(getfem::old_asm_H1_semi_norm(mim2, mf_u, U)));

      SCAL_TEST_2("Id(meshdim)*Grad_u:Grad_u", mim2);

      if (N == 2) {
        SCAL_TEST_2("Grad_u(1,:).Grad_u(1,:) + Grad_u(2,:).Grad_u(2,:)", mim2);
        SCAL_TEST_2("Grad_u(:,1).Grad_u(:,1) + Grad_u(:,2).Grad_u(:,2)", mim2);
        SCAL_TEST_2("Grad_u(1,1)*Grad_u(1,1) + Grad_u(1,2)*Grad_u(1,2)"
                    "+ Grad_u(2,1)*Grad_u(2,1) + Grad_u(2,2)*Grad_u(2,2)",
                    mim2);
      }
      
      if (N == 3) {
        SCAL_TEST_2("Grad_u(1,:).Grad_u(1,:) + Grad_u(2,:).Grad_u(2,:) +"
                    "Grad_u(3,:).Grad_u(3,:)", mim2);
        SCAL_TEST_2("Grad_u(:,1).Grad_u(:,1) + Grad_u(:,2).Grad_u(:,2) +"
                    "Grad_u(:,3).Grad_u(:,3)", mim2);
      }
    }


    if (all) {
      VEC_TEST_1("Test for source term", ndofu, "u.Test_u", mim, size_type(-1),
                 Iu, getfem::old_asm_source_term(V, mim, mf_u, mf_u, U));

    }

    if (all) {

      {VEC_TEST_1("Test for Neumann term", ndofu, "u.Test_u",
                  mim, NEUMANN_BOUNDARY_NUM,
                  Iu, getfem::old_asm_source_term(V, mim, mf_u, mf_u,
                                              U, NEUMANN_BOUNDARY_NUM));}

      {VEC_TEST_1("Test for Neumann term", ndofu,
                  "(((Reshape(A,meshdim,meshdim))')*Normal).Test_u",
                  mim, NEUMANN_BOUNDARY_NUM,
                  Iu, getfem::old_asm_normal_source_term(V, mim, mf_u, mf_u,
                                              A, NEUMANN_BOUNDARY_NUM));}
      
      if (N == 2)
      {VEC_TEST_1("Test for Neumann term", ndofu,
                  "(A'*Normal).Test_u", mim,
                  NEUMANN_BOUNDARY_NUM,
                  Iu, getfem::old_asm_normal_source_term(V, mim, mf_u, mf_u,
                                                 A, NEUMANN_BOUNDARY_NUM));}
      if (N == 3)
      {VEC_TEST_1("Test for Neumann term", ndofu,
                  "(A'*Normal).Test_u", mim, NEUMANN_BOUNDARY_NUM,
                  Iu, getfem::old_asm_normal_source_term(V, mim, mf_u, mf_u,
                                                 A, NEUMANN_BOUNDARY_NUM));}
    }

    if (all) {
      {VEC_TEST_1("Test for Neumann term with reduced fem", ndofchi,
                  "p*Test_chi", mim, DIRICHLET_BOUNDARY_NUM,
                  Ichi, getfem::old_asm_source_term(V, mim, mf_chi, mf_p,
                                                P, DIRICHLET_BOUNDARY_NUM));}
    }





    if (all) {
      MAT_TEST_1("Test for Mass matrix", ndofu, ndofu, "Test_u.Test2_u", mim,
                 Iu, Iu,  getfem::old_asm_mass_matrix(K, mim, mf_u));
    }

    if (all) {
      MAT_TEST_1("Test for Laplacian stiffness matrix", ndofp, ndofp,
                 "Grad_Test_p:Grad_Test2_p", mim2, Ip, Ip,
                 getfem::old_asm_stiffness_matrix_for_homogeneous_laplacian
                 (K, mim2, mf_p));
      MAT_TEST_2(ndofp, ndofp, "(Grad_p:Grad_p)/2", mim2, Ip, Ip);
      MAT_TEST_2(ndofp, ndofp, "sqr(Norm(Grad_p))/2", mim2, Ip, Ip);
      MAT_TEST_2(ndofp, ndofp, "Norm_sqr(Grad_p)/2", mim2, Ip, Ip);
      if (N == 2) {
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
      if (N == 3) {
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

    if (all) {
      base_vector lambda(1); lambda[0] = 3.0;
      workspace.add_fixed_size_constant("lambda", lambda);
      base_vector mu(1); mu[0] = 2.0;
      workspace.add_fixed_size_constant("mu", mu);
      
      MAT_TEST_1("Test for linear homogeneous elasticity stiffness matrix",
                 ndofu, ndofu, "(lambda*Trace(Grad_Test_u)*Id(qdim(u)) "
                 "+ mu*(Grad_Test_u'+Grad_Test_u)):Grad_Test2_u", mim2,
                 Iu, Iu,
                 getfem::old_asm_stiffness_matrix_for_homogeneous_linear_elasticity
                 (K, mim2, mf_u, lambda, mu));
      MAT_TEST_2(ndofu, ndofu, "lambda*Div_Test_u*Div_Test2_u "
                 "+ mu*(Grad_Test_u'+Grad_Test_u):Grad_Test2_u", mim2, Iu, Iu);
      
      MAT_TEST_2(ndofu, ndofu,
                 "lambda*((Grad_Test2_u@Grad_Test_u):Id(meshdim))"
                 ":Id(meshdim) + mu*(Grad_Test_u'+Grad_Test_u):Grad_Test2_u",
                 mim2, Iu, Iu);
      
      MAT_TEST_2(ndofu, ndofu,
                 "lambda*Id(meshdim)@Id(meshdim)*Grad_Test_u"
                 ":Grad_Test2_u + mu*(Grad_Test_u'+Grad_Test_u):Grad_Test2_u",
                 mim2, Iu, Iu);
      
      MAT_TEST_2(ndofu, ndofu,
                 "lambda*(Id(meshdim)*Id(meshdim))@Id(meshdim)"
                 "*Grad_Test_u:Grad_Test2_u"
                 "+ mu*(Grad_Test_u'+Grad_Test_u):Grad_Test2_u",
                 mim2, Iu, Iu);

      if (N == 2) {
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
      if (N == 3) {
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

    if (all) {
      base_vector lambda2(ndofp, 3.0);
      workspace.add_fem_constant("lambda2", mf_p, lambda2);
      base_vector mu2(ndofp, 2.0);
      workspace.add_fem_constant("mu2", mf_p, mu2);

      MAT_TEST_1("Test for linear non homogeneous elasticity stiffness matrix",
                 ndofu, ndofu, "(lambda2*Trace(Grad_Test_u)*Id(meshdim) "
                 "+ mu2*(Grad_Test_u'+Grad_Test_u)):Grad_Test2_u",
                 mim2, Iu, Iu,
                 getfem::old_asm_stiffness_matrix_for_linear_elasticity
                 (K, mim2, mf_u, mf_p, lambda2, mu2));
    }

}




int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.
  
  test_new_assembly(2, 25, 2);
  test_new_assembly(3, 7, 2);


  // testbug();
  
  param.init(argc,argv);
  std::vector<bool> tests(NB_TESTS, true);
  if (param.do_what >=0 && param.do_what < NB_TESTS) {
    std::fill(tests.begin(),tests.end(),false);
    tests[param.do_what]=true;
  }
  
  
  cerr << "\n\n-----------------------------PERFORMANCE TESTS------------"
       << "---------\n\n";   
  
  cout << "failures: " << fail_cnt << endl;
  return fail_cnt; 
}
