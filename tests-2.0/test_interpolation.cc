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
#include <getfem/getfem_export.h>
#include <getfem/getfem_export.h>
#include <getfem/getfem_regular_meshes.h>
#ifdef GETFEM_HAVE_SYS_TIMES
#  include <sys/times.h>
#endif
#include <unistd.h>
#include <iomanip>


using getfem::scalar_type;
using getfem::size_type;
using getfem::dim_type;
using getfem::mesh;
using getfem::mesh_fem;
using getfem::pfem;
using getfem::base_node;
using bgeot::base_small_vector;
using std::setw;

bool quick = false;

#ifdef GETFEM_HAVE_SYS_TIMES
struct chrono {
  struct ::tms t;
  ::clock_t t_elapsed;
  float cpu_, elapsed_, system_;
  float nbclocktk;
public:
  chrono() { nbclocktk = ::sysconf(_SC_CLK_TCK); init(); }
  chrono& init() { elapsed_=0; cpu_=0; system_ =0; return *this; }
  void tic() { t_elapsed = ::times(&t); }
  chrono& toc() { 
    struct tms t2; ::clock_t t2_elapsed = ::times(&t2); 
    elapsed_ += (t2_elapsed - t_elapsed) / nbclocktk;
    cpu_     += (t2.tms_utime - t.tms_utime) / nbclocktk;
    system_  += (t2.tms_stime - t.tms_stime) / nbclocktk;
    memcpy(&t, &t2, sizeof(struct tms));
    return *this;
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
  chrono& init() { cpu_=0; return *this; }
  void tic() { t = ::clock()/float(CLOCKS_PER_SEC); }
  chrono& toc() {
    float t2 = ::clock()/float(CLOCKS_PER_SEC);
    cpu_ += t2 - t; t = t2; return *this;
  }
  float cpu() const { return cpu_; }
  float elapsed() const { return cpu_; }
  float system() const { return 0.; }
};
#endif

scalar_type func(const base_node& x) {
  return sin(x[0])*cos(x[1]+x[0]/3.);
}

/* deformation inside a square */
base_node shake_func(const base_node& x) {
  base_node z(x.size());
  scalar_type c1 = 1., c2 = 1.;
  for (size_type i=0; i < x.size(); ++i) {
    c1*=(x[i]*(1.-x[i]));
    c2*=(.5 - gmm::abs(x[i]-.5));
  }
  z[0] = x[0] + c1;
  for (size_type i=1; i < x.size(); ++i) {
    z[i] = x[i] + c2/10.;
  }
  return z;
}

void build_mesh(mesh& m, int MESH_TYPE, size_type dim, size_type N, size_type NX, size_type K, bool noised) {
  mesh msh;
  base_node org(N); gmm::clear(org);
  std::vector<base_small_vector> vtab(N);
  std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NX);
  for (dim_type i = 0; i < N; i++) { 
    vtab[i] = base_small_vector(N); gmm::clear(vtab[i]);
    (vtab[i])[i] = 1. / scalar_type(NX) * 1.;
  }
  switch (MESH_TYPE) {
  case 0 : getfem::parallelepiped_regular_simplex_mesh
      (msh, N, org, vtab.begin(), ref.begin()); break;
  case 1 : getfem::parallelepiped_regular_mesh
      (msh, N, org, vtab.begin(), ref.begin()); break;
  case 2 : getfem::parallelepiped_regular_prism_mesh
      (msh, N, org, vtab.begin(), ref.begin()); break;
  default: GMM_THROW(dal::failure_error, "invalid mesh type\n");
  }
  msh.optimize_structure();
  m.clear();
  /* build a mesh with a geotrans of degree K */
  {
    bgeot::pgeometric_trans pgt;
    switch (MESH_TYPE) {
    case 0: pgt = bgeot::simplex_geotrans(N,K); break;
    case 1: pgt = bgeot::parallelepiped_geotrans(N,K); break;
    case 2: pgt = bgeot::prism_geotrans(N,K); break;
    default: assert(0);
    }
    for (dal::bv_visitor cv(msh.convex_index()); !cv.finished(); ++cv) {
      if (K == 1) { 
	m.add_convex_by_points(msh.trans_of_convex(cv), msh.points_of_convex(cv).begin()); 
      } else {
	std::vector<base_node> pts(pgt->nb_points());
	for (size_type i=0; i < pgt->nb_points(); ++i) {
	  pts[i] = msh.trans_of_convex(cv)->transform(pgt->convex_ref()->points()[i], 
						       msh.points_of_convex(cv));
	}
	m.add_convex_by_points(pgt, pts.begin());
      }
    }
  }

  /* apply a continuous deformation + some noise */
  if (noised) {
    for (dal::bv_visitor ip(m.points().index()); !ip.finished(); ++ip) {
      bool is_border = false;
      base_node& P = m.points()[ip];
      for (size_type i=0; i < N; ++i) { if (gmm::abs(P[i]) < 1e-10 || gmm::abs(P[i]-1.) < 1e-10) is_border = true; }
      if (!is_border) { 
	P = shake_func(P); 
	for (size_type i=0; i < N; ++i) P[i] += 0.05*(1./NX)*gmm::random(double());
      }
    }
  }
  /* add other dimensions to the points */
  assert(dim >= N);
  if (dim > N) {
    getfem::base_matrix T(dim,N);
    for (unsigned i=0; i < N; ++i) T(i,i) = 1;
    for (unsigned i=N; i < dim; ++i) T(i, i%N) = -.5;
    m.transformation(T);
    assert(m.dim() == dim);
  }
}

typedef gmm::col_matrix<gmm::rsvector<scalar_type> > rsc_matrix;
typedef gmm::row_matrix<gmm::rsvector<scalar_type> > rsr_matrix;
typedef gmm::col_matrix<gmm::wsvector<scalar_type> > wsc_matrix;
typedef gmm::row_matrix<gmm::wsvector<scalar_type> > wsr_matrix;

scalar_type interpolate_check(const mesh_fem &mf1, const mesh_fem& mf2, int i, int mat_version) {  
  static std::auto_ptr<rsc_matrix> rsc12, rsc21;
  static std::auto_ptr<rsr_matrix> rsr12, rsr21;
  static std::auto_ptr<wsr_matrix> wsr12, wsr21;
  static std::auto_ptr<wsc_matrix> wsc12, wsc21;

  if (i == 0) {
    switch (mat_version) {
      case 0: break;
      case 1: 
        rsc12.reset(new rsc_matrix(mf2.nb_dof(), mf1.nb_dof())); 
        rsc21.reset(new rsc_matrix(mf1.nb_dof(), mf2.nb_dof())); 
        getfem::interpolation(mf1, mf2, *rsc12); getfem::interpolation(mf2, mf1, *rsc21);
        return 0.;
      case 2: 
        rsr12.reset(new rsr_matrix(mf2.nb_dof(), mf1.nb_dof())); 
        rsr21.reset(new rsr_matrix(mf1.nb_dof(), mf2.nb_dof())); 
        getfem::interpolation(mf1, mf2, *rsr12); getfem::interpolation(mf2, mf1, *rsr21);
        return 0.;
      case 3: 
        wsc12.reset(new wsc_matrix(mf2.nb_dof(), mf1.nb_dof())); 
        wsc21.reset(new wsc_matrix(mf1.nb_dof(), mf2.nb_dof())); 
        getfem::interpolation(mf1, mf2, *wsc12); getfem::interpolation(mf2, mf1, *wsc21);
        return 0.;
      case 4: 
        wsr12.reset(new wsr_matrix(mf2.nb_dof(), mf1.nb_dof())); 
        wsr21.reset(new wsr_matrix(mf1.nb_dof(), mf2.nb_dof())); 
        getfem::interpolation(mf1, mf2, *wsr12); getfem::interpolation(mf2, mf1, *wsr21);
        return 0.;
      default: assert(0);
    }
  }
  std::vector<scalar_type> U(mf1.nb_dof()), U2(mf1.nb_dof());
  std::vector<scalar_type> V(mf2.nb_dof());
  for (size_type d=0; d < mf1.nb_dof(); ++d) U[d] = func(mf1.point_of_dof(d));
  switch (mat_version) {
    case 0: getfem::interpolation(mf1,mf2,U,V); getfem::interpolation(mf2,mf1,V,U2);
      break;
    case 1: gmm::mult(*(rsc12.get()), U, V);  gmm::mult(*(rsc21.get()), V, U2); break;
    case 2: gmm::mult(*(rsr12.get()), U, V);  gmm::mult(*(rsr21.get()), V, U2); break;
    case 3: gmm::mult(*(wsc12.get()), U, V);  gmm::mult(*(wsc21.get()), V, U2); break;
    case 4: gmm::mult(*(wsr12.get()), U, V);  gmm::mult(*(wsr21.get()), V, U2); break;
  }
  gmm::add(gmm::scaled(U,-1.),U2);
  return gmm::vect_norminf(U2)/gmm::vect_norminf(U);
}

void test_same_mesh(int mat_version, size_type N, size_type NX, size_type K, size_type Qdim1=1, size_type Qdim2=1) {
  chrono c;
  cout << "  Same simplex mesh,  N=" << N << ", NX=" << setw(3) << NX 
       << ", P" << K << "<->P" << K+1 << ":"; cout.flush();
  mesh m;
  build_mesh(m, 0, N, N, NX, K, false);
  mesh_fem mf1(m,Qdim1); mf1.set_finite_element(getfem::PK_fem(N,K));
  mesh_fem mf2(m,Qdim2); mf2.set_finite_element(getfem::PK_fem(N,K+1));
  /* force evaluation of a number of things which are not part of interpolation */
  size_type d = mf1.nb_dof(); d -= mf2.nb_dof();
  double err = 0.;  
  for (int i=0; i<3; ++i) { /* pour amortir/cout de la construction du maillage, et de divers trucs
			       (ça a un gros impact) */
    c.init().tic();
    double err2 = interpolate_check(mf1, mf2, i, mat_version); 
    if (i==0 || (mat_version > 0 && i == 1)) err = err2;
    else if (err != err2) GMM_ASSERT1(false, "");
    printf(" %5.1f ", c.toc().cpu()*1000.); //cout << " " << setw(4) << c.toc().cpu(); 
    cout.flush();
  }
  cout << " ms/interpolation -- rel.err = " << err << "\n";
  assert(err < 1e-8);
}


void test_different_mesh(int mat_version, size_type dim, size_type N, size_type NX, size_type K) {
  chrono c; c.init();
  cout << "  Different meshes, dim=" << dim << " N=" << N << ", NX=" << setw(3) << NX 
       << ", P" << K << ":"; cout.flush();
  mesh m1, m2;
  size_type gK=1;
  build_mesh(m1, 0, dim, N, NX, gK, true);
  build_mesh(m2, 0, dim, N, NX, gK, true);
  mesh_fem mf1(m1); mf1.set_finite_element(getfem::PK_fem(N,K));
  mesh_fem mf2(m2); mf2.set_finite_element(getfem::PK_fem(N,K));
  /* force evaluation of a number of things which are not part of interpolation */
  size_type d = mf1.nb_dof(); d -= mf2.nb_dof();
  double err = 0;
  for (int i=0; i<3; ++i) { /* pour amortir/cout de la construction du maillage, et de divers trucs
			       (ça a un gros impact) */
    c.init().tic();
    double err2 = interpolate_check(mf1, mf2, i, mat_version); 
    if (i==0 || (mat_version > 0 && i == 1)) err = err2;
    else GMM_ASSERT1(err == err2, "");
    printf(" %5.1f ", c.toc().cpu()*1000.); //cout << " " << setw(4) << c.toc().cpu(); 
    cout.flush();
  }
  cout << " ms/interpolation -- rel.err = " << err << "\n";
  //mf1.write_to_file("toto.mf",true);
}

void test0() {
  mesh m1, m2;
  std::stringstream ss1("BEGIN POINTS LIST\n"
		       "  POINT  0  2.5  0.6\n"
		       "  POINT  1  5  0\n"
		       "  POINT  2  2.5  1.8\n"
		       "  POINT  3  3.2  1.5\n"
		       "  POINT  4  2.1  1.8\n"
		       "  POINT  5  2.1  3\n"
		       "  POINT  6  2.7  2.4\n"
		       "END POINTS LIST\n\n"
		       "BEGIN MESH STRUCTURE DESCRIPTION\n"
		       "CONVEX 0    \'GT_QK(2,1)\'      0  1  2  3\n"
		       "CONVEX 1    \'GT_QK(2,1)\'      4  2  5  6\n"
		       "END MESH STRUCTURE DESCRIPTION");
  m1.read_from_file(ss1);
  std::stringstream ss2("BEGIN POINTS LIST\n"
			"  POINT  0  3.809523809523809  0.2857142857142857\n"
			"  POINT  1  5  0\n"
			"  POINT  2  3.2  1.5\n"
			"  POINT  3  3.092307692307692  1.361538461538462\n"
			"  POINT  4  2.92  1.62\n"
			"  POINT  5  2.52  2.22\n"
			"  POINT  6  2.6  2.1\n"
			"  POINT  7  2.7  2.4\n"
			"  POINT  8  2.1  2.85\n"
			"  POINT  9  2.1  3\n"
			"END POINTS LIST\n"
			"BEGIN MESH STRUCTURE DESCRIPTION\n"
			"CONVEX 0    \'GT_PK(2,1)\'      0  1  2\n"
			"CONVEX 1    \'GT_PK(2,1)\'      3  0  2\n"
			"CONVEX 2    \'GT_PK(2,1)\'      4  3  2\n"
			"CONVEX 3    \'GT_PK(2,1)\'      5  6  7\n"
			"CONVEX 4    \'GT_PK(2,1)\'      8  5  7\n"
			"CONVEX 5    \'GT_PK(2,1)\'      9  8  7\n"
			"END MESH STRUCTURE DESCRIPTION\n");
  m2.read_from_file(ss2);
  mesh_fem mf1(m1,1), mf2(m2,1);
  mf1.set_finite_element(getfem::fem_descriptor("FEM_QK(2,1)"));
  mf2.set_finite_element(getfem::fem_descriptor("FEM_PK(2,1)"));
  rsc_matrix M(mf2.nb_dof(), mf1.nb_dof());
  getfem::interpolation(mf1, mf2, M);
}

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.


  if (argc == 2 && strcmp(argv[1],"-quick")==0) quick = true;
  try {
    test0();
    for (int mat_version = 0; mat_version < 5; ++mat_version) {
      const char *msg[] = {"Testing interpolation", 
			   "Testing stored interpolator in rsc matrix",
			   "Testing stored interpolator in rsr matrix",
			   "Testing stored interpolator in wsc matrix",
			   "Testing stored interpolator in wsr matrix"};
      cout << msg[mat_version] << "..\n";
      test_same_mesh(mat_version, 2,quick ? 17 : 80,1);
      test_same_mesh(mat_version, 2,quick ? 8 : 20,4);
      test_same_mesh(mat_version, 3,quick ? 5 : 15,1);
      test_different_mesh(mat_version, 2, 2, quick ? 17 : 80,1);
      test_different_mesh(mat_version, 3, 3, quick ? 6 : 15,1);
      if (mat_version == 0) {
	test_different_mesh(0, 2, 1, 100, 2);
	test_different_mesh(0, 3, 1, 500, 1);
	test_different_mesh(0, 3, 2, quick ? 8 : 50, 2);
      }
    }
  }  
  GMM_STANDARD_CATCH_ERROR;
}
