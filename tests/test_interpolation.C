#include <getfem_export.h>
#include <getfem_regular_meshes.h>
#include <ftool.h>
#ifdef GETFEM_HAVE_SYS_TIMES
#  include <sys/times.h>
#endif
#include <unistd.h>
#include <iomanip>
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif

using getfem::scalar_type;
using getfem::size_type;
using getfem::dim_type;
using getfem::getfem_mesh;
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
    c2*=(.5 - dal::abs(x[i]-.5));
  }
  z[0] = x[0] + c1;
  for (size_type i=1; i < x.size(); ++i) {
    z[i] = x[i] + c2/10.;
  }
  return z;
}

void build_mesh(getfem_mesh& m, int MESH_TYPE, size_type N, size_type NX, size_type K, bool noised) {
  getfem_mesh mesh;
  base_node org(N); org.fill(0.0);
  std::vector<base_small_vector> vtab(N);
  std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NX);
  for (dim_type i = 0; i < N; i++) { 
    vtab[i] = base_small_vector(N); vtab[i].fill(0.0);
    (vtab[i])[i] = 1. / scalar_type(NX) * 1.;
  }
  switch (MESH_TYPE) {
  case 0 : getfem::parallelepiped_regular_simplex_mesh
      (mesh, N, org, vtab.begin(), ref.begin()); break;
  case 1 : getfem::parallelepiped_regular_mesh
      (mesh, N, org, vtab.begin(), ref.begin()); break;
  case 2 : getfem::parallelepiped_regular_prism_mesh
      (mesh, N, org, vtab.begin(), ref.begin()); break;
  default: DAL_THROW(dal::failure_error, "invalid mesh type\n");
  }
  mesh.optimize_structure();
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
    for (dal::bv_visitor cv(mesh.convex_index()); !cv.finished(); ++cv) {
      if (K == 1) { 
	m.add_convex_by_points(mesh.trans_of_convex(cv), mesh.points_of_convex(cv).begin()); 
      } else {
	std::vector<base_node> pts(pgt->nb_points());
	for (size_type i=0; i < pgt->nb_points(); ++i) {
	  pts[i] = mesh.trans_of_convex(cv)->transform(pgt->convex_ref()->points()[i], 
						       mesh.points_of_convex(cv));
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
      for (size_type i=0; i < N; ++i) { if (dal::abs(P[i]) < 1e-10 || dal::abs(P[i]-1.) < 1e-10) is_border = true; }
      if (!is_border) { 
	P = shake_func(P); 
	for (size_type i=0; i < N; ++i) P[i] += 0.05*(1./NX)*dal::random(double());
      }
    }
  }
}

scalar_type interpolate_check(const mesh_fem &mf1, const mesh_fem& mf2) {
  std::vector<scalar_type> U(mf1.nb_dof()), U2(mf1.nb_dof());
  std::vector<scalar_type> V(mf2.nb_dof());
  for (size_type d=0; d < mf1.nb_dof(); ++d) U[d] = func(mf1.point_of_dof(d));
  getfem::interpolation_solution(mf1,mf2,U,V);
  getfem::interpolation_solution(mf2,mf1,V,U2);
  gmm::add(gmm::scaled(U,-1.),U2);
  return gmm::vect_norminf(U2)/gmm::vect_norminf(U);
}

void test_same_mesh(size_type N, size_type NX, size_type K, size_type Qdim1=1, size_type Qdim2=1) {
  chrono c;
  cout << "  Interpolation on same simplex mesh,  N=" << N << ", NX=" << setw(3) << NX 
       << ", P" << K << "<->P" << K+1 << ":"; cout.flush();
  getfem_mesh m;
  build_mesh(m, 0, N, NX, K, false);
  mesh_fem mf1(m,Qdim1); mf1.set_finite_element(m.convex_index(), getfem::PK_fem(N,K), getfem::exact_simplex_im(N));
  mesh_fem mf2(m,Qdim2); mf2.set_finite_element(m.convex_index(), getfem::PK_fem(N,K+1), getfem::exact_simplex_im(N));
  /* force evaluation of a number of things which are not part of interpolation */
  size_type d = mf1.nb_dof(); d -= mf2.nb_dof();
  double err = 0.;  
  for (int i=0; i<3; ++i) { /* pour amortir/cout de la construction du maillage, et de divers trucs
			       (ça a un gros impact) */
    c.init().tic();
    double err2 = interpolate_check(mf1, mf2); 
    if (i==0) err = err2;
    else if (err != err2) DAL_INTERNAL_ERROR("");
    cout << " " << setw(4) << c.toc().cpu(); cout.flush();
  }
  cout << " seconds/interpolation -- rel.err = " << err << "\n";
  assert(err < 1e-8);
}

void test_different_mesh(size_type N, size_type NX, size_type K) {
  chrono c; c.init();
  cout << "  Interpolation on different linear meshes, N=" << N << ", NX=" << setw(3) << NX 
       << ", P" << K << ":"; cout.flush();
  getfem_mesh m1, m2;
  size_type gK=1;
  build_mesh(m1, 0, N, NX, gK, true);
  build_mesh(m2, 0, N, NX, gK, true);
  mesh_fem mf1(m1); mf1.set_finite_element(m1.convex_index(), getfem::PK_fem(N,K), getfem::exact_simplex_im(N));
  mesh_fem mf2(m2); mf2.set_finite_element(m2.convex_index(), getfem::PK_fem(N,K), getfem::exact_simplex_im(N));    
  /* force evaluation of a number of things which are not part of interpolation */
  size_type d = mf1.nb_dof(); d -= mf2.nb_dof();
  double err = 0;
  for (int i=0; i<3; ++i) { /* pour amortir/cout de la construction du maillage, et de divers trucs
			       (ça a un gros impact) */
    c.init().tic();
    double err2 = interpolate_check(mf1, mf2); 
    if (i==0) err = err2;
    else if (err != err2) DAL_INTERNAL_ERROR("");
    cout << " " << setw(4) << c.toc().cpu(); cout.flush();
  }
  cout << " seconds/interpolation -- rel.err = " << err << "\n";
  //mf1.write_to_file("toto.mf",true);
}

int main(int argc, char *argv[]) {
#ifdef GETFEM_HAVE_FEENABLEEXCEPT /* trap SIGFPE */
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
  if (argc == 2 && strcmp(argv[1],"-quick")==0) quick = true;
  cout << "Testing interpolation..\n";
  test_same_mesh(2,quick ? 20 : 80,1);
  test_same_mesh(2,8,1,2,2);
  test_same_mesh(2,quick ? 10 : 20,4);
  test_same_mesh(3,quick ? 5 : 15,1);
  test_different_mesh(2,quick ? 20 : 80,1);
  test_different_mesh(3,quick ? 8 : 15,1);
}
