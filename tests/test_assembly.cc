/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2007-2012 Yves Renard, Julien Pommier.
 
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
 
===========================================================================*/
#include "getfem/getfem_assembling.h"
#include "getfem/getfem_generic_assembly.h"
#include "getfem/getfem_export.h"
#include "getfem/getfem_regular_meshes.h"
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

static void classical_mesh_fem(getfem::mesh_fem& mf, getfem::short_type K) {
  for (dal::bv_visitor cv(mf.linked_mesh().convex_index()); !cv.finished();
       ++cv) {
    bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
    mf.set_finite_element(cv, getfem::classical_fem(pgt,K));
  }
  //mf.set_classical_finite_element(K,2*K);
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

namespace getfem {
  template<class VECT1, class VECT2>
  void old_asm_Neumann_condition(VECT1 &B, const mesh_im &mim,
				 const mesh_fem &mf,
				 size_type boundary, const mesh_fem &mfdata,
				 const VECT2 &F, dim_type N) {
    size_type cv, nbd1, nbd2, f;
    dal::bit_vector nn = mf.convex_index(), nf;
    base_tensor t;
    pfem pf1, pf2, pf1prec = NULL, pf2prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      GMM_ASSERT1(false,
		  "This assembling procedure only works on a single mesh");
  
    for (cv << nn; cv != ST_NIL; cv << nn) {
      nf =
        dal::bit_vector(mf.linked_mesh().region(boundary).faces_of_convex(cv));
      if (nf.card() > 0) {
	pf1 =     mf.fem_of_element(cv); nbd1 = pf1->nb_dof(cv);
	pf2 = mfdata.fem_of_element(cv); nbd2 = pf2->nb_dof(cv);
	pgt = mf.linked_mesh().trans_of_convex(cv);
	pim = mim.int_method_of_element(cv);
	if (pf1prec != pf1 || pf2prec != pf2 || pgtprec!=pgt || pimprec!=pim) {
	  pme = mat_elem_product(mat_elem_base(pf1), mat_elem_base(pf2));
	  pmec = mat_elem(pme, pim, pgt);
	  pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
	}
	for (f << nf; f != ST_NIL; f << nf) {
	  pmec->gen_compute_on_face(t,mf.linked_mesh().points_of_convex(cv),
                                    f, cv);
	  base_tensor::iterator p = t.begin();
	  for (size_type i = 0; i < nbd2; i++)
	    {
	      size_type dof2 = mfdata.ind_basic_dof_of_element(cv)[i];
	      for (size_type j = 0; j < nbd1; j++, ++p)
		{
		  size_type dof1 = mf.ind_basic_dof_of_element(cv)[j];
		  for (size_type k = 0; k < N; k++) {
		    B[dof1*N + k] += F[dof2*N+k]*(*p);
		  }
		}
	    }
	  if (p != t.end()) GMM_ASSERT1(false, "internal error"); 
	}
      }
    }
  }

  template<class VECT1, class VECT2>
  void old_asm_volumic_source_term(VECT1 &B, const mesh_im &mim,
                                   const mesh_fem &mf,
				   const mesh_fem &mfdata,
                                   const VECT2 &F, dim_type N)
  {
    size_type cv, nbd1, nbd2;
    dal::bit_vector nn = mf.convex_index();
    base_tensor t;
    pfem pf1, pf2, pf1prec = NULL, pf2prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      GMM_ASSERT1(false,
		  "This assembling procedure only works on a single mesh");

    for (cv << nn; cv != ST_NIL; cv << nn)
      {
	pf1 =     mf.fem_of_element(cv); nbd1 = pf1->nb_dof(cv);
	pf2 = mfdata.fem_of_element(cv); nbd2 = pf2->nb_dof(cv);
	pgt = mf.linked_mesh().trans_of_convex(cv);
	pim = mim.int_method_of_element(cv);
	if (pf1prec != pf1 || pf2prec != pf2 || pgtprec != pgt ||
            pimprec != pim) {
          pme = mat_elem_product(mat_elem_base(pf1), mat_elem_base(pf2));
          pmec = mat_elem(pme, pim, pgt);
          pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
        }
	pmec->gen_compute(t, mf.linked_mesh().points_of_convex(cv), cv);
	base_tensor::iterator p = t.begin();
	for (size_type i = 0; i < nbd2; i++) {
          size_type dof2 = mfdata.ind_basic_dof_of_element(cv)[i];
          for (size_type j = 0; j < nbd1; j++, ++p) {
            size_type dof1 = mf.ind_basic_dof_of_element(cv)[j];
            for (size_type k = 0; k < N; k++)
              B[dof1*N + k] += F[dof2*N+k]*(*p);
          }
        }
	if (p != t.end()) GMM_ASSERT1(false, "internal error"); 
      }
  }

  template<class MATRM, class MESH_FEM>
  void old_asm_mass_matrix(MATRM &M, const mesh_im &mim, const MESH_FEM &mf1,
			   const MESH_FEM &mf2, dim_type N) {
    size_type cv, nbd1, nbd2;
    dal::bit_vector nn = mf1.convex_index();
    base_tensor t;
    pfem pf1, pf1prec = 0, pf2, pf2prec = 0;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf1.linked_mesh()) != &(mf2.linked_mesh()))
      GMM_ASSERT1(false,
		  "This assembling procedure only works on a single mesh");

    for (cv << nn; cv != ST_NIL; cv << nn) {
	pf1 = mf1.fem_of_element(cv); nbd1 = pf1->nb_dof(cv);
	pf2 = mf2.fem_of_element(cv); nbd2 = pf2->nb_dof(cv);
	pgt = mf1.linked_mesh().trans_of_convex(cv);
	pim = mim.int_method_of_element(cv);
	if (pf1prec != pf1 || pf2prec != pf2 || pgtprec != pgt ||
            pimprec != pim)
	  {
	    pme = mat_elem_product(mat_elem_base(pf1), mat_elem_base(pf2));
	    pmec = mat_elem(pme, pim, pgt);
	    pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
	  }
	pmec->gen_compute(t, mf1.linked_mesh().points_of_convex(cv), cv);

	// cout << "t = " << t << endl;
      
	base_tensor::iterator p = t.begin();
	for (size_type i = 0; i < nbd2; i++) {
	    size_type dof2 = mf2.ind_basic_dof_of_element(cv)[i];
	    // cout << "cv = " << cv << " dof2 = " << dof2 << endl;
	    for (size_type j = 0; j < nbd1; j++, ++p) {
		size_type dof1 = mf1.ind_basic_dof_of_element(cv)[j];
		// cout << "dof1 = " << dof1 << " dof2 = " << dof2 << endl;
		for (size_type k = 0; k < N; k++)
		  M(dof1*N + k, dof2*N + k) += (*p);
	      }
	  }
	if (p != t.end()) GMM_ASSERT1(false, "internal error"); 
      }
  }

  template<class MAT, class VECT>
  void old_asm_stiffness_matrix_for_linear_elasticity
  (MAT &RM, const mesh_im &mim, const mesh_fem &mf, 
   const mesh_fem &mfdata, const VECT &LAMBDA, const VECT &MU) {

    size_type cv, nbd2, N = mf.linked_mesh().dim();
    dal::bit_vector nn = mf.convex_index();
    base_tensor t;
    pfem pf1, pf2, pf1prec = NULL, pf2prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      GMM_ASSERT1(false,
		  "This assembling procedure only works on a single mesh");
  
    for (cv << nn; cv != ST_NIL; cv << nn) {
      pf1 =     mf.fem_of_element(cv); 
      pf2 = mfdata.fem_of_element(cv); nbd2 = pf2->nb_dof(cv);
      pgt = mf.linked_mesh().trans_of_convex(cv);
      pim = mim.int_method_of_element(cv);
      if (pf1prec != pf1 || pf2prec != pf2 || pgtprec != pgt || pimprec != pim)
	{
	  pme = mat_elem_product(mat_elem_product(mat_elem_grad(pf1),
						  mat_elem_grad(pf1)), 
				 mat_elem_base(pf2));
	  pmec = mat_elem(pme, pim, pgt);
	  pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
	}
      pmec->gen_compute(t, mf.linked_mesh().points_of_convex(cv), cv);
      base_tensor::iterator p = t.begin();
      
      size_type nbd = mf.nb_basic_dof_of_element(cv);
      
      for (size_type r = 0; r < nbd2; r++) {
	size_type dof3 = mfdata.ind_basic_dof_of_element(cv)[r];
	for (dim_type l = 0; l < N; l++)
	  for (size_type j = 0; j < nbd; j++) {
	    size_type dof2 = mf.ind_basic_dof_of_element(cv)[j];
	    
	    for (dim_type k = 0; k < N; k++)
	      for (size_type i = 0; i < nbd; i++, ++p) {
		size_type dof1 = mf.ind_basic_dof_of_element(cv)[i];
		
		if (dof1*N + k >= dof2*N + l) {
		  RM(dof1*N + k, dof2*N + l) += LAMBDA[dof3] * (*p);
		  RM(dof2*N + l, dof1*N + k) = RM(dof1*N + k, dof2*N + l);
		}
		
		if (dof1*N + l >= dof2*N + k) {
		  RM(dof1*N + l, dof2*N + k) += MU[dof3] * (*p);
		  RM(dof2*N + k, dof1*N + l) = RM(dof1*N + l, dof2*N + k);
		}
		
		if (l == k && dof1 >= dof2)
		  for (size_type n = 0; n < N; ++n) {
		    RM(dof1*N + n, dof2*N + n) += MU[dof3] * (*p);
		    RM(dof2*N + n, dof1*N + n) = RM(dof1*N + n, dof2*N + n);
		  }
		
	      }
	  }
      }
      if (p != t.end()) GMM_ASSERT1(false, "internal error"); 
    }
  }

} /* namespace getfem */


static void gen_mesh(getfem::mesh& mesh) {
  cout << "Mesh generation, N=" << param.NX << " Ndim=" << param.Ndim << endl;
  base_node org(param.Ndim); gmm::clear(org);
  std::vector<base_small_vector> vtab(param.Ndim);
  std::vector<size_type> ref(param.Ndim);
  std::fill(ref.begin(), ref.end(), param.NX);
  for (size_type i = 0; i < param.Ndim; i++) { 
    vtab[i] = base_small_vector(param.Ndim); gmm::clear(vtab[i]);
    (vtab[i])[i] = 1. / scalar_type(param.NX);
  }
  switch (param.mesh_type) {
  case 0: getfem::parallelepiped_regular_simplex_mesh
      (mesh, dim_type(param.Ndim), org,vtab.begin(), ref.begin()); 
    cerr << mesh.convex_index().card() << " " << param.Ndim
         << "D simplexes generated\n";
    break;
  case 1 : getfem::parallelepiped_regular_mesh
      (mesh, dim_type(param.Ndim), org, vtab.begin(), ref.begin()); 
    cerr << mesh.convex_index().card() << " " << param.Ndim
         << "D parallelepipeds generated\n";
    break;
  case 2 : getfem::parallelepiped_regular_prism_mesh
      (mesh, dim_type(param.Ndim), org, vtab.begin(), ref.begin()); 
    cerr << mesh.convex_index().card() << " " << param.Ndim
         << "D prisms generated\n";
    break;
  default : GMM_ASSERT1(false, "Unknown type of mesh");
  }

  assert(param.NX>2);
  /* un ptit trou dans la liste des convexes ne fait pas de mal */
  mesh.sup_convex(param.NX/2);
  mesh.sup_convex(param.NX/2 + 1);
  mesh.optimize_structure();

  /* bouge un peu les noeuds */
  /*  for (size_type i=0; i < mesh.points().size(); ++i) {
    for (size_type j=0; j < param.Ndim; ++j) {
      float d = ((rand() % 100)-50)/(500.*param.NX);
      mesh.points()[i][j] += d;
    }
    }*/
  for (unsigned cv=0; cv < std::min(mesh.convex_index().card(),
			   param.NX*param.Ndim*param.Ndim*10); cv += 2) {
    mesh.region(1).add(cv, (cv/4) % (param.Ndim > 1 ? 3 : 2)); 
  }
  mesh.region(1).add(0,0);
}

static void init_mesh_fem(getfem::mesh_fem &mf, bool datamf) {
  if (datamf)
    mf.set_classical_finite_element(dim_type(param.Kdata));
  else {
    dal::bit_vector cvlst = mf.linked_mesh().convex_index();
    for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
      if ((cv+1) % 100) {
	mf.set_finite_element(cv,
                              getfem::classical_fem(pgt,short_type(param.K)));
      } else {
	mf.set_finite_element(cv,
                              getfem::classical_fem(pgt,short_type(param.K2)));
      }
    }
  }
}

static void init_mesh_im(getfem::mesh_im &mim, bool use_exact_im=true) {
  size_type cv;
  dal::bit_vector cvlst = mim.linked_mesh().convex_index();
  for (cv << cvlst; cv != size_type(-1); cv << cvlst) {      
    bgeot::pgeometric_trans pgt = mim.linked_mesh().trans_of_convex(cv);
    if ((cv+1) % 100) {
      mim.set_integration_method(cv, 
	(use_exact_im && (rand() % 10)==0) ? getfem::classical_exact_im(pgt) : 
		      getfem::classical_approx_im(pgt,dim_type(param.K*3)));
    } else {
      mim.set_integration_method(cv, 
      (use_exact_im && (rand() % 10)==0)  ? getfem::classical_exact_im(pgt) : 
		       getfem::classical_approx_im(pgt,dim_type(param.K2*3)));
    }
  }
}

static void comp_mat(const sparse_matrix_type& M1,
                     const sparse_matrix_type& M2) {
  scalar_type d = 0;
  scalar_type mx = 1e-200; /* avoid triggering an FPE for bound assembly */
                           /* when there is no boundary                  */
  sparse_vector_type r(gmm::mat_ncols(M1));
  for (size_type i = 0; i < gmm::mat_nrows(M1); ++i) {
    mx = std::max(mx,gmm::vect_norminf(gmm::mat_const_row(M1,i)));
    mx = std::max(mx,gmm::vect_norminf(gmm::mat_const_row(M2,i)));
    /*    int r = gmm::add(gmm::scaled(gmm::mat_const_row(M1,i), -1.0),
	  gmm::mat_row(M2,i));*/
    gmm::copy(gmm::mat_const_row(M2,i),r);
    gmm::add(gmm::scaled(gmm::mat_const_row(M1,i), -1.0),r);
    scalar_type d2 = gmm::vect_norminf(r);
    d = std::max(d,d2);
    if (mx > 1e-10 && d/mx > 1e-6) {
      sparse_vector_type r1(gmm::mat_ncols(M1));
      sparse_vector_type r2(gmm::mat_ncols(M2));
      gmm::copy(gmm::mat_const_row(M1,i),r1);
      gmm::copy(gmm::mat_const_row(M2,i),r2);    
      cout << "\nrow(" << i+1 << "),\nM1=" << r1 << "\nM2=" << r2 << endl;
      fail_cnt++;
      GMM_ASSERT1(false, "Failed ! ");
      break;
    }
  }
  assert(mx!=0.);
  cout << " ---> difference between assemblies: " << d / mx << "\n\n";
}

static void comp_vec(const base_vector& V1, const base_vector& V2) {
  scalar_type mx = std::max(gmm::vect_norminf(V1),gmm::vect_norminf(V2));
  base_vector dv = V2;
  gmm::add(gmm::scaled(V1, -1.0),dv);
  scalar_type d = gmm::vect_norminf(dv);
  if (mx != 0. && d/mx > 1e-6) {
    fail_cnt++;
    cout << " FAILED !";
  }
  assert(mx!=0.);
  cout << " ---> difference between assemblies: " << d / mx << "\n\n";
}



static double nrand() { return (::rand() % 10000) / 10000. + 0.01; }


static void run_tests(getfem::mesh_im &mim, 
               getfem::mesh_fem& mf, getfem::mesh_fem& mfq,
               getfem::mesh_fem& mfd, getfem::mesh_fem& mfdq,
               bool do_new, bool do_old, const std::vector<bool>& do_what,
               unsigned nloop, unsigned nloop_bound) {
  size_type Ndim = mf.linked_mesh().dim();
  base_vector V1q(Ndim*mf.nb_dof()), V2q(mfq.nb_dof());
  base_vector V1(mf.nb_dof()), V2(mf.nb_dof());
  sparse_matrix_type M1(mfq.nb_dof(),mfq.nb_dof());
  sparse_matrix_type M2(mfq.nb_dof(),mfq.nb_dof());

  chrono c;
    

  cout << "mf.nb_dof=" << mf.nb_dof() << " mfq=" << mfq.nb_dof() << endl;
  cout << "mfd.nb_dof=" << mfd.nb_dof() << " mfdq=" << mfdq.nb_dof() << endl;

  base_vector A(mfd.nb_dof()); std::generate(A.begin(), A.end(), nrand);
  base_vector A2(mfd.nb_dof()); std::generate(A2.begin(), A2.end(), nrand);
  base_vector Aq(mfdq.nb_dof()); std::generate(Aq.begin(), Aq.end(), nrand);


  /* --- SCALAR VOLUMIC SOURCE --- */
  if (do_what[DO_SCAL_VOLUMIC_SOURCE]) {
    if (do_old) {
      cout << "volumic source, Q=" << 1 << ", old way [" << nloop_bound
           << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop_bound; ++cnt) {
	gmm::clear(V1); c.tic(); //gmm::resize(M1, mfq.nb_dof(),mfq.nb_dof());
	getfem::old_asm_volumic_source_term(V1, mim, mf, mfd, A, 1u);
	c.toc(); cout << "#" << flushy;
      }
      cout << "done " << c << endl;
    }
    if (do_new) {
      cout << "volumic source, Q=" << 1 << ", new way [" << nloop_bound
           << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop_bound; ++cnt) {
	gmm::clear(V2); c.tic();
	getfem::asm_source_term(V2, mim, mf, mfd, A);
	c.toc(); cout << "#" << flushy;
      }
      cout << "done " << c << endl;
    }
    if (do_old && do_new) comp_vec(V1,V2);
  }
  //  cerr << "V1(old)=" << V1 << endl;
  //cerr << "V2(new)=" << V2 << endl;

  /* --- VECTOR VOLUMIC SOURCE --- */
  if (do_what[DO_VEC_VOLUMIC_SOURCE]) {
    if (do_old) {
      cout << "volumic source, Q=" << Ndim << ", old way ["
           << nloop_bound << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop_bound; ++cnt) {
	gmm::clear(V1q); c.tic(); //gmm::resize(M1, mfq.nb_dof(),mfq.nb_dof());
	getfem::old_asm_volumic_source_term(V1q, mim, mf, mfd,
                                            Aq, dim_type(Ndim));
	c.toc(); cout << "#" << flushy;
      }
      cout << "done " << c << endl;
    }
    if (do_new) {
      cout << "volumic source, Q=" << Ndim << ", new way ["
           << nloop_bound << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop_bound; ++cnt) {
	gmm::clear(V2q); c.tic();
	getfem::asm_source_term(V2q, mim, mfq, mfd, Aq);
	c.toc(); cout << "#" << flushy;
      }
      cout << "done " << c << endl;
    }  
    if (do_old && do_new) comp_vec(V1q,V2q);
  }

  /* --- SCALAR MASS MATRIX --- */
  if (do_what[DO_SCAL_MASS_MATRIX]) {
    if (do_old) {
      cout << "mass matrix, Q=" << 1 << ", old way ["
           << nloop << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop; ++cnt) {
	gmm::clear(M1); c.tic(); //gmm::resize(M1, mfq.nb_dof(),mfq.nb_dof());
	getfem::old_asm_mass_matrix(M1, mim, mf, mfd, 1);
      c.toc(); cout << "#" << flushy;
    }
    cout << "done " << c << endl;
  }
  if (do_new) {
    cout << "mass matrix, Q=" << 1 << ", new way ["
         << nloop << " times] .." << flushy;
    c.init();
    for (size_type cnt = 0; cnt < nloop; ++cnt) {
      gmm::clear(M2); c.tic(); 
      getfem::asm_mass_matrix(M2, mim, mf, mfd);
      c.toc(); cout << "#" << flushy;
    }
    cout << "done " << c << endl;
  }
  if (do_old && do_new) comp_mat(M1,M2);
  }

  /* --- VECTOR MASS MATRIX --- */
  if (do_what[DO_VEC_MASS_MATRIX]) {
  if (do_old) {
    cout << "mass matrix, Q=" << Ndim << ", old way [" << nloop
         << " times] .." << flushy;
    c.init();
    for (size_type cnt = 0; cnt < nloop; ++cnt) {
      gmm::clear(M1); c.tic(); //gmm::resize(M1, mfq.nb_dof(),mfq.nb_dof());
      getfem::old_asm_mass_matrix(M1, mim, mf, mfd, dim_type(Ndim));
      c.toc(); cout << "#" << flushy;
    }
    cout << "done " << c << endl;
  }
  if (do_new) {
    cout << "mass matrix, Q=" << Ndim << ", new way [" << nloop
         << " times] .." << flushy;
    c.init();
    for (size_type cnt = 0; cnt < nloop; ++cnt) {
      gmm::clear(M2); c.tic();
      getfem::asm_mass_matrix(M2, mim, mfq, mfdq);
      c.toc(); cout << "#" << flushy;
    }
    cout << "done " << c << endl;
  }
  if (do_old && do_new) comp_mat(M1,M2);
  }

  /* ---- LINEAR ELASTICITY ---- */
  if (do_what[DO_LIN_ELAST]) {
  if (do_old) {
    cout << "linear elasticity, Q=" << Ndim<<", old way ["
         << nloop << " times] .." << flushy;
    c.init();
    for (size_type cnt = 0; cnt < nloop; ++cnt) {
      gmm::clear(M1); c.tic();
      getfem::old_asm_stiffness_matrix_for_linear_elasticity(M1, mim, mf,
                                                             mfd, A, A2);
      c.toc(); cout << "#" << flushy;
    }
    cout << "done " << c << endl;
  }
  if (do_new) {
    cout << "linear elasticity, Q=" << Ndim << ", new way ["
         << nloop << " times] .." << flushy;
    c.init();
    for (size_type cnt = 0; cnt < nloop; ++cnt) {
      gmm::clear(M2); c.tic();
      getfem::asm_stiffness_matrix_for_linear_elasticity(M2, mim, mfq,
                                                         mfd, A, A2);
      c.toc(); cout << "#" << flushy;
    }
    cout << "done " << c << endl;

  }
  if (do_old && do_new) comp_mat(M1,M2);
  }
}


struct dummy_nonlin : public getfem::nonlinear_elem_term {
  unsigned i,j;
  bgeot::multi_index sizes_;
  dummy_nonlin(size_type N) : sizes_(2)
  { sizes_[0] = sizes_[1] = short_type(N); }
  const bgeot::multi_index &sizes(size_type) const { return sizes_; }
  virtual void compute(getfem::fem_interpolation_context& /*ctx*/,
		       bgeot::base_tensor &t) {
    t.adjust_sizes(sizes_); std::fill(t.begin(), t.end(), 0.);
    t[j*sizes_[0]+i] = 1.0;
  }
};

static void test_nonlin(const getfem::mesh_im &mim, const getfem::mesh_fem &mf)
{
  size_type N = mf.linked_mesh().dim();
  dummy_nonlin bidon(N);
  cerr << "testing assembly of nonlinear terms\n";
  for (bidon.i=0; bidon.i < N; ++bidon.i) {
    for (bidon.j=0; bidon.j < N; ++bidon.j) {
      std::vector<scalar_type> V1(mf.nb_dof()), V2(mf.nb_dof());
      char s[512]; sprintf(s,"t=comp(NonLin(#1).vGrad(#1));"
			   "V$1(#1) += t(i,j,:,i,j); "
			   "V$2(#1) += comp(vGrad(#1))(:,%d,%d)",
                           bidon.i+1, bidon.j+1);
      cout << s << "\n";
      getfem::generic_assembly assem(s);
      assem.push_mi(mim);
      assem.push_mf(mf);
      assem.push_nonlinear_term(&bidon);
      assem.push_vec(V1);
      assem.push_vec(V2);
      assem.assembly();
      gmm::add(gmm::scaled(V2,-1.),V1);
      scalar_type err = gmm::vect_norm2(V1);
      cout << "i=" << bidon.i << ", j=" << bidon.j << " |V1-V2| = "
           << err << "\n";
      assert(err < 1e-10);      
    }
  }
}

template<typename VECT1> class shape_der_nonlinear_term 
  : public getfem::nonlinear_elem_term {
  
  const getfem::mesh_fem &mf;
  const VECT1 &U;
  size_type N;
  base_vector coeff;
  base_matrix gradU, E, Sigma;
  bgeot::multi_index sizes_;
  scalar_type lambda, mu;
  
public:
  shape_der_nonlinear_term(const getfem::mesh_fem &mf_, const VECT1 &U_,
			  scalar_type lambda_, scalar_type mu_) 
    : mf(mf_), U(U_),
      N(mf_.get_qdim()),
      gradU(N, N), E(N, N), Sigma(N,N), sizes_(N,N),
      lambda(lambda_), mu(mu_) { }
  
  const bgeot::multi_index &sizes(size_type) const { return sizes_; }
  
  virtual void compute(getfem::fem_interpolation_context& ,
		       bgeot::base_tensor &t) {
    assert(t.size() == N*N);
    for (size_type i = 0; i < N; ++i) 
      for (size_type j = 0; j < N; ++j)
	t(i,j) = 0.0;
    
  }
};

void testbug() {
  std::vector<size_type> nsubdiv(3); nsubdiv[0] = nsubdiv[1] = nsubdiv[2] = 3;
  getfem::mesh m; 
  getfem::regular_unit_mesh(m, nsubdiv, bgeot::simplex_geotrans(3,1));
  getfem::mesh_fem mf1(m,3), mf2(m,3);
  mf1.set_classical_finite_element(m.convex_index(), 1);
  mf2.set_classical_finite_element(m.convex_index(), 1);
  getfem::mesh_im mim(m); mim.set_integration_method(m.convex_index(), 5);
  std::vector<scalar_type> U(mf1.nb_dof()), SD(mf2.nb_dof());
  gmm::fill_random(U);

  shape_der_nonlinear_term<std::vector<scalar_type> > nl(mf1, U, 0, 1);
  
  // trigers a real bug: printing an empty (because of the "vectorization"
  // of vBase) subtensor will crash
  getfem::generic_assembly assem2
    ("t=comp(vBase(#1).vBase(#2));"
     "print(t(:,:,2,3)); ");

  assem2.push_mi(mim);
  assem2.push_mf(mf1);
  assem2.push_nonlinear_term(&nl);
  assem2.push_mf(mf2);
  assem2.push_vec(SD);
  
  double t0 = gmm::uclock_sec();
  assem2.assembly();
  cerr << " done : " << gmm::uclock_sec() - t0 << "\n"; exit(1);
  exit(1);
}



static void test_new_assembly(void) {

    // std::string expr="([1,2;3,4]@[1,2;1,2])(:,2,1,1)(1)+ [1,2;3,4](1,:)(2)"; // should give 4
    // std::string expr="[1,2;3,4]@[1,2;1,2]*[2,3;2,1]/4 + [1,2;3,1]*[1;1](1)"; // should give [4, 8; 12, 13]
    // std::string expr="[1,2;3,a](2,:) + b(:)"; // should give [6, 9]
    // std::string expr="[1,1;1,2,,1,1;1,2;;1,1;1,2,,1,1;1,3](:,:,:,2)";
    // std::string expr="sin([0;pi;2*pi])";
    // std::string expr="Id(meshdim(u)+qdim(u))";
    // std::string expr="[sin(pi);-2] + Derivative_Norm(Grad_u) + Derivative_Norm(b) + Derivative_sin(pi)*[0;2]";
    // std::string expr="Trace([1,2;3,5;5,6]')";
    // std::string expr = "([1,2;3,4]@[1,2;1,2]).[1;2]";
    // std::string expr = "[u.u; u(1); (u./u)(1); a*Norm(u); c]";
    // std::string expr = "(3*(1*Grad_u)).Grad_Test_u*2 + 0*[1;2].Grad_Test_u + c*Grad_Test_u(1) + [u;1](1)*Test_u";
    // std::string expr = "-(4+(2*3)+2*(1+2))/-(-3+5)"; // should give 8
    // std::string expr="[1,2;3,4]@[1,2;1,2]*(Grad_u@Grad_u)/4 + [1,2;3,1]*[1;1](1)";
    std::string expr = "Test_u.Test2_u";

    getfem::ga_workspace workspace;

    base_vector a(1); a[0] = 3.0;
    workspace.add_fixed_size_constant("a", a);
    base_vector b(2); b[0] = 3.0; b[1] = 6.0;
    workspace.add_fixed_size_constant("b", b);
    // base_vector c(1); c[0] = 1.0;
    // workspace.add_fixed_size_variable("c", gmm::sub_interval(0, 1), c);
    
    getfem::mesh m;

    bgeot::pgeometric_trans pgt =
      bgeot::geometric_trans_descriptor("GT_PK(2,1)");
    size_type N = pgt->dim();
    std::vector<size_type> nsubdiv(N);
    std::fill(nsubdiv.begin(),nsubdiv.end(), 200);
    getfem::regular_unit_mesh(m, nsubdiv, pgt);

    getfem::mesh_fem mf_u(m);
    getfem::pfem pf_u = getfem::fem_descriptor("FEM_PK(2,2)");
    mf_u.set_finite_element(m.convex_index(), pf_u);
    mf_u.set_qdim(dim_type(N));

    getfem::mesh_fem mf_p(m);
    getfem::pfem pf_p = getfem::fem_descriptor("FEM_PK(2,2)");
    mf_p.set_finite_element(m.convex_index(), pf_p);
    // mf_p.set_qdim(dim_type(N));

    getfem::mesh_im mim(m);
    mim.set_integration_method(m.convex_index(), 4);
    
    getfem::mesh_im mim2(m);
    mim2.set_integration_method(m.convex_index(), 2);

    std::vector<scalar_type> U(mf_u.nb_dof());
    std::vector<scalar_type> P(mf_p.nb_dof(), 1.);
    size_type ndofu = mf_u.nb_dof(), ndofp = mf_p.nb_dof();
    cout << "ndofu = " << ndofu << " ndofp = " << ndofp << endl;
    
    gmm::sub_interval Iu(0, ndofu);
    gmm::sub_interval Ip(ndofu, ndofp);
    
    workspace.add_fem_variable("u", mf_u, Iu, U);
    workspace.add_fem_variable("p", mf_p, Ip, P);
    
    chrono ch;

    if (0) {
      cout << "\nTest for Mass matrix" << endl;
      workspace.clear_expressions();
      workspace.add_expression("Test_u.Test2_u", mim);
      
      ch.init(); ch.tic();
      workspace.assembly(2);
      ch.toc();
      cout << "Elapsed time for new assembly " << ch.elapsed() << endl;
      getfem::model_real_sparse_matrix K(ndofu, ndofu);
      ch.init(); ch.tic();
      getfem::asm_mass_matrix(K, mim, mf_u);
      ch.toc();
      cout << "Elapsed time for old assembly " << ch.elapsed() << endl;
      gmm::add(gmm::scaled(gmm::sub_matrix(workspace.assembled_matrix(),Iu,Iu),
                           scalar_type(-1)), K);
      scalar_type norm_error = gmm::mat_norminf(K);
      cout << "Error : " << norm_error << endl;
      GMM_ASSERT1(norm_error < 1E-11,
                  "Error in high or low level generic assembly");
    }
    

    if (1) {
      cout << "\nTest for Laplacian stiffness matrix" << endl;
      workspace.clear_expressions();
      workspace.add_expression("Grad_Test_p:Grad_Test2_p", mim2);
      // workspace.add_expression("(Grad_p:Grad_p)/2", mim2);
      // workspace.add_expression("sqr(Norm(Grad_p))/2", mim2);

      ch.init(); ch.tic();
      workspace.assembly(2);
      ch.toc();
      cout << "Elapsed time for new assembly " << ch.elapsed() << endl;
      getfem::model_real_sparse_matrix K(ndofp, ndofp);
      ch.init(); ch.tic();
      getfem::asm_stiffness_matrix_for_homogeneous_laplacian(K, mim2, mf_p);
      ch.toc();
      cout << "Elapsed time for old assembly " << ch.elapsed() << endl;
      gmm::add(gmm::scaled(gmm::sub_matrix(workspace.assembled_matrix(),Ip,Ip),
                           scalar_type(-1)), K);
      scalar_type norm_error = gmm::mat_norminf(K);
      cout << "Error : " << norm_error << endl;
      GMM_ASSERT1(norm_error < 1E-11,
                  "Error in high or low level generic assembly");
    }

    if (1) {
      cout << "\nTest for linear homogeneous elasticity stiffness matrix"
           << endl;
      workspace.clear_expressions();
      base_vector lambda(1); lambda[0] = 3.0;
      workspace.add_fixed_size_constant("lambda", lambda);
      base_vector mu(1); mu[0] = 2.0;
      workspace.add_fixed_size_constant("mu", mu);

//       workspace.add_expression
//         ("lambda*Trace(Grad_Test_u)*Trace(Grad_Test2_u) "
//          "+ mu*(Grad_Test_u'+Grad_Test_u):Grad_Test2_u", mim2);

      workspace.add_expression // Slighly better
        ("(lambda*Trace(Grad_Test_u)*Id(meshdim(u)) "
         "+ mu*(Grad_Test_u'+Grad_Test_u)):Grad_Test2_u", mim2);

      
      ch.init(); ch.tic();
      workspace.assembly(2);
      ch.toc();
      cout << "Elapsed time for new assembly " << ch.elapsed() << endl;
      getfem::model_real_sparse_matrix K(ndofu, ndofu);
      ch.init(); ch.tic();
      getfem::asm_stiffness_matrix_for_homogeneous_linear_elasticity
        (K, mim2, mf_u, lambda, mu);
      ch.toc();
      cout << "Elapsed time for old assembly " << ch.elapsed() << endl;
      gmm::add(gmm::scaled(gmm::sub_matrix(workspace.assembled_matrix(),Iu,Iu),
                           scalar_type(-1)), K);
      scalar_type norm_error = gmm::mat_norminf(K);
      cout << "Error : " << norm_error << endl;
      GMM_ASSERT1(norm_error < 1E-10,
                  "Error in high or low level generic assembly");
    }

    if (1) {
      cout << "\nTest for linear non homogeneous elasticity stiffness matrix"
           << endl;
      workspace.clear_expressions();
      base_vector lambda2(ndofp, 3.0);
      workspace.add_fem_constant("lambda2", mf_p, lambda2);
      base_vector mu2(ndofp, 2.0);
      workspace.add_fem_constant("mu2", mf_p, mu2);

      workspace.add_expression
        ("(lambda2*Trace(Grad_Test_u)*Id(meshdim(u)) "
         "+ mu2*(Grad_Test_u'+Grad_Test_u)):Grad_Test2_u", mim2);

      
      ch.init(); ch.tic();
      workspace.assembly(2);
      ch.toc();
      cout << "Elapsed time for new assembly " << ch.elapsed() << endl;
      getfem::model_real_sparse_matrix K(ndofu, ndofu);
      ch.init(); ch.tic();
      getfem::asm_stiffness_matrix_for_linear_elasticity
        (K, mim2, mf_u, mf_p, lambda2, mu2);
      ch.toc();
      cout << "Elapsed time for old assembly " << ch.elapsed() << endl;
      gmm::add(gmm::scaled(gmm::sub_matrix(workspace.assembled_matrix(),Iu,Iu),
                           scalar_type(-1)), K);
      scalar_type norm_error = gmm::mat_norminf(K);
      cout << "Error : " << norm_error << endl;
      GMM_ASSERT1(norm_error < 1E-10,
                  "Error in high or low level generic assembly");
    }



    // + faire un test sur un bord avec un fem réduit ...

}




int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.
  
  test_new_assembly();
  return 0;
  


  // testbug();
  
  param.init(argc,argv);
  std::vector<bool> tests(NB_TESTS, true);
  if (param.do_what >=0 && param.do_what < NB_TESTS) {
    std::fill(tests.begin(),tests.end(),false);
    tests[param.do_what]=true;
  }
  
  
  cerr << "\n\n-----------------------------PERFORMANCE TESTS------------"
       << "---------\n\n";   
  {
    getfem::mesh m; 
    gen_mesh(m);
    
    getfem::mesh_fem mf(m); 
    init_mesh_fem(mf,false);
    
    getfem::mesh_im mim(m);
    init_mesh_im(mim, false);
    
    getfem::mesh_fem mfq(m); 
    mfq.set_qdim(m.dim());
    init_mesh_fem(mfq,false);
    
    getfem::mesh_fem mfqne(m);
    init_mesh_fem(mfqne,false);
    
    test_nonlin(mim,mfq);
    
    getfem::mesh_fem mfd(m); 
    init_mesh_fem(mfd,true);
    
    getfem::mesh_fem mfdq(m); 
    mfdq.set_qdim(m.dim());
    init_mesh_fem(mfdq,true);
    
    run_tests(mim,mf,mfq,mfd,mfdq,param.do_new,param.do_old,tests,1,1);
  }
  
  cout << "failures: " << fail_cnt << endl;
  return fail_cnt; 
}
