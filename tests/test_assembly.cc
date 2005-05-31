#include <getfem_assembling.h>
#include <getfem_export.h>
#include <getfem_regular_meshes.h>
#include <getfem_mat_elem.h>
#include <gmm.h>
#ifdef GETFEM_HAVE_SYS_TIMES
# include <sys/times.h>
#endif
#include <unistd.h>
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif

using bgeot::base_vector;
using bgeot::base_matrix;
using bgeot::base_small_vector;
using bgeot::base_node;
using bgeot::scalar_type;
using bgeot::size_type;
using bgeot::dim_type;

typedef gmm::wsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef std::vector<scalar_type> linalg_vector;


using std::flush;
#define flushy flush

int fail_cnt = 0;

void classical_mesh_fem(getfem::mesh_fem& mf, getfem::short_type K) {
  for (dal::bv_visitor cv(mf.linked_mesh().convex_index()); !cv.finished();
       ++cv) {
    bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
    mf.set_finite_element(cv, getfem::classical_fem(pgt,K));
  }
  //mf.set_classical_finite_element(K,2*K);
}
#define ASSEMBLY_CHECK

#ifdef ASSEMBLY_CHECK

typedef enum {DO_BOUNDARY_MASS,
      DO_SCAL_VOLUMIC_SOURCE,
      DO_VEC_VOLUMIC_SOURCE,
      DO_SCAL_MASS_MATRIX,
      DO_VEC_MASS_MATRIX,
      DO_SCAL_LAPLACIAN,
      DO_SCAL_L2_NORM,
      DO_VECT_H1_NORM,
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
  void tic() { t = ::clock()/float(CLOCKS_PER_SEC); }
  void toc() {
    float t2 = ::clock()/float(CLOCKS_PER_SEC);
    cpu_ += t2 - t; t = t2;
  }
  float cpu() const { return cpu_; }
  float elapsed() const { return cpu_; }
  float system() const { return 0.; }
};
#endif

std::ostream& operator<<(std::ostream& o, const chrono& c) {
  o << "[elapsed=" << int(c.elapsed()*1000) << "ms, cpu=" << int(c.cpu()*1000) << "ms, system=" << int(c.system()*1000) << "ms]";
  return o;
}

struct g_params {
  ftool::md_param PARAM;

  size_type NX,Ndim;
  int mesh_type;
  int K, K2, Kdata;
  bool do_new, do_old;
  int do_what;
  void init(int argc, char *argv[]);
};

g_params param;

void g_params::init(int argc, char *argv[]) {
  PARAM.add_int_param("NX", 50);
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
  mesh_type = PARAM.int_value("MESH_TYPE", "Mesh type ");
  K = PARAM.int_value("K", "Finite element degree");
  K2 = PARAM.int_value("K", "Finite element degree");
  Kdata = PARAM.int_value("KDATA", "Finite element degree for data meshfem");
  do_new = PARAM.int_value("BENCH_NEW", "bench new assembly routines");
  do_old = PARAM.int_value("BENCH_OLD", "bench old assembly routines");
  do_what = PARAM.int_value("BENCH_WHAT", "which test do you want to run?");
}

namespace getfem {
  template<class VECT1, class VECT2>
  void old_asm_Neumann_condition(VECT1 &B, const mesh_im &mim, const mesh_fem &mf,
				 size_type boundary, const mesh_fem &mfdata, const VECT2 &F, dim_type N)
  {
    size_type cv, nbd1, nbd2, f;
    dal::bit_vector nn = mf.convex_index(), nf;
    base_tensor t;
    pfem pf1, pf2, pf1prec = NULL, pf2prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");
  
    for (cv << nn; cv != ST_NIL; cv << nn) {
      nf = dal::bit_vector(mf.linked_mesh().region(boundary).faces_of_convex(cv));
      if (nf.card() > 0) {
	pf1 =     mf.fem_of_element(cv); nbd1 = pf1->nb_dof(cv);
	pf2 = mfdata.fem_of_element(cv); nbd2 = pf2->nb_dof(cv);
	pgt = mf.linked_mesh().trans_of_convex(cv);
	pim = mim.int_method_of_element(cv);
	if (pf1prec != pf1 || pf2prec != pf2 || pgtprec!=pgt || pimprec != pim) {
	  pme = mat_elem_product(mat_elem_base(pf1), mat_elem_base(pf2));
	  pmec = mat_elem(pme, pim, pgt);
	  pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
	}
	for (f << nf; f != ST_NIL; f << nf) {
	  pmec->gen_compute_on_face(t,mf.linked_mesh().points_of_convex(cv),f, cv);
	  base_tensor::iterator p = t.begin();
	  for (size_type i = 0; i < nbd2; i++)
	    {
	      size_type dof2 = mfdata.ind_dof_of_element(cv)[i];
	      for (size_type j = 0; j < nbd1; j++, ++p)
		{
		  size_type dof1 = mf.ind_dof_of_element(cv)[j];
		  for (size_type k = 0; k < N; k++) {
		    B[dof1*N + k] += F[dof2*N+k]*(*p);
		  }
		}
	    }
	  if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
	}
      }
    }
  }

  template<class VECT1, class VECT2>
  void old_asm_volumic_source_term(VECT1 &B, const mesh_im &mim, const mesh_fem &mf,
				   const mesh_fem &mfdata, const VECT2 &F, dim_type N)
  {
    size_type cv, nbd1, nbd2;
    dal::bit_vector nn = mf.convex_index();
    base_tensor t;
    pfem pf1, pf2, pf1prec = NULL, pf2prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    for (cv << nn; cv != ST_NIL; cv << nn)
      {
	pf1 =     mf.fem_of_element(cv); nbd1 = pf1->nb_dof(cv);
	pf2 = mfdata.fem_of_element(cv); nbd2 = pf2->nb_dof(cv);
	pgt = mf.linked_mesh().trans_of_convex(cv);
	pim = mim.int_method_of_element(cv);
	if (pf1prec != pf1 || pf2prec != pf2 || pgtprec != pgt || pimprec != pim)
	  {
	    pme = mat_elem_product(mat_elem_base(pf1), mat_elem_base(pf2));
	    pmec = mat_elem(pme, pim, pgt);
	    pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
	  }
	pmec->gen_compute(t, mf.linked_mesh().points_of_convex(cv), cv);
	base_tensor::iterator p = t.begin();
	for (size_type i = 0; i < nbd2; i++)
	  {
	    size_type dof2 = mfdata.ind_dof_of_element(cv)[i];
	    for (size_type j = 0; j < nbd1; j++, ++p)
	      {
		size_type dof1 = mf.ind_dof_of_element(cv)[j];
		for (size_type k = 0; k < N; k++) B[dof1*N + k] += F[dof2*N+k]*(*p);
	      }
	  }
	if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
      }
  }

  template<class MATRM, class MESH_FEM>
  void old_asm_mass_matrix_on_boundary(MATRM &M, const mesh_im &mim, const MESH_FEM &mf1,
				       const MESH_FEM &mf2, size_type boundary, dim_type N)
  {
    size_type cv, nbd1, nbd2, f;
    dal::bit_vector nn = mf1.convex_index();
    getfem::mesh_region::face_bitset nf;
    base_tensor t;
    pfem pf1, pf1prec = 0, pf2, pf2prec = 0;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;
    // M(0,0) = 1.0;  ??

    if (&(mf1.linked_mesh()) != &(mf2.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    for (cv << nn; cv != ST_NIL; cv << nn)
      {
	nf = mf1.linked_mesh().region(boundary).faces_of_convex(cv);
	if (nf.count() > 0) {
	  pf1 = mf1.fem_of_element(cv); nbd1 = pf1->nb_dof(cv);
	  pf2 = mf2.fem_of_element(cv); nbd2 = pf2->nb_dof(cv);
	  pgt = mf1.linked_mesh().trans_of_convex(cv);
	  pim = mim.int_method_of_element(cv);
	  if (pf1prec != pf1 || pf2prec != pf2 || pgtprec != pgt 
	      || pimprec != pim) {
	    pme = mat_elem_product(mat_elem_base(pf1), mat_elem_base(pf2));
	    pmec = mat_elem(pme, pim, pgt);
	    pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
	  }

	  for (f = 0; f < MAX_FACES_PER_CV; ++f) if (nf[f]) {

	    pmec->gen_compute_on_face(t,
				      mf1.linked_mesh().points_of_convex(cv),
				      f, cv);
	    
	    base_tensor::iterator p = t.begin();
	    for (size_type i = 0; i < nbd2; i++) {
	      size_type dof2 = mf2.ind_dof_of_element(cv)[i];
	      // cout << "cv = " << cv << " dof2 = " << dof2 << endl;
	      for (size_type j = 0; j < nbd1; j++, ++p) {
		size_type dof1 = mf1.ind_dof_of_element(cv)[j];
		// cout << "dof1 = " << dof1 << " dof2 = " << dof2 << endl;
		for (size_type k = 0; k < N; k++)
		  M(dof1*N + k, dof2*N + k) += (*p);
	      }
	    }
	    if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
	  }
	}
      }
  }

  template<class MATRM, class MESH_FEM>
  void old_asm_mass_matrix(MATRM &M, const mesh_im &mim, const MESH_FEM &mf1, const MESH_FEM &mf2, dim_type N)
  {
    size_type cv, nbd1, nbd2;
    dal::bit_vector nn = mf1.convex_index();
    base_tensor t;
    pfem pf1, pf1prec = 0, pf2, pf2prec = 0;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;
    // M(0,0) = 1.0;  ??

    if (&(mf1.linked_mesh()) != &(mf2.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    for (cv << nn; cv != ST_NIL; cv << nn)
      {
	pf1 = mf1.fem_of_element(cv); nbd1 = pf1->nb_dof(cv);
	pf2 = mf2.fem_of_element(cv); nbd2 = pf2->nb_dof(cv);
	pgt = mf1.linked_mesh().trans_of_convex(cv);
	pim = mim.int_method_of_element(cv);
	if (pf1prec != pf1 || pf2prec != pf2 || pgtprec != pgt || pimprec != pim)
	  {
	    pme = mat_elem_product(mat_elem_base(pf1), mat_elem_base(pf2));
	    pmec = mat_elem(pme, pim, pgt);
	    pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
	  }
	pmec->gen_compute(t, mf1.linked_mesh().points_of_convex(cv), cv);

	// cout << "t = " << t << endl;
      
	base_tensor::iterator p = t.begin();
	for (size_type i = 0; i < nbd2; i++)
	  {
	    size_type dof2 = mf2.ind_dof_of_element(cv)[i];
	    // cout << "cv = " << cv << " dof2 = " << dof2 << endl;
	    for (size_type j = 0; j < nbd1; j++, ++p)
	      {
		size_type dof1 = mf1.ind_dof_of_element(cv)[j];
		// cout << "dof1 = " << dof1 << " dof2 = " << dof2 << endl;
		for (size_type k = 0; k < N; k++)
		  M(dof1*N + k, dof2*N + k) += (*p);
	      }
	  }
	if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
      }
  }

  template<class MAT, class VECT>
  void old_asm_boundary_qu_term(MAT &M, 
				const mesh_im &mim, 
				const mesh_fem &mf_u, size_type boundary, 
				const mesh_fem &mf_d, const VECT &Q, dim_type N)
  {
    size_type cv;
    dal::bit_vector nn = mf_u.convex_index(), nf;
    base_tensor t;
    pfem pf_u, pf_d, pf_u_prec = NULL, pf_d_prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf_u.linked_mesh()) != &(mf_d.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    for (cv << nn; cv != ST_NIL; cv << nn)
      {
	nf = dal::bit_vector(mf_u.linked_mesh().region(boundary).faces_of_convex(cv));
	if (nf.card() > 0)
	  {
	    size_type f, nbdof_u, nbdof_d;

	    pf_u = mf_u.fem_of_element(cv); nbdof_u = pf_u->nb_dof(cv);
	    pf_d = mf_d.fem_of_element(cv); nbdof_d = pf_d->nb_dof(cv);
	    pgt = mf_u.linked_mesh().trans_of_convex(cv);
	    pim = mim.int_method_of_element(cv);
	    if (pf_u_prec != pf_u || pf_d_prec != pf_d || pgtprec!=pgt 
		|| pimprec != pim)
	      {
		pme = mat_elem_product(mat_elem_base(pf_d), 
				       mat_elem_product(mat_elem_base(pf_u),
							mat_elem_base(pf_u)));
		pmec = mat_elem(pme, pim, pgt);
		pf_u_prec = pf_u; pf_d_prec = pf_d; pgtprec = pgt; pimprec = pim;
	      }
	    for (f << nf; f != ST_NIL; f << nf)
	      {
		pmec->gen_compute_on_face(t,
				       mf_u.linked_mesh().points_of_convex(cv),
					  f, cv);
		base_tensor::iterator p = t.begin();
		scalar_type vmax = gmm::vect_norminf(base_vector(t));

		for (size_type j = 0; j < nbdof_u; j++) {
		  size_type dof_j = mf_u.ind_dof_of_element(cv)[j];
		  for (size_type i = 0; i < nbdof_u; i++) {
		    size_type dof_i = mf_u.ind_dof_of_element(cv)[i];
		    for (size_type id = 0; id < nbdof_d; id++) {

		      size_type dof_d = mf_d.ind_dof_of_element(cv)[id];

		      /* for every element of the matrix Q */
		      for (int ii=0; ii < N; ii++) {
			for (int jj=0; jj < N; jj++) {
			  /* get Q[ii][jj] for the degree of freedom 'dof_d' */
			  scalar_type data = Q[(jj*N+ii) + N*N*(dof_d)];

			  /* we filter out noise since this matrix can be used 
			     as a constraints matrix for dirichlet conditions,
			     noise may lead to 'fictive' dirichlet condition
			     (this is the case for ex. with laplace/PK(1,4)) 

			     NON NON ET NON !!
			     finaly we DON'T FILTER NOISE since it breaks 
			     the assembling of dirichlet conditions against
			     hierarchical FEMS ...
			  */
			  if (data != 0.) {// && vmax != .0 && (*p)/vmax > 1e-5) {
			    /*
			      cerr << "QU : adding " << data << "*" << (*p) << " at dof_i=" << 
			      dof_i << "*" << N << "+" << ii << ", dof_j=" << dof_i << "*" << 
			      N << "+" << ii << endl;
			    */
			    M(dof_i*N+ii, dof_j*N+jj) += data* (*p);
			  }
			}
		      }
		      p++;
		    }
		  }
		}
		if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
	      }
	  }
      }
  }

  template<class MAT, class VECT>
  void old_asm_stiffness_matrix_for_linear_elasticity(MAT &RM,
						      const mesh_im &mim, 
						     const mesh_fem &mf, 
						     const mesh_fem &mfdata, 
						     const VECT &LAMBDA, const VECT &MU)
  { // à verifier

    size_type cv, nbd2, N = mf.linked_mesh().dim();
    dal::bit_vector nn = mf.convex_index();
    base_tensor t;
    pfem pf1, pf2, pf1prec = NULL, pf2prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");
  
    for (cv << nn; cv != ST_NIL; cv << nn)
      {
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
      
	size_type nbd = mf.nb_dof_of_element(cv);

	for (size_type r = 0; r < nbd2; r++)
	  {
	    size_type dof3 = mfdata.ind_dof_of_element(cv)[r];
	    for (dim_type l = 0; l < N; l++)
	      for (size_type j = 0; j < nbd; j++)
		{
		  size_type dof2 = mf.ind_dof_of_element(cv)[j];
	    
		  for (dim_type k = 0; k < N; k++)
		    for (size_type i = 0; i < nbd; i++, ++p)
		      {
			size_type dof1 = mf.ind_dof_of_element(cv)[i];
		
			if (dof1*N + k >= dof2*N + l)
			  {
			    RM(dof1*N + k, dof2*N + l) += LAMBDA[dof3] * (*p);
			    RM(dof2*N + l, dof1*N + k) = RM(dof1*N + k, dof2*N + l);
			  }
		
			if (dof1*N + l >= dof2*N + k)
			  {
			    RM(dof1*N + l, dof2*N + k) += MU[dof3] * (*p);
			    RM(dof2*N + k, dof1*N + l) = RM(dof1*N + l, dof2*N + k);
			  }

			// cout << "matr elem : " << int(l) << " " << int(j) << " " << int(k) << " " << int(i) << " : " << *p << endl; getchar();
	      
			if (l == k && dof1 >= dof2)
			  for (size_type n = 0; n < N; ++n)
			    {
			      RM(dof1*N + n, dof2*N + n) += MU[dof3] * (*p);
			      RM(dof2*N + n, dof1*N + n) = RM(dof1*N + n, dof2*N + n);
			    }
		
		      }
		}
	  }
	if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
      }
  }

  template<class MAT, class VECT>
  void old_asm_mixed_pressure_term(MAT &B,
				   const mesh_im &mim, 
				   const mesh_fem &mf_u,
				   const mesh_fem &mf_p,
				   const mesh_fem &mf_d,
				   const VECT &DATA) {
    size_type cv;
    dal::bit_vector nn = mf_u.convex_index();

    base_tensor t;

    pmat_elem_computation pmec = 0;

    pfem pf_u, pf_p, pf_d;
    pfem pf_u_prec = NULL, pf_p_prec = NULL, pf_d_prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    size_type nbdof_u, nbdof_p, nbdof_d;
    size_type N = mf_u.linked_mesh().dim();

    if (&(mf_u.linked_mesh()) != &(mf_p.linked_mesh())
	|| &(mf_u.linked_mesh()) != &(mf_d.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    /* loop over all convexes */
    for (cv << nn; cv != ST_NIL; cv << nn) {

      pf_u = mf_u.fem_of_element(cv); nbdof_u = pf_u->nb_dof(cv);
      pf_p = mf_p.fem_of_element(cv); nbdof_p = pf_p->nb_dof(cv);
      pf_d = mf_d.fem_of_element(cv); nbdof_d = pf_d->nb_dof(cv);
      pgt = mf_u.linked_mesh().trans_of_convex(cv);
      pim = mim.int_method_of_element(cv);

      /* avoids recomputation of already known pmat_elem_computation */
      if (pf_u_prec != pf_u || pf_p_prec != pf_p || pf_d_prec != pf_d 
	  || pgtprec != pgt || pimprec != pim) {
	pmec = mat_elem(mat_elem_product(mat_elem_product(mat_elem_grad(pf_u),
							  mat_elem_base(pf_p)),
					 mat_elem_base(pf_d)), pim, pgt);
	pf_u_prec = pf_u;
	pf_p_prec = pf_p;
	pf_d_prec = pf_d; pgtprec = pgt; pimprec = pim;
      }
      pmec->gen_compute(t, mf_u.linked_mesh().points_of_convex(cv), cv);
      
      base_tensor::iterator p = t.begin();
      for (size_type i = 0; i < nbdof_d; i++) {
	size_type dof_d = mf_d.ind_dof_of_element(cv)[i];
	for (size_type j = 0; j < nbdof_p; j++) {
	  size_type dof_p = mf_p.ind_dof_of_element(cv)[j];
	  for (size_type l = 0; l < N; l++) {
	    // loop over derivation directions (d/dx, d/dy ..)
	    //	    for (size_type m = 0; m < N; m++) {
	    // loop over vector base function components (phi_x, phi_y ...)
	    for (size_type k = 0; k < nbdof_u; k++) {
	      //		if (m == l) {
	      /*
		ssert(finite(DATA[dof_d])); 
		  
		ssert(p < t.end());
		  
		ssert(finite(*p));
	      */
	      size_type dof_u = mf_u.ind_dof_of_element(cv)[k];
	      B(dof_u*N+l, dof_p) += DATA[dof_d]*(*p);
	      //		}
	      p++;
	    }
	    //	    }
	  } 
	}
      }
      if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
    }
  }

  template<class MAT, class VECT>
  void old_asm_stiffness_matrix_for_laplacian(MAT &RM, const mesh_im &mim, const mesh_fem &mf,
					     const mesh_fem &mfdata, const VECT &A)
  { // optimisable
    size_type cv, nbd1, nbd2, N = mf.linked_mesh().dim();
    dal::bit_vector nn = mf.convex_index();
    base_tensor t;
    pfem pf1, pf2, pf1prec = 0, pf2prec = 0;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = 0;
    pmat_elem_computation pmec = 0;

    if (&(mf.linked_mesh()) != &(mfdata.linked_mesh()))
      DAL_THROW(std::invalid_argument,
		"This assembling procedure only works on a single mesh");

    for (cv << nn; cv != ST_NIL; cv << nn)
      {
	pf1 =     mf.fem_of_element(cv); nbd1 = pf1->nb_dof(cv);
	pf2 = mfdata.fem_of_element(cv); nbd2 = pf2->nb_dof(cv);
	pgt = mf.linked_mesh().trans_of_convex(cv);
	pim = mim.int_method_of_element(cv);
	if (pf1prec != pf1 || pf2prec != pf2 || pgtprec != pgt || pimprec != pim)
	  {
	    pmat_elem_type pme = mat_elem_product(mat_elem_product(
								   mat_elem_grad(pf1), mat_elem_grad(pf1)), mat_elem_base(pf2));
	    pmec = mat_elem(pme, pim, pgt);
	    pf1prec = pf1; pf2prec = pf2; pgtprec = pgt; pimprec = pim;
	  }
	pmec->gen_compute(t, mf.linked_mesh().points_of_convex(cv), cv);
	// cout << "elem matrix " << t << endl;
	base_tensor::iterator p = t.begin();
	for (size_type r = 0; r < nbd2; r++) {
	  size_type dof3 = mfdata.ind_dof_of_element(cv)[r];
	  for (size_type l = 0; l < N; l++) {
	    for (size_type i = 0; i < nbd1; i++) {
	      size_type dof2 = mf.ind_dof_of_element(cv)[i];
	      p += l * nbd1;
	      for (size_type j = 0; j < nbd1; j++, ++p) {
		size_type dof1 = mf.ind_dof_of_element(cv)[j];
		if (dof1 >= dof2) { 
		  RM(dof1, dof2) += A[dof3]*(*p);
		  RM(dof2, dof1) = RM(dof1, dof2);
		}
	      }
	      p += (N-l-1) * nbd1;
	    }
	  }
	}
	if (p != t.end()) DAL_THROW(dal::internal_error, "internal error"); 
      }
  }

  template<class MESH_FEM, class VECT>
  scalar_type old_L2_norm(const mesh_im &mim, MESH_FEM &mf, const VECT &U, size_type N, const dal::bit_vector &cvlst)
  { /* optimisable */
    size_type cv;
    scalar_type no = 0.0;
    dal::bit_vector nn = cvlst;
    dal::dynamic_array<base_vector, 2> vval;
    base_tensor t;
    pfem pf1, pf1prec = NULL;
    pintegration_method pim, pimprec = 0;

    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;
    
    for (cv << nn; cv != ST_NIL; cv << nn)
      {
	pf1 =     mf.fem_of_element(cv);
	pgt = mf.linked_mesh().trans_of_convex(cv);
	pim = mim.int_method_of_element(cv);
	size_type nbd = mf.nb_dof_of_element(cv);
	if (pf1prec != pf1 || pgtprec != pgt || pimprec != pim)
	  {
	    pme = mat_elem_product(mat_elem_base(pf1), mat_elem_base(pf1));
	    pmec = mat_elem(pme, pim, pgt);
	    pf1prec = pf1; pgtprec = pgt; pimprec = pim;
	  }
	pmec->gen_compute(t, mf.linked_mesh().points_of_convex(cv), cv);
	base_tensor::iterator p = t.begin();

	for (size_type i = 0; i < nbd; i++)
	  { 
	    size_type dof1 = mf.ind_dof_of_element(cv)[i];
	    if (vval[i].size() != N) vval[i] = base_vector(N); 
	    for (size_type k = 0; k < N; k++) (vval[i])[k] = U[dof1*N+k];
	  }

	for (size_type i = 0; i < nbd; i++)
	  for (size_type j = 0; j < nbd; j++, ++p)
	    no += bgeot::vect_sp(vval[i], vval[j]) * (*p);
      
      }
    return sqrt(no);
  }

  template<class MESH_FEM, class VECT>
  scalar_type old_H1_semi_norm(const mesh_im &mim, MESH_FEM &mf, const VECT &U, size_type N, const dal::bit_vector& cvlst)
  { /* optimisable */
    size_type cv, NN = mf.linked_mesh().dim();
    scalar_type no = 0.0;
    dal::bit_vector nn = cvlst;
    dal::dynamic_array<base_vector, 2> vval;
    base_tensor t;
    pfem pf1, pf1prec = NULL;
    pintegration_method pim, pimprec = 0;
    bgeot::pgeometric_trans pgt, pgtprec = NULL;
    pmat_elem_type pme; pmat_elem_computation pmec = 0;

    for (cv << nn; cv != ST_NIL; cv << nn)
      {
	pf1 =     mf.fem_of_element(cv);
	pgt = mf.linked_mesh().trans_of_convex(cv);
	pim = mim.int_method_of_element(cv);
	size_type nbd = mf.nb_dof_of_element(cv);
	if (pf1prec != pf1 || pgtprec != pgt || pimprec != pim)
	  {
	    pme = mat_elem_product(mat_elem_grad(pf1), mat_elem_grad(pf1));
	    pmec = mat_elem(pme, pim, pgt);
	    pf1prec = pf1; pgtprec = pgt; pimprec = pim;
	  }
	pmec->gen_compute(t, mf.linked_mesh().points_of_convex(cv), cv);
	base_tensor::iterator p = t.begin();
	for (size_type i = 0; i < nbd; i++)
	  { 
	    size_type dof1 = mf.ind_dof_of_element(cv)[i];
	    if (vval[i].size() != N) vval[i] = base_vector(N); 
	    for (size_type k = 0; k < N; k++) (vval[i])[k] = U[dof1*N+k];
	  }
	for (size_type l = 0; l < NN; l++)
	  for (size_type i = 0; i < nbd; i++)
	    for (size_type k = 0; k < NN; k++)
	      for (size_type j = 0; j < nbd; j++, ++p)
		if (k == l)
		  no += (*p) * bgeot::vect_sp(vval[i], vval[j]);
      }
    return sqrt(no);
  }

  template<class MESH_FEM, class VECT>
  scalar_type old_H1_norm(const mesh_im &mim, MESH_FEM &mf, const VECT &U, size_type N, const dal::bit_vector& cvlst) {
    return sqrt( gmm::sqr(old_L2_norm(mim, mf, U, N, cvlst)) 
		 + gmm::sqr(old_H1_semi_norm(mim, mf, U, N, cvlst)));
  }


  /* old2 *******************************************************************
     inline reduction tests */

  template<class MAT, class VECT>
    void old2_asm_stiffness_matrix_for_linear_elasticity(const MAT &RM_,
							 const mesh_im &mim, 
                                                         const mesh_fem &mf,
                                                         const mesh_fem &mfdata,
                                                         const VECT &LAMBDA,const VECT &MU) {
    MAT &RM = const_cast<MAT &>(RM_);
    if (mfdata.get_qdim() != 1)
      DAL_THROW(invalid_argument, "invalid data mesh fem (Qdim=1 required)");
    
    if (mf.get_qdim() != mf.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");
    generic_assembly assem("lambda=data$1(#2); mu=data$2(#2);"
			   "M(#1,#1)+= sym(comp(vGrad(#1)(:,i,j).vGrad(#1)(:,i,j).Base(#2)(k).mu(k)) +"
			   "               comp(vGrad(#1)(:,j,i).vGrad(#1)(:,i,j).Base(#2)(k).mu(k)) +"
			   "               comp(vGrad(#1)(:,i,i).vGrad(#1)(:,j,j).Base(#2)(k).lambda(k)));");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_data(LAMBDA);
    assem.push_data(MU);
    assem.push_mat(RM);
    assem.assembly();
  }

  template<typename VEC>
  scalar_type old2_asm_L2_norm(const mesh_im &mim, const mesh_fem &mf, const VEC &U) {
    generic_assembly assem;    
    if (mf.get_qdim() == 1)
      assem.set("u=data(#1); v=u; V()+=comp(Base(#1)(i).Base(#1)(j).u(i).v(j))");
    else
      assem.set("u=data(#1); v=u;"
		"V()+=comp(vBase(#1)(i,k).vBase(#1)(j,k).u(i).v(j))");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(U);
    bgeot::vsvector<scalar_type> v(1);
    assem.push_vec(v);
    assem.assembly();
    return sqrt(v[0]);
  }

  template<typename VEC>
  scalar_type old2_asm_H1_norm(const mesh_im &mim, const mesh_fem &mf, const VEC &U) {
    generic_assembly assem;    
    if (mf.get_qdim() == 1)
      assem.set("u=data(#1); V()+=comp(Grad(#1)(i,d).Grad(#1)(j,d).u(i).u(j))");
    else
      assem.set("u=data(#1);"
		"V()+=comp(vGrad(#1)(i,k,d).vGrad(#1)(j,k,d).u(i).u(j))");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(U);
    bgeot::vsvector<scalar_type> v(1);
    assem.push_vec(v);
    assem.assembly();
    return sqrt(v[0] + gmm::sqr(old2_asm_L2_norm(mim,mf,U)));
  }

} /* namespace getfem */


void gen_mesh(getfem::getfem_mesh& mesh) {
  cout << "Mesh generation, N=" << param.NX << " Ndim=" << param.Ndim << endl;
  base_node org(param.Ndim); org.fill(0.0);
  std::vector<base_small_vector> vtab(param.Ndim);
  std::vector<size_type> ref(param.Ndim); std::fill(ref.begin(), ref.end(), param.NX);
  for (size_type i = 0; i < param.Ndim; i++) { 
    vtab[i] = base_small_vector(param.Ndim); vtab[i].fill(0.0);
    (vtab[i])[i] = 1. / scalar_type(param.NX);
  }
  switch (param.mesh_type) {
  case 0: getfem::parallelepiped_regular_simplex_mesh(mesh, param.Ndim, org,
						      vtab.begin(), ref.begin()); 
    cerr << mesh.convex_index().card() << " " << param.Ndim << "D simplexes generated\n";
    break;
  case 1 : getfem::parallelepiped_regular_mesh
		   (mesh, param.Ndim, org, vtab.begin(), ref.begin()); 
    cerr << mesh.convex_index().card() << " " << param.Ndim << "D parallelepipeds generated\n";
    break;
  case 2 : getfem::parallelepiped_regular_prism_mesh
		     (mesh, param.Ndim, org, vtab.begin(), ref.begin()); 
    cerr << mesh.convex_index().card() << " " << param.Ndim << "D prisms generated\n";
    break;
  default : DAL_THROW(dal::internal_error, "Unknown type of mesh");
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

void init_mesh_fem(getfem::mesh_fem &mf, bool datamf) {
  if (datamf)
    mf.set_classical_finite_element(param.Kdata);
  else {
    dal::bit_vector cvlst = mf.linked_mesh().convex_index();
    for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
      if ((cv+1) % 100) {
	mf.set_finite_element(cv, getfem::classical_fem(pgt,param.K));
      } else {
	mf.set_finite_element(cv, getfem::classical_fem(pgt,param.K2));
      }
    }
  }
}

void init_mesh_im(getfem::mesh_im &mim, bool use_exact_im=true) {
  size_type cv;
  dal::bit_vector cvlst = mim.linked_mesh().convex_index();
  for (cv << cvlst; cv != size_type(-1); cv << cvlst) {      
    bgeot::pgeometric_trans pgt = mim.linked_mesh().trans_of_convex(cv);
    if ((cv+1) % 100) {
      mim.set_integration_method(cv, 
				 (use_exact_im && (rand() % 10)==0) ? getfem::classical_exact_im(pgt) : 
				 getfem::classical_approx_im(pgt,param.K*3));
    } else {
      mim.set_integration_method(cv, 
				 (use_exact_im && (rand() % 10)==0)  ? getfem::classical_exact_im(pgt) : 
				 getfem::classical_approx_im(pgt,param.K2*3));
    }
  }
}

void comp_mat(const sparse_matrix_type& M1, const sparse_matrix_type& M2)
{
  scalar_type d = 0;
  scalar_type mx = 1e-200; /* avoid triggering an FPE for bound assembly when there is no boundary */
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
      cout << " FAILED !";
      break;
    }
  }
  assert(mx!=0.);
  cout << " ---> difference between assemblies: " << d / mx << "\n\n";
}

void comp_vec(const base_vector& V1, const base_vector& V2)
{
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

void comp_scal(scalar_type a, scalar_type b) {
  scalar_type d = gmm::abs(a-b)/std::max(gmm::abs(a),gmm::abs(b));
  if (d > 1e-10) {
    fail_cnt++;
    cout << " FAILED !";
  }
  cout << " ---> difference between assemblies: " << d << "\n\n";
}

#endif /* ASSEMBLY_CHECK */

base_node mknode(double a,double b,double c) {
  base_node n(3); n[0]=a; n[1]=b; n[2]=c; return n;
}
base_node mknode(double a,double b) {
  base_node n(2); n[0]=a; n[1]=b; return n;
}

void tensor_shape_check() {
  bgeot::tensor_ranges r1(5); r1[0] = 2; r1[1]=1; r1[2]=2; r1[3]=3; r1[4]=2;
  bgeot::tensor_ranges r2(4); r2[0] = 2; r2[1]=2; r2[2]=3; r2[3]=2;
  
  bgeot::tensor_shape m1(r1),m2(r2);

    cerr << "m1=\n" << m1 << endl;

  cerr << "m2=\n" << m2 << endl;

  cerr << "slice_shape(m1, 1, 0)=" << m1.slice_shape(bgeot::tensor_mask::Slice(1, 0)) << endl;
    

  bgeot::tensor_shape m3(m1.slice_shape(bgeot::tensor_mask::Slice(1, 0)));

  m3.remove_mask(1);m3.shift_dim_num_ge(1,-1);
  m3.set_ndim_noclean(4);m3.update_idx2mask();
  cerr << "m3=slice(m1)=\n" << m3 << endl;

  bgeot::tensor_shape m4(m2.diag_shape(bgeot::tensor_mask::Diagonal(0,1)));
  cerr << "m4=diag(m2)=\n" << m4 << endl;

  bgeot::tensor_shape m5(m4.diag_shape(bgeot::tensor_mask::Diagonal(1,3)));
  cerr << "m5=diag(m4)=\n" << m5 << endl;

  assert(m5.card()==6);

  bgeot::tensor_shape m6(m3); 
  m6.merge(m5, true);
  cerr << "m6=m5.and.m3=\n" << m6 << endl;
  assert(m6.card()==6);


  bgeot::tensor_shape m7(m3); 
  m7.merge(m5, false);
  //cerr << "m7=m5.or.m3=\n" << m7 << endl;
  assert(m7.card() == 24);
}

void tensor_ref_check() {
  scalar_type s1_[] = {1.0,2.0,3.0, 4.,5.,6., 7.,8.,9., 10.,11.,12.,13.,14.,15.,16.};
  scalar_type *s1 = s1_;
  bgeot::tensor_ranges r1(3); r1[0]=3; r1[1]=2; r1[2]=2;
  bgeot::tensor_ref tr1(r1,&s1);
  cerr << "tr1=" << tr1 << endl;
  
  bgeot::tensor_ref tr2(tr1, bgeot::tensor_mask::Slice(0,1));
  cerr << "tr2=tr1(1,:,:)=" << tr2 << endl;

  bgeot::tensor_ref tr20(tr1, bgeot::tensor_mask::Slice(1,1));
  cerr << "tr20=tr1(:,1,:)=" << tr20 << endl;

  bgeot::tensor_ref tr21(tr1, bgeot::tensor_mask::Slice(2,1));
  cerr << "tr21=tr1(:,:,1)=" << tr21 << endl;

  bgeot::tensor_ref tr22(tr21, bgeot::tensor_mask::Slice(1,1));
  cerr << "tr22=tr1(:,1,1)=" << tr22 << endl;

  bgeot::tensor_ref tr3(tr2, bgeot::tensor_mask::Diagonal(0,1));
  cerr << "tr3=tr2(i,i)=[2 0;0 11]=" << tr3 << endl;

  bgeot::tensor_ranges r2(4); r2[0] = 2; r2[1]=2; r2[2]=2; r2[3]=2;
  bgeot::tensor_ref tr4(r2, &s1);
  cerr << "tr4=" << tr4 << endl;

  bgeot::tensor_ref tr5_(bgeot::tensor_ref(tr4, bgeot::tensor_mask::Diagonal(2,3)));
  
  bgeot::tensor_ref tr5(bgeot::tensor_ref(tr5_, bgeot::tensor_mask::Diagonal(0,1)));

  cerr << "@@@@@@@@@@@@@-------------------------------------------\n" 
       << "tr5=tr4(i,i,j,j)=" << tr5 << endl;
  
  cerr << "-------------------------------------------reduction en cours...\n"; 
  bgeot::tensor_reduction red; red.insert(tr5, " i i");
  red.prepare(NULL); 
  red.do_reduction();
 bgeot::tensor_ref tr6; red.result(tr6);
  cerr << "@@@@@@@@@@@@@@@@@-------------------------------------------\n"; 
  cerr << "tr6=sum(i,tr5(:,i,:,i))=[1 0; 0 16]=" << tr6 << endl;
  assert(tr6.base()[0]==1 && 
	 tr6.base()[1]==0 && 
	 tr6.base()[2]==0 && 
	 tr6.base()[3]==16);
  cerr << "-------------------------------------------reduction 2 en cours...\n" ;
  red.clear();
  red.insert(tr5, "  kl");
  red.prepare(NULL); 
  red.do_reduction();
  red.result(tr6);
  cerr << "@@@@@@@@@@@@@@@@@-------------------------------------------\n"; 
  cerr << "sum(i,tr5(:,:,i,j))=[14 0; 0 20]=" << tr6 << endl;
  assert(tr6.base()[0]==14 && 
	 tr6.base()[1]==0 && 
	 tr6.base()[2]==0 && 
	 tr6.base()[3]==20);
  

  red.clear();
  red.insert(tr3, "ij");
  red.insert(tr5, "ijkl");

  red.prepare(NULL); 
  red.do_reduction();
  bgeot::tensor_ref tr7; red.result(tr7);
  cerr << "@@@@@@@@@@@@@@@@@-------------------------------------------\n"; 
  cerr << "tr7=sum(tr3(i,j)*tr5(i,j,l,k))=" << tr7 << endl;
  assert(tr7.base()[0] == 248.);
  

  red.clear();
  red.insert(tr4, "ijij");
  red.prepare(NULL); 
  red.do_reduction();
  bgeot::tensor_ref tr8; red.result(tr8);
  cerr << "@@@@@@@@@@@@@@@@@-------------------------------------------\n"; 
  cerr << "tr8=sum(i,j,tr4(i,j,i,j))=" << tr8 << endl;
  assert(tr8.base()[0] == 34.);
  
}

void tensor_ref_check2() {
  cout << "more checks with strange strides..\n";
  /* NOTE : ALTHOUGH THIS TEST PASSES, THERE IS SOMETHING BROKEN WHEN
     STRIDES.SIZE() != CARD() .. do not use that for emulation symmetricy
  */
  scalar_type s1_[] = {1.0,2.0,3.0,4.,5.,6., 7.,8.,9., 10.,11.,12.,13.,14.,15.,16.,17.,18.};
  std::vector<scalar_type> s1(s1_,  s1_+ sizeof s1_/sizeof(scalar_type));

  bgeot::tensor_ref tr1;
  tr1.set_ndim_noclean(3);
  bgeot::tensor_mask m0; m0.set_full(0,4);
  bgeot::tensor_strides strd0(4); std::fill(strd0.begin(),strd0.end(),0); strd0[1] = 1;
  tr1.push_mask(m0); tr1.strides().push_back(strd0);
  bgeot::tensor_mask m1; m1.set_full(1,3);
  bgeot::tensor_strides strd1(3); strd1[0] = 0; strd1[1] = 2; strd1[2] = 4;
  tr1.push_mask(m1); tr1.strides().push_back(strd1);
  bgeot::tensor_mask m2; m2.set_full(2,3);
  bgeot::tensor_strides strd2(3); strd2[0] = 6; strd2[1] = 0; strd2[2] = 6;
  tr1.push_mask(m2); tr1.strides().push_back(strd2);
  
  scalar_type *pbase = &s1[0];
  tr1.set_base(pbase);
  
  cerr << "tr1=" << tr1 << endl;

  scalar_type s2_[] = {1.0,2.0,3.0, 4.,5.,6., 7.,8.,9., 10.,11.,12.,13.,14.,15.,16.};
  scalar_type *s2 = s2_;
  bgeot::tensor_ranges r2(3); r2[0]=3; r2[1]=2; r2[2]=2;
  bgeot::tensor_ref tr2(r2,&s2);

  bgeot::tensor_reduction red; 
  red.insert(tr1, " ii");
  //red.insert(tr2, "i  ");
  red.prepare(NULL); 
  red.do_reduction();
  bgeot::tensor_ref tr9; red.result(tr9);
  
  cerr << "tr9 = " << tr9 << endl;
}  

bgeot::tensor_ref tr_from_matrix(const base_matrix &A, scalar_type *&p) {
  bgeot::tensor_ranges r(2); r[0] = A.nrows(); r[1] = A.ncols();
  bgeot::tensor_ref tr(r); tr.set_base(p);
  return tr;
}

void tensor_ref_check3(unsigned n1, unsigned n2, unsigned n3, unsigned n4, unsigned n5) {
  base_matrix A(n1,n2); base_matrix B(n2,n3); base_matrix C(n3,n4); base_matrix D(n4,n5);
  gmm::fill_random(A); gmm::fill_random(B); gmm::fill_random(C); gmm::fill_random(D);
  base_matrix AB(n1,n3); gmm::mult(A,B,AB);
  base_matrix ABC(n1,n4); gmm::mult(AB,C,ABC);
  base_matrix BC(n2,n4); gmm::mult(B,C,BC);
  base_matrix ABCD(n1,n5); gmm::mult(ABC,D,ABCD);
  scalar_type *pA = &A[0], *pB = &B[0], *pC = &C[0], *pD = &D[0];
  bgeot::tensor_ref trA = tr_from_matrix(A,pA);
  bgeot::tensor_ref trB = tr_from_matrix(B,pB);
  bgeot::tensor_ref trC = tr_from_matrix(C,pC);
  bgeot::tensor_ref trD = tr_from_matrix(D,pD);
  bgeot::tensor_reduction red;
  red.insert(trC, "jk");
  red.insert(trB, "ij");
  red.insert(trA, " i");
  red.insert(trD, "k ");
  cout << "A=" << A << "\nBC = " << BC << "\n";
  red.prepare(NULL); red.do_reduction();
  bgeot::tensor_ref trABCD; red.result(trABCD);
  cout << "Final: " << trABCD << "\n";
  cout << "ABCD = " << ABCD << "\n =?= " << "\n";
  for (unsigned j=0; j < n5; ++j)
    for (unsigned i=0; i < n1; ++i) {
      if (gmm::abs(trABCD.base()[i+j*n1] - ABCD(i,j)) > 1e-10) {
	cerr << "FAILED : " << i << ", " << j << ", " << trABCD.base()[i+j*n1] << "!=" << ABCD(i,j) << "\n";
	DAL_INTERNAL_ERROR("");
      }
    }
}

void
do_general_check() {
  tensor_shape_check();
  tensor_ref_check();  
  cerr << "Basic check OK..\n";
  tensor_ref_check3(4,2,3,2,5);
  cerr << "Advanced tensor check OK..\n";
}


double nrand() { return (::rand() % 10000) / 10000. + 0.01; }

#ifdef ASSEMBLY_CHECK

void
run_tests(getfem::mesh_im &mim, 
	  getfem::mesh_fem& mf, getfem::mesh_fem& mfq,
	  getfem::mesh_fem& mfd, getfem::mesh_fem& mfdq,
	  bool do_new, bool do_old, const std::vector<bool>& do_what, unsigned nloop, unsigned nloop_bound) {
  size_type Ndim = mf.linked_mesh().dim();
  bool do_new2 = do_new;
  base_vector V1q(Ndim*mf.nb_dof()), V2q(mfq.nb_dof());
  base_vector V1(mf.nb_dof()), V2(mf.nb_dof());
  sparse_matrix_type M1(mfq.nb_dof(),mfq.nb_dof()), M2(mfq.nb_dof(),mfq.nb_dof());

  chrono c;
    

  cout << "mf.nb_dof=" << mf.nb_dof() << " mfq=" << mfq.nb_dof() << endl;
  cout << "mfd.nb_dof=" << mfd.nb_dof() << " mfdq=" << mfdq.nb_dof() << endl;

  base_vector A(mfd.nb_dof()); std::generate(A.begin(), A.end(), nrand);
  base_vector A2(mfd.nb_dof()); std::generate(A2.begin(), A2.end(), nrand);
  base_vector Aq(mfdq.nb_dof()); std::generate(Aq.begin(), Aq.end(), nrand);
  scalar_type s1=0,s2=0;

  /* --- BOUNDARY MASS MATRIX --- */
  if (do_what[DO_BOUNDARY_MASS]) {
    if (do_old) {
      cout << "boundary mass matrix, old way [" << nloop_bound << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop_bound; ++cnt) {
	gmm::clear(M1); c.tic();
	getfem::old_asm_mass_matrix_on_boundary(M1, mim, mf, mfd, 1, Ndim);
	c.toc(); cout << "#" << flushy; 
      }
      cout << "done " << c << endl;
    }
    if (do_new) {
      cout << "boundary mass matrix, new way [" << nloop_bound << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop_bound; ++cnt) {
	gmm::clear(M2); c.tic();
	getfem::asm_mass_matrix(M2, mim, mfq, mfdq, size_type(1));
	c.toc(); cout << "#" << flushy;
      }
      cout << "done " << c << endl;    
    }
    if (do_old && do_new) comp_mat(M1,M2);
  }


  /* --- SCALAR VOLUMIC SOURCE --- */
  if (do_what[DO_SCAL_VOLUMIC_SOURCE]) {
    if (do_old) {
      cout << "volumic source, Q=" << 1 << ", old way [" << nloop_bound << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop_bound; ++cnt) {
	gmm::clear(V1); c.tic(); //gmm::resize(M1, mfq.nb_dof(),mfq.nb_dof());
	getfem::old_asm_volumic_source_term(V1, mim, mf, mfd, A, 1u);
	c.toc(); cout << "#" << flushy;
      }
      cout << "done " << c << endl;
    }
    if (do_new) {
      cout << "volumic source, Q=" << 1 << ", new way [" << nloop_bound << " times] .." << flushy;
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
      cout << "volumic source, Q=" << Ndim << ", old way [" << nloop_bound << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop_bound; ++cnt) {
	gmm::clear(V1q); c.tic(); //gmm::resize(M1, mfq.nb_dof(),mfq.nb_dof());
	getfem::old_asm_volumic_source_term(V1q, mim, mf, mfd, Aq, Ndim);
	c.toc(); cout << "#" << flushy;
      }
      cout << "done " << c << endl;
    }
    if (do_new) {
      cout << "volumic source, Q=" << Ndim << ", new way [" << nloop_bound << " times] .." << flushy;
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
      cout << "mass matrix, Q=" << 1 << ", old way [" << nloop << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop; ++cnt) {
	gmm::clear(M1); c.tic(); //gmm::resize(M1, mfq.nb_dof(),mfq.nb_dof());
	getfem::old_asm_mass_matrix(M1, mim, mf, mfd, 1);
      c.toc(); cout << "#" << flushy;
    }
    cout << "done " << c << endl;
  }
  if (do_new) {
    cout << "mass matrix, Q=" << 1 << ", new way [" << nloop << " times] .." << flushy;
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
    cout << "mass matrix, Q=" << Ndim << ", old way [" << nloop << " times] .." << flushy;
    c.init();
    for (size_type cnt = 0; cnt < nloop; ++cnt) {
      gmm::clear(M1); c.tic(); //gmm::resize(M1, mfq.nb_dof(),mfq.nb_dof());
      getfem::old_asm_mass_matrix(M1, mim, mf, mfd, Ndim);
      c.toc(); cout << "#" << flushy;
    }
    cout << "done " << c << endl;
  }
  if (do_new) {
    cout << "mass matrix, Q=" << Ndim << ", new way [" << nloop << " times] .." << flushy;
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



  /* ---- SCALAR LAPLACIAN ---- */
  if (do_what[DO_SCAL_LAPLACIAN]) {
  if (do_old) {
    cout << "laplacian, Q=1, old way [" << nloop << " times] .." << flushy;
    c.init();
    for (size_type cnt = 0; cnt < nloop; ++cnt) {
      gmm::clear(M1); c.tic();
      getfem::old_asm_stiffness_matrix_for_laplacian(M1, mim, mf, mfd, A);
      c.toc(); cout << "#" << flushy;
    }
    cout << "done " << c << endl;
  }
  if (do_new) {
    cout << "laplacian, Q=1, new way [" << nloop << " times] .." << flushy;
    c.init();
    for (size_type cnt = 0; cnt < nloop; ++cnt) {
      gmm::clear(M2); c.tic();
      getfem::asm_stiffness_matrix_for_laplacian(M2, mim, mf, mfd, A); 
      c.toc(); cout << "#" << flushy;
    }
    cout << "done " << c << endl;
  }
  if (do_old && do_new) comp_mat(M1,M2);
  }

  /* ---- L2 NORM, Q=1 ---- */
  if (do_what[DO_SCAL_L2_NORM]) {
    for (size_type i=0; i < mf.nb_dof(); ++i) { V1[i] = V2[i] = float(rand())/RAND_MAX; }
    if (do_old) {
      cout << "L2 norm, Q=" << 1 << ", old way [" << nloop << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop; ++cnt) {
	c.tic();
	s1 = getfem::old_L2_norm(mim, mf, V1, 1, mf.convex_index());
	c.toc(); cout << "#" << flushy;
      }
      cout << "done " << c << endl;
    }
    if (do_new) {
      cout << "L2 norm, Q=" << 1 << ", new way [" << nloop << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop; ++cnt) {
	c.tic();
	s2 = getfem::asm_L2_norm(mim, mf, V2);
	c.toc(); cout << "#" << flushy;
      }
      cout << "done " << c << endl;
    }
    if (do_old && do_new) comp_scal(s1,s2);
    if (do_new2) {
      cout << "L2 norm, Q=" << 1 << ", new way2 [" << nloop << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop; ++cnt) {
	c.tic();
	s2 = getfem::old2_asm_L2_norm(mim, mf, V2);
	c.toc(); cout << "#" << flushy;
      }
      cout << "done " << c << endl;
    }
    if (do_old && do_new2) comp_scal(s1,s2);
  }

  /* ---- VECT H1 NORM ---- */
  if (do_what[DO_VECT_H1_NORM]) {
    for (size_type i=0; i < mfq.nb_dof(); ++i) { V1q[i] = V2q[i] = float(rand())/RAND_MAX; }
    if (do_old) {
      cout << "H1 norm, Q=" << Ndim << ", old way [" << nloop << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop; ++cnt) {
	c.tic();
	s1 = getfem::old_H1_norm(mim, mf, V1q, Ndim, mf.convex_index());
	c.toc(); cout << "#" << flushy;
      }
      cout << "done " << c << endl;
    }
    if (do_new) {
      cout << "H1 norm, Q=" << Ndim << ", new way [" << nloop << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop; ++cnt) {
	c.tic();
	s2 = getfem::asm_H1_norm(mim, mfq, V2q);
	c.toc(); cout << "#" << flushy;
      }
      cout << "done " << c << endl;
    }
    if (do_old && do_new) comp_vec(V1q,V2q);

    if (do_new2) {
      cout << "H1 norm, Q=" << Ndim << ", new way2 [" << nloop << " times] .." << flushy;
      c.init();
      for (size_type cnt = 0; cnt < nloop; ++cnt) {
	c.tic();
	s2 = getfem::old2_asm_H1_norm(mim, mfq, V2q);
	c.toc(); cout << "#" << flushy;
      }
      cout << "done " << c << endl;
    }
    if (do_old && do_new2) comp_vec(V1q,V2q);
  }


  /* ---- LINEAR ELASTICITY ---- */
  if (do_what[DO_LIN_ELAST]) {
  if (do_old) {
    cout << "linear elasticity, Q=" << Ndim<<", old way [" << nloop << " times] .." << flushy;
    c.init();
    for (size_type cnt = 0; cnt < nloop; ++cnt) {
      gmm::clear(M1); c.tic();
      getfem::old_asm_stiffness_matrix_for_linear_elasticity(M1, mim, mf, mfd, A, A2);
      c.toc(); cout << "#" << flushy;
    }
    cout << "done " << c << endl;
  }
  if (do_new) {
    cout << "linear elasticity, Q=" << Ndim<<", new way [" << nloop << " times] .." << flushy;
    c.init();
    for (size_type cnt = 0; cnt < nloop; ++cnt) {
      gmm::clear(M2); c.tic();
      getfem::asm_stiffness_matrix_for_linear_elasticity(M2, mim, mfq, mfd, A, A2);
      c.toc(); cout << "#" << flushy;
    }
    cout << "done " << c << endl;

    cout << "linear elasticity, Q=" << Ndim<<", new way2 [" << nloop << " times] .." << flushy;
    c.init();
    for (size_type cnt = 0; cnt < nloop; ++cnt) {
      gmm::clear(M2); c.tic();
      getfem::old2_asm_stiffness_matrix_for_linear_elasticity(M2, mim, mfq, mfd, A, A2);
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
  dummy_nonlin(size_type N) : sizes_(2) { sizes_[0] = sizes_[1] = N; }
  const bgeot::multi_index &sizes() const { return sizes_; }
  virtual void compute(getfem::fem_interpolation_context& /*ctx*/,
		       bgeot::base_tensor &t) {
    t.adjust_sizes(sizes_); std::fill(t.begin(), t.end(), 0.);
    t[j*sizes_[0]+i] = 1.0;
  }
};

void test_nonlin(const getfem::mesh_im &mim, const getfem::mesh_fem &mf)
{
  size_type N = mf.linked_mesh().dim();
  dummy_nonlin bidon(N);
  cerr << "testing assembly of nonlinear terms\n";
  for (bidon.i=0; bidon.i < N; ++bidon.i) {
    for (bidon.j=0; bidon.j < N; ++bidon.j) {
      std::vector<scalar_type> V1(mf.nb_dof()), V2(mf.nb_dof());
      char s[512]; sprintf(s,"t=comp(NonLin(#1).vGrad(#1));"
			   "V$1(#1) += t(i,j,:,i,j); "
			   "V$2(#1) += comp(vGrad(#1))(:,%d,%d)", bidon.i+1, bidon.j+1);
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
      cout << "i=" << bidon.i << ", j=" << bidon.j << " |V1-V2| = " << err << "\n";
      assert(err < 1e-10);      
    }
  }
}

void inline_red_test(const getfem::mesh_im &mim, const getfem::mesh_fem &mf1, const getfem::mesh_fem &mf2) {
  std::vector<scalar_type> U(mf1.nb_dof()); gmm::fill_random(U);
  std::vector<scalar_type> V(mf2.nb_dof()); gmm::fill_random(V);
  getfem::generic_assembly assem;    
  cout << "INLINE RED\n";
  assem.set("u=data(#1);v=data$2(#2);"
            "V()+=u(i).v(j).comp(Grad(#1)(:,d).Grad(#2)(:,d))(i,j); print comp(Grad(#1)(:,d).Grad(#2)(:,d))");
  assem.push_mi(mim);
  assem.push_mf(mf1);
  assem.push_mf(mf2);
  assem.push_data(U);
  assem.push_data(V);
  std::vector<scalar_type> v1(1);
  assem.push_vec(v1);
  assem.assembly();
  
  cout << "OLD SCHOOL\n";
  getfem::generic_assembly assem2;
  assem2.set("u=data(#1);v=data$2(#2);"
            "V()+=u(i).v(j).comp(Grad(#1).Grad(#2))(i,d,j,d); print comp(Grad(#1).Grad(#2))(:,d,:,d)");
  assem2.push_mi(mim);
  assem2.push_mf(mf1);
  assem2.push_mf(mf2);
  assem2.push_data(U);
  assem2.push_data(V);
  std::vector<scalar_type> v2(1);
  assem2.push_vec(v2);
  assem2.assembly();

  cout << "v1 = " << v1 << ", v2 = " << v2 << endl;
  assert(gmm::abs(v1[0]-v2[0]) < 1e-14);
}

#endif /* ASSEMBLY_CHECK */

int main(int argc, char *argv[])
{
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
  try {

    // getfem::pfem pf = getfem::fem_descriptor("FEM_PK_PRISM_HIERARCHICAL(3,3)");
    

    do_general_check();

#ifdef ASSEMBLY_CHECK
    bool do_old=true, do_new=true;
    param.init(argc,argv);
    std::vector<bool> tests(NB_TESTS, true);
    if (param.do_what >=0 && param.do_what < NB_TESTS) {
      std::fill(tests.begin(),tests.end(),false);
      tests[param.do_what]=true;
    }

   cerr << "\n\n-----------------------------SIMPLE MESH TESTS---------------------\n\n";
   {
     getfem::getfem_mesh m(2);
     m.add_triangle_by_points(mknode(0.,0.),mknode(1.2,0.),mknode(0.1,1.5));     
     m.add_triangle_by_points(mknode(0.,0.),mknode(-1.2,0.),mknode(0.1,1.5));
     m.region(1).add(0, 0);
     m.region(1).add(0, 1);
     getfem::mesh_fem mf(m);
     classical_mesh_fem(mf, 2);     
     getfem::mesh_fem mfq(m); 
     mfq.set_qdim(m.dim());
     classical_mesh_fem(mfq, 2);
     getfem::mesh_fem mfd(m); 
     mfd.set_classical_finite_element(1);
     getfem::mesh_fem mfdq(m); 
     mfdq.set_classical_finite_element(1);
     mfdq.set_qdim(m.dim());     
     
     getfem::mesh_im mim(m);
     init_mesh_im(mim, true);
     inline_red_test(mim,mf,mfd);

     run_tests(mim,mf,mfq,mfd,mfdq,do_new,do_old,tests,1,1);
   }

   cerr << "\n\n-----------------------------PERFORMANCE TESTS---------------------\n\n";   
   {
     getfem::getfem_mesh m; 
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
     
     //mf.write_to_file("toto1.mf",true);
     //mfq.write_to_file("totoq.mf",true);

     run_tests(mim,mf,mfq,mfd,mfdq,param.do_new,param.do_old,tests,1,1);
   }
#endif /* ASSEMBLY_CHECK */
  }
  DAL_STANDARD_CATCH_ERROR;
  cout << "failures: " << fail_cnt << endl;
  return fail_cnt; 
}
