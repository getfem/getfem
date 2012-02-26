/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2002-2012 Yves Renard, Julien Pommier.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/

/**
   @file mixed_elastostatic.cc
   @brief Linear Elastostatic problem. A dummy linear
   elastotatic problem is solved on a regular mesh, and is compared to
   the analytical solution.

   This program is used to check that getfem++ is working. This is also 
   a good example of use of Getfem++.

   @see laplacian.cc
   @see nonlinear_elastostatic.cc
*/

#include "getfem/getfem_config.h"
#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "gmm/gmm.h"
#include "getfem/getfem_interpolation.h"
#include "getfem/getfem_import.h"

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::short_type;
using bgeot::dim_type;
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. 
 * default types of getfem_model_solvers.h
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/**************************************************************************/
/*  Brique mixed elasticity.                                              */
/**************************************************************************/


namespace getfem {

  template<class MAT, class VECT>
  void asm_mixed_stiffness_matrix_for_linear_elasticity
  (const MAT &RM_, const mesh_im &mim, const mesh_fem &mf,
   const mesh_fem &mf_data, const VECT &LAMBDA_INV, const VECT &MU_INV,
   const mesh_region &rg = mesh_region::all_convexes()) {
    MAT &RM = const_cast<MAT &>(RM_);
    GMM_ASSERT1(mf_data.get_qdim() == 1,
		"invalid data mesh fem (Qdim=1 required)");
    
    GMM_ASSERT1(mf.get_qdim() == gmm::sqr(mf.linked_mesh().dim()),
		"wrong qdim for the mesh_fem");

    generic_assembly assem("lambda=data$1(#2); mu=data$2(#2);"
			   "t=comp(mBase(#1).mBase(#1).Base(#2));"
                           "M(#1,#1)+= sym(t(:,i,j,:,i,j,k).mu(k)*2"
			   "+ t(:,i,i,:,j,j,k).lambda(k))");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mf_data);
    assem.push_data(LAMBDA_INV);
    assem.push_data(MU_INV);
    assem.push_mat(RM);
    assem.assembly(rg);
  }


  template<class MAT>
  void asm_mixed_thetasitau
  (const MAT &RM_, const mesh_im &mim, const mesh_fem &mf_sigma,
   const mesh_fem &mf_theta,
   const mesh_region &rg = mesh_region::all_convexes()) {
    MAT &RM = const_cast<MAT &>(RM_);
    size_type N = mf_sigma.linked_mesh().dim();

    GMM_ASSERT1(mf_sigma.get_qdim() == N*N &&
		mf_theta.get_qdim() == (N*(N-1)/2),
		"wrong qdim for a mesh_fem");
    
    base_matrix A(N*(N-1)/2, N*N);
    size_type k = 0;
    for (size_type i = 0; i < N; ++i)
      for (size_type j = i + 1; j < N; ++j, ++k) {
	A(k, i*N+j) = 1.0;
	A(k, j*N+i) = -1.0;
      }

    generic_assembly assem("A=data(qdim(#1), qdim(#2));"
			   "t=comp(vBase(#1).vBase(#2));"
                           "M(#1,#2)+= t(:,i,:,j).A(i,j)");
    assem.push_mi(mim);
    assem.push_mf(mf_theta);
    assem.push_mf(mf_sigma);
    assem.push_mat(RM);
    assem.push_data(gmm::array1D_reference<scalar_type *>(&A(0,0), A.size()));
    assem.assembly(rg);
  }

  template<class MAT>
  void asm_mixed_udivtau
  (const MAT &RM_, const mesh_im &mim, const mesh_fem &mf_sigma,
   const mesh_fem &mf_u, const mesh_region &rg = mesh_region::all_convexes()) {
    MAT &RM = const_cast<MAT &>(RM_);
    size_type N = mf_sigma.linked_mesh().dim();

    GMM_ASSERT1(mf_sigma.get_qdim() == N*N && mf_u.get_qdim() == N,
		"wrong qdim for a mesh_fem");

    generic_assembly assem("t=comp(vBase(#1).mGrad(#2));"
                           "M(#1,#2)+= t(:,i,:,i,j,j)");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_sigma);
    assem.push_mat(RM);
    assem.assembly(rg);
  }




# define MDBRICK_MIXED_ELASTICITY 43224677

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_mixed_elasticity : public mdbrick_abstract<MODEL_STATE> {

    TYPEDEF_MODEL_STATE_TYPES;

    const mesh_im &mim;
    const mesh_fem &mf_sigma, &mf_u, &mf_theta;
    mdbrick_parameter<VECTOR> lambda_inv_, mu_inv_;
    T_MATRIX K;
    bool K_uptodate;
    size_type nbdof;

    virtual void proper_update(void) {
      K_uptodate = false;
      nbdof = mf_sigma.nb_dof() + mf_u.nb_dof() + mf_theta.nb_dof();
    }

  public :

    const T_MATRIX &get_K(void) {
      this->context_check(); 
      if (!K_uptodate || this->parameters_is_any_modified()) {
	GMM_ASSERT1(&lambda_inv_.mf() == &mu_inv_.mf(), 
		    "Lame coefficients should share the same mesh_fem");
	gmm::resize(K, nbdof, nbdof);
	gmm::clear(K);
	gmm::sub_interval I1(0, mf_sigma.nb_dof());
	gmm::sub_interval I2(mf_sigma.nb_dof(), mf_u.nb_dof());
	gmm::sub_interval I3(mf_sigma.nb_dof()+mf_u.nb_dof(),
			     mf_theta.nb_dof());
	VECTOR vlambda(lambda_inv_.get()), vmu(mu_inv_.get());

	asm_mixed_stiffness_matrix_for_linear_elasticity
	  (gmm::sub_matrix(K, I1), mim, mf_sigma, lambda_inv_.mf(), vlambda,
	   vmu, mf_u.linked_mesh().get_mpi_region());
	
	T_MATRIX B1(mf_u.nb_dof(), mf_sigma.nb_dof());
	asm_mixed_udivtau
	  (B1, mim, mf_sigma, mf_u, mf_u.linked_mesh().get_mpi_region());

	gmm::copy(gmm::transposed(B1), gmm::sub_matrix(K, I1, I2));	
	gmm::copy(B1, gmm::sub_matrix(K, I2, I1));
	
	T_MATRIX B2(mf_theta.nb_dof(), mf_sigma.nb_dof());
	asm_mixed_thetasitau
	  (B2, mim, mf_sigma, mf_theta, mf_u.linked_mesh().get_mpi_region());
	gmm::copy(gmm::transposed(B2), gmm::sub_matrix(K, I1, I3));	
	gmm::copy(B2, gmm::sub_matrix(K, I3, I1));

	K_uptodate = true;
	this->parameters_set_uptodate();
      }
      return K;
    }

    mdbrick_parameter<VECTOR> &lambda_inv(void) { return lambda_inv_; }
    const mdbrick_parameter<VECTOR> &lambda_inv(void) const
    { return lambda_inv_; }
    mdbrick_parameter<VECTOR> &mu_inv(void) { return mu_inv_; }
    const mdbrick_parameter<VECTOR> &mu_inv(void) const { return mu_inv_; }


    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0, nbdof);
      gmm::copy(get_K(), gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				size_type) {
      gmm::sub_interval SUBI(i0, nbdof);
      gmm::mult(get_K(), gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residual(), SUBI));
    }

    SUBVECTOR get_solution(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index(), nbdof);
      return gmm::sub_vector(MS.state(), SUBU);
    }
    SUBVECTOR get_U(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index()+mf_sigma.nb_dof(), mf_u.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    void init_(void) {
      size_type N = mf_sigma.linked_mesh().dim();
      GMM_ASSERT1(mf_sigma.get_qdim() == N*N,
		  "Qdim of mf_sigma should be " << N*N << ".");
      GMM_ASSERT1(mf_u.get_qdim() == N,
		  "Qdim of mf_u should be " << N << ".");
      GMM_ASSERT1(mf_theta.get_qdim() == (N*(N-1)/2),
		  "Qdim of mf_theta should be " << (N*(N-1)/2) << ".");
      this->add_proper_mesh_im(mim);
      this->add_proper_mesh_fem(mf_sigma, MDBRICK_MIXED_ELASTICITY);
      this->add_proper_mesh_fem(mf_u, MDBRICK_MIXED_ELASTICITY);
      this->add_proper_mesh_fem(mf_theta, MDBRICK_MIXED_ELASTICITY);
      this->force_update();
    }

    /* constructor for a homogeneous material (constant lambda and mu).
     * @param epsilon the thickness of the plate.
     */
    mdbrick_mixed_elasticity
    (const mesh_im &mim_, const mesh_fem &mf_sigma_, const mesh_fem &mf_u_,
     const mesh_fem &mf_theta_, value_type lambdai, value_type mui)
      : mim(mim_), mf_sigma(mf_sigma_), mf_u(mf_u_),
	mf_theta(mf_theta_), lambda_inv_("lambda", mf_u_.linked_mesh(), this),
	mu_inv_("mu", mf_u_.linked_mesh(), this) {
      size_type N = mf_u_.linked_mesh().dim();
      cout << "lambda_inv = " << -lambdai/(2*mui*(double(N)*lambdai+2*mui))
	   << endl;
      cout << "mu_inv = " << 1. / (4. * mui) << endl;
      lambda_inv_.set(-lambdai/(2*mui*(double(N)*lambdai+2*mui)));
      mu_inv_.set(1. / (4. * mui));
      init_();
    }
 
  };


   template<typename VECT1, typename VECT2>
   void asm_source_term_normal(VECT1 &B, const mesh_im &mim, const mesh_fem &mf,
			       const mesh_fem &mf_data, const VECT2 &F,
			       const mesh_region &rg) {
     GMM_ASSERT1(mf_data.get_qdim() == 1, "invalid data mesh_fem");

    const char *st;
    if (mf.get_qdim_n() == 1)
      st = "F=data(#2);"
	"V(#1)+=comp(vBase(#1).Base(#2).Normal())(:,j,k,j).F(k);";
    else
      st = "F=data(mdim(#1),#2);" // a corriger pour les tenseurs non carres
	"V(#1)+=comp(mBase(#1).Base(#2).Normal())(:,i,j,k,j).F(i,k);";

    asm_real_or_complex_1_param(B, mim, mf, mf_data, F, rg, st);
  }

  
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_source_term_normal : public mdbrick_abstract<MODEL_STATE>  {

    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_parameter<VECTOR> B_;
    VECTOR F_;
    bool F_uptodate;
    size_type boundary, num_fem, i1, nbd;

    void proper_update(void) {
      const mesh_fem &mf_u = this->get_mesh_fem(num_fem);
      i1 = this->mesh_fem_positions[num_fem];
      nbd = mf_u.nb_dof();

      gmm::resize(F_, mf_u.nb_dof());
      gmm::clear(F_);
      F_uptodate = false;
    }

  public :

    mdbrick_parameter<VECTOR> &source_term(void) {
      const mesh_fem &mf_u = this->get_mesh_fem(num_fem);
      B_.reshape(mf_u.get_qdim() / mf_u.linked_mesh().dim());
      return B_;
    }
    const mdbrick_parameter<VECTOR> &source_term(void) const { return B_; }

    // gives the right hand side of the linear system.
    const VECTOR &get_F(void) { 
      this->context_check();
      if (!F_uptodate || this->parameters_is_any_modified()) {
	const mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
	F_uptodate = true;
	GMM_TRACE2("Assembling a source term");
	asm_source_term_normal(F_, *(this->mesh_ims[0]), mf_u, B_.mf(), B_.get(),
			       mf_u.linked_mesh().get_mpi_sub_region(boundary));
	this->parameters_set_uptodate();
      }
      return F_;
    }

    virtual void do_compute_tangent_matrix(MODEL_STATE &, size_type,
					   size_type) { }
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::add(gmm::scaled(get_F(), value_type(-1)),
	       gmm::sub_vector(MS.residual(), gmm::sub_interval(i0+i1, nbd)));
    }

    /* Constructor not defining the rhs
	@param problem the sub-problem to which this brick applies.
	@param bound the mesh boundary number on which the source term is applied 
	(by default, it is a volumic source term as the whole mesh is taken).
	@param num_fem_ the mesh_fem number on which this brick is is applied.
    */
    mdbrick_source_term_normal(mdbrick_abstract<MODEL_STATE> &problem,
			size_type bound = size_type(-1), size_type num_fem_=0)
      : B_("source_term", this), boundary(bound),
	num_fem(num_fem_) {
      this->add_sub_brick(problem);
      if (bound != size_type(-1))
	this->add_proper_boundary_info(num_fem, bound, MDBRICK_NEUMANN);
      this->force_update();
      source_term();
    }
  };



}

/**************************************************************************/
/*  Exact solution.                                                       */
/**************************************************************************/

gmm::row_matrix<base_small_vector> sol_K;
static scalar_type sol_lambda, sol_mu, alph = 0.3;

base_small_vector sol_u(const base_node &x) {
  int N = x.size(); base_small_vector res(N);
  for (int i = 0; i < N; ++i)
    res[i] = alph * sin(gmm::vect_sp(sol_K.row(i), x));
  return res;
}

base_small_vector sol_f(const base_node &x) {
  int N = x.size();
  base_small_vector res(N);
  for (int i = 0; i < N; i++) {
    res[i] = alph * ( sol_mu * gmm::vect_sp(sol_K.row(i), sol_K.row(i)) )
      * sin(gmm::vect_sp(sol_K.row(i), x));
    for (int j = 0; j < N; j++)
      res[i] += alph * ( (sol_lambda + sol_mu) * sol_K(j,j) * sol_K(j,i))
	* sin(gmm::vect_sp(sol_K.row(j), x));
  }
  return res;
}

base_matrix sol_sigma(const base_node &x) {
  int N = x.size();
  base_matrix res(N,N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j <= i; j++) {
      res(j,i) = res(i,j) = alph * sol_mu *
	( sol_K(i,j) * cos(gmm::vect_sp(sol_K.row(i), x))
	  +  sol_K(j,i) * cos(gmm::vect_sp(sol_K.row(j), x))
	  );
      if (i == j)
	for (int k = 0; k < N; k++)
	  res(i,j) += alph * sol_lambda * sol_K(k,k)
	    * cos(gmm::vect_sp(sol_K.row(k), x));
    }
  return res;
}

/*
  structure for the elastostatic problem
*/
struct elastostatic_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1};
  getfem::mesh mesh;         /* the mesh */
  getfem::mesh_im mim;       /* the integration methods.                     */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the elastostatic solution */
  getfem::mesh_fem mf_sigma; /* mesh_fem for the stress tensor.              */
  getfem::mesh_fem mf_theta; /* mesh_fem for the multiplier for the symmetry *
			      * of the stress tensor.                        */
  getfem::mesh_fem mf_mult;  /* mesh_fem for the Dirichlet condition.        */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  scalar_type lambda, mu;    /* Lamé coefficients.                           */

  scalar_type residual;       /* max residual for iterative solvers          */
  getfem::constraints_type dirichlet_version;

  std::string datafilename;
  bgeot::md_param PARAM;

  bool solve(plain_vector &U);
  void init(void);
  void compute_error(plain_vector &U);
  elastostatic_problem(void) : mim(mesh),mf_u(mesh), mf_sigma(mesh),
			       mf_theta(mesh),mf_mult(mesh),
			       mf_rhs(mesh) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void elastostatic_problem::init(void) {
  std::string MESH_FILE = PARAM.string_value("MESH_FILE", "Mesh file");
  std::string FEM_TYPE_U  = PARAM.string_value("FEM_TYPE_U","FEM name");
  std::string FEM_TYPE_SIGMA = PARAM.string_value("FEM_TYPE_SIGMA","FEM name");
  std::string FEM_TYPE_THETA = PARAM.string_value("FEM_TYPE_THETA","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");

  cout << "MESH_FILE=" << MESH_FILE << "\n";
  cout << "FEM_TYPE_U="  << FEM_TYPE_U << "\n";
  cout << "FEM_TYPE_SIGMA="  << FEM_TYPE_SIGMA << "\n";
  cout << "FEM_TYPE_THETA="  << FEM_TYPE_THETA << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  size_type NX = PARAM.int_value("NX");
  size_type N = PARAM.int_value("N");
  std::stringstream filename; filename << MESH_FILE;
  if ((MESH_FILE.compare(0,11,"structured:") == 0) && NX > 0) {
    filename << ";NSUBDIV=[" << NX;
    for (size_type i = 1; i < N; ++i) filename << "," << NX;
    filename << "];";
  }
  getfem::import_mesh(filename.str(), mesh);
  
  GMM_ASSERT1(N == mesh.dim(), "The mesh has not the right dimension");

  dirichlet_version
    = getfem::constraints_type(PARAM.int_value("DIRICHLET_VERSION",
					       "Dirichlet version"));
  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  scalar_type FT = PARAM.real_value("FT", "parameter for exact solution");
  residual = PARAM.real_value("RESIDUAL");
  if (residual == 0.) residual = 1e-10;
  gmm::resize(sol_K, N, N);
  for (size_type i = 0; i < N; i++)
    for (size_type j = 0; j < N; j++)
      sol_K(i,j) = (i == j) ? 0.0 : FT;

  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  sol_lambda = lambda; sol_mu = mu;
  mf_u.set_qdim(dim_type(N));
  mf_mult.set_qdim(dim_type(N));
  mf_sigma.set_qdim_mn(dim_type(N), dim_type(N));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE_U);
  getfem::pfem pf_sigma = getfem::fem_descriptor(FEM_TYPE_SIGMA);
  getfem::pfem pf_theta = getfem::fem_descriptor(FEM_TYPE_THETA);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_sigma.set_finite_element(mesh.convex_index(), pf_sigma);
  mf_theta.set_finite_element(mesh.convex_index(), pf_theta);

  std::string dirichlet_fem_name
    = PARAM.string_value("DIRICHLET_FEM_TYPE", "Fem for dirichlet condition");
  cout << "DIRICHLET_FEM_TYPE="  << dirichlet_fem_name << "\n";
  mf_mult.set_finite_element(mesh.convex_index(), 
			     getfem::fem_descriptor(dirichlet_fem_name));

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }
  
  /* set boundary conditions
   * (Neuman on the upper face, Dirichlet elsewhere) */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
    base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);
    if (gmm::abs(un[N-1] - 1.0) < 0.5) { // new Neumann face
      mesh.region(NEUMANN_BOUNDARY_NUM).add(i.cv(), i.f());
    } else {
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
    }
  }
}

/* compute the error with respect to the exact solution */
void elastostatic_problem::compute_error(plain_vector &U) {
  size_type N = mesh.dim();
  std::vector<scalar_type> V(mf_rhs.nb_dof()*N), W(mf_rhs.nb_dof()*N);
  getfem::interpolation(mf_u, mf_rhs, U, V);
  getfem::interpolation_function(mf_rhs, W, sol_u);
  gmm::add(gmm::scaled(W, -1.0), V);

  cout.precision(16);
  mf_rhs.set_qdim(dim_type(N));
  scalar_type l2 = getfem::asm_L2_norm(mim, mf_rhs, V);
  scalar_type h1 = getfem::asm_H1_norm(mim, mf_rhs, V);

  cout << "L2 error = " << l2 << endl
       << "H1 error = " << h1 << endl
       << "Linfty error = " << gmm::vect_norminf(V) << endl;
  
  getfem::vtk_export exp(datafilename + "_err.vtk",
			 PARAM.int_value("VTK_EXPORT")==1);
  exp.exporting(mf_rhs); 
  exp.write_point_data(mf_rhs, V, "elastostatic_displacement");

  mf_rhs.set_qdim(1);
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/


bool elastostatic_problem::solve(plain_vector &U) {

  size_type N = mesh.dim();

  cout << "Number of dof for u: " << mf_u.nb_dof() << endl;
  cout << "Number of dof for sigma: " << mf_sigma.nb_dof() << endl;

  // Linearized elasticity brick.
  getfem::mdbrick_mixed_elasticity<>
    ELAS(mim, mf_sigma, mf_u, mf_theta, lambda, mu);
  
  // Volumic source term brick.
  getfem::mdbrick_source_term<> VOL_F(ELAS, size_type(-1), 1);
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  plain_vector F(nb_dof_rhs * N);
  getfem::interpolation_function(mf_rhs, F, sol_f);
  VOL_F.source_term().set(mf_rhs, gmm::scaled(F, -1.0));


  getfem::mdbrick_source_term_normal<>
    NEUMANN(VOL_F, NEUMANN_BOUNDARY_NUM);
  gmm::resize(F, nb_dof_rhs * N);
  getfem::interpolation_function(mf_rhs, F, sol_u, NEUMANN_BOUNDARY_NUM);
  NEUMANN.source_term().set(mf_rhs, F);


  // Dirichlet condition brick.
  getfem::mdbrick_normal_component_Dirichlet<>
    final_model(NEUMANN, DIRICHLET_BOUNDARY_NUM, mf_mult);
  final_model.set_constraints_type(dirichlet_version);
  final_model.set_coeff_dimension(2);
  
  gmm::resize(F, nb_dof_rhs * N * N);
  getfem::interpolation_function(mf_rhs, F, sol_sigma, DIRICHLET_BOUNDARY_NUM);
  final_model.rhs().set(mf_rhs, F);

  cout << "Total number of variables : " << final_model.nb_dof() << endl;
  getfem::standard_model_state MS(final_model);
  gmm::iteration iter(residual, 1, 40000);
  
  iter.init();
  getfem::standard_solve(MS, final_model, iter);
  gmm::resize(U, mf_u.nb_dof());
  gmm::copy(ELAS.get_U(MS), U);
 
  return (iter.converged());
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {

    elastostatic_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh");

    plain_vector U;

    if (!p.solve(U)) GMM_ASSERT1(false, "Solve has failed");

    p.compute_error(U);

    if (p.PARAM.int_value("VTK_EXPORT")) {
      cout << "export to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk",
			     p.PARAM.int_value("VTK_EXPORT")==1);
      exp.exporting(p.mf_u); 
      exp.write_point_data(p.mf_u, U, "elastostatic_displacement");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d " << p.datafilename << ".vtk -f ExtractVectorNorm -f "
	"WarpVector -m BandedSurfaceMap -m Outline\n";
    }

  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
