/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2005 Yves Renard, Julien Pommier.                    */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */

/**
 * Solve the Signorini problem with static Coulomb friction.
 * 
 * Research program.
 */


#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_import.h>
#include <getfem_regular_meshes.h>
#include <getfem_Coulomb_friction.h>
#include <getfem_derivatives.h>
#include <gmm.h>
#include <fstream>
/* try to enable the SIGFPE if something evaluates to a Not-a-number
 * of infinity during computations
 */
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. 
 * default types of getfem_modeling.h
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/*
 * structure for the friction problem
 */
struct friction_problem {

  enum { DIRICHLET_BOUNDARY, CONTACT_BOUNDARY, NEUMANN_BOUNDARY };
  getfem::getfem_mesh mesh;  /* the mesh */
  getfem::mesh_im  mim;      /* integration methods.                         */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the friction solution     */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_coef;  /* mesh_fem used to represent pde coefficients  */
  getfem::mesh_fem mf_vm;    /* mesh_fem used for the VonMises stress        */
  scalar_type lambda, mu;    /* Lamé coefficients.                           */
  scalar_type rho, PG;       /* density, and gravity                         */
  scalar_type friction_coef; /* friction coefficient.                        */

  scalar_type residue;       /* max residue for the iterative solvers        */
  
  size_type N, noisy, method, population, Dirichlet, Neumann;
  scalar_type r;
  scalar_type dtexport;
  scalar_type Dirichlet_ratio, Neumann_intensity;
  bool dxexport;

  std::string datafilename;
  ftool::md_param PARAM;

  void solve(void);
  void init(void);
  friction_problem(void)
    : mim(mesh), mf_u(mesh), mf_rhs(mesh), mf_coef(mesh), mf_vm(mesh) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void friction_problem::init(void) {
  const char *MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  const char *FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  const char *INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  std::string meshname
    (PARAM.string_value("MESHNAME", "Nom du fichier de maillage"));

  /* First step : build the mesh */
  bgeot::pgeometric_trans pgt = bgeot::geometric_trans_descriptor(MESH_TYPE);
  N = pgt->dim();
  if (meshname.compare(0,5, "splx:")==0) {
    std::vector<size_type> nsubdiv(N);
    std::fill(nsubdiv.begin(),nsubdiv.end(),
	      PARAM.int_value("NX", "Number of space steps "));
    getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			      PARAM.int_value("MESH_NOISED") != 0);
  
    bgeot::base_matrix M(N,N);
    for (size_type i=0; i < N; ++i) {
      static const char *t[] = {"LX","LY","LZ"};
      M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
    }
    mesh.transformation(M);
  }
  else getfem::import_mesh(meshname, mesh);
  mesh.optimize_structure();

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residue = PARAM.real_value("RESIDUE"); if (residue == 0.) residue = 1e-10;

  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  rho = PARAM.real_value("RHO", "Density");
  PG = PARAM.real_value("PG", "Gravity constant");
  friction_coef = PARAM.real_value("FRICTION_COEF", "Friction coefficient");

  Dirichlet = PARAM.int_value("DIRICHLET","Dirichlet condition or not");
  Dirichlet_ratio = PARAM.real_value("DIRICHLET_RATIO",
				     "parameter for Dirichlet condition");
  Neumann = PARAM.int_value("NEUMANN","Neumann condition or not");
  Neumann_intensity = PARAM.real_value("NEUMANN_INTENSITY",
				       "Intensity of the applied force");
  dxexport = (PARAM.int_value("DX_EXPORT", "Exporting on OpenDX format")
	      != 0);
  r = PARAM.real_value("R", "augmentation parameter");
  method = PARAM.int_value("METHOD", "solve method");
  population = PARAM.int_value("POPULATION", "genetic population");
  noisy = (PARAM.int_value("NOISY", "verbosity of iterative methods") != 0);
  mf_u.set_qdim(N);

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_vm.set_classical_discontinuous_finite_element(1);
  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  const char *data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name == 0) {
    if (!pf_u->is_lagrange()) {
      DAL_THROW(dal::failure_error, "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    }
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }
  
  mf_coef.set_finite_element(mesh.convex_index(),
			     getfem::classical_fem(pgt,0));

  /* set boundary conditions */
  base_node center(0.,0.,20.);
  std::cout << "Reperage des bord de contact et Dirichlet\n";  
  for (dal::bv_visitor cv(mesh.convex_index()); !cv.finished(); ++cv) {
    size_type nf = mesh.structure_of_convex(cv)->nb_faces();
    for (size_type f = 0; f < nf; f++) {
      if (bgeot::neighbour_of_convex(mesh, cv, f).empty()) {
	base_small_vector un = mesh.normal_of_face_of_convex(cv, f);
	un /= gmm::vect_norm2(un);	
	base_node pt = mesh.points_of_face_of_convex(cv,f)[0];
	if (un[N-1] < -0.000001 && (N != 3 || (bgeot::vect_dist2(pt, center)
			   > .99*sqrt(25. + 15*15) && pt[N-1] < 20.1)))
	  mesh.region(CONTACT_BOUNDARY).add(cv, f); 
	if ((un[N-1] > 0.1 && (Dirichlet == 1))
	    || (un[0] < -0.5 && (Dirichlet == 2)))
	  mesh.region(DIRICHLET_BOUNDARY).add(cv, f);
	if (un[N-1] > 0.1 && (Neumann == 1))
	  mesh.region(NEUMANN_BOUNDARY).add(cv, f);
      }
    }
  }
}

template <typename VEC1, typename VEC2>
void calcul_von_mises(const getfem::mesh_fem &mf_u, const VEC1 &U,
		      const getfem::mesh_fem &mf_vm, VEC2 &VM,
		      scalar_type mu=1) {
  assert(mf_vm.get_qdim() == 1); 
  unsigned N = mf_u.linked_mesh().dim(); assert(N == mf_u.get_qdim());
  std::vector<scalar_type> DU(mf_vm.nb_dof() * N * N);

  getfem::compute_gradient(mf_u, mf_vm, U, DU);
  
  gmm::resize(VM, mf_vm.nb_dof());
  scalar_type vm_min, vm_max;
  for (size_type i=0; i < mf_vm.nb_dof(); ++i) {
    VM[i] = 0;
    scalar_type sdiag = 0.;
    for (unsigned j=0; j < N; ++j) {
      sdiag += DU[i*N*N + j*N + j];
      for (unsigned k=0; k < N; ++k) {
	scalar_type e = .5*(DU[i*N*N + j*N + k] + DU[i*N*N + k*N + j]);
	VM[i] += e*e;	
      }      
    }
    VM[i] -= 1./N * sdiag * sdiag;
    vm_min = (i == 0 ? VM[0] : std::min(vm_min, VM[i]));
    vm_max = (i == 0 ? VM[0] : std::max(vm_max, VM[i]));
    assert(VM[i] > -1e-6);
  }
  cout << "Von Mises : min=" << 4*mu*mu*vm_min << ", max="
       << 4*mu*mu*vm_max << "\n";
  gmm::scale(VM, 4*mu*mu);
}

/**************************************************************************/
/*  Friction model brick for genetic algorithm.                           */
/**************************************************************************/

namespace getfem {

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_genetic_Coulomb_friction
    : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    size_type num_fem;

    T_MATRIX BN, BT;
    VECTOR gap, friction_coef;
    std::vector<int> situation;
    value_type r;
    size_type d, nbc;

    mesh_fem *mf_u;
    gmm::sub_interval SUBU, SUBN, SUBT;

    void proper_update(void) {
      mf_u = this->mesh_fems[num_fem];
      d = mf_u->linked_mesh().dim();
      r = value_type(1);
      gmm::resize(BN, nbc, mf_u->nb_dof());
      gmm::resize(BT, nbc*(d-1), mf_u->nb_dof());
      gmm::resize(gap, nbc); gmm::resize(friction_coef, nbc);
      this->proper_additional_dof = gmm::mat_nrows(BN) + gmm::mat_nrows(BT);
      this->proper_mixed_variables.clear();
      this->proper_mixed_variables.add(sub_problem.nb_dof(),
				       this->proper_additional_dof);
    }

    void precomp(MODEL_STATE &, size_type i0) {
      size_type i1 = this->mesh_fem_positions[num_fem];
      SUBU = gmm::sub_interval(i0 + i1, mf_u->nb_dof());
      SUBN = gmm::sub_interval(i0 + sub_problem.nb_dof(), gmm::mat_nrows(BN));
      SUBT = gmm::sub_interval(i0 + sub_problem.nb_dof() + gmm::mat_nrows(BN),
			       gmm::mat_nrows(BT));
    }

  public :

    inline size_type nb_contact_nodes(void) const
    { return gmm::mat_nrows(BN); }

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      precomp(MS, i0);
      gmm::copy(gmm::scaled(gmm::transposed(BN), value_type(-1)),
		gmm::sub_matrix(MS.tangent_matrix(), SUBU, SUBN));
      gmm::copy(gmm::scaled(gmm::transposed(BT), value_type(-1)), 
		gmm::sub_matrix(MS.tangent_matrix(), SUBU, SUBT));
      gmm::copy(gmm::scaled(BN, value_type(-1)),
		gmm::sub_matrix(MS.tangent_matrix(), SUBN, SUBU));
      gmm::copy(gmm::scaled(BT, value_type(-1)), 
		gmm::sub_matrix(MS.tangent_matrix(), SUBT, SUBU));
      gmm::clear(gmm::sub_matrix(MS.tangent_matrix(),
				 gmm::sub_interval(SUBN.first(),
				 gmm::mat_nrows(BN) + gmm::mat_nrows(BT))));
      for (size_type i=0; i < nb_contact_nodes(); ++i) {
	switch (situation[i]) {
	case 0 : // no contact
	  gmm::clear(gmm::sub_matrix(MS.tangent_matrix(),
		     gmm::sub_interval(SUBN.first()+i,1), SUBU));
	  gmm::clear(gmm::sub_matrix(MS.tangent_matrix(),
		     gmm::sub_interval(SUBT.first()+i,1), SUBU));
	  MS.tangent_matrix()(SUBN.first()+i, SUBN.first()+i)=-value_type(1)/r;
	  MS.tangent_matrix()(SUBT.first()+i, SUBT.first()+i)=-value_type(1)/r;
	  break;
	case 1 : // stick
	  break;
	case 2 : // slip direction 1
	  gmm::clear(gmm::sub_matrix(MS.tangent_matrix(),
		     gmm::sub_interval(SUBT.first()+i,1), SUBU));
	  MS.tangent_matrix()(SUBT.first()+i, SUBN.first()+i)
	    = -friction_coef[i]/r;
	  MS.tangent_matrix()(SUBT.first()+i, SUBT.first()+i)
	    = -value_type(1)/r;
	  break;
	case 3 : // slip direction 2
	  gmm::clear(gmm::sub_matrix(MS.tangent_matrix(),
		     gmm::sub_interval(SUBT.first()+i,1), SUBU));
	  MS.tangent_matrix()(SUBT.first()+i, SUBN.first()+i)
	    = -friction_coef[i]/r;
	  MS.tangent_matrix()(SUBT.first()+i, SUBT.first()+i)
	    = value_type(1)/r;
	  break;
	}
      }
    }
    
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0, size_type) {
      precomp(MS, i0);
      value_type c1(1);
      gmm::clear(gmm::sub_vector(MS.residu(), SUBN));
      gmm::clear(gmm::sub_vector(MS.residu(), SUBT));
      gmm::mult(BN, gmm::scaled(gmm::sub_vector(MS.state(), SUBU), -c1), gap,
		gmm::sub_vector(MS.residu(), SUBN));
      gmm::mult(BT, gmm::scaled(gmm::sub_vector(MS.state(), SUBU), -c1),
		gmm::sub_vector(MS.residu(), SUBT));
      gmm::mult_add(gmm::transposed(BN),
		    gmm::scaled(gmm::sub_vector(MS.state(), SUBN),-c1),
		    gmm::sub_vector(MS.residu(), SUBU));
      gmm::mult_add(gmm::transposed(BT),
		    gmm::scaled(gmm::sub_vector(MS.state(), SUBT),-c1),
		    gmm::sub_vector(MS.residu(), SUBU));
      for (size_type i=0; i < nb_contact_nodes(); ++i) {
	switch(situation[i]) {
	case 0 :
	  MS.residu()[SUBN.first()+i] = -(MS.state()[SUBN.first()+i])/r;
	  MS.residu()[SUBT.first()+i] = -(MS.state()[SUBT.first()+i])/r;
	  break;
	case 1 :
	  break;
	case 2 :
	  MS.residu()[SUBT.first()+i] = - (MS.state()[SUBT.first()+i])/r
	    - friction_coef[i]*(MS.state()[SUBN.first()+i])/r;
	  break;
	case 3 :
	  MS.residu()[SUBT.first()+i] = + (MS.state()[SUBT.first()+i])/r
	    - friction_coef[i]*(MS.state()[SUBN.first()+i])/r;
	  break;
	default : assert(false);
	}
      }
    }

    void init(void) {
      this->add_sub_brick(sub_problem);
      this->proper_is_linear_ = true; 
      this->proper_is_coercive_ = false;
      this->proper_is_symmetric_ =  false;
      this->update_from_context();
    }

    void change_situation(std::vector<int> &sit) 
    { gmm::copy(sit, situation); }
    void set_r(value_type r_) { r = r_; }
    value_type get_r(void) const { return r; }

    VECTOR &get_gap(void) { return gap; }
    const VECTOR &get_gap(void) const { return gap; }

    SUBVECTOR get_LN(MODEL_STATE &MS) {
      SUBN = gmm::sub_interval(this->first_index() + sub_problem.nb_dof(),
			       gmm::mat_nrows(BN));
      return gmm::sub_vector(MS.state(), SUBN);
    }

    SUBVECTOR get_LT(MODEL_STATE &MS) {
      SUBT = gmm::sub_interval(this->first_index() + sub_problem.nb_dof()
			       + gmm::mat_nrows(BN),  gmm::mat_nrows(BT));
      return gmm::sub_vector(MS.state(), SUBT);
    }

    template <class MAT, class VEC> mdbrick_genetic_Coulomb_friction
    (mdbrick_abstract<MODEL_STATE> &problem, const MAT &BN_, const VEC &gap_,
     scalar_type FC_, const MAT &BT_, size_type num_fem_=0)
      : sub_problem(problem), num_fem(num_fem_) {
      nbc = gmm::mat_nrows(BN_);
      situation.resize(nbc);
      init();
      gmm::copy(BN_, BN); gmm::copy(BT_, BT); gmm::copy(gap_, gap);
      std::fill(friction_coef.begin(), friction_coef.end(), FC_);
    }

  };
}



/**************************************************************************/
/*  Function determining the situation of a solution.                     */
/**************************************************************************/

template <class Matrix, class Vector1, class Vector2>
void situation_of(const Vector1 &U, const Vector2 &gap,
		  const Matrix &BN, const Matrix &BT,
		  std::vector<int> &situation) {
  size_type nbc = gmm::mat_nrows(BN);
  plain_vector RLN(nbc), RLT(nbc);
  gmm::mult(BN, U, gmm::scaled(gap, -1.0), RLN);
  gmm::mult(BT, U, RLT);
  for (size_type j = 0; j < nbc; ++j) {
    if (RLN[j] < -1E-10) situation[j] = 0;
    else if (gmm::abs(RLT[j]) < 1E-10) situation[j] = 1;
    else if (RLT[j] < 0) situation[j] = 2;
    else situation[j] = 3;
  }
}


/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

struct fellow {
  scalar_type residue;
  std::vector<int> *situation_;
  std::vector<int> &situation(void) { return *situation_; }
  const std::vector<int> &situation(void) const { return *situation_; }
  fellow(void) : situation_(0) {}
};

bool  operator < (const fellow &a, const fellow &b)
{ return (a.residue < b.residue); }


void friction_problem::solve(void) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  N = mesh.dim();
  cout << "Number of dof for u: " << mf_u.nb_dof() << endl;

  // Linearized elasticity brick.
  getfem::mdbrick_isotropic_linearized_elasticity<>
    ELAS(mim, mf_u, mf_coef, lambda, mu);

  // Defining the volumic source term brick.
  plain_vector F(nb_dof_rhs * N);
  plain_vector f(N);
  if (Dirichlet == 2) f[0] = -rho*PG; else f[N-1] = -rho*PG;
  for (size_type i = 0; i < nb_dof_rhs; ++i)
      gmm::copy(f,gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  getfem::mdbrick_source_term<> VOL_F(ELAS, mf_rhs, F);
  
  // Defining the applied force source term brick
  gmm::clear(f);
  if (Neumann == 1) f[N-1] = Neumann_intensity;
  for (size_type i = 0; i < nb_dof_rhs; ++i)
      gmm::copy(f,gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  getfem::mdbrick_source_term<> NEUMANN_F(VOL_F, mf_rhs, F, NEUMANN_BOUNDARY);
 

  // Dirichlet condition brick.
  gmm::clear(F);
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    F[(i+1)*N-1] = Dirichlet_ratio * mf_rhs.point_of_dof(i)[N-1];
  getfem::mdbrick_Dirichlet<> DIRICHLET(NEUMANN_F, mf_rhs, F,
					DIRICHLET_BOUNDARY);
  
  // contact condition for Lagrange elements
  dal::bit_vector cn = mf_u.dof_on_set(CONTACT_BOUNDARY);
  dal::bit_vector dn = mf_u.dof_on_set(DIRICHLET_BOUNDARY);
  cn.setminus(dn);
  size_type nbc = cn.card()/N;
  sparse_matrix BN(nbc, mf_u.nb_dof());
  sparse_matrix BT((N-1)*nbc, mf_u.nb_dof());
  plain_vector gap(nbc);
  size_type jj = 0;
  for (dal::bv_visitor i(cn); !i.finished(); ++i)
    if (i % N == 0) {
      BN(jj, i+N-1) = -1.;
      gap[jj] = mf_u.point_of_dof(i)[N-1];
      for (size_type k = 0; k < N-1; ++k) BT((N-1)*jj+k, i+k) = 1.;
      ++jj;
    }

  getfem::mdbrick_Coulomb_friction<>
    FRICTION(DIRICHLET, BN, gap, friction_coef, BT);
  
  cout << "Total number of variables: " << FRICTION.nb_dof() << endl;
  getfem::standard_model_state MS(FRICTION);
 

  std::vector<scalar_type> U(mf_u.nb_dof());

  switch (method) {
  case 0 : {
      gmm::iteration iter(residue, noisy, 40000);
      getfem::standard_solve(MS, FRICTION, iter);
      gmm::copy(ELAS.get_solution(MS), U);
    }
    break;
  case 1 : {
      assert(N==2);
      std::vector<fellow> people(population);
      
      // Specific brick for the linear system
      getfem::mdbrick_genetic_Coulomb_friction<>
	FRICTION2(DIRICHLET, BN, gap, friction_coef, BT);
      getfem::standard_model_state MS2(FRICTION2);
      FRICTION.compute_tangent_matrix(MS);
      size_type first_computed = 0, nbiter = 0;

      for (;;) {

	if (nbiter == 0) {
	  // initial random population
	  for (size_type i = 0; i < population; ++i) {
	    people[i].situation_ = new std::vector<int>(nbc);
	    for (size_type j = 0; j < nbc; ++j)
	      people[i].situation()[j] = (rand() & 3);
	    // cout << "fellow " << i <<" : " << people[i].situation() << endl;
	  }
	  first_computed = 0;
	}

	// compute residue for the new fellows
	for (size_type i = first_computed; i < population; ++i) {
	  cout << "computing solution " << i << " on " << population << endl;
	  FRICTION2.change_situation(people[i].situation());
	  gmm::iteration iter(residue, noisy, 40000);
	  getfem::standard_solve(MS2, FRICTION2, iter);
	  gmm::copy(MS2.state(), MS.state());
	  FRICTION.compute_residu(MS);
	  MS.compute_reduced_system();
	  people[i].residue = MS.reduced_residu_norm();
	}
	
	// elimination of the 10% worst
	first_computed = (population * 90) / 100;
	std::sort(people.begin(), people.end());
	for (size_type i = 0; i < population; ++i) {
	  cout << "fellow " << i << " : ";
	  for (size_type j = 0; j < nbc; ++j) cout << people[i].situation()[j];
	  cout << " residue : " << people[i].residue << endl;
	}

	// new generation
	for (size_type i = first_computed; i < population; ++i) {
	  int n = (rand() & 255);
	  if (n < 230) { // Evolution by mixing
	    size_type a = (rand() % first_computed);
	    size_type b = (rand() % first_computed);
	    while (a == b) b = (rand() % first_computed);
	    // cout << "Mixte between " << a << " and " << b << endl;
	    for (size_type j = 0; j < nbc; ++j)
	      people[i].situation()[j] = (rand() & 1) ?
		people[a].situation()[j] :
		people[b].situation()[j];
 	    if (n < 80) // Mutations
 	      for (size_type j = 0; j < nbc; ++j) {
		// cout << "mutation" << endl;
		if ((rand() & 15) == 0) people[i].situation()[j] = (rand()&3);
	      }
	  } else { // Evolution by Newton iterations
	    cout << "Newton\n";
	    size_type a = (rand() % first_computed);
	    FRICTION2.change_situation(people[a].situation());
	    gmm::iteration iter(residue, noisy, 40000);
	    getfem::standard_solve(MS2, FRICTION2, iter);  
	    gmm::copy(MS2.state(), MS.state());
	    gmm::iteration iterbis(residue, noisy, 2);
	    getfem::standard_solve(MS, FRICTION, iterbis);
	    situation_of(gmm::sub_vector(MS.state(),
					 gmm::sub_interval(0, mf_u.nb_dof())),
			 gap, BN, BT, people[i].situation());
	  }
	  
	}
	
	++nbiter;
	// if (nbiter == 200) nbiter = 0;
	cout << "iter = " << nbiter << endl;
	// getchar();
	
	
      }

      // free memory
      for (size_type i = 0; i < population; ++i)
	delete people[i].situation_;

    }
    break;
  }

  std::vector<int> situation1(nbc);
  situation_of(gmm::sub_vector(MS.state(),
			       gmm::sub_interval(0, mf_u.nb_dof())),
	       gap, BN, BT, situation1);
  cout << "situation of solution : ";
  for (size_type j = 0; j < nbc; ++j) cout << situation1[j];
  cout << endl;
  cout << "Final residu : " <<  MS.reduced_residu_norm() << endl;
  cout << "Norm of solution : " << gmm::vect_norm2(U) << endl;

  // gmm::copy(FRICTION.get_LN(MS), LN1);
  // gmm::copy(FRICTION.get_LT(MS), LT1);
    
  if (dxexport) {
    getfem::dx_export exp(datafilename + ".dx", false);
    getfem::stored_mesh_slice sl;
    if (N <= 2)
      sl.build(mesh, getfem::slicer_none(),4);
    else
      sl.build(mesh, getfem::slicer_boundary(mesh),4);
    exp.exporting(sl,true);
    exp.exporting_mesh_edges();
    exp.write_point_data(mf_u, U);
    exp.serie_add_object("deformationsteps");
    std::vector<scalar_type> VM;
    calcul_von_mises(mf_u, U, mf_vm, VM, mu);
    exp.write_point_data(mf_vm, VM); 
    exp.serie_add_object("vonmisessteps");    
  }
}

  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {
  dal::exception_callback_debug cb;
  dal::exception_callback::set_exception_callback(&cb); // in order to debug

#ifdef GETFEM_HAVE_FEENABLEEXCEPT /* trap SIGFPE */
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif

  try {    
    friction_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.solve();
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0; 
}
