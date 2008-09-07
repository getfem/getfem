// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2008 Yves Renard.
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

/**
 * Dynamic friction in linear elasticity.
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
 */


#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_import.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_Coulomb_friction.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_derivatives.h"
#include "gmm/gmm.h"
#include <fstream>

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

/*
 * structure for the friction problem
 */
struct friction_problem {

  enum {
    DIRICHLET_BOUNDARY, CONTACT_BOUNDARY
  };
  getfem::mesh mesh;  /* the mesh */
  getfem::mesh_im  mim;      /* integration methods.                         */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the friction solution     */
  getfem::mesh_fem mf_v;     /* main mesh_fem, for the friction solution     */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_vm;    /* mesh_fem used for the VonMises stress        */
  scalar_type lambda, mu;    /* Lamé coefficients.                           */
  scalar_type rho, PG;       /* density, and gravity                         */
  scalar_type friction_coef; /* friction coefficient.                        */

  scalar_type residual;      /* max residual for the iterative solvers       */
  size_type N, noisy;
  scalar_type T, dt, r;
  scalar_type init_vert_pos, init_vert_speed, hspeed, dtexport;
  scalar_type Dirichlet_ratio;
  bool  dxexport, Dirichlet;

  std::string datafilename;
  bgeot::md_param PARAM;

  void solve(void);
  void init(void);
  friction_problem(void) : mim(mesh), mf_u(mesh), mf_v(mesh),
			   mf_rhs(mesh), mf_vm(mesh) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void friction_problem::init(void) {
  std::string FEM_TYPE_U  = PARAM.string_value("FEM_TYPE_U","FEM name");
  std::string FEM_TYPE_V  = PARAM.string_value("FEM_TYPE_V","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "FEM_TYPE_U = "  << FEM_TYPE_U << "\n";
  cout << "FEM_TYPE_V = "  << FEM_TYPE_V << "\n";
  cout << "INTEGRATION = " << INTEGRATION << "\n";

  std::string meshname
    (PARAM.string_value("MESHNAME", "Nom du fichier de maillage"));

  getfem::import_mesh(meshname, mesh);
  N = mesh.dim();
  mesh.optimize_structure();

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL");
  if (residual == 0.) residual = 1e-10;

  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  rho = PARAM.real_value("RHO", "Density");
  PG = PARAM.real_value("PG", "Gravity constant");
  friction_coef = PARAM.real_value("FRICTION_COEF", "Friction coefficient");
  GMM_ASSERT1(friction_coef == 0., "Sorry, friction desactived");

  Dirichlet = PARAM.int_value("DIRICHLET","Dirichlet condition or not");
  Dirichlet_ratio = PARAM.real_value("DIRICHLET_RATIO",
				     "parameter for Dirichlet condition");
  T = PARAM.real_value("T", "from [0,T] the time interval");
  dt = PARAM.real_value("DT", "time step");
  dxexport = (PARAM.int_value("DX_EXPORT", "Exporting on OpenDX format")
	      != 0);
  dtexport = PARAM.real_value("DT_EXPORT", "time step for the export");
  dtexport = dt * double(int(dtexport / dt + 0.5));
 
  r = PARAM.real_value("R", "augmentation parameter");
  noisy = (PARAM.int_value("NOISY", "verbosity of iterative methods") != 0);
  init_vert_pos = PARAM.real_value("INIT_VERT_POS", "initial position");
  init_vert_speed = PARAM.real_value("INIT_VERT_SPEED","initial speed");
  hspeed = PARAM.real_value("FOUNDATION_HSPEED","initial speed");

  mf_u.set_qdim(dim_type(N));
  mf_v.set_qdim(dim_type(N));
  

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE_U);
  getfem::pfem pf_v = getfem::fem_descriptor(FEM_TYPE_V);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_v.set_finite_element(mesh.convex_index(), pf_v);
  mf_vm.set_classical_discontinuous_finite_element(1);
  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    if (!pf_u->is_lagrange()) {
      GMM_ASSERT1(false, "You are using a non-lagrange FEM. "
		  << "In that case you need to set "
		  << "DATA_FEM_TYPE in the .param file");
    }
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }
  
  /* set boundary conditions */
  base_node center(0.,0.,20.);
  std::cout << "Reperage des bord de contact et Dirichlet\n";  
  for (dal::bv_visitor cv(mesh.convex_index()); !cv.finished(); ++cv) {
    short_type nf = mesh.structure_of_convex(cv)->nb_faces();
    for (short_type f = 0; f < nf; f++) {
      if (!mesh.is_convex_having_neighbour(cv, f)) {
	base_small_vector un = mesh.normal_of_face_of_convex(cv, f);
	un /= gmm::vect_norm2(un);	
	base_node pt = mesh.points_of_face_of_convex(cv,f)[0];
	if (un[N-1] < -0.000001 && (N != 3 || (gmm::vect_dist2(pt, center)
			   > .99*sqrt(25. + 15*15) && pt[N-1] < 20.1)))
	  mesh.region(CONTACT_BOUNDARY).add(cv, f);
	if (un[N-1] > 0.1 && Dirichlet)
	  mesh.region(DIRICHLET_BOUNDARY).add(cv, f);
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
/*  Brique dynamique mixte.                                               */
/**************************************************************************/

namespace getfem {

# define MDBRICK_MIXED_DYNAMIC 738053

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_mixed_dynamic : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    const mesh_fem *mf_u;
     const mesh_fem &mf_v;
    mdbrick_parameter<VECTOR> RHO_;
    VECTOR DFU, DFV;
    T_MATRIX B_, C_;
    size_type num_fem;
    value_type Bcoef, Kcoef, Ccoef;
    std::set<size_type> boundary_sup;
    bool M_uptodate;

    virtual void proper_update(void) {
      mf_u = this->mesh_fems[num_fem];
      M_uptodate = false;
    }

    void proper_update_M(void) {
      GMM_TRACE2("Assembling mass matrices for mdbrick_mixed_dynamic");
      gmm::clear(B_);
      gmm::resize(B_, mf_v.nb_dof(), mf_u->nb_dof());
      asm_mass_matrix_param(B_, *(this->mesh_ims[0]), mf_v, *mf_u,
			    RHO_.mf(), RHO_.get());
      gmm::clear(C_);
      gmm::resize(C_, mf_v.nb_dof(), mf_v.nb_dof());
      asm_mass_matrix_param(C_, *(this->mesh_ims[0]), mf_v,
			    RHO_.mf(), RHO_.get());
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), mf_v.nb_dof());
      gmm::sub_interval SUBJ(i0+this->mesh_fem_positions[num_fem],
			     mf_u->nb_dof());

      if (Kcoef != value_type(1)) gmm::scale(MS.tangent_matrix(), Kcoef);
      gmm::copy(gmm::scaled(get_B(), Bcoef),
		gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
      gmm::copy(gmm::transposed(gmm::scaled(get_B(), Bcoef)),
		gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
      gmm::copy(gmm::scaled(get_C(), -Ccoef),
		gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
    }
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), mf_v.nb_dof());
      gmm::sub_interval SUBJ(i0+this->mesh_fem_positions[num_fem],
			     mf_u->nb_dof());

      if (Kcoef != value_type(1))  gmm::scale(MS.residual(), Kcoef);
      gmm::mult(get_B(), gmm::scaled(gmm::sub_vector(MS.state(), SUBJ), Bcoef),
		gmm::sub_vector(MS.residual(), SUBI));
      gmm::mult_add(gmm::transposed(get_B()),
		    gmm::scaled(gmm::sub_vector(MS.state(), SUBI), Bcoef),
		    gmm::sub_vector(MS.residual(), SUBJ));
      gmm::mult_add(get_C(),
		    gmm::scaled(gmm::sub_vector(MS.state(), SUBI), -Ccoef),
		    gmm::sub_vector(MS.residual(), SUBI));
      if (gmm::vect_size(DFU) > 0)
	gmm::add(gmm::scaled(DFU, -value_type(1)),
		 gmm::sub_vector(MS.residual(), SUBJ));
      if (gmm::vect_size(DFV) > 0)
	gmm::add(gmm::scaled(DFV, -value_type(1)),
		 gmm::sub_vector(MS.residual(), SUBI));
    }

    void set_dynamic_coeff(value_type a, value_type b, value_type c)
    { Bcoef=a; Kcoef=b; Ccoef = c; }
    template <class VEC> void set_DFU(const VEC &DF_)
    { gmm::resize(DFU, gmm::vect_size(DF_)); gmm::copy(DF_, DFU); }
    template <class VEC> void set_DFV(const VEC &DF_)
    { gmm::resize(DFV, gmm::vect_size(DF_)); gmm::copy(DF_, DFV); }

    const T_MATRIX &get_B(void) {
      this->context_check();
      if (!M_uptodate || this->parameters_is_any_modified()) {
	proper_update_M();
	M_uptodate = true;
	this->parameters_set_uptodate();
      }
      return B_; 
    }

    const T_MATRIX &get_C(void) {
      this->context_check();
      if (!M_uptodate || this->parameters_is_any_modified()) {
	proper_update_M();
	M_uptodate = true;
	this->parameters_set_uptodate();
      }
      return C_; 
    }

    /** provide access to the value of the solution corresponding to
	the local mesh_fem mf_v. */
    SUBVECTOR get_V(MODEL_STATE &MS) {
      gmm::sub_interval SUBU = gmm::sub_interval(this->first_index()
				     + sub_problem.nb_dof(), mf_v.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    /**
       @param num_fem_ the mesh_fem number on which this brick is is applied.
    */
    mdbrick_mixed_dynamic(mdbrick_abstract<MODEL_STATE> &problem,
			  const mesh_fem &mf_v_, value_type RHO__, 
			  size_type num_fem_=0)
      : sub_problem(problem), mf_v(mf_v_), RHO_("rho", this),
	num_fem(num_fem_) {
      
      Bcoef = Kcoef = Ccoef = value_type(1);
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = false;
      this->add_proper_mesh_fem(mf_v, MDBRICK_MIXED_DYNAMIC);
      this->force_update();
      
      RHO_.set(classical_mesh_fem(mf_u->linked_mesh(), 0), RHO__);
    }
  };



}



/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

void friction_problem::solve(void) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  N = mesh.dim();
  cout << "Number of dof for u: " << mf_u.nb_dof() << endl;
  cout << "Number of dof for v: " << mf_v.nb_dof() << endl;

  size_type ref_dof = 0;
  for (size_type i = 1; i < mf_u.nb_dof(); ++i)
    if (mf_u.point_of_dof(i)[N-1] < mf_u.point_of_dof(ref_dof)[N-1])
      ref_dof = i;

  // Linearized elasticity brick.
  getfem::mdbrick_isotropic_linearized_elasticity<>
    ELAS(mim, mf_u, lambda, mu);

  scalar_type h = mesh.minimal_convex_radius_estimate();
  cout << "minimal convex radius estimate : " << h << endl;
  cout << "CFL estimate = "
       << h / gmm::sqrt((lambda + 2.0 * mu) / rho) << endl;

  // Defining the volumic source term.
  plain_vector F(nb_dof_rhs * N);
  plain_vector f(N); f[N-1] = -rho*PG;
  for (size_type i = 0; i < nb_dof_rhs; ++i)
      gmm::copy(f,gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  
  // Volumic source term brick.
  getfem::mdbrick_source_term<> VOL_F(ELAS, mf_rhs, F);

  // Dirichlet condition brick.
  gmm::clear(F);
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    F[(i+1)*N-1] = Dirichlet_ratio * mf_rhs.point_of_dof(i)[N-1];
  getfem::mdbrick_Dirichlet<> DIRICHLET(VOL_F, DIRICHLET_BOUNDARY);
  DIRICHLET.rhs().set(mf_rhs, F);
  
  // contact condition for Lagrange elements
  dal::bit_vector cn = mf_u.dof_on_set(CONTACT_BOUNDARY);
  cout << "cn = " << cn << endl;
  sparse_matrix BN(cn.card()/N, mf_u.nb_dof());
  sparse_matrix BT((N-1)*cn.card()/N, mf_u.nb_dof());
  plain_vector gap(cn.card()/N);
  size_type jj = 0;
  for (dal::bv_visitor i(cn); !i.finished(); ++i)
    if (i % N == 0) {
      BN(jj, i+N-1) = -1.;
      gap[jj] = mf_u.point_of_dof(i)[N-1];
      for (size_type k = 0; k < N-1; ++k) BT((N-1)*jj+k, i+k) = 1.;
      ++jj;
    }

  getfem::mdbrick_Coulomb_friction<>
    // FRICTION(DIRICHLET, BN, gap, friction_coef, BT);
    FRICTION(DIRICHLET, BN, gap);
  FRICTION.set_r(r);

  // Dynamic brick.
  getfem::mdbrick_mixed_dynamic<> DYNAMIC(FRICTION, mf_v, rho);
  
  cout << "Total number of variables: " << DYNAMIC.nb_dof() << endl;
  getfem::standard_model_state MS(DYNAMIC);

  // cout << "C = " << DYNAMIC.get_C() << endl;

  plain_vector WT(mf_u.nb_dof()), WN(mf_u.nb_dof());
  plain_vector DFU(mf_u.nb_dof()), DFV(mf_v.nb_dof());
  plain_vector HSPEED(mf_u.nb_dof());
  plain_vector U0(mf_u.nb_dof()), V0(mf_v.nb_dof());
  plain_vector U1(mf_u.nb_dof()), V1(mf_v.nb_dof());
  plain_vector Udemi(mf_u.nb_dof()), Vdemi(mf_v.nb_dof());
  plain_vector LT0(gmm::mat_nrows(BT)), LN0(gmm::mat_nrows(BN));
  plain_vector LT1(gmm::mat_nrows(BT)), LN1(gmm::mat_nrows(BN));
  scalar_type t(0), t_export(dtexport), beta_(0);
  scalar_type alpha_(0);

  // Initial conditions (U0, V0)
  gmm::clear(U0); gmm::clear(V0); gmm::clear(LT0);
  for (size_type i=0; i < mf_u.nb_dof(); ++i)
    if ((i % N) == 0) { 
      U0[i+N-1] = Dirichlet ? (Dirichlet_ratio * mf_u.point_of_dof(i)[N-1])
	: init_vert_pos;
      HSPEED[i] = hspeed;
    }

  for (size_type i=0; i < mf_v.nb_dof(); ++i)
    if ((i % N) == 0) V0[i+N-1] = Dirichlet ? 0.0 : init_vert_speed;
  
  gmm::iteration iter(residual, 0, 40000);
  iter.set_noisy(int(noisy));

  scalar_type J0 = 0.5*gmm::vect_sp(ELAS.get_K(), U0, U0)
    + 0.5 * gmm::vect_sp(DYNAMIC.get_C(), V0, V0)
    - gmm::vect_sp(VOL_F.get_F(), U0);

  std::auto_ptr<getfem::dx_export> exp;
  getfem::stored_mesh_slice sl;
  if (dxexport) {
    exp.reset(new getfem::dx_export(datafilename + ".dx", false));
    if (N <= 2)
      sl.build(mesh, getfem::slicer_none(),4);
    else
      sl.build(mesh, getfem::slicer_boundary(mesh),4);
    exp->exporting(sl,true);
    exp->exporting_mesh_edges();
    exp->write_point_data(mf_u, U0, "stepinit"); 
    exp->serie_add_object("deformationsteps");
    std::vector<scalar_type> VM;
    calcul_von_mises(mf_u, U0, mf_vm, VM, mu);
    exp->write_point_data(mf_vm, VM, "vonmises"); 
    exp->serie_add_object("vonmisessteps");
  }
  
  std::ofstream fileout1("time", std::ios::out);   
  std::ofstream fileout2("energy", std::ios::out);
 
  scalar_type Einit = (0.5*gmm::vect_sp(ELAS.get_K(), U0, U0)
		       + 0.5 * gmm::vect_sp(DYNAMIC.get_C(), V0, V0)
		       - gmm::vect_sp(VOL_F.get_F(), U0));
  cout << "t=0, initial energy: " << Einit << endl;
  
  while (t <= T) {

    beta_ = 2./dt; alpha_ = 1.;
    gmm::mult(gmm::transposed(DYNAMIC.get_B()), gmm::scaled(V0, 2./dt), DFU);
    gmm::mult(DYNAMIC.get_B(), gmm::scaled(U0, 2./dt), DFV);
    gmm::copy(gmm::scaled(U0, -1.), WT);
    gmm::clear(WN);
    gmm::add(gmm::scaled(HSPEED, -1./beta_), WT);

    FRICTION.set_WN(WN); FRICTION.set_r(r); 
    // FRICTION.set_WT(WT); FRICTION.set_beta(beta_);
    FRICTION.set_alpha(alpha_); 
    DYNAMIC.set_dynamic_coeff(2./dt, 1., 1.);
    DYNAMIC.set_DFU(DFU);
    DYNAMIC.set_DFV(DFV);
    
    iter.init();
    gmm::default_newton_line_search ls(size_type(-1), 4.0/3.0,
				       1.0/20.0, 9.0/10.0, 1.1);
    getfem::standard_solve(MS, DYNAMIC, iter,
			   getfem::default_linear_solver(DYNAMIC), ls);
    
    
    gmm::copy(ELAS.get_solution(MS), Udemi);
    gmm::copy(DYNAMIC.get_V(MS), Vdemi);
    gmm::copy(FRICTION.get_LN(MS), LN1);
    gmm::copy(FRICTION.get_LT(MS), LT1);

    gmm::add(gmm::scaled(V0, -1.), gmm::scaled(Vdemi, 2.), V1);
    gmm::add(gmm::scaled(U0, -1.), gmm::scaled(Udemi, 2.), U1);
    
    scalar_type J1 =  0.5*gmm::vect_sp(ELAS.get_K(), U1, U1)
      + 0.5 * gmm::vect_sp(DYNAMIC.get_C(), V1, V1)
      - gmm::vect_sp(VOL_F.get_F(), U1);

    scalar_type Jdemi = 0.5*gmm::vect_sp(ELAS.get_K(), Udemi, Udemi)
      + 0.5 * gmm::vect_sp(DYNAMIC.get_C(), Vdemi, Vdemi)
      - gmm::vect_sp(VOL_F.get_F(), Udemi);

    t += dt;
    
    size_type nbsl= 0, nbst = 0;
    for (size_type i = 0; i < gmm::mat_nrows(BN); ++i) {
      if (LN1[i] < -1E-12) {
	if (gmm::vect_norm2(gmm::sub_vector(LT1,
					    gmm::sub_interval(i*(N-1), N-1)))
	    < -friction_coef * double(LN1[i]) * 0.99999)
	  nbst++; else nbsl++;
      }
    }
    
    cout << "t = " << t << " J1 = " << J1 << " Jdemi = " << Jdemi;
    cout << " (st " << nbst << ", sl " << nbsl << ")" << endl;
    
    gmm::copy(U1, U0); gmm::copy(V1, V0);  J0 = J1;
    gmm::copy(LN1, LN0); gmm::copy(LT1, LT0);
    if (dxexport && t >= t_export-dt/20.0) {
      
      fileout1 << t << "\n";
      fileout2 << J1   << "\n";	
      
      exp->write_point_data(mf_u, U0);
      exp->serie_add_object("deformationsteps");
      std::vector<scalar_type> VM;
      calcul_von_mises(mf_u, U0, mf_vm, VM, mu);
      exp->write_point_data(mf_vm, VM); 
      exp->serie_add_object("vonmisessteps");
      
      t_export += dtexport;
    }
  }
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {    
    friction_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.solve();
  }
  GMM_STANDARD_CATCH_ERROR;

  cout << "To see the simulation, you have to set DX_EXPORT to 1 in "
    "dynamic_friction.param and DT_EXPORT to a suitable value (for "
    "instance equal to DT). Then you can use Open_DX (type just \"dx\" "
    "if it is installed on your system) with the Visual Program "
    "dynamic_friction.net (use for instance \"Edit Visual Programs ...\" "
    "with dynamic_friction.net, then \"execute once\" in Execute menu and "
    "use the sequencer to start the animation).\n";

  return 0; 
}
