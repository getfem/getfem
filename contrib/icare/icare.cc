// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2008 Michel Fournié, Julien Pommier,
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

/**@file icare.cc
   @brief Fluid flow (Navier-Stokes) around an obstacle.
*/

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_import.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_Navier_Stokes.h"
#include "icare.h"

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrices. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

enum {
  DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM,
  NORMAL_PART_DIRICHLET_BOUNDARY_NUM,
  NONREFLECTIVE_BOUNDARY_NUM
};


struct problem_definition;

/*
 * structure for the navier_stokes problem
 */
struct navier_stokes_problem {

  getfem::mesh mesh;         /* the mesh */
  base_node BBmin, BBmax;    /* bounding box of the mesh */
  getfem::mesh_im  mim;      /* integration methods.                         */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the velocity              */
  getfem::mesh_fem mf_p;     /* mesh_fem for the pressure                    */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_mult;  /* mesh_fem for Dirichlet condition.            */
  scalar_type Re;            /* Reynolds number */
  scalar_type nu;            /* 1/Re */
  scalar_type dt, T, dt_export;
  unsigned N;
  scalar_type residual;      /* max residual for the iterative solvers       */
  int noisy;
  int export_to_opendx;

  std::auto_ptr<getfem::dx_export> exp;
  getfem::stored_mesh_slice sl;
  bool first_export;
  scalar_type t_export;
  void do_export(scalar_type t);

  int option;

  std::auto_ptr<problem_definition> pdef;

  std::string datafilename;
  bgeot::md_param PARAM;

  plain_vector Un1, Un0, Pn1, Pn0; /* U_{n+1}, U_{n}, P_{n+1} and P_{n} */

  base_small_vector sol_f(const base_small_vector &P, scalar_type t);
  base_small_vector Dir_cond(const base_small_vector &P, scalar_type t);

  void solve(void);
  void solve_METHOD_SPLITTING(bool);
  void solve_FULLY_CONSERVATIVE();
  void solve_PREDICTION_CORRECTION();
  void solve_PREDICTION_CORRECTION2();
  void init(void);
  navier_stokes_problem(void) : mim(mesh), mf_u(mesh), mf_p(mesh),
				mf_rhs(mesh), mf_mult(mesh) {}
};

struct problem_definition {
  virtual void choose_boundaries(navier_stokes_problem &p) {
    getfem::outer_faces_of_mesh(p.mesh, p.mesh.region(DIRICHLET_BOUNDARY_NUM));
  }
  virtual void validate_solution(navier_stokes_problem &p, scalar_type t) {
    plain_vector R; dirichlet_condition(p, t, R);
    p.mf_rhs.set_qdim(p.N);
    scalar_type err = getfem::asm_L2_dist(p.mim, 
					  p.mf_u, p.Un1,
					  p.mf_rhs, R);
    p.mf_rhs.set_qdim(1);
    cout << " L2 error(t=" << t <<") : " << err << "\n";
  }
  virtual void dirichlet_condition(navier_stokes_problem &p,
				   const base_node &x, scalar_type t,
				   base_small_vector &r) = 0;
  virtual void source_term(navier_stokes_problem &p,
			   const base_node &x, scalar_type t,
			   base_small_vector &F) = 0;
  virtual scalar_type initial_pressure(navier_stokes_problem &, const base_node &) {
    return 0.;
  }
  virtual base_small_vector initial_velocity(navier_stokes_problem &p, const base_node &P) {
    base_small_vector r; dirichlet_condition(p,P,0,r); return r;
  }
  void dirichlet_condition(navier_stokes_problem &p, scalar_type t, 
			   plain_vector &R) {
    unsigned N = p.mesh.dim();
    gmm::resize(R, N*p.mf_rhs.nb_dof());
    base_small_vector r; 
    for (unsigned i=0; i < p.mf_rhs.nb_dof(); ++i) {
      dirichlet_condition(p, p.mf_rhs.point_of_dof(i), t, r);
      gmm::copy(r, gmm::sub_vector(R, gmm::sub_interval(i*N, N)));
    }
  }
  void source_term(navier_stokes_problem &p, scalar_type t, plain_vector &F) {
    unsigned N = p.mesh.dim();
    gmm::resize(F, N*p.mf_rhs.nb_dof());
    base_small_vector f;
    for (unsigned i=0; i < p.mf_rhs.nb_dof(); ++i) {
      source_term(p, p.mf_rhs.point_of_dof(i), t, f);
      gmm::copy(f, gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
    }
  }
  
  void initial_condition_u(navier_stokes_problem &p, plain_vector &U0) {
    plain_vector R(p.N*p.mf_rhs.nb_dof()), F(p.mf_u.nb_dof());
    for (unsigned i=0; i < p.mf_rhs.nb_dof(); ++i) {
      base_small_vector r = initial_velocity(p, p.mf_rhs.point_of_dof(i));
      gmm::copy(r, gmm::sub_vector(R, gmm::sub_interval(i*p.N, p.N)));
    }
    /* L2 projection from mf_rhs onto mf_u (we cannot interpolate directly onto mf_u
     since it can be non-lagrangian) */
    sparse_matrix M(p.mf_u.nb_dof(), p.mf_u.nb_dof());
    getfem::asm_mass_matrix(M, p.mim, p.mf_u);
    getfem::asm_source_term(F, p.mim, p.mf_u, p.mf_rhs, R);
    gmm::iteration iter(1E-13);
    gmm::cg(M, U0, F, gmm::identity_matrix(), iter);
  }
  void initial_condition_p(navier_stokes_problem &p, plain_vector &P0) {
    plain_vector PP(p.mf_rhs.nb_dof());

    for (unsigned i=0; i < p.mf_rhs.nb_dof(); ++i)
      PP[i] = initial_pressure(p, p.mf_rhs.point_of_dof(i));

    /* L2 projection from mf_rhs onto mf_p */
    plain_vector F(p.mf_p.nb_dof());
    sparse_matrix M(p.mf_p.nb_dof(), p.mf_p.nb_dof());
    getfem::asm_mass_matrix(M, p.mim, p.mf_p);
    getfem::asm_source_term(F, p.mim, p.mf_p, p.mf_rhs, PP);
    gmm::iteration iter(1E-13);
    gmm::cg(M, P0, F, gmm::identity_matrix(), iter);
  }
  virtual ~problem_definition() {}
};

struct problem_definition_Stokes_analytic : public problem_definition {
  virtual void dirichlet_condition(navier_stokes_problem &p,
				   const base_node &P, scalar_type t,
				   base_small_vector &r) {
    r = base_small_vector(p.N);
    scalar_type x = P[0], y = P[1];
    r[0] =  2.*(2.*y-1.)*(1.-1.*gmm::sqr(2.*x-1.))*t;
    r[1] = -2.*(2.*x-1.)*(1.-1.*gmm::sqr(2.*y-1.))*t;
  }
  virtual void source_term(navier_stokes_problem &p,
			   const base_node &P, scalar_type t,
			   base_small_vector &F) {
    scalar_type x = P[0], y = P[1];
    F = base_small_vector(p.N);
    F[0] = -16.*y*x*x+16.*y*x+8.*x*x-8.*x+32.*p.nu*t*y-16.*p.nu*t+8.*t*x*x-8.*t*x;
    F[1] =  16.*x*y*y-16.*y*x-8.*y*y+8.*y-32.*p.nu*t*x+16.*p.nu*t+8.*t*y*y-8.*t*y;
  }
};

struct problem_definition_Green_Taylor_analytic : public problem_definition {
  virtual void dirichlet_condition(navier_stokes_problem &p,
				   const base_node &P, scalar_type t,
				   base_small_vector &r) {
    r = base_small_vector(p.N);
    scalar_type x = P[0], y = P[1];
    r[0] =  -cos(x)*sin(y)*exp(-2*t*p.nu);
    r[1] =  +sin(x)*cos(y)*exp(-2*t*p.nu);
  }
  virtual void source_term(navier_stokes_problem &p,
			   const base_node &P, scalar_type t,
			   base_small_vector &F) {
    scalar_type x = P[0], y = P[1];
    F = base_small_vector(p.N);
    F[0] = -exp(-4.*t*p.nu)*sin(2.*x);
    F[1] = -exp(-4.*t*p.nu)*sin(2.*y);
  }
  virtual scalar_type initial_pressure(navier_stokes_problem &p,
				       const base_node &P) {
    scalar_type x = P[0], y = P[1], t = 0;
    return -1./4 * (cos(2*x)+cos(2*y)) * exp(-4*p.nu * t);
  }
};

struct problem_rotating_cylinder : public problem_definition {
  scalar_type alpha;
  
  virtual void choose_boundaries(navier_stokes_problem &p) {
    getfem::mesh_region r; 
    getfem::outer_faces_of_mesh(p.mesh, r);
    for (getfem::mr_visitor i(r); !i.finished(); ++i) {
      base_node G = gmm::mean_value(p.mesh.points_of_face_of_convex(i.cv(),i.f()));
      /*if (gmm::abs(G[0] - p.BBmax[0]) < 1e-7)
	p.mesh.region(NONREFLECTIVE_BOUNDARY_NUM).add(i.cv(),i.f());
      else if (gmm::abs(G[1] - p.BBmax[1]) < 1e-7
	       || gmm::abs(G[1] - p.BBmin[1]) < 1e-7)
	p.mesh.region(NORMAL_PART_DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());
	else */
	p.mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());
    }
  }
  virtual void dirichlet_condition(navier_stokes_problem &p,
				   const base_node &P, scalar_type /*t*/,
				   base_small_vector &r) {
    r = base_small_vector(p.N);
    scalar_type x = P[0], y = P[1];
    bool on_cylinder = true;
    if (gmm::abs(x-p.BBmin[0]) < 1e-7) on_cylinder = false;
    if (gmm::abs(y-p.BBmin[1]) < 1e-7) on_cylinder = false;
    if (gmm::abs(x-p.BBmax[0]) < 1e-7) on_cylinder = false;
    if (gmm::abs(y-p.BBmax[1]) < 1e-7) on_cylinder = false;
    if (!on_cylinder) {
      //if (!(gmm::abs(y- p.BBmin[1]) < 1e-7 || gmm::abs(y - p.BBmax[1])< 1e-7))
      r[0] = 1;
    } else {
      r[0] = -2.*alpha*y; /* HYPOTHESIS: cylinder centered at (0,0) */
      r[1] = 2.*alpha*x;
    }
  }

   virtual void source_term(navier_stokes_problem &p,
			    const base_node &, scalar_type /*t*/,
			   base_small_vector &F) {
    F = base_small_vector(p.N);
  }
  
  void validate_solution(navier_stokes_problem &p, scalar_type t) {
    cout << "Validate_solution : t = " << t << " , |u| = " 
	 << gmm::vect_norm2(p.Un1) << ", |p|="
	 << gmm::vect_norm2(p.Pn1) << "\n";
  }
  virtual base_small_vector initial_velocity(navier_stokes_problem &,
					     const base_node &) {
    base_small_vector r(2); r[0] = 1; r[1] = 0; return r;
  }
  virtual scalar_type initial_pressure(navier_stokes_problem &,
				       const base_node &) {
    return 0.;
  }
  problem_rotating_cylinder(scalar_type aa) : alpha(aa) {}
};


/*
 * Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void navier_stokes_problem::init(void) {
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string FEM_TYPE_P  = PARAM.string_value("FEM_TYPE_P","FEM name P");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  std::string meshname
    (PARAM.string_value("MESHNAME", "Nom du fichier de maillage"));

  getfem::import_mesh(meshname, mesh);
  N = mesh.dim();

  mesh.bounding_box(BBmin, BBmax);
  cout << "mesh bounding box: " << BBmin << " ... " << BBmax << "\n";
  datafilename = PARAM.string_value("ROOTFILENAME","Data files base name.");
  residual = PARAM.real_value("RESIDUAL"); 
  if (residual == 0.) residual = 1e-10;

  nu = PARAM.real_value("NU", "Viscosity");
  dt = PARAM.real_value("DT", "Time step");
  T = PARAM.real_value("T", "Final time");
  dt_export = PARAM.real_value("DT_EXPORT", "Time step for export");
  noisy = PARAM.int_value("NOISY", "");
  option = PARAM.int_value("OPTION", "option");

  int prob = PARAM.int_value("PROBLEM", "the problem");
  switch (prob) {
    case 1: pdef.reset(new problem_definition_Stokes_analytic); break;
    case 2: pdef.reset(new problem_definition_Green_Taylor_analytic); break;
    case 3: pdef.reset(new problem_rotating_cylinder(PARAM.real_value("CYL_ROT_SPEED"))); break;
    default: GMM_ASSERT1(false, "wrong PROBLEM value");
  }

  export_to_opendx = PARAM.int_value("DX_EXPORT", "");
  first_export = true;

  Re = 1 / nu;
  mf_u.set_qdim(N);

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION); 

  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_p.set_finite_element(mesh.convex_index(),
			  getfem::fem_descriptor(FEM_TYPE_P));

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM "
		<< FEM_TYPE << ". In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }

  std::string mult_fem_name = PARAM.string_value("MULTIPLIER_FEM_TYPE");
  if (mult_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM "
		<< FEM_TYPE << ". In that case you need to set "
		<< "MULTIPLIER_FEM_TYPE in the .param file");
    mf_mult.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_mult.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(mult_fem_name));
  }

  /* set boundary conditions */
  cout << "Choosing boundaries\n";
  pdef->choose_boundaries(*this);
}


void navier_stokes_problem::solve() {
  cout << "Number of dof for u : " << mf_u.nb_dof() << endl;
  cout << "Number of dof for p : " << mf_p.nb_dof() << endl;
  gmm::resize(Un0, mf_u.nb_dof());
  gmm::resize(Un1, mf_u.nb_dof());
  gmm::resize(Pn0, mf_p.nb_dof());
  gmm::resize(Pn1, mf_p.nb_dof());
  switch (option) {
  case 0 : solve_METHOD_SPLITTING(true); break;
  case 1 : solve_METHOD_SPLITTING(false); break;
  case 2 : solve_FULLY_CONSERVATIVE(); break;
  case 3 : solve_PREDICTION_CORRECTION(); break;
  case 4 : solve_PREDICTION_CORRECTION2(); break;
  default: GMM_ASSERT1(false, "unknown method");
  }
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

void navier_stokes_problem::solve_METHOD_SPLITTING(bool stokes_only) {
  // in the splitting method, there are TWO problems:

  //
  // definition of the first problem 
  //

  // Velocity brick.
  getfem::mdbrick_abstract<> *previous;
  getfem::mdbrick_generic_elliptic<> velocity(mim, mf_u, nu);
  previous = &velocity;

  
  std::auto_ptr<getfem::mdbrick_NS_uuT<> > velocity_nonlin;
  if (!stokes_only) {
    velocity_nonlin.reset(new getfem::mdbrick_NS_uuT<>(velocity));
    previous = velocity_nonlin.get();
  }
  
  // Volumic source term
  getfem::mdbrick_source_term<> velocity_f(*previous, 
					   mf_rhs, 
					   plain_vector(mf_rhs.nb_dof()*N));

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<> velocity_dir(velocity_f, DIRICHLET_BOUNDARY_NUM);
 
  // Normal part Dirichlet condition brick.
  getfem::mdbrick_normal_component_Dirichlet<>
    velocity_dirnp(velocity_dir, NORMAL_PART_DIRICHLET_BOUNDARY_NUM, mf_rhs);

  // Non-reflective condition brick
  getfem::mdbrick_NS_nonref1<> nonreflective(velocity_dirnp,
					     NONREFLECTIVE_BOUNDARY_NUM, dt);
 
  // Dynamic brick.
  getfem::mdbrick_dynamic<> velocity_dyn(nonreflective, 1.);
  velocity_dyn.set_dynamic_coeff(1.0/dt, 1.0);

  // 
  // definition of the second problem
  //

  getfem::mdbrick_mass_matrix<> mixed(mim, mf_u, 1./dt);
  
  // Pressure term
  getfem::mdbrick_linear_incomp<> mixed_p(mixed, mf_p);

  // Condition on the pressure
  sparse_matrix G(1, mf_p.nb_dof());
  G(0,0) = 1.;
  plain_vector gr(1);
  getfem::mdbrick_constraint<> set_pressure(mixed_p, 1);
  set_pressure.set_constraints(G, gr);
  set_pressure.set_constraints_type(getfem::AUGMENTED_CONSTRAINTS);
    
  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<> mixed_dir(set_pressure, DIRICHLET_BOUNDARY_NUM);
  
  // Normal part Dirichlet condition brick.
  getfem::mdbrick_normal_component_Dirichlet<>
    mixed_dirnp(mixed_dir, NORMAL_PART_DIRICHLET_BOUNDARY_NUM, mf_rhs);

  // Dynamic brick.
  getfem::mdbrick_dynamic<> mixed_dyn(mixed_dirnp, 1.);
  mixed_dyn.set_dynamic_coeff(0.0, 1.0);

  // 
  // dynamic problem
  //

  plain_vector DF(mf_u.nb_dof()), F(mf_rhs.nb_dof());

  pdef->initial_condition_u(*this, Un0);
  
  gmm::iteration iter(residual, noisy);
  getfem::standard_model_state MSL(velocity_dyn), MSM(mixed_dyn);
  
  do_export(0);
  for (scalar_type t = dt; t <= T; t += dt) {

    if (!stokes_only) velocity_nonlin->set_U0(Un0);
    pdef->source_term(*this, t, F);
    velocity_f.source_term().set(F);
    pdef->dirichlet_condition(*this, t, F);
    velocity_dir.rhs().set(mf_rhs, F);
    nonreflective.set_Un(Un0);
    
    gmm::mult(velocity_dyn.get_M(), gmm::scaled(Un0, 1./dt), DF);
    velocity_dyn.set_DF(DF);
    iter.init();
    getfem::standard_solve(MSL, velocity_dyn, iter);

    gmm::copy(velocity.get_solution(MSL), Un1);
    gmm::mult(velocity_dyn.get_M(), gmm::scaled(Un1, 1./dt), DF);
    mixed_dyn.set_DF(DF);
    mixed_dir.rhs().set(mf_rhs, F);
    iter.init();
    getfem::standard_solve(MSM, mixed_dyn, iter);
    gmm::copy(mixed.get_solution(MSM), Un1);
    gmm::copy(mixed_p.get_pressure(MSM), Pn1);

    pdef->validate_solution(*this, t);

    gmm::copy(Un1, Un0); gmm::copy(Pn1, Pn0);
    do_export(t);
  }
}



void navier_stokes_problem::solve_FULLY_CONSERVATIVE() {

  // Velocity brick.    
  getfem::mdbrick_navier_stokes<> velocity(mim, mf_u, mf_p, nu);
  // Condition on the pressure
  sparse_matrix G(1, mf_p.nb_dof());
  G(0,0) = 1.;
  plain_vector gr(1);
  getfem::mdbrick_constraint<> velocity_ctr(velocity);
  velocity_ctr.set_constraints(G, gr);
  velocity_ctr.set_constraints_type(getfem::AUGMENTED_CONSTRAINTS);


  // Volumic source term
  getfem::mdbrick_source_term<> velocity_f(velocity_ctr, 
					   mf_rhs, 
					   plain_vector(mf_rhs.nb_dof()*N));

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<> velocity_dir(velocity_f, DIRICHLET_BOUNDARY_NUM);
  
  // Normal part Dirichlet condition brick.
  getfem::mdbrick_normal_component_Dirichlet<>
    velocity_dirnp(velocity_dir, NORMAL_PART_DIRICHLET_BOUNDARY_NUM, mf_rhs);

  // Non-reflective condition brick
  getfem::mdbrick_NS_nonref1<>
    nonreflective(velocity_dirnp, NONREFLECTIVE_BOUNDARY_NUM, dt);

  // Dynamic brick.
  getfem::mdbrick_dynamic<> velocity_dyn(nonreflective, 1.);
  velocity_dyn.set_dynamic_coeff(1.0/dt, 1.0);

  // 
  // dynamic problem
  //

  plain_vector DF(mf_u.nb_dof()), F(mf_rhs.nb_dof());

  pdef->initial_condition_u(*this, Un0);
  
  gmm::iteration iter(residual, noisy);
  getfem::standard_model_state MSL(velocity_dyn);
  
  do_export(0);
  for (scalar_type t = dt; t <= T; t += dt) {

  
    pdef->source_term(*this, t, F);
    velocity_f.source_term().set(F);
    pdef->dirichlet_condition(*this, t, F);
    velocity_dir.rhs().set(mf_rhs, F);
    nonreflective.set_Un(Un0);

    gmm::mult(velocity_dyn.get_M(), gmm::scaled(Un0, 1./dt), DF);
    velocity_dyn.set_DF(DF);
    iter.init();
    getfem::standard_solve(MSL, velocity_dyn, iter);

    gmm::copy(velocity.get_velocity(MSL), Un1);
    gmm::copy(velocity.get_pressure(MSL), Pn1);
   
    pdef->validate_solution(*this, t);     

    gmm::copy(Un1, Un0); gmm::copy(Pn1, Pn0);
    do_export(t);
  }
}

/************************************************************/
void navier_stokes_problem::solve_PREDICTION_CORRECTION() {
  // Velocity brick.  
  getfem::mdbrick_generic_elliptic<> velocity(mim, mf_u, nu);
  getfem::mdbrick_NS_uuT<> velocity_nonlin(velocity);

  // Volumic source term
  getfem::mdbrick_source_term<> velocity_f(velocity_nonlin, 
					   mf_rhs, 
					   plain_vector(mf_rhs.nb_dof()*N));

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<> velocity_dir(velocity_f, DIRICHLET_BOUNDARY_NUM);
  
  // Normal part Dirichlet condition brick.
  getfem::mdbrick_normal_component_Dirichlet<>
    velocity_dirnp(velocity_dir, NORMAL_PART_DIRICHLET_BOUNDARY_NUM, mf_rhs);
  
  // Non-reflective condition brick
  getfem::mdbrick_NS_nonref1<> nonreflective( velocity_dirnp, 
					      NONREFLECTIVE_BOUNDARY_NUM, dt);
  

  // Dynamic brick.
    getfem::mdbrick_dynamic<> velocity_dyn(nonreflective, 1.);
  velocity_dyn.set_dynamic_coeff(1.0/dt, 1.0);

  // 
  // Poisson problem for prediction correction scheme
  //
  
  getfem::mdbrick_generic_elliptic<> poisson(mim, mf_p, 1.0);
  sparse_matrix G(1, mf_p.nb_dof()); G(0,0) = 1.;  
  getfem::mdbrick_constraint<> poisson_setonedof(poisson);
  poisson_setonedof.set_constraints(G, plain_vector(1));
  poisson_setonedof.set_constraints_type(getfem::ELIMINATED_CONSTRAINTS);

  getfem::mdbrick_source_term<>
    poisson_source(poisson_setonedof, mf_rhs, plain_vector(mf_rhs.nb_dof()));

  sparse_matrix B(mf_p.nb_dof(), mf_u.nb_dof());
  asm_stokes_B(B, mim, mf_u, mf_p);

  // 
  // dynamic problem
  //
  plain_vector DF(mf_u.nb_dof()), F(mf_rhs.nb_dof()), 
    USTAR(mf_u.nb_dof()), USTARbis(mf_u.nb_dof());

  pdef->initial_condition_u(*this, Un0);
  pdef->initial_condition_p(*this, Pn0);

  gmm::iteration iter(residual, noisy);
  getfem::standard_model_state MSL(velocity_dyn), MSM(poisson_source);
  
  do_export(0);
  for (scalar_type t = dt; t <= T; t += dt) {

    velocity_nonlin.set_U0(Un0);
    nonreflective.set_Un(Un0);
    pdef->source_term(*this, t, F);
    velocity_f.source_term().set(F);
    pdef->dirichlet_condition(*this, t, F);
    velocity_dir.rhs().set(mf_rhs, F);
    
    gmm::mult(velocity_dyn.get_M(), gmm::scaled(Un0, 1./dt), DF);
    gmm::mult_add(gmm::transposed(B), gmm::scaled(Pn0, -1.), DF);
    velocity_dyn.set_DF(DF);
    iter.init();
    getfem::standard_solve(MSL, velocity_dyn, iter);

    gmm::copy(velocity.get_solution(MSL), USTAR);
    gmm::mult(B, USTAR, Pn1);
 
    cout << "PN1 : " << gmm::vect_norm2(Pn1) << endl;
    cout << "U1 - USTAR = " << gmm::vect_dist2(Un0, USTAR);

    poisson_source.set_auxF(Pn1);
    iter.init();
    getfem::standard_solve(MSM, poisson_source, iter);
    gmm::mult(gmm::transposed(B), gmm::scaled(poisson.get_solution(MSM), -1.),
	      USTARbis);

    gmm::iteration iter2 = iter; iter2.reduce_noisy(); iter2.init();
    gmm::cg(velocity_dyn.get_M(), Un1, USTARbis,
	    gmm::identity_matrix(), iter2);

    gmm::add(USTAR, Un1);
    gmm::add(gmm::scaled(poisson.get_solution(MSM), 1./dt), Pn0, Pn1);
    
    pdef->validate_solution(*this, t);
    cout << "U1 - Un1 = " << gmm::vect_dist2(Un1, Un0);
    
    gmm::copy(Un1, Un0); gmm::copy(Pn1, Pn0);
    do_export(t);

 }
}


/************************************************************/
void navier_stokes_problem::solve_PREDICTION_CORRECTION2() {

  size_type nbdof_u = mf_u.nb_dof(), nbdof_p = mf_p.nb_dof();
  size_type nbdof_rhs = mf_rhs.nb_dof();
  getfem::mesh_region mpirg = mf_u.linked_mesh().get_mpi_region();
  gmm::sub_interval I1(0, nbdof_u);

  cout << "nbdof rhs = " << nbdof_rhs << endl;

  // Discretization of laplace operator for u
  sparse_matrix K1(nbdof_u, nbdof_u);
  asm_stiffness_matrix_for_homogeneous_laplacian_componentwise
    (K1, mim, mf_u, mpirg);
  gmm::scale(K1, nu);
  
  // Mass Matrix
  sparse_matrix M(nbdof_u, nbdof_u);
  asm_mass_matrix(M, mim, mf_u, mpirg);
  


  // Matrix p div u
  sparse_matrix B(nbdof_p, nbdof_u);
  asm_stokes_B(B, mim, mf_u, mf_p, mpirg);
  
  mf_mult.set_qdim(N);
  dal::bit_vector dofon_Dirichlet = mf_mult.dof_on_region(DIRICHLET_BOUNDARY_NUM);
  dal::bit_vector dofon_nonref =mf_mult.dof_on_region(NONREFLECTIVE_BOUNDARY_NUM);
  dofon_Dirichlet.setminus(dofon_nonref);

  // Normal part Dirichlet condition
  mf_mult.set_qdim(1);
  dal::bit_vector dofon_NDirichlet
    = mf_mult.dof_on_region(NORMAL_PART_DIRICHLET_BOUNDARY_NUM);
  std::vector<size_type> ind_ct_ndir;
  for (dal::bv_visitor i(dofon_NDirichlet); !i.finished(); ++i) {
    if (dofon_Dirichlet.is_in(i*N) || dofon_nonref.is_in(i*N)) 
      dofon_NDirichlet.sup(i);  // Suppress i because it is on the
    else                        // Dirichlet or non reflective boundary.
      ind_ct_ndir.push_back(i);
  }
  size_type nbdof_NDir = dofon_NDirichlet.card();
  gmm::sub_index SUB_CT_NDIR(ind_ct_ndir);
  gmm::sub_interval I2(nbdof_u, nbdof_NDir);
  getfem::mesh_region mpindirrg
    =mf_u.linked_mesh().get_mpi_sub_region(NORMAL_PART_DIRICHLET_BOUNDARY_NUM);

  sparse_matrix HND(nbdof_NDir, nbdof_u);
  {
    sparse_matrix A(mf_mult.nb_dof(), nbdof_u); 
    getfem::generic_assembly assem;
    assem.set("M(#2,#1)+=comp(Base(#2).vBase(#1).Normal())(:,:,i,i);");
    assem.push_mi(mim); assem.push_mf(mf_u); assem.push_mf(mf_mult);
    assem.push_mat(A); assem.assembly(mpindirrg);
    gmm::copy(gmm::sub_matrix(A, SUB_CT_NDIR, I1), HND);
  }
  cout << "Nb of Normal par Dirichlet constraints : " << nbdof_NDir << endl;

  // Dirichlet condition
  mf_mult.set_qdim(N);
  size_type nbdof_Dir = dofon_Dirichlet.card();
  std::vector<size_type> ind_ct_dir;
  for (dal::bv_visitor i(dofon_Dirichlet); !i.finished(); ++i) 
    ind_ct_dir.push_back(i);
  gmm::sub_index SUB_CT_DIR(ind_ct_dir);
  gmm::sub_interval I3(nbdof_u+nbdof_NDir, nbdof_Dir);
  getfem::mesh_region mpidirrg
    = mf_u.linked_mesh().get_mpi_sub_region(DIRICHLET_BOUNDARY_NUM);

  sparse_matrix HD(nbdof_Dir, nbdof_u);
  {
    sparse_matrix A(mf_mult.nb_dof(), nbdof_u); 
    getfem::generic_assembly assem;
    assem.set("M(#2,#1)+=comp(vBase(#2).vBase(#1))(:,i,:,i);");
    assem.push_mi(mim); assem.push_mf(mf_u); assem.push_mf(mf_mult);
    assem.push_mat(A); assem.assembly(mpidirrg);
    gmm::copy(gmm::sub_matrix(A, SUB_CT_DIR, I1), HD);
  }
  cout << "Nb of Dirichlet constraints : " << nbdof_Dir << endl;


  // Non reflective condition
  size_type nbdof_nonref = dofon_nonref.card();
  std::vector<size_type> ind_ct_nonref;
  for (dal::bv_visitor i(dofon_nonref); !i.finished(); ++i) 
    ind_ct_nonref.push_back(i);
  gmm::sub_index SUB_CT_NONREF(ind_ct_nonref);
  gmm::sub_interval I4(nbdof_u+nbdof_NDir+nbdof_Dir, nbdof_nonref);
  getfem::mesh_region mpinonrefrg
    = mf_u.linked_mesh().get_mpi_sub_region(NONREFLECTIVE_BOUNDARY_NUM);

  sparse_matrix HNR(nbdof_nonref, nbdof_u);
  {
    sparse_matrix A(mf_mult.nb_dof(), nbdof_u); 
    getfem::generic_assembly assem;
    assem.set("M(#2,#1)+=comp(vBase(#2).vBase(#1))(:,i,:,i);");
    assem.push_mi(mim); assem.push_mf(mf_u); assem.push_mf(mf_mult);
    assem.push_mat(A); assem.assembly(mpinonrefrg);
    gmm::copy(gmm::sub_matrix(A, SUB_CT_NONREF, I1), HNR);

  
  }
  cout << "Nb on Non reflective condition: " << nbdof_nonref << endl;


  // Discretization of laplace operator for p
  sparse_matrix K2(nbdof_p+1, nbdof_p+1);
  gmm::sub_interval IP(0, nbdof_p);
  asm_stiffness_matrix_for_homogeneous_laplacian(gmm::sub_matrix(K2, IP),
						 mim, mf_p, mpirg);
  K2(nbdof_p, 0) = K2(0, nbdof_p) = 1.0; // set the first pressure dof to 0
  

  // 
  // dynamic problem
  //
  plain_vector DF(nbdof_u), F(nbdof_rhs), USTAR(nbdof_u), USTARbis(nbdof_u);
  plain_vector Phi(nbdof_p);

  pdef->initial_condition_u(*this, Un0);
  pdef->initial_condition_p(*this, Pn0);

  //gmm::clear(Un0); gmm::clear(Pn0);
  gmm::vecsave("Un0", Un0);
  gmm::vecsave("Pn0", Pn0);
  
  do_export(0);
  for (scalar_type t = dt; t <= T; t += dt) {

    //
    // Assembly of the first linear system
    //
    size_type sizelsystem = nbdof_u + nbdof_NDir + nbdof_Dir + nbdof_nonref;
    sparse_matrix A1(sizelsystem, sizelsystem);
   
    sparse_matrix A2(sizelsystem, sizelsystem);  // sera utilisée pour reappliquer la CNR
    
    plain_vector Y(sizelsystem), YY(nbdof_u);
    // laplace operator
    gmm::copy(K1, gmm::sub_matrix(A1, I1));

    // Nonlinear part
    getfem::asm_NS_uuT(gmm::sub_matrix(A1, I1), mim, mf_u, Un0, mpirg);
    // Dynamic part
    gmm::add(gmm::scaled(M, 1./dt), gmm::sub_matrix(A1, I1));
        
     gmm::add(M, gmm::sub_matrix(A2, I1));

    
    gmm::mult(M, gmm::scaled(Un0, 1./dt), gmm::sub_vector(Y, I1));

    //    plain_vector subY1(nbdof_u);
    // gmm::mult(M, gmm::scaled(Un0, 1./dt), subY1);
    /*
    plain_vector subY1(nbdof_u);
    gmm::mult(M, gmm::scaled(Un0, 1./dt), subY1);

    // Volumic source term
    pdef->source_term(*this, t, F);
    getfem::asm_source_term(subY1, mim, mf_u, mf_rhs, F,
			    mpirg);
    */
    // Normal Dirichlet condition
    gmm::copy(HND, gmm::sub_matrix(A1, I2, I1));
    gmm::copy(gmm::transposed(HND), gmm::sub_matrix(A1, I1, I2));
    
    gmm::copy(HND, gmm::sub_matrix(A2, I2, I1));
    gmm::copy(gmm::transposed(HND), gmm::sub_matrix(A2, I1, I2));
    
    // Dirichlet condition
    gmm::copy(HD, gmm::sub_matrix(A1, I3, I1));
    gmm::copy(gmm::transposed(HD), gmm::sub_matrix(A1, I1, I3));
    
    gmm::copy(HD, gmm::sub_matrix(A2, I3, I1));
    gmm::copy(gmm::transposed(HD), gmm::sub_matrix(A2, I1, I3));
    
    gmm::resize(F, N * nbdof_rhs);
    pdef->dirichlet_condition(*this, t, F);
    {
      plain_vector VV(mf_mult.nb_dof());
      getfem::asm_source_term(VV, mim, mf_mult, mf_rhs, F, mpidirrg);
      gmm::copy(gmm::sub_vector(VV, SUB_CT_DIR), gmm::sub_vector(Y, I3));
    }
    // Non reflective condition
    gmm::copy(HNR, gmm::sub_matrix(A1, I4, I1));
    gmm::copy(gmm::transposed(HNR), gmm::sub_matrix(A1, I1, I4));
    
    gmm::copy(HNR, gmm::sub_matrix(A2, I4, I1));
    gmm::copy(gmm::transposed(HNR), gmm::sub_matrix(A2, I1, I4));
    {
      plain_vector VV(mf_mult.nb_dof());
    /*  if (t < 0.2)
	getfem::asm_source_term(VV, mim, mf_mult, mf_rhs, F, mpinonrefrg);
      else {
	*/
	//getfem::asm_source_term(VV, mim, mf_mult, mf_rhs, Un0, mpinonrefrg);
	
      //}
	getfem::generic_assembly assem;
    // construction du terme de droite dans [M]*Unp1=F
    
    // mise en place de Un + Un.N*(dUn/dn).N

    std::stringstream ss;
    ss << "u=data$1(#1); "
      "V(#1)+="
       << -dt << "*"
      " comp(vBase(#1).Normal().vGrad(#1).Normal().vBase(#2))"
      "(l,i,i,m,j,k,j,:,k).u(l).u(m)+"
      "comp(vBase(#1).vBase(#2))(j,i,:,i).u(j)";
 
    assem.set(ss.str());
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_mult);  
    assem.push_data(Un0);
    assem.push_vec(VV);
    assem.assembly(mpinonrefrg);
      gmm::copy(gmm::sub_vector(VV, SUB_CT_NONREF), gmm::sub_vector(Y, I4));
    }

    // pressure part
    gmm::mult(gmm::transposed(B), gmm::scaled(Pn0, -1.0), YY);
    gmm::add(YY, gmm::sub_vector(Y, I1));
    
    //
    // Solving the first linear system
    //
    {
      double rcond;
      plain_vector X(sizelsystem);
      SuperLU_solve(A1, X, Y, rcond);
      // if (noisy) cout << "condition number: " << 1.0/rcond << endl;
      gmm::copy(gmm::sub_vector(X, I1), USTAR);
      
      gmm::HarwellBoeing_IO::write("A1.hb", A1);
      gmm::vecsave("Y", Y);
      gmm::vecsave("USTAR", USTAR); exit(1);
    }

    cout << "U* - Un0 = " << gmm::vect_dist2(USTAR, Un0) << endl;

    //
    // Solving the second linear system
    //
    {
      double rcond;
      plain_vector X(nbdof_p+1), Z(nbdof_p+1);
      gmm::mult(B, USTAR, gmm::sub_vector(Z, IP));
      SuperLU_solve(K2, X, Z, rcond);
      // if (noisy) cout << "condition number: " << 1.0/rcond << endl;
      gmm::copy(gmm::sub_vector(X, IP), Phi);
    }

    gmm::mult(M, USTAR, USTARbis);
    gmm::mult(gmm::transposed(B), gmm::scaled(Phi, -1.), USTARbis, USTARbis);
    gmm::copy(USTARbis, gmm::sub_vector(Y, I1));
    //gmm::copy(M, gmm::sub_matrix(A1, I1));
    
    {
      double rcond;
      plain_vector X(sizelsystem);
      SuperLU_solve(A2, X, Y, rcond);     // A2 contient la matrice de masse et les blocs des CL
      gmm::copy(gmm::sub_vector(X, I1), Un1);

    }
   

    // gmm::add(USTAR, Un1);
    gmm::add(gmm::scaled(Phi, 1./dt), Pn0, Pn1);
    
    pdef->validate_solution(*this, t);
    cout << "Un1 - Un0 = " << gmm::vect_dist2(Un1, Un0) << ", t=" << t << ", T=" << T << endl;
    
    gmm::copy(Un1, Un0); gmm::copy(Pn1, Pn0);
    do_export(t);
  }
}



void navier_stokes_problem::do_export(scalar_type t) {
  if (!export_to_opendx) return;
  if (first_export) {
    mf_u.write_to_file("icare.mf_u", true);
    mf_p.write_to_file("icare.mf_p", true);

    exp.reset(new getfem::dx_export(datafilename + ".dx", false));
    if (N <= 2)
      sl.build(mesh, getfem::slicer_none(),2);
    else
      sl.build(mesh, getfem::slicer_boundary(mesh),2);
    exp->exporting(sl,true);
    exp->exporting_mesh_edges();
    t_export = 0;
    first_export = false;
  }
  if (t >= t_export-dt/20.0) {
    t_export += dt_export;
    
    static int cnt = 0;
    char s[128]; sprintf(s, "icare.U%d", cnt++);
    gmm::vecsave(s, Un0);
    
    exp->write_point_data(mf_u, Un0);
    exp->serie_add_object("velocity");
    cout << "Saving Pressure, |p| = " << gmm::vect_norm2(Pn1) << "\n";
    exp->write_point_data(mf_p, Pn1);
    exp->serie_add_object("pressure");

        
    static int cntp=0;
    char sp[128]; sprintf(sp, "icare.P%d", cntp++);
    gmm::vecsave(sp, Pn0);
    
    if (N == 2) {
      plain_vector DU(mf_rhs.nb_dof() * N * N);
      plain_vector Rot(mf_rhs.nb_dof());
      compute_gradient(mf_u, mf_rhs, Un0, DU);
      for (unsigned i=0; i < mf_rhs.nb_dof(); ++i) {
	Rot[i] = DU[i*N*N + 3] - DU[i*N*N + 2];
	if ((Rot[i]*Rot[i])<=1.5){Rot[i]=0;}
      }
      cout << "Saving Rot, |rot| = " << gmm::vect_norm2(Rot) << "\n";
      exp->write_point_data(mf_rhs, Rot);
      exp->serie_add_object("rot");
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
    navier_stokes_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh");
    p.solve();
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
