// -*- c++ -*- (enables emacs c++ mode)
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2005 Michel Fournié, Julien Pommier,                 */
/*                         Yves Renard, Nicolas Roux.                      */
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

#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_regular_meshes.h>
#include <getfem_modeling.h>
#include <getfem_Navier_Stokes.h>
#include "navier_stokes.h"


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
using bgeot::base_matrix; /* small dense matrices. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1 };

struct problem_definition {
  virtual void choose_boundaries(navier_stokes_problem &p) {
    getfem::outer_faces_of_mesh(p.mesh, mesh.region(DIRICHLET_BOUNDARY_NUM));
  }
  virtual void validate_solution(navier_stokes_problem &p, scalar_type t) {
    plain_vector H,R; dirichlet_condition(p, t, H, R);
    scalar_type err = getfem::asm_L2_dist(p.mim, 
					  p.mf_u, p.Un1,
					  p.mf_rhs, R);
    cout << " L2 error(t=" << t <<") : " << err << "\n";
  }
  virtual void dirichlet_condition(navier_stokes_problem &p,
				   const base_node &x, scalar_type t,
				   base_small_vector &r) = 0;
  virtual void dirichlet_condition_h(navier_stokes_problem &p,
				   const base_node &x, scalar_type t,
				     base_matrix &h) {
    h.resize(p.N, p.N); gmm::copy(gmm::identity_matrix(), h);
  }
  virtual void source_term(navier_stokes_problem &p,
			   const base_node &x, scalar_type t,
			   base_small_vector &F) = 0;
  virtual void initial_pressure(navier_stokes_problem &, const base_node &) {
    return 0.;
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
  void dirichlet_condition_h(navier_stokes_problem &p, scalar_type t, 
			     plain_vector &H) {
    unsigned N = p.mesh.dim();
    gmm::resize(H, N*N*p.mf_rhs.nb_dof());
    base_matrix h;
    for (unsigned i=0; i < p.mf_rhs.nb_dof(); ++i) {
      dirichlet_condition_h(p, p.mf_rhs.point_of_dof(i), t, h);
      gmm::copy(base_vector(h), gmm::sub_vector(H, gmm::sub_interval(i*N*N, N*N)));
    }
  }
  void source_term(navier_stokes_problem &p, scalar_type t, plain_vector &F) {
    unsigned N = p.mesh.dim();
    gmm::resize(F, N*p.mf_rhs.nb_dof());
    for (unsigned i=0; i < p.mf_rhs.nb_dof(); ++i)
      gmm::copy(source_term(p.mf_rhs.point_of_dof(i), t),
		gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  }
  
  void initial_condition_u(navier_stokes_problem &p, plain_vector &U0) {
    plain_vector R, F(p.mf_u.nb_dof());
    dirichlet_condition(p, 0., R);
    /* L2 projection from mf_rhs onto mf_u */
    sparse_matrix M;
    getfem::asm_mass_matrix(M, p.mim, p.mf_u);
    getfem::asm_source_term(F, p.mim, p.mf_u, p.mf_rhs, R);
    gmm::iteration iter(1E-13);
    gmm::cg(M, U0, F, gmm::identity_matrix(), iter);
  }
  void initial_condition_p(navier_stokes_problem &p, plain_vector &P0) {
    plain_vector PP(mf_rhs.nb_dof());

    for (unsigned i=0; i < p.mf_rhs.nb_dof(); ++i)
      PP[i] = initial_pressure(p, p.mf_rhs.point_of_dof(i));

    /* L2 projection from mf_rhs onto mf_p */
    plain_vector F(p.mf_p.nb_dof());
    sparse_matrix M;
    getfem::asm_mass_matrix(M, p.mim, p.mf_p);
    getfem::asm_source_term(F, p.mim, p.mf_p, p.mf_rhs, PP);
    gmm::iteration iter(1E-13);
    gmm::cg(M, P0, F, gmm::identity_matrix(), iter);
  }

};

struct problem_definition_Stokes_analytic : public problem_definition {
  virtual void dirichlet_condition(navier_stokes_problem &,
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
  virtual void dirichlet_condition(navier_stokes_problem &,
				   const base_node &P, scalar_type t,
				   base_small_vector &r) {
    r = base_small_vector(p.N);
    scalar_type x = P[0], y = P[1];
    r[0] =  -cos(x)*sin(y)*exp(-2*t*p.nu);
    r[1] =  +sin(x)*cos(y)*exp(-2*t*p.ny);
  }
  virtual void source_term(navier_stokes_problem &p,
			   const base_node &P, scalar_type t,
			   base_small_vector &F) {
    scalar_type x = P[0], y = P[1];
    F = base_small_vector(p.N);
    F[0] = -exp(-4.*t*nu)*sin(2.*x);
    F[1] = -exp(-4.*t*nu)*sin(2.*y);
  }
  virtual scalar_type initial_pressure(navier_stokes_problem &p,
				       const base_node &P) {
    scalar_type x = P[0], y = P[1], t = 0;
    return -1./4 * (cos(2*x)+cos(2*y)) * exp(-4*p.nu * t);
  }
};

struct problem_rotating_cylinder : public problem_definition {
  scalar_type alpha;
  
  virtual void dirichlet_condition(navier_stokes_problem &p,
				   const base_node &P, scalar_type t,
				   base_small_vector &r) {
    r = base_small_vector(p.N);
    scalar_type x = P[0], y = P[1];
    bool on_cylinder = true;
    if (gmm::abs(x-p.BBmin[0]) < 1e-7) on_cylinder = false;
    if (gmm::abs(y-p.BBmin[1]) < 1e-7) on_cylinder = false;
    if (gmm::abs(x-p.BBmax[0]) < 1e-7) on_cylinder = false;
    if (gmm::abs(y-p.BBmax[1]) < 1e-7) on_cylinder = false;
    if (!on_cylinder) {
      r[0] = 1; r[1] = 0;
    } else {
      r[0] = -2.*alpha*y; /* HYPOTHESIS: cylinder centered at (0,0) */
      r[1] = 2.*alpha*x;
    }
  }

  virtual void dirichlet_condition_h(navier_stokes_problem &p,
				     const base_node &P, scalar_type t,
				     base_matrix &h) {
    h.resize(p.N,p.N);
    scalar_type y = P[1];   
    gmm::copy(gmm::identity_matrix(), h);
    if (gmm::abs(y- BBmin[1]) < 1e-7 || gmm::abs(y - BBmax[1]) < 1e-7)
      h(0,0) = 0;
  }

   virtual void source_term(navier_stokes_problem &p,
			   const base_node &P, scalar_type t,
			   base_small_vector &F) {
    F = base_small_vector(p.N);
  }
  

  problem_rotating_cylinder(scalar_type aa) : a(aa) {}
};

/*
 * structure for the navier_stokes problem
 */
struct navier_stokes_problem {

  getfem::getfem_mesh mesh;  /* the mesh */
  getfem::mesh_im  mim;      /* integration methods.                         */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the velocity              */
  getfem::mesh_fem mf_p;     /* mesh_fem for the pressure                    */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_coef;  /* mesh_fem used to represent pde coefficients  */
  scalar_type Re;            /* Reynolds number */
  scalar_type nu;            /* 1/Re */
  scalar_type dt, T, dt_export;
  unsigned N;
  scalar_type residue;       /* max residue for the iterative solvers        */
  int noisy;
  int export_to_opendx;

  std::auto_ptr<getfem::dx_export> exp;
  getfem::stored_mesh_slice sl;
  bool first_export;
  scalar_type t_export;
  void export();

  enum { METHOD_SPLITTING_STOKES = 0,      /* cf. works of Michel Fournie */
	 METHOD_SPLITTING = 1,             /* cf. works of Michel Fournie */
	 METHOD_FULLY_CONSERVATIVE = 2, 
	 METHOD_PREDICTION_CORRECTION = 3, /* cf. works of Marianna Braza */
  } option;

  std::auto_ptr<problem_definition> pdef;

  std::string datafilename;
  ftool::md_param PARAM;

  plain_vector Un1, Un0, Pn1, Pn0; /* U_{n+1}, U_{n}, P_{n+1} and P_{n} */

  base_small_vector sol_f(const base_small_vector &P, scalar_type t);
  base_small_vector Dir_cond(const base_small_vector &P, scalar_type t);

  bool solve(void);
  void init(void);
  navier_stokes_problem(void) : mim(mesh), mf_u(mesh), mf_p(mesh),
				mf_rhs(mesh), mf_coef(mesh) {}
};

/*
 * Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void navier_stokes_problem::init(void) {
  const char *FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  const char *FEM_TYPE_P  = PARAM.string_value("FEM_TYPE_P","FEM name P");
  const char *INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  std::string meshname
    (PARAM.string_value("MESHNAME", "Nom du fichier de maillage"));

  /* First step : build the mesh */
  if (meshname.compare(0,5, "regular:")==0) {
    const char *MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
    bgeot::pgeometric_trans pgt = 
      bgeot::geometric_trans_descriptor(MESH_TYPE);
    N = pgt->dim();
    std::vector<size_type> nsubdiv(N);
    std::fill(nsubdiv.begin(),nsubdiv.end(),
	      PARAM.int_value("NX", "Nomber of space steps "));
    getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			      PARAM.int_value("MESH_NOISED") != 0);
    bgeot::base_matrix M(N,N);
    for (size_type i=0; i < N; ++i) {
      static const char *t[] = {"LX","LY","LZ"};
      M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
    }
    mesh.transformation(M);
  } else {
    getfem::import_mesh(meshname, mesh);
    N = mesh.dim();
  }
  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residue = PARAM.real_value("RESIDUE"); if (residue == 0.) residue = 1e-10;

  nu = PARAM.real_value("NU", "Viscosity");
  dt = PARAM.real_value("DT", "Time step");
  T = PARAM.real_value("T", "Final time");
  dtexport = PARAM.real_value("DT_EXPORT", "Time step for export");
  noisy = PARAM.int_value("NOISY", "");
  option = PARAM.int_value("OPTION", "option");

  int prob = PARAM.int_value("PROBLEM", "the problem");
  switch (prob) {
    case 1: pdef.reset(new problem_definition_Stokes_analytic); break;
    case 2: pdef.reset(new problem_definition_Green_Taylor_analytic); break;
    case 3: pdef.reset(new problem_rotating_cylinder(PARAM.real_value("CYL_ROT_SPEED"))); break;
    default: DAL_THROW(dal::failure_error, "wrong PROBLEM value");
  }

  export_to_opendx = PARAM.int_value("DX_EXPORT", "");
  first_export = false;

  Re = 1 / nu;
  mf_u.set_qdim(N);

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION); 

  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_p.set_finite_element(mesh.convex_index(),
			  getfem::fem_descriptor(FEM_TYPE_P));

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  const char *data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name == 0) {
    if (!pf_u->is_lagrange()) {
      DAL_THROW(dal::failure_error, "You are using a non-lagrange FEM "
		<< data_fem_name << ". In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    }
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }
  
  /* set the finite element on mf_coef. Here we use a very simple element
   *  since the only function that need to be interpolated on the mesh_fem 
   * is f(x)=1 ... */
  mf_coef.set_classical_finite_element(0);

  /* set boundary conditions */
  cout << "Choosing boundaries\n";
  pdef->choose_boundaries(*this);
}


void solve() {
  cout << "Number of dof for u : " << mf_u.nb_dof() << endl;
  cout << "Number of dof for p : " << mf_p.nb_dof() << endl;
  gmm::resize(Un0, mf_u.nb_dof());
  gmm::resize(Un1, mf_u.nb_dof());
  gmm::resize(Pn0, mf_p.nb_dof());
  gmm::resize(Pn1, mf_p.nb_dof());
  switch (option) {
    case METHOD_SPLITTING_STOKES: solve_METHOD_SPLITTING(true); break;
    case METHOD_SPLITTING: solve_METHOD_SPLITTING(false); break;
    case METHOD_FULLY_CONSERVATIVE: solve_FULLY_CONSERVATIVE(); break;
    case METHOD_PREDICTION_CORRECTION: solve_PREDICTION_CORRECTION(); break;
    default: DAL_THROW(dal::failure_error, "unknown method");
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
  getfem::mdbrick_scalar_elliptic<> velocity(mim, mf_u, mf_coef, nu); previous = &velocity;

  
  std::auto_ptr<getfem::mdbrick_NS_uuT<> > velocity_nonlin;
  if (!stokes_only) {
    velocity_nonlin.reset(new getfem::mdbrick_NS_uuT<>(velocity));
    previous = &velocity_nonlin;
  }
  
  // Volumic source term
  getfem::mdbrick_source_term<> velocity_f(*previous, 
					   mf_rhs, 
					   plain_vector(mf_rhs.nb_dof()*N));

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<> velocity_dir(velocity_f, mf_rhs, DIRICHLET_BOUNDARY_NUM);
  
  // Dynamic brick.
  getfem::mdbrick_dynamic<> velocity_dyn(velocity_dir, mf_coef, 1.);
  velocity_dyn.set_dynamic_coeff(1.0/dt, 1.0);

  // 
  // definition of the second problem
  //

  getfem::mdbrick_mass_matrix<> mixed(mim, mf_u, mf_coef, 1./dt);
  
  // Pressure term
  getfem::mdbrick_linear_incomp<>mixed_p(mixed, mf_p);

  // Condition on the pressure
  sparse_matrix G(1, mf_p.nb_dof());
  G(0,0) = 1.;
  plain_vector gr(1);
  getfem::mdbrick_constraint<> set_pressure(mixed_p, G, gr, 1);
    
  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<> mixed_dir(set_pressure, mf_rhs, DIRICHLET_BOUNDARY_NUM);
  {
    plain_vector H(mf_rhs.nb_dof()*N*N);
    pdef->dirichlet_condition_h(*this, 0, H);
    velocity_dir.set_H(H);
    mixed_dir.set_H(H);
  }

  // Dynamic brick.
  getfem::mdbrick_dynamic<> mixed_dyn(mixed_dir, mf_coef, 1.);
  mixed_dyn.set_dynamic_coeff(0.0, 1.0);

  // 
  // dynamic problem
  //

  plain_vector DF(mf_u.nb_dof()), F(mf_rhs.nb_dof());

  pdef->initial_condition_u(*this, Un0);
  
  gmm::iteration iter(residue, noisy);
  getfem::standard_model_state MSL(velocity_dyn);
  getfem::standard_model_state MSM(mixed_dyn);

  for (scalar_type t = dt; t <= T; t += dt) {

    iter.init();
    if (!stokes_only) velocity_NS_uuT->set_U0(Un0);
    pdef->source_term(*this, t, F);
    velocity_f.set_rhs(F);
    pdef->dirichlet_condition(*this, t, F);
    velocity_dir.set_rhs(F);
    
    gmm::mult(velocity_dyn.mass_matrix(), gmm::scaled(Un0, 1./dt), DF);
    velocity_dyn.set_DF(DF);
    getfem::standard_solve(MSL, velocity_dyn, iter);

    gmm::copy(basic_velocity->get_solution(MSL), Un1);
    cout << " gmm::vect_dist2(Un0, USTAR) = " << gmm::vect_dist2(Un0, Un1) << endl;
      
    iter.init();
    gmm::mult(velocity_dyn.mass_matrix(), gmm::scaled(Un1, 1./dt), DF);
    mixed_dyn.set_DF(DF);
    mixed_dir.set_rhs(F);
    getfem::standard_solve(*MSM, *mixed_dyn, iter);
    gmm::copy(mixed.get_solution(*MSM), Un1);
    gmm::copy(mixed_p.get_pressure(*MSM), Pn1);

    pdef->validate_solution(*this);

    gmm::copy(Un1, Un0);
  }
}



bool navier_stokes_problem::solve_FULLY_CONSERVATIVE() {

  // Velocity brick.  
  getfem::mdbrick_abstract<>  *velocity_nonlin, *velocity;
  getfem::mdbrick_NS_uuT<> *velocity_NS_uuT = 0;
  getfem::mdbrick_scalar_elliptic<> *basic_velocity = 0;
  getfem::mdbrick_navier_stokes<> *global_velocity = 0;
  sparse_matrix B;
  
  getfem::mdbrick_navier_stokes<> velocity(mim, mf_u, mf_p, nu);
  // Condition on the pressure
  sparse_matrix G(1, mf_p.nb_dof());
  G(0,0) = 1.;
  plain_vector gr(1);
  getfem::mdbrick_constraint<> velocity_ctr(velocity, G, gr, 1);
  

  // Volumic source term
  getfem::mdbrick_source_term<> velocity_f(velocity_ctr, 
					   mf_rhs, 
					   plain_vector(mf_rhs.nb_dof()*N));

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<> velocity_dir(velocity_f, mf_rhs, DIRICHLET_BOUNDARY_NUM);
  
  // Dynamic brick.
  getfem::mdbrick_dynamic<> velocity_dyn(velocity_dir, mf_coef, 1.);
  velocity_dyn.set_dynamic_coeff(1.0/dt, 1.0);

  // 
  // dynamic problem
  //

  plain_vector DF(mf_u.nb_dof()), F(mf_rhs.nb_dof());

  pdef->initial_condition_u(*this, Un0);
  
  gmm::iteration iter(residue, noisy);
  getfem::standard_model_state MSL(velocity_dyn);
  
  for (scalar_type t = dt; t <= T; t += dt) {

    iter.init();
  
    pdef->source_term(*this, t, F);
    velocity_f.set_rhs(F);
    pdef->dirichlet_condition(*this, t, F);
    velocity_dir.set_rhs(F);

    gmm::mult(velocity_dyn.mass_matrix(), gmm::scaled(Un0, 1./dt), DF);
    velocity_dyn.set_DF(DF);
    getfem::standard_solve(MSL, velocity_dyn, iter);

    gmm::copy(global_velocity->get_velocity(MSL), Un1);
   
    pdef->validate_solution(*this);
    
    gmm::copy(Un1, Un0);
  }
}


bool navier_stokes_problem::solve_PREDICTION_CORRECTION() {
  // Velocity brick.  
  sparse_matrix B;
  getfem::mdbrick_scalar_elliptic<> velocity(mim, mf_u, mf_coef, nu);
  getfem::mdbrick_NS_uuT<> velocity_nonlin(velocity);

  gmm::resize(B, mf_p.nb_dof(), mf_u.nb_dof());
  asm_stokes_B(B, mim, mf_u, mf_p);

  // Volumic source term
  getfem::mdbrick_source_term<> velocity_f(velocity_nonlin, 
					   mf_rhs, 
					   plain_vector(mf_rhs.nb_dof()*N));

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<> velocity_dir(velocity_f, mf_rhs, DIRICHLET_BOUNDARY_NUM);
  
  // Dynamic brick.
  getfem::mdbrick_dynamic<> velocity_dyn(velocity_dir, mf_coef, 1.);
  velocity_dyn.set_dynamic_coeff(1.0/dt, 1.0);

  // 
  // Poisson problem for prediction correction scheme
  //
  
  getfem::mdbrick_scalar_elliptic<> poisson(mim, mf_p, mf_coef, 1.0);
  sparse_matrix G(1, mf_p.nb_dof()); G(0,0) = 1.;  
  getfem::mdbrick_constraint<> poisson_setonedof(poisson, G, plain_vector(1));
  getfem::mdbrick_source_term<> poisson_source(poisson_setonedof,
					       mf_rhs, plain_vector(mf_rhs.nb_dof()));

  // 
  // dynamic problem
  //
  plain_vector DF(mf_u.nb_dof()), F(mf_rhs.nb_dof()), USTAR(mf_rhs.nb_dof());

  pdef->initial_condition_u(*this, Un0);
  pdef->initial_condition_p(*this, Pn0);

  gmm::iteration iter(residue, noisy);
  getfem::standard_model_state MSL(velocity_dyn);
  getfem::standard_model_state MSM(poisson_source);
  

  for (scalar_type t = dt; t <= T; t += dt) {

    iter.init();

    velocity_NS_uuT.set_U0(Un0);
    pdef->source_term(*this, t, F);
    velocity_f.set_rhs(F);
    pdef->dirichlet_condition(*this, t, F);
    velocity_dir.set_rhs(F);
    
    gmm::mult(velocity_dyn.mass_matrix(), gmm::scaled(Un0, 1./dt), DF);
    gmm::mult_add(gmm::transposed(B), gmm::scaled(Pn0, -1.), DF);
    velocity_dyn.set_DF(DF);
    getfem::standard_solve(MSL, velocity_dyn, iter);

    gmm::copy(basic_velocity->get_solution(MSL), USTAR);
    cout << " gmm::vect_dist2(Un0, USTAR) = " << gmm::vect_dist2(Un0, USTAR) << endl;
    gmm::mult(B, USTAR, Pn1);
    cout << "div ustar = " << gmm::vect_norm2(Pn1) << endl;
    poisson_source->set_auxF(Pn1);
    iter.init();
    getfem::standard_solve(*MSM, *poisson_source, iter);
    gmm::mult(gmm::transposed(B), gmm::scaled(poisson->get_solution(*MSM), -1.), USTARbis);

    gmm::iteration iter2 = iter; iter2.reduce_noisy(); iter2.init();
    gmm::cg(velocity_dyn.mass_matrix(), Un1, USTARbis,
	    gmm::identity_matrix(), iter2);

    gmm::add(USTAR, Un1);
    gmm::add(gmm::scaled(poisson->get_solution(*MSM), 1./dt), Pn0, Pn1);
    
    pdef->validate_solution(*this);
    
    gmm::copy(Un1, Un0); gmm::copy(Pn1, Pn0);
  }
}


void navier_stokes_problem::export() {
  if (!export_to_opendx) return;
  if (first_export) {
    std::auto_ptr<getfem::dx_export> exp;
    getfem::stored_mesh_slice sl;
    exp.reset(new getfem::dx_export(datafilename + ".dx", false));
    if (N <= 2)
      sl.build(mesh, getfem::slicer_none(),4);
    else
      sl.build(mesh, getfem::slicer_boundary(mesh),4);
    exp->exporting(sl,true);
    exp->exporting_mesh_edges();
    exp->write_point_data(mf_u, Un0, "stepinit"); 
    exp->serie_add_object("deformationsteps");
    t_export = 0;
    first_export = false;
  } else if (t >= t_export-dt/20.0) {
    exp->write_point_data(mf_u, Un0);
    exp->serie_add_object("deformationsteps");
    t_export += dtexport;
  }
}

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {
  dal::exception_callback_debug cb;
  dal::exception_callback::set_exception_callback(&cb); // to debug ...

#ifdef GETFEM_HAVE_FEENABLEEXCEPT /* trap SIGFPE */
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif

  try {    
    navier_stokes_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh");
    p.solve();
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0; 
}
