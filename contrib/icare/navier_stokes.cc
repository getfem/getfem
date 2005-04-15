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

/*
 * structure for the navier_stokes problem
 */
struct navier_stokes_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1 };
  getfem::getfem_mesh mesh;  /* the mesh */
  getfem::mesh_im  mim;      /* integration methods.                         */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the velocity              */
  getfem::mesh_fem mf_p;     /* mesh_fem for the pressure                    */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_coef;  /* mesh_fem used to represent pde coefficients  */
  scalar_type nu, dt, T, dtexport;

  scalar_type residu;        /* max residu for the iterative solvers         */
  int noisy, dxexport, option;

  std::string datafilename;
  ftool::md_param PARAM;

  bool solve(void);
  void init(void);
  navier_stokes_problem(void) : mim(mesh), mf_u(mesh), mf_p(mesh),
				mf_rhs(mesh), mf_coef(mesh) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void navier_stokes_problem::init(void) {
  const char *MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  const char *FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  const char *FEM_TYPE_P  = PARAM.string_value("FEM_TYPE_P","FEM name P");
  const char *INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  /* First step : build the mesh */
  bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Number of space steps "));
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			    PARAM.int_value("MESH_NOISED") != 0);
  
  /* scale the unit mesh to [LX,LY,..] and incline it */
   bgeot::base_matrix M(N,N);
  for (size_type i=0; i < N; ++i) {
    static const char *t[] = {"LX","LY","LZ"};
    M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
  }
  mesh.transformation(M);

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residu = PARAM.real_value("RESIDU"); if (residu == 0.) residu = 1e-10;

  nu = PARAM.real_value("NU", "Viscosity");
  dt = PARAM.real_value("DT", "Time step");
  T = PARAM.real_value("T", "Final time");
  dtexport = PARAM.real_value("DT_EXPORT", "Final time");
  noisy = PARAM.int_value("NOISY", "");
  option = PARAM.int_value("OPTION", "option");
  dxexport = PARAM.int_value("DX_EXPORT", "");
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
  mf_coef.set_finite_element(mesh.convex_index(),
			     getfem::classical_fem(pgt,0));

  /* set boundary conditions
   * (Neuman on the upper face, Dirichlet elsewhere) */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::convex_face_ct border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::convex_face_ct::const_iterator it = border_faces.begin();
       it != border_faces.end(); ++it) {
    assert(it->f != size_type(-1));
    base_node un = mesh.normal_of_face_of_convex(it->cv, it->f);
    un /= gmm::vect_norm2(un);
    if (0) {
    // if (gmm::abs(un[N-1] - 1.0) < 1.0E-7) { // new Neumann face
      mesh.add_face_to_set(NEUMANN_BOUNDARY_NUM, it->cv, it->f);
    } else {
      mesh.add_face_to_set(DIRICHLET_BOUNDARY_NUM, it->cv, it->f);
    }
  }
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

scalar_type NU_;

base_small_vector sol_f(const base_small_vector &P, scalar_type t) {
  base_small_vector res(P.size());
  scalar_type x = P[0], y = P[1];
  res[0] = -16.*y*x*x+16.*y*x+8.*x*x-8.*x+32.*NU_*t*y
    -16.*NU_*t+8.*t*x*x-8.*t*x;
  res[1] =  16.*x*y*y-16.*y*x-8.*y*y+8.*y-32.*NU_*t*x
    +16.*NU_*t+8.*t*y*y-8.*t*y;
  // gmm::clear(res);
  return res;
}

base_small_vector Dir_cond(const base_small_vector &P, scalar_type t) {
  base_small_vector res(P.size());
  scalar_type x = P[0], y = P[1];
  res[0] =  2.*(2.*y-1.)*(1.-1.*gmm::sqr(2.*x-1.))*t;
  res[1] = -2.*(2.*x-1.)*(1.-1.*gmm::sqr(2.*y-1.))*t;
  // if (x > .5) res[0] = res[1] = 0;
  return res;
}

bool navier_stokes_problem::solve() {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();

  cout << "Number of dof for u : " << mf_u.nb_dof() << endl;
  cout << "Number of dof for p : " << mf_p.nb_dof() << endl;
  NU_ = nu;
  // 
  // definition of the Laplacian problem
  //

  // Velocity brick.  
  std::auto_ptr<getfem::mdbrick_abstract<> > velocity_nonlin, velocity;
  getfem::mdbrick_NS_uuT<> *velocity_NS_uuT = 0;
  getfem::mdbrick_scalar_elliptic<> *basic_velocity = 0;
  getfem::mdbrick_navier_stokes<> *global_velocity = 0;
  switch (option) {
    case 1 : case 2 :
      basic_velocity
	= new getfem::mdbrick_scalar_elliptic<>(mim, mf_u, mf_coef, nu, true);
      if (option == 2) {
	velocity.reset(basic_velocity);
	velocity_nonlin.reset(velocity_NS_uuT
			      = new getfem::mdbrick_NS_uuT<>(*velocity));
      }
      else velocity_nonlin.reset(basic_velocity);
      break;
    case 3 :
      global_velocity 
	= new getfem::mdbrick_navier_stokes<>(mim, mf_u, mf_p, nu);
      velocity.reset(global_velocity);
      { // Condition on the pressure
	sparse_matrix G(1, mf_p.nb_dof());
	G(0,0) = 1.;
	plain_vector gr(1);
	velocity_nonlin.reset(new getfem::mdbrick_constraint<>
			      (*global_velocity, G, gr, 1));
      }
      break;
    default : DAL_THROW(dal::failure_error, "This option does not exist");
  }

  // Volumic source term
  plain_vector F(nb_dof_rhs * N);
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    gmm::copy(sol_f(mf_rhs.point_of_dof(i), 0.),
	      gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  getfem::mdbrick_source_term<> velocity_f(*velocity_nonlin, mf_rhs, F);


  // Dirichlet condition brick.
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    gmm::copy(Dir_cond(mf_rhs.point_of_dof(i), 0.),
	      gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  getfem::mdbrick_Dirichlet<> velocity_dir(velocity_f, mf_rhs,
					  F, DIRICHLET_BOUNDARY_NUM);
  
  // Dynamic brick.
  getfem::mdbrick_dynamic<> velocity_dyn(velocity_dir, mf_coef, 1.);
  velocity_dyn.set_dynamic_coeff(1.0/dt, 1.0);

  // 
  // definition of the mixed problem
  //

  getfem::mdbrick_mass_matrix<> mixed(mim, mf_u, mf_coef, 1./dt, true);
  
  // Pressure term
  getfem::mdbrick_linear_incomp<> mixed_p(mixed, mf_p);

  // Condition on the pressure
  sparse_matrix G(1, mf_p.nb_dof());
  G(0,0) = 1.;
  plain_vector gr(1);
  getfem::mdbrick_constraint<> set_pressure(mixed_p, G, gr, 1);
  
  // Dirichlet condition brick.
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    gmm::copy(Dir_cond(mf_rhs.point_of_dof(i), 0.),
	      gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  getfem::mdbrick_Dirichlet<> mixed_dir(set_pressure, mf_rhs,
					F, DIRICHLET_BOUNDARY_NUM);

  // Dynamic brick.
  getfem::mdbrick_dynamic<> mixed_dyn(mixed_dir, mf_coef, 1.);
  mixed_dyn.set_dynamic_coeff(0.0, 1.0);

  // 
  // dynamic problem
  //

  plain_vector DF(mf_u.nb_dof()), U0(mf_u.nb_dof()), USTAR(mf_u.nb_dof()),
    USTARbis(mf_u.nb_dof());
  
  gmm::iteration iter(residu, noisy);
  getfem::standard_model_state MSL(velocity_dyn);
  getfem::standard_model_state MSM(mixed_dyn);
  
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
  }

  scalar_type t_export(dtexport);
  for (scalar_type t = dt; t <= T; t += dt) {

    iter.init();
    if (option == 2) {
      velocity_NS_uuT->set_U0(U0);
    }
    for (size_type i = 0; i < nb_dof_rhs; ++i)
      gmm::copy(sol_f(mf_rhs.point_of_dof(i), t),
		gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
    velocity_f.set_rhs(F);
    for (size_type i = 0; i < nb_dof_rhs; ++i)
      gmm::copy(Dir_cond(mf_rhs.point_of_dof(i), t),
		gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
    velocity_dir.set_rhs(F);

    
    gmm::mult(velocity_dyn.mass_matrix(), gmm::scaled(U0, 1./dt), DF);
    velocity_dyn.set_DF(DF);
    getfem::standard_solve(MSL, velocity_dyn, iter);

    if (option == 1 || option == 2) {
      gmm::copy(basic_velocity->get_solution(MSL), USTAR);
    
      iter.init();
      gmm::mult(velocity_dyn.mass_matrix(), gmm::scaled(USTAR, 1./dt), DF);
      mixed_dyn.set_DF(DF);
      for (size_type i = 0; i < nb_dof_rhs; ++i)
	gmm::copy(Dir_cond(mf_rhs.point_of_dof(i), t),
		gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
      mixed_dir.set_rhs(F);
      getfem::standard_solve(MSM, mixed_dyn, iter);
      gmm::copy(mixed.get_solution(MSM), U0);
    }
    if (option == 3) {
      gmm::copy(global_velocity->get_velocity(MSL), U0);
    }
    
    cout << "Kinetic energy : " <<
      0.5 * gmm::vect_sp(velocity_dyn.mass_matrix(), U0, U0) << endl;
    
    cout << "error = " << gmm::vect_dist2(U0, F) << endl;

    if (dxexport && t >= t_export-dt/20.0) {
      exp->write_point_data(mf_u, U0);
      exp->serie_add_object("deformationsteps");
      t_export += dtexport;
    }

  }
  return 1;
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
    if (!p.solve()) DAL_THROW(dal::failure_error,"Solve has failed");
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0; 
}
