/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2005 Yves Renard, Julien Pommier, Houari Khenous.    */
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
 * Dynamic friction in linear elasticity.
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
*/

#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_import.h>
#include <getfem_regular_meshes.h>
#include <getfem_Coulomb_friction.h>
#include <gmm.h>
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

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/*
  structure for the friction problem
*/
struct friction_problem {

  enum {
    DIRICHLET_BOUNDARY, CONTACT_BOUNDARY, PERIODIC_BOUNDARY1,
    PERIODIC_BOUNDARY2 
  };
  getfem::getfem_mesh mesh;  /* the mesh */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the friction solution */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_coef;  /* mesh_fem used to represent pde coefficients  */
  scalar_type lambda, mu;    /* Lamé coefficients.                           */
  scalar_type rho, PG;       /* density, and gravity                         */
  scalar_type friction_coef; /* friction coefficient.                        */

  scalar_type residu;        /* max residu for the iterative solvers         */
  
  int scheme;  /* 0 = theta method, 1 = Newmark, 2 = middle point,           */
               /* 3 = middle point with separated contact forces.            */
  size_type N, noisy, nocontact_mass;
  scalar_type beta, theta, gamma;
  plain_vector Dirichlet_displacement;
  scalar_type T, dt, r;
  scalar_type init_vert_pos, init_vert_speed, hspeed;
  bool dt_adapt, periodic, dxexport;

  std::string datafilename;
  ftool::md_param PARAM;

  void solve(void);
  void init(void);
  friction_problem(void) : mf_u(mesh), mf_rhs(mesh), mf_coef(mesh) {}
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
  bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(MESH_TYPE);
  N = pgt->dim();
  if (meshname.compare(0,5, "splx:")==0) {
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
  }
  else getfem::import_mesh(meshname, mesh);
  mesh.optimize_structure();

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residu = PARAM.real_value("RESIDU"); if (residu == 0.) residu = 1e-10;

  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  rho = PARAM.real_value("RHO", "Density");
  PG = PARAM.real_value("PG", "Gravity constant");
  friction_coef = PARAM.real_value("FRICTION_COEF", "Friction coefficient");

  scheme = PARAM.int_value("SCHEME", "Time scheme");
  theta = PARAM.real_value("THETA", "Parameter for the theta-method"); 
  beta = PARAM.real_value("BETA", "Parameter beta for Newmark");
  gamma = PARAM.real_value("GAMMA", "Parameter gamma for Newmark");

  bool Dirichlet = PARAM.int_value("DIRICHLET","Dirichlet condition or not");
  Dirichlet_displacement.resize(N);
  Dirichlet_displacement[N-1] = PARAM.real_value("DIRICHLET_VAL",
				     "parameter for Dirichlet condition");
  T = PARAM.real_value("T", "from [0,T] the time interval");
  dt = PARAM.real_value("DT", "time step");
  dt_adapt = (PARAM.int_value("DT_ADAPT", "time step adaptation") != 0);
  periodic = (PARAM.int_value("PERIODICITY", "peridiodic condition or not")
	      != 0);
  dxexport = (PARAM.int_value("DX_EXPORT", "Exporting on OpenDX format")
	      != 0);
  nocontact_mass = PARAM.int_value("NOCONTACT_MASS", "Suppress the mass "
				   "of contact nodes");
  r = PARAM.real_value("R", "augmentation parameter");
  noisy = (PARAM.int_value("NOISY", "verbosity of iterative methods") != 0);
  init_vert_pos = PARAM.real_value("INIT_VERT_POS", "initial position");
  init_vert_speed = PARAM.real_value("INIT_VERT_SPEED","initial speed");
  hspeed = PARAM.real_value("FOUNDATION_HSPEED","initial speed");
  mf_u.set_qdim(N);

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mf_u.set_finite_element(mesh.convex_index(), pf_u, ppi);

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  const char *data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name == 0) {
    if (!pf_u->is_lagrange()) {
      DAL_THROW(dal::failure_error, "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    }
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u, ppi);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name), ppi);
  }
  
  mf_coef.set_finite_element(mesh.convex_index(),
			     getfem::classical_fem(pgt,0), ppi);

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
	if (un[N-1] < 0.0001 && (N != 3 || (bgeot::vect_dist2(pt, center)
			   > .99*sqrt(25. + 15*15) && pt[N-1] < 20.1)))
	  mf_u.add_boundary_elt(CONTACT_BOUNDARY, cv, f); 
	if (un[0] > 0.98) mf_u.add_boundary_elt(PERIODIC_BOUNDARY1, cv, f); 
	if (un[0] < -0.98) mf_u.add_boundary_elt(PERIODIC_BOUNDARY2, cv, f); 
	if (un[N-1] > 0.1 && Dirichlet)
	  mf_u.add_boundary_elt(DIRICHLET_BOUNDARY, cv, f);
      }
    }
  }
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

void friction_problem::solve(void) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  N = mesh.dim();
  cout << "Number of dof for u: " << mf_u.nb_dof() << endl;

  // Linearized elasticity brick.
  getfem::mdbrick_isotropic_linearized_elasticity<>
    ELAS(mf_u, mf_coef, lambda, mu, true);

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
    gmm::copy(Dirichlet_displacement,
	      gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  getfem::mdbrick_Dirichlet<> DIRICHLET(VOL_F, mf_rhs, F, DIRICHLET_BOUNDARY);

  // Eventual periodic condition (lagrange element only).
  getfem::mdbrick_abstract<> *PERIODIC;
  if (periodic) {
    dal::bit_vector b1 = mf_u.dof_on_boundary(PERIODIC_BOUNDARY1);
    dal::bit_vector b2 = mf_u.dof_on_boundary(PERIODIC_BOUNDARY2);
    sparse_matrix BP(b1.card(), mf_u.nb_dof());
    size_type k =0;
    for (dal::bv_visitor i(b1); !i.finished(); ++i, ++k)
      if (i % N == 0) {
	for (dal::bv_visitor j(b2); !j.finished(); ++j) 
	  if (j % N == 0) {
	    base_node pt = mf_u.point_of_dof(i) - mf_u.point_of_dof(j);
	    pt[0] = 0.;
	    if (gmm::vect_norm2(pt) < 1E-4) {
	      for (size_type l = 0; l < N; ++l)
		{ BP(k+l, i+l) = 1.; BP(k+l, j+l) = -1.; }
	      break; 
	    }
	  }
      }
    gmm::resize(F, b1.card()); gmm::clear(F);
    PERIODIC = new getfem::mdbrick_constraint<>(DIRICHLET, BP, F);
  }
  else PERIODIC = &DIRICHLET;
  
  // Dynamic brick.
  getfem::mdbrick_dynamic<> DYNAMIC(*PERIODIC, mf_coef, rho);
  if (nocontact_mass) DYNAMIC.no_mass_on_boundary(CONTACT_BOUNDARY);

  // contact condition for Lagrange elements
  dal::bit_vector cn = mf_u.dof_on_boundary(CONTACT_BOUNDARY);
  sparse_matrix BN(cn.card()/N, mf_u.nb_dof());
  sparse_matrix BT((N-1)*cn.card()/N, mf_u.nb_dof());
  plain_vector gap(cn.card()/N);
  size_type j = 0;
  for (dal::bv_visitor i(cn); !i.finished(); ++i)
    if (i % N == 0) {
      BN(j, i+N-1) = -1.;
      gap[j] = mf_u.point_of_dof(i)[N-1];
      for (size_type k = 0; k < N-1; ++k) BT((N-1)*j+k, i+k) = 1.;
      ++j;
    }
  getfem::mdbrick_Coulomb_friction<> FRICTION(DYNAMIC, BN, gap,
					      friction_coef, BT);
  cout << "BT = " << BT << endl;

  cout << "Total number of variables: " << FRICTION.nb_dof() << endl;
  getfem::standard_model_state MS(FRICTION);

  plain_vector WT(mf_u.nb_dof()), DF(mf_u.nb_dof()), HSPEED(mf_u.nb_dof());
  plain_vector U0(mf_u.nb_dof()), V0(mf_u.nb_dof()), A0(mf_u.nb_dof());
  plain_vector U1(mf_u.nb_dof()), V1(mf_u.nb_dof()), A1(mf_u.nb_dof());
  scalar_type a(1), b(1), dt0 = dt, t(0);

  // Initial conditions (U0, V0, M A0 = F)
  gmm::clear(U0); gmm::clear(V0);
  for (size_type i=0; i < mf_u.nb_dof(); ++i)
    if ((i % N) == 0) { 
      U0[i+N-1] = init_vert_pos;
      V0[i+N-1] = init_vert_speed;
      HSPEED[i] = hspeed;
    }
  gmm::iteration iter(residu, 0, 40000);
  if ((scheme == 0 || scheme == 1) && !nocontact_mass)
    gmm::cg(DYNAMIC.mass_matrix(), A0, VOL_F.source_term(),
	    gmm::identity_matrix(), iter);
  iter.set_noisy(noisy);

  scalar_type J0 = 0.5*gmm::vect_sp(ELAS.stiffness_matrix(), U0, U0)
    + 0.5 * gmm::vect_sp(DYNAMIC.mass_matrix(), V0, V0)
    - gmm::vect_sp(VOL_F.source_term(), U0);

  std::auto_ptr<getfem::dx_export> exp;
  getfem::stored_mesh_slice sl;
  if (dxexport) {
    exp.reset(new getfem::dx_export(datafilename + ".dx", true));
    sl.build(mesh, getfem::slicer_boundary(mesh),4);
    exp->exporting(sl,true); exp->exporting_mesh_edges();
    exp->write_point_data(mf_u, U0, "stepinit"); 
    exp->serie_add_object("deformationsteps");
  }
  

  while (t <= T) {

    switch (scheme) { // complementary left hand side and velocity complement
    case 0 :
      gmm::add(U0, gmm::scaled(V0, dt), U1);
      gmm::add(gmm::scaled(A0, dt*dt*theta*(1.-theta)), U1);
      gmm::mult(DYNAMIC.mass_matrix(), U1, DF);
      gmm::add(gmm::scaled(U0, -1.), gmm::scaled(V0, -dt*(1.-theta)), WT);
      gmm::add(gmm::scaled(HSPEED, -theta*dt), WT);
      a = 1.; b = dt*dt*theta*theta;
      break;
    case 1 :
      gmm::add(U0, gmm::scaled(V0, dt), U1);
      gmm::add(gmm::scaled(A0, dt*dt*(1.-beta)*0.5), U1);
      gmm::mult(DYNAMIC.mass_matrix(), U1, DF);
      gmm::add(gmm::scaled(U0, -1.),
	       gmm::scaled(V0, dt*(beta*0.5/gamma -1.)), WT);
      gmm::add(gmm::scaled(A0, dt*dt*0.5*(beta-gamma)/gamma), WT);
      gmm::add(gmm::scaled(HSPEED, -0.5*dt*beta/gamma), WT);
      a = 1.; b = dt*dt*beta*0.5;
      break;
    case 2 :
      gmm::add(U0, gmm::scaled(V0, dt*0.5), U1);
      gmm::mult(DYNAMIC.mass_matrix(), U1, DF);
      gmm::copy(gmm::scaled(U0, -1.), WT);
      gmm::add(gmm::scaled(HSPEED, -dt*0.5), WT);
      a = 1.; b = dt*dt*0.25;
      break;
    case 3 : // for the friction, it should be better to take the average 
      // for the contact forces to define the friction threshold
      gmm::add(U0, gmm::scaled(V0, dt), U1);
      gmm::mult(DYNAMIC.mass_matrix(), U1, DF);
      gmm::add(gmm::scaled(A0, (1.-theta)/theta), DF);
      gmm::copy(gmm::scaled(U0, -1.), WT);
      gmm::add(gmm::scaled(HSPEED, -dt), WT);
      a = 1.; b = dt*dt*0.5; 
      break;
    }

    FRICTION.set_WT(WT); FRICTION.set_r(r); 
    DYNAMIC.set_dynamic_coeff(a, b);
    DYNAMIC.set_DF(DF);
    
    iter.init();
    getfem::standard_solve(MS, FRICTION, iter);
    gmm::copy(ELAS.get_solution(MS), U1); 

//     {
//       plain_vector w(gmm::mat_nrows(BN));
//       gmm::mult(BN, U1, gmm::scaled(gap, -1.), w);
//       cout << "Normal dep : " << w << endl;
//       cout << "Contact pressure : " << FRICTION.get_LN(MS) << endl;
//     }

    switch (scheme) { // computation of U^{n+1}, V^{n+1}, A^{n+1}
    case 0 :
      gmm::add(gmm::scaled(U1, 1./dt), gmm::scaled(U0, -1./dt), V1);
      gmm::add(gmm::scaled(V0, -(1.-theta)), V1);
      gmm::scale(V1, 1./theta);
      gmm::add(gmm::scaled(V1, 1./dt), gmm::scaled(V0, -1./dt), A1);
      gmm::add(gmm::scaled(A0, -(1.-theta)), A1);
      gmm::scale(A1, 1./theta);
      break;
    case 1 :
      gmm::add(gmm::scaled(U1, 2./(beta*dt*dt)),
	       gmm::scaled(U0, -2./(beta*dt*dt)), A1);
      gmm::add(gmm::scaled(V0, -2./(beta*dt)), A1);
      gmm::add(gmm::scaled(A0, -(1. - beta)/beta), A1);
      gmm::add(gmm::scaled(A0, (1.-gamma)*dt), gmm::scaled(A1, gamma*dt), V1);
      gmm::add(V0, V1);
      break;
    case 2 :
      gmm::copy(U1, V1);
      gmm::add(gmm::scaled(V1, 2.), gmm::scaled(U0, -1.), U1);
      gmm::add(gmm::scaled(U1, 2./dt), gmm::scaled(U0, -2./dt), V1);
      gmm::add(gmm::scaled(V0, -1), V1);
      break;
    case 3 :
      gmm::mult(gmm::transposed(BN), FRICTION.get_LN(MS), A1);
      gmm::add(gmm::scaled(U1, 2./dt), gmm::scaled(U0, -2./dt), V1);
      gmm::add(gmm::scaled(V0, -1), V1);
      break;
    }

    scalar_type J1 = 0.5*gmm::vect_sp(ELAS.stiffness_matrix(), U1, U1)
      + 0.5 * gmm::vect_sp(DYNAMIC.mass_matrix(), V1, V1)
      - gmm::vect_sp(VOL_F.source_term(), U1);

    if (dt_adapt && gmm::abs(J0-J1) > 1E-4 && dt > 1E-5) {
      dt /= 2.;
      gmm::clear(MS.state());
      cout << "Trying with dt = " << dt << endl;
    }
    else {
      t += dt;    cout << "t = " << t << " dt = " << dt;
      cout << " total energy : " << J1 << endl;
      dt = std::min(2.*dt, dt0);

      gmm::copy(U1, U0); gmm::copy(V1, V0); gmm::copy(A1, A0); J0 = J1;
       if (dxexport) {
	 exp->write_point_data(mf_u, U0);
	 exp->serie_add_object("deformationsteps");
       }
    }
    
  }
  if (periodic) delete PERIODIC;
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
