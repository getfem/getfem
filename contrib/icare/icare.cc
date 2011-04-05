// -*- c++ -*- (enables emacs c++ mode)
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2009 Michel FourniŽ, Julien Pommier,                 */
/*                         Yves Renard, Nicolas Renon, Nicolas Roux.       */
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
/* You  should  have received a copy of the GNU Lesser General Public      */
/* License along  with  this program;  if not, write to the Free Software  */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */

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

#if GETFEM_PARA_LEVEL > 1
  #include "/usr/include/mpi.h"
#endif

//#include <iostream>
//#include "gmm/gmm_MUMPS_interface.h"


/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
//using bgeot::dim_type;
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrices. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

enum {
  DIRICHLET_BOUNDARY_NUM = 0, 
  NONREFLECTIVE_BOUNDARY_NUM,
  NORMAL_PART_DIRICHLET_BOUNDARY_NUM, 
  ON_CYLINDER_BOUNDARY_NUM,  
  NEUMANN_BOUNDARY_NUM};

#if GETFEM_PARA_LEVEL > 1
#ifdef GMM_USES_MPI
  template <typename VECT> inline void MPI_SUM_VECTOR2(VECT &V) {
    typedef typename gmm::linalg_traits<VECT>::value_type T;
    std::vector<T> W(gmm::vect_size(V));
//    cout<<"MPI ALLREDUCE"<<endl;
    MPI_Allreduce(&(V[0]), &(W[0]), gmm::vect_size(V), gmm::mpi_type(T()),
		  MPI_SUM, MPI_COMM_WORLD);
 gmm::copy(W, V);
  }
#endif
#endif

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
  getfem::mesh_fem mf_mult;  /* mesh_fem for Dirichlet condition ...            */
  scalar_type Re;            /* Reynolds number */
  scalar_type nu;            /* 1/Re */
  //  scalar_type R;              /*Radius of the cylinder */
  scalar_type dt, T, Tinitial, dt_export;
  unsigned N;
  scalar_type residual;      /* max residual for the iterative solvers       */
  int noisy;
  int export_to_opendx;
  int non_reflective_bc;

  std::auto_ptr<getfem::dx_export> exp;
  getfem::stored_mesh_slice sl;
  bool first_export;
  scalar_type t_export;
  void do_export(scalar_type t);

  int option, time_order;

  std::auto_ptr<problem_definition> pdef;

  std::string datafilename;
  ftool::md_param PARAM;

  plain_vector Un1, Un0, Pn1, Pn0, lambda, tmp; /* U_{n+1}, U_{n}, P_{n+1} and P_{n} */
  plain_vector Unm1; /*U_{n-1} */

  base_small_vector sol_f(const base_small_vector &P, scalar_type t);
  base_small_vector Dir_cond(const base_small_vector &P, scalar_type t);

  void solve(void);
  //void solve_METHOD_SPLITTING(bool);
  //void solve_FULLY_CONSERVATIVE();
  //void solve_PREDICTION_CORRECTION();
  void solve_PREDICTION_CORRECTION2();
  //void solve_PREDICTION_CORRECTION_ORDER2();
  //  void solve_FULLY_EXPLICIT_ORDER2();

  void init(void);
  navier_stokes_problem(void) : mim(mesh), mf_u(mesh), mf_p(mesh),
				mf_rhs(mesh), mf_mult(mesh)  {}
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
    return 0.5;
  }
  virtual base_small_vector initial_velocity(navier_stokes_problem &p, const base_node &P) {
    base_small_vector r; dirichlet_condition(p,P,0,r); return r;
  }
  void dirichlet_condition(navier_stokes_problem &p, scalar_type t, 
			   plain_vector &R) {
    unsigned N = p.mesh.dim();
    gmm::resize(R, N*p.mf_rhs.nb_dof());
    base_small_vector r; 
    GMM_ASSERT1(!p.mf_rhs.is_reduced(), "To be adapted");
    for (unsigned i=0; i < p.mf_rhs.nb_dof(); ++i) {
      dirichlet_condition(p, p.mf_rhs.point_of_basic_dof(i), t, r);
      gmm::copy(r, gmm::sub_vector(R, gmm::sub_interval(i*N, N)));
    }
  }
  void source_term(navier_stokes_problem &p, scalar_type t, plain_vector &F) {
    unsigned N = p.mesh.dim();
    gmm::resize(F, N*p.mf_rhs.nb_dof());
    GMM_ASSERT1(!p.mf_rhs.is_reduced(), "To be adapted");
    base_small_vector f;
    for (unsigned i=0; i < p.mf_rhs.nb_dof(); ++i) {
      source_term(p, p.mf_rhs.point_of_basic_dof(i), t, f);
      gmm::copy(f, gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
    }
  }
  
  void initial_condition_u(navier_stokes_problem &p, plain_vector &U0) {
    GMM_ASSERT1(!p.mf_rhs.is_reduced(), "To be adapted");
    plain_vector R(p.N*p.mf_rhs.nb_dof());
    for (unsigned i=0; i < p.mf_rhs.nb_dof(); ++i) {
      base_small_vector
	r = initial_velocity(p, p.mf_rhs.point_of_basic_dof(i));
      gmm::copy(r, gmm::sub_vector(R, gmm::sub_interval(i*p.N, p.N)));
    }
    gmm::copy(R,U0);
  }
//  void initial_condition_u(navier_stokes_problem &p, plain_vector &U0) {
//     GMM_ASSERT1(!p.mf_rhs.is_reduced(), "To be adapted");
//     plain_vector R(p.N*p.mf_rhs.nb_dof()), F(p.mf_u.nb_dof());
//     for (unsigned i=0; i < p.mf_rhs.nb_dof(); ++i) {
//       base_small_vector
// 	r = initial_velocity(p, p.mf_rhs.point_of_basic_dof(i));
//       gmm::copy(r, gmm::sub_vector(R, gmm::sub_interval(i*p.N, p.N)));
//     }
//     /* L2 projection from mf_rhs onto mf_u (we cannot interpolate directly
//        onto mf_u since it can be non-lagrangian) */
//     sparse_matrix M(p.mf_u.nb_dof(), p.mf_u.nb_dof());
//     getfem::asm_mass_matrix(M, p.mim, p.mf_u);
//     getfem::asm_source_term(F, p.mim, p.mf_u, p.mf_rhs, R);
//     gmm::iteration iter(1E-13);
//     gmm::cg(M, U0, F, gmm::identity_matrix(), iter);
//   }

  void initial_condition_p(navier_stokes_problem &p, plain_vector &P0) {
    GMM_ASSERT1(!p.mf_p.is_reduced(), "To be adapted");
    for (unsigned i=0; i < p.mf_p.nb_dof(); ++i)
      P0[i] = initial_pressure(p, p.mf_p.point_of_basic_dof(i));
  }

//   void initial_condition_p(navier_stokes_problem &p, plain_vector &P0) {
//     plain_vector PP(p.mf_rhs.nb_dof());
//     GMM_ASSERT1(!p.mf_rhs.is_reduced(), "To be adapted");
//     for (unsigned i=0; i < p.mf_rhs.nb_dof(); ++i)
//       PP[i] = initial_pressure(p, p.mf_rhs.point_of_basic_dof(i));

//     /* L2 projection from mf_rhs onto mf_p */
//     plain_vector F(p.mf_p.nb_dof());
//     sparse_matrix M(p.mf_p.nb_dof(), p.mf_p.nb_dof());
//     getfem::asm_mass_matrix(M, p.mim, p.mf_p);
//     getfem::asm_source_term(F, p.mim, p.mf_p, p.mf_rhs, PP);
//     gmm::iteration iter(1E-13);
//     gmm::cg(M, P0, F, gmm::identity_matrix(), iter);
//   }
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
			   const base_node &/* P */, scalar_type /* t */,
			   base_small_vector &F) {
    //    scalar_type x = P[0], y = P[1];
    F = base_small_vector(p.N);
    //F[0] = -exp(-4.*t*p.nu)*sin(2.*x);
    //F[1] = -exp(-4.*t*p.nu)*sin(2.*y);
    F[0] =0;
    F[1] =0;
  }
  virtual scalar_type initial_pressure(navier_stokes_problem &p,
				       const base_node &P) {
    scalar_type x = P[0], y = P[1], t = 0;
    return -1./4 * (cos(2*x)+cos(2*y)) * exp(-4*p.nu * t);
  }
};

struct problem_rotating_cylinder : public problem_definition {
  scalar_type alpha;
  
  virtual bool is_on_west_face(base_node BBmin,base_node G){
    bool res = false;
    if (gmm::abs(G[0] - BBmin[0]) < 1e-7) 
      res = true;
    return res;
  }

  virtual bool is_on_est_face(base_node BBmax, base_node G){
    bool res = false;
    if (gmm::abs(G[0] - BBmax[0]) < 1e-7) 
      res = true;
    return res;
  }

 virtual bool is_on_nord_face(base_node BBmax, base_node G){
    bool res = false;
    if (gmm::abs(G[1] - BBmax[1]) < 1e-7) 
      res = true;
    //if (is_on_west_face(BBmin,G))
    //		   res = false;
    //if (is_on_est_face(BBmax,G))
    //		       res = false;
    return res;
 }

 virtual bool is_on_sud_face(base_node BBmin, base_node G){
    bool res = false;
    if (gmm::abs(G[1] - BBmin[1]) < 1e-7) 
      res = true;
    //if (is_on_west_face(BBmin,G))
    //		   res = false;
    //if (is_on_est_face(BBmax,G))
    //		    res = false;
    return res;
  }

 virtual bool is_on_down_face(base_node BBmin, base_node G){
    bool res = false;
    if (gmm::abs(G[2] - BBmin[2]) < 1e-7) 
      res = true;
    // if (is_on_est_face(BBmax,G)) 
    //		   res = false;
    //if (is_on_west_face(BBmin,G)) 
    //		   res = false;
    //if (is_on_nord_face(BBmax,G)) 
    //			res = false;
    //if (is_on_sud_face(BBmin,G)) 
    //		       res = false;
    return res;
  }

 virtual bool is_on_up_face(base_node BBmax, base_node G){
    bool res = false;
    if (gmm::abs(G[2] - BBmax[2]) < 1e-7) 
      res = true;
    //if (is_on_est_face(BBmax,G)) 
    //		       res = false;
    //if (is_on_west_face(BBmin,G)) 
    //			res = false;
    //if (is_on_nord_face(BBmax,G)) 
    //			res = false;
    //if (is_on_sud_face(BBmin,G)) 
    //		       res = false;
    return res;
 }



  virtual void choose_boundaries(navier_stokes_problem &p) {
    getfem::mesh_region r; 
    getfem::outer_faces_of_mesh(p.mesh, r);
    unsigned N = p.mesh.dim();

    if (N==2){
      for (getfem::mr_visitor i(r); !i.finished(); ++i) {
	base_node G = gmm ::mean_value(p.mesh.points_of_face_of_convex(i.cv(),i.f()));
   	//cout << "x=" << G[0] << "     y="<< G[1] << "\n";
	bool on_cyl = true;

       	if  (is_on_west_face(p.BBmin,G) )
	  {on_cyl = false;
	    p.mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());};
	
	//if  (is_on_west_face(p.BBmin,G)||is_on_nord_face(p.BBmax,G) || is_on_sud_face(p.BBmin,G) )
	//  {on_cyl = false;
	//    p.mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());};
	
    	if (is_on_est_face(p.BBmax,G))
	  {on_cyl = false;
	    p.mesh.region(NONREFLECTIVE_BOUNDARY_NUM).add(i.cv(),i.f());};

	if (is_on_nord_face(p.BBmax,G) || is_on_sud_face(p.BBmin,G))
	  {on_cyl = false;
	    p.mesh.region(NORMAL_PART_DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());};

	//scalar_type h = p.mesh.convex_radius_estimate(i.cv());
	//if (gmm::abs(sqrt(G[0]*G[0] + G[1]*G[1]) - p.R) < h)
    	//  {p.mesh.region(ON_CYLINDER_BOUNDARY_NUM).add(i.cv(),i.f());};	
	if (on_cyl)
	  {p.mesh.region(ON_CYLINDER_BOUNDARY_NUM).add(i.cv(),i.f());};
      }
    }
    
    if (N==3){
      for (getfem::mr_visitor i(r); !i.finished(); ++i) {
	base_node G = gmm ::mean_value(p.mesh.points_of_face_of_convex(i.cv(),i.f()));
   	//cout << "x=" << G[0] << "     y="<< G[1]<< "     z="<< G[2] << "\n";
	bool on_cyl = true;

	if  (is_on_west_face(p.BBmin,G) ) 
	  {on_cyl = false;
	    p.mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());};
	
	//if  (is_on_west_face(p.BBmin,G) ||  is_on_nord_face(p.BBmax,G) || is_on_sud_face(p.BBmin,G) ) 
	//  {on_cyl = false;
	//    p.mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());};
	
	if (is_on_est_face(p.BBmax,G))
	  { on_cyl = false;
	    p.mesh.region(NONREFLECTIVE_BOUNDARY_NUM).add(i.cv(),i.f());};
	
	if (is_on_nord_face(p.BBmax,G) || is_on_sud_face(p.BBmin,G))
	  { on_cyl = false;
	    p.mesh.region(NORMAL_PART_DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());};
	
	if  (is_on_down_face(p.BBmin,G) || is_on_up_face(p.BBmax,G))
	  { on_cyl = false;
	    p.mesh.region(NEUMANN_BOUNDARY_NUM).add(i.cv(),i.f());};
	
	//scalar_type h = p.mesh.convex_radius_estimate(i.cv());
	if (gmm::abs(sqrt(G[0]*G[0] + G[1]*G[1])) < 0.5) 
	  {p.mesh.region(ON_CYLINDER_BOUNDARY_NUM).add(i.cv(),i.f());};
	
	//if (on_cyl)
	//  {p.mesh.region(ON_CYLINDER_BOUNDARY_NUM).add(i.cv(),i.f());};
      }
    }
  }
  virtual void dirichlet_condition(navier_stokes_problem &p,
				   const base_node &P, scalar_type /*t*/,
				   base_small_vector &r) {
    r = base_small_vector(p.N); 

    if (p.N==2){
      scalar_type x = P[0], y = P[1];
      //  if (gmm::abs(x - p.BBmax[0]) < 1e-7)
      //       {}
      //       else if ((gmm::abs(y - p.BBmax[1]) < 1e-7
      // 		|| gmm::abs(y - p.BBmin[1]) < 1e-7 ) && !(gmm::abs(x - p.BBmin[0]) < 1e-7) ){}
      //       else{ 
      //         r[0] = 1.0;
      //         r[1] = 0.0;
      //         if (gmm::sqr(x*x + y*y) < 1.2){
      //     	    r[0] = -2.*alpha*y; /* HYPOTHESIS: cylinder centered at (0,0) */
      // 	    r[1] = 2.*alpha*x;
      // 	}
      //       }
      r[0] = 1.0;
      r[1] = 0.0;
      if (gmm::sqrt(x*x + y*y) < 1.0 ){
	r[0] = -2.*alpha*y; /* HYPOTHESIS: cylinder centered at (0,0) */
	r[1] = 2.*alpha*x;
      };
      if (gmm::abs(x - p.BBmin[0]) < 1e-7){
	r[0] = 1.0;
	r[1] = 0.0;
      };
    }
    
    if (p.N==3){
      scalar_type x = P[0], y = P[1]; //, z = P[2];
      //   if (gmm::abs(x - p.BBmax[0]) < 1e-7)
      // 	{}
      //       e2lse if ((gmm::abs(y - p.BBmax[1]) < 1e-7
      // 		|| gmm::abs(y - p.BBmin[1]) < 1e-7 ) && !(gmm::abs(x - p.BBmin[0]) < 1e-7) ){}
      //       else if ((gmm::abs(z - p.BBmax[2] < 1e-7) 
      //                 || gmm::abs(z - p.BBmin[2] < 1e-7) ) && 
      //                 !(gmm::abs(y - p.BBmax[1]) < 1e-7 || gmm::abs(y - p.BBmin[1]) < 1e-7 ) && 
      //                 !(gmm::abs(x - p.BBmin[0]) < 1e-7)){}
      //       else{ 
      //         r[0] = 1.0;
      //         r[1] = 0.0;
      //         r[2] = 0.0;
      //         if (gmm::sqr(x*x + y*y) < 1.2){
      // 	  r[0] = -2.*alpha*y; 
      // 	  r[1] = 2.*alpha*x;
      // 	  r[2] = 0.0;
      // 	}
      //       }
      r[0] = 1.0; 
      r[1] = 0.0;
      r[2] = 0.0;
      if (gmm::sqrt(x*x + y*y) < 1.0 ){
	r[0] = -2.0*alpha*y; /* HYPOTHESIS: cylinder centered at (0,0) */
	r[1] =  2.0*alpha*x;
	r[2] =  0.0;
      };
      if  (gmm::abs(x - p.BBmin[0]) < 1e-7){
	r[0] = 1.0;
	r[1] = 0.0;
	r[2] = 0.0;
      };
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
  virtual base_small_vector initial_velocity(navier_stokes_problem &p,const base_node &) { 
    //unsigned N = p.mesh.dim();
    unsigned N = p.N;

    base_small_vector r(N);
    switch(N) {
    case 1 : {
      r[0] = 1; 
      break;}
    case 2 : {
      r[0] = 1.0; 
      r[1] = 0.0;
      break;}
    case 3 : {
      r[0] = 1.0;
      r[1] = 0.0;
      r[2] = 0.0;
      break;}
    }
    return r;
  }


  virtual scalar_type initial_pressure(navier_stokes_problem &,
				       const base_node &) {
    return 0.0; // the mean value = 0 
  }

  problem_rotating_cylinder(scalar_type aa) : alpha(aa) {}
};


/*
 * Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void navier_stokes_problem::init(void) {
  cout << "-----------------------------------------------------------" << endl;
  cout << "Reading the icare.param file" << endl;
  cout << "-----------------------------------------------------------" << endl;

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
  T  = PARAM.real_value("T", "Final time");
  Tinitial = PARAM.real_value("Tinitial","Initial Time");
  dt_export = PARAM.real_value("DT_EXPORT", "Time step for export");
  noisy = PARAM.int_value("NOISY", "");
  time_order = PARAM.int_value("TIME_ORDER", "Discretization time order");
  option = PARAM.int_value("OPTION", "option");

  //  R = PARAM.real_value("RADIUS","Radius of the cylinder");

  int prob = PARAM.int_value("PROBLEM", "the problem");
  switch (prob) {
    case 1: pdef.reset(new problem_definition_Stokes_analytic); break;
    case 2: pdef.reset(new problem_definition_Green_Taylor_analytic); break;
    case 3: pdef.reset(new problem_rotating_cylinder(PARAM.real_value("CYL_ROT_SPEED"))); break;
  default: GMM_ASSERT1(false, "wrong PROBLEM value");
  }

  non_reflective_bc = PARAM.int_value("NON_REFLECTIVE_BC", "the type of non-reflective boundary condition");
  GMM_ASSERT1(non_reflective_bc >= 0 && non_reflective_bc <= 2, "arg wrong bc");

  export_to_opendx = PARAM.int_value("DX_EXPORT", "");
  first_export = true;

  Re = 1 / nu;
  mf_u.set_qdim(N);

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  getfem::pfem pf_p = getfem::fem_descriptor(FEM_TYPE_P);

  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION); 

  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_p.set_finite_element(mesh.convex_index(), pf_p);

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
  cout << "-----------------------------------------------------------" << endl;
  cout << "Choosing boundaries\n";
  cout << "-----------------------------------------------------------" << endl;

  pdef->choose_boundaries(*this);
}


void navier_stokes_problem::solve() {
  cout << "Number of dof for u : " << mf_u.nb_dof() << endl;
  cout << "Number of dof for p : " << mf_p.nb_dof() << endl;
  gmm::resize(Un0, mf_u.nb_dof());
  gmm::resize(Un1, mf_u.nb_dof());
  gmm::resize(Pn0, mf_p.nb_dof());
  gmm::resize(Pn1, mf_p.nb_dof());
 
  if (time_order==2){
    gmm::resize(Unm1, mf_u.nb_dof());
  }

  switch (option) {
  //case 0 : solve_METHOD_SPLITTING(true); break;
  //case 1 : solve_METHOD_SPLITTING(false); break;
  //case 2 : solve_FULLY_CONSERVATIVE(); break;
  //case 3 : solve_PREDICTION_CORRECTION(); break;
  case 4 : solve_PREDICTION_CORRECTION2(); break;
    //case 5 : solve_PREDICTION_CORRECTION_ORDER2(); break;
    //case 6 : solve_FULLY_EXPLICIT_ORDER2();break;
  default: GMM_ASSERT1(false, "unknown method");
  }
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

// void navier_stokes_problem::solve_METHOD_SPLITTING(bool stokes_only) {
//   // in the splitting method, there are TWO problems:

//   //
//   // definition of the first problem 
//   //

//   // Velocity brick.
//   getfem::mdbrick_abstract<> *previous;
//   getfem::mdbrick_generic_elliptic<> velocity(mim, mf_u, nu);
//   previous = &velocity;

  
//   std::auto_ptr<getfem::mdbrick_NS_uuT<> > velocity_nonlin;
//   if (!stokes_only) {
//     velocity_nonlin.reset(new getfem::mdbrick_NS_uuT<>(velocity));
//     previous = velocity_nonlin.get();
//   }
  
//   // Volumic source term
//   getfem::mdbrick_source_term<> velocity_f(*previous, 
// 					   mf_rhs, 
// 					   plain_vector(mf_rhs.nb_dof()*N));

//   // Dirichlet condition brick.
//   getfem::mdbrick_Dirichlet<> velocity_dir(velocity_f, DIRICHLET_BOUNDARY_NUM);
 
//   // Normal part Dirichlet condition brick.
//   getfem::mdbrick_normal_component_Dirichlet<>
//     velocity_dirnp(velocity_dir, NORMAL_PART_DIRICHLET_BOUNDARY_NUM, mf_rhs);

//   // Non-reflective condition brick
//   getfem::mdbrick_NS_nonref1<> nonreflective(velocity_dirnp,
// 					     NONREFLECTIVE_BOUNDARY_NUM, dt);
 
//   // Dynamic brick.
//   getfem::mdbrick_dynamic<> velocity_dyn(nonreflective, 1.);
//   velocity_dyn.set_dynamic_coeff(1.0/dt, 1.0);

//   // 
//   // definition of the second problem
//   //

//   getfem::mdbrick_mass_matrix<> mixed(mim, mf_u, 1./dt);
  
//   // Pressure term
//   getfem::mdbrick_linear_incomp<> mixed_p(mixed, mf_p);

//   // Condition on the pressure
//   sparse_matrix G(1, mf_p.nb_dof());
//   G(0,0) = 1.;
//   plain_vector gr(1);
//   getfem::mdbrick_constraint<> set_pressure(mixed_p, 1);
//   set_pressure.set_constraints(G, gr);
//   set_pressure.set_constraints_type(getfem::AUGMENTED_CONSTRAINTS);
    
//   // Dirichlet condition brick.
//   getfem::mdbrick_Dirichlet<> mixed_dir(set_pressure, DIRICHLET_BOUNDARY_NUM);
  
//   // Normal part Dirichlet condition brick.
//   getfem::mdbrick_normal_component_Dirichlet<>
//     mixed_dirnp(mixed_dir, NORMAL_PART_DIRICHLET_BOUNDARY_NUM, mf_rhs);

//   // Dynamic brick.
//   getfem::mdbrick_dynamic<> mixed_dyn(mixed_dirnp, 1.);
//   mixed_dyn.set_dynamic_coeff(0.0, 1.0);

//   // 
//   // dynamic problem
//   //

//   plain_vector DF(mf_u.nb_dof()), F(mf_rhs.nb_dof());

//   pdef->initial_condition_u(*this, Un0);
  
//   gmm::iteration iter(residual, noisy);
//   getfem::standard_model_state MSL(velocity_dyn), MSM(mixed_dyn);
  
//   do_export(0);
  
//   for (scalar_type t = dt; t <= T; t += dt) {

//     if (!stokes_only) velocity_nonlin->set_U0(Un0);
//     pdef->source_term(*this, t, F);
//     velocity_f.source_term().set(F);
//     pdef->dirichlet_condition(*this, t, F);
//     velocity_dir.rhs().set(mf_rhs, F);
//     nonreflective.set_Un(Un0);
    
//     gmm::mult(velocity_dyn.get_M(), gmm::scaled(Un0, 1./dt), DF);
//     velocity_dyn.set_DF(DF);
//     iter.init();
//     getfem::standard_solve(MSL, velocity_dyn, iter);

//     gmm::copy(velocity.get_solution(MSL), Un1);
//     gmm::mult(velocity_dyn.get_M(), gmm::scaled(Un1, 1./dt), DF);
//     mixed_dyn.set_DF(DF);
//     mixed_dir.rhs().set(mf_rhs, F);
//     iter.init();
//     getfem::standard_solve(MSM, mixed_dyn, iter);
//     gmm::copy(mixed.get_solution(MSM), Un1);
//     gmm::copy(mixed_p.get_pressure(MSM), Pn1);

//     pdef->validate_solution(*this, t);

//     gmm::copy(Un1, Un0); gmm::copy(Pn1, Pn0);
//     do_export(t);
//   }
// }



// void navier_stokes_problem::solve_FULLY_CONSERVATIVE() {

//   // Velocity brick.    
//   getfem::mdbrick_navier_stokes<> velocity(mim, mf_u, mf_p, nu);
//   // Condition on the pressure
//   sparse_matrix G(1, mf_p.nb_dof());
//   G(0,0) = 1.;
//   plain_vector gr(1);
//   getfem::mdbrick_constraint<> velocity_ctr(velocity);
//   velocity_ctr.set_constraints(G, gr);
//   velocity_ctr.set_constraints_type(getfem::AUGMENTED_CONSTRAINTS);


//   // Volumic source term
//   getfem::mdbrick_source_term<> velocity_f(velocity_ctr, 
// 					   mf_rhs, 
// 					   plain_vector(mf_rhs.nb_dof()*N));

//   // Dirichlet condition brick.
//   getfem::mdbrick_Dirichlet<> velocity_dir(velocity_f, DIRICHLET_BOUNDARY_NUM);
  
//   // Normal part Dirichlet condition brick.
//   getfem::mdbrick_normal_component_Dirichlet<>
//     velocity_dirnp(velocity_dir, NORMAL_PART_DIRICHLET_BOUNDARY_NUM, mf_rhs);

//   // Non-reflective condition brick
//   getfem::mdbrick_NS_nonref1<>
//     nonreflective(velocity_dirnp, NONREFLECTIVE_BOUNDARY_NUM, dt);

//   // Dynamic brick.
//   getfem::mdbrick_dynamic<> velocity_dyn(nonreflective, 1.);
//   velocity_dyn.set_dynamic_coeff(1.0/dt, 1.0);

//   // 
//   // dynamic problem
//   //

//   plain_vector DF(mf_u.nb_dof()), F(mf_rhs.nb_dof());

//   pdef->initial_condition_u(*this, Un0);
  
//   gmm::iteration iter(residual, noisy);
//   getfem::standard_model_state MSL(velocity_dyn);
  
//   do_export(0);
//   for (scalar_type t = dt; t <= T; t += dt) {

  
//     pdef->source_term(*this, t, F);
//     velocity_f.source_term().set(F);
//     pdef->dirichlet_condition(*this, t, F);
//     velocity_dir.rhs().set(mf_rhs, F);
//     nonreflective.set_Un(Un0);

//     gmm::mult(velocity_dyn.get_M(), gmm::scaled(Un0, 1./dt), DF);
//     velocity_dyn.set_DF(DF);
//     iter.init();
//     getfem::standard_solve(MSL, velocity_dyn, iter);

//     gmm::copy(velocity.get_velocity(MSL), Un1);
//     gmm::copy(velocity.get_pressure(MSL), Pn1);
   
//     pdef->validate_solution(*this, t);     

//     gmm::copy(Un1, Un0); gmm::copy(Pn1, Pn0);
//     do_export(t);
//   }
// }

// /************************************************************/
// void navier_stokes_problem::solve_PREDICTION_CORRECTION() {
//   // Velocity brick.  
//   getfem::mdbrick_generic_elliptic<> velocity(mim, mf_u, nu);
//   getfem::mdbrick_NS_uuT<> velocity_nonlin(velocity);

//   // Volumic source term
//   getfem::mdbrick_source_term<> velocity_f(velocity_nonlin, 
// 					   mf_rhs, 
// 					   plain_vector(mf_rhs.nb_dof()*N));

//   // Dirichlet condition brick.
//   getfem::mdbrick_Dirichlet<> velocity_dir(velocity_f, DIRICHLET_BOUNDARY_NUM);
  
//   // Normal part Dirichlet condition brick.
//   getfem::mdbrick_normal_component_Dirichlet<>
//     velocity_dirnp(velocity_dir, NORMAL_PART_DIRICHLET_BOUNDARY_NUM, mf_rhs);
  
//   // Non-reflective condition brick
//   getfem::mdbrick_NS_nonref1<> nonreflective( velocity_dirnp, 
// 					      NONREFLECTIVE_BOUNDARY_NUM, dt);
  

//   // Dynamic brick.
//     getfem::mdbrick_dynamic<> velocity_dyn(nonreflective, 1.);
//   velocity_dyn.set_dynamic_coeff(1.0/dt, 1.0);

//   // 
//   // Poisson problem for prediction correction scheme
//   //
  
//   getfem::mdbrick_generic_elliptic<> poisson(mim, mf_p, 1.0);
//   sparse_matrix G(1, mf_p.nb_dof()); G(0,0) = 1.;  
//   getfem::mdbrick_constraint<> poisson_setonedof(poisson);
//   poisson_setonedof.set_constraints(G, plain_vector(1));
//   poisson_setonedof.set_constraints_type(getfem::ELIMINATED_CONSTRAINTS);

//   getfem::mdbrick_source_term<>
//     poisson_source(poisson_setonedof, mf_rhs, plain_vector(mf_rhs.nb_dof()));

//   sparse_matrix B(mf_p.nb_dof(), mf_u.nb_dof());
//   asm_stokes_B(B, mim, mf_u, mf_p);

//   // 
//   // dynamic problem
//   //
//   plain_vector DF(mf_u.nb_dof()), F(mf_rhs.nb_dof()), 
//     USTAR(mf_u.nb_dof()), USTARbis(mf_u.nb_dof());

//   /////// REPRISE DES CALCULS EVENTUELS //////
//   /////// utiliser gmm::vecload("sortie.U100",Un0) 
//   //////  gmm::vecload("sortie.P100",Pn0)
//   ////// initialiser t et eviter les 2 initialisations suivantes

//   pdef->initial_condition_u(*this, Un0);
//   pdef->initial_condition_p(*this, Pn0);

//   gmm::iteration iter(residual, noisy);
//   getfem::standard_model_state MSL(velocity_dyn), MSM(poisson_source);
  
//   do_export(0);
//   for (scalar_type t = dt; t <= T; t += dt) {

//     velocity_nonlin.set_U0(Un0);
//     nonreflective.set_Un(Un0);
//     pdef->source_term(*this, t, F);
//     velocity_f.source_term().set(F);
//     pdef->dirichlet_condition(*this, t, F);
//     velocity_dir.rhs().set(mf_rhs, F);
    
//     gmm::mult(velocity_dyn.get_M(), gmm::scaled(Un0, 1./dt), DF);
//     gmm::mult_add(gmm::transposed(B), gmm::scaled(Pn0, -1.), DF);
//     velocity_dyn.set_DF(DF);
//     iter.init();
//     getfem::standard_solve(MSL, velocity_dyn, iter);

//     gmm::copy(velocity.get_solution(MSL), USTAR);
//     gmm::mult(B, USTAR, Pn1);
 
//     cout << "PN1 : " << gmm::vect_norm2(Pn1) << endl;
//     cout << "U1 - USTAR = " << gmm::vect_dist2(Un0, USTAR);

//     poisson_source.set_auxF(Pn1);
//     iter.init();
//     getfem::standard_solve(MSM, poisson_source, iter);
//     gmm::mult(gmm::transposed(B), gmm::scaled(poisson.get_solution(MSM), -1.),
// 	      USTARbis);

//     gmm::iteration iter2 = iter; iter2.reduce_noisy(); iter2.init();
//     gmm::cg(velocity_dyn.get_M(), Un1, USTARbis,
// 	    gmm::identity_matrix(), iter2);

//     gmm::add(USTAR, Un1);
//     gmm::add(gmm::scaled(poisson.get_solution(MSM), 1./dt), Pn0, Pn1);
    
//     pdef->validate_solution(*this, t);
//     cout << "U1 - Un1 = " << gmm::vect_dist2(Un1, Un0);
    
//     gmm::copy(Un1, Un0); gmm::copy(Pn1, Pn0);
//     do_export(t);

//  }
// }


/******************************************************************************************/
void navier_stokes_problem::solve_PREDICTION_CORRECTION2() {

  size_type nbdof_u = mf_u.nb_dof(), nbdof_p = mf_p.nb_dof();
  size_type nbdof_rhs = mf_rhs.nb_dof();
  getfem::mesh_region mpirg = mf_u.linked_mesh().get_mpi_region();
  getfem::mesh_region mpirgp = mf_p.linked_mesh().get_mpi_region();
  gmm::sub_interval I1(0, nbdof_u);

  cout << "nbdof rhs = " << nbdof_rhs << endl;
	cout << "h = "<< mesh.minimal_convex_radius_estimate()<< endl;

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

  // Boundary terms p v.n
  sparse_matrix Bbc(nbdof_p, nbdof_u);
  asm_B_boundary(Bbc, mim, mf_u, mf_p,mf_u.linked_mesh().get_mpi_sub_region(DIRICHLET_BOUNDARY_NUM));
  asm_B_boundary(Bbc, mim, mf_u, mf_p, mf_u.linked_mesh().get_mpi_sub_region(NONREFLECTIVE_BOUNDARY_NUM));
  asm_B_boundary(Bbc, mim, mf_u, mf_p, mf_u.linked_mesh().get_mpi_sub_region(NORMAL_PART_DIRICHLET_BOUNDARY_NUM));
  asm_B_boundary(Bbc, mim, mf_u, mf_p,mf_u.linked_mesh().get_mpi_sub_region(ON_CYLINDER_BOUNDARY_NUM));

  if (N==3){
    asm_B_boundary(Bbc, mim, mf_u, mf_p,mf_u.linked_mesh().get_mpi_sub_region(NEUMANN_BOUNDARY_NUM)); 
  }

//////////////////////////////////////////////////////////////////////////
// To take into account the BOUNDAY CONDITIONS (Lagrange multipliers)
//////////////////////////////////////////////////////////////////////////
  mf_mult.set_qdim(N);
  GMM_ASSERT1(!mf_rhs.is_reduced(), "To be adapted");
  
  dal::bit_vector dofon_nonref = mf_mult.basic_dof_on_region(NONREFLECTIVE_BOUNDARY_NUM);
  
  dal::bit_vector dofon_Dirichlet_Out_Cylinder = mf_mult.basic_dof_on_region(DIRICHLET_BOUNDARY_NUM);
  
	
	
	
//(NORMAL_PART_DIRICHLET_BOUNDARY_NUM);  
//dofon_NDirichlet.setminus(dofon_nonref);
//dofon_NDirichlet.setminus(dofon_Dirichlet_Out_Cylinder);
//(NEUMANN_BOUNDARY_NUM);
//dofon_Neumann.setminus(dofon_Dirichlet_On_cylinder);
//dofon_Neumann.setminus(dofon_nonref);
//dofon_Neumann.setminus(dofon_Dirichlet_Out_Cylinder);
	
//*******************************************************************************************************  
// Normal part Dirichlet condition -- sur v en 2D  -- sur v et w en 3D
//*******************************************************************************************************  
//cout << dofon_Dirichlet_Out_Cylinder << endl;
//cout << dofon_nonref << endl;
//
//  mf_mult.set_qdim(1);
  dal::bit_vector dofon_NDirichlet = mf_mult.basic_dof_on_region(NORMAL_PART_DIRICHLET_BOUNDARY_NUM);
  dofon_NDirichlet.setminus(dofon_nonref);
  dofon_NDirichlet.setminus(dofon_Dirichlet_Out_Cylinder);
	
  std::vector<size_type> ind_ct_ndir;
  for (dal::bv_visitor i(dofon_NDirichlet); !i.finished(); ++i) {
    if (dofon_Dirichlet_Out_Cylinder.is_in(i) || dofon_nonref.is_in(i))   
      dofon_NDirichlet.sup(i);  // Suppress i because it is on the
    else                        // Dirichlet or non reflective boundary.
      ind_ct_ndir.push_back(i);
  } 
  size_type nbdof_NDir = dofon_NDirichlet.card();
  gmm::sub_index SUB_CT_NDIR(ind_ct_ndir);

  gmm::sub_interval I2(nbdof_u, nbdof_NDir);
  getfem::mesh_region mpindirrg
    =mf_u.linked_mesh().get_mpi_sub_region(NORMAL_PART_DIRICHLET_BOUNDARY_NUM);

// sparse_matrix HND1(nbdof_NDir, nbdof_u);
//  {
//    sparse_matrix A(mf_mult.nb_dof(), nbdof_u); 
//    getfem::generic_assembly assem;
//    // assem.set("M(#2,#1)+=comp(Base(#2).vBase(#1).Normal())(:,:,i,i);"); - Good in 2D / Bad in 3D 
//    assem.set("M(#2,#1)+=comp(Base(#2).vBase(#1))(:,:,2);"); // 2ieme composante V = (u,v,w) : v = 0
//    assem.push_mi(mim); assem.push_mf(mf_u); assem.push_mf(mf_mult);
//    assem.push_mat(A); assem.assembly(mpindirrg);
//    gmm::copy(gmm::sub_matrix(A, SUB_CT_NDIR, I1), HND1);
//  }
//  cout << "Nb of Normal part Dirichlet constraints (qdim=1): " << nbdof_NDir << endl;
// 
//  sparse_matrix HND;
  if (N==2){
    //gmm::resize(HND,nbdof_NDir, nbdof_u);
    //gmm::copy(HND1,HND);
    I2 = gmm::sub_interval(nbdof_u, nbdof_NDir);
  }
//  else { //(N==3){
//    // Normal part Dirichlet condition Left - Right (in 3D only)
//    sparse_matrix HND2(nbdof_NDir, nbdof_u);
//    {
//      sparse_matrix A(mf_mult.nb_dof(), nbdof_u); 
//      getfem::generic_assembly assem;
//      assem.set("M(#2,#1)+=comp(Base(#2).vBase(#1))(:,:,3);"); // 3ieme composante V = (u,v,w) : w = 0 
//      assem.push_mi(mim); assem.push_mf(mf_u); assem.push_mf(mf_mult);
//      assem.push_mat(A); assem.assembly(mpindirrg);
//      gmm::copy(gmm::sub_matrix(A, SUB_CT_NDIR, I1), HND2);
//    }
//    gmm::resize(HND,2*nbdof_NDir, nbdof_u);
//    gmm::copy(HND1,gmm::sub_matrix(HND,gmm::sub_interval(0,nbdof_NDir),I1) );
//
//    gmm::copy(HND2,gmm::sub_matrix(HND,gmm::sub_interval(nbdof_NDir,nbdof_NDir),I1) );
//
//    nbdof_NDir = 2*nbdof_NDir;
//    I2 = gmm::sub_interval(nbdof_u, nbdof_NDir);
//  }
// 
//	
//  //M gmm::dense_matrix<double> DM(nbdof_NDir,nbdof_u);
//  //M gmm:: copy(HND,DM);
//  //M for(unsigned i=0;i<nbdof_NDir;++i){
//  //M    for(unsigned j=0;j<nbdof_u;++j){
//  //M      if (DM(i,j)!=0.0) cout << i<<" , "<<j<<" = "<<DM(i,j)<< endl;
//  //M   }
//  //M  }
//
//  // Dirichlet condition except on cylinder
//  mf_mult.set_qdim(N);
//  GMM_ASSERT1(!mf_rhs.is_reduced(), "To be adapted");
	//size_type nbdof_NDir = 0;
//*******************************************************************************************************  
  
  //dal::bit_vector dofon_Dirichlet_Out_Cylinder 
  //  = mf_mult.basic_dof_on_region(DIRICHLET_BOUNDARY_NUM);
  
  std::vector<size_type> ind_ct_dir_out_cyl;
  for (dal::bv_visitor i(dofon_Dirichlet_Out_Cylinder); !i.finished(); ++i) {
    //if (dofon_nonref.is_in(i)) 
    //  dofon_Dirichlet_Out_Cylinder.sup(i);  
    //else  
      ind_ct_dir_out_cyl.push_back(i);
  }
  size_type nbdof_Dir_Out_Cylinder = dofon_Dirichlet_Out_Cylinder.card();
  gmm::sub_index SUB_CT_DIR(ind_ct_dir_out_cyl);
  gmm::sub_interval I3(nbdof_u+nbdof_NDir, nbdof_Dir_Out_Cylinder);
    
  getfem::mesh_region mpidirrg
    = mf_u.linked_mesh().get_mpi_sub_region(DIRICHLET_BOUNDARY_NUM);
  
  sparse_matrix HD(nbdof_Dir_Out_Cylinder, nbdof_u);
  {
    sparse_matrix A(mf_mult.nb_dof(), nbdof_u); 
    getfem::generic_assembly assem;
    assem.set("M(#2,#1)+=comp(vBase(#2).vBase(#1))(:,i,:,i);");
    assem.push_mi(mim); assem.push_mf(mf_u); assem.push_mf(mf_mult);
    assem.push_mat(A); assem.assembly(mpidirrg);
    gmm::copy(gmm::sub_matrix(A, SUB_CT_DIR, I1), HD);
  }
  cout << "Nb of Dirichlet constraints without cylinder : " << nbdof_Dir_Out_Cylinder << endl;


 // Dirichlet condition on cylinder
  mf_mult.set_qdim(N);
  GMM_ASSERT1(!mf_rhs.is_reduced(), "To be adapted");

 dal::bit_vector dofon_Dirichlet_On_Cylinder 
    = mf_mult.basic_dof_on_region(ON_CYLINDER_BOUNDARY_NUM);

  std::vector<size_type> ind_ct_dir_on_cylinder;
  for (dal::bv_visitor i(dofon_Dirichlet_On_Cylinder); !i.finished(); ++i) 
    ind_ct_dir_on_cylinder.push_back(i);

  size_type nbdof_Dir_On_Cylinder = dofon_Dirichlet_On_Cylinder.card();
  gmm::sub_index SUB_CT_DIR_CYL(ind_ct_dir_on_cylinder);
  gmm::sub_interval I3C(nbdof_u+nbdof_NDir+nbdof_Dir_Out_Cylinder, nbdof_Dir_On_Cylinder);
  
  getfem::mesh_region mpidircylrg
    = mf_u.linked_mesh().get_mpi_sub_region(ON_CYLINDER_BOUNDARY_NUM);

  sparse_matrix HDC(nbdof_Dir_On_Cylinder, nbdof_u);
  {
    sparse_matrix A(mf_mult.nb_dof(), nbdof_u); 
    getfem::generic_assembly assem;
    assem.set("M(#2,#1)+=comp(vBase(#2).vBase(#1))(:,i,:,i);");
    assem.push_mi(mim); assem.push_mf(mf_u); assem.push_mf(mf_mult);
    assem.push_mat(A); assem.assembly(mpidircylrg);
    gmm::copy(gmm::sub_matrix(A, SUB_CT_DIR_CYL, I1), HDC);
  }
  cout << "Nb of Dirichlet constraints on cylinder : " << nbdof_Dir_On_Cylinder << endl;
  //??//gmm :: resize(lambda,  nbdof_Dir_On_Cylinder);


  // Non reflective condition
  mf_mult.set_qdim(N);
  GMM_ASSERT1(!mf_rhs.is_reduced(), "To be adapted");
  //   dal::bit_vector dofon_nonref
  //  = mf_mult.basic_dof_on_region(NONREFLECTIVE_BOUNDARY_NUM);

  std::vector<size_type> ind_ct_nonref;
  for (dal::bv_visitor i(dofon_nonref); !i.finished(); ++i) {
    //if (dofon_NDirichlet.is_in(i)) 
    //  dofon_nonref.sup(i); 
    //else  
      ind_ct_nonref.push_back(i);}
  size_type nbdof_nonref = dofon_nonref.card();
  gmm::sub_index SUB_CT_NONREF(ind_ct_nonref);
  gmm::sub_interval I4(nbdof_u+nbdof_NDir+nbdof_Dir_Out_Cylinder+nbdof_Dir_On_Cylinder, nbdof_nonref);
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


  // In 3D - NEUMANN_BOUNDARY_NUM
  sparse_matrix HN;
  size_type nbdof_neumann=0;
  gmm::sub_interval I5;
  
  dofon_nonref= mf_mult.basic_dof_on_region(NONREFLECTIVE_BOUNDARY_NUM);
  dofon_Dirichlet_Out_Cylinder = mf_mult.basic_dof_on_region(DIRICHLET_BOUNDARY_NUM);
  //dofon_NDirichlet = mf_mult.basic_dof_on_region(NORMAL_PART_DIRICHLET_BOUNDARY_NUM);
  dofon_Dirichlet_On_Cylinder = mf_mult.basic_dof_on_region(ON_CYLINDER_BOUNDARY_NUM);

  //if (N==3){
//    mf_mult.set_qdim(1);
//
//    dal::bit_vector dofon_neumann = mf_mult.basic_dof_on_region(NEUMANN_BOUNDARY_NUM);
//
//    std::vector<size_type> ind_ct_neumann;
//    for (dal::bv_visitor i(dofon_neumann); !i.finished(); ++i) {
//      if (dofon_Dirichlet_Out_Cylinder.is_in(i*N) || dofon_nonref.is_in(i*N) ||
//	  dofon_NDirichlet.is_in(i*N) || dofon_Dirichlet_On_Cylinder.is_in(i*N) )   
//	dofon_neumann.sup(i); 
//      else                     
//	ind_ct_neumann.push_back(i);
//    } 
//    //cout << mf_mult.point_of_basic_dof(i)<< endl;
//    
//    nbdof_neumann = dofon_neumann.card();
//    gmm::sub_index SUB_CT_NEUMANN(ind_ct_neumann);
//    
//    I5 = gmm::sub_interval(nbdof_u+nbdof_NDir+nbdof_Dir_Out_Cylinder+nbdof_Dir_On_Cylinder+nbdof_nonref, nbdof_neumann);
//    getfem::mesh_region mpineumannrg = mf_u.linked_mesh().get_mpi_sub_region(NEUMANN_BOUNDARY_NUM);
//        
//    gmm::resize(HN,nbdof_neumann, nbdof_u);
//    {
//      sparse_matrix A(mf_mult.nb_dof(), nbdof_u); 
//      getfem::generic_assembly assem;
//      assem.set("M(#2,#1)+=comp(Base(#2).vBase(#1))(:,:,3);"); // 3ieme composante V = (u,v,w) : w = 0
//      assem.push_mi(mim); assem.push_mf(mf_u); assem.push_mf(mf_mult);
//      assem.push_mat(A); assem.assembly(mpineumannrg);
//      gmm::copy(gmm::sub_matrix(A, SUB_CT_NEUMANN, I1), HN);
//    }
//    cout << "Nb of Neumann part in 3D  (qdim=1): " << nbdof_neumann << endl;
//    
//    mf_mult.set_qdim(N);
//
//
//  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
// Operators to solve NSE
/////////////////////////////////////////////////////////////////////////////////////////////////

  // Discretization of laplace operator for p
  sparse_matrix K2(nbdof_p+1, nbdof_p+1);
  gmm::sub_interval IP(0, nbdof_p);
  asm_stiffness_matrix_for_homogeneous_laplacian(gmm::sub_matrix(K2, IP),
						 mim, mf_p, mpirgp);


  /*
  dal::bit_vector dofon_cylp 
    = mf_p.basic_dof_on_region(ON_CYLINDER_BOUNDARY_NUM);

  getfem::mesh_region mpicylp
    = mf_p.linked_mesh().get_mpi_sub_region(ON_CYLINDER_BOUNDARY_NUM);

  std::vector<size_type> ind_cylp;
  for (dal::bv_visitor i(dofon_cylp); !i.finished(); ++i) 
    ind_cylp.push_back(i);

  //size_type nbdof_cylp = dofon_cylp.card();
  gmm::sub_index SUB_CYLP(ind_cylp);
  


  if (N>=2)
  {
    sparse_matrix A(nbdof_p, nbdof_p); 
    getfem::generic_assembly assem;
    assem.set("M(#1,#1)+=comp(Grad(#1).Normal().Base(#1))(:,1,1,:);");
    assem.push_mi(mim); assem.push_mf(mf_p);
    assem.push_mat(A); assem.assembly(mpicylp);
    //gmm::add(gmm::sub_matrix(A,SUB_CYLP,IP),gmm::sub_matrix(K2,SUB_CYLP,IP));
    gmm::add(A,gmm::sub_matrix(K2,IP,IP));
  }

  if (N==3)
  {
    sparse_matrix A(nbdof_p, nbdof_p); 
    getfem::generic_assembly assem;
    assem.set("M(#1,#1)+=comp(Grad(#1).Normal().Base(#1))(:,3,3,:);");
    assem.push_mi(mim); assem.push_mf(mf_p);
    assem.push_mat(A); assem.assembly(mpicylp);
    //    gmm::add(gmm::sub_matrix(A,SUB_CYLP,IP),gmm::sub_matrix(K2,SUB_CYLP,IP));
    gmm::add(A,gmm::sub_matrix(K2,IP,IP));

  }

  //Prise en compte de façon faible
  plain_vector HP(nbdof_p);
  { plain_vector A(nbdof_p); 
    getfem::generic_assembly assem;
    assem.set("V(#1)+=comp(Base(#1))(:);");
    assem.push_mi(mim); assem.push_mf(mf_p);
    assem.push_vec(A); assem.assembly(mpirgp);
    gmm::copy(A,HP);
  }
 
  // for (unsigned i=0;i<=nbdof_p; ++i) {
  //  HP[i];} // mean value of the pressure = 0
  //K2(nbdof_p,nbdof_p) = 0.0;  
  */

  //Prise en compte de façon forte
  for (unsigned i=0;i<=nbdof_p; ++i) {
    K2(nbdof_p,i) = K2(i,nbdof_p) = 1.0;} // mean value of the pressure = 0
  //K2(nbdof_p, 0) = K2(0, nbdof_p) = 1.0;} // set the first pressure dof to 0
  K2(nbdof_p,nbdof_p) = 0.0;  

  gmm::SuperLU_factor<double> SLUsys2, SLUsys3,SLUsys1;
  SLUsys2.build_with(K2);

	
///////////////////////////////////////////////////////////////////////////////////////////////// 
// dynamic problem
/////////////////////////////////////////////////////////////////////////////////////////////////

  plain_vector DF(nbdof_u), F(nbdof_rhs), USTAR(nbdof_u), USTARbis(nbdof_u);
  plain_vector Phi(nbdof_p);

  /////// REPRISE DES CALCULS EVENTUELS //////
  std::string initfile_Um1 = PARAM.string_value("INIT_FILE_Um1");
  
  std::string initfile_U = PARAM.string_value("INIT_FILE_U");
  std::string initfile_P = PARAM.string_value("INIT_FILE_P");
  
  if (initfile_U.size() == 0 || initfile_P.size() == 0) {
    pdef->initial_condition_u(*this, Un0);
    pdef->initial_condition_p(*this, Pn0);
  }  else {
    gmm::vecload(initfile_U,Un0); 
    gmm::vecload(initfile_P,Pn0);
  }
  
  if (initfile_Um1.size() != 0) {
    gmm::vecload(initfile_Um1,Unm1); 
  }
  
  ////// initialiser t et eviter les 2 initialisations suivantes

  std::ofstream coeffTP("coeffTP.data");
  std::ofstream ptPartData("ptPart.data");

  // Recherche d'un point de réference pour le calcul des coeff de trainée ...

  scalar_type BoxXmin =  PARAM.real_value("BOXXmin", "Particular Point xMin");
  scalar_type BoxXmax =  PARAM.real_value("BOXXmax", "Particular Point xMax");
  scalar_type BoxYmin =  PARAM.real_value("BOXYmin", "Particular Point yMin");
  scalar_type BoxYmax =  PARAM.real_value("BOXYmax", "Particular Point yMax");
  scalar_type BoxZmin, BoxZmax;

  if (N==3) {
    BoxZmin =  PARAM.real_value("BOXZmin", "Particular Point zMin");
    BoxZmax =  PARAM.real_value("BOXZmax", "Particular Point zMax");
  }
	
  base_node ptPartU(N+1), ptPartP(N+1);

	
  GMM_ASSERT1(!mf_p.is_reduced(), "To be adapted");


 if (N==2) {	 
	 for (unsigned i=0; i< mf_p.nb_dof(); ++i){
      bgeot :: base_node BN =  mf_p.point_of_basic_dof(i);
      if (BN[0]<BoxXmax && BN[0] > BoxXmin && BN[1]< BoxYmax && BN[1] > BoxYmin ) {
		  cout << "Point Part in Box -- i on mf_p= " << i <<",x="<<BN[0]<<",y="<<BN[1]<< endl;
		  ptPartP[0] = BN[0];
		  ptPartP[1] = BN[1];
		  ptPartP[2] = i;
		  break;
      }
    } 
	 
	 
	 for (unsigned i=0; i< mf_u.nb_dof(); ++i){
		 bgeot :: base_node BN =  mf_u.point_of_dof(i);
		 if (BN[0]==ptPartP[0] && BN[1]==ptPartP[1] ) {
			 cout << "Point Part in Box -- i on mf_u= " << i <<",x="<<BN[0]<<",y="<<BN[1]<< endl;
			 // Attention c'est vectoriel => en sortie i <-> pt sur la dernière composante de la vitesse
			 // Vitesse =(U,V,W) alors en 2D --> sur V, en 3D --> sur W
			 // D'ou la modif dans ptPartU[2]
			 ptPartU[0] = BN[0];
			 ptPartU[1] = BN[1]; 
			 ptPartU[2] = i - N + 1;
		     break;

		 } 
	 }
	 ptPartData << ptPartP[0] << " " << ptPartP[1] << endl; 
	 
  }
	if (N==3) {	 
		for (unsigned i=0; i< mf_p.nb_dof(); ++i){
			bgeot :: base_node BN =  mf_p.point_of_basic_dof(i);
			if (BN[0]<BoxXmax && BN[0] > BoxXmin && BN[1]< BoxYmax && BN[1] > BoxYmin  && BN[2]< BoxZmax && BN[2] > BoxZmin ) {
				cout << "Point Part in Box -- i on mf_p= " << i <<",x="<<BN[0]<<",y="<<BN[1]<< ",z="<<BN[2]<< endl;
				ptPartP[0] = BN[0];
				ptPartP[1] = BN[1];
				ptPartP[2] = BN[2];
				ptPartP[3] = i;
				break;
			}
		} 
		
		
		for (unsigned i=0; i< mf_u.nb_dof(); ++i){
			bgeot :: base_node BN =  mf_u.point_of_dof(i);
			if (BN[0]==ptPartP[0] && BN[1]==ptPartP[1]&& BN[2]==ptPartP[2] ) {
				cout << "Point Part in Box -- i on mf_u= " << i <<",x="<<BN[0]<<",y="<<BN[1]<<",z="<<BN[2]<< endl;
				// Attention c'est vectoriel => en sortie i <-> pt sur la dernière composante de la vitesse
				// Vitesse =(U,V,W) alors en 2D --> sur V, en 3D --> sur W
				// D'ou la modif dans ptPartU[3]
				ptPartU[0] = BN[0];
				ptPartU[1] = BN[1]; 
				ptPartU[2] = BN[2];
				ptPartU[3] = i - N + 1;
				break;
				
			} 
		}
		ptPartData << ptPartP[0] << " " << ptPartP[1] << " " << ptPartP[2] << endl; 
		
	} 

  size_type sizelsystem = 0;
  size_type sizelsystemP = nbdof_p+1;

  if (N==2) 
    sizelsystem = nbdof_u + nbdof_NDir + nbdof_Dir_Out_Cylinder + nbdof_Dir_On_Cylinder + nbdof_nonref;
  if (N==3)
    sizelsystem = nbdof_u + nbdof_NDir + nbdof_Dir_Out_Cylinder + nbdof_Dir_On_Cylinder + nbdof_nonref;// + nbdof_neumann;

  sparse_matrix A1(sizelsystem, sizelsystem);
  sparse_matrix A2(sizelsystem, sizelsystem); 
  sparse_matrix A2u(sizelsystem/2, sizelsystem/2); 
  sparse_matrix A2v(sizelsystem/2, sizelsystem/2); 
  sparse_matrix A1u(sizelsystem/2, sizelsystem/2); 
  sparse_matrix A1v(sizelsystem/2, sizelsystem/2); 
	
	
  plain_vector    Y(sizelsystem), YY(nbdof_u), Ytmp(nbdof_u);

  do_export(0);
  cout << "-----------------------------------------------------------" << endl;
  cout << "Dynamic problem implementation" << endl;
  cout << "-----------------------------------------------------------" << endl;

	
gmm :: sub_slice SUB_CT_Vu(0,sizelsystem/N,N);
gmm :: sub_slice SUB_CT_Vv(1,sizelsystem/N,N);
gmm :: sub_interval SUB_CT_V(0,nbdof_u);
gmm :: sub_interval SUB_CT_P(0,nbdof_p);
	
// idŽees 	
//gmm::copy(gmm::sub_vector(F1generic,SUB_CT_V1full),F1full);
//gmm::copy(gmm::sub_vector(X1full,SUB_CT_V),X1);
//gmm::copy(gmm::sub_vector(X2full,SUB_CT_V),X2);
//gmm::copy(X1,gmm::sub_vector(USTAR,gmm::sub_slice(0,nbdof_u,2)));
//gmm::copy(X2,gmm::sub_vector(USTAR,gmm::sub_slice(1,nbdof_u,2)));
	
  for (scalar_type t = Tinitial + dt; t <= T; t += dt) {

    //******* Construction of the Matrix for the 3rd system (to obtain velocity) **********//
    if (t==Tinitial+dt){ //time_order = 1 or first iterations with time_order = 2, computed once
      gmm::clear(A2);
      gmm::add(M, gmm::sub_matrix(A2, I1));
      //// Normal Dirichlet condition
      //gmm::copy(HND, gmm::sub_matrix(A2, I2, I1));
      //gmm::copy(gmm::transposed(HND), gmm::sub_matrix(A2, I1, I2));
      
      // Dirichlet condition except on cylinder
      gmm::copy(HD, gmm::sub_matrix(A2, I3, I1));
      gmm::copy(gmm::transposed(HD), gmm::sub_matrix(A2, I1, I3));
      
      // Dirichlet condition on cylinder
      gmm::copy(HDC, gmm::sub_matrix(A2, I3C, I1));
      gmm::copy(gmm::transposed(HDC), gmm::sub_matrix(A2, I1, I3C));
      
      // Non reflective condition
      gmm::copy(HNR, gmm::sub_matrix(A2, I4, I1));
      gmm::copy(gmm::transposed(HNR), gmm::sub_matrix(A2, I1, I4));
      
      // Factorization LU
      // SLUsys3.build_with(A2); // pour le systeme couple (u,v)
		
		
		
		gmm::clear(A2u);
		//gmm::clear(A2v);

		
		gmm::copy(gmm::sub_matrix(A2,SUB_CT_Vu,SUB_CT_Vu),A2u);
		//gmm::copy(gmm::sub_matrix(A2,SUB_CT_Vv,SUB_CT_Vv),A2v);

	
		//gmm::add(gmm::scaled(A2u,-1.0),A2v);
		//cout <<"A2 "<< gmm::mat_norminf(A2v) << endl;
		
	SLUsys3.build_with(A2u);
		
    }

    if ((time_order==2)&&(t==Tinitial+2*dt)){ // computed once
      gmm::clear(A2);
      gmm::add(gmm::scaled(M,1.5), gmm::sub_matrix(A2, I1));
      //// Normal Dirichlet condition
      //gmm::copy(HND, gmm::sub_matrix(A2, I2, I1));
      //gmm::copy(gmm::transposed(HND), gmm::sub_matrix(A2, I1, I2));
      
      // Dirichlet condition except on cylinder
      gmm::copy(HD, gmm::sub_matrix(A2, I3, I1));
      gmm::copy(gmm::transposed(HD), gmm::sub_matrix(A2, I1, I3));
      
      // Dirichlet condition on cylinder
      gmm::copy(HDC, gmm::sub_matrix(A2, I3C, I1));
      gmm::copy(gmm::transposed(HDC), gmm::sub_matrix(A2, I1, I3C));
      
      // Non reflective condition
      gmm::copy(HNR, gmm::sub_matrix(A2, I4, I1));
      gmm::copy(gmm::transposed(HNR), gmm::sub_matrix(A2, I1, I4));
      
      // Factorization LU
      //SLUsys3.build_with(A2); // pb en (u,v) couplŽ
      
		gmm::clear(A2u);
		//gmm::clear(A2v);
		
		
		gmm::copy(gmm::sub_matrix(A2,SUB_CT_Vu,SUB_CT_Vu),A2u);
		//gmm::copy(gmm::sub_matrix(A2,SUB_CT_Vv,SUB_CT_Vv),A2v);
		
		SLUsys3.build_with(A2u);
	
		
    }
    
  
    //******* The matrix of the 3rd system is constructed and factorized
   
  /*
    //  cout<<"cylinder"<<gmm::sub_vector(Un0,SUB_CT_DIR_CYL) << endl;

    //cout<<"dir"<<gmm::sub_vector(Un0,SUB_CT_DIR) << endl;
  
  dal::bit_vector dofon_NDirichlet_tmp
    = mf_mult.basic_dof_on_region(NORMAL_PART_DIRICHLET_BOUNDARY_NUM);
  std::vector<size_type> ind_ct_ndir_tmp;
  for (dal::bv_visitor i(dofon_NDirichlet_tmp); !i.finished(); ++i) {
      ind_ct_ndir_tmp.push_back(i);
  } 
  gmm::sub_index SUB_CT_NDIR_tmp(ind_ct_ndir_tmp);
    cout<<"NdirV"<<gmm::sub_vector(Un0,SUB_CT_NDIR_tmp) << endl;
  

    // if (N==3){cout<<"NdirW"<<gmm::sub_vector(Un0,SUB_CT_NDIR_W) << endl;};

    // cout<<"Nonref"<<gmm::sub_vector(Un0,SUB_CT_NONREF) << endl;
    */

    // Assembly of the first linear system
   
    if ((time_order==1)||(t==Tinitial+dt)){ //time_order = 1 or first iterations with time_order = 2, computed once
      gmm::clear(A1);
      gmm::clear(Y);
      // laplace operator
      gmm::copy(K1, gmm::sub_matrix(A1, I1));     
      // Nonlinear part
      //cout << "avant  asm_NS_uuT"<< endl;
      getfem::asm_NS_uuT(gmm::sub_matrix(A1, I1), mim, mf_u, Un0, mpirg);   
      //cout << "apres  asm_NS_uuT"<< endl;
      // Dynamic part 
      gmm::add(gmm::scaled(M, 1./dt), gmm::sub_matrix(A1, I1));
		
      gmm::mult(M, gmm::scaled(Un0, 1.0/dt), gmm::sub_vector(Y, I1));
	
	}

    if ((time_order==2)&&(t>=Tinitial+2*dt)){
      gmm::clear(A1);
      gmm::clear(Y);
      gmm::clear(Ytmp);
      // laplace operator
      gmm::copy(K1, gmm::sub_matrix(A1, I1));     
      // Nonlinear part
      //cout << "avant  asm_NS_uuT"<< endl;
      getfem::asm_NS_uuT(gmm::sub_matrix(A1, I1), mim, mf_u, Un0, mpirg);  
      //cout << "apres  asm_NS_uuT"<< endl;
      // Dynamic part 
      gmm::add(gmm::scaled(M, 1.5/dt), gmm::sub_matrix(A1, I1));
      gmm::add(gmm::scaled(Un0, 2.0/dt),gmm::scaled(Unm1,-0.5/dt),Ytmp);
      gmm::mult(M, Ytmp, gmm::sub_vector(Y, I1));
		
    }


#if GETFEM_PARA_LEVEL > 1
    MPI_SUM_VECTOR2(Y);
#endif
   
    // Volumic source term -- inutile d'assambler car F = 0
    //pdef->source_term(*this, t, F);
    //getfem::asm_source_term(Y, mim, mf_u, mf_rhs, F,mpirg); // inutile (F=0)
   
    // // Normal Dirichlet condition
    //gmm::copy(HND, gmm::sub_matrix(A1, I2, I1));
    //gmm::copy(gmm::transposed(HND), gmm::sub_matrix(A1, I1, I2));
   
    // Dirichlet condition except on cylinder
    gmm::copy(HD, gmm::sub_matrix(A1, I3, I1));
    gmm::copy(gmm::transposed(HD), gmm::sub_matrix(A1, I1, I3));

    // Dirichlet condition on cylinder
    gmm::copy(HDC, gmm::sub_matrix(A1, I3C, I1));
    gmm::copy(gmm::transposed(HDC), gmm::sub_matrix(A1, I1, I3C));

    // Non reflective condition
    gmm::copy(HNR, gmm::sub_matrix(A1, I4, I1));
    gmm::copy(gmm::transposed(HNR), gmm::sub_matrix(A1, I1, I4));
		
    
   //  if (N==3){
//     // Neumann condition in 3D
//       gmm::copy(HN, gmm::sub_matrix(A1, I5, I1));
//       gmm::copy(gmm::transposed(HN), gmm::sub_matrix(A1, I1, I5));
      
//       gmm::copy(HN, gmm::sub_matrix(A2, I5, I1));
//       gmm::copy(gmm::transposed(HN), gmm::sub_matrix(A2, I1, I5));
//     }
    
    // Right hand side with Dirichlet boundary conditions

    gmm::resize(F, N * nbdof_rhs);
    pdef->dirichlet_condition(*this, t, F); // time independent => exit loop

    {
      plain_vector VV(mf_mult.nb_dof());
      gmm::clear(VV);
      getfem::asm_source_term(VV, mim, mf_mult, mf_rhs, F, mpidirrg); // a optimiser indept time
#if GETFEM_PARA_LEVEL > 1
      MPI_SUM_VECTOR2(VV);
#endif
      gmm::copy(gmm::sub_vector(VV, SUB_CT_DIR), gmm::sub_vector(Y, I3));
    }
    {
      plain_vector VV(mf_mult.nb_dof());
      gmm::clear(VV);
      getfem::asm_source_term(VV, mim, mf_mult, mf_rhs, F, mpidircylrg); // a optimiser  indept time
#if GETFEM_PARA_LEVEL > 1
      MPI_SUM_VECTOR2(VV);
#endif
      gmm::copy(gmm::sub_vector(VV, SUB_CT_DIR_CYL), gmm::sub_vector(Y, I3C));
    }

    /*{
      plain_vector VV(mf_mult.nb_dof());
      getfem::asm_source_term(VV, mim, mf_mult, mf_rhs, F, mpindirrg); // a optimiser  indept time
      #if GETFEM_PARA_LEVEL > 1
      MPI_SUM_VECTOR2(VV);
      #endif
      gmm::copy(gmm::sub_vector(VV, SUB_CT_NDIR), gmm::sub_vector(Y, I2));       
      }*/

       {
      plain_vector VV(mf_mult.nb_dof());
      gmm::clear(VV);
      // getfem::asm_source_term(VV, mim, mf_mult, mf_rhs, F, mpinonrefrg); // CL Périodiques
      if (non_reflective_bc == 0) {
 	asm_basic_non_reflective_bc(VV, mim, mf_u, Un0, mf_mult, dt, mpinonrefrg);
      } else {
 	asm_improved_non_reflective_bc(VV, mim, mf_u, Un0, mf_mult, dt, nu, mpinonrefrg); 
      }

#if GETFEM_PARA_LEVEL > 1
	MPI_SUM_VECTOR2(VV);
#endif   
      gmm::copy(gmm::sub_vector(VV, SUB_CT_NONREF), gmm::sub_vector(Y, I4));
    }
    
    // pressure part
    gmm::clear(YY);
    gmm::mult(gmm::transposed(B), gmm::scaled(Pn0, -1.0), YY);

#if GETFEM_PARA_LEVEL > 1
    MPI_SUM_VECTOR2(YY);
#endif
    
    gmm::add(YY, gmm::sub_vector(Y, I1));
    
    gmm::clear(YY);
    gmm::mult(gmm::transposed(Bbc), Pn0, YY);
 

#if GETFEM_PARA_LEVEL > 1
    MPI_SUM_VECTOR2(YY);
#endif
    gmm::add(YY, gmm::sub_vector(Y, I1));
    
	  
	  gmm::clear(A1u);
	 // gmm::clear(A1v);
	  
	  
	  gmm::copy(gmm::sub_matrix(A1,SUB_CT_Vu,SUB_CT_Vu),A1u);
	 // gmm::copy(gmm::sub_matrix(A1,SUB_CT_Vv,SUB_CT_Vv),A1v);
	  
	  
	  
    //
    // Solving the first linear system
    //
	  
    {
      double rcond;
      plain_vector X(sizelsystem);

		//plain_vector Xu(sizelsystem/2);
		//plain_vector Xv(sizelsystem/2);
		//plain_vector Yu(sizelsystem/2);
		//plain_vector Yv(sizelsystem/2);

		
#if  (GETFEM_PARA_LEVEL > 1 && GETFEM_PARA_SOLVER == MUMPS_PARA_SOLVER)
      MUMPS_distributed_matrix_solve(A1,X,Y);
#elif (GETFEM_PARA_LEVEL==1 && GMM_USES_MUMPS)
      MUMPS_solve(A1,X,Y);
      //#elif (GETFEM_PARA_LEVEL==0 && GMM_USES_MUMPS)
      //MUMPS_solve(A1,X,Y);
#elif (GETFEM_PARA_LEVEL==0)
      // SuperLU_solve(A1, X, Y, rcond);
	  
		//gmm::copy(gmm::sub_vector(Y,SUB_CT_Vu),Yu);
	    //gmm::copy(gmm::sub_vector(Y,SUB_CT_Vv),Yv);

		
		//SuperLU_solve(A1u, Xu, Yu, rcond);
		//SuperLU_solve(A1u, Xv, Yv, rcond);
		
		// Factorisation LU
		SLUsys1.build_with(A1u);
		//SLUsys1.solve(Xu,Yu);
		//SLUsys1.solve(Xv,Yv);
		SLUsys1.solve(gmm::sub_vector(X,SUB_CT_Vu),gmm::sub_vector(Y,SUB_CT_Vu));
		SLUsys1.solve(gmm::sub_vector(X,SUB_CT_Vv),gmm::sub_vector(Y,SUB_CT_Vv));
		
		
		//gmm::copy(gmm::sub_vector(Un0,gmm::sub_slice(0,nbdof_u/2,2)),gmm::sub_vector(Xu,gmm::sub_interval(0,nbdof_u/2)));
	    //gmm::copy(gmm::sub_vector(Un0,gmm::sub_slice(1,nbdof_u/2,2)),gmm::sub_vector(Xv,gmm::sub_interval(0,nbdof_u/2)));
		//gmm::iteration iter(1E-8);
		//gmm::bicgstab(A1u,Xu,Yu,gmm::identity_matrix(),iter);				
		//gmm::bicgstab(A1u,Xv,Yv,gmm::identity_matrix(),iter);
		
		
		
		//gmm::copy(Xu,gmm::sub_vector(X,SUB_CT_Vu));
		//gmm::copy(Xv,gmm::sub_vector(X,SUB_CT_Vv));


		//gmm::add(gmm::scaled(Xfull,-1.0),X);
		//cout << gmm::vect_norminf(X) << endl;
		
      // gmm::iteration iter(1E-8);
      // gmm::gmres(A1,X,Y,gmm::identity_matrix(),10,iter);
      // gmm::bicgstab(A1,X,Y,gmm::identity_matrix(),iter);
      // ?? gmm::bicgstab(A1,X,Y,gmm::diagonal_precond<sparse_matrix>(A1),iter);
      // if (noisy) cout << "condition number: " << 1.0/rcond << endl;
#endif
      gmm::copy(gmm::sub_vector(X, I1), USTAR);
      //??//  gmm::copy(gmm::sub_vector(X, I3C), lambda);
		
		// Relation de compatibilitŽ int_domaine div(ustar)=0 //
		scalar_type delta_in;
		sparse_matrix Bbc_flux_in(nbdof_p, nbdof_u);
		asm_B_boundary(Bbc_flux_in, mim, mf_u, mf_p,mf_u.linked_mesh().get_mpi_sub_region(DIRICHLET_BOUNDARY_NUM));
		
		gmm :: resize(lambda, nbdof_p); gmm :: clear(lambda); gmm::fill(lambda,1.0);
		
		gmm :: resize(tmp, nbdof_Dir_Out_Cylinder+nbdof_Dir_On_Cylinder); gmm :: clear(tmp);
		
		gmm :: mult(Bbc_flux_in,lambda,tmp);
		gmm :: resize(lambda, nbdof_Dir_On_Cylinder+nbdof_Dir_Out_Cylinder); gmm :: clear(lambda); 
		gmm :: copy(gmm::sub_vector(USTAR,I3),gmm::sub_vector(lambda,gmm::sub_interval(0,nbdof_Dir_Out_Cylinder)));
		gmm :: copy(gmm::sub_vector(USTAR,I3C),gmm::sub_vector(lambda,gmm::sub_interval(nbdof_Dir_Out_Cylinder,nbdof_Dir_On_Cylinder)));
		delta_in = gmm :: vect_sp(tmp,lambda);
		
		scalar_type delta_out;

		sparse_matrix Bbc_flux_out(nbdof_p, nbdof_u);
		asm_B_boundary(Bbc_flux_out, mim, mf_u, mf_p, mf_u.linked_mesh().get_mpi_sub_region(NONREFLECTIVE_BOUNDARY_NUM));
		asm_B_boundary(Bbc_flux_out, mim, mf_u, mf_p, mf_u.linked_mesh().get_mpi_sub_region(NORMAL_PART_DIRICHLET_BOUNDARY_NUM));
		//asm_B_boundary(Bbc_flux_out, mim, mf_u, mf_p,mf_u.linked_mesh().get_mpi_sub_region(ON_CYLINDER_BOUNDARY_NUM));
		//asm_B_boundary(Bbc_flux_in, mim, mf_u, mf_p,mf_u.linked_mesh().get_mpi_sub_region(DIRICHLET_BOUNDARY_NUM));
		gmm :: resize(lambda, nbdof_p); gmm :: clear(lambda); gmm::fill(lambda,1.0);
		
		gmm :: resize(tmp, nbdof_NDir+nbdof_nonref); gmm :: clear(tmp);
		
		gmm :: mult(Bbc_flux_out,lambda,tmp);
		gmm :: resize(lambda, nbdof_NDir+nbdof_nonref); gmm :: clear(lambda); 
		gmm :: copy(gmm::sub_vector(USTAR,I2),gmm::sub_vector(lambda,gmm::sub_interval(0,nbdof_NDir)));
		gmm :: copy(gmm::sub_vector(USTAR,I4),gmm::sub_vector(lambda,gmm::sub_interval(nbdof_NDir,nbdof_nonref)));
		delta_out = gmm :: vect_sp(tmp,lambda);
		
		cout<<"compatibilite" << delta_in << "  "<< delta_out << "  "<< delta_in-delta_out << endl;
		gmm::scaled(gmm::sub_vector(USTAR,I2),delta_in/delta_out); 
		gmm::scaled(gmm::sub_vector(USTAR,I4),delta_in/delta_out); 
		gmm::scaled(gmm::sub_vector(USTAR,I2),delta_in/delta_out); 

		if (N==3){
			asm_B_boundary(Bbc_flux_out, mim, mf_u, mf_p,mf_u.linked_mesh().get_mpi_sub_region(NEUMANN_BOUNDARY_NUM)); 
		}
		
		
		//delta = 0.0;
		//for (int ii = 0; ii <= nbdof_nonref; ii += 2) {
		//   delta = delta+lambda[ii];
		//}
	    //delta=delta/nbdof_nonref;
		//for (int ii = 0; ii <= nbdof_nonref; ii += 2) {
		//	lambda[ii] = 1.0-delta;
		//}
		//cout << "delta" << delta << endl;
		//gmm::add(lambda,gmm::sub_vector(X, I4));
		//gmm::copy(gmm::sub_vector(X, I1), USTAR);
		
    }

    cout << "U* - Un0 = " << gmm::vect_dist2(USTAR, Un0) << endl;


 if (N==3){
    plain_vector DU(mf_rhs.nb_dof() * N * N);
    plain_vector DIV(mf_rhs.nb_dof());
    compute_gradient(mf_u, mf_rhs, USTAR, DU);
    for (unsigned i=0; i < mf_rhs.nb_dof(); ++i) {
      DIV[i] = DU[i*N*N] +DU[i*N*N+4] +DU[i*N*N+8];
    }
    if (PARAM.int_value("VTK_EXPORT")) {
      static int cnta=0;
      char sa[128]; sprintf(sa, "SolIcare/DivUstar%d.vtk", cnta++);
      getfem::vtk_export tata(sa, PARAM.int_value("VTK_EXPORT")==1);
      tata.exporting(mf_rhs); 
      tata.write_point_data(mf_rhs, DIV, "divergence");

      char sb[128]; sprintf(sb, "SolIcare/Ustar%d.vtk", cnta++);
      getfem::vtk_export tbtb(sb, PARAM.int_value("VTK_EXPORT")==1);
      tbtb.exporting(mf_u); 
      tbtb.write_point_data(mf_u, USTAR, "Ustar");

    }

    }

    //
    // Solving the second linear system
    //
    {
      //double rcond;
      plain_vector X(sizelsystemP), Z(sizelsystemP);
      gmm::mult(B, USTAR, gmm::sub_vector(Z, IP));
#if GETFEM_PARA_LEVEL > 1
      MPI_SUM_VECTOR2(Z);
#endif

#if  (GETFEM_PARA_LEVEL > 1 && GETFEM_PARA_SOLVER == MUMPS_PARA_SOLVER)
      MUMPS_distributed_matrix_solve(K2,X,Z);
#elif (GETFEM_PARA_LEVEL==1 && GMM_USES_MUMPS)
      MUMPS_solve(K2,X,Z);
      //#elif (GETFEM_PARA_LEVEL==0 && GMM_USES_MUMPS)
      //MUMPS_solve(K2,X,Z);
#elif (GETFEM_PARA_LEVEL==0)
      //      SuperLU_solve(K2,X,Z,rcond);

      if ((time_order==1)||(t==Tinitial+dt)){ //time_order = 1 or first iterations with time_order = 2
	SLUsys2.solve(X, Z);
      }
      if ((time_order==2)&&(t>=Tinitial+2*dt)){
	//	SLUsys2.solve(X, gmm::scaled(gmm::sub_vector(Z,IP),1.5));
	SLUsys2.solve(X, gmm::scaled(Z,1.5));
      }


      // gmm::iteration iter(1E-8);
      // gmm::gmres(K2,X,Z,gmm::identity_matrix(),10,iter);
      // gmm::bicgstab(K2,X,Z,gmm::identity_matrix(),iter);
      // ?? gmm::bicgstab(K2,X,Z,gmm::diagonal_precond<sparse_matrix>(K2),iter);
      // if (noisy) cout << "condition number: " << 1.0/rcond << endl;
#endif
      gmm::copy(gmm::sub_vector(X, IP), Phi);
    }

    
    
    /*
      plain_vector toto(nbdof_p);
      for (unsigned i=0; i <= nbdof_p; ++i) {
      toto[i]=1.0;}
      cout << "Mean value of phi " << (gmm::vect_sp(Phi,toto))/nbdof_p<< endl;
    */

    if ((time_order==1)||(t==Tinitial+dt)){ //time_order = 1 or first iterations with time_order = 2
      gmm::mult(M, USTAR, USTARbis);
    }
    if ((time_order==2)&&(t>=Tinitial+2*dt)){
      gmm::mult(gmm::scaled(M,1.5), USTAR, USTARbis);
    }

    gmm::mult(gmm::transposed(B), gmm::scaled(Phi, -1.0), USTARbis, USTARbis); // -B^t*Phi + USTARbis -> USTARbis
#if GETFEM_PARA_LEVEL > 1
    MPI_SUM_VECTOR2(USTARbis);
#endif 

    gmm::copy(USTARbis, gmm::sub_vector(Y, I1)); //gmm::copy(M, gmm::sub_matrix(A1, I1));
 
  
	  //
	  // Solving the third linear system
	  //
    gmm::clear(YY);
    gmm::mult(gmm::transposed(Bbc), Phi, YY);
    
#if GETFEM_PARA_LEVEL > 1
    MPI_SUM_VECTOR2(YY);
#endif
    gmm::add(YY, gmm::sub_vector(Y, I1));
 
    {
      double rcond;
      plain_vector X(sizelsystem);
	  //plain_vector Xu(sizelsystem/2);
	  //plain_vector Xv(sizelsystem/2);
	  //plain_vector Yu(sizelsystem/2);
	  //plain_vector Yv(sizelsystem/2);
		
#if  (GETFEM_PARA_LEVEL > 1 && GETFEM_PARA_SOLVER == MUMPS_PARA_SOLVER)
                MUMPS_distributed_matrix_solve(A2,X,Y);
#elif (GETFEM_PARA_LEVEL==1 && GMM_USES_MUMPS)
      MUMPS_solve(A2,X,Y);
      //#elif (GETFEM_PARA_LEVEL==0 && GMM_USES_MUMPS)
      //MUMPS_solve(A2,X,Y);
#elif (GETFEM_PARA_LEVEL==0)
      //SuperLU_solve(A2, X, Y, rcond);
	  //SLUsys3.solve(X, Y);
		
	   
		//gmm::copy(gmm::sub_vector(Y,SUB_CT_Vu),Yu);
		//gmm::copy(gmm::sub_vector(Y,SUB_CT_Vv),Yv);
		
		//SuperLU_solve(A2u, Xu, Yu, rcond);
		//SuperLU_solve(A2u, Xv, Yv, rcond);
		//SLUsys3.solve(Xu, Yu);
		//SLUsys3.solve(Xv, Yv);
	
		SLUsys3.solve(gmm::sub_vector(X,SUB_CT_Vu),gmm::sub_vector(Y,SUB_CT_Vu));
		SLUsys3.solve(gmm::sub_vector(X,SUB_CT_Vv),gmm::sub_vector(Y,SUB_CT_Vv));
		
		
      //gmm::iteration iter(1E-8);
      //gmm::gmres(A2,X,Y,gmm::identity_matrix(),10,iter);
      //gmm::bicgstab(A2,X,Y,gmm::identity_matrix(),iter);
      // ?? gmm::bicgstab(A2,X,Y,gmm::diagonal_precond<sparse_matrix>(A2),iter);
      // if (noisy) cout << "condition number: " << 1.0/rcond << endl;
#endif		  

      //gmm:: clear(Un1);
      gmm::copy(gmm::sub_vector(X, I1), Un1);
		
		
		// Relation de compatibilitŽ int_domaine div(ustar)=0 //
		scalar_type delta;
		gmm :: resize(lambda,  nbdof_nonref);
		gmm:: clear(lambda);
		gmm::copy(gmm::sub_vector(X, I4),lambda);
		delta = 0.0;
		for (int ii = 0; ii <= nbdof_nonref; ii += 2) {
			delta = delta+lambda[ii];
		}
	    delta=delta/nbdof_nonref;
		for (int ii = 0; ii <= nbdof_nonref; ii += 2) {
			lambda[ii] = 1.0-delta;
		}
		cout << "delta" << delta << endl;
		gmm::add(lambda,gmm::sub_vector(X, I4));
		gmm::copy(gmm::sub_vector(X, I1), Un1);
		
		
		
		
      
//??//      if ((time_order==1)||(t==Tinitial+dt)){
//??//	gmm::add(gmm::scaled(gmm::sub_vector(X, I3C),1.0/dt), lambda);
//??//      }
//??//      if ((time_order==2)&&(t>=Tinitial+2*dt)){
//??//	gmm::add(gmm::scaled(gmm::sub_vector(X, I3C),1.5/dt), lambda);
//??//      }
      
    }


    //////////////////////////////////////////////////////////////////////////
    // INFORMATION POUR VERIFS 
    //cout << "I2   "<< gmm::sub_vector(Y, I2) << endl;

    // for (dal::bv_visitor i(dofon_nonref); !i.finished(); ++i) {
    //  bgeot :: base_node BN =  mf_u.point_of_basic_dof(i);
    //  cout << i  << "  " << BN  <<" Un1 =  " << Un1[i]  <<   endl;
    //}

 
    // gmm::sub_index SUB_CT_NDIR_tmp(ind_ct_ndir_tmp);
    //cout<<"NdirV ----- "<<gmm::sub_vector(Un1,SUB_CT_NDIR_tmp) << endl;
  
    //////////////////////////////////////////////////////////////////////////


   
    // gmm::add(USTAR, Un1);
    gmm::add(gmm::scaled(Phi, 1./dt), Pn0, Pn1);
    
    pdef->validate_solution(*this, t);
    cout << "Un1 - Un0 = " << gmm::vect_dist2(Un1, Un0) << endl;
    
    if (time_order==2){
      gmm::copy(Un0, Unm1);
    }
    
    gmm::copy(Un1, Un0); gmm::copy(Pn1, Pn0);
    
    if (N==3){
      plain_vector DU(mf_rhs.nb_dof() * N * N);
      plain_vector DIV(mf_rhs.nb_dof());
      compute_gradient(mf_u, mf_rhs, Un0, DU);
      for (unsigned i=0; i < mf_rhs.nb_dof(); ++i) {
	DIV[i] = DU[i*N*N] +DU[i*N*N+4] +DU[i*N*N+8];
      }
    if (PARAM.int_value("VTK_EXPORT")) {
      static int cnta=0;
      char sa[128]; sprintf(sa, "SolIcare/DivU%d.vtk", cnta++);
      getfem::vtk_export tata(sa, PARAM.int_value("VTK_EXPORT")==1);
      tata.exporting(mf_rhs); 
      tata.write_point_data(mf_rhs, DIV, "divergence");
    }

    }

    do_export(t);
 
    //
    // SORTIES : Coefficient de TrainŽ (Cd) et de Portance (Cl)
    //
 
	  
	  
	  if (N==2) {
		  std::vector<scalar_type> Cxn(1), Cxp(1), Cyn(1), Cyp(1);
		  getfem :: mesh_region mpioncylinder = mf_u.linked_mesh().get_mpi_sub_region(ON_CYLINDER_BOUNDARY_NUM);
		  //Cxn[0] = 0; Cxp[0] = 0; Cyn[0] = 0; Cyp[0] = 0;
		  getfem :: ClCd2D(Cxn,Cxp,Cyn,Cyp,mim,mf_u,mf_p,Un0,Pn1,mpioncylinder);
		  coeffTP<<t<<" "<<nu*Cxn[0]<<" "<<Cxp[0]<<" "<<nu*Cxn[0]+Cxp[0]<<" "<<nu*Cyn[0]<<" "<<Cyp[0]<<" "<<nu*Cyn[0]+Cyp[0]<<" "<<endl;
		  ptPartData << t << " " << Un1[ptPartU[2]] << " " << Un1[ptPartU[2] + 1] << " " << Pn1[ptPartP[2]] << endl;
	  }
	  
	  if (N==3) {
		  std::vector<scalar_type> Cxn(1), Cxp(1), Cyn(1), Cyp(1);
		  getfem :: mesh_region mpioncylinder = mf_u.linked_mesh().get_mpi_sub_region(ON_CYLINDER_BOUNDARY_NUM);
		  getfem :: ClCd3D(Cxn,Cxp,Cyn,Cyp,mim,mf_u,mf_p,Un0,Pn1,mpioncylinder);
		  coeffTP <<t<<" "<<Cxn[0]<<" "<<Cxp[0]<<" "<<Cxn[0]+Cxp[0]<<" "<<Cyn[0]<<" "<<Cyp[0]<<" "<<Cyn[0]+Cyp[0]<<" " << endl;
		  ptPartData <<t<<" "<<Un1[ptPartU[3]]<<" "<<Un1[ptPartU[3]+1]<<" " <<Un1[ptPartU[3]+2]<<" "<<Pn1[ptPartP[3]]<< endl;
		}
	  
  }
	  
  
  coeffTP.close();
  ptPartData.close();
}

  void navier_stokes_problem::do_export(scalar_type t) {
    /*dal :: bit_vector  ddd = mf_p.points_index();
      cout << ddd  << endl;*/ 
  if (!export_to_opendx) return;
  if (first_export) {
    mf_u.write_to_file("SolIcare/icare.mf_u", true);
    mf_p.write_to_file("SolIcare/icare.mf_p", true);
    mim.write_to_file("SolIcare/icare.mim",true);

     std::ofstream Uformat("Uformat.data");
   
     for(unsigned i=0;i<mf_u.nb_dof();++i){
       base_node G = mf_u.point_of_basic_dof(i);
       if (N==2){
 	Uformat << G[0] << "  " << G[1] << endl;
       }
     }
//     if (N==3){
//       for(unsigned i=0;i<mf_u.nb_dof();++i){
// 	base_node G = mf_u.point_of_basic_dof(i);
// 	Uformat << G[0] << "  " << G[1] << "  " << G[2] << endl;
//       }
//     }
     Uformat.close();
    
     std::ofstream Pformat("Pformat.data");
    
     for(unsigned i=0;i<mf_p.nb_dof();++i){
		 base_node G = mf_p.point_of_basic_dof(i);
       if (N==2){
 	Pformat << G[0] << "  " << G[1]  << endl;
       }
    }
    
//     if (N==3){
//       for(unsigned i=0;i<mf_p.nb_dof();++i){
// 	base_node G = mf_p.point_of_basic_dof(i);
// 	Pformat << G[0] << "  " << G[1] << "  " << G[2]  << endl;
//       }
//     }
     Pformat.close();
    
    
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
  if ((t >= t_export-dt/20.0)||(t>=99)) {
    t_export += dt_export;
    
    static int cnt = 0;
    char s[128]; sprintf(s, "SolIcare/icare.U%d", cnt++);
    gmm::vecsave(s, Un0);
    exp->write_point_data(mf_u, Un0);
    exp->serie_add_object("velocity");
    
    cout << "Saving Pressure, |p| = " << gmm::vect_norm2(Pn1) << "\n";
    exp->write_point_data(mf_p, Pn1);
    exp->serie_add_object("pressure");
  
    static int cntp=0;
    char sp[128]; sprintf(sp, "SolIcare/icare.P%d", cntp++);
    gmm::vecsave(sp, Pn0);
    
    if (PARAM.int_value("TIME_ORDER")==2){
      static int cntm1 = 0;
      char sm1[128]; sprintf(sm1, "SolIcare/icare.Um%d", cntm1++);
      gmm::vecsave(sm1, Unm1);
      exp->write_point_data(mf_u, Unm1);
      exp->serie_add_object("velocity");
    }


    //if (N == 2) {
    //plain_vector DU(mf_rhs.nb_dof() * N * N);
    //plain_vector Rot(mf_rhs.nb_dof());
    //compute_gradient(mf_u, mf_rhs, Un0, DU);
    //for (unsigned i=0; i < mf_rhs.nb_dof(); ++i) {
    //Rot[i] = DU[i*N*N + 3] - DU[i*N*N + 2];
    //if ((Rot[i]*Rot[i])<=1.5){Rot[i]=0;}
    //  }
    //  cout << "Saving Rot, |rot| = " << gmm::vect_norm2(Rot) << "\n";
    //  exp->write_point_data(mf_rhs, Rot);
    //  exp->serie_add_object("rot");
    //}
    if (PARAM.int_value("VTK_EXPORT")) {
      cout << "export to " << datafilename + "U.vtk" << "..\n";
      static int cnta=0;
      char sa[128]; sprintf(sa, "SolIcare/icareU%d.vtk", cnta++);
      getfem::vtk_export tata( sa,PARAM.int_value("VTK_EXPORT")==1);
      tata.exporting(mf_u); 
      tata.write_point_data(mf_u, Un0, "vitesse");
      
      if (PARAM.int_value("TIME_ORDER")==2){
	static int cntam1=0;
	char sam1[128]; sprintf(sam1, "SolIcare/icareUm%d.vtk", cntam1++);
	getfem::vtk_export tamtam( sam1,PARAM.int_value("VTK_EXPORT")==1);
	tamtam.exporting(mf_u); 
	tamtam.write_point_data(mf_u, Unm1, "vitesse");
      }
      
      static int cnte=0;
      char se[128]; sprintf(se, "SolIcare/icareP%d.vtk", cnte++);
      getfem::vtk_export tete( se,PARAM.int_value("VTK_EXPORT")==1);
      tete.exporting(mf_p);
      tete.write_point_data(mf_p, Pn0, "pression");

      //static int cnti=0;
      //char si[128]; sprintf(si, "SolIcare/icareRot%d.vtk", cnti++);
      //getfem::vtk_export titi( si, PARAM.int_value("VTK_EXPORT")==1);
      //titi.exporting(mf_rhs);
      //titi.write_point_data(mf_rhs, Rot, "rotationnel");
      }
    
    //else if (N == 3) {
      // plain_vector DU(mf_rhs.nb_dof() * N * N);
      //plain_vector RotX(mf_rhs.nb_dof());
      //plain_vector RotY(mf_rhs.nb_dof());
      //plain_vector RotZ(mf_rhs.nb_dof());
      //compute_gradient(mf_u, mf_rhs, Un0, DU);
      //for (unsigned i=0; i < mf_rhs.nb_dof(); ++i) {
      //	RotX[i] = DU[i*N*N + 7] - DU[i*N*N + 5];
      // RotY[i] = DU[i*N*N + 2] - DU[i*N*N + 6];
      // RotZ[i] = DU[i*N*N + 3] - DU[i*N*N + 1];
      //
      //static int cntt=0;
      //char st[128]; sprintf(st, "SolIcare/icareRotX%d.vtk", cntt++);
      //getfem::vtk_export titi( st, PARAM.int_value("VTK_EXPORT")==1);
      //titi.exporting(mf_rhs);
      //titi.write_point_data(mf_rhs, RotX, "rotationnelX");
      //static int cntu=0;
      //char su[128]; sprintf(su, "SolIcare/icareRotY%d.vtk", cntu++);
      //getfem::vtk_export tete( su, PARAM.int_value("VTK_EXPORT")==1);
      //tete.exporting(mf_rhs);
      //tete.write_point_data(mf_rhs, RotY, "rotationnelY");
      //static int cntv=0;
      //char sv[128]; sprintf(sv, "SolIcare/icareRotZ%d.vtk", cntv++);
      //getfem::vtk_export tyty( sv, PARAM.int_value("VTK_EXPORT")==1);
      //tyty.exporting(mf_rhs);
      //tyty.write_point_data(mf_rhs, RotZ, "rotationnelZ");
    //}
  }
  }

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {
    GETFEM_MPI_INIT(argc,argv);// For parallelized version
  // DAL_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {    
    navier_stokes_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh");
    p.solve();
   }
  GMM_STANDARD_CATCH_ERROR;
  GETFEM_MPI_FINALIZE;

  return 0; 
}
