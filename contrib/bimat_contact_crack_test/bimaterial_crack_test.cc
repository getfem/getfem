// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2002-2006 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================
  
/**
 * Linear Elastostatic problem with a crack.
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
*/

#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_derivatives.h>
#include <getfem_regular_meshes.h>
#include <getfem_model_solvers.h>
#include <getfem_mesh_im_level_set.h>
#include <getfem_mesh_fem_level_set.h>
#include <getfem_mesh_fem_product.h>
#include <getfem_mesh_fem_global_function.h>
#include <getfem_spider_fem.h>
#include <getfem_mesh_fem_sum.h>
#include <gmm.h>
#include <getfem_error_estimate.h>

#include <getfem_interpolated_fem.h>


/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are builtmayavi -d crack.vtk -f WarpVector -m BandedSurfaceMap -m Outline

 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;


base_small_vector sol_f(const base_node &x) {
  int N = x.size();
  base_small_vector res(N); res[N-1] = x[N-1];
  return res;
}

/**************************************************************************/
/*  Structure for the crack problem.                                      */
/**************************************************************************/

struct crack_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM_down = 1, NEUMANN_BOUNDARY_NUM_up=2, NEUMANN_BOUNDARY_NUM_right=3};
  getfem::mesh mesh;  /* the mesh */
  getfem::mesh_level_set mls;       /* the integration methods.              */
  getfem::mesh_im_level_set mim;    /* the integration methods.              */
  getfem::mesh_fem mf_pre_u;
  getfem::mesh_fem mf_mult;
  getfem::mesh_fem_level_set mfls_u; 
  getfem::mesh_fem_sum mf_u_sum;
  
  getfem::mesh_fem& mf_u() { return mf_u_sum; }
  
  scalar_type lambda, mu;    /* Lame coefficients.                */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  
  
  double lx,ly;             /* size of the mesh */
  int bimaterial;           /* For bimaterial interface fracture */
  bool all_dirichlet;
  double F11,F12,F21,F22,F31,F32,F41,F42;       /* NEUMANN forces */
  double lambda_up, lambda_down, mu_up, mu_down;  /*Lame coeff for bimaterial case*/
  getfem::level_set ls;      /* The two level sets defining the crack.       */
  
  base_small_vector translation;

  scalar_type residual;       /* max residual for the iterative solvers        */
  scalar_type conv_max;
  unsigned dir_with_mult;
  
  std::string datafilename;
  ftool::md_param PARAM;
  
  bool adapted_refine, solve(plain_vector &U);
  
  void init(void);
  crack_problem(void) : mls(mesh), mim(mls), mf_pre_u(mesh), mf_mult(mesh),
			mfls_u(mls, mf_pre_u),
			
			mf_u_sum(mesh), mf_rhs(mesh), 

			ls(mesh, 1, true) {}

};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void crack_problem::init(void) {
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");

  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
					 "Name of simplex integration method");
  std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION");

  
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";


  /* First step : build the mesh */
  bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Nomber of space steps "));
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			    PARAM.int_value("MESH_NOISED") != 0);
  
  
  lx = PARAM.real_value("LX","length x'ox");
  ly = PARAM.real_value("LY","length y'oy");
  
  bgeot::base_matrix M(2,2);
  M(0,0) = lx;   
  M(1,1) = ly;
  mesh.transformation(M);
  
  base_small_vector tt(N); tt[0] = tt[1] = -(lx/2.);
  mesh.translation(tt); 

  conv_max = PARAM.int_value("CONV_MAX","Maximal number of convexes in the mesh");
  adapted_refine = PARAM.int_value("ADAPTED_REFINE","Adapted Refinement");
  
  scalar_type refinement_radius;
  refinement_radius = PARAM.real_value("REFINEMENT_RADIUS","Refinement Radius");
  size_type refinement_process;
  refinement_process = PARAM.int_value("REFINEMENT_PROCESS","Refinement process");
  dal::bit_vector conv_to_refine;
  size_type ref = 0;  
  if (refinement_radius > 0){
    while(ref <= refinement_process){
      conv_to_refine.clear();
      for(dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i){
      for(size_type j=0; j < 3; ++j)
	if(fabs(mesh.points()[mesh.ind_points_of_convex(i)[j]][0])<refinement_radius 
	   && fabs(mesh.points()[mesh.ind_points_of_convex(i)[j]][1])<refinement_radius){
	  conv_to_refine.add(i);
	}
      }
      mesh.Bank_refine(conv_to_refine);
      ref = ref + 1;
      refinement_radius = refinement_radius/2.;
      if(refinement_radius > 1e-16)
	cout<<"refining process step " << ref << "... refining "<< conv_to_refine.size() <<" convexes..." << endl ; 
    }
  cout<<"refining process complete." << endl ;
  }
  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL");
  if (residual == 0.) residual = 1e-10;
  
  bimaterial = PARAM.int_value("BIMATERIAL", "Bimaterial interface crack");
  all_dirichlet = PARAM.int_value("all_dirichlet", "Dirichlet condition");
  F11 = PARAM.real_value("F11","F11");
  F12 = PARAM.real_value("F12","F12");
  F21 = PARAM.real_value("F21","F21");
  F22 = PARAM.real_value("F22","F22");
  F31 = PARAM.real_value("F31","F31");
  F32 = PARAM.real_value("F32","F32");
  F41 = PARAM.real_value("F41","F41");
  F42 = PARAM.real_value("F42","F42");
 
  if (bimaterial == 1){
    mu = PARAM.real_value("MU", "Lame coefficient mu"); 
    lambda_up = PARAM.int_value("LAMBDA_UP", "Lame Coef");
    lambda_down = PARAM.int_value("LAMBDA_DOWN", "Lame Coef");
    lambda = PARAM.real_value("LAMBDA", "Lame coefficient lambda");
    mu_up = PARAM.int_value("MU_UP", "Lame Coef");
    mu_down = PARAM.int_value("MU_DOWN", "Lame Coef");
    
  }
  else{
    
    mu = PARAM.real_value("MU", "Lame coefficient mu");
    lambda = PARAM.real_value("LAMBDA", "Lame coefficient lambda");
  }
  

  mf_u().set_qdim(N);


  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
 
 
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);
  getfem::pintegration_method simp_ppi = 
    getfem::int_method_descriptor(SIMPLEX_INTEGRATION);
  getfem::pintegration_method sing_ppi = (SINGULAR_INTEGRATION.size() ? getfem::int_method_descriptor(SINGULAR_INTEGRATION) : 0);
  
  mim.set_integration_method(mesh.convex_index(), ppi);
  mls.add_level_set(ls);


  mim.set_simplex_im(simp_ppi, sing_ppi);
  mf_pre_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_mult.set_finite_element(mesh.convex_index(), pf_u);
  mf_mult.set_qdim(N);


  dir_with_mult = PARAM.int_value("DIRICHLET_VERSINO");
 
  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
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
  
  /* set boundary conditions
   * (Neuman on the upper face, Dirichlet elsewhere) */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
    
    base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);
    if(all_dirichlet) 
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
    else {
      if (un[0]  > 1.0E-7 ) { // new Neumann face
	mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
      } else {
	if (un[1]  > 1.0E-7 ) {
	  //cout << "normal = " << un << endl;
	  mesh.region(NEUMANN_BOUNDARY_NUM_up).add(i.cv(), i.f());
	}
	else {
	  if (un[1]  < -1.0E-7 ) {
	    //cout << "normal = " << un << endl;
	    mesh.region(NEUMANN_BOUNDARY_NUM_down).add(i.cv(), i.f());
	  }
	  else {
	    
	    if (un[0]  < -1.0E-7 ) {
	      
	      if(mesh.points_of_convex(i.cv())[mesh.structure_of_convex(i.cv())->ind_points_of_face(i.f())[0]][1] > 1.0E-16 ||  mesh.points_of_convex(i.cv())[mesh.structure_of_convex(i.cv())->ind_points_of_face(i.f())[1]][1] > 1.0E-16) {
		
		// cout << "normal = " << un << endl;
		mesh.region(NEUMANN_BOUNDARY_NUM_right).add(i.cv(), i.f());
	      }
	      else
		mesh.region(NEUMANN_BOUNDARY_NUM_right+1).add(i.cv(), i.f());
	    }
	  }
	}
      }
    }
  }

}


base_small_vector ls_function(const base_node P, int num = 0) {
  
  scalar_type x = P[0], y = P[1];
  base_small_vector res(2);
  switch (num) {
    case 0: {
      res[0] =  y;
      res[1] =  x;
    } break;
    case 1: {
      res[0] = gmm::vect_dist2(P, base_node(0.5, 0.)) - .25;
      res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.27;
    } break;
    case 2: {
      res[0] = x - 0.25;
      res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
    } break;
    default: assert(0);
  }
  return res;
}

bool crack_problem::solve(plain_vector &U) {
  size_type N = mesh.dim();
  
  
  // Linearized elasticity brick.
  dal::bit_vector conv_to_refine;
  bool iteration;
  
  do {
    cout << "Number of convexes : "<<  mesh.convex_index().size() <<endl;
    size_type nb_dof_rhs = mf_rhs.nb_dof();
    ls.reinit();  
    cout << "ls.get_mesh_fem().nb_dof() = " << ls.get_mesh_fem().nb_dof() << "\n";
    for (size_type d = 0; d < ls.get_mesh_fem().nb_dof(); ++d) {
      ls.values(0)[d] = ls_function(ls.get_mesh_fem().point_of_dof(d), 0)[0];
      ls.values(1)[d] = ls_function(ls.get_mesh_fem().point_of_dof(d), 0)[1];
    }
    ls.touch();
    mls.adapt();
    mim.adapt();
    mfls_u.adapt();
 
    mf_u_sum.set_mesh_fems(mfls_u);
  
    U.resize(mf_u().nb_dof());

    conv_to_refine.clear();
    getfem::mdbrick_isotropic_linearized_elasticity<>
      ELAS(mim, mf_u(), lambda, mu);
     
     
    if(bimaterial == 1){
      cout<<"______________________________________________________________________________"<<endl;
      cout<<"CASE OF BIMATERIAL CRACK  with lambda_up = "<<lambda_up<<" and lambda_down = "<<lambda_down<<endl;
      cout<<"______________________________________________________________________________"<<endl;
      std::vector<float> bi_lambda(ELAS.lambda().mf().nb_dof());
      std::vector<float> bi_mu(ELAS.lambda().mf().nb_dof());
    
      cout<<"ELAS.lambda().mf().nb_dof()==="<<ELAS.lambda().mf().nb_dof()<<endl;
    
      for (size_type ite = 0; ite < ELAS.lambda().mf().nb_dof();ite++) {
	if (ELAS.lambda().mf().point_of_dof(ite)[1] > 0){
	  bi_lambda[ite] = lambda_up;
	  bi_mu[ite] = mu_up;
	}
	else{
	  bi_lambda[ite] = lambda_down;
	  bi_mu[ite] = mu_down;
	}
      } 
      //cout<<"bi_lambda.size() = "<<bi_lambda.size()<<endl;
      // cout<<"ELAS.lambda().mf().nb_dof()==="<<ELAS.lambda().mf().nb_dof()<<endl;
    
      ELAS.lambda().set(bi_lambda);
      ELAS.mu().set(bi_mu);
    }
  
 
    // Defining the volumic source term.
    plain_vector F(nb_dof_rhs * N);
    for (size_type i = 0; i < nb_dof_rhs; ++i)
      gmm::copy(sol_f(mf_rhs.point_of_dof(i)),
		gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  
    // Volumic source term brick.
    getfem::mdbrick_source_term<> VOL_F(ELAS, mf_rhs, F);

    // Defining the Neumann condition right hand side.
    
    // Neumann condition brick.
  
    getfem::mdbrick_abstract<> *pNEUMANN;

    gmm::clear(F);
    for(size_type i = 0; i<F.size(); i=i+2) 
      {F[i] = F41; F[i+1] = F42;}
  
    getfem::mdbrick_source_term<> NEUMANN_down(VOL_F, mf_rhs, F,NEUMANN_BOUNDARY_NUM_down);
  
    gmm::clear(F);
    for(size_type i = 0; i<F.size(); i=i+2) 
      {F[i] = F21; F[i+1] = F22;}
  
    getfem::mdbrick_source_term<> NEUMANN_right_up(NEUMANN_down, mf_rhs, F,NEUMANN_BOUNDARY_NUM_right);
  
  
    gmm::clear(F);
    for(size_type i = 0; i<F.size(); i=i+2) 
      {F[i] = F31; F[i+1] = F32;}
  
    getfem::mdbrick_source_term<> NEUMANN_right_down(NEUMANN_right_up, mf_rhs, F,NEUMANN_BOUNDARY_NUM_right+1);
  
    gmm::clear(F);
    for(size_type i = 0; i<F.size(); i=i+2) 
      {F[i] = F11; F[i+1] = F12;}
  
    getfem::mdbrick_source_term<> NEUMANN_up(NEUMANN_right_down, mf_rhs, F,NEUMANN_BOUNDARY_NUM_up);
  
    if (all_dirichlet)
      pNEUMANN = & VOL_F; 
    else 
      pNEUMANN = & NEUMANN_up; 
    
    //toto_solution toto(mf_rhs.linked_mesh()); toto.init();
    //assert(toto.mf.nb_dof() == 1);
  
    // Dirichlet condition brick.
    getfem::mdbrick_Dirichlet<> final_model(*pNEUMANN, DIRICHLET_BOUNDARY_NUM,
					    mf_mult);
    if(all_dirichlet){

    } else {

    }
    final_model.set_constraints_type(getfem::constraints_type(dir_with_mult));
  
    // Generic solve.
    cout << "Total number of variables : " << final_model.nb_dof() << endl;
    getfem::standard_model_state MS(final_model);
    gmm::iteration iter(residual, 1, 40000);
  
    getfem::standard_solve(MS, final_model, iter);
  
    // Solution extraction
    gmm::copy(ELAS.get_solution(MS), U);
    iteration = iter.converged();  

    // Adapted Refinement (suivant une erreur a posteriori)
    if (adapted_refine) {
      if (mesh.convex_index().size() < conv_max){
      plain_vector ERR(mesh.convex_index().last_true()+1);
      getfem::error_estimate(mim, mf_u(), U, ERR);
    
      cout << "max = " << gmm::vect_norminf(ERR) << endl;
      scalar_type threshold = 0.01, min_ = 1e18;
      conv_to_refine.clear();
      for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i) {
	if (ERR[i] > threshold) conv_to_refine.add(i);
	min_ = std::min(min_, ERR[i]);
      }
      cout << "min = " << min_ << endl;
      cout << "Refining " <<  conv_to_refine.card() << " convexes..."<< endl;  
      mesh.Bank_refine(conv_to_refine);
    }    
      else
	break;
    }
  } while(adapted_refine && conv_to_refine.size() > 0);
  mesh.write_to_file(datafilename + ".meshh");
  cout << "Refining process complete. The mesh contains now " <<  mesh.convex_index().size() << " convexes "<<endl;
  
  dal::bit_vector blocked_dof = mf_u().dof_on_set(5);
  getfem::mesh_fem mf_printed(mesh, N);
  std::string FEM_DISC = PARAM.string_value("FEM_DISC","fem disc ");
  mf_printed.set_finite_element(mesh.convex_index(),
				getfem::fem_descriptor(FEM_DISC));
  plain_vector W(mf_printed.nb_dof());
  getfem::interpolation(mf_u(), mf_printed, U, W);
  mf_printed.write_to_file(datafilename + ".meshfem", true);
  gmm::vecsave(datafilename + ".U", W);

  return (iteration);
}

  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  DAL_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  //getfem::getfem_mesh_level_set_noisy();


  try {
    crack_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh");
    plain_vector U(p.mf_u().nb_dof());
     if (!p.solve(U)) DAL_THROW(dal::failure_error,"Solve has failed");
   
    {
      getfem::mesh mcut;
      p.mls.global_cut_mesh(mcut);
      unsigned Q = p.mf_u().get_qdim();
      getfem::mesh_fem mf(mcut, Q);
      mf.set_classical_discontinuous_finite_element(2, 0.001);
      // mf.set_finite_element
      //	(getfem::fem_descriptor("FEM_PK_DISCONTINUOUS(2, 2, 0.0001)"));
      plain_vector V(mf.nb_dof());

      getfem::interpolation(p.mf_u(), mf, U, V);

      getfem::stored_mesh_slice sl;
      getfem::mesh mcut_refined;

      unsigned NX = p.PARAM.int_value("NX"), nn;
      if (NX < 6) nn = 24;
      else if (NX < 12) nn = 8;
      else if (NX < 30) nn = 3;
      else nn = 1;

      /* choose an adequate slice refinement based on the distance to the crack tip */
      std::vector<bgeot::short_type> nrefine(mcut.convex_index().last_true()+1);
      for (dal::bv_visitor cv(mcut.convex_index()); !cv.finished(); ++cv) {
	scalar_type dmin=0, d;
	base_node Pmin,P;
	for (unsigned i=0; i < mcut.nb_points_of_convex(cv); ++i) {
	  P = mcut.points_of_convex(cv)[i];
	  d = gmm::vect_norm2(ls_function(P));
	  if (d < dmin || i == 0) { dmin = d; Pmin = P; }
	}

	if (dmin < 1e-5)
	  nrefine[cv] = nn*8;
	else if (dmin < .1) 
	  nrefine[cv] = nn*2;
	else nrefine[cv] = nn;
	if (dmin < .01)
	  cout << "cv: "<< cv << ", dmin = " << dmin << "Pmin=" << Pmin << " " << nrefine[cv] << "\n";
      }

      {
	getfem::mesh_slicer slicer(mcut); 
	getfem::slicer_build_mesh bmesh(mcut_refined);
	slicer.push_back_action(bmesh);
	slicer.exec(nrefine, getfem::mesh_region::all_convexes());
      }
      /*
      sl.build(mcut, 
      getfem::slicer_build_mesh(mcut_refined), nrefine);*/

      getfem::mesh_im mim_refined(mcut_refined); 
      mim_refined.set_integration_method(getfem::int_method_descriptor
					 ("IM_TRIANGLE(6)"));

      getfem::mesh_fem mf_refined(mcut_refined, Q);
      mf_refined.set_classical_discontinuous_finite_element(2, 0.0001);
      plain_vector W(mf_refined.nb_dof());

      getfem::interpolation(p.mf_u(), mf_refined, U, W);



      if (p.PARAM.int_value("VTK_EXPORT")) {
	getfem::mesh_fem mf_refined_vm(mcut_refined, 1);
	mf_refined_vm.set_classical_discontinuous_finite_element(1, 0.0001);
	cerr << "mf_refined_vm.nb_dof=" << mf_refined_vm.nb_dof() << "\n";
	plain_vector VM(mf_refined_vm.nb_dof());

	cout << "computing von mises\n";
	getfem::interpolation_von_mises(mf_refined, mf_refined_vm, W, VM);

	plain_vector D(mf_refined_vm.nb_dof() * Q), 
	  DN(mf_refined_vm.nb_dof());
	

	cout << "export to " << p.datafilename + ".vtk" << "..\n";
	getfem::vtk_export exp(p.datafilename + ".vtk",
			       p.PARAM.int_value("VTK_EXPORT")==1);
	exp.exporting(mf_refined); 
	//exp.write_point_data(mf_refined_vm, DN, "error");
	exp.write_point_data(mf_refined_vm, VM, "von mises stress");
	exp.write_point_data(mf_refined, W, "elastostatic_displacement");
      


	cout << "export done, you can view the data file with (for example)\n"
	  "mayavi -d " << p.datafilename << ".vtk -f "
	  "WarpVector -m BandedSurfaceMap -m Outline\n";
      }
    }

  }
  DAL_STANDARD_CATCH_ERROR;

  return 0; 
}



