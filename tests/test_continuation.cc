/*===========================================================================
 
 Copyright (C) 2011-2015 Yves Renard, Tomas Ligursky.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/

/**
   @file test_continuation.cc
   
   Continuation of solutions to the following state problem parameterized by
   lambda:
   -u'' +  u = lambda * exp(u) in (0, 1), u'(0) = u'(1) = 0.

   This program is used to check that getfem++ is working. This is also 
   a good example of use of Getfem++.
*/

#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_continuation.h" /* import continuation method */
#include "gmm/gmm_inoutput.h"

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

/* some Getfem++ types that we will be using */
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */

/* definition of a vector type; this one is built using the predefined types
   in Gmm++ */
typedef getfem::modeling_standard_plain_vector plain_vector;

/*
  structure for the state problem
*/
struct state_problem {

  getfem::mesh mesh;         /* the mesh */
  getfem::mesh_fem mf_u;     /* the FEM-method */
  getfem::mesh_im mim;       /* the integration method */
  scalar_type lambda;

  bgeot::md_param PARAM;

  bool cont(plain_vector &U);
  void init(void);
  state_problem(void) : mf_u(mesh), mim(mesh) {}
};


/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods.
 */
void state_problem::init(void) {
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";
  
  /* First step: build the mesh */
  std::vector<getfem::size_type> nsubdiv(1);
  nsubdiv[0] = PARAM.int_value("NX", "Number of space steps ");
  regular_unit_mesh(mesh, nsubdiv, bgeot::simplex_geotrans(1, 1));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mf_u.set_finite_element(pf_u);
  mim.set_integration_method(ppi);

}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool state_problem::cont(plain_vector &U) {
  
  //Define the model
  getfem::model model;
  model.add_fem_variable("u", mf_u);
  getfem::add_Laplacian_brick(model, mim, "u");
  std::string f = "u-lambda*exp(u)", dfdu = "1-lambda*exp(u)";
  model.add_fixed_size_data("lambda", 1);
  getfem::add_nonlinear_generic_assembly_brick(model, mim,
					       "(u-lambda*exp(u))*Test_u");


  // Initialise the continuation
  getfem::rmodel_plsolver_type ls =
    getfem::default_linear_solver<getfem::model_real_sparse_matrix,
                                  getfem::model_real_plain_vector>(model);
  size_type nb_dof = mf_u.nb_dof();
  scalar_type scfac = 1./ scalar_type(nb_dof);
  size_type nb_step = int(PARAM.int_value("NBSTEP",
					  "Number of continuation steps"));
  int singularities = (int) PARAM.int_value("SINGULARITIES", 
				      "Deal with singularities?");
  scalar_type  h_init = PARAM.real_value("H_INIT", "h_init"),
    h_max = PARAM.real_value("H_MAX", "h_max"),
    h_min = PARAM.real_value("H_MIN", "h_min"),
    h_inc = PARAM.real_value("H_INC", "h_inc"),
    h_dec = PARAM.real_value("H_DEC", "h_dec");
  size_type  maxit = PARAM.int_value("MAXITER", "maxit"),
    thrit = PARAM.int_value("THR_ITER", "thrit");
  scalar_type maxres = PARAM.real_value("RESIDUAL", "maxres"),
    maxdiff = PARAM.real_value("DIFFERENCE", "maxdiff"),
    mincos = PARAM.real_value("COS", "mincos"),
    maxres_solve = PARAM.real_value("RESIDUAL_SOLVE", "maxres_solve");
  int noisy = (int) PARAM.int_value("NOISY", "noisy");
  std::string datapath = PARAM.string_value("DATAPATH",
					    "Directory of data files");
  gmm::set_traces_level(noisy - 1);
  getfem::cont_struct_getfem_model
    S(model, "lambda", scfac, ls, h_init, h_max, h_min, h_inc, h_dec, maxit,
      thrit, maxres, maxdiff, mincos, maxres_solve, noisy, singularities);

  std::string bp_rootfilename = PARAM.string_value("BP_ROOTFILENAME").size()
    ? PARAM.string_value("BP_ROOTFILENAME") : "";
  scalar_type direction = PARAM.real_value("DIRECTION", "Initial direction"),
    h, T_lambda;
  plain_vector T_U(U), Y(nb_dof + 1);
  
  if (bp_rootfilename.size() > 0) {
    gmm::vecload(datapath + bp_rootfilename + ".Y", Y);
    gmm::copy(gmm::sub_vector(Y, gmm::sub_interval(0, nb_dof)), U);
    lambda = Y[nb_dof];
    char s[100];
    sprintf(s, ".T_Y%d", (int) PARAM.int_value("IND_BRANCH", "Branch"));
    gmm::vecload(datapath + bp_rootfilename + s, Y);
    gmm::copy(gmm::scaled(gmm::sub_vector(Y, gmm::sub_interval(0, nb_dof)),
			  direction), T_U);
    T_lambda = direction * Y[nb_dof];
    h = S.h_init();
  } else {
    lambda = PARAM.real_value("LAMBDA0", "lambda0");
    model.set_real_variable("lambda")[0] = lambda;   
    if (noisy > 0) cout << "Starting computing an initial point" << endl;
    gmm::iteration iter(maxres_solve, noisy - 1, 40000);
    getfem::standard_solve(model, iter);
    gmm::copy(model.real_variable("u"), U);
    T_lambda = direction;

    S.init_Moore_Penrose_continuation(U, lambda, T_U, T_lambda, h);
  }

//   cout << "U = " << U << endl;
//   cout << "lambda - u * exp(-u) = " << lambda - U[0] * exp(-U[0]) << endl;

  // Continuation
  std::string sing_label;
  char s1[100], s2[100];
  std::vector<std::string> sing_out;
  for (size_type step = 0; step < nb_step; ++step) {
    cout << endl << "Beginning of step " << step + 1 << endl;
   
    S.Moore_Penrose_continuation(U, lambda, T_U, T_lambda, h);
    if (h == 0) break;

//     cout << "U = " << U << endl;
//     cout << "lambda = " << lambda << endl;
//     cout << "lambda - U[0] * exp(-U[0]) = "
// 	 << lambda - U[0] * exp(-U[0]) << endl;

    sing_label = S.get_sing_label();
    if (sing_label.size() > 0) {
      if (sing_label == "limit point")
	sprintf(s1, "Step %lu: %s", step + 1, sing_label.c_str());
      else if (sing_label == "smooth bifurcation point") {
	gmm::copy(S.get_x_sing(),
		  gmm::sub_vector(Y, gmm::sub_interval(0, nb_dof)));
	Y[nb_dof] = S.get_gamma_sing();
	sprintf(s1, "continuation_step_%lu", step + 1);
	gmm::vecsave(datapath + s1 + "_bp.Y", Y);
	for (size_type i = 0; i < S.nb_tangent_sing(); i++) {
	  gmm::copy(S.get_tx_sing(i),
		    gmm::sub_vector(Y, gmm::sub_interval(0, nb_dof)));
	  Y[nb_dof] = S.get_tgamma_sing(i);
	  sprintf(s2, "_bp.T_Y%lu", i + 1);
	  gmm::vecsave(datapath + s1 + s2, Y);
	}
	sprintf(s1, "Step %lu: %s, %u branch(es) located", step + 1,
		sing_label.c_str(), (unsigned int) S.nb_tangent_sing());
      }
      sing_out.push_back(s1);
    }
    cout << "End of Step nº " << step + 1 << " / " << nb_step << endl;
  }

  if (sing_out.size() > 0) {
    cout << endl
	 << "------------------------------" 
	 << endl
	 << "   Detected singular points"
	 << endl
	 << "------------------------------"
	 << endl;
    for (size_type i = 0; i < sing_out.size(); i++)
      cout << sing_out[i] << endl;
    cout << endl;
  }

  return (h > 0);
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {    
    state_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    plain_vector U(p.mf_u.nb_dof());
    p.cont(U); 
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
