// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2008 Yves Renard, Julien Pommier.
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
 * Linear Elastostatic problem with a crack.
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
 */

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_nonlinear_elasticity.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_mesh_fem_product.h"
#include "getfem/getfem_mesh_fem_global_function.h"
#include "getfem/getfem_spider_fem.h"
#include "getfem/getfem_mesh_fem_sum.h"
#include "getfem/getfem_superlu.h"
#include "getfem_nonlinear_elastoptim.h"
#include "gmm/gmm.h"
#include "gmm/gmm_inoutput.h"

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::base_vector;
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::dim_type; 
using bgeot::short_type;
using bgeot::base_vector;
using bgeot::base_matrix; /* small dense matrix. */



/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;




/**************************************************************************/
/*                                                                        */
/* Brick definition.                                                      */
/*                                                                        */
/**************************************************************************/

namespace getfem {

struct nonlinear_elasticity_optim_brick : public virtual_brick {
  
  const abstract_hyperelastic_law &AHL;
  const getfem::level_set &ls;
  
  virtual void asm_real_tangent_terms(const model &md, size_type /* ib */,
				      const model::varnamelist &vl,
				      const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {
      GMM_ASSERT1(mims.size() == 1,
		  "Nonlinear elasticity brick need a single mesh_im");
      GMM_ASSERT1(vl.size() == 4,
		  "Nonlinear elasticity brick need a single variable");
      GMM_ASSERT1(dl.size() == 1,
		  "Wrong number of data for nonlinear elasticity brick, "
                  << dl.size() << " should be 1 (vector).");
      GMM_ASSERT1(matl.size() == 7,  "Wrong number of terms for nonlinear "
		  "elasticity brick");

      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const mesh_fem &mf_u = *(md.pmesh_fem_of_variable(vl[0]));

      const model_real_plain_vector &p = md.real_variable(vl[1]);
      const mesh_fem &mf_p = *(md.pmesh_fem_of_variable(vl[1]));


      const model_real_plain_vector &alpha = md.real_variable(vl[2]);
      const model_real_plain_vector &beta = md.real_variable(vl[3]);

      const mesh_fem *mf_params = md.pmesh_fem_of_variable(dl[0]);
      const model_real_plain_vector &params = md.real_variable(dl[0]);
      const mesh_im &mim = *mims[0];

      size_type sl = gmm::vect_size(params);
      if (mf_params) sl = sl * mf_params->get_qdim() / mf_params->nb_dof();
      GMM_ASSERT1(sl == AHL.nb_params(), "Wrong number of coefficients for the "
		  "nonlinear constitutive elastic law");

      mesh_region rg(region);
      mf_u.linked_mesh().intersect_with_mpi_region(rg);

      
      std::vector<scalar_type> P_ls(gmm::vect_size(p));
      
      dim_type d = mf_u.linked_mesh().dim();
      dal::bit_vector p_enriched_dof;
      for (dal::bv_visitor cv(mf_p.linked_mesh().convex_index());
	   !cv.finished(); ++cv) {
	pfem pf = mf_p.fem_of_element(cv);
	for (size_type j = 0; j < pf->nb_dof(cv); ++j)
	  if (pf->dof_types()[j] == global_dof(d)) {
	    size_type dof = mf_p.ind_basic_dof_of_element(cv)[j];
	    p_enriched_dof.add(dof);
	    P_ls[dof] = p[dof];
	    
	  }
      }
      
      
      std::vector<scalar_type> U_ls(gmm::vect_size(u));
      
      dal::bit_vector u_enriched_dof;
      for (dal::bv_visitor cv(mf_u.linked_mesh().convex_index());
	   !cv.finished(); ++cv) {
	pfem pf = mf_u.fem_of_element(cv);
	for (size_type j = 0; j< pf->nb_dof(cv); ++j)
	  if (pf->dof_types()[j] == global_dof(d)) {
	    for (size_type k = 0; k < d; ++k) {
	      size_type dof = mf_u.ind_basic_dof_of_element(cv)[j*d+k];
	      u_enriched_dof.add(dof);
	      U_ls[dof] = u[dof];
	    }
	  }
      }
      

      if (version & model::BUILD_MATRIX) {
	std::vector<scalar_type> V(gmm::vect_size(u));
	GMM_TRACE2("Nonlinear elasticity stiffness matrix assembly");
	asm_nonlinear_elasticity_optim_tangent_matrix
	  (V, mim, mf_u, u, U_ls, u_enriched_dof,
	   mf_p, p, P_ls, p_enriched_dof, alpha, beta, mf_params,
	   params, AHL, ls,rg);
	gmm::copy(gmm::row_vector(V), matl[0]);
      }


      if (version & model::BUILD_RHS) {
	asm_nonlinear_elasticity_optim_rhs(vecl[0], vecl[3], mim, 
					   mf_u, u, U_ls, u_enriched_dof,
					   mf_p, p, P_ls, p_enriched_dof,
					   alpha, beta,
					   mf_params, params, AHL, ls, rg);
	gmm::scale(vecl[0], scalar_type(-1));
      }

    }


    nonlinear_elasticity_optim_brick(const abstract_hyperelastic_law &AHL_,
				     const getfem::level_set &ls_)
      : AHL(AHL_), ls(ls_) {
      set_flags("Nonlinear elasticity brick", false /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
		true /* is real */, false /* is complex */);
    }

  };
  
  //=========================================================================
  //  Add a nonlinear elasticity brick.  
  //=========================================================================

  size_type add_nonlinear_elasticity_optim_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &varname_p,  const std::string &varname_alpha,
   const std::string &varname_beta,
   const abstract_hyperelastic_law &AHL, const std::string &dataname_law,
   const getfem::level_set &ls, size_type region = size_type(-1)) {
    pbrick pbr = new nonlinear_elasticity_optim_brick(AHL, ls);

    model::termlist tl;
    tl.push_back(model::term_description(varname_alpha, varname_u, true));
    tl.push_back(model::term_description(varname_alpha, varname_p, true));
    tl.push_back(model::term_description(varname_alpha, varname_alpha, true));
    tl.push_back(model::term_description(varname_beta, varname_u, true));
    tl.push_back(model::term_description(varname_beta, varname_p, true));
    tl.push_back(model::term_description(varname_beta, varname_alpha, true));
    tl.push_back(model::term_description(varname_beta, varname_beta, true));
    model::varnamelist dl(1, dataname_law);
    model::varnamelist vl(1, varname_u);
    vl.push_back(varname_p);
    vl.push_back(varname_alpha);
    vl.push_back(varname_beta);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1,&mim), region);
  }







}






/**************************************************************************/
/*                                                                        */
/*  Structure for the crack problem.                                      */
/*                                                                        */
/**************************************************************************/

struct crack_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN1_BOUNDARY_NUM = 1, NEUMANN2_BOUNDARY_NUM=2, NEUMANN3_BOUNDARY_NUM=3, NEUMANN4_BOUNDARY_NUM=4, MORTAR_BOUNDARY_IN=42, MORTAR_BOUNDARY_OUT=43};
  getfem::mesh mesh;  /* the mesh */
  getfem::level_set ls;      /* The two level sets defining the crack.       */
  getfem::mesh_level_set mls;       /* the integration methods for cutted element.    */
  getfem::mesh_im_level_set mim;    /* the integration methods.              */
  getfem::mesh_fem mf_pre_u, mf_pre_mortar;
  getfem::mesh_fem mf_mult, mf_mult_p;
  getfem::mesh_fem_level_set mfls_u,mfls_mortar;
  getfem::mesh_fem_global_function mf_sing_u;
  
  getfem::mesh_fem mf_partition_of_unity;
  getfem::mesh_fem_product mf_product;
  getfem::mesh_fem_sum mf_u_sum;
  
  getfem::mesh_fem mf_pre_p; /* mesh_fem for the pressure for mixed form     */
  getfem::mesh_fem_level_set mfls_p;   /* mesh_fem for the pressure enriched with H.   */
  getfem::mesh_fem_global_function mf_sing_p;
  getfem::mesh_fem_product mf_product_p;
  getfem::mesh_fem_sum mf_p_sum;

  scalar_type pr1, pr2, pr3, AMP_LOAD_X, AMP_LOAD_Y, nb_step;   /* elastic coefficients,Amplitude of load, number of step loadings*/

  base_small_vector cracktip;



  struct spider_param {
    getfem::spider_fem *fem;
    scalar_type theta0;
    scalar_type radius;
    unsigned Nr;
    unsigned Ntheta;
    int K;
    int bimat_enrichment;
    scalar_type epsilon;
  };

  spider_param spider;

  getfem::mesh_fem mf_us;
  //getfem::mesh_fem mf_p;

  getfem::mesh_fem& mf_u() { return mf_u_sum; }
  //getfem::mesh_fem& mfls_p() { return mf_p_sum; }
  getfem::mesh_fem& mf_pe() { return mf_p_sum; }
  // getfem::mesh_fem& mf_u() { return mf_us; }
  
  //scalar_type mu;            /* Lame coefficients.                */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  
  
  int bimaterial;           /* For bimaterial interface fracture */
  //double mu_up, mu_down;  /*Lame coeff for bimaterial case*/
 
  scalar_type residual;      /* max residual for the iterative solvers      */
  bool mixed_pressure;
  unsigned dir_with_mult;
  scalar_type cutoff_radius, cutoff_radius1, cutoff_radius0, enr_area_radius;
  
  size_type cutoff_func;

  typedef enum { NO_ENRICHMENT=0, 
		 FIXED_ZONE=1, 
		 GLOBAL_WITH_MORTAR=2,
		 GLOBAL_WITH_CUTOFF=3,
		 SPIDER_FEM_ALONE=4,
		 SPIDER_FEM_ENRICHMENT=5 } enrichment_option_enum;
  enrichment_option_enum enrichment_option;

  std::string datafilename;
  
  int reference_test;
  std::string GLOBAL_FUNCTION_MF, GLOBAL_FUNCTION_U,  GLOBAL_FUNCTION_P;

  bgeot::md_param PARAM;

  bool solve(plain_vector &U, plain_vector &P);
  void init(void);
  crack_problem(void) : ls(mesh, 1, true), mls(mesh), mim(mls), 
			mf_pre_u(mesh), mf_pre_mortar(mesh), mf_mult(mesh), mf_mult_p(mesh),
			mfls_u(mls, mf_pre_u), mfls_mortar(mls, mf_pre_mortar),	
			mf_sing_u(mesh),
			mf_partition_of_unity(mesh),
                        mf_product(mf_partition_of_unity, mf_sing_u),
			mf_u_sum(mesh), mf_pre_p(mesh), mfls_p(mls, mf_pre_p), 
                        mf_sing_p(mesh), mf_product_p(mf_partition_of_unity, mf_sing_p), 
                	mf_p_sum(mesh),
			mf_us(mesh), /*mf_pe(mesh),*/ mf_rhs(mesh)
 {}

};

std::string name_of_dof(getfem::pdof_description dof) {
  char s[200];
  sprintf(s, "UnknownDof[%p]", (void*)dof);
  for (dim_type d = 0; d < 4; ++d) {
    if (dof == getfem::lagrange_dof(d)) {
      sprintf(s, "Lagrange[%d]", d); goto found;
    }
    if (dof == getfem::normal_derivative_dof(d)) {
      sprintf(s, "D_n[%d]", d); goto found;
    }
    if (dof == getfem::global_dof(d)) {
      sprintf(s, "GlobalDof[%d]", d);
    }
    if (dof == getfem::mean_value_dof(d)) {
      sprintf(s, "MeanValue[%d]", d);
    }
    if (getfem::dof_xfem_index(dof) != 0) {
      sprintf(s, "Xfem[idx:%d]", int(dof_xfem_index(dof)));
    }
    
    for (dim_type r = 0; r < d; ++r) {
      if (dof == getfem::derivative_dof(d, r)) {
	sprintf(s, "D_%c[%d]", "xyzuvw"[r], d); goto found;
      }
      for (dim_type t = 0; t < d; ++t) {
	if (dof == getfem::second_derivative_dof(d, r, t)) {
	  sprintf(s, "D2%c%c[%d]", "xyzuvw"[r], "xyzuvw"[t], d); 
	  goto found;
	}
      }
    }
  }
 found:
  return s;
}


 /* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void crack_problem::init(void) {
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string FEM_TYPE_P  = PARAM.string_value("FEM_TYPE_P","FEM name p");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
					 "Name of simplex integration method");
  std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION");

  enrichment_option = enrichment_option_enum(PARAM.int_value("ENRICHMENT_OPTION",
							     "Enrichment option"));
 
 /* Affichage */

  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "FFEM_TYPE_P  =" << FEM_TYPE_P <<"\n";  
  cout << "INTEGRATION=" << INTEGRATION << "\n";
  pr1 = PARAM.real_value("P1", "First Elastic coefficient");
  pr2 = PARAM.real_value("P2", "Second Elastic coefficient");
  pr3 = PARAM.real_value("P3", "Third Elastic coefficient");
  AMP_LOAD_X = PARAM.real_value("AMP_LOAD_X", "Amp load x");
  AMP_LOAD_Y = PARAM.real_value("AMP_LOAD_Y", "Amp load y");
  nb_step  = PARAM.real_value("nb_step", "nb_step");
  reference_test = int(PARAM.int_value("REFERENCE_TEST", "Reference test"));
  cout << "AMP_LOAD_X = " << AMP_LOAD_X << "\n";
  cout << "AMP_LOAD_Y = " << AMP_LOAD_Y << "\n";
  
  /* First step : build the mesh */
 
  bgeot::pgeometric_trans pgt =  bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Nomber of space steps "));

  long int  NX_aff=PARAM.int_value("NX", "Nomber of space steps ");
  cout << "Number of space step NX= " << NX_aff << "\n";
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			    PARAM.int_value("MESH_NOISED") != 0);
  base_small_vector tt(N); tt[1] = -0.5;
  mesh.translation(tt); 
  
  cracktip.resize(2); // Coordonn�e du fond de fissure
  cracktip[0] = 0.9;
  cracktip[1] = 0.;

  scalar_type refinement_radius;
  refinement_radius = PARAM.real_value("REFINEMENT_RADIUS", "Refinement Radius");
  cout << "refinement_radius= " << refinement_radius << "\n";
  size_type refinement_process;
  refinement_process = PARAM.int_value("REFINEMENT_PROCESS", "Refinement process");
  cout << "refinement_process= " << refinement_process << "\n";
  
  if (refinement_radius > 0) {
    for (size_type ref = 0; ref < refinement_process; ++ref){
      dal::bit_vector conv_to_refine;
      for(dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i){
	for(size_type j=0; j < 3; ++j)
	  if(gmm::vect_dist2(mesh.points_of_convex(i)[j],cracktip)
	     < refinement_radius )
	    conv_to_refine.add(i);
      }
      mesh.Bank_refine(conv_to_refine);
      
      refinement_radius = refinement_radius/3.;
      cout <<"refining process step " << ref << " ... refining "
	   << conv_to_refine.size() <<" convexes..." << endl;
    }
    cout << "refinement process completed." << endl ;
  }

  mesh.write_to_file("toto.mesh");


  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL");
  if (residual == 0.) residual = 1e-10;
  enr_area_radius = PARAM.real_value("RADIUS_ENR_AREA",
				     "radius of the enrichment area");
  
  GLOBAL_FUNCTION_MF = PARAM.string_value("GLOBAL_FUNCTION_MF");
  GLOBAL_FUNCTION_U = PARAM.string_value("GLOBAL_FUNCTION_U");
  GLOBAL_FUNCTION_P = PARAM.string_value("GLOBAL_FUNCTION_P");

  


  cutoff_func = PARAM.int_value("CUTOFF_FUNC", "cutoff function");

  cutoff_radius = PARAM.real_value("CUTOFF", "Cutoff");
  cutoff_radius1 = PARAM.real_value("CUTOFF1", "Cutoff1");
  cutoff_radius0 = PARAM.real_value("CUTOFF0", "Cutoff0");
  mf_u().set_qdim(dim_type(N));
  

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);
  getfem::pintegration_method simp_ppi = 
    getfem::int_method_descriptor(SIMPLEX_INTEGRATION);
  getfem::pintegration_method sing_ppi = (SINGULAR_INTEGRATION.size() ?
		getfem::int_method_descriptor(SINGULAR_INTEGRATION) : 0);
  
  mim.set_integration_method(mesh.convex_index(), ppi);
  mls.add_level_set(ls);
  
  mim.set_simplex_im(simp_ppi, sing_ppi);
  mf_pre_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_pre_mortar.set_finite_element(mesh.convex_index(), 
				   getfem::fem_descriptor(PARAM.string_value("MORTAR_FEM_TYPE")));
  mf_mult.set_finite_element(mesh.convex_index(), pf_u);
  mf_mult.set_qdim(dim_type(N));
  mf_partition_of_unity.set_classical_finite_element(1);

  mixed_pressure =
    (PARAM.int_value("MIXED_PRESSURE","Mixed version or not.") != 0);
  dir_with_mult = unsigned(PARAM.int_value("DIRICHLET_VERSION",
					   "Version of Dirichlet"));
  if (mixed_pressure) {
        getfem::pfem pf_p = 
        getfem::fem_descriptor(FEM_TYPE_P);

  mf_pre_p.set_finite_element(mesh.convex_index(), pf_p);
  mf_mult_p.set_finite_element(mesh.convex_index(), pf_p);

  // mf_p.set_finite_element(mesh.convex_index(), pf_p);
  
    
  }

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
    
    if (un[0]  > 0.5) mesh.region(NEUMANN1_BOUNDARY_NUM).add(i.cv(), i.f());
    if (un[1]  > 0.5) mesh.region(NEUMANN2_BOUNDARY_NUM).add(i.cv(), i.f());
    if (un[0]  < -0.5) mesh.region(NEUMANN3_BOUNDARY_NUM).add(i.cv(), i.f());
    if (un[1]  < -0.5) mesh.region(NEUMANN4_BOUNDARY_NUM).add(i.cv(), i.f());
  }
  
  
}


base_small_vector ls_function(const base_node P, int num = 0) {
  scalar_type x = P[0], y = P[1];
  base_small_vector res(2);
  switch (num) {
    case 0: {
      res[0] = y;
      res[1] = -.5 + x;
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



struct matrix_G {

  const sparse_matrix &B;
  const sparse_matrix &S;
  mutable plain_vector W1, W2;

  gmm::SuperLU_factor<scalar_type> SLUF;

  matrix_G(const sparse_matrix &BB, const sparse_matrix &SS)
    : B(BB), S(SS), W1(gmm::mat_nrows(SS)), W2(gmm::mat_nrows(SS)) {
    SLUF.build_with(SS);
  }
  
};

template <typename vector1, typename vector2>
void mult(const matrix_G &G, const vector1 &X, vector2 &Y) {
  gmm::mult(gmm::transposed(G.B), X, G.W1);
  // gmm::iteration it(1E-6, 0);
  // gmm::cg(G.S, G.W2, G.W1,  gmm::identity_matrix(), it);
  G.SLUF.solve(G.W2, G.W1);
  gmm::mult(G.B, G.W2, Y);
}

template <typename vector1, typename vector2>
void mult(const matrix_G &G, const vector1 &X, const vector2 &b, vector2 &Y)
{ mult(G, X, Y); gmm::add(b, Y); }


scalar_type smallest_eigen_value(const sparse_matrix &B,
				 const sparse_matrix &M,
				 const sparse_matrix &S) {
  cout << "matrice B = "<< B;
  size_type n = gmm::mat_nrows(M);
  scalar_type lambda;
  plain_vector V(n), W(n), V2(n);
  gmm::fill_random(V2);
  matrix_G G(B, S);
  
  do {
    gmm::copy(V2, V);
    gmm::scale(V, 1./gmm::vect_norm2(V));
    gmm::mult(M, V, W);
    
    gmm::iteration it(1E-3, 0);
    gmm::cg(G, V2, W,  gmm::identity_matrix(), it);    
    lambda = gmm::vect_norm2(V2);

//  compute the Rayleigh quotient
//     mult(G, V2, W);
//     scalar_type lambda2 = gmm::vect_sp(V2, W);
//     gmm::mult(M, V2, W);
//     lambda2 /= gmm::vect_sp(V2, W);
//     cout << "lambda2 = " << sqrt(lambda2) << endl;

    cout << "lambda = " << sqrt(1./lambda) << endl;
    cout << "residu = " << gmm::vect_dist2(V2, gmm::scaled(V, lambda)) << endl;
    
  } while (gmm::vect_dist2(V2, gmm::scaled(V, lambda)) > 1E-3);
  
  return sqrt(1./lambda);
}
  


/*****************************************************************/
/*   Model.                                                      */
/*****************************************************************/

bool crack_problem::solve(plain_vector &U, plain_vector &P) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();
  ls.reinit();
  size_type law_num = PARAM.int_value("LAW");
  //size_type newton_version = PARAM.int_value("newton_version");
  base_vector pr(3); pr[0] = pr1; pr[1] = pr2; pr[2] = pr3;

  for (size_type d = 0; d < ls.get_mesh_fem().nb_basic_dof(); ++d) {
    ls.values(0)[d] = ls_function(ls.get_mesh_fem().point_of_basic_dof(d), 0)[0];
    ls.values(1)[d] = ls_function(ls.get_mesh_fem().point_of_basic_dof(d), 0)[1];
  }
  ls.touch();

  mls.adapt();
  mim.adapt();
  mfls_u.adapt();
  mfls_p.adapt();
  mfls_mortar.adapt(); mfls_mortar.set_qdim(2);

  bool load_global_fun = GLOBAL_FUNCTION_MF.size() != 0;


  cout << "Setting up the singular functions for the enrichment\n";

  std::vector<getfem::pglobal_function> vfunc(2);
  std::vector<getfem::pglobal_function> vfunc_p(2);
  if (!load_global_fun) {
  std::cout << "Using default singular functions\n";
    for (unsigned i = 0; i < vfunc.size(); ++i){
      /* use the singularity */
     
      getfem::abstract_xy_function *s = 
	new getfem::crack_singular_xy_function(i+12);
      
      if (enrichment_option != FIXED_ZONE && 
	  enrichment_option != GLOBAL_WITH_MORTAR) {
	/* use the product of the singularity function
	   with a cutoff */
	getfem::abstract_xy_function *c = 
	  new getfem::cutoff_xy_function(int(cutoff_func),
					 cutoff_radius, 
					 cutoff_radius1,cutoff_radius0);
	s  = new getfem::product_of_xy_functions(*s, *c);
        
      }
      vfunc[i]=getfem::global_function_on_level_set(ls, *s);           
    }

     for (unsigned i = 0; i < vfunc_p.size() ;++i){
       /* use the singularity */
      
       getfem::abstract_xy_function *sp = 
 	new getfem::crack_singular_xy_function(i+14);      
       if (enrichment_option != FIXED_ZONE && 
 	  enrichment_option != GLOBAL_WITH_MORTAR) {
 	/* use the product of the singularity function
 	   with a cutoff */
 	getfem::abstract_xy_function *cp = 
	  new getfem::cutoff_xy_function(int(cutoff_func),
 					 cutoff_radius, 
 					 cutoff_radius1,cutoff_radius0);
 	sp  = new getfem::product_of_xy_functions(*sp, *cp);
        
       }
       vfunc_p[i]=getfem::global_function_on_level_set(ls, *sp);           
    }
  } else {
    cout << "Load singular functions from " << GLOBAL_FUNCTION_MF << " and " << GLOBAL_FUNCTION_U <<" and " << GLOBAL_FUNCTION_P << "\n";
    getfem::mesh *m = new getfem::mesh(); 
    m->read_from_file(GLOBAL_FUNCTION_MF);
    getfem::mesh_fem *mf_c = new getfem::mesh_fem(*m); 
    mf_c->read_from_file(GLOBAL_FUNCTION_MF);
    std::fstream f(GLOBAL_FUNCTION_U.c_str(), std::ios::in);
    std::fstream fp(GLOBAL_FUNCTION_P.c_str(), std::ios::in);
    plain_vector W(mf_c->nb_dof());
    plain_vector WP(mf_c->nb_dof());


  
    for (unsigned i=0; i < mf_c->nb_dof(); ++i) {
      f >> W[i]; GMM_ASSERT1(f.good(), "problem while reading " << GLOBAL_FUNCTION_U);
      
      //cout << "The precalculated dof " << i << " of coordinates " << mf_c->point_of_dof(i) << " is "<< W[i] <<endl; 
      /*scalar_type x = pow(mf_c->point_of_dof(i)[0],2); scalar_type y = pow(mf_c->point_of_dof(i)[1],2);
	scalar_type r = std::sqrt(pow(x,2) + pow(y,2));
	scalar_type sgny = (y < 0 ? -1.0 : 1.0);
	scalar_type sin2 = sqrt(gmm::abs(.5-x/(2*r))) * sgny;
	scalar_type cos2 = sqrt(gmm::abs(.5+x/(2*r)));
	W[i] = std::sqrt(r) * sin2;
      */
    }
    unsigned nb_func = mf_c->get_qdim();
    cout << "read " << nb_func << " global functions OK.\n";
    vfunc.resize(nb_func);
    getfem::interpolator_on_mesh_fem *global_interp = 
      new getfem::interpolator_on_mesh_fem(*mf_c, W);
    for (size_type i=0; i < nb_func; ++i) {
      /* use the precalculated function for the enrichment*/
      //getfem::abstract_xy_function *s = new getfem::crack_singular_xy_function(i);
      getfem::abstract_xy_function *s = new getfem::interpolated_xy_function(*global_interp,i);

      if (enrichment_option != FIXED_ZONE && 
	  enrichment_option != GLOBAL_WITH_MORTAR) {

	/* use the product of the enrichment function
	   with a cutoff */
	getfem::abstract_xy_function *c = 
	  new getfem::cutoff_xy_function(int(cutoff_func),
					 cutoff_radius, 
					 cutoff_radius1,cutoff_radius0);
	s = new getfem::product_of_xy_functions(*s, *c);
      }    
      vfunc[i] = getfem::global_function_on_level_set(ls, *s);
    }    
  }
  
  mf_sing_u.set_functions(vfunc);
  mf_sing_p.set_functions(vfunc_p);  

  if (enrichment_option == SPIDER_FEM_ALONE || 
      enrichment_option == SPIDER_FEM_ENRICHMENT) {
    spider.fem = new getfem::spider_fem(spider.radius, mim, spider.Nr,
					spider.Ntheta, spider.K, cracktip,
					spider.theta0, spider.bimat_enrichment,
					spider.epsilon);
    mf_us.set_finite_element(mesh.convex_index(),spider.fem->get_pfem());
    for (dal::bv_visitor_c i(mf_us.convex_index()); !i.finished(); ++i) {
      if (mf_us.fem_of_element(i)->nb_dof(i) == 0) {
	mf_us.set_finite_element(i,0);
      }
    }
    spider.fem->check();
  }

  switch (enrichment_option) {

  case FIXED_ZONE :
    {
      dal::bit_vector enriched_dofs;
      plain_vector X(mf_partition_of_unity.nb_dof());
      plain_vector Y(mf_partition_of_unity.nb_dof());
      getfem::interpolation(ls.get_mesh_fem(), mf_partition_of_unity,
			    ls.values(1), X);
      getfem::interpolation(ls.get_mesh_fem(), mf_partition_of_unity,
			    ls.values(0), Y);
      for (size_type j = 0; j < mf_partition_of_unity.nb_dof(); ++j) {
	if (gmm::sqr(X[j]) + gmm::sqr(Y[j]) <= gmm::sqr(enr_area_radius))
	  enriched_dofs.add(j);
      }
      if (enriched_dofs.card() < 3)
	GMM_WARNING0("There is " << enriched_dofs.card() <<
		     " enriched dofs for the crack tip");
      mf_product.set_enrichment(enriched_dofs);
      mf_u_sum.set_mesh_fems(mf_product, mfls_u);
      mf_product_p.set_enrichment(enriched_dofs);
      mf_p_sum.set_mesh_fems(mf_product_p, mfls_p);
    }
    break;
  

    case GLOBAL_WITH_MORTAR: {
      // Selecting the element in the enriched domain

      dal::bit_vector cvlist_in_area;
      dal::bit_vector cvlist_out_area;
      for (dal::bv_visitor cv(mesh.convex_index()); 
	   !cv.finished(); ++cv) {
	bool in_area = true;
	/* For each element, we test all of its nodes. 
	   If all the nodes are inside the enrichment area,
	   then the element is completly inside the area too */ 
	for (unsigned j=0; j < mesh.nb_points_of_convex(cv); ++j) {
	  if (gmm::sqr(mesh.points_of_convex(cv)[j][0] - cracktip[0]) + 
	      gmm::sqr(mesh.points_of_convex(cv)[j][1] - cracktip[1]) > 
	      gmm::sqr(enr_area_radius)) {
	    in_area = false; break;
	  }
	}

	/* "remove" the global function on convexes outside the enrichment
	   area */
	if (!in_area) {
	  cvlist_out_area.add(cv);
	  mf_sing_u.set_finite_element(cv, 0);
	  mf_u().set_dof_partition(cv, 1);
	  //mf_sing_p.set_finite_element(cv, 0);
	  //mf_p().set_dof_partition(cv, 1);
	  
	} else cvlist_in_area.add(cv);
      }

      /* extract the boundary of the enrichment area, from the
	 "inside" point-of-view, and from the "outside"
	 point-of-view */
      getfem::mesh_region r_border, r_enr_out;
      getfem::outer_faces_of_mesh(mesh, r_border);

      getfem::outer_faces_of_mesh(mesh, cvlist_in_area, 
				  mesh.region(MORTAR_BOUNDARY_IN));
      getfem::outer_faces_of_mesh(mesh, cvlist_out_area, 
				  mesh.region(MORTAR_BOUNDARY_OUT));
      for (getfem::mr_visitor v(r_border); !v.finished(); ++v) {
	mesh.region(MORTAR_BOUNDARY_OUT).sup(v.cv(), v.f());
      }
      mf_u_sum.set_mesh_fems(mf_sing_u, mfls_u);
      mf_p_sum.set_mesh_fems(mf_sing_p, mfls_p);
    } break;

  case GLOBAL_WITH_CUTOFF :{
    if(cutoff_func == 0)
      cout<<"Using exponential Cutoff..."<<endl;
    else
      cout<<"Using Polynomial Cutoff..."<<endl;
      mf_u_sum.set_mesh_fems(mf_sing_u, mfls_u);
      mf_p_sum.set_mesh_fems(mf_sing_p, mfls_p);
  } break;

  case SPIDER_FEM_ALONE : { 
    mf_u_sum.set_mesh_fems(mf_us); 
  } break;
    
  case SPIDER_FEM_ENRICHMENT : {
    mf_u_sum.set_mesh_fems(mf_us, mfls_u); 
  } break;
    
  case NO_ENRICHMENT: {
    cout<<"No enrichment..."<<endl;
    mf_u_sum.set_mesh_fems(mfls_u);
    mf_p_sum.set_mesh_fems(mfls_p);
  } break;
  
  }
  

  U.resize(mf_u().nb_basic_dof());
  P.resize(mf_pe().nb_basic_dof());
  


  if (mixed_pressure)
    //cout << "Number of dof for P mfls: " << mfls_p().nb_basic_dof() << endl;
    cout << "Number of dof for P mf_pe: " << mf_pe().nb_basic_dof() << endl;
    cout << "Number of dof for u: " << mf_u().nb_basic_dof() << endl;

  unsigned Q = mf_u().get_qdim();
  if (0) {
    for (unsigned d=0; d < mf_u().nb_dof(); d += Q) {
      printf("dof %4d @ %+6.2f:%+6.2f: ", d, 
	     mf_u().point_of_basic_dof(d)[0], mf_u().point_of_basic_dof(d)[1]);

      const getfem::mesh::ind_cv_ct cvs = mf_u().convex_to_basic_dof(d);
      for (unsigned i=0; i < cvs.size(); ++i) {
	size_type cv = cvs[i];
	//if (pm_cvlist.is_in(cv)) flag1 = true; else flag2 = true;

	getfem::pfem pf = mf_u().fem_of_element(cv);
	unsigned ld = unsigned(-1);
	for (unsigned dd = 0; dd < mf_u().nb_basic_dof_of_element(cv); dd += Q) {
	  if (mf_u().ind_basic_dof_of_element(cv)[dd] == d) {
	    ld = dd/Q; break;
	  }
	}
	if (ld == unsigned(-1)) {
	  cout << "DOF " << d << "NOT FOUND in " << cv << " BUG BUG\n";
	} else {
	  printf(" %3d:%.16s", int(cv), name_of_dof(pf->dof_types().at(ld)).c_str());
	}
      }
      printf("\n");
    }
  }
  
  ///////////////////////////////////////////////////////////////
  // find the dofs on the upper right and lower right corners  //
  ///////////////////////////////////////////////////////////////
  cout << "Find the dofs on the upper right and lower right corners" << endl;
  scalar_type d1 = 1.0, d2 = 1.0;
  size_type icorner1 = size_type(-1), icorner2 = size_type(-1);
  base_node corner1 = base_node(1.0, -0.5);
  base_node corner2 = base_node(1.0, 0.5);
  GMM_ASSERT1(!(mf_u().is_reduced()), "To be adapted for reduced fems");
  for (size_type i = 0; i < mf_u().nb_basic_dof(); i+=N) {
    scalar_type dd1 = gmm::vect_dist2(mf_u().point_of_basic_dof(i), corner1);
    if (dd1 < d1) { icorner1 = i; d1 = dd1; }
    scalar_type dd2 = gmm::vect_dist2(mf_u().point_of_basic_dof(i), corner2);
    if (dd2 < d2) { icorner2 = i; d2 = dd2; }
  }

  GMM_ASSERT1(((d1 < 1E-8) && (d2 < 1E-8)),
	      "Upper right or lower right corners not found d1 = "
	      << d1 << " d2 = " << d2);

  /***********************************************/
  /*                                             */   
  /* Add Nonlinear-elasticity brick.             */
  /*                                             */
  /***********************************************/
  // Choose the material law.
 
 //  getfem::abstract_hyperelastic_law *pl = 0;
 //     switch (law_num) {
 //       case 0: pl = new getfem::Mooney_Rivlin_hyperelastic_law(); break;
 //       case 1: pl = new getfem::SaintVenant_Kirchhoff_hyperelastic_law(); break;
 //       case 3: pl = new getfem::Ciarlet_Geymonat_hyperelastic_law(); break;
 //       default: GMM_ASSERT1(false, "no such law");
 //     }

  // if (mixed_pressure) 
  //   cout << "Number of dof for P: " << mf_pe.nb_dof() << endl;
  //   cout << "Number of dof for u: " << mf_u.nb_dof() << endl;
  //pr.resize(pl->nb_params());
  
  getfem::model model;
  size_type alpha,beta;
  model.add_fem_variable("u", mf_u());
  //model.add_fixed_size_variable("alpha", alpha);
  //model.add_fixed_size_variable("beta", beta);
  
  model.add_initialized_fixed_size_data("Mooney_Riv_coefficient", pr);
  
  getfem::add_nonlinear_elasticity_brick(model, mim, "u",
					 getfem::Mooney_Rivlin_hyperelastic_law(), "Mooney_Riv_coefficient");
  // Incompressibility
  if (mixed_pressure) {
    model.add_fem_variable("p", mf_pe());
    getfem::add_nonlinear_incompressibility(model, mim, "u" , "p", size_type(-1));
    }
  // Add nonlinear elasticity optimazition brick

  add_nonlinear_elasticity_optim_brick(model, mim, "u" , "p", "alpha", "beta",
	   getfem::Mooney_Rivlin_hyperelastic_law(), "Mooney_Riv_coefficient", ls);

 //  getfem::mdbrick_nonlinear_elasticity<>  ELAS(*pl, mim, mf_u(), pr);
 //   getfem::mdbrick_nonlinear_incomp<> INCOMP(ELAS, mf_pe());







  // Defining the Neumann condition right hand side.
  plain_vector F_Neumann1(nb_dof_rhs * N);
  plain_vector F_Neumann2(nb_dof_rhs * N);
  plain_vector F_Neumann3(nb_dof_rhs * N);
  plain_vector F_Neumann4(nb_dof_rhs * N);
  // Neumann condition brick.
   
  for(size_type i = 0; i < F_Neumann1.size(); i=i+N) F_Neumann1[i] = AMP_LOAD_X;
  for(size_type i = 1; i < F_Neumann1.size(); i=i+N) F_Neumann1[i] = AMP_LOAD_Y;
  for(size_type i = 0; i < F_Neumann2.size(); i=i+N) F_Neumann2[i] = AMP_LOAD_X;
  for(size_type i = 1; i < F_Neumann2.size(); i=i+N) F_Neumann2[i] = AMP_LOAD_Y;
  for(size_type i = 0; i < F_Neumann3.size(); i=i+N) F_Neumann3[i] = -AMP_LOAD_X;
  for(size_type i = 1; i < F_Neumann3.size(); i=i+N) F_Neumann3[i] = -AMP_LOAD_Y;
  for(size_type i = 0; i < F_Neumann4.size(); i=i+N) F_Neumann4[i] = -AMP_LOAD_X;
  for(size_type i = 1; i < F_Neumann4.size(); i=i+N) F_Neumann4[i] = -AMP_LOAD_Y;
   
   model.add_initialized_fem_data("NeumannData1", mf_rhs,F_Neumann1 );
  getfem::add_normal_source_term_brick
    (model, mim, "u", "NeumannData1", NEUMANN1_BOUNDARY_NUM);

   model.add_initialized_fem_data("NeumannData2", mf_rhs,F_Neumann2 );
  getfem::add_normal_source_term_brick
    (model, mim, "u", "NeumannData2", NEUMANN2_BOUNDARY_NUM);

   model.add_initialized_fem_data("NeumannData3", mf_rhs,F_Neumann3 );
   getfem::add_normal_source_term_brick
    (model, mim, "u", "NeumannData3", NEUMANN3_BOUNDARY_NUM);

   model.add_initialized_fem_data("NeumannData4", mf_rhs,F_Neumann3 );
   getfem::add_normal_source_term_brick
    (model, mim, "u", "NeumannData4", NEUMANN4_BOUNDARY_NUM);
   
  
  //  getfem::mdbrick_source_term<> NEUMANN1(INCOMP, mf_rhs, F_Neumann1, NEUMANN1_BOUNDARY_NUM);
  //  getfem::mdbrick_source_term<> NEUMANN2(NEUMANN1, mf_rhs, F_Neumann2, NEUMANN2_BOUNDARY_NUM);
  //  gmm::scale(F, -1.0);
  //  getfem::mdbrick_source_term<> NEUMANN3(NEUMANN2, mf_rhs, F_Neumann3, NEUMANN3_BOUNDARY_NUM);
  //  getfem::mdbrick_source_term<> NEUMANN4(NEUMANN3, mf_rhs, F_Neumann4, NEUMANN4_BOUNDARY_NUM);
  //  getfem::mdbrick_constraint<> KILL_RIGID_MOTIONS(NEUMANN4);

  GMM_ASSERT1(N==2, "To be corrected for 3D computation");
  sparse_matrix BB(4, mf_u().nb_dof());
  BB(0, icorner1) = 1.0;
  BB(1, icorner1+1) = 1.0;
  BB(2, icorner2) = 1.0;
  BB(3, icorner2+1) = 1.0;
   size_type size(3);
   std::vector<scalar_type> LRH(size);
   model.add_fixed_size_variable("dir", size);
   getfem::add_constraint_with_multipliers(model, "u", "dir", BB, LRH);
   
   // KILL_RIGID_MOTIONS.set_constraints(BB, plain_vector(4));
   // KILL_RIGID_MOTIONS.set_constraints_type(getfem::constraints_type(dir_with_mult));

   // Dirichlet condition brick.
   // getfem::mdbrick_Dirichlet<> DIRICHLET (KILL_RIGID_MOTIONS, DIRICHLET_BOUNDARY_NUM, mf_mult);
    
   // DIRICHLET.set_constraints_type(getfem::constraints_type(PARAM.int_value("DIRICHLET_VERSION")));

   // getfem::mdbrick_abstract<> *final_model = &DIRICHLET;

/*************************************/
/*       Generic solve.              */
/*************************************/

   // getfem::standard_model_state MS(*final_model);
   // size_type maxit = PARAM.int_value("MAXITER"); 
   // gmm::iteration iter;
  
   //size_type stnst = PARAM.int_value("stnst");
  
  //if (stnst==1) {
 
 /**********************/
 /*   Step loading     */
 /**********************/
  // cout << "By step amigoo###################################<<<<<>>>>>>######" << endl;
  // cout << "Nb de Step " << nb_step <<endl;     
  // for (int step = 0; step < nb_step; ++step) {
  // plain_vector DF(F);

  //gmm::copy(gmm::scaled(F, (step+1.)/(scalar_type)nb_step), DF);
  //NEUMANN4.source_term().set(DF);

       /************************************/
       /* increment  imposed displacement  */
       /************************************/ 
     
  //KILL_RIGID_MOTIONS.set_constraints(BB, plain_vector(3));

    /********************************************************/
    /* let the default non-linear solve (Newton) do its job */
    /********************************************************/ 

    // cout << "step " << step << ", number of variables : " << final_model->nb_dof() << endl;
    // cout << "DF " << DF <<  endl;

    // iter = gmm::iteration(residual, int(PARAM.int_value("NOISY", "Noisy = ")), maxit ? maxit : 40000);
    
    // gmm::abstract_newton_line_search alnrs;
    // gmm::simplest_newton_line_search silnrs;
    // gmm::default_newton_line_search dlnrs;
    // gmm::systematic_newton_line_search sylnrs;
    
    // getfem::standard_solve(MS,*final_model, iter, getfem::default_linear_solver(*final_model), dlnrs);

    // pl->reset_unvalid_flag();
    //   final_model->compute_residual(MS);
    // if (pl->get_unvalid_flag())
    //  GMM_WARNING1("The solution is not completely valid, the determinant "
    //               "of the transformation is negative on "
    //               << pl->get_unvalid_flag() << " gauss points");

  //}   

  //}
  //else{
  //cout << "Resolution direct =============D=I=R=E=C=T====================> " << endl;
 /****************/
 /*Sans iteration*/
 /****************/
 // iter = gmm::iteration(residual, int(PARAM.int_value("NOISY", "Noisy = ")),
// 			  maxit ? maxit : 40000);
    
    
//  // gmm::abstract_newton_line_search alnrs;
//     gmm::simplest_newton_line_search silnrs;
//     gmm::default_newton_line_search dlnrs;
//     gmm::systematic_newton_line_search sylnrs;
    
//     switch (newton_version){
    
//  // case 0: getfem::standard_solve(MS,*final_model, iter, getfem::default_linear_solver(*final_model), alnrs);
//     case 1:{
//       getfem::standard_solve(MS,*final_model, iter, getfem::default_linear_solver(*final_model), silnrs);
//       cout << "============================ " << endl;
//       cout << "=:simplest_newton_line_search= " << endl;
//       cout << "============================ " << endl;
//     }break;
//     case 2:{
//       getfem::standard_solve(MS,*final_model, iter, getfem::default_linear_solver(*final_model), dlnrs);
//       cout << "============================ " << endl;
//       cout << "=default_newton_line_search= " << endl;
//       cout << "============================ " << endl;
//     }break;
    

//     case 3: {
//       getfem::standard_solve(MS,*final_model, iter, getfem::default_linear_solver(*final_model), sylnrs);
//       cout << "=============================== " << endl;
//       cout << "=systematic_newton_line_search= " << endl;
//       cout << "=============================== " << endl;
//     }break;
//     default: GMM_ASSERT1(false, "No such newton");
//     }
//     pl->reset_unvalid_flag();
//     final_model->compute_residual(MS);
//     if (pl->get_unvalid_flag()) 
//       GMM_WARNING1("The solution is not completely valid, the determinant "
// 		   "of the transformation is negative on "
// 		   << pl->get_unvalid_flag() << " gauss points");
//   }
  gmm::iteration iter(residual, 1, 40000);
  cout << "Solving..." << endl;
  iter.init();
  getfem::standard_solve(model, iter);
  gmm::resize(U, mf_u().nb_dof());
  gmm::copy(model.real_variable("u"), U);



//  /*************************************/
//  /*  Computation of the inf-sup bound */ 
//  /*************************************/

   // if (PARAM.int_value("INF_SUP_COMP") && mixed_pressure) {
    
//      cout << "Sparse matrices computation for the test of inf-sup condition"
//  	 << endl;

//     sparse_matrix Mis(mf_pe().nb_dof(), mf_pe().nb_dof());
//     sparse_matrix Sis(mf_u().nb_dof(), mf_u().nb_dof());
//     sparse_matrix Bis(mf_pe().nb_dof(),mf_u().nb_dof());
//     getfem::asm_mass_matrix(Mis, mim, mf_pe());
//     getfem::asm_stiffness_matrix_for_homogeneous_laplacian_componentwise(Sis, mim, mf_u());
//     getfem::asm_mass_matrix(Sis, mim, mf_u());    
//     cout << "Inf-sup condition test" << endl;
//     INCOMP.get_B(MS,Bis);
//     scalar_type lambda = smallest_eigen_value(Bis, Mis, Sis);
//     cout << "The inf-sup test gives " << lambda << endl;
//   }

   
//   gmm::iteration iter(residual, 1, 40000);
//   cout << "Solving..." << endl;
//   iter.init();
//   getfem::standard_solve(model, iter);
//   gmm::resize(U, mf_u().nb_dof());
//   gmm::copy(model.real_variable("u"), U);
 
 
  

  return (iter.converged());




  // Solution extraction
    
  //   gmm::copy(ELAS.get_solution(MS), U);
//     gmm::copy(INCOMP.get_pressure(MS), P);
//     return (iter.converged());

      
      
      
 if (reference_test == 1) {
      
      cout << "Exporting reference solution...";
      dal::bit_vector blocked_dof = mf_u().basic_dof_on_region(5);
      getfem::mesh_fem mf_refined(mesh, dim_type(N));
      std::string FEM_DISC = PARAM.string_value("FEM_DISC","fem disc ");
      mf_refined.set_finite_element(mesh.convex_index(),
				    getfem::fem_descriptor(FEM_DISC));
      
      plain_vector W(mf_refined.nb_dof());
      getfem::interpolation(mf_u(), mf_refined, U, W);
      
      
      mf_refined.write_to_file(datafilename + "_refined_test.meshfem_refined", true);
      gmm::vecsave(datafilename + "_refined_test.U_refined", W);

      
      mf_refined.set_qdim(1);
      plain_vector PP(mf_refined.nb_dof());
      getfem::interpolation(mf_pe(), mf_refined, P, PP);
      mf_refined.write_to_file(datafilename + "_refined_test.p_meshfem_refined", true);
      gmm::vecsave(datafilename + "_refined_test.P_refined", PP);
      cout << "done" << endl;
  }
  }

void export_interpolated_on_line(const getfem::mesh_fem &mf,
				 const getfem::base_vector &U,
				 const base_node &x0,
				 const base_small_vector &dir,
				 const int nb_points, 
				 const std::string &filename) {
  getfem::mesh_trans_inv mti(mf.linked_mesh());
  scalar_type h = 1.0/(2*nb_points);
  for (int i=-nb_points; i <= nb_points; ++i) {
    mti.add_point(x0 + 2*(i*h)*dir);
  }

  getfem::base_vector V(mti.nb_points() * mf.get_qdim());
  getfem::base_matrix M;
  getfem::interpolation(mf, mti, U, V, M, 0, false);
  
  std::ofstream f(filename.c_str()); f.precision(16);
  
  for (size_type i=0; i < mti.nb_points(); ++i) {
    for (unsigned q=0; q < mf.get_qdim(); ++q) {
      f << V[i*mf.get_qdim()+q] << " ";
    }
    f << "\n";
  }
}
				 

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  //getfem::getfem_mesh_level_set_noisy();


  crack_problem p;
  p.PARAM.read_command_line(argc, argv);
  p.init();
  p.mesh.write_to_file(p.datafilename + ".mesh");
  
  plain_vector U, P;
  if (!p.solve(U, P)) GMM_ASSERT1(false,"Solve has failed");
  
  cout << "Saving the solution" << endl;
  getfem::mesh mcut;
  p.mls.global_cut_mesh(mcut);
  unsigned Q = p.mf_u().get_qdim();
  
  getfem::mesh_fem mf(mcut, dim_type(Q));
  mf.set_classical_discontinuous_finite_element(2, 1E-7);
  
  plain_vector V(mf.nb_dof());
  getfem::interpolation(p.mf_u(), mf, U, V);
  mf.write_to_file(p.datafilename + ".meshfem", true);
  gmm::vecsave(p.datafilename + ".U", V);
  
  getfem::mesh_fem mf_p(mcut);
  mf_p.set_classical_discontinuous_finite_element(2, 1E-7);
  
  plain_vector PP(mf_p.nb_dof());
  getfem::interpolation(p.mf_pe(), mf_p, P, PP);
  mf_p.write_to_file(p.datafilename + ".p_meshfem", true);
  gmm::vecsave(p.datafilename + ".P", PP);
      
  cout << "Interpolating solution for the drawing" << endl;
  getfem::stored_mesh_slice sl;
  getfem::mesh mcut_refined;
  //getfem::mesh mcut_refined_p;
  
  unsigned NX = unsigned(p.PARAM.int_value("NX")), nn;
  if (NX < 6) nn = 24;
  else if (NX < 12) nn = 6;
  else if (NX < 30) nn = 3;
  else nn = 3;
  
  /* choose an adequate slice refinement based on the distance
     to the crack tip */
  std::vector<bgeot::short_type> nrefine(mcut.convex_index().last_true()+1);
  for (dal::bv_visitor cv(mcut.convex_index()); !cv.finished(); ++cv) {
    scalar_type dmin=0, d;
    base_node Pmin,Pp;
    for (unsigned i=0; i < mcut.nb_points_of_convex(cv); ++i) {
      Pp = mcut.points_of_convex(cv)[i];
      d = gmm::vect_norm2(ls_function(Pp));
      if (d < dmin || i == 0) { dmin = d; Pmin = Pp; }
    }
    
    if (dmin < 1e-5)
      nrefine[cv] = short_type(nn*8);
    else if (dmin < .1) 
      nrefine[cv] = short_type(nn*2);
    else nrefine[cv] = short_type(nn);
    // if (dmin < .01)
    //  cout << "cv: "<< cv << ", dmin = " << dmin << "Pmin=" << Pmin
    //       << " " << nrefine[cv] << "\n";
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
  
  // getfem::mesh_im mim_refined_p(mcut_refined_p); 
  // mim_refined.set_integration_method(getfem::int_method_descriptor
  //				     ("IM_TRIANGLE(6)"));  
  mcut.write_to_file(p.datafilename + ".meshvm");
  
  getfem::mesh_fem mf_refined(mcut_refined, dim_type(Q));
  mf_refined.set_classical_discontinuous_finite_element(2,0.001);
  plain_vector W(mf_refined.nb_dof());
  
  //getfem::mesh_fem mf_refined_p(mcut_refined_p, dim_type(Q));
  //mf_refined_p.set_classical_discontinuous_finite_element(2,0.001);
  //plain_vector PPW(mf_refined.nb_dof());
  //getfem::interpolation(p.mf_pe(), mf_refined, P, PPW);
  
  getfem::interpolation(p.mf_u(), mf_refined, U, W);
  mf_refined.write_to_file(p.datafilename + ".meshfemuvm", true);
  gmm::vecsave(p.datafilename + ".Uvm", W);
  

  
  getfem::mesh_fem mf_vm(mcut_refined,  1);
  mf_vm.set_classical_discontinuous_finite_element(2, 0.001);
  plain_vector Vm(mf_vm.nb_dof());
 
  cout << "compute Von_mises" << endl;
  getfem::interpolation_von_mises(mf_refined, mf_vm, W, Vm);
  
  gmm::vecsave(p.datafilename + ".VM",Vm);
  mf_vm.write_to_file(p.datafilename + ".meshfemvm", true); 
  
  if (p.PARAM.int_value("VTK_EXPORT")) {
    
    getfem::mesh_fem mf_refined_vm(mcut_refined, 1);
    getfem::mesh_fem mf_refined_p(mcut_refined, 1);

    mf_refined_vm.set_classical_discontinuous_finite_element(2, 0.001);
    mf_refined_p.set_classical_discontinuous_finite_element(2, 0.001);

    cerr << "mf_refined_vm.nb_dof=" << mf_refined_vm.nb_dof() << "\n";
    
    plain_vector VM(mf_refined_vm.nb_dof());
    plain_vector PPW(mf_refined_p.nb_dof());
   
    cout << "computing von mises\n";
    
    getfem::interpolation_von_mises(mf_refined, mf_refined_vm, W, VM);
    getfem::interpolation(p.mf_pe(), mf_refined_p, P, PPW);
    plain_vector D(mf_refined_vm.nb_dof() * Q), 
      DN(mf_refined_vm.nb_dof());
    
    cout << "export to " << p.datafilename + ".vtk" << "..\n";
    getfem::vtk_export exp(p.datafilename + ".vtk",
			   p.PARAM.int_value("VTK_EXPORT")==1);
    
    exp.exporting(mf_refined); 
    
    //exp.write_point_data(mf_refined_vm, DN, "error");
    getfem::interpolation(p.mf_pe(), mf_refined_p, P, PPW);
    exp.write_point_data(mf_refined_vm, VM, "von_mises_stress");
    exp.write_point_data(mf_refined, W, "elastostatic_displacement");
    exp.write_point_data(mf_refined_vm, PPW, "Pressure");

    base_node line_x0(0.70001,0);
    base_small_vector line_dir(0, 0.5001);
    unsigned line_nb_points = 1000;
    export_interpolated_on_line(mf_refined_vm, VM, 
				line_x0, line_dir, line_nb_points,
				"von_mises_on_line.data");
    
    
    cout << "export done, you can view the data file with (for example)\n"
      "mayavi -d " << p.datafilename << ".vtk -f  "
      "WarpVector -m BandedSurfaceMap -m Outline\n";
  }
    
  if(p.PARAM.int_value("ERROR_TO_REF_SOL") == 1){
    cout << "Computing error with respect to a reference solution..." << endl;
    
    std::string REFERENCE_MF = "Solution_160_ref_raffine_P1BP1.meshfem";
    std::string REFERENCE_U = "Solution_160_ref_raffine_P1BP1.U";
    std::string REFERENCE_MFP = "Solution_160_ref_raffine_P1BP1.p_meshfem";
    std::string REFERENCE_P = "Solution_160_ref_raffine_P1BP1.P";
    
    cout << "Load reference displacement from "
	 << REFERENCE_MF << " and " << REFERENCE_U << "\n";
    getfem::mesh ref_m; 
    ref_m.read_from_file(REFERENCE_MF);
    getfem::mesh_fem ref_mf(ref_m); 
    ref_mf.read_from_file(REFERENCE_MF);
    plain_vector ref_U(ref_mf.nb_dof());
    gmm::vecload(REFERENCE_U, ref_U);
    
    cout << "Load reference pressure from "
	 << REFERENCE_MFP << " and " << REFERENCE_P << "\n";
    getfem::mesh_fem ref_mfp(ref_m); 
    ref_mfp.read_from_file(REFERENCE_MFP);
    plain_vector ref_P(ref_mfp.nb_dof());
    gmm::vecload(REFERENCE_P, ref_P);
    
    getfem::mesh_im ref_mim(ref_m);
    getfem::pintegration_method ppi = 
      getfem::int_method_descriptor("IM_TRIANGLE(6)");
    ref_mim.set_integration_method(ref_m.convex_index(), ppi);
    
    plain_vector interp_U(ref_mf.nb_dof());
    getfem::interpolation(p.mf_u(), ref_mf, U, interp_U);


    plain_vector interp_U_error(ref_mf.nb_dof());
    gmm::add(interp_U, gmm::scaled(ref_U, -1.), interp_U_error);
    gmm::vecsave(p.datafilename+".U_map_error", interp_U_error);

    cout << "To ref L2 ERROR on U:"
	 << getfem::asm_L2_dist(ref_mim, ref_mf, interp_U,
				ref_mf, ref_U) << endl;
    
    cout << "To ref H1 ERROR on U:"
	 << getfem::asm_H1_dist(ref_mim, ref_mf, interp_U,
				ref_mf, ref_U) << endl;
    
    plain_vector interp_P(ref_mfp.nb_dof());
    getfem::interpolation(p.mf_pe(), ref_mfp, P, interp_P);

    plain_vector interp_P_error(ref_mfp.nb_dof());
    gmm::add(interp_P, gmm::scaled(ref_P, -1.), interp_P_error);
    gmm::vecsave(p.datafilename+".P_map_error", interp_P_error);
    
    cout << "To ref L2 ERROR on P:"
	 << getfem::asm_L2_dist(ref_mim, ref_mfp, interp_P,
				ref_mfp, ref_P) << endl;
    
    gmm::add(gmm::scaled(interp_U, -1.), ref_U);
    gmm::vecsave(p.datafilename + ".diff_ref", ref_U);
    
  }
  
  return 0; 
}
