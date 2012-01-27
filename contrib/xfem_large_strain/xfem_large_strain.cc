// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2010 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

/**
  Nonlinear elastostatic  crack problem
  This is an example of use of Getfem++.
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
#include "getfem_nonlinear_elastoptim.h" /*Optimization procedure to evaluate the order of singularities*/
#include "gmm/gmm.h"
#include "gmm/gmm_solver_Newton.h"
#include "gmm/gmm_inoutput.h"

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::base_vector;
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::dim_type;
using bgeot::short_type;
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/**************************************************************************/
/*                                                                        */
/*                            Enrichment.                                 */
/*                                                                        */
/**************************************************************************/

scalar_type alpha_md = 0.5;
scalar_type beta_md = 0.5;

struct generic_u_singular_xy_function : public getfem::abstract_xy_function {
    int n;
    virtual scalar_type val(scalar_type x, scalar_type y) const;
    virtual base_small_vector grad(scalar_type x, scalar_type y) const;
    virtual base_matrix hess(scalar_type x, scalar_type y) const;
    generic_u_singular_xy_function(int n_) : n(n_) {}
  };


scalar_type 
generic_u_singular_xy_function::val(scalar_type x, scalar_type y) const {
  scalar_type r = sqrt(x*x + y*y);  
  if (r < 1e-10)  return 0;
  scalar_type theta = atan2(y, x);
  
  if (n <= 0)
    return pow(r, alpha_md-n) * cos(scalar_type(n) * theta * 0.5);
  else
    return pow(r, alpha_md+n) * sin(scalar_type(n) * theta * 0.5);
}


base_small_vector
generic_u_singular_xy_function::grad(scalar_type x, scalar_type y) const {
  scalar_type r = sqrt(x*x + y*y);
  if (r < 1e-10) {
    GMM_WARNING0("Warning, point close to the singularity (r=" << r << ")");
  }
  scalar_type theta = atan2(y,x);
  base_small_vector res(2);
  scalar_type cos_n_2 = cos(scalar_type(n) * theta * 0.5);
  scalar_type sin_n_2 = sin(scalar_type(n) * theta * 0.5);

  if (n <= 0) {
    res[0] = alpha_md * x * cos_n_2 - scalar_type(n) * 0.5 * y * sin_n_2;
    res[1] = alpha_md  * y * cos_n_2 + scalar_type(n) * 0.5 * x * sin_n_2;
  } else {
    res[0] = alpha_md  * x * sin_n_2 + scalar_type(n) * 0.5 * y * cos_n_2;
    res[1] = alpha_md  * y * sin_n_2 - scalar_type(n) * 0.5 * x * cos_n_2;
  }

  gmm::scale(res, pow(r, alpha_md - 2));
  return res;
}


base_matrix generic_u_singular_xy_function::hess(scalar_type, scalar_type)
  const {
  GMM_ASSERT1(false, "To be done !");
}

struct generic_p_singular_xy_function : public getfem::abstract_xy_function {
  int n;
  virtual scalar_type val(scalar_type x, scalar_type y) const;
  virtual base_small_vector grad(scalar_type x, scalar_type y) const;
  virtual base_matrix hess(scalar_type x, scalar_type y) const;
  generic_p_singular_xy_function(int n_) : n(n_) {}
};

scalar_type
generic_p_singular_xy_function::val(scalar_type x, scalar_type y) const {
  scalar_type r = sqrt(x*x + y*y);  
  if (r < 1e-10)  return 0;
  scalar_type theta = atan2(y, x);
  
  if (n <= 0)
    return pow(r, beta_md) * cos(scalar_type(n) * theta * 0.5);
  else
    return pow(r, beta_md) * sin(scalar_type(n) * theta * 0.5);
}


base_small_vector
generic_p_singular_xy_function::grad(scalar_type x, scalar_type y) const {
  scalar_type r = sqrt(x*x + y*y);
  if (r < 1e-10) {
    GMM_WARNING0("Warning, point close to the singularity (r=" << r << ")");
  }
  scalar_type theta = atan2(y, x);
  base_small_vector res(2);
  scalar_type cos_n_2 = cos(scalar_type(n) * theta * 0.5);
  scalar_type sin_n_2 = sin(scalar_type(n) * theta * 0.5);

  if (n <= 0) {
    res[0] = beta_md * x * cos_n_2 - scalar_type(n) * 0.5 * y * sin_n_2;
    res[1] = beta_md  * y * cos_n_2 + scalar_type(n) * 0.5 * x * sin_n_2;
  } else {
    res[0] = beta_md  * x * sin_n_2 + scalar_type(n) * 0.5 * y * cos_n_2;
    res[1] = beta_md  * y * sin_n_2 - scalar_type(n) * 0.5 * x * cos_n_2;
  }

  gmm::scale(res, pow(r, beta_md - 2));
  return res;
}


base_matrix generic_p_singular_xy_function::hess(scalar_type, scalar_type)
  const {
  GMM_ASSERT1(false, "To be done !");
}


///////////*/*/*/*////////////////////////////////////////////////////////////////////////
//                                                                                      //
//    Displacement Singular functions for crack problem  from Stephenson article 1982   //
//                                                                                      //
///////////*/*/*/*////////////////////////////////////////////////////////////////////////
 
struct steph_u_singular_xy_function : public getfem::abstract_xy_function {
    unsigned l;
    virtual scalar_type val(scalar_type x, scalar_type y) const;
    virtual base_small_vector grad(scalar_type x, scalar_type y) const;
    virtual base_matrix hess(scalar_type x, scalar_type y) const;
    steph_u_singular_xy_function(unsigned l_) : l(l_) {}
  };


scalar_type steph_u_singular_xy_function::val(scalar_type x, scalar_type y) const {
  
  scalar_type r = sqrt(x*x + y*y);  
  if (r < 1e-10)  return 0;
  scalar_type theta = atan2(y, x);
  scalar_type res = 0;
    switch(l){
      case 0 : res = r*sin(theta*0.5)*sin(theta*0.5) ; break;
      case 1 : res = sqrt(r)*sin(theta*0.5)          ; break;
      default: GMM_ASSERT2(false, "arg");
    }
  return res;
}


base_small_vector steph_u_singular_xy_function::grad(scalar_type x, scalar_type y) const {
  scalar_type r = sqrt(x*x + y*y);
  
  scalar_type theta = atan2(y, x);
  base_small_vector res(2);
  if (r < 1e-10) {
    GMM_WARNING0("Warning, point close to the singularity (r=" << r << ")");
  }
  
    switch(l){
    case 0 :
      res[0] = sin(theta*0.5)*sin(theta*0.5);
      res[1] = sin(theta*0.5)*cos(theta*0.5);
      break;
    case 1 :
      res[0] = (-1/sqrt(r))*0.5*sin(theta*0.5) ;
      res[1] = 0.5*sqrt(r)*sin(theta*0.5);
      break;
      
      default: GMM_ASSERT2(false, "arg");
    }

  return res;
}


base_matrix steph_u_singular_xy_function::hess(scalar_type, scalar_type)
  const {
  GMM_ASSERT1(false, "To be done !");
}


///////////*/*/*/*////////////////////////////////////////////////////////////////////
//                                                                                  //
//    Pressure Singular functions for crack problem  from Stephenson article 1982   //
//                                                                                  //
///////////*/*/*/*////////////////////////////////////////////////////////////////////
 
struct steph_p_singular_xy_function : public getfem::abstract_xy_function {
    unsigned l;
    virtual scalar_type val(scalar_type x, scalar_type y) const;
    virtual base_small_vector grad(scalar_type x, scalar_type y) const;
    virtual base_matrix hess(scalar_type x, scalar_type y) const;
    steph_p_singular_xy_function(unsigned l_) : l(l_) {}
  };


scalar_type steph_p_singular_xy_function::val(scalar_type x, scalar_type y) const {
  
  scalar_type r = sqrt(x*x + y*y);  
  if (r < 1e-10)  return 0;
  scalar_type theta = atan2(y, x);
  scalar_type res = 0;
    switch(l){
      case 0 : res = 1/r ; break;
      case 1 : res = sin(theta*0.5)/sqrt(r) ; break;
      default: GMM_ASSERT2(false, "arg");
    }
  return res;
}


base_small_vector steph_p_singular_xy_function::grad(scalar_type x, scalar_type y) const {
  scalar_type r = sqrt(x*x + y*y);
  scalar_type theta = atan2(y, x);
  base_small_vector res(2);
  if (r < 1e-10) {
    GMM_WARNING0("Warning, point close to the singularity (r=" << r << ")");
  }
  
    switch(l){
    case 0 :
      res[0] =- 1/(r*r) ;
      res[1] = 0 ;
      break;
    case 1 :
      res[0] = (-0.5*sin(theta*0.5))/(r*sqrt(r)) ;
      res[1] = (0.5*sin(theta*0.5))/(r*sqrt(r));
      break;
      
      default: GMM_ASSERT2(false, "arg");
    }

  return res;
}


base_matrix steph_p_singular_xy_function::hess(scalar_type, scalar_type)
  const {
  GMM_ASSERT1(false, "To be done !");
}





/*
  structure for the elastostatic problem
*/
struct cr_nl_elastostatic_problem {

  enum { BOUNDARY_NUM0 = 0, BOUNDARY_NUM1 = 1, BOUNDARY_NUM2 = 2, BOUNDARY_NUM3 = 3, BOUNDARY_NUM4 = 4, PART_CAL_ERROR = 5,DIRICHLET_BOUND=6, MORTAR_BOUNDARY_IN=42, MORTAR_BOUNDARY_OUT=43};
  getfem::mesh mesh;         /* the mesh */
  getfem::level_set ls;      /* The two level sets defining the crack.       */
  getfem::mesh_level_set mls;       /* the integration methods for cutted element.    */
  getfem::mesh_im_level_set mim;
  getfem::mesh_fem mf_pre_u, mf_pre_mortar;
  getfem::mesh_fem mf_mult, mf_mult_p;

  getfem::mesh_fem_level_set mfls_u,mfls_mortar;
  getfem::mesh_fem_global_function mf_sing_u;
  
  getfem::mesh_fem mf_partition_of_unity;
  getfem::mesh_fem_product mf_product;
  getfem::mesh_fem_sum mf_u_sum/*, mf_u_sumCE*/;
  
  getfem::mesh_fem mf_pre_p; /* mesh_fem for the pressure for mixed form     */
  getfem::mesh_fem_level_set mfls_p;   /* mesh_fem for the pressure enriched with H.   */
  getfem::mesh_fem_global_function mf_sing_p;
  getfem::mesh_fem_product mf_product_p;
  getfem::mesh_fem_sum mf_p_sum /*, mf_p_sumCE*/;

  base_small_vector cracktip;
  
  scalar_type pr1, pr2, pr3, AMP_LOAD_X, AMP_LOAD_Y, nb_step;
 


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
  getfem::mesh_fem& mf_u() { return mf_u_sum; }
  getfem::mesh_fem& mf_pe() { return mf_p_sum; }
  getfem::mesh_fem mf_rhs;
  
 
  scalar_type residual;        /* max residual for the iterative solvers         */
  bool mixed_pressure;
  //bool sing_search;
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
  cr_nl_elastostatic_problem(void) : ls(mesh, 1, true), mls(mesh), mim(mls), 
				     mf_pre_u(mesh), mf_pre_mortar(mesh), mf_mult(mesh), mf_mult_p(mesh),
				     mfls_u(mls, mf_pre_u), mfls_mortar(mls, mf_pre_mortar),	
				     mf_sing_u(mesh),
				     mf_partition_of_unity(mesh),
				     mf_product(mf_partition_of_unity, mf_sing_u),
				     mf_u_sum(mesh),
				     mf_pre_p(mesh), mfls_p(mls, mf_pre_p), 
				     mf_sing_p(mesh), mf_product_p(mf_partition_of_unity, mf_sing_p), 
				     mf_p_sum(mesh),
				     mf_us(mesh),  mf_rhs(mesh) {}
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
   and integration methods and selects the boundaries.

   (this is boilerplate code, not very interesting)
 */
void cr_nl_elastostatic_problem::init(void) {

  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE", "Mesh type");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE", "FEM name");
  std::string FEM_TYPE_P  = PARAM.string_value("FEM_TYPE_P", "FEM name for the pressure");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION", "Name of integration method");
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION", "Name of simplex integration method");
  std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION", "Singular integration");
  enrichment_option = enrichment_option_enum(PARAM.int_value("ENRICHMENT_OPTION", "Enrichment option"));
  /* Lecture des parametres */
  pr1 = PARAM.real_value("P1", "First Elastic coefficient");
  pr2 = PARAM.real_value("P2", "Second Elastic coefficient");
  pr3 = PARAM.real_value("P3", "Third Elastic coefficient");
  AMP_LOAD_X = PARAM.real_value("AMP_LOAD_X", "Amp load x");
  AMP_LOAD_Y = PARAM.real_value("AMP_LOAD_Y", "Amp load y");
  nb_step  = PARAM.real_value("nb_step", "nb_step");
  reference_test = int(PARAM.int_value("REFERENCE_TEST", "Reference test"));
  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;
  
 /* Affichage */

  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "FEM_TYPE_P  =" << FEM_TYPE_P <<"\n";  
  cout << "INTEGRATION=" << INTEGRATION << "\n";
  cout << "AMP_LOAD_X = " << AMP_LOAD_X << "\n";
  cout << "AMP_LOAD_Y = " << AMP_LOAD_Y << "\n";
  cout << "Reference_test " << AMP_LOAD_Y << "\n";

  /* First step : build the mesh */

  bgeot::pgeometric_trans pgt = bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  std::vector<size_type> nsubdiv(N);
  
  long int  NX_aff = PARAM.int_value("NX", "Nomber of space steps ");
  std::fill(nsubdiv.begin(),nsubdiv.end(),NX_aff);
  cout << "Number of space step NX= " << NX_aff << "\n";
  
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt, PARAM.int_value("MESH_NOISED") != 0);
  
  base_small_vector tt(N); tt[1] = -0.5;
  mesh.translation(tt);
  
  cracktip.resize(2); // Coordonnee du fond de fissure
  cracktip[0] = 0.5;
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

  mesh.write_to_file("le_maillage.mesh");
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

  /* set the finite element and integration method on the singular mesh */

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

  mixed_pressure = (PARAM.int_value("MIXED_PRESSURE","Mixed version or not.") != 0);
  //sing_search = (PARAM.int_value("SINGULAR_SEARCH","Singular search or not.") != 0);
  dir_with_mult = unsigned(PARAM.int_value("DIRICHLET_VERSION", "Version of Dirichlet"));

  if (mixed_pressure) {
        getfem::pfem pf_p = 
        getfem::fem_descriptor(FEM_TYPE_P);
	mf_pre_p.set_finite_element(mesh.convex_index(), pf_p);
	mf_mult_p.set_finite_element(mesh.convex_index(), pf_p);
    
  }  




  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM"
		". In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }
  
  /* set the finite element on mf_coef. Here we use a very simple element
   *  since the only function that need to be interpolated on the mesh_fem 
   * is f(x)=1 ... */
  

  /* set boundary conditions
   * (Neuman on the upper face, Dirichlet elsewhere) */

/*******************************************************************************/
/*         Select the part of convex where the error will be calculation       */
/*******************************************************************************/


  //for(getfem::mr_visitor it_err(mesh.convex_index()); !it_err.finished(); ++it_err){   
  
  //base_node un1_cal_err = mesh.normal_of_face_of_convex(it_err.cv(), it_err.f());
  //un1_cal_err /= gmm::vect_norm2(un1_cal_err);
  //if (un1_cal_err[1]> -0.5) mesh.region(PART_CAL_ERROR).add(it_err.cv());   
  //
  //    }



  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::mr_visitor it(border_faces); !it.finished(); ++it) {
    assert(it.is_face());
    base_node un = mesh.normal_of_face_of_convex(it.cv(), it.f());
    un /= gmm::vect_norm2(un);

    if (un[0]  > 0.5) mesh.region(BOUNDARY_NUM1).add(it.cv(), it.f());
    if (un[1]  > 0.5) mesh.region(BOUNDARY_NUM2).add(it.cv(), it.f());
    if (un[0]  < -0.5) mesh.region(BOUNDARY_NUM3).add(it.cv(), it.f());
    if (un[1]  < -0.5) mesh.region(BOUNDARY_NUM4).add(it.cv(), it.f());   

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
/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool cr_nl_elastostatic_problem::solve(plain_vector &U, plain_vector &P) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  // size_type nb_dof_mult = mf_mult.nb_dof();
  size_type N = mesh.dim();
  ls.reinit();
  size_type law_num = PARAM.int_value("LAW");
  size_type line_search_version = PARAM.int_value("line_search_version");
  size_type Pseudo_Potential = PARAM.int_value("Pseudo_Potential");
  std::cout<<"law num  "<< law_num << endl;
  
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

  //  bool load_global_fun = GLOBAL_FUNCTION_MF.size() != 0;


  cout << "Setting up the singular functions for the enrichment\n";

  // size_type nb_enr_func_u = size_type(PARAM.int_value("NB_ENR_FUNC_U", "Number of Enriched function for u"));
  // size_type nb_enr_func_p = size_type(PARAM.int_value("NB_ENR_FUNC_P", "Number of Enriched function for p"));

    //std::vector<getfem::pglobal_function> vfunc(2*nb_enr_func_u);
    //std::vector<getfem::pglobal_function> vfunc_p(2*nb_enr_func_p);
  std::vector<getfem::pglobal_function> vfunc(2);
  std::vector<getfem::pglobal_function> vfunc_p(2);

  std::cout << "Using default singular functions\n";
  for (unsigned i = 0; i < vfunc.size(); ++i){
    
    // getfem::abstract_xy_function *s;
//     if (i < nb_enr_func_u){
//      s = new generic_u_singular_xy_function(i+1);
//     }
//     else{
//      s = new generic_u_singular_xy_function(-(i+1));
//       }

    getfem::abstract_xy_function *s = 
        new steph_u_singular_xy_function(i);

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
    

   // getfem::abstract_xy_function *sp;
   //        if (i < nb_enr_func_p)
   //          sp = new generic_p_singular_xy_function(i+1);
   //        else
   //          sp = new generic_p_singular_xy_function(-(i+1));
    getfem::abstract_xy_function *sp = 
        new steph_p_singular_xy_function(i);
                 
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
  
  //////////////////////////////////////////////////////////////////////
  //                                                                  //
  //     find the dofs on the upper right and lower right corners     //
  //                                                                  //
  //////////////////////////////////////////////////////////////////////

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

  /*******************************************/
  /*                                         */
  /*         choose the material law         */
  /*                                         */
  /*******************************************/

  getfem::abstract_hyperelastic_law *pl1 = 0, *pl = 0;
  switch (law_num) {
    case 0:
    case 1: pl1 = new getfem::SaintVenant_Kirchhoff_hyperelastic_law(); break;
    case 2: pl1 = new getfem::Ciarlet_Geymonat_hyperelastic_law(); break;
    case 3: pl1 = new getfem::Mooney_Rivlin_hyperelastic_law(); break;
    default: GMM_ASSERT1(false, "no such law");
  }

  if (N == 2) {
   cout << "2D plane strain hyper-elasticity\n";
   pl = new getfem::plane_strain_hyperelastic_law(pl1);
   } else pl = pl1;

  pr.resize(pl->nb_params());
  cout << "parametre du loi de comportements "<< pr << endl ;

  //pl1->test_derivatives(3, 5e-9, pr);


  getfem::model model;
  
  // Main unknown of the problem (displacement).

   model.add_fem_variable("u", mf_u());

  // model.add_fem_variable("u",  error_cal);

  // Nonlinear elasticity brick
  model.add_initialized_fixed_size_data("params", pr);
  getfem::add_nonlinear_elasticity_brick(model, mim, "u", *pl, "params");

  // Incompressibility

  if (mixed_pressure && (law_num == 1 || law_num == 3)) {
    cout << "mixed pressure <|-------------------------|>" << endl;
    model.add_fem_variable("p", mf_pe());
    getfem::add_nonlinear_incompressibility_brick(model, mim, "u", "p", size_type(-1));
  }

 // Defining the Neumann condition right hand side.
  plain_vector F_Neumann(nb_dof_rhs * N);
 // Neumann condition brick.
   
  for(size_type i = 0; i < F_Neumann.size(); i=i+N) F_Neumann[i] = AMP_LOAD_X;
  for(size_type i = 1; i < F_Neumann.size(); i=i+N) F_Neumann[i] = AMP_LOAD_Y;
  
   
   model.add_initialized_fem_data("NeumannData", mf_rhs, F_Neumann );

  getfem::add_source_term_brick (model, mim, "u", "NeumannData", BOUNDARY_NUM2);

    //Dirichlet condition brick

    plain_vector F_Diri(nb_dof_rhs * N);
    
    // for(size_type i = 0; i < F_Diri.size(); i=i+N) F_Diri[i] = 0;
    // for(size_type i = 1; i < F_Diri.size(); i=i+N) F_Diri[i] = 0;
  
    // gmm::resize(F_Diri, nb_dof_mult);
    
    model.add_initialized_fem_data("Dirichletdata", mf_rhs, F_Diri);
    if (PARAM.int_value("DIRICHLET_VERSION") == 0)
    add_Dirichlet_condition_with_multipliers (model, mim, "u", mf_u(), BOUNDARY_NUM4, "Dirichletdata");
    else
    add_Dirichlet_condition_with_penalization (model, mim, "u", 1E15, BOUNDARY_NUM4, "Dirichletdata");

//     ////////////////////////////////////////// 
//     //                                      //
//     //     Symetrie condition               // 
//     //                                      //
//     //////////////////////////////////////////

    model.add_initialized_fixed_size_data("Dirichletsymdata",
					  plain_vector(N, 0.0));
    plain_vector HH(N*N); HH[0] = 1.0;
    model.add_initialized_fixed_size_data("Hdata", HH); 

    getfem::add_generalized_Dirichlet_condition_with_multipliers(model, mim, "u", 1, BOUNDARY_NUM3, "Dirichletsymdata", "Hdata");


  


//     getfem::partial_mesh_fem symetrie_NUM3(mf_u());
//     dal::bit_vector symetrie_dofs_NUM3 = mf_u().basic_dof_on_region(BOUNDARY_NUM3);
//     cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
//     cout << "Size of vector  :" << symetrie_dofs_NUM3.card() << endl;
//     cout << "Basic dof on region BOUNDARY_NUM3  :" << symetrie_dofs_NUM3 << endl;
//     cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
//     cout << "test  :" << symetrie_dofs_NUM3[1] << endl;
//     cout << "test  :" << symetrie_dofs_NUM3[2] << endl;
//     GMM_ASSERT1(N==2, "To be corrected for 3D computation");

//     sparse_matrix BB((symetrie_dofs_NUM3.size()/2), mf_u().nb_dof());
    
   
//     // for (size_type i = 0; i < mf_u().nb_basic_dof(); i+=N) 
//     //    {}
//     for (size_type j = 0; j < symetrie_dofs_NUM3.size()+1 ; j+=N ){
//     for (size_type i = 0; i < (symetrie_dofs_NUM3.size()/2)+1 ; ++i )
//       {
//         BB(j,symetrie_dofs_NUM3[i])=1.0;
// 	  }}
   
//     // BB (0, icorner1)   = 1.0;
//     // BB (1, icorner1+1) = 1.0;
//     // BB (2, icorner2)   = 1.0;

//   // cout << "matrice BB   :" << BB << endl;
//   // std::vector<scalar_type> LRH(3);
//   // model.add_fixed_size_variable("dir", 3);
//   // getfem::add_constraint_with_multipliers(model, "u", "dir", BB, LRH);


 
  gmm::iteration iter(residual, 1, 40000);

  /* prepare the export routine for OpenDX (all time steps will be exported) 
     (can be viewed with "dx -edit nonlinear_elastostatic.net")
  */
  getfem::dx_export exp(datafilename + ".dx",
			PARAM.int_value("VTK_EXPORT")==1);
  getfem::stored_mesh_slice sl; sl.build(mesh, getfem::slicer_boundary(mesh),8); 
  exp.exporting(sl,true); exp.exporting_mesh_edges();
  //exp.begin_series("deformationsteps");
  exp.write_point_data(mf_u(), U, "stepinit"); 
  exp.serie_add_object("deformationsteps");

  GMM_ASSERT1(!mf_rhs.is_reduced(), "To be adapted for reduced mesh_fems");
    
  getfem::simplest_newton_line_search simls;
  getfem::default_newton_line_search dlnrs;
  getfem::systematic_newton_line_search sylnrs;
  getfem::quadratic_newton_line_search qdlnrs;

//  simplest_newton_line_search 1 *** default_newton_line_search 2 *** systematic_newton_line_search 3
 cout << "line search value" <<line_search_version<< endl;
 switch (line_search_version){
     
    case 1:{
      if (Pseudo_Potential == 1){
	bool with_pseudo_potential = true;
        //gmm::simplest_newton_line_search simls;
	getfem::standard_solve(model, iter, getfem::default_linear_solver<getfem::model_real_sparse_matrix,
			       getfem::model_real_plain_vector>(model), simls, with_pseudo_potential);
	cout << "/******************/" << endl;
	cout << "/* Enery Criteria */" << endl;
	cout << "/******************/" << endl;
      }
      getfem::standard_solve(model, iter,getfem::default_linear_solver<getfem::model_real_sparse_matrix,
			       getfem::model_real_plain_vector>(model), simls );
      cout << "=============================== " << endl;
      cout << "= Simplest_newton_line_search = " << endl;
      cout << "=============================== " << endl;
      cout << "alpha_md valeur de l'ordre de singularite =  " << alpha_md << endl;
    }break;

    case 2:{
      if (Pseudo_Potential == 1){
	bool with_pseudo_potential = true;
	getfem::standard_solve(model, iter, getfem::default_linear_solver<getfem::model_real_sparse_matrix,
			       getfem::model_real_plain_vector>(model), dlnrs, with_pseudo_potential);
	cout << "/******************/" << endl;
	cout << "/* Enery Criteria */" << endl;
	cout << "/******************/" << endl;
      }
      getfem::standard_solve(model, iter, getfem::default_linear_solver<getfem::model_real_sparse_matrix,
			       getfem::model_real_plain_vector>(model), dlnrs);
      cout << "================================ " << endl;
      cout << "=  Default_newton_line_search  = " << endl;
      cout << "================================ " << endl;   
      cout << "alpha_md valeur de l'ordre de singularite =  " << alpha_md << endl;
    }break;
    

    case 3: {
      if (Pseudo_Potential == 1){
	bool with_pseudo_potential = true;
	getfem::standard_solve(model, iter,getfem::default_linear_solver<getfem::model_real_sparse_matrix,
			       getfem::model_real_plain_vector>(model), sylnrs,  with_pseudo_potential);
	cout << "/******************/" << endl;
	cout << "/* Enery Criteria */" << endl;
	cout << "/******************/" << endl;
      }
      getfem::standard_solve(model, iter, getfem::default_linear_solver<getfem::model_real_sparse_matrix,
			       getfem::model_real_plain_vector>(model), sylnrs); 
      cout << "=================================== " << endl;
      cout << "=  Systematic_newton_line_search  = " << endl;
      cout << "=================================== " << endl;
      cout << "alpha_md valeur de l'ordre de singularite =  " << alpha_md << endl;
    }break;

    case 4: {
      if (Pseudo_Potential == 1){
	bool with_pseudo_potential = true;
	getfem::standard_solve(model, iter,getfem::default_linear_solver<getfem::model_real_sparse_matrix,
			       getfem::model_real_plain_vector>(model), qdlnrs,  with_pseudo_potential);
	cout << "/******************/" << endl;
	cout << "/* Enery Criteria */" << endl;
	cout << "/******************/" << endl;
      }
      getfem::standard_solve(model, iter, getfem::default_linear_solver<getfem::model_real_sparse_matrix,
			       getfem::model_real_plain_vector>(model), qdlnrs); 
      cout << "================================== " << endl;
      cout << "=  Quadratic_newton_line_search  = " << endl;
      cout << "================================== " << endl;
      cout << "alpha_md valeur de l'ordre de singularite =  " << alpha_md << endl;
    }break;




    default: GMM_ASSERT1(false, "No such line search");
    }

  
  

   // Solution extraction
   gmm::copy(model.real_variable("u"), U);
   if (mixed_pressure && (law_num == 1 || law_num == 3)) {
     gmm::copy(model.real_variable("p"), P);}

  
  return (iter.converged());

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

  cr_nl_elastostatic_problem p;
  p.PARAM.read_command_line(argc, argv);
  p.init();
  p.mesh.write_to_file(p.datafilename + ".mesh");
  //getfem::mesh_region &rg_error_calc = p.mesh.region(1);
  //getfem::mesh_region &rg2 = p.mesh.region(2);
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
  
  
  getfem::mesh_fem mf_refined(mcut_refined, dim_type(Q));
  mf_refined.set_classical_discontinuous_finite_element(2,0.001);
  plain_vector W(mf_refined.nb_dof());
 
  
  if (p.PARAM.int_value("VTK_EXPORT")) {
    
     getfem::mesh_fem mf_refined_p(mcut_refined, 1);
    
    mf_refined_p.set_classical_discontinuous_finite_element(2, 0.001);
  
    plain_vector PPW(mf_refined_p.nb_dof());
   
   
    getfem::interpolation(p.mf_pe(), mf_refined_p, P, PPW);

    
    cout << "export to " << p.datafilename + ".vtk" << "..\n";
    getfem::vtk_export exp(p.datafilename + ".vtk",
			   p.PARAM.int_value("VTK_EXPORT")==1);
    
    exp.exporting(mf_refined); 
    
    //exp.write_point_data(mf_refined_vm, DN, "error");
    getfem::interpolation(p.mf_pe(), mf_refined_p, P, PPW);
    
    exp.write_point_data(mf_refined, W, "elastostatic_displacement");
    exp.write_point_data(mf_refined_p, PPW, "Pressure");

    base_node line_x0(0.70001,0);
    base_small_vector line_dir(0, 0.5001);

    
    cout << "export done, you can view the data file with (for example)\n"
      "mayavi -d " << p.datafilename << ".vtk -f  "
      "WarpVector -m BandedSurfaceMap -m Outline\n";
  }
    
  if(p.PARAM.int_value("ERROR_TO_REF_SOL") == 1){
    cout << "Computing error with respect to a reference solution..." << endl;

   
    std::string REFERENCE_MF  = "TestEnrichemen**tPressureStephNX40.meshfem";
    std::string REFERENCE_U   = "TestEnrichemen**tPressureStephNX40.U";
    std::string REFERENCE_MFP = "TestEnrichemen**tPressureStephNX40.p_meshfem";
    std::string REFERENCE_P   = "TestEnrichemen**tPressureStephNX40.P";
                                                                  
    cout << "Load reference displacement from "
	 << REFERENCE_MF << " and " << REFERENCE_U << "\n";

    getfem::mesh ref_m ; 
    ref_m.read_from_file(REFERENCE_MF);
        
    
    enum { Part_cal_error = 82};
      for(dal::bv_visitor i(ref_m.convex_index()); !i.finished(); ++i){
	for(size_type j=0; j < 3; ++j){
        base_node cord_err = ref_m.points_of_convex(i)[j];
	//cout << "coordonne-erreur " << cord_err[1] << endl;
	 if (cord_err[1] > -0.39) ref_m.region(Part_cal_error).add(i);
	}   
        }


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
	 << getfem::asm_L2_dist(ref_mim, ref_mf, interp_U, ref_mf, ref_U, Part_cal_error) << endl;
    
    cout << "To ref H1 ERROR on U:"
	 << getfem::asm_H1_dist(ref_mim, ref_mf, interp_U,
				ref_mf, ref_U, Part_cal_error ) << endl;
    
    plain_vector interp_P(ref_mfp.nb_dof());
    getfem::interpolation(p.mf_pe(), ref_mfp, P, interp_P, Part_cal_error);

    plain_vector interp_P_error(ref_mfp.nb_dof());
    gmm::add(interp_P, gmm::scaled(ref_P, -1.), interp_P_error);
    gmm::vecsave(p.datafilename+".P_map_error", interp_P_error);
    
    cout << "To ref L2 ERROR on P:"
	 << getfem::asm_L2_dist(ref_mim, ref_mfp, interp_P,
				ref_mfp, ref_P, Part_cal_error) << endl;
    cout << "/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/ " << endl;
    cout << "/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/ " << endl;
    cout << "                                                                                            " << endl;

    cout << "Norme of displacement error vector %%% asm_L2_norm %%% : " 
	 << getfem::asm_L2_norm(ref_mim, ref_mf, interp_U_error, Part_cal_error)
         << endl;
    cout << "Norme of displacement error vector %%% asm_H1_norm %%% : " 
	 << getfem::asm_H1_norm(ref_mim, ref_mf, interp_U_error, Part_cal_error)
         << endl;
    cout << "Norme of pressure error vector     %%% asm_L2_norm %%% : " 
	 << getfem::asm_L2_norm(ref_mim, ref_mfp, interp_P_error, Part_cal_error)
         << endl;
    gmm::add(gmm::scaled(interp_U, -1.), ref_U);
    gmm::vecsave(p.datafilename + ".diff_ref", ref_U);
    
  }


  return 0; 
}

