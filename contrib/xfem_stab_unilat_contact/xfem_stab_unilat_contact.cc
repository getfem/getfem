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
 * Goal : stabilization of unilateral contact with Xfem.
 *
 * Research program.
 */

#include "getfem/getfem_assembling.h" /* import assembly methods      */
#include "getfem/getfem_export.h"     /* export functions             */
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_mesh_fem_product.h"
#include "getfem/getfem_mesh_fem_global_function.h"
#include "getfem/getfem_partial_mesh_fem.h"
#include "getfem/getfem_Coulomb_friction.h"
#include "getfem/getfem_import.h"
#include "getfem/getfem_inter_element.h"
#include "gmm/gmm.h"
#include "getfem/getfem_mesh_fem_sum.h"
#include "gmm/gmm_inoutput.h"


/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_vector;
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::short_type;  /* = short */
using bgeot::size_type;   /* = unsigned long */
using bgeot::dim_type;
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

typedef gmm::row_matrix<sparse_vector> sparse_row_matrix;


/* level set unit normal */

template<typename VECT1> class level_set_unit_normal 
  : public getfem::nonlinear_elem_term {
  const getfem::mesh_fem &mf;
  std::vector<scalar_type> U;
  size_type N;
  base_matrix gradU;
  bgeot::base_vector coeff;
  bgeot::multi_index sizes_;
public:
  level_set_unit_normal(const getfem::mesh_fem &mf_, const VECT1 &U_) 
    : mf(mf_), U(mf_.nb_basic_dof()), N(mf_.linked_mesh().dim()),
      gradU(1, N) {
    sizes_.resize(1); sizes_[0] = short_type(N);
    mf.extend_vector(U_, U);
  }
  const bgeot::multi_index &sizes() const {  return sizes_; }
  virtual void compute(getfem::fem_interpolation_context& ctx,
		       bgeot::base_tensor &t) {
    size_type cv = ctx.convex_num();
    coeff.resize(mf.nb_basic_dof_of_element(cv));
    gmm::copy
      (gmm::sub_vector(U,gmm::sub_index(mf.ind_basic_dof_of_element(cv))),
       coeff);
    ctx.pf()->interpolation_grad(ctx, coeff, gradU, 1);
    scalar_type norm = gmm::vect_norm2(gmm::mat_row(gradU, 0));
    for (size_type i = 0; i < N; ++i) t[i] = gradU(0, i) / norm;
  }
};

/***********************************/
/* asembling stabilised mixed term */
/***********************************/

template<class MAT, class VECT>
void asm_stabilization_mixed_term
(const MAT &RM_, const getfem::mesh_im &mim, const getfem::mesh_fem &mf, 
 const getfem::mesh_fem &mf_mult, getfem::level_set &ls,const VECT &LAMBDA, const VECT &MU,
 const getfem::mesh_region &rg = getfem::mesh_region::all_convexes()) {
  MAT &RM = const_cast<MAT &>(RM_);

  level_set_unit_normal<std::vector<scalar_type> >
    nterm(ls.get_mesh_fem(), ls.values());

  getfem::generic_assembly assem("lambda=data$1(1); mu=data$2(1);"
                                 //"t=comp(Base(#2).vGrad(#1).NonLin(#3));"
                                 "t=comp(Base(#2).NonLin(#3).vGrad(#1).NonLin(#3))"
				 "M(#2, #1)+= t(:,i,:,j,j,i).lambda(1)"
				 "+t(:,i,:,i,j,j).2.mu(1)");
  assem.push_mi(mim);
  assem.push_mf(mf);
  assem.push_mf(mf_mult);
  assem.push_mf(ls.get_mesh_fem());
  assem.push_data(LAMBDA);
  assem.push_data(MU);
  assem.push_mat(RM);
  assem.push_nonlinear_term(&nterm);



 ls.set_shift(1e-7);
 assem.assembly(rg);
 ls.set_shift(-1e-7);
 assem.assembly(rg);
 ls.set_shift(0.);
}

/***********************************************/
/* asembling stabilization symetric term       */
/***********************************************/

template<class MAT, class VECT>
void asm_stabilization_symm_term
(const MAT &RM_, const getfem::mesh_im &mim, const getfem::mesh_fem &mf, 
 getfem::level_set &ls, const VECT &LAMBDA, const VECT &MU,
 const getfem::mesh_region &rg = getfem::mesh_region::all_convexes()) {
  MAT &RM = const_cast<MAT &>(RM_);

  level_set_unit_normal<std::vector<scalar_type> >
    nterm(ls.get_mesh_fem(), ls.values());

  getfem::generic_assembly
    assem("lambda=data$1(1); mu=data$2(1);"
          "t=comp(NonLin(#3).vGrad(#1).NonLin(#3).NonLin(#3).vGrad(#1).NonLin(#3));"
	  "M(#1, #1)+= sym(t(i,:,j,j,i,k,:,l,l,k).lambda(1).lambda(1)"
          "+t(i,:,j,j,i,k,:,k,l,l).2.lambda(1).mu(1)"
	   "+t(i,:,i,j,j,k,:,l,l,k).2.lambda(1).mu(1)"
          "+t(i,:,i,j,j,k,:,k,l,l).4.mu(1).mu(1)");
  assem.push_mi(mim);
  assem.push_mf(mf);
  assem.push_mf(ls.get_mesh_fem());
  assem.push_data(LAMBDA);
  assem.push_data(MU);
  assem.push_mat(RM);
  assem.push_nonlinear_term(&nterm);

 ls.set_shift(1e-7);
 assem.assembly(rg);
 ls.set_shift(-1e-7);
 assem.assembly(rg);
 ls.set_shift(0.);
}

/**************************************************************************/
/*structure of unilateral contact proble                                  */
/**************************************************************************/

struct unilateral_contact_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN1_BOUNDARY_NUM = 1, NEUMANN2_BOUNDARY_NUM=2, NEUMANN3_BOUNDARY_NUM=3, NEUMANN4_BOUNDARY_NUM=4};
  getfem::mesh mesh;  /* the mesh */
  getfem::level_set ls;      /* The two level sets defining the crack.       */
  getfem::mesh_level_set mls, mls_bound;       /* two meshs level set                    */
  getfem::mesh_im_level_set mim, mimbound;    /* the integration methods.              */
  getfem::mesh_fem mf_pre_u;
  getfem::mesh_fem mf_mult_dir, mf_mult_cont;
  getfem::mesh_fem_level_set mfls_u;
  getfem::mesh_fem_global_function mf_sing_u;
  
  getfem::mesh_fem mf_partition_of_unity;
  getfem::mesh_fem_product mf_product;
  getfem::mesh_fem_sum mf_u_sum;
  
  getfem::mesh_fem mf_pre_cont; /* mesh_fem for the contact multiplier     */
  getfem::mesh_fem_level_set mfls_cont;   /* mesh_fem for the multiplier contact enriched with H.   */
  getfem::mesh_fem_sum mf_cont_sum;


  base_small_vector cracktip;


  getfem::mesh_fem& mf_u() { return mf_u_sum; }
  getfem::mesh_fem& mf_cont() { return mf_cont_sum; }
  
  scalar_type mu, lambda;    /* Lame coeff                   */

  int dgr ;          /* Order of enrichement for u */
 

  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  
  
  int bimaterial;           /* For bimaterial interface unilateral contact problem */

  double lambda_up, lambda_down, mu_up, mu_down;  /*Lame coeff for bimaterial case*/
 
  scalar_type residual;      /* max residual for the iterative solvers      */
  bool stabilized_problem;
  unsigned dir_with_mult;
  scalar_type cutoff_radius, cutoff_radius1, cutoff_radius0, enr_area_radius;
  
  size_type cutoff_func;

  typedef enum { NO_ENRICHMENT=0, 
		 FIXED_ZONE=1, 
		 GLOBAL_WITH_CUTOFF=2 } enrichment_option_enum;
  enrichment_option_enum enrichment_option;

  scalar_type h, cont_gamma0, R; // mesh parameter
  

  std::string datafilename;
  
  int reference_test;

  bgeot::md_param PARAM;

  bool solve(plain_vector &U, plain_vector &LAMBDA);
  void init(void);
  unilateral_contact_problem(void) : ls(mesh, 1, true), mls(mesh), mls_bound(mesh), mim(mls), mimbound(mls_bound),
		                     mf_pre_u(mesh), mf_mult_dir(mesh), mf_mult_cont(mesh),
				     mfls_u(mls, mf_pre_u),
				     mf_sing_u(mesh),
				     mf_partition_of_unity(mesh),
				     mf_product(mf_partition_of_unity, mf_sing_u),
				     mf_u_sum(mesh), mf_pre_cont(mesh), mfls_cont(mls, mf_pre_cont), 				      
				     mf_cont_sum(mesh),
				     /*mf_pe(mesh),*/ mf_rhs(mesh)
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


void  unilateral_contact_problem::init(void) {
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string FEM_TYPE_CONT  = PARAM.string_value("FEM_TYPE_cont","FEM name mult contact viriable");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
					 "Name of simplex integration method");
  std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION");

  enrichment_option = enrichment_option_enum(PARAM.int_value("ENRICHMENT_OPTION",
							     "Enrichment option"));
  cout << "MESH_TYPE="      << MESH_TYPE         << "\n";
  cout << "FEM_TYPE="       << FEM_TYPE          << "\n";
  cout << "FEM_TYPE_CONT="  << FEM_TYPE_CONT     << "\n";
  cout << "INTEGRATION="    << INTEGRATION       << "\n";


 dgr = int(PARAM.int_value("dgr", "Degr%G�%@ d'enrichissement en u"));


 //  bimaterial = int(PARAM.int_value("BIMATERIAL", "bimaterial interface crack"));
  
  mu    = PARAM.real_value("MU", "Lame coefficient mu"); 
  lambda =PARAM.real_value("Lamda", "Lame coefficient lambda");

  // if (bimaterial == 1){
  //   mu_up = PARAM.real_value("MU_UP", "Lame coefficient mu"); 
  // mu_down = PARAM.real_value("MU_DOWN", "Lame coefficient mu"); 
  //Lambda_up=PARAM.real_value("LAMBDA_UP", "Lame coefficient Lambda_up");
  //Lambda_down=PARAM.real_value("LAMBDA_DOWN", "Lame coefficient Lambda_down");
  // }
  

  /* First step : build the mesh */
  bgeot::pgeometric_trans pgt = 
  bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Nomber of space steps "));
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			    PARAM.int_value("MESH_NOISED") != 0);
  base_small_vector tt(N); tt[1] = -0.5;
  mesh.translation(tt); 
  
  cracktip.resize(2); // Cordinate of cracktip
  cracktip[0] = 0.5;
  cracktip[1] = 0.;

  scalar_type refinement_radius;
  refinement_radius
    = PARAM.real_value("REFINEMENT_RADIUS", "Refinement Radius");
  size_type refinement_process;
  refinement_process
    = PARAM.int_value("REFINEMENT_PROCESS", "Refinement process");
  
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

  h = mesh.minimal_convex_radius_estimate();
  cout << "h = " << h << endl;

 cont_gamma0 = PARAM.real_value("CONTACT_GAMMA0",
				  "Stabilization parameter for "
				  "contact condition");
 R = PARAM.real_value( "R", "augmented parameter");	  

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL");
  if (residual == 0.) residual = 1e-10;
  enr_area_radius = PARAM.real_value("RADIUS_ENR_AREA",
				     "radius of the enrichment area");

  reference_test = int(PARAM.int_value("REFERENCE_TEST", "Reference test")); 


  cutoff_func = PARAM.int_value("CUTOFF_FUNC", "cutoff function");

  cutoff_radius = PARAM.real_value("CUTOFF", "Cutoff");
  cutoff_radius1 = PARAM.real_value("CUTOFF1", "Cutoff1");
  cutoff_radius0 = PARAM.real_value("CUTOFF0", "Cutoff0");
  mf_u().set_qdim(dim_type(N));
  

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pfem pf_mult_cont = 
    getfem::fem_descriptor(FEM_TYPE_CONT);
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
  mf_mult_dir.set_finite_element(mesh.convex_index(), pf_u);
  mf_mult_dir.set_qdim(dim_type(N));
  mf_partition_of_unity.set_classical_finite_element(1);
  
  /*******************************************************************************/
// Integration method on the boudary


  int intbound = getfem::mesh_im_level_set::INTEGRATE_BOUNDARY;
  mimbound.set_integration_method(intbound, ppi);
  mimbound.set_simplex_im(simp_ppi, sing_ppi);
 // mimbound(mls, intbound, getfem::int_method_descriptor(SIMPLEX_INTEGRATION));
 // mimbound.set_integration_method(mesh.convex_index(),
//				      getfem::int_method_descriptor(INTEGRATION);
  mimbound.adapt();


   stabilized_problem =
    (PARAM.int_value("STABILIZED_PROBLEM"," stabilized_problem or not.") != 0);
  dir_with_mult = unsigned(PARAM.int_value("DIRICHLET_VERSION",
					   "Version of Dirichlet"));
  // if ( stabilized_problem ) {
        //getfem::pfem pf_mult_cont = getfem::fem_descriptor(FEM_TYPE_CONT);

  mf_pre_cont.set_finite_element(mesh.convex_index(), pf_mult_cont);
  mf_mult_cont.set_finite_element(mesh.convex_index(), pf_mult_cont);
   
  

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
    
}//end intialisation


/*****************************************************************/

base_small_vector ls_function(const base_node P, int num = 0) {
  scalar_type x = P[0], y = P[1], x0=0, y0=0, phi, ll;
  base_small_vector res(2), cracktip(2), t(2), xx(2);
  cracktip[0]=0.5; cracktip[1]= 0.;
  switch (num) {

    case 0: {
      xx[0]= P[0]; xx[1]=P[1];
      t[0]= cracktip[0] - x0;
      t[1]= cracktip[1] - y0;
      t= t * (1./ gmm::vect_norm2(t));
      ll= gmm::sqr(t[0]* t[0] + t[1]*t[1]);
      phi= t[0]*y - t[1]*x + x0*cracktip[1] - y0*cracktip[0];
      res[0]= phi/ll;
      res[1]= gmm::vect_sp(xx-cracktip, t);
    } break;

    case 1: {
      res[0] = y;
      res[1] = -0.5 + x;
    } break;
    case 2: {
      res[0] = x - 0.25;
      res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
    } break;
    default: assert(0);
  }
  return res;
}//end ls_function

/*************************************************************************************/
/*     solv  unilateral_contact_problem                                                                              */
/*************************************************************************************/

bool  unilateral_contact_problem::solve(plain_vector &U, plain_vector &LAMBDA) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();
  ls.reinit();  
  
  for (size_type d = 0; d < ls.get_mesh_fem().nb_basic_dof(); ++d) {
    ls.values(0)[d] = ls_function(ls.get_mesh_fem().point_of_basic_dof(d), 0)[0];
    ls.values(1)[d] = ls_function(ls.get_mesh_fem().point_of_basic_dof(d), 0)[1];
  }
  ls.touch();

  mls.adapt();
  mim.adapt();
  mfls_u.adapt();
  mfls_cont.adapt();

  bool load_global_fun =  0;


  cout << "Setting up the singular functions for the enrichment\n";
  cout << dgr << endl;

  std::vector<getfem::pglobal_function> vfunc(dgr*4);
  //std::vector<getfem::pglobal_function> vfunc_u(6);
  //std::vector<getfem::pglobal_function> vfunc_p(dgrp*2);
  if (!load_global_fun) {
  std::cout << "Using default singular functions\n";
    for (unsigned i = 0; i < vfunc.size(); ++i){
      /* use the singularity */
      
      getfem::abstract_xy_function *s = 
	new getfem::crack_singular_xy_function(i);
      
      if (enrichment_option != FIXED_ZONE ) {
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

           
    }
  
  mf_sing_u.set_functions(vfunc);
  
  //vfunc_p.resize(2);
  // mf_sing_p.set_functions(vfunc_p);  

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
     // mf_product_cont.set_enrichment(enriched_dofs);
      mf_cont_sum.set_mesh_fems( mfls_cont);
    }
    break;

  case GLOBAL_WITH_CUTOFF :{
    if(cutoff_func == 0)
      cout<<"Using exponential Cutoff..."<<endl;
    else
      cout<<"Using Polynomial Cutoff..."<<endl;
      mf_u_sum.set_mesh_fems(mf_sing_u, mfls_u);
      mf_cont_sum.set_mesh_fems(mfls_cont);
  } break;

  case NO_ENRICHMENT: {
    mf_u_sum.set_mesh_fems(mfls_u);
    mf_cont_sum.set_mesh_fems(mfls_cont);
  } break;
  
  }// end switch
  
  
  /************************************************************************************/
  /* Range_basis call                                                                 */
  /************************************************************************************/
  std::set<size_type> cols;
  cols.clear();
  sparse_matrix BRBB(mf_u().nb_basic_dof(), mf_cont().nb_basic_dof());
  getfem::asm_mass_matrix(BRBB, mimbound, mf_u(), mf_cont());
  cout << "Selecting dofs for the multiplier" << endl;
  cout << "nb_dof_mult = " << mf_cont().nb_dof() << endl;
  gmm::range_basis(BRBB, cols);
  mf_cont().reduce_to_basic_dof(cols);



  size_type nb_dof = mf_u().nb_basic_dof();
  cout << "nb_dof = " << nb_dof << endl;

  size_type nb_dof_cont = mf_cont().nb_basic_dof();
  cout << "nb_dof_cont = " << nb_dof_cont << endl;

  U.resize(nb_dof);
  LAMBDA.resize(nb_dof_cont);

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
  }//end if zero
  

  // find the dofs on the upper right and lower right corners

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
  }//end for

  GMM_ASSERT1(((d1 < 1E-8) && (d2 < 1E-8)),
	      "Upper right or lower right corners not found d1 = "
	      << d1 << " d2 = " << d2);

  //asembling mixed term for contact problem

 getfem::CONTACT_B_MATRIX BN(nb_dof_cont, nb_dof);
 ls.set_shift(1e-7);
 getfem::asm_mass_matrix(BN, mimbound, mf_cont(), mf_u());
 ls.set_shift(-1e-7);
 getfem::asm_mass_matrix(BN, mimbound, mf_cont(), mf_u());
 ls.set_shift(0.);

 getfem::CONTACT_B_MATRIX CA(nb_dof_cont, nb_dof);
 asm_stabilization_mixed_term(CA, mimbound, mf_u(), mf_cont(), ls, lambda, mu);
 gmm::scale(CA, -cont_gamma0 * h);

 
 // assembling symetrique term for contact problem


sparse_matrix KA(nb_dof, nb_dof);

 asm_stabilization_symm_term(KA, mimbound, mf_u(), ls, lambda, mu);
 gmm::scale(KA, -cont_gamma0 * h);


 //assembling mass term for stabilized problem
 
 sparse_matrix MA(nb_dof_cont, nb_dof_cont);

 ls.set_shift(1e-7);
 getfem::asm_mass_matrix(MA, mimbound, mf_cont());
 ls.set_shift(-1e-7);
 getfem::asm_mass_matrix(MA, mimbound, mf_cont());
 ls.set_shift(0.);
 gmm::scale(MA, -cont_gamma0 * h);



  getfem::model model;

  // Main unknown of the problem.
  model.add_fem_variable("u", mf_u());
  model.add_fem_variable("Lambda", mf_cont()); // Adding contact variable


 // Linearized elasticity brick.
  model.add_initialized_scalar_data("lambda", lambda);
  model.add_initialized_scalar_data("mu", mu);
  getfem::add_isotropic_linearized_elasticity_brick
    (model, mim, "u", "lambda", "mu");
    
  getfem::add_explicit_matrix(model, "u", "u", KA);

  model.add_initialized_scalar_data("augmentation_parameter", R);

 if (stabilized_problem) {
    gmm::add(CA, BN);
  }

 // Defining the contact condition.

 getfem::add_basic_contact_brick(model, "u", "Lambda", "augmentation_parameter", BN,"","", false);

  
  // Defining the Neumann condition right hand side.

  std::vector<scalar_type> F(nb_dof_rhs * N);

  // Neumann condition brick.
  
  //down side
  for(size_type i = 1; i < F.size(); i=i+N) F[i] = 1.;
  for(size_type i = 0; i < F.size(); i=i+N) F[i] = 1.;
  
 
  model.add_initialized_fem_data("NeumannData", mf_rhs, F);
  getfem::add_normal_source_term_brick
	(model, mim, "u", "NeumannData", NEUMANN1_BOUNDARY_NUM);
  getfem::add_normal_source_term_brick
	(model, mim, "u", "NeumannData", NEUMANN2_BOUNDARY_NUM);

  gmm::scale(F, -1.0);
  model.add_initialized_fem_data("NeumannData2", mf_rhs, F);
  getfem::add_normal_source_term_brick
        (model, mim, "u", "NeumannData2", NEUMANN3_BOUNDARY_NUM);
  getfem::add_normal_source_term_brick
	(model, mim, "u", "NeumannData2", NEUMANN4_BOUNDARY_NUM);




GMM_ASSERT1(N==2, "To be corrected for 3D computation");
sparse_matrix BB(3, mf_u().nb_dof());
BB(0, icorner1) = 1.0;
BB(1, icorner1+1) = 1.0;
BB(2, icorner2) = 1.0;

std::vector<scalar_type> LRH(3);
getfem::add_constraint_with_multipliers(model, "u", "mf_mult_dir", BB, LRH);


// Generic solve.

    gmm::iteration iter(residual, 1, 40000);
    cout << "Solving..." << endl;
    iter.init();
    getfem::standard_solve(model, iter);
    gmm::resize(U, mf_u().nb_dof());
    gmm::copy(model.real_variable("u"), U);
    gmm::resize(LAMBDA, mf_cont().nb_dof());
    gmm::copy(model.real_variable("Lambda"), LAMBDA);
    
 //  Exporting reference solution  
    
if (reference_test) {
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
      getfem::interpolation(mf_cont(), mf_refined, LAMBDA, PP);
      mf_refined.write_to_file(datafilename + "_refined_test.cont_meshfem_refined", true);
      gmm::vecsave(datafilename + "_refined_test.cont_refined", PP);
      
      cout << "done" << endl;
    }//end reference_test
    
  return (iter.converged());


}//end solv






/* ************************************************************************************/
/* Main program                                                                       */
/**************************************************************************************/

int main(int argc, char *argv[]) {


 FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  //getfem::getfem_mesh_level_set_noisy();


  unilateral_contact_problem p;
  p.PARAM.read_command_line(argc, argv);
  p.init();
  p.mesh.write_to_file(p.datafilename + ".mesh");
  
  plain_vector U, P;
  if (!p.solve(U, P)) GMM_ASSERT1(false,"Solve has failed");

  cout << "Saving the solution" << endl;
  getfem::mesh mcut;
  p.mls.global_cut_mesh(mcut);
  unsigned Q = p.mf_u().get_qdim();

  cout << p.mf_u().nb_basic_dof() << endl;
  cout << p.mf_cont().nb_basic_dof() << endl;

 
  
}


