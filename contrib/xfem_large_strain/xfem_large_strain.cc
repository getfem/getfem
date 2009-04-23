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
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_mesh_fem_product.h"
#include "getfem/getfem_mesh_fem_global_function.h"
#include "getfem/getfem_spider_fem.h"
#include "getfem/getfem_mesh_fem_sum.h"
#include "gmm/gmm.h"
#include "gmm/gmm_inoutput.h"

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
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
/*  Structure for the crack problem.                                      */
/**************************************************************************/

struct crack_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN1_BOUNDARY_NUM = 1, NEUMANN2_BOUNDARY_NUM=2, NEUMANN3_BOUNDARY_NUM=3, NEUMANN4_BOUNDARY_NUM=4, MORTAR_BOUNDARY_IN=42, MORTAR_BOUNDARY_OUT=43};
  getfem::mesh mesh;  /* the mesh */
  getfem::level_set ls;      /* The two level sets defining the crack.       */
  getfem::mesh_level_set mls;       /* the integration methods.              */
  getfem::mesh_im_level_set mim;    /* the integration methods.              */
  getfem::mesh_fem mf_pre_u, mf_pre_mortar;
  getfem::mesh_fem mf_mult;
  getfem::mesh_fem_level_set mfls_u, mfls_mortar;
  getfem::mesh_fem_global_function mf_sing_u;
  getfem::mesh_fem mf_partition_of_unity;
  getfem::mesh_fem_product mf_product;
  getfem::mesh_fem_sum mf_u_sum;
  
  getfem::mesh_fem mf_pre_p; /* mesh_fem for the pressure for mixed form     */
  getfem::mesh_fem_level_set mfls_p;   /* mesh_fem for the pressure enriched with H.   */
 


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

  getfem::mesh_fem& mf_u() { return mf_u_sum; }
  // getfem::mesh_fem& mf_u() { return mf_us; }
  
  scalar_type mu;            /* Lame coefficients.                */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  
  
  int bimaterial;           /* For bimaterial interface fracture */
  double mu_up, mu_down;  /*Lame coeff for bimaterial case*/
 
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
  std::string GLOBAL_FUNCTION_MF, GLOBAL_FUNCTION_U;

  bgeot::md_param PARAM;

  bool solve(plain_vector &U, plain_vector &P);
  void init(void);
  crack_problem(void) : ls(mesh, 1, true), mls(mesh), mim(mls), 
			mf_pre_u(mesh), mf_pre_mortar(mesh), mf_mult(mesh),
			mfls_u(mls, mf_pre_u), mfls_mortar(mls, mf_pre_mortar),	
			mf_sing_u(mesh),
			mf_partition_of_unity(mesh),
			mf_product(mf_partition_of_unity, mf_sing_u),

			mf_u_sum(mesh), mf_pre_p(mesh), mfls_p(mls, mf_pre_p),
			mf_us(mesh), mf_rhs(mesh)
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
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
					 "Name of simplex integration method");
  std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION");

  enrichment_option = enrichment_option_enum(PARAM.int_value("ENRICHMENT_OPTION",
							     "Enrichment option"));
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  bimaterial = int(PARAM.int_value("BIMATERIAL", "bimaterial interface crack"));
  
  mu = PARAM.real_value("MU", "Lame coefficient mu"); 
  if (bimaterial == 1){
    mu_up = PARAM.real_value("MU_UP", "Lame coefficient mu"); 
    mu_down = PARAM.real_value("MU_DOWN", "Lame coefficient mu"); 
  }
  
  
  spider.radius = PARAM.real_value("SPIDER_RADIUS","spider_radius");
  spider.Nr = unsigned(PARAM.int_value("SPIDER_NR","Spider_Nr "));
  spider.Ntheta = unsigned(PARAM.int_value("SPIDER_NTHETA","Ntheta "));
  spider.K = int(PARAM.int_value("SPIDER_K","K "));
  spider.theta0 = 0;
  spider.bimat_enrichment = int(PARAM.int_value("SPIDER_BIMAT_ENRICHMENT","spider_bimat_enrichment"));

  /* The following constants are taken from the Chang and Xu paper of 2007 
     titled The singular stress field and stress intensity factors of a crack 
     terminating at a bimaterial interface published in IJMECSCI. epsilon 
     being the constant present in the asymptotic displacement analytic solution: 
     r^{1/2} * cos (\epsilon Ln(r)) * f(\theta) and  r^{1/2} * sin (\epsilon Ln(r)) * f(\theta).
  */
  
  if (spider.bimat_enrichment == 1){
    scalar_type nu1 = 10;
    scalar_type nu2 = 3;
    scalar_type kappa1 = 3. - 4. * nu1;
    scalar_type kappa2 = 3. - 4. * nu2;
        
    scalar_type beta = (mu_up*(kappa2-1.) - mu_down*(kappa1-1.)) / (mu_up*(kappa2+1.) - mu_down*(kappa1+1.));
      spider.epsilon = 1./(2.*M_PI) * log( (1.-beta) / (1.+beta) );
  
    
    //spider.epsilon = PARAM.real_value("SPIDER_BIMAT_ENRICHMENT","spider_bimat_enrichment");
  }
  

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
  
  cracktip.resize(2); // Coordonn�e du fond de fissure
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


  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL");
  if (residual == 0.) residual = 1e-10;
  enr_area_radius = PARAM.real_value("RADIUS_ENR_AREA",
				     "radius of the enrichment area");
  
  GLOBAL_FUNCTION_MF = PARAM.string_value("GLOBAL_FUNCTION_MF");
  GLOBAL_FUNCTION_U = PARAM.string_value("GLOBAL_FUNCTION_U");

  reference_test = int(PARAM.int_value("REFERENCE_TEST", "Reference test")); 


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
  
//   if (enrichment_option == 3 || enrichment_option == 4) {
//     spider = new getfem::spider_fem(spider_radius, mim, spider_Nr,
// 				    spider_Ntheta, spider_K, cracktip,
// 				    theta0);
//     mf_us.set_finite_element(mesh.convex_index(),spider->get_pfem());
//     for (dal::bv_visitor_c i(mf_us.convex_index()); !i.finished(); ++i) {
//       if (mf_us.fem_of_element(i)->nb_dof(i) == 0) {
// 	mf_us.set_finite_element(i,0);
//       }
//     }
//   }

  mixed_pressure =
    (PARAM.int_value("MIXED_PRESSURE","Mixed version or not.") != 0);
  dir_with_mult = unsigned(PARAM.int_value("DIRICHLET_VERSION",
					   "Version of Dirichlet"));
  if (mixed_pressure) {
    std::string FEM_TYPE_P  = PARAM.string_value("FEM_TYPE_P","FEM name P");
    mf_pre_p.set_finite_element(mesh.convex_index(),
			    getfem::fem_descriptor(FEM_TYPE_P));
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


bool crack_problem::solve(plain_vector &U, plain_vector &P) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();
  ls.reinit();  
  
  for (size_type d = 0; d < ls.get_mesh_fem().nb_dof(); ++d) {
    ls.values(0)[d] = ls_function(ls.get_mesh_fem().point_of_dof(d), 0)[0];
    ls.values(1)[d] = ls_function(ls.get_mesh_fem().point_of_dof(d), 0)[1];
  }
  ls.touch();

  mls.adapt();
  mim.adapt();
  mfls_u.adapt();
  mfls_p.adapt();
  mfls_mortar.adapt(); mfls_mortar.set_qdim(2);

  bool load_global_fun = GLOBAL_FUNCTION_MF.size() != 0;


  cout << "Setting up the singular functions for the enrichment\n";

  std::vector<getfem::pglobal_function> vfunc(4);
  if (!load_global_fun) {
    cout << "Using default singular functions\n";
    for (unsigned i = 0; i < vfunc.size(); ++i){
      /* use the singularity */
      getfem::abstract_xy_function *s = 
	new getfem::crack_singular_xy_function(i);
      if (enrichment_option != FIXED_ZONE && 
	  enrichment_option != GLOBAL_WITH_MORTAR) {
	/* use the product of the singularity function
	   with a cutoff */
	getfem::abstract_xy_function *c = 
	  new getfem::cutoff_xy_function(int(cutoff_func),
					 cutoff_radius, 
					 cutoff_radius1,cutoff_radius0);
	s = new getfem::product_of_xy_functions(*s, *c);
      }
      vfunc[i] = getfem::global_function_on_level_set(ls, *s);
    }
  } else {
    cout << "Load singular functions from " << GLOBAL_FUNCTION_MF << " and " << GLOBAL_FUNCTION_U << "\n";
    getfem::mesh *m = new getfem::mesh(); 
    m->read_from_file(GLOBAL_FUNCTION_MF);
    getfem::mesh_fem *mf_c = new getfem::mesh_fem(*m); 
    mf_c->read_from_file(GLOBAL_FUNCTION_MF);
    std::fstream f(GLOBAL_FUNCTION_U.c_str(), std::ios::in);
    plain_vector W(mf_c->nb_dof());


  
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
    } break;

  case GLOBAL_WITH_CUTOFF :{
    if(cutoff_func == 0)
      cout<<"Using exponential Cutoff..."<<endl;
    else
      cout<<"Using Polynomial Cutoff..."<<endl;
    mf_u_sum.set_mesh_fems(mf_sing_u, mfls_u);
  } break;

  case SPIDER_FEM_ALONE : { 
    mf_u_sum.set_mesh_fems(mf_us); 
  } break;
    
  case SPIDER_FEM_ENRICHMENT : {
    mf_u_sum.set_mesh_fems(mf_us, mfls_u); 
  } break;
    
  case NO_ENRICHMENT: {
    mf_u_sum.set_mesh_fems(mfls_u);
  } break;
  
  }
  

  U.resize(mf_u().nb_dof());
  P.resize(mfls_p.nb_dof());


  if (mixed_pressure)
    cout << "Number of dof for P: " << mfls_p.nb_dof() << endl;
  cout << "Number of dof for u: " << mf_u().nb_dof() << endl;

  unsigned Q = mf_u().get_qdim();
  if (0) {
    for (unsigned d=0; d < mf_u().nb_dof(); d += Q) {
      printf("dof %4d @ %+6.2f:%+6.2f: ", d, 
	     mf_u().point_of_dof(d)[0], mf_u().point_of_dof(d)[1]);

      const getfem::mesh::ind_cv_ct cvs = mf_u().convex_to_dof(d);
      for (unsigned i=0; i < cvs.size(); ++i) {
	size_type cv = cvs[i];
	//if (pm_cvlist.is_in(cv)) flag1 = true; else flag2 = true;

	getfem::pfem pf = mf_u().fem_of_element(cv);
	unsigned ld = unsigned(-1);
	for (unsigned dd = 0; dd < mf_u().nb_dof_of_element(cv); dd += Q) {
	  if (mf_u().ind_dof_of_element(cv)[dd] == d) {
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
  

  // find the dofs on the upper right and lower right corners

  scalar_type d1 = 1.0, d2 = 1.0;
  size_type icorner1 = size_type(-1), icorner2 = size_type(-1);
  base_node corner1 = base_node(1.0, -0.5);
  base_node corner2 = base_node(1.0, 0.5);
  for (size_type i = 0; i < mf_u().nb_dof(); i+=N) {
    scalar_type dd1 = gmm::vect_dist2(mf_u().point_of_dof(i), corner1);
    if (dd1 < d1) { icorner1 = i; d1 = dd1; }
    scalar_type dd2 = gmm::vect_dist2(mf_u().point_of_dof(i), corner2);
    if (dd2 < d2) { icorner2 = i; d2 = dd2; }
  }

  GMM_ASSERT1(((d1 < 1E-8) && (d2 < 1E-8)),
	      "Upper right or lower right corners not found d1 = "
	      << d1 << " d2 = " << d2);

  // Linearized elasticity brick.
  
  
  getfem::mdbrick_isotropic_linearized_elasticity<>
    ELAS(mim, mf_u(), mixed_pressure ? 0.0 : 0.0, mu);

  
  if(bimaterial == 1){
    cout<<"______________________________________________________________________________"<<endl;
    cout<<"CASE OF BIMATERIAL CRACK  with mu_up = "<<mu_up<<" and mu_down = "<<mu_down<<endl;
    cout<<"______________________________________________________________________________"<<endl;
    std::vector<double> bi_lambda(ELAS.lambda().mf().nb_dof());
    std::vector<double> bi_mu(ELAS.lambda().mf().nb_dof());
    
    cout<<"ELAS.lambda().mf().nb_dof()==="<<ELAS.lambda().mf().nb_dof()<<endl;
    
    for (size_type ite = 0; ite < ELAS.lambda().mf().nb_dof();ite++) {
      if (ELAS.lambda().mf().point_of_dof(ite)[1] > 0){
	bi_lambda[ite] = 0;
	bi_mu[ite] = mu_up;
      }
	else{
	  bi_lambda[ite] = 0;
	  bi_mu[ite] = mu_down;
	}
    }
    
    //cout<<"bi_lambda.size() = "<<bi_lambda.size()<<endl;
    // cout<<"ELAS.lambda().mf().nb_dof()==="<<ELAS.lambda().mf().nb_dof()<<endl;
    
    ELAS.lambda().set(bi_lambda);
    ELAS.mu().set(bi_mu);
  }
  

  getfem::mdbrick_abstract<> *pINCOMP;
  getfem::mdbrick_linear_incomp<> *incomp;
  if (mixed_pressure) {
    incomp = new getfem::mdbrick_linear_incomp<>(ELAS, mfls_p);
    // incomp->penalization_coeff().set(1.0/lambda);
    pINCOMP = incomp;
  } else { pINCOMP = &ELAS; incomp = 0; }

  // Defining the Neumann condition right hand side.
  plain_vector F(nb_dof_rhs * N);

  // Neumann condition brick.
  
  //down side
  for(size_type i = 1; i < F.size(); i=i+N) F[i] = 1;
  for(size_type i = 0; i < F.size(); i=i+2) F[i] = 1;
  
  getfem::mdbrick_source_term<> NEUMANN1(*pINCOMP, mf_rhs, F,
					 NEUMANN1_BOUNDARY_NUM);
  getfem::mdbrick_source_term<> NEUMANN2(NEUMANN1, mf_rhs, F,
					 NEUMANN2_BOUNDARY_NUM);
  gmm::scale(F, -1.0);
  getfem::mdbrick_source_term<> NEUMANN3(NEUMANN2, mf_rhs, F,
					 NEUMANN3_BOUNDARY_NUM);
  getfem::mdbrick_source_term<> NEUMANN4(NEUMANN3, mf_rhs, F,
					 NEUMANN4_BOUNDARY_NUM);
  

  getfem::mdbrick_constraint<> KILL_RIGID_MOTIONS(NEUMANN4);
  GMM_ASSERT1(N==2, "To be corrected for 3D computation");
  sparse_matrix BB(3, mf_u().nb_dof());
  BB(0, icorner1) = 1.0;
  BB(1, icorner1+1) = 1.0;
  BB(2, icorner2) = 1.0;
  KILL_RIGID_MOTIONS.set_constraints(BB, plain_vector(3));
  KILL_RIGID_MOTIONS.set_constraints_type(getfem::constraints_type(dir_with_mult));

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<> DIRICHLET(KILL_RIGID_MOTIONS, DIRICHLET_BOUNDARY_NUM, mf_mult);
    
  DIRICHLET.set_constraints_type(getfem::constraints_type(dir_with_mult));

  getfem::mdbrick_abstract<> *final_model = &DIRICHLET;

 if (enrichment_option == GLOBAL_WITH_MORTAR) {
    /* add a constraint brick for the mortar junction between
       the enriched area and the rest of the mesh */
    /* we use mfls_u as the space of lagrange multipliers */
    getfem::mesh_fem &mf_mortar = mfls_mortar; 
    /* adjust its qdim.. this is just evil and dangerous
       since mf_u() is built upon mfls_u.. it would be better
       to use a copy. */
    mf_mortar.set_qdim(2); // EVIL 
    getfem::mdbrick_constraint<> &mortar = 
      *(new getfem::mdbrick_constraint<>(DIRICHLET,0));

    cout << "Handling mortar junction\n";

    /* list of dof of mf_mortar for the mortar condition */
    std::vector<size_type> ind_mortar;
    /* unfortunately , dof_on_region sometimes returns too much dof
       when mf_mortar is an enriched one so we have to filter them */
    sparse_matrix M(mf_mortar.nb_dof(), mf_mortar.nb_dof());
    getfem::asm_mass_matrix(M, mim, mf_mortar, MORTAR_BOUNDARY_OUT);
    for (dal::bv_visitor_c d(mf_mortar.dof_on_region(MORTAR_BOUNDARY_OUT)); 
	 !d.finished(); ++d) {
      if (M(d,d) > 1e-8) ind_mortar.push_back(d);
      else cout << "  removing non mortar dof" << d << "\n";
    }

    cout << ind_mortar.size() << " dof for the lagrange multiplier)\n";

    sparse_matrix H0(mf_mortar.nb_dof(), mf_u().nb_dof()), 
      H(ind_mortar.size(), mf_u().nb_dof());


    gmm::sub_index sub_i(ind_mortar);
    gmm::sub_interval sub_j(0, mf_u().nb_dof());
    
    /* build the mortar constraint matrix -- note that the integration
       method is conformal to the crack
     */
    getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u(), 
			    MORTAR_BOUNDARY_OUT);
    gmm::copy(gmm::sub_matrix(H0, sub_i, sub_j), H);

    gmm::clear(H0);
    getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u(), 
			    MORTAR_BOUNDARY_IN);
    gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i, sub_j), -1), H);


    /* because of the discontinuous partition of mf_u(), some levelset
       enriched functions do not contribute any more to the
       mass-matrix (the ones which are null on one side of the
       levelset, when split in two by the mortar partition, may create
       a "null" dof whose base function is all zero..
    */
    sparse_matrix M2(mf_u().nb_dof(), mf_u().nb_dof());
    getfem::asm_mass_matrix(M2, mim, mf_u(), mf_u());

    for (size_type d = 0; d < mf_u().nb_dof(); ++d) {
      if (M2(d,d) < 1e-10) {
	cout << "  removing null mf_u() dof " << d << "\n";
	size_type n = gmm::mat_nrows(H);
	gmm::resize(H, n+1, gmm::mat_ncols(H));
	H(n, d) = 1;
      }
    }

    
    
    getfem::base_vector R(gmm::mat_nrows(H));
    mortar.set_constraints(H,R);

    final_model = &mortar;
  }



  // Computation of the inf-sup bound 
  if (PARAM.int_value("INF_SUP_COMP") && mixed_pressure) {
    
    cout << "Sparse matrices computation for the test of inf-sup condition"
	 << endl;

    sparse_matrix Mis(mfls_p.nb_dof(), mfls_p.nb_dof());
    sparse_matrix Sis(mf_u().nb_dof(), mf_u().nb_dof());

    getfem::asm_mass_matrix(Mis, mim, mfls_p);
    getfem::asm_stiffness_matrix_for_homogeneous_laplacian_componentwise
      (Sis, mim, mf_u());
    getfem::asm_mass_matrix(Sis, mim, mf_u());
    
    cout << "Inf-sup condition test" << endl;
    scalar_type lambda = smallest_eigen_value(incomp->get_B(), Mis, Sis);
    cout << "The inf-sup test gives " << lambda << endl;
  }

  // Generic solve.
  cout << "Total number of variables : " << final_model->nb_dof() << endl;
  getfem::standard_model_state MS(*final_model);
  gmm::iteration iter(residual, 2, 40000);
  cout << "Solving..." << endl;
  getfem::standard_solve(MS, *final_model, iter,
			 getfem::select_linear_solver(*final_model,"superlu"));
  cout << "Solving... done" << endl;
  // Solution extraction
  gmm::copy(ELAS.get_solution(MS), U);
  if (mixed_pressure)  gmm::copy(incomp->get_pressure(MS), P);

  
  if (reference_test) {
      cout << "Exporting reference solution...";
      dal::bit_vector blocked_dof = mf_u().dof_on_set(5);
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
      getfem::interpolation(mfls_p, mf_refined, P, PP);
      mf_refined.write_to_file(datafilename + "_refined_test.p_meshfem_refined", true);
      gmm::vecsave(datafilename + "_refined_test.P_refined", PP);
     
      
      cout << "done" << endl;
    }
  
  
  return (iter.converged());
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
  
  getfem::interpolation(p.mfls_p, mf_p, P, PP);
  
  mf_p.write_to_file(p.datafilename + ".p_meshfem", true);
  gmm::vecsave(p.datafilename + ".P", PP);
      
  cout << "Interpolating solution for the drawing" << endl;
  getfem::stored_mesh_slice sl;
  getfem::mesh mcut_refined;
  
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
  
  getfem::mesh_fem mf_refined(mcut_refined, dim_type(Q));
  mf_refined.set_classical_discontinuous_finite_element(2, 0.001);
  plain_vector W(mf_refined.nb_dof());
  
  getfem::interpolation(p.mf_u(), mf_refined, U, W);
  
  
  if (p.PARAM.int_value("VTK_EXPORT")) {
    getfem::mesh_fem mf_refined_vm(mcut_refined, 1);
    mf_refined_vm.set_classical_discontinuous_finite_element(1, 0.001);
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
    std::string REFERENCE_MF = "large_strain_refined_test.meshfem_refined";
    std::string REFERENCE_U = "large_strain_refined_test.U_refined";
    std::string REFERENCE_MFP = "large_strain_refined_test.p_meshfem_refined";
    std::string REFERENCE_P = "large_strain_refined_test.P_refined";
    
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
    getfem::interpolation(p.mfls_p, ref_mfp, P, interp_P);
    
    cout << "To ref L2 ERROR on P:"
	 << getfem::asm_L2_dist(ref_mim, ref_mfp, interp_P,
				ref_mfp, ref_P) << endl;
    
    gmm::add(gmm::scaled(interp_U, -1.), ref_U);
    gmm::vecsave(p.datafilename + ".diff_ref", ref_U);
    
  }
  
  return 0; 
}
