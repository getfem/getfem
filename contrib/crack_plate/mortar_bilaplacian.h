// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2008 Yves Renard, Julien Pommier.
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

#include "getfem/getfem_config.h"
#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_fourth_order.h"
#include "getfem/getfem_model_solvers.h"
#include "gmm/gmm.h"
#include "getfem/getfem_superlu.h"
#include "getfem/getfem_derivatives.h"
#include "gmm/gmm_inoutput.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_mesh_fem_product.h"
#include "getfem/getfem_mesh_fem_global_function.h"
#include "getfem/getfem_mesh_fem_sum.h"

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::base_vector; /* dense vector. */
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* dense matrix. */

/* definition of some matrix/vector types. 
 * default types of getfem_model_solvers.h
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/**************************************************************************/
/*  Exact solution.                                                       */
/**************************************************************************/

scalar_type FT = 10.0;
scalar_type D  = 1.0;
scalar_type nu = 0.0;
scalar_type pressure  = 1.0 ;

#if 1

scalar_type sol_u(const base_node &x)
{ return sin(FT*std::accumulate(x.begin(), x.end(), 0.0)); }
scalar_type sol_lapl_u(const base_node &x)
{ return -FT*FT*sol_u(x) * x.size(); }
scalar_type sol_f(const base_node &x)
{ return FT*FT*FT*FT*sol_u(x)*gmm::sqr(x.size()); }
base_small_vector sol_du(const base_node &x) {
  base_small_vector res(x.size());
  std::fill(res.begin(), res.end(),
	    FT * cos(FT*std::accumulate(x.begin(), x.end(), 0.0)));
  return res;
}
base_small_vector neumann_val(const base_node &x)
{ return -FT*FT*sol_du(x) * scalar_type(x.size()); }

base_matrix sol_hessian(const base_node &x) {
  base_matrix m(x.size(), x.size());
  std::fill(m.begin(), m.end(), -FT*FT*sol_u(x));
  return m;
}

#else

scalar_type sol_u(const base_node &x) {
 return pressure/(2. * D) * x[0] * x[0] * ( x[0] * x[0] / 12. - x[0] / 3. + 1./2. ); }
 
scalar_type sol_lapl_u(const base_node &) { 
return pressure/(2. * D) * ( x[0] * x[0] - 2. * x[0] + 1. ); }

scalar_type sol_f(const base_node &) { return pressure ; }
base_small_vector sol_du(const base_node &x) {
  base_small_vector res(x.size()); 
  res[0] = pressure/(2. * D) * x[0] * ( x[0] * x[0] / 3. - x[0] + 1. ); 
  res[1] = 0.; 
  return res;
}
base_small_vector neumann_val(const base_node &x)
{ cout << "erreur : KL = 1 doit etre selectionne \n";
  return 0 ; /*base_small_vector(x.size()); */ }

base_matrix sol_hessian(const base_node &x)
{ base_matrix m(x.size(), x.size()); 
  m(0,0) = sol_lapl_u(x);
  return m; }


// scalar_type sol_u(const base_node &x) { return x[0]*x[1]; }
// scalar_type sol_lapl_u(const base_node &) { return 0.0; }
// scalar_type sol_f(const base_node &) { return 0.0; }
// base_small_vector sol_du(const base_node &x) {
//   base_small_vector res(x.size()); res[0] = x[1]; res[1] = x[0]; 
//   return res;
// }
// base_small_vector neumann_val(const base_node &x)
// { return base_small_vector(x.size());  }
// 
// base_matrix sol_hessian(const base_node &x)
// { base_matrix m(x.size(), x.size()); m(1,0) = m(0,1) = 1.0; return m; }

#endif

base_matrix sol_mtensor(const base_node &x) { 
  base_matrix m = sol_hessian(x), mm(x.size(), x.size());
  scalar_type l = sol_lapl_u(x);
  for (size_type i = 0; i < x.size(); ++i) mm(i,i) = l * nu;
  gmm::scale(m, (1-nu));
  gmm::add(mm, m);
  gmm::scale(m, -D);
  return m;
}

base_small_vector sol_bf(const base_node &x)
{ return -D * neumann_val(x); }




namespace getfem {
// function for assembling the constraints of the integral matching 
  template<typename MAT, typename VECT1, typename VECT2>
  void asm_normal_derivative_dirichlet_constraints_bis 
 (MAT &H, VECT1 &R, const mesh_im &mim, const mesh_fem &mf_u,
   const mesh_fem &mf_mult, const mesh_fem &mf_r,
   const VECT2 &r_data, const mesh_region &rg, bool R_must_be_derivated, 
   int version) {
    typedef typename gmm::linalg_traits<VECT1>::value_type value_type;
    typedef typename gmm::number_traits<value_type>::magnitude_type magn_type;
    
    rg.from_mesh(mim.linked_mesh()).error_if_not_faces();
    GMM_ASSERT1(mf_r.get_qdim() == 1, 
		"invalid data mesh fem (Qdim=1 required)");
    if (version & ASMDIR_BUILDH) {
      const char *s;
      if (mf_u.get_qdim() == 1 && mf_mult.get_qdim() == 1)
	s = "M(#1,#2)+=comp(Grad(#1).Normal().Grad(#2).Normal())(:,i,i,:,j,j)";
      else
	s = "M(#1,#2)+=comp(vGrad(#1).Normal().vGrad(#2).Normal())(:,j,i,i,:,j,k,k);";
      
      getfem::generic_assembly assem(s);
      assem.push_mi(mim);
      assem.push_mf(mf_mult);
      assem.push_mf(mf_u);
      assem.push_mat(H);
      assem.assembly(rg);
      gmm::clean(H, gmm::default_tol(magn_type())
		 * gmm::mat_maxnorm(H) * magn_type(1000));
    }
    if (version & ASMDIR_BUILDR) {
      if (!R_must_be_derivated) {
	asm_normal_source_term(R, mim, mf_mult, mf_r, r_data, rg);
      } else {
	asm_real_or_complex_1_param
	  (R, mim, mf_mult, mf_r, r_data, rg,
	   "R=data(#2); V(#1)+=comp(Grad(#1).Normal().Grad(#2).Normal())(i,j,k,k).R(j)");
      }
    }
  }
  
  template<typename MAT>
void asm_constraint_gradient_vectorial_mult 
(MAT &H, const getfem::mesh_im &mim, const getfem::mesh_fem &mf_u,
 const getfem::mesh_fem &mf_mult, 
 const getfem::mesh_region &rg, int version) {
  typedef typename gmm::linalg_traits<MAT>::value_type value_type;
  typedef typename gmm::number_traits<value_type>::magnitude_type magn_type;
    
  rg.from_mesh(mim.linked_mesh()).error_if_not_faces();
  if (version & getfem::ASMDIR_BUILDH) {
    assert(mf_u.get_qdim() == 1 && mf_mult.get_qdim() == 2);
    const char *s = "M(#1,#2)+=comp(vBase(#1).Grad(#2))(:,i,:,i)";
    getfem::generic_assembly assem(s);
    assem.push_mi(mim);
    assem.push_mf(mf_mult);
    assem.push_mf(mf_u);
    assem.push_mat(H);
    assem.assembly(rg);
    gmm::clean(H, gmm::default_tol(magn_type())
               * gmm::mat_maxnorm(H) * magn_type(1000));
  }
}

} // end getfem

/******** Struct for the Bilaplacian problem *******************************/


struct bilaplacian_mortar_problem {

  enum { CLAMPED_BOUNDARY_NUM = 0, SIMPLE_SUPPORT_BOUNDARY_NUM = 1,
	 FORCE_BOUNDARY_NUM = 2, MOMENTUM_BOUNDARY_NUM = 3,
	 MORTAR_BOUNDARY_IN = 40, MORTAR_BOUNDARY_OUT = 41};
  
  getfem::mesh mesh;        /* the mesh */
  getfem::mesh_im mim;    /* the integration methods.              */
  getfem::mesh_fem mf_u, mf_mortar, mf_mortar_deriv ;  
  
  getfem::mesh_fem mf_rhs;  /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_mult; /* mesh_fem for the Dirichlet condition.        */
  getfem::mesh_fem mf_mult_d; /* mesh_fem for the normal_derivative Dirichlet condition. */
  
  
  scalar_type residual;     /* max residual for the iterative solvers       */
  getfem::constraints_type dirichlet_version, mortar_version;
  
  scalar_type radius_mortar;
  size_type NX;
  
  std::string datafilename;
  bgeot::md_param PARAM;

  bool KL;

  dal::bit_vector pm_convexes; /* convexes inside the enrichment 
				  area when point-wise matching is used.*/
  
  scalar_type epsilon ;      /* half-plate thickness */
  

    
  bool solve(plain_vector &U);
  void init(void);
  void compute_error(plain_vector &U);
  void compute_error_without_assembling(plain_vector &U, plain_vector &V);
  
  bilaplacian_mortar_problem(void) : mim(mesh), mf_u(mesh), 
				    mf_mortar(mesh), mf_mortar_deriv(mesh), 
				    mf_rhs(mesh), mf_mult(mesh), mf_mult_d(mesh)  
				    { KL = true; } 
};


/*                                                          */
/*****  Methods for class bilaplacian_mortar_problem  ********/
/*                                                          */





void bilaplacian_mortar_problem::init(void) {
  std::string MESH_FILE = PARAM.string_value("MESH_FILE");
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  size_type N = 2 ;

    cout << "MESH_TYPE=" << MESH_TYPE << "\n";
    

    radius_mortar = PARAM.real_value("RADIUS_MORTAR",
				     "radius of the enrichment area");
    
    /* First step : build the mesh */
    if (!MESH_FILE.empty()) {
    mesh.read_from_file(MESH_FILE);
    base_small_vector tt(N); 
    tt[0] = PARAM.real_value("TRANSLAT_X") ; //0.02 ; 
    tt[1] = PARAM.real_value("TRANSLAT_Y") ;// 0.04 ; 
    cout << "TRANSLAT_X = " << tt[0] << " ; TRANSLAT_Y = " << tt[1] << "\n" ;
    mesh.translation(tt); 
    MESH_TYPE = bgeot::name_of_geometric_trans
      (mesh.trans_of_convex(mesh.convex_index().first_true()));
    bgeot::pgeometric_trans pgt = 
      bgeot::geometric_trans_descriptor(MESH_TYPE);
    cout << "MESH_TYPE=" << MESH_TYPE << "\n";
    N = mesh.dim();
    } else {
    bgeot::pgeometric_trans pgt = 
      bgeot::geometric_trans_descriptor(MESH_TYPE);
    N = pgt->dim();
    GMM_ASSERT1(N == 2, "For a plate problem, N should be 2");
    std::vector<size_type> nsubdiv(N);
    NX = PARAM.int_value("NX", "Number of space steps ") ;
    std::fill(nsubdiv.begin(),nsubdiv.end(), NX);
    getfem::regular_unit_mesh(mesh, nsubdiv, pgt, PARAM.int_value("MESH_NOISED") != 0);

    bgeot::base_matrix M(N,N);
    for (size_type i=0; i < N; ++i) {
      static const char *t[] = {"LX","LY","LZ"};
      M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
    }
    /* scale the unit mesh to [LX,LY,..] and incline it */
    mesh.transformation(M);
    
    base_small_vector tt(N); 
    tt[0] = -0.5 ;
    tt[1] = -0.5;
    mesh.translation(tt); 
  }
    
    
   /* read the parameters   */
  epsilon = PARAM.real_value("EPSILON", "thickness") ;
  int dv = PARAM.int_value("DIRICHLET_VERSION", "Dirichlet version");
  int mv = PARAM.int_value("MORTAR_VERSION", "Mortar version");
  dirichlet_version = getfem::constraints_type(dv);
  mortar_version = getfem::constraints_type(mv);
  datafilename=PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual=PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;


  // Setting the integration methods 
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);    
  mim.set_integration_method(mesh.convex_index(), ppi);
  
  // Setting the finite element on the mf_u  
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  mf_u.set_finite_element(mesh.convex_index(), pf_u); 
  mf_mortar.set_finite_element(mesh.convex_index(),
             getfem::fem_descriptor(PARAM.string_value("MORTAR_FEM_TYPE")));
  mf_mortar_deriv.set_finite_element(mesh.convex_index(),
             getfem::fem_descriptor(PARAM.string_value("MORTAR_DERIV_FEM_TYPE")));
  
  // set the mesh_fem of the multipliers (for the dirichlet condition)    
  std::string dirichlet_fem_name = PARAM.string_value("DIRICHLET_FEM_TYPE");
  std::string dirichlet_derivative_fem_name = PARAM.string_value("DIRICHLET_DER_FEM_TYPE");
  
  if (dirichlet_fem_name.size() == 0){
    mf_mult.set_finite_element(mesh.convex_index(), pf_u);
    mf_mult_d.set_finite_element(mesh.convex_index(), pf_u);}
  else {
    cout << "DIRICHLET_FEM_TYPE="  << dirichlet_fem_name << "\n";
    mf_mult.set_finite_element(mesh.convex_index(), 
			       getfem::fem_descriptor(dirichlet_fem_name));
    cout << "DIRICHLET_DER_FEM_TYPE="  << dirichlet_derivative_fem_name << "\n";
    mf_mult_d.set_finite_element(mesh.convex_index(), 
			       getfem::fem_descriptor(dirichlet_derivative_fem_name));
  }

  /* set the finite element on mf_rhs (same as mf_u if DATA_FEM_TYPE is
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
    mesh.region(SIMPLE_SUPPORT_BOUNDARY_NUM).add(i.cv(), i.f());
    mesh.region(CLAMPED_BOUNDARY_NUM).add(i.cv(), i.f()); 
  }



}

/* compute the relative error with respect to the exact solution */
void bilaplacian_mortar_problem::compute_error(plain_vector &U) {
  std::vector<scalar_type> V(mf_rhs.nb_dof());
  getfem::interpolation(mf_u, mf_rhs, U, V);
  for (size_type i = 0; i < mf_rhs.nb_dof(); ++i)
    V[i] -= sol_u(mf_rhs.point_of_dof(i));
  cout.precision(16);
//   cout  << "L2 error = " << getfem::asm_L2_norm(mim, mf_rhs, V)  << endl
//         << "H1 error = " << getfem::asm_H1_norm(mim, mf_rhs, V)  << endl
//         << "H2 error = " << getfem::asm_H2_norm(mim, mf_rhs, V)  << endl
//         /*<< "Linfty error = " << gmm::vect_norminf(V)  << endl*/; 
//   cout  << "semi-norme H1 = " << getfem::asm_H1_semi_norm(mim, mf_rhs, V)  << endl 
//         << "semi-norme H2 = " << getfem::asm_H2_semi_norm(mim, mf_rhs, V)  << endl ;
cout << "erreur L2 / errur H1 / erreur H2 / semi-H1 / semi-H2 :\n" << getfem::asm_L2_norm(mim, mf_rhs, V) << " ";
cout << getfem::asm_H1_norm(mim, mf_rhs, V) << " " << getfem::asm_H2_norm(mim, mf_rhs, V) << " ";
cout << getfem::asm_H1_semi_norm(mim, mf_rhs, V) << " " <<getfem::asm_H2_semi_norm(mim, mf_rhs, V) << endl ;
       
}

/* The aim of this function is to provide the errors without 
   assembling the norms, as it is some expensive computations.
   The error is computed after interpolation on a lagrangian 
   mesh, so the error can be calculated by a single substraction. 
   This method is only useful for plotting the 2D graph of the error, 
   when we don't need the L2 or H1 norm evaluation of the error. */
void bilaplacian_mortar_problem::compute_error_without_assembling(plain_vector &U, plain_vector &V) {
  gmm::resize( V, mf_rhs.nb_dof());
  getfem::interpolation(mf_u, mf_rhs, U, V);
  for (size_type i = 0; i < mf_rhs.nb_dof(); ++i)
    V[i] -= sol_u(mf_rhs.point_of_dof(i));

}



/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool bilaplacian_mortar_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();
   

    
    cout << "\nIntegral Matching (Mortar)\n" ;
    
    dal::bit_vector cvlist_in_area;
    dal::bit_vector cvlist_out_area;
    for (dal::bv_visitor cv(mesh.convex_index()); 
	   !cv.finished(); ++cv) {
	bool in_area = true;
	/* For each element, we test all of its nodes. 
	   If all the nodes are inside the enrichment area,
	   then the element is completly inside the area too */ 
	for (unsigned j=0; j < mesh.nb_points_of_convex(cv); ++j) {
	  if (gmm::sqr(mesh.points_of_convex(cv)[j][0] ) + 
	      gmm::sqr(mesh.points_of_convex(cv)[j][1] ) > 
	      gmm::sqr(radius_mortar)) {
	    in_area = false; break;
	  }
	}


	if (!in_area) {
	  cvlist_out_area.add(cv);
	  mf_u.set_dof_partition(cv, 1);
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
     

      cout << "cvlist_in_area: " << cvlist_in_area << "\n";
      cout << "MORTAR_BOUNDARY_IN: " << mesh.region(MORTAR_BOUNDARY_IN) << "\n";
      cout << "MORTAR_BOUNDARY_OUT: " << mesh.region(MORTAR_BOUNDARY_OUT) << "\n";
    

  
  // the model bricks ---------------------------------------------------------------------------
  
  cout << "Number of dof for u: " << mf_u.nb_dof() << endl;

  // Bilaplacian brick.
  getfem::mdbrick_bilaplacian<> BIL(mim, mf_u);
  BIL.D().set(D);
  if (KL) { BIL.set_to_KL(); BIL.nu().set(nu); }
  
  // Volumic source term
  plain_vector F(nb_dof_rhs);
  getfem::interpolation_function(mf_rhs, F, sol_f);
  getfem::mdbrick_source_term<> VOL_F(BIL, mf_rhs, F);
  
  //Defining the normal derivative Dirichlet condition value.
  
  F.resize(nb_dof_rhs*N);
  getfem::interpolation_function(mf_rhs, F, sol_du, CLAMPED_BOUNDARY_NUM);    
  getfem::mdbrick_normal_derivative_Dirichlet<> NDER_DIRICHLET(VOL_F, CLAMPED_BOUNDARY_NUM, mf_mult_d); 
  NDER_DIRICHLET.set_constraints_type(dirichlet_version);
  NDER_DIRICHLET.rhs().set(mf_rhs, F);
  
  
//   F.resize(nb_dof_rhs);
//   getfem::interpolation_function(mf_rhs, F, sol_u, CLAMPED_BOUNDARY_NUM);    
//   getfem::mdbrick_normal_derivative_Dirichlet<> NDER_DIRICHLET(VOL_F, CLAMPED_BOUNDARY_NUM, mf_mult_d); 
//   NDER_DIRICHLET.set_constraints_type(dirichlet_version);
//   NDER_DIRICHLET.R_must_be_derivated(); // hence we give the exact solution , and its gradient will be taken
//   NDER_DIRICHLET.rhs().set(mf_rhs, F);
  
  
  // Defining the Dirichlet condition value.
  gmm::resize(F, nb_dof_rhs);
  getfem::interpolation_function(mf_rhs, F, sol_u,SIMPLE_SUPPORT_BOUNDARY_NUM);

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<>
    DIRICHLET(NDER_DIRICHLET, SIMPLE_SUPPORT_BOUNDARY_NUM, mf_mult);
  DIRICHLET.rhs().set(mf_rhs, F);
  DIRICHLET.set_constraints_type(getfem::constraints_type(dirichlet_version));  
  if (dirichlet_version == getfem::PENALIZED_CONSTRAINTS)
    DIRICHLET.set_penalization_parameter(PARAM.real_value("EPS_DIRICHLET_PENAL")) ;
 
     /* add a constraint brick for the mortar junction between
       the enriched area and the rest of the mesh */
    getfem::mdbrick_constraint<> MORTAR(DIRICHLET, 0);
    MORTAR.set_constraints_type(getfem::constraints_type(mortar_version));  
       if (mortar_version == getfem::PENALIZED_CONSTRAINTS)
    MORTAR.set_penalization_parameter(PARAM.real_value("EPS_MORTAR_PENAL")) ;

   sparse_matrix H0, H ;
   getfem::base_vector R; 
   
   if (PARAM.int_value("MORTAR_MATCHING") == 1)	{
    /* older version of integral matching (jan-feb 2007) :
     * \int_\Gamma (u-v) \lambda + \partial_n (u+v) \partial_n \lambda = 0, for all \lambda \in \Lambda
     */
          
    /* build the list of dof for the "(u-v) lambda" and for the  
    "\partial_n(u+v) \partial_n lambda" term in the mortar condition */  
    dal::bit_vector bv_mortar;
    dal::bit_vector bv_deriv;
    dal::bit_vector bv_union;
    sparse_matrix MM(mf_mortar.nb_dof(), mf_mortar.nb_dof());
    sparse_matrix MD(mf_mortar.nb_dof(), mf_mortar.nb_dof());
    std::vector<size_type> ind_mortar;
    std::vector<size_type> ind_deriv;
    std::vector<size_type> ind_union;
    getfem::asm_mass_matrix(MM, mim, mf_mortar, MORTAR_BOUNDARY_OUT);

    gmm::resize(R, gmm::mat_nrows(MD));
    getfem::asm_normal_derivative_dirichlet_constraints_bis(MD, R, mim, mf_mortar,
	     		mf_mortar, mf_mortar, R, MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;

    for (size_type i=0; i < mf_mortar.nb_dof(); ++i) {
      if ( (MM(i,i) > 1e-15) & (MD(i,i) > 1e-15) ) { 
	bv_mortar.add(i);
	bv_deriv.add(i);
	bv_union.add(i);
      }
      if ( (MM(i,i) > 1e-15) & (MD(i,i) <= 1e-15) ) { 
	bv_mortar.add(i);
	bv_union.add(i);
      }
      if ( (MM(i,i) <= 1e-15) & (MD(i,i) > 1e-15) ) { 
	bv_deriv.add(i);
	bv_union.add(i);
      }
    }

    //cout << "matrice des dérivées normales du mesh_fem_mortar : \n" << MD << "\n" ;

    for (dal::bv_visitor d(bv_mortar); !d.finished(); ++d)
      ind_mortar.push_back(d);
    for (dal::bv_visitor d(bv_deriv); !d.finished(); ++d)
      ind_deriv.push_back(d);
    for (dal::bv_visitor d(bv_union); !d.finished(); ++d)
      ind_union.push_back(d);

  /* This command was not suitable here :
    dal::bit_vector bv_mortar = 
      mf_mortar.dof_on_region(MORTAR_BOUNDARY_OUT);  
      The reason for that is that unfortunately, the method "dof_on_region"
      sometimes return too much dof when the mesh_fem is enriched. */
  



    /* building matrices */
    cout << "Handling mortar junction (" << ind_union.size() << 
      " dof for the lagrange multiplier)\n";
      
    gmm::resize(H0, mf_mortar.nb_dof(), mf_u.nb_dof()) ;
    gmm::resize(H , ind_union.size(),   mf_u.nb_dof());

    cout << "bv_mortar = " << bv_mortar << "\n";
    cout << "bv_deriv = " << bv_deriv << "\n" ;
    cout << "bv_union = " << bv_union << "\n" ;
    
    gmm::sub_index sub_i(ind_union);
    gmm::sub_index sub_i1(ind_mortar);
    gmm::sub_index sub_i2(ind_deriv);
    gmm::sub_interval sub_j(0, mf_u.nb_dof());
    // build sub_indices of dofs which are either value or derivatives in the matrix of constraints H
    std::vector<size_type> ind_val_H, ind_deriv_H ;
    for (unsigned i=0; i< ind_union.size(); ++i) {
      if ( bv_mortar[ind_union[i]] ) ind_val_H.push_back(i) ;
      if ( bv_deriv[ind_union[i]] ) ind_deriv_H.push_back(i) ;
    }   
    gmm::sub_index sub_val_H(ind_val_H) ;
    gmm::sub_index sub_deriv_H(ind_deriv_H) ;

    /* build the mortar constraint matrix -- note that the integration
       method is conformal to the crack
     */
    getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u, MORTAR_BOUNDARY_OUT);
    gmm::copy(gmm::sub_matrix(H0, sub_i1, sub_j), gmm::sub_matrix(H, sub_val_H, sub_j) );

    gmm::clear(H0);
    getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u, MORTAR_BOUNDARY_IN);
    gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i1, sub_j), -1), gmm::sub_matrix(H, sub_val_H, sub_j));
    
   
    gmm::clear(H0) ;
    gmm::resize(R, mat_nrows(H));
    getfem::asm_normal_derivative_dirichlet_constraints_bis(H0, R, mim, mf_u,
	     		mf_mortar, mf_u, R, MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;
    gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i2, sub_j), 1. ), gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
    
    gmm::clear(H0) ;
    getfem::asm_normal_derivative_dirichlet_constraints_bis(H0, R, mim, mf_u,
			mf_mortar, mf_u, R, MORTAR_BOUNDARY_IN, 0, getfem::ASMDIR_BUILDH) ;
    gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i2, sub_j), -1.), gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
    
    // -----------------------------------
    }
    
    if (PARAM.int_value("MORTAR_MATCHING") == 2 ){
    /* Another version of the integral matching :
     *     \int_Gamma        (u-v) \lambda  = 0, for all \lambda in \Lambda
     *     \int_Gamma \nabla (u-v).\mu  = 0, for all \mu in M 
     */

    // Be carefull : the multiplier is vectorial.
    mf_mortar_deriv.set_qdim(2) ;
    
    // selecting nodes indices on the two meth. mult.
    // assembling matrices, selecting the lines corresponding to those who are holding some information

   /* build the list of dof for the "(u-v) lambda" condition  */  
    dal::bit_vector bv_mortar;
    dal::bit_vector bv_deriv;
    sparse_matrix MM(mf_mortar.nb_dof(), mf_mortar.nb_dof());
    sparse_matrix MD(mf_mortar_deriv.nb_dof(), mf_mortar_deriv.nb_dof());
    std::vector<size_type> ind_mortar;
    std::vector<size_type> ind_deriv;
    getfem::asm_mass_matrix(MM, mim, mf_mortar, MORTAR_BOUNDARY_OUT);
    getfem::asm_mass_matrix(MD, mim, mf_mortar_deriv, MORTAR_BOUNDARY_OUT);
    //  getfem::base_vector R( mf_mortar_deriv.nb_dof() );
    //  getfem::asm_normal_derivative_dirichlet_constraints_bis(MD, R, mim, mf_mortar_deriv,
    // 	  		mf_mortar_deriv, mf_u, R, MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;
    for (size_type i=0; i < mf_mortar.nb_dof(); ++i)
      if (MM(i,i) > 1e-15) bv_mortar.add(i);
    for (size_type i=0; i < mf_mortar_deriv.nb_dof(); ++i)
      if (MD(i,i) > 1e-15) bv_deriv.add(i);
    
    for (dal::bv_visitor d(bv_mortar); !d.finished(); ++d)
      ind_mortar.push_back(d);
    for (dal::bv_visitor d(bv_deriv); !d.finished(); ++d)
      ind_deriv.push_back(d);

    // building matrices

    cout << "Handling mortar junction (" << ind_mortar.size() << 
      " dof for the lagrange multiplier of the displacement, " <<
      ind_deriv.size() << " dof for the lagrange multiplier of the derivative)\n";
      
    gmm::resize(H0, mf_mortar.nb_dof(), mf_u.nb_dof()) ;
    gmm::resize(H,  ind_mortar.size() + ind_deriv.size(), mf_u.nb_dof()) ; 

    // Defining sub_indexes of the matrices calculated with the 
    // complete set of dofs.
    gmm::sub_index sub_i(ind_mortar);
    gmm::sub_index sub_i1(ind_deriv);
    gmm::sub_interval sub_j(0, mf_u.nb_dof());

    gmm::sub_interval sub_val_H(0, ind_mortar.size()) ;
    gmm::sub_interval sub_deriv_H(ind_mortar.size(), ind_deriv.size()) ;

    cout << "sub_indexes built\n" ;
    /* build the mortar constraint matrix -- note that the integration
       method is conformal to the crack
     */
    getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u, MORTAR_BOUNDARY_OUT);
    gmm::copy(gmm::sub_matrix(H0, sub_i, sub_j), gmm::sub_matrix(H, sub_val_H, sub_j) );

    gmm::clear(H0);
    getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u, MORTAR_BOUNDARY_IN);
    gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i, sub_j), -1), gmm::sub_matrix(H, sub_val_H, sub_j) );


    cout << "first contraint asm\n" ;
    
    gmm::clear(H0);
    gmm::resize(H0, mf_mortar_deriv.nb_dof(), mf_u.nb_dof() ) ;
    getfem::asm_constraint_gradient_vectorial_mult
	(H0, mim, mf_u, mf_mortar_deriv,
	 MORTAR_BOUNDARY_OUT, getfem::ASMDIR_BUILDH) ;       
    gmm::add(gmm::sub_matrix(H0, sub_i1, sub_j),
             gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
   
    cout << "first step \n" ;
    gmm::clear(H0);
    getfem::asm_constraint_gradient_vectorial_mult
	(H0, mim, mf_u, mf_mortar_deriv,
	 MORTAR_BOUNDARY_IN, getfem::ASMDIR_BUILDH) ;
    gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i1, sub_j), -1),
             gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
    
    }
    
    if (PARAM.int_value("MORTAR_MATCHING") == 3 ){
       
    /* Other version of the integral matching :
     *     \int_Gamma        (u-v) \lambda  = 0, for all \lambda in \Lambda
     *     \int_Gamma \partial_n (u+v)\mu  = 0, for all \mu in M 
     */

// selecting nodes indices on the two meth. mult.

// assembling matrices, selecting the lines corresponding to those who are holding some information

   /* build the list of dof for the "(u-v) lambda" condition  */  
    dal::bit_vector bv_mortar;
    dal::bit_vector bv_deriv;
    sparse_matrix MM(mf_mortar.nb_dof(), mf_mortar.nb_dof());
    sparse_matrix MD(mf_mortar_deriv.nb_dof(), mf_mortar_deriv.nb_dof());
    std::vector<size_type> ind_mortar;
    std::vector<size_type> ind_deriv;
    getfem::asm_mass_matrix(MM, mim, mf_mortar, MORTAR_BOUNDARY_OUT);
    getfem::asm_mass_matrix(MD, mim, mf_mortar_deriv, MORTAR_BOUNDARY_OUT);
    //  getfem::base_vector R( mf_mortar_deriv.nb_dof() );
    //  getfem::asm_normal_derivative_dirichlet_constraints_bis(MD, R, mim, mf_mortar_deriv,
    // 	  		mf_mortar_deriv, mf_u, R, MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;
    for (size_type i=0; i < mf_mortar.nb_dof(); ++i)
      if (MM(i,i) > 1e-15) bv_mortar.add(i);
    for (size_type i=0; i < mf_mortar_deriv.nb_dof(); ++i)
      if (MD(i,i) > 1e-15) bv_deriv.add(i);
    
    for (dal::bv_visitor d(bv_mortar); !d.finished(); ++d)
      ind_mortar.push_back(d);
    for (dal::bv_visitor d(bv_deriv); !d.finished(); ++d)
      ind_deriv.push_back(d);

    // building matrices

    cout << "Handling mortar junction (" << ind_mortar.size() << 
      " dof for the lagrange multiplier of the displacement, " <<
      ind_deriv.size() << " dof for the lagrange multiplier of the derivative)\n";
      
    gmm::resize(H0, mf_mortar.nb_dof(), mf_u.nb_dof()) ;
    gmm::resize(H,  ind_mortar.size() + ind_deriv.size(), mf_u.nb_dof()) ; 

    // Defining sub_indexes of the matrices calculated with the 
    // complete set of dofs.
    gmm::sub_index sub_i(ind_mortar);
    gmm::sub_index sub_i1(ind_deriv);
    gmm::sub_interval sub_j(0, mf_u.nb_dof());

    gmm::sub_interval sub_val_H(0, ind_mortar.size()) ;
    gmm::sub_interval sub_deriv_H(ind_mortar.size(), ind_deriv.size()) ;

    cout << "sub_indexes built\n" ;
    /* build the mortar constraint matrix -- note that the integration
       method is conformal to the crack
     */
    getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u, MORTAR_BOUNDARY_OUT);
    gmm::copy(gmm::sub_matrix(H0, sub_i, sub_j), gmm::sub_matrix(H, sub_val_H, sub_j) );

    gmm::clear(H0);
    getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u, MORTAR_BOUNDARY_IN);
    gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i, sub_j), -1), 
             gmm::sub_matrix(H, sub_val_H, sub_j) );


    cout << "first contraint asm\n" ;
    
    gmm::clear(H0);
    gmm::resize(H0, mf_mortar_deriv.nb_dof(), mf_u.nb_dof() ) ;
    getfem::asm_normal_derivative_dirichlet_constraints
      (H0, R, mim, mf_u, mf_mortar_deriv, mf_u, R,
       MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;
    gmm::add(gmm::sub_matrix(H0, sub_i1, sub_j),
             gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
   
    cout << "first step \n" ;
    gmm::clear(H0);
    getfem::asm_normal_derivative_dirichlet_constraints
      (H0, R, mim, mf_u, mf_mortar_deriv, mf_u, R,
       MORTAR_BOUNDARY_IN, 0, getfem::ASMDIR_BUILDH) ;
    gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i1, sub_j), 1),
             gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
    
    }
    
    // ------------------------------------------ end of new version

    /* because of the discontinuous partition of mf_u, some levelset 
       enriched functions do not contribute any more to the
       mass-matrix (the ones which are null on one side of the
       levelset, when split in two by the mortar partition, may create
       a "null" dof whose base function is all zero.. */
    sparse_matrix M2(mf_u.nb_dof(), mf_u.nb_dof());
    getfem::asm_mass_matrix(M2, mim, mf_u, mf_u);
    //gmm::HarwellBoeing_IO::write("M2.hb", M2);
    for (size_type d = 0; d < mf_u.nb_dof(); ++d) {
      if (M2(d,d) < PARAM.real_value("SEUIL")) {
	cout << "  removing null mf_u dof " << d << " @ " << 
	  mf_u.point_of_dof(d) << "\n";	
	unsigned n = gmm::mat_nrows(H);
	gmm::resize(H, n+1, gmm::mat_ncols(H));
	H(n, d) = 1;
      }
    }  
    gmm::resize(R, gmm::mat_nrows(H)); 
    MORTAR.set_constraints(H,R);
//    gmm::HarwellBoeing_IO::write("H.hb", H);        
    //------------------------------------------------------------------
    // Matching of the normal derivative
//     getfem::mdbrick_constraint<> &mortar_derivative = 
//       *(new getfem::mdbrick_constraint<>(mortar,0));
//       
//     gmm::clear(H0) ;
//     gmm::clear(H) ;
//     gmm::resize(H, ind_mortar.size(), mf_u.nb_dof());
//     getfem::asm_normal_derivative_dirichlet_constraints(H0, R, mim, mf_u,
// 	     		mf_mortar_deriv, mf_u, R, MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;
//     gmm::copy(gmm::sub_matrix(H0, sub_i, sub_j), H) ;
//     
//     gmm::clear(H0) ;
//     getfem::asm_normal_derivative_dirichlet_constraints(H0, R, mim, mf_u,
// 			mf_mortar_deriv, mf_u, R, MORTAR_BOUNDARY_IN, 0, getfem::ASMDIR_BUILDH) ;
//     gmm::add(gmm::sub_matrix(H0, sub_i, sub_j), H) ;
//     
//     //cout << "matrice de contraintes : \n" << H << "\n" ;
//     cout << "matching of derivative values done.\n" ;  
//     gmm::resize(R, ind_mortar.size());
//     mortar_derivative.set_constraints(H,R);
//     final_model = &mortar_derivative ;
  
  

  // Generic solve.
  cout << "Total number of variables : " << MORTAR.nb_dof() << endl;
  getfem::standard_model_state MS(MORTAR);
  gmm::iteration iter(residual, 1, 40000);
  getfem::standard_solve(MS, MORTAR, iter);

  // Solution extraction
  gmm::resize(U, mf_u.nb_dof());
  gmm::copy(BIL.get_solution(MS), U);
  return (iter.converged());
}






