/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2006-2012 Yves Renard, Julien Pommier, Jeremie Lasry.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/

#include "crack_bilaplacian.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_fourth_order.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_superlu.h"



scalar_type moment = - 1. ;
scalar_type bbeta = M_PI / 4. ;
//scalar_type epsilon = 1e-1 ;

base_matrix sol_moment(const base_node &x)
{ base_matrix m(x.size(), x.size()); 
  m(1,1) = moment;
  return m; }
  
base_small_vector sol_ff(const base_node &x)
{ base_small_vector res(x.size());
  if ( x[1] > 0.4999){
     res[1] = moment ; //* (-1./3.) * epsilon * epsilon;
  }
  if ( x[1] < -0.4999){
     res[1] = - moment ; //* (-1./3.) * epsilon * epsilon;
  }  
  res[0] = 0. ;
  return res ;}

scalar_type sol_beta(const base_node &x){
base_matrix m(2,2) ;
m(0,0) = m(1,1) = cos(bbeta) ;
m(0,1) = - sin(bbeta) ;
m(1,0) = sin(bbeta) ;
base_node z(2) ;
gmm::mult(m, x, z) ;
return - moment * cos(bbeta) * sin(bbeta) * z[0] * z[1] / 0.7 ; }

base_matrix moment_beta(const base_node &x)
{ base_matrix m(x.size(), x.size()); 
  m(1,0) = m(0,1) = moment * sin(bbeta) * cos(bbeta) ;
  return m; }

  
/***********************************/

size_type is_lagrange_dof_type(getfem::pdof_description dof){
  size_type displ_dof = 0 ;
  for (dim_type k = 0; k < 4; ++k) {
      if (dof == getfem::lagrange_dof(k)) 
         displ_dof = 1 ;
    }
return displ_dof ;
}

size_type is_dx_dof_type(getfem::pdof_description dof){
  char s[200];
  size_type dx_dof = 0 ;
    for (unsigned r = 0; r < 2; ++r) {
      if (dof == getfem::derivative_dof(2,0)){
         sprintf(s, "D_%c[%d]", "xyzuvw"[r], 2);
         dx_dof = 1 ;
	 }
  }
return dx_dof ;
}

size_type is_dy_dof_type(getfem::pdof_description dof){
  char s[200];
  size_type dy_dof = 0 ;
    for (unsigned r = 0; r < 2; ++r) {
      if (dof == getfem::derivative_dof(2,1)){
         sprintf(s, "D_%c[%d]", "xyzuvw"[r], 2);
         dy_dof = 1 ;
	 cout << "passage dans is_dy_dof_type OK \n" ;
	 }
    }
return dy_dof ;
}

bool bilaplacian_crack_problem::solve_moment(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();

  // setting singularities 
  cout << "setting singularities \n" ;
  if (PARAM.int_value("SING_BASE_TYPE") == 0){
	std::vector<getfem::pglobal_function> ufunc(4);
	for (size_type i = 0 ; i < ufunc.size() ; ++i) {
	ufunc[i] = bilaplacian_crack_singular(i, ls, nu, 0.);
	}
  mf_sing_u.set_functions(ufunc);
  }
  if (PARAM.int_value("SING_BASE_TYPE") == 1){
	std::vector<getfem::pglobal_function> ufunc(2);
	for (size_type i = 0 ; i < ufunc.size() ; ++i) {
	ufunc[i] = bilaplacian_crack_singular(i + 4, ls, nu, 0.);
	}
  mf_sing_u.set_functions(ufunc);
  }
  
  
  // Setting the enrichment --------------------------------------------/
   
  switch(enrichment_option) {
  case 0 :  // No enrichment
    mf_u_sum.set_mesh_fems(mfls_u);
    break ;
  case 1 : 
    {
      cout << "\npointwise matching\n";
     /* first : selecting the convexes that are completly included in the enrichment area */ 
     for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i) {
       pm_convexes.add(i) ;
	/* For each element, we test all of its nodes. 
	   If all the nodes are inside the enrichment area,
	   then the element is completly inside the area too */ 
       for (unsigned j=0; j < mesh.nb_points_of_convex(i); ++j) {
	 if (gmm::sqr(mesh.points_of_convex(i)[j][0]) + 
	     gmm::sqr(mesh.points_of_convex(i)[j][1]) > 
	     gmm::sqr(enr_area_radius)) 
	   pm_convexes.sup(i); break;
       }
     }
      
      for (dal::bv_visitor cv(mf_sing_u.convex_index()); !cv.finished(); ++cv) {
	if (!pm_convexes.is_in(cv))
	  mf_sing_u.set_finite_element(cv, 0);
      }
      cout << "mf_sing_u: convex_index() = " << mf_sing_u.convex_index().card() << " convexes\n";

      //mf_u_sum.set_mesh_fems(mfls_u_ext, mf_pre_u); //_ext, mf_sing_u);
      mf_u_sum.set_smart_global_dof_linking(true);
      mf_u_sum.set_mesh_fems(mf_pre_u, mf_sing_u);


      cout << "mf_u_sum.nb_dof = " << mf_u_sum.nb_dof() << "\n";
      cout << "mfls_u.convex_index = " << mfls_u.convex_index() << "\nmf_sing_u: " << mf_sing_u.convex_index() << "\n";
      
    } break ;
  case 2 :  // standard XFEM on a fixed zone
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
      //cout << "enriched_dofs: " << enriched_dofs << "\n";
      if (enriched_dofs.card() < 3)
            GMM_WARNING0("There is " << enriched_dofs.card() <<
		   " enriched dofs for the crack tip");
      mf_u_product.set_enrichment(enriched_dofs);
      mf_u_sum.set_mesh_fems(mf_u_product, mfls_u);
      cout << "enrichment done \n" ;}
      break ;
  case 3 : // Integral matching (mortar)
    {
    cout << "\nIntegral Matching (Mortar)\n" ;    
    dal::bit_vector cvlist_in_area, cvlist_out_area;
    bool in_area = true;
    for (dal::bv_visitor cv(mesh.convex_index()); 
	   !cv.finished(); ++cv) {
	in_area = true;
	/* For each element, we test all of its nodes. 
	   If all the nodes are inside the enrichment area,
	   then the element is completly inside the area too */ 
	for (unsigned j=0; j < mesh.nb_points_of_convex(cv); ++j) {
	  if (gmm::sqr(mesh.points_of_convex(cv)[j][0] ) + 
	      gmm::sqr(mesh.points_of_convex(cv)[j][1] ) > 
	      gmm::sqr(enr_area_radius)) {
	          in_area = false; 
		  break;
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
      if (PARAM.int_value("MORTAR_WITHOUT_SINGUL"))
         mf_u_sum.set_mesh_fems(mfls_u);
      else
         mf_u_sum.set_mesh_fems(mf_sing_u, mfls_u);
      
      //cout << "cvlist_in_area: " << cvlist_in_area << "\n";
      cout << "mfls_u.nb_dof: " << mfls_u.nb_dof() << "\n";
      cout << "mf_u_sum.nb_dof: " << mf_u_sum.nb_dof() << "\n";
      //cout << "MORTAR_BOUNDARY_IN: " << mesh.region(MORTAR_BOUNDARY_IN) << "\n";
      //cout << "MORTAR_BOUNDARY_OUT: " << mesh.region(MORTAR_BOUNDARY_OUT) << "\n";
      
//       // an optional treatment : creating a representation of the enrichment area     
//       getfem::mesh_fem mf_enrich(mesh);
//       getfem::pfem pf_mef = getfem::classical_fem(mesh.trans_of_convex(mesh.convex_index().first_true()), 1 );
//       mf_enrich.set_finite_element(mesh.convex_index(), pf_mef) ;
//       std::vector<scalar_type> UU(mf_enrich.nb_dof()) ;
//       std::fill(UU.begin(), UU.end() ,0.) ;
//       cout << "exporting the enrichment zone: \n" ;
//       for (dal::bv_visitor i(cvlist_in_area) ; !i.finished() ; ++i){ 
// 	  for (unsigned int j = 0 ; j < mf_enrich.ind_dof_of_element(i).size() ; ++j )  
// 	  UU[mf_enrich.ind_dof_of_element(i)[j]] = 1. ;         
//       }
//       
//       cout << "exporting enrichment to " << "enrichment_zone.vtk" << "..\n";
//       getfem::vtk_export exp("enrichment_zone.vtk", false);
//       exp.exporting(mf_enrich); 
//       exp.write_point_data(mf_enrich, UU, "enrichment");
//       cout << "export done, you can view the data file with (for example)\n"
// 	"mayavi -d enrichment_zone.vtk -f "
// 	"WarpScalar -m BandedSurfaceMap -m Outline\n";
	
//       // Another optional treatment :
//       // Searching the elements that are both crossed by the crack
//       // and with one of their faces which constitutes a part of the 
//       // boundary between the enriched zone and the rest of the domain.
//       getfem::mesh_region &boundary = mesh.region(MORTAR_BOUNDARY_IN);
//       unsigned int cpt = 0 ;
//       for (dal::bv_visitor i(cvlist_in_area); !i.finished(); ++i) {
//          if (mls.is_convex_cut(i)){
// 	    // Among the faces of the convex, we search if some are
// 	    // part of the boundary
// 	    cpt = 0 ;
// 	    for (unsigned j=0; j < mesh.structure_of_convex(i) ->nb_faces(); ++j) {
// 	        if (boundary.is_in(i,j))
// 		   cpt += 1;
// 	    }
// 	    if (cpt) {
//                cout << "\n The convex number " << i << " is crossed by the crack :\n" ;
// 	       cout << "  it has : " << cpt << " face(s) among the boundary.\n \n " ;
// 	    }
// 	 }
//       }
 
    }
    break ;
  default : 
	GMM_ASSERT1(false, "Enrichment_option parameter is undefined");
	break ;  
	}
  mesh.write_to_file("toto.mesh");

  cout << "Number of dof for u: " << mf_u().nb_dof() << endl;

  // Bilaplacian brick.
  getfem::mdbrick_bilaplacian<> BIL(mim, mf_u());
  BIL.D().set(D);
  if (KL) { BIL.set_to_KL(); BIL.nu().set(nu); }

  size_type N = mesh.dim();
  plain_vector F(nb_dof_rhs);

  // Defining the prescribed momentum.
  gmm::resize(F, nb_dof_rhs*N*N);
  if (PARAM.int_value("SOL_REF")==1){
     getfem::interpolation_function(mf_rhs, F, sol_moment, MOMENTUM_BOUNDARY_NUM);
  }
  if (PARAM.int_value("SOL_REF")==2){
     getfem::interpolation_function(mf_rhs, F, moment_beta, MOMENTUM_BOUNDARY_NUM);
  }
//   getfem::mdbrick_normal_derivative_source_term<>
//   MOMENTUM(VOL_F, mf_rhs, F, MOMENTUM_BOUNDARY_NUM);
getfem::mdbrick_normal_derivative_source_term<>
  MOMENTUM(BIL, mf_rhs, F, MOMENTUM_BOUNDARY_NUM);

  // Normal derivative Dirichlet condition brick
  gmm::resize(F, nb_dof_rhs*N);
  gmm::clear(F);
  getfem::mdbrick_normal_derivative_Dirichlet<> 
	NDER_DIRICHLET(MOMENTUM, CLAMPED_BOUNDARY_NUM, mf_mult_d);
  NDER_DIRICHLET.set_constraints_type(dirichlet_version);
  NDER_DIRICHLET.rhs().set(mf_rhs, F);

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<>
    DIRICHLET(NDER_DIRICHLET, SIMPLE_SUPPORT_BOUNDARY_NUM, mf_mult);
  DIRICHLET.set_constraints_type(dirichlet_version);
  if (dirichlet_version == getfem::PENALIZED_CONSTRAINTS)
    DIRICHLET.set_penalization_parameter(PARAM.real_value("EPS_DIRICHLET_PENAL")) ;
  gmm::resize(F, nb_dof_rhs);
  gmm::clear(F) ;
  DIRICHLET.rhs().set(mf_rhs, F);

//   // model brick constraint for the integral matching
//   getfem::mdbrick_constraint<> &mortar = 
//       *(new getfem::mdbrick_constraint<>(NDER_DIRICHLET, 0));

  getfem::mdbrick_abstract<> *final_model = &DIRICHLET ;

  if (enrichment_option == 3 ) {
     /* add a constraint brick for the mortar junction between
       the enriched area and the rest of the mesh */

    getfem::mdbrick_constraint<> &mortar = 
      *(new getfem::mdbrick_constraint<>(DIRICHLET, 0));
    mortar.set_constraints_type(getfem::constraints_type(mortar_version));
    if (mortar_version == getfem::PENALIZED_CONSTRAINTS)
      mortar.set_penalization_parameter(PARAM.real_value("EPS_MORTAR_PENAL")) ;

    plain_vector R(1) ;
    sparse_matrix H(1, mf_u().nb_dof());
    (*this).set_matrix_mortar(H) ;

    /* because of the discontinuous partition of mf_u(), some levelset 
       enriched functions do not contribute any more to the
       mass-matrix (the ones which are null on one side of the
       levelset, when split in two by the mortar partition, may create
       a "null" dof whose base function is all zero.. */
    sparse_matrix M2(mf_u().nb_dof(), mf_u().nb_dof());
    getfem::asm_mass_matrix(M2, mim, mf_u(), mf_u());
    //gmm::HarwellBoeing_IO::write("M2.hb", M2);
    cout << "PARAM.real_value(\"SEUIL\") : " << PARAM.real_value("SEUIL") << "\n" ;
    for (size_type d = 0; d < mf_u().nb_dof(); ++d) {
      //       if (M2(d,d) < 1e-7) cout << "  weak mf_u() dof " << d << " @ " << 
      // 	  mf_u().point_of_dof(d) << " M2(d,d) = " << M2(d,d) << "\n";
      if (M2(d,d) < PARAM.real_value("SEUIL")) {
	cout << "removed\n";	
	size_type n = gmm::mat_nrows(H);
	gmm::resize(H, n+1, gmm::mat_ncols(H));
	H(n, d) = 1;
      }
    }  
    gmm::resize(R, gmm::mat_nrows(H)); 
    mortar.set_constraints(H,R);
    final_model = &mortar;
    gmm::Harwell_Boeing_save("H.hb", H);

  }

    cout << "Total number of variables : " << final_model->nb_dof() << endl;
    getfem::standard_model_state MS(*final_model);
    gmm::iteration iter(residual, 1, 40000);
    getfem::standard_solve(MS, *final_model, iter /* , p*/);  
    // getfem::useful_types<getfem::standard_model_state>::plsolver_type p;
    // p.reset(new getfem::linear_solver_cg_preconditioned_ildlt<sparse_matrix,plain_vector>);
    // Solution extraction
    gmm::resize(U, mf_u().nb_dof());
    gmm::copy(BIL.get_solution(MS), U); 

  return true;  
}

void bilaplacian_crack_problem::compute_error_beta(plain_vector &U) {
  GMM_ASSERT1(!mf_rhs.is_reduced(), "To be adapted");
  std::vector<scalar_type> V(mf_rhs.nb_dof());
  getfem::interpolation(mf_u(), mf_rhs, U, V);
  for (size_type i = 0; i < mf_rhs.nb_dof(); ++i)
    V[i] -= sol_beta(mf_rhs.point_of_basic_dof(i));
  cout.precision(8);
//   cout  << "L2 error = " << getfem::asm_L2_norm(mim, mf_rhs, V)  << endl
//         << "H1 error = " << getfem::asm_H1_norm(mim, mf_rhs, V)  << endl
//         << "H2 error = " << getfem::asm_H2_norm(mim, mf_rhs, V)  << endl
//         /*<< "Linfty error = " << gmm::vect_norminf(V)  << endl*/; 
//   cout  << "semi-norme H1 = " << getfem::asm_H1_semi_norm(mim, mf_rhs, V)  << endl 
//         << "semi-norme H2 = " << getfem::asm_H2_semi_norm(mim, mf_rhs, V)  << endl ;
cout << "erreur L2 / erreur H1 / erreur H2 / semi-H1 / semi-H2 :\n" << getfem::asm_L2_norm(mim, mf_rhs, V) << " ";
cout << getfem::asm_H1_norm(mim, mf_rhs, V) << " " << getfem::asm_H2_norm(mim, mf_rhs, V) << " ";
cout << getfem::asm_H1_semi_norm(mim, mf_rhs, V) << " " <<getfem::asm_H2_semi_norm(mim, mf_rhs, V) << endl ;
  if(PARAM.real_value("NORM_EXACT") != 0. ){
    for (size_type i = 0; i < mf_rhs.nb_dof(); ++i)
        V[i] = sol_beta(mf_rhs.point_of_basic_dof(i));
    cout << "display exact solution: \n" ; 
    cout << "L2 norm: " <<  getfem::asm_L2_norm(mim, mf_rhs, V) << endl ;
    cout << "H1 norm: " <<  getfem::asm_H1_norm(mim, mf_rhs, V) << endl ;
    cout << "H2 norm: " <<  getfem::asm_H2_norm(mim, mf_rhs, V) << endl ;
    }
}
