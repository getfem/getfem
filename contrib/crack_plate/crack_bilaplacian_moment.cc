#include "crack_bilaplacian.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_fourth_order.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_superlu.h"



scalar_type moment = 1. ;
scalar_type epsilon = 1e-2 ;

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

size_type is_lagrange_dof_type(getfem::pdof_description dof){
  size_type displ_dof = 0 ;
  for (unsigned k = 0; k < 4; ++k) {
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
  
  // Setting the level-set
  ls.reinit();  
  scalar_type a = PARAM.real_value("CRACK_SEMI_LENGTH") ;
  for (size_type d = 0; d < ls.get_mesh_fem().nb_dof(); ++d) {
    scalar_type x = ls.get_mesh_fem().point_of_dof(d)[0];
    scalar_type y = ls.get_mesh_fem().point_of_dof(d)[1];
    if (PARAM.int_value("SOL_REF") == 0){
    ls.values(0)[d] = y  ; // + 1/4.*(x + .25);
    ls.values(1)[d] = x;}
    if (PARAM.int_value("SOL_REF") == 1){
     ls.values(0)[d] = y  ;
     ls.values(1)[d] = x - a ; //x * x - a * a ;
    }
  }
  //ls.simplify(0.5);
  ls.touch();  
  mls.adapt();
  mim.adapt();
  mfls_u.adapt();
  mfls_mortar.adapt();
  mfls_mortar_deriv.adapt();
  cout << "mfls_u.nb_dof()=" << mfls_u.nb_dof() << "\n";
  
  // setting singularities 
  cout << "setting singularities \n" ;
  if (PARAM.int_value("SING_BASE_TYPE") == 0){
	std::vector<getfem::pglobal_function> ufunc(4);
	for (size_type i = 0 ; i < ufunc.size() ; ++i) {                              
	ufunc[i] = bilaplacian_crack_singular(i, ls, nu);
	}
  mf_sing_u.set_functions(ufunc);
  }
  if (PARAM.int_value("SING_BASE_TYPE") == 1){
	std::vector<getfem::pglobal_function> ufunc(2);
	for (size_type i = 0 ; i < ufunc.size() ; ++i) {                              
	ufunc[i] = bilaplacian_crack_singular(i + 4, ls, nu);
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
      
      // an optional treatment : creating a representation of the enrichment area     
      getfem::mesh_fem mf_enrich(mesh);
      getfem::pfem pf_mef = getfem::classical_fem(mesh.trans_of_convex(mesh.convex_index().first_true()), 1 );
      mf_enrich.set_finite_element(mesh.convex_index(), pf_mef) ;
      std::vector<scalar_type> UU(mf_enrich.nb_dof()) ;
      std::fill(UU.begin(), UU.end() ,0.) ;
      cout << "exporting the enrichment zone: \n" ;
      for (dal::bv_visitor i(cvlist_in_area) ; !i.finished() ; ++i){ 
	  for (unsigned int j = 0 ; j < mf_enrich.ind_dof_of_element(i).size() ; ++j )  
	  UU[mf_enrich.ind_dof_of_element(i)[j]] = 1. ;         
      }
      
      cout << "exporting enrichment to " << "enrichment_zone.vtk" << "..\n";
      getfem::vtk_export exp("enrichment_zone.vtk", false);
      exp.exporting(mf_enrich); 
      exp.write_point_data(mf_enrich, UU, "enrichment");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d enrichment_zone.vtk -f "
	"WarpScalar -m BandedSurfaceMap -m Outline\n";
	
      // Another optional treatment :
      // Searching the elements that are both crossed by the crack
      // and with one of their faces which constitutes a part of the 
      // boundary between the enriched zone and the rest of the domain.
      getfem::mesh_region &boundary = mesh.region(MORTAR_BOUNDARY_IN);
      unsigned int cpt = 0 ;
      for (dal::bv_visitor i(cvlist_in_area); !i.finished(); ++i) {
         if (mls.is_convex_cut(i)){
	    // Among the faces of the convex, we search if some are
	    // part of the boundary
	    cpt = 0 ;
	    for (unsigned j=0; j < mesh.structure_of_convex(i) ->nb_faces(); ++j) {
	        if (boundary.is_in(i,j))
		   cpt += 1;
	    }
	    if (cpt) {
               cout << "\n The convex number " << i << " is crossed by the crack :\n" ;
	       cout << "  it has : " << cpt << " face(s) among the boundary.\n \n " ;
	    }
	 }
      }
 
    }  // end of "enrichment_option = 3"
    break ;
  default : 
	GMM_ASSERT1(false, "Enrichment_option parameter is undefined");
	break ;  
	}
  mesh.write_to_file("toto.mesh");
  
  if (0) {  // printing the type of each dof
    unsigned Q = mf_u().get_qdim();
    for (unsigned d=0; d < mf_u().nb_dof(); d += Q) {
      printf("dof %4d @ %+6.2f:%+6.2f: ", d, 
             mf_u().point_of_dof(d)[0], mf_u().point_of_dof(d)[1]);
      
      
      const getfem::mesh::ind_cv_ct cvs = mf_u().convex_to_dof(d);
      for (unsigned i=0; i < cvs.size(); ++i) {
        unsigned cv = cvs[i];
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
          printf(" %3d:%.16s", cv, name_of_dof(pf->dof_types().at(ld)).c_str());
        }
      }
      printf("\n");
    }
  }
  
  cout << "Number of dof for u: " << mf_u().nb_dof() << endl;

  // Bilaplacian brick.
  getfem::mdbrick_bilaplacian<> BIL(mim, mf_u());
  BIL.D().set(D);
  if (KL) { BIL.set_to_KL(); BIL.nu().set(nu); }
  
  // Volumic source term brick.
  size_type N = mesh.dim();
  plain_vector F(nb_dof_rhs);
  getfem::mdbrick_source_term<> VOL_F(BIL, mf_rhs, F);  
  
  // Defining the prescribed momentum.
  gmm::resize(F, nb_dof_rhs*N*N);
  getfem::mdbrick_normal_derivative_source_term<>
  MOMENTUM(VOL_F, mf_rhs, F, MOMENTUM_BOUNDARY_NUM);
  
  // Defining the Neumann condition right hand side.
  plain_vector HH(nb_dof_rhs*N*N);
  gmm::resize(F, nb_dof_rhs*N);
  getfem::interpolation_function(mf_rhs, F, sol_ff, FORCE_BOUNDARY_NUM);

  // Neumann condition brick.
  getfem::mdbrick_neumann_KL_term<> NEUMANN(MOMENTUM, mf_rhs, HH, F, FORCE_BOUNDARY_NUM);
 

  // Normal derivative Dirichlet condition brick
  gmm::resize(F, nb_dof_rhs*N);
  gmm::clear(F);
  getfem::mdbrick_normal_derivative_Dirichlet<> 
	NDER_DIRICHLET(NEUMANN, CLAMPED_BOUNDARY_NUM, mf_mult_d);
  NDER_DIRICHLET.set_constraints_type(dirichlet_version);
  NDER_DIRICHLET.rhs().set(mf_rhs, F);

  getfem::mdbrick_constraint<> CLOSING_MODEL(NDER_DIRICHLET, 0);
  
  
  CLOSING_MODEL.set_constraints_type(getfem::constraints_type(closing_version));  
    if (closing_version == getfem::PENALIZED_CONSTRAINTS)
      CLOSING_MODEL.set_penalization_parameter(PARAM.real_value("EPS_CLOSING_PENAL")) ;
  
  cout << "building constraints matrix : \n" ;
  base_matrix C(1, mf_u().nb_dof() );
  int nb_line_C = 0 ;
  
  unsigned q = mf_u().get_qdim();
  unsigned cv, ld = unsigned(-1);
  dal::bit_vector dofs_up, dofs_down ;
  for (size_type i=0; i< mf_u().nb_dof(); ++i){
      cv = mf_u().first_convex_of_dof(i) ;
      getfem::pfem pf = mf_u().fem_of_element(cv);
      ld = unsigned(-1) ;
      for (unsigned dd = 0; dd < mf_u().nb_dof_of_element(cv); dd += q) {
	     if (mf_u().ind_dof_of_element(cv)[dd] == i) {
	        ld = dd/q; break;
	     }
         }  
//       if ( ( gmm::abs(mf_u().point_of_dof(i)[1]) > (0.5-1e-6) )
//       //  && ( gmm::abs(mf_u().point_of_dof(i)[0]) > (0.5-1e-6) )
//         && ( is_lagrange_dof_type(pf->dof_types().at(ld)) ) ){
// 	    nb_line_C += 1 ;
// 	    gmm::resize(C, nb_line_C, mf_u().nb_dof()) ;
// 	    C(nb_line_C - 1, i) = 1 ;     
//       }

//       if ( ( gmm::abs(mf_u().point_of_dof(i)[1]) > (0.5-1e-6) ) 
//         && ( gmm::abs(mf_u().point_of_dof(i)[0]) > (0.5-1e-6) ) 
//         && ( is_dx_dof_type(pf->dof_types().at(ld)) ) ){
// 	    cout << name_of_dof(pf->dof_types().at(ld)) << "\n" ;
// 	    nb_line_C += 1 ;
// 	    gmm::resize(C, nb_line_C, mf_u().nb_dof()) ;
// 	    C(nb_line_C - 1, i) = 1 ;
//       }
      if ( ( gmm::abs(mf_u().point_of_dof(i)[1]) < 1e-6 )
        && ( mf_u().point_of_dof(i)[0] > (0.5-1e-6) )
        && ( is_lagrange_dof_type(pf->dof_types().at(ld)) ) ){
	    cout << "valeur \n" ;
	    nb_line_C += 1 ;
	    gmm::resize(C, nb_line_C, mf_u().nb_dof()) ;
	    C(nb_line_C - 1, i) = 1 ;
      }
      if ( ( gmm::abs(mf_u().point_of_dof(i)[1]) < 1e-6 ) 
        && ( mf_u().point_of_dof(i)[0] > (0.5-1e-6) ) 
        && ( is_dy_dof_type(pf->dof_types().at(ld)) ) ){
	    cout << "dérivée par rapport à y \n" ;
	    cout << name_of_dof(pf->dof_types().at(ld)) << "\n" ;
	    nb_line_C += 1 ;
	    gmm::resize(C, nb_line_C, mf_u().nb_dof()) ;
	    C(nb_line_C - 1, i) = 1 ;
      }
      
      if ( mf_u().point_of_dof(i)[1] > (0.5 - 1e-6) ){
         dofs_up.add(i) ;      }
      if ( mf_u().point_of_dof(i)[1] < ( -0.5 + 1e-6) ){
         dofs_down.add(i) ;    }
  }

//   size_type cpt = 0 ;
//   for (unsigned i = 0 ; i < mf_u().nb_dof() ; ++i){  
//      if (dofs_up[i]){
//         if (cpt){
//         nb_line_C += 1 ;
// 	gmm::resize(C, nb_line_C, mf_u().nb_dof()) ;
// 	C(nb_line_C - 1, i  ) = 1 ;
// 	C(nb_line_C - 1, i-1) = -1 ;
//         }
//      else cpt +=1 ;
//      }
//   }


  cout << "nbre de lignes de C : " << gmm::mat_nrows(C) << "\n" ;
  plain_vector R(gmm::mat_nrows(C)) ;
  CLOSING_MODEL.set_constraints(C, R);
  
  getfem::mdbrick_constraint<> &mortar = 
      *(new getfem::mdbrick_constraint<>(CLOSING_MODEL, 0));    
  if (enrichment_option == 3 ) {
     /* add a constraint brick for the mortar junction between
       the enriched area and the rest of the mesh */
       int mult_with_H = PARAM.int_value("MULT_WITH_H") ;
       mortar_type = PARAM.int_value("MORTAR_TYPE") ;
       getfem::mesh_fem &mf_mortar = (mult_with_H == 1) ? mfls_mortar : mf_pre_mortar;
       getfem::mesh_fem &mf_mortar_deriv = (mult_with_H == 1) ? mfls_mortar_deriv : mf_pre_mortar_deriv;     
 
    mortar.set_constraints_type(getfem::constraints_type(mortar_version));  
    if (mortar_version == getfem::PENALIZED_CONSTRAINTS)
      mortar.set_penalization_parameter(PARAM.real_value("EPS_MORTAR_PENAL")) ;

    sparse_matrix H0(mf_mortar.nb_dof(), mf_u().nb_dof()), 
       H(mf_mortar_deriv.nb_dof(), mf_u().nb_dof());
    gmm::resize(R, mf_mortar.nb_dof());
  /* other version of the integral matching :
       *     \int_Gamma        (u-v) \lambda  = 0, for all \lambda in \Lambda
       *     \int_Gamma \partial_n (u-v)\mu  = 0, for all \mu in M
      */
      
      /* build the list of dof for the "(u-v) lambda" condition  */  
      dal::bit_vector bv_mortar;
      dal::bit_vector bv_deriv;
      sparse_matrix MM(mf_mortar.nb_dof(), mf_mortar.nb_dof());
      sparse_matrix MD(mf_mortar_deriv.nb_dof(), mf_mortar_deriv.nb_dof());
      std::vector<size_type> ind_mortar;
      std::vector<size_type> ind_deriv;
      getfem::asm_mass_matrix(MM, mim, mf_mortar, MORTAR_BOUNDARY_OUT);
      getfem::asm_mass_matrix(MD, mim, mf_mortar_deriv, MORTAR_BOUNDARY_OUT);

      for (size_type i=0; i < mf_mortar.nb_dof(); ++i)
	if (MM(i,i) > 1e-15) { 
	  // cout << "dof " << i << " MM = " << MM(i,i) << endl;
	  bv_mortar.add(i); ind_mortar.push_back(i);
	}
      for (size_type i=0; i < mf_mortar_deriv.nb_dof(); ++i)
	if (MD(i,i) > 1e-15) { bv_deriv.add(i); ind_deriv.push_back(i); }
      
      // building matrices
      
      cout << "Handling mortar junction (" << ind_mortar.size() << 
	" dof for the lagrange multiplier of the displacement, " <<
	ind_deriv.size() << " dof for the lagrange multiplier of the derivative)\n";
      
      gmm::resize(H0, mf_mortar.nb_dof(), mf_u().nb_dof()) ;
      gmm::resize(H,  ind_mortar.size() + ind_deriv.size(), mf_u().nb_dof()) ; 
      
      // Defining sub_indexes of the matrices calculated with the 
      // complete set of dofs.
      gmm::sub_index sub_i(ind_mortar);
      gmm::sub_index sub_i1(ind_deriv);
      gmm::sub_interval sub_j(0, mf_u().nb_dof());
      
      gmm::sub_interval sub_val_H(0, ind_mortar.size()) ;
      gmm::sub_interval sub_deriv_H(ind_mortar.size(), ind_deriv.size()) ;
      
      cout << "sub_indexes built\n" ;
      /* build the mortar constraint matrix -- note that the integration
	 method is conformal to the crack
      */
      gmm::clear(H0);
      getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u(), MORTAR_BOUNDARY_OUT);
      gmm::copy(gmm::sub_matrix(H0, sub_i, sub_j), gmm::sub_matrix(H, sub_val_H, sub_j) );
    
      gmm::clear(H0);
      getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u(), MORTAR_BOUNDARY_IN);
      gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i, sub_j), -1), gmm::sub_matrix(H, sub_val_H, sub_j) );
      
      
      cout << "first contraint asm\n" ;
      gmm::resize(R, ind_deriv.size());
      gmm::clear(H0);
      gmm::resize(H0, mf_mortar_deriv.nb_dof(), mf_u().nb_dof() ) ;

      getfem::asm_normal_derivative_dirichlet_constraints
	(H0, R, mim, mf_u(),
	 mf_mortar_deriv, mf_pre_u, R, MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;

      gmm::add(gmm::sub_matrix(H0, sub_i1, sub_j),
	       gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
      
      cout << "first step \n" ;
      gmm::clear(H0);

      getfem::asm_normal_derivative_dirichlet_constraints
	(H0, R, mim, mf_u(),
	 mf_mortar_deriv, mf_pre_u, R, MORTAR_BOUNDARY_IN, 0, getfem::ASMDIR_BUILDH) ;
      gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i1, sub_j), +1.0),
	       gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
	       
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
	unsigned n = gmm::mat_nrows(H);
	gmm::resize(H, n+1, gmm::mat_ncols(H));
	H(n, d) = 1;
      }
    }  
    gmm::resize(R, gmm::mat_nrows(H)); 
    mortar.set_constraints(H,R);
    cout << "Total number of variables : " << mortar.nb_dof() << endl;
    getfem::standard_model_state MS(mortar);
    gmm::iteration iter(residual, 1, 40000); 
    getfem::standard_solve(MS, mortar, iter /* , p*/);  
    // getfem::useful_types<getfem::standard_model_state>::plsolver_type p;
    // p.reset(new getfem::linear_solver_cg_preconditioned_ildlt<sparse_matrix,plain_vector>);
  // Solution extraction
  gmm::resize(U, mf_u().nb_dof());
  gmm::copy(BIL.get_solution(MS), U); 
  }
  else{
    cout << "Total number of variables : " << CLOSING_MODEL.nb_dof() << endl;
    getfem::standard_model_state MS(CLOSING_MODEL);
    gmm::iteration iter(residual, 1, 40000);
    getfem::standard_solve(MS, CLOSING_MODEL, iter /* , p*/);  
    // getfem::useful_types<getfem::standard_model_state>::plsolver_type p;
    // p.reset(new getfem::linear_solver_cg_preconditioned_ildlt<sparse_matrix,plain_vector>);
    // Solution extraction
    gmm::resize(U, mf_u().nb_dof());
    gmm::copy(BIL.get_solution(MS), U); 
    
  }


  return true;  
}
