#include "crack_bilaplacian.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_fourth_order.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_superlu.h"




scalar_type crack_level(base_node P) {
return P[1];
}

scalar_type crack_tip_level(base_node P) {
return P[0];
}

// functions for assembling the constraints of the integral matching 
  
template<typename MAT, typename VECT1, typename VECT2>
void asm_normal_derivative_dirichlet_constraints_bis
(MAT &H, VECT1 &R, const getfem::mesh_im &mim, const getfem::mesh_fem &mf_u,
 const getfem::mesh_fem &mf_mult, const getfem::mesh_fem &mf_r,
 const VECT2 &r_data, const getfem::mesh_region &rg, bool R_must_be_derivated, 
 int version) {
  typedef typename gmm::linalg_traits<VECT1>::value_type value_type;
  typedef typename gmm::number_traits<value_type>::magnitude_type magn_type;
    
  rg.from_mesh(mim.linked_mesh()).error_if_not_faces();
  GMM_ASSERT1(mf_r.get_qdim() == 1, 
              "invalid data mesh fem (Qdim=1 required)");
  if (version & getfem::ASMDIR_BUILDH) {
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
  if (version & getfem::ASMDIR_BUILDR) {
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

void bilaplacian_crack_problem::set_matrix_mortar(sparse_matrix &H){

  int mult_with_H = int(PARAM.int_value("MULT_WITH_H")) ;
  mortar_type = int(PARAM.int_value("MORTAR_TYPE")) ;
       // MODIFICATION : trying to define 1-Dimensionnal 
       // multipliers on the boundary of enrichment zone.
//   cout << "Re-initialisation des mesh_fem mortar : " ;
//   mf_pre_mortar.set_finite_element(MORTAR_BOUNDARY_IN,
//              getfem::fem_descriptor(PARAM.string_value("MORTAR_FEM_TYPE")));
//   mf_pre_mortar_deriv.set_finite_element(MORTAR_BOUNDARY_IN,
//              getfem::fem_descriptor(PARAM.string_value("MORTAR_DERIV_FEM_TYPE")));
//   cout << " OK !! \n" ;
// 
//   mfls_mortar.adapt();
//   mfls_mortar_deriv.adapt();
       // END MODIFICATION
       getfem::mesh_fem &mf_mortar = (mult_with_H == 1) ? mfls_mortar : mf_pre_mortar;
       getfem::mesh_fem &mf_mortar_deriv = (mult_with_H == 1) ? mfls_mortar_deriv : mf_pre_mortar_deriv;

sparse_matrix H0(mf_mortar.nb_dof(), mf_u().nb_dof()) ;
    getfem::base_vector R(mf_mortar.nb_dof());

gmm::resize( H, mf_mortar_deriv.nb_dof(), mf_u().nb_dof());

 if (mortar_type == 1) {
      /* older version of integral matching (jan-feb 2007) :
       *
       *  \int_Gamma (u-v) \lambda + \partial_n ( u + v) \partial_n \lambda = 0, for all \lambda in \Lambda      
      */
 
      /* build the list of dof for the "(u-v) lambda" and for the  "\partial_n(u-v) \partial_n lambda" term in the mortar condition */  
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
      asm_normal_derivative_dirichlet_constraints_bis(MD, R, mim, mf_mortar,
            mf_mortar, mf_pre_mortar, R, MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;
      
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
      gmm::resize(H0, mf_mortar.nb_dof(), mf_u().nb_dof()) ;
      gmm::resize(H,  ind_union.size(),   mf_u().nb_dof()) ;
//       cout << "bv_mortar = " << bv_mortar << "\n";
//       cout << "bv_deriv = " << bv_deriv << "\n" ;
//       cout << "bv_union = " << bv_union << "\n" ;
      
      gmm::sub_index sub_i(ind_union);
      gmm::sub_index sub_i1(ind_mortar);
      gmm::sub_index sub_i2(ind_deriv);
      gmm::sub_interval sub_j(0, mf_u().nb_dof());
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
      getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u(), MORTAR_BOUNDARY_OUT);
      gmm::copy(gmm::sub_matrix(H0, sub_i1, sub_j), gmm::sub_matrix(H, sub_val_H, sub_j) );
      
      gmm::clear(H0);
      getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u(), MORTAR_BOUNDARY_IN);
      gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i1, sub_j), -1), gmm::sub_matrix(H, sub_val_H, sub_j));
      
      
      gmm::clear(H0) ;
      gmm::resize(R, mat_nrows(H));
      asm_normal_derivative_dirichlet_constraints_bis(H0, R, mim, mf_u(),
                                                      mf_mortar, mf_pre_u, R, MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;
      gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i2, sub_j), 1. ), gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
      
      gmm::clear(H0) ;
      asm_normal_derivative_dirichlet_constraints_bis(H0, R, mim, mf_u(),
                                                      mf_mortar, mf_pre_u, R, MORTAR_BOUNDARY_IN, 0, getfem::ASMDIR_BUILDH) ;
      gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i2, sub_j), -1.), gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
      
      // -----------------------------------
    }
    if (mortar_type == 2 || mortar_type == 3) {
      /* other version of the integral matching.
       * version 2 :
       *     \int_Gamma        (u-v) \lambda  = 0, for all \lambda in \Lambda
       *     \int_Gamma \nabla (u-v).\mu      = 0, for all \mu in M    (be carefull : \mu is vectorial.
       * version 3 : only second constraint is different.   
       *     \int_Gamma \partial_n (u-v)\mu  = 0, for all \mu in M
      */
      if (mortar_type == 2){
          mf_mortar_deriv.set_qdim(2) ; //  \mu is vectorial !
      }
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
	  bv_mortar.add(i); 
	  ind_mortar.push_back(i);
	}
      for (size_type i=0; i < mf_mortar_deriv.nb_dof(); ++i)
	if (MD(i,i) > 1e-15) { 
	  bv_deriv.add(i); 
	  ind_deriv.push_back(i); 
        }
      
      // building matrices
      
      cout << "Handling mortar junction (" << ind_mortar.size() << 
	" dof for the lagrange multiplier of the displacement, " <<
	ind_deriv.size() << " dof for the lagrange multiplier of the derivative)\n";
      gmm::clear(H0);
      gmm::clear(H);
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
      getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u(), MORTAR_BOUNDARY_OUT);
      gmm::copy(gmm::sub_matrix(H0, sub_i, sub_j), gmm::sub_matrix(H, sub_val_H, sub_j) );
    
      gmm::clear(H0);
      getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u(), MORTAR_BOUNDARY_IN);
      gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i, sub_j), -1), gmm::sub_matrix(H, sub_val_H, sub_j) );
      
      
      cout << "first contraint asm\n" ;
      gmm::resize(R, ind_deriv.size());
      gmm::clear(H0);
      gmm::resize(H0, mf_mortar_deriv.nb_dof(), mf_u().nb_dof() ) ;
      if (mortar_type == 2)
         asm_constraint_gradient_vectorial_mult
           (H0, mim, mf_u(), mf_mortar_deriv,
           MORTAR_BOUNDARY_OUT, getfem::ASMDIR_BUILDH) ;
      if (mortar_type == 3)
         getfem::asm_normal_derivative_dirichlet_constraints
	   (H0, R, mim, mf_u(), mf_mortar_deriv, mf_pre_u, R,
	    MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;

      gmm::add(gmm::sub_matrix(H0, sub_i1, sub_j),
	       gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
      
      cout << "first step \n" ;
      gmm::clear(H0);
      if (mortar_type == 2) {
         asm_constraint_gradient_vectorial_mult
	   (H0, mim, mf_u(), mf_mortar_deriv, 
	   MORTAR_BOUNDARY_IN, getfem::ASMDIR_BUILDH) ;
      gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i1, sub_j), - 1.0),
	       gmm::sub_matrix(H, sub_deriv_H, sub_j)) ; 
      }
      if (mortar_type == 3) {
         getfem::asm_normal_derivative_dirichlet_constraints
	   (H0, R, mim, mf_u(), mf_mortar_deriv, mf_pre_u, R,
	    MORTAR_BOUNDARY_IN, 0, getfem::ASMDIR_BUILDH) ;
      gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i1, sub_j), + 1.0),
	       gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
      }


    }

}

void bilaplacian_crack_problem::init_mixed_elements(void) {
  std::string MESH_FILE = PARAM.string_value("MESH_FILE");
  std::string TRI_MESH_TYPE = PARAM.string_value("TRI_MESH_TYPE","Mesh type ");
  std::string QUAD_MESH_TYPE = PARAM.string_value("QUAD_MESH_TYPE","Mesh type ");
  std::string TRI_FEM_TYPE  = PARAM.string_value("TRI_FEM_TYPE","FEM name");
  std::string QUAD_FEM_TYPE  = PARAM.string_value("QUAD_FEM_TYPE","FEM name");
  std::string TRI_INTEGRATION = PARAM.string_value("TRI_INTEGRATION",
                                               "Name of integration method");
  std::string QUAD_INTEGRATION = PARAM.string_value("QUAD_INTEGRATION",
                                               "Name of integration method");
  
  cout << "TRI_MESH_TYPE=" << TRI_MESH_TYPE << "\n";
  cout << "TRI_FEM_TYPE="  << TRI_FEM_TYPE << "\n";
  cout << "TRI_INTEGRATION=" << TRI_INTEGRATION << "\n";
  cout << "QUAD_MESH_TYPE=" << QUAD_MESH_TYPE << "\n";
  cout << "QUAD_FEM_TYPE="  << QUAD_FEM_TYPE << "\n";
  cout << "QUAD_INTEGRATION=" << QUAD_INTEGRATION << "\n";
  
  size_type /*N_tri = 2, N_quad = 2,*/ N = 2 ;
    
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
                                         "Name of simplex integration method");
  std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION");
  enrichment_option = int(PARAM.int_value("ENRICHMENT_OPTION",
					  "Enrichment option"));
    enr_area_radius = PARAM.real_value("RADIUS_ENR_AREA",
                                     "radius of the enrichment area");

sol_ref = PARAM.int_value("SOL_REF") ;
  
  bgeot::pgeometric_trans pgt_tri = 
      bgeot::geometric_trans_descriptor(TRI_MESH_TYPE);
  bgeot::pgeometric_trans pgt_quad = 
      bgeot::geometric_trans_descriptor(QUAD_MESH_TYPE);
    
    /* First step : build the mesh */

 if (!MESH_FILE.empty()) {
    mesh.read_from_file(MESH_FILE);
    cout << "WARNING : THE MESH FILES MUST DESIGN A MESH MADE OF TRIANGLES \n" ;
    base_small_vector tt(N); 
    tt[0] = PARAM.real_value("TRANSLAT_X") ;  
    tt[1] = PARAM.real_value("TRANSLAT_Y") ; 
    cout << "TRANSLAT_X = " << tt[0] << " ; TRANSLAT_Y = " << tt[1] << "\n" ;
    mesh.translation(tt); 
    if (mesh.points_of_convex(mesh.convex_index().first_true()).size() == 3){
       TRI_MESH_TYPE = bgeot::name_of_geometric_trans
          (mesh.trans_of_convex(mesh.convex_index().first_true()));
       cout << "TRI_MESH_TYPE=" << TRI_MESH_TYPE << "\n";
      }
    else if (mesh.points_of_convex(mesh.convex_index().first_true()).size() == 4){
       QUAD_MESH_TYPE = bgeot::name_of_geometric_trans
          (mesh.trans_of_convex(mesh.convex_index().first_true()));
       cout << "QUAD_MESH_TYPE=" << QUAD_MESH_TYPE << "\n";
      }
    else cout << "ERROR : elements are not triangles or quadrangles \n" ;
    pgt_tri = bgeot::geometric_trans_descriptor(TRI_MESH_TYPE);
    pgt_quad = bgeot::geometric_trans_descriptor(QUAD_MESH_TYPE);
    N = mesh.dim();
 } else {
    GMM_ASSERT1(N == 2, "For a plate problem, N should be 2");
    std::vector<size_type> nsubdiv(N);
    NX = PARAM.int_value("NX", "Number of space steps ") ;
    std::fill(nsubdiv.begin(),nsubdiv.end(), NX);
    if (PARAM.int_value("QUAD") == 0)
       getfem::regular_unit_mesh(mesh, nsubdiv, pgt_tri, PARAM.int_value("MESH_NOISED") != 0);
    else
       getfem::regular_unit_mesh(mesh, nsubdiv, pgt_quad, PARAM.int_value("MESH_NOISED") != 0);
    bgeot::base_matrix M(N,N);
    for (size_type i=0; i < N; ++i) {
      static const char *t[] = {"LX","LY","LZ"};
      M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
    }
    /* scale the unit mesh to [LX,LY,..] and incline it */
    mesh.transformation(M);
    base_small_vector tt(N); 
    tt[0] = - 0.5 + PARAM.real_value("TRANSLAT_X") ; 
    tt[1] = - 0.5 + PARAM.real_value("TRANSLAT_Y") ; 
    mesh.translation(tt); 
    
    /* MODIFICATION RELATIVE TO THE USUAL WAY OF PROGRAMMING :
        In the other xfem programs, the level-set is initialized 
        after the mesh_fems. However, to define the mesh_fem, here, 
        we need to select the cut quadrangles. So, it is necessary to
        switch the orders here :
		- first, the level-set is initialized
		- second, the mesh_fem are initialized too.
    */

   // Setting the level-set
  ls.reinit();
  scalar_type a = PARAM.real_value("CRACK_SEMI_LENGTH") ; 
  for (size_type d = 0; d < ls.get_mesh_fem().nb_dof(); ++d) {
    scalar_type x = ls.get_mesh_fem().point_of_dof(d)[0];
    scalar_type y = ls.get_mesh_fem().point_of_dof(d)[1];
    if (sol_ref == 0){
       ls.values(0)[d] = y  ; // + 1/4.*(x + .25);
       ls.values(1)[d] = x;}
    if (sol_ref == 1){
     ls.values(0)[d] = y  ;
     ls.values(1)[d] = x - a ; //x * x - a * a ;
    }
    if (sol_ref == 2){
     ls.values(0)[d] = y - x ;
     ls.values(1)[d] = gmm::abs(x + y) - a ; //x * x - a * a ;
    }
  }
  //ls.simplify(0.5);
  ls.touch(); 
  mls.add_level_set(ls); 
  mls.adapt();

    if (PARAM.int_value("MOVE_NODES")){
       cout << "d�placement des noeuds \n" ;
       //size_type nb_x_pos, nb_y_pos = 0 ;
       scalar_type seuil_select = PARAM.real_value("SEUIL_SELECT") ;
       //scalar_type seuil_move = PARAM.real_value("SEUIL_MOVE") ;


// for(dal::bv_visitor i(mesh.convex_index()) ; !i.finished() ; ++i){
//    nb_x_pos = 0 ; // nombre de noeuds d'abscisse positive
//    nb_y_pos = 0 ; // nombre de noeuds d'ordonn�e positive
//    for (int j=0; j<4 ; ++j){
//        if (mesh.points_of_convex(i)[j][0] > 0.) 
//            nb_x_pos += 1 ;
//        if (mesh.points_of_convex(i)[j][1] > 0.) 
//            nb_y_pos += 1 ;
//    }
//  
//    if (nb_x_pos < 4){ // juste pour �viter de traiter des �l�ments inutiles
//        if ( nb_y_pos == 1){
// 	  for (int j=0; j<4 ; ++j){
// 	      if ( (mesh.points_of_convex(i)[j][1] > 0.)
//                &&  (mesh.points_of_convex(i)[j][1] < seuil_select )  // 1./5./NX
// 	       ){
// 	         bgeot::base_node Q = mesh.points_of_convex(i)[j] ;
// 	         for (dal::bv_visitor ip(mesh.points().index()); !ip.finished(); ++ip) {
//                     bgeot::base_node& P = mesh.points()[ip];
// 	            if( gmm::vect_dist2(P, Q) < 1e-8){
// 		      cout << "d�plac� de (" << P[0] << " ; " << P[1] << ") � : " ;
// 		      //P[1] = 1./2./NX ;
// 		      //P[1] = 1./2./NX ;
// 		      P[1] = seuil_move ;
// 		      cout << P[1] << "\n" ;
// 	            }  
// 	        }
//               }
// 	  }
//        } 
//     if ( nb_y_pos == 3){
// 	  for (int j=0; j<4 ; ++j){
// 	      if ( (mesh.points_of_convex(i)[j][1] < 0.)
//                &&  (mesh.points_of_convex(i)[j][1] > -1.*seuil_select ) // -1./5./NX 
// 	       ){
// 	         bgeot::base_node Q = mesh.points_of_convex(i)[j] ;
// 	         for (dal::bv_visitor ip(mesh.points().index()); !ip.finished(); ++ip) {
//                     bgeot::base_node& P = mesh.points()[ip];
// 	            if( gmm::vect_dist2(P, Q) < 1e-8){
// 		      cout << "d�plac� de (" << P[0] << " ; " << P[1] << ") � : " ;
// 		      //P[1] = - 1./2./NX ;
// 		      //P[1] = -1./2./NX ;
// 		      P[1] = -1. *seuil_move ;
// 		      cout << P[1] << "\n" ;
// 	            }  
// 	        }
//               }
// 	  }
//        } 
// 
//        }
//    }
    
  
   for (dal::bv_visitor ip(mesh.points().index()); !ip.finished(); ++ip) {
                    bgeot::base_node& P = mesh.points()[ip];
	            if( gmm::abs(P[1]) < seuil_select){
		      cout << "d�plac� de (" << P[0] << " ; " << P[1] << ") � : " ;
		      P[1] = 0. ;
		      cout << P[1] << "\n" ;
	            }
    } 
  } // end of moving nodes		     
 } // end of else
 
    /* Separation of the crossed quadrangles in two triangles */
    

    if (PARAM.real_value("QUAD") == 1 ) {
       cout << "transforming quad in two triangles\n";
       // Selecting crossed quadrangles :
       dal::bit_vector quad_crossed ;
       for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i) {
           if (mls.is_convex_cut(i) && mesh.points_of_convex(i).size() == 4){
	      quad_crossed.add(i) ;
           }
       }
       cout << "quadrangles travers�s par la fissure : \n" << quad_crossed << "\n" ;
       // Replace quadrangles by two triangles
       for (dal::bv_visitor i(quad_crossed) ; !i.finished() ; ++i){
              //cout << "Convex " << i << "to split ; values of the level_set : " << values_ls << "\n" ;
	      // First : save the global indexes of the 4 nodes
              std::vector<size_type> ind_sommets_cvx(4) ; 
              for (unsigned j = 0 ; j < (mesh.ind_points_of_convex(i)).size() ; ++j)
                  ind_sommets_cvx[j] = mesh.ind_points_of_convex(i)[j] ; 
              std::vector<bgeot::base_node> P(4) ;
	      // Second : save the coordinates, in order to calculate the length of the two diagonals
              for (unsigned node = 0 ; node < (mesh.points_of_convex(i)).size() ; ++node)
                  P[node] = mesh.points_of_convex(i)[node] ;
              mesh.sup_convex(i) ; // Delete convex now, we won't need it beyond this line.
              scalar_type dist1 = 0., dist2 = 0. ;
              for(int k = 0 ; k < 2 ; k++){ 
                  dist1 += (P[0][k] - P[3][k]) * (P[0][k] - P[3][k]) ;
                  dist2 += (P[2][k] - P[1][k]) * (P[2][k] - P[1][k]) ; 
              }
              // Cut the quadrangle along the shortest diagonal
	      std::vector<bgeot::size_type> ind(3) ;
              if (dist1 > dist2){ //then cut along the [s1 s2] diagonal
                 ind[0] = ind_sommets_cvx[0] ;
                 ind[1] = ind_sommets_cvx[1] ;
                 ind[2] = ind_sommets_cvx[2] ;
              }
              else{  // cut along the [s0 s3] diagonal
                 ind[0] = ind_sommets_cvx[0] ;
                 ind[1] = ind_sommets_cvx[2] ;
                 ind[2] = ind_sommets_cvx[3] ;
              }
              mesh.add_convex(pgt_tri, ind.begin());  // add first triangle
              if ( dist1 > dist2)
                 ind[0] = ind_sommets_cvx[3] ;
              else
                 ind[1] = ind_sommets_cvx[1] ;
              mesh.add_convex(pgt_tri, ind.begin());  // add second triangle 
       }
       // The procedure used beyond was made possible to use thanks to the
       // structure of the containers returned : the elements are ordered 
       // in accordance with the numerotation of the vertex and faces in a quadrangle:
       //         2____3        _2__
       //         |    |     1 |    |  0
       //         |____|       |____|
       //         0    1         3     
    }
    // Holding those informations is mandatory for the following
    dal::bit_vector quad_among_cvx, tri_among_cvx ;
    for (dal::bv_visitor i(mesh.convex_index()) ; !i.finished() ; ++i){
       if (mesh.points_of_convex(i).size() == 3)
          tri_among_cvx.add(i) ;
       else if (mesh.points_of_convex(i).size() == 4)
          quad_among_cvx.add(i) ;
       else cout << "WARNING : an element has nor 3 or 4 nodes ! \n" ;
    }
    
    scalar_type quality = 1.0, avg_area = 0. , min_area = 1. , max_area = 0., area ;
    scalar_type radius, avg_radius = 0., min_radius = 1., max_radius = 0. ;
    size_type cpt = 0 ; 
    for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i){
    	quality = std::min(quality, mesh.convex_quality_estimate(i));
	area = mesh.convex_area_estimate(i) ;
	radius = mesh.convex_radius_estimate(i) ;
	avg_area += area ;
	avg_radius += radius ;
	min_radius = std::min(radius, min_radius) ;
	max_radius = std::max(radius, max_radius) ;	
	min_area = std::min(min_area, area) ; 
	max_area = std::max(max_area, area) ;
	cpt++ ;
    }
    avg_area /= scalar_type(cpt) ;
    avg_radius /= scalar_type(cpt) ; 
    cout << "quality of mesh : " << quality << endl;
    cout << "average radius : " << avg_radius << endl;
    cout << "radius min : " << min_radius << " ; radius max : " << max_radius << endl;
    cout << "average area : " << avg_area << endl ;
    cout << "area min : " << min_area << " ; area max : " << max_area << endl;        
    
   /* read the parameters   */
  epsilon = PARAM.real_value("EPSILON", "thickness") ;
  nu = PARAM.real_value("NU", "nu") ;
  D = PARAM.real_value("D", "Flexure modulus") ;
  int dv = int(PARAM.int_value("DIRICHLET_VERSION", "Dirichlet version"));
  int mv = int(PARAM.int_value("MORTAR_VERSION", "Mortar version"));
  int cv = int(PARAM.int_value("CLOSING_VERSION"));  
  dirichlet_version = getfem::constraints_type(dv);
  mortar_version = getfem::constraints_type(mv);
  closing_version = getfem::constraints_type(cv);
  datafilename=PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual=PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;
  KL = PARAM.int_value("KL", "Kirchhoff-Love model or not") != 0;
  D = PARAM.real_value("D", "Flexion modulus");
  if (KL) nu = PARAM.real_value("NU", "Poisson ratio");

  // Setting the integration methods
  getfem::pfem pf_u_tri = getfem::fem_descriptor(TRI_FEM_TYPE);
  getfem::pfem pf_u_quad = getfem::fem_descriptor(QUAD_FEM_TYPE);
  getfem::pintegration_method ppi_tri = 
    getfem::int_method_descriptor(TRI_INTEGRATION);
  getfem::pintegration_method ppi_quad = 
    getfem::int_method_descriptor(QUAD_INTEGRATION);
  getfem::pintegration_method sppi = 
    getfem::int_method_descriptor(SIMPLEX_INTEGRATION);
  getfem::pintegration_method sing_ppi = (SINGULAR_INTEGRATION.size() ?
		getfem::int_method_descriptor(SINGULAR_INTEGRATION) : 0);
  
  mim.set_integration_method(tri_among_cvx, ppi_tri);
  mim.set_integration_method(quad_among_cvx, ppi_quad);
  mim.set_simplex_im(sppi, sing_ppi);
  
    /* Setting the finite element on the mf_u */  
  mf_pre_u.set_finite_element(tri_among_cvx, pf_u_tri); 
  mf_pre_u.set_finite_element(quad_among_cvx, pf_u_quad);  
  getfem::pfem pf_partition_of_unity_tri = getfem::fem_descriptor(PARAM.string_value("TRI_PARTITION_OF_UNITY_FEM_TYPE")) ;
  getfem::pfem pf_partition_of_unity_quad = getfem::fem_descriptor(PARAM.string_value("QUAD_PARTITION_OF_UNITY_FEM_TYPE")) ; 
  mf_partition_of_unity.set_finite_element(tri_among_cvx, pf_partition_of_unity_tri);
  mf_partition_of_unity.set_finite_element(quad_among_cvx, pf_partition_of_unity_quad);

  mf_pre_mortar.set_finite_element(tri_among_cvx,
             getfem::fem_descriptor(PARAM.string_value("TRI_MORTAR_FEM_TYPE")));
  mf_pre_mortar_deriv.set_finite_element(tri_among_cvx,
             getfem::fem_descriptor(PARAM.string_value("TRI_MORTAR_DERIV_FEM_TYPE")));
  mf_pre_mortar.set_finite_element(quad_among_cvx,
             getfem::fem_descriptor(PARAM.string_value("QUAD_MORTAR_FEM_TYPE")));
  mf_pre_mortar_deriv.set_finite_element(quad_among_cvx,
             getfem::fem_descriptor(PARAM.string_value("QUAD_MORTAR_DERIV_FEM_TYPE")));  
             
  // set the mesh_fem of the multipliers (for the dirichlet condition)    
  std::string dirichlet_fem_name_tri = PARAM.string_value("TRI_DIRICHLET_FEM_TYPE");
  std::string dirichlet_fem_name_quad = PARAM.string_value("QUAD_DIRICHLET_FEM_TYPE");
  if (dirichlet_fem_name_tri.size() == 0 || dirichlet_fem_name_quad.size() == 0){
    mf_mult.set_finite_element(tri_among_cvx, pf_u_tri);
    mf_mult.set_finite_element(quad_among_cvx, pf_u_quad);    
    }
  else {
    cout << "TRI_DIRICHLET_FEM_TYPE="  << dirichlet_fem_name_tri << "\n";
    cout << "QUAD_DIRICHLET_FEM_TYPE="  << dirichlet_fem_name_quad << "\n";
    mf_mult.set_finite_element(tri_among_cvx, 
                               getfem::fem_descriptor(dirichlet_fem_name_tri));
    mf_mult.set_finite_element(quad_among_cvx, 
                               getfem::fem_descriptor(dirichlet_fem_name_quad));
  }
  std::string dirichlet_der_fem_name_tri
    = PARAM.string_value("TRI_DIRICHLET_DER_FEM_TYPE", "");
  std::string dirichlet_der_fem_name_quad
    = PARAM.string_value("QUAD_DIRICHLET_DER_FEM_TYPE", "");
  if (dirichlet_der_fem_name_tri.size() == 0 || dirichlet_der_fem_name_quad.size() == 0){
    mf_mult_d.set_finite_element(tri_among_cvx, pf_u_tri);
    mf_mult_d.set_finite_element(quad_among_cvx, pf_u_quad);
  }
  else {
    cout << "TRI_DIRICHLET_DER_FEM_TYPE="  << dirichlet_der_fem_name_tri << "\n";
    cout << "QUAD_DIRICHLET_DER_FEM_TYPE="  << dirichlet_der_fem_name_quad << "\n";
    mf_mult_d.set_finite_element(tri_among_cvx, 
                             getfem::fem_descriptor(dirichlet_der_fem_name_tri));
    mf_mult_d.set_finite_element(quad_among_cvx, 
                             getfem::fem_descriptor(dirichlet_der_fem_name_quad));
  }

  /* set the finite element on mf_rhs (same as mf_u if DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name_tri = PARAM.string_value("TRI_DATA_FEM_TYPE");
  std::string data_fem_name_quad = PARAM.string_value("QUAD_DATA_FEM_TYPE");
  if (data_fem_name_tri.size() == 0 || data_fem_name_quad.size() == 0) {
    GMM_ASSERT1(pf_u_tri->is_lagrange(), "You are using a non-lagrange FEM. "
                << "In that case you need to set "
                << "DATA_FEM_TYPE in the .param file");
    mf_rhs.set_finite_element(tri_among_cvx, pf_u_tri);
    mf_rhs.set_finite_element(quad_among_cvx, pf_u_quad);
  } else {
    mf_rhs.set_finite_element(tri_among_cvx, 
                              getfem::fem_descriptor(data_fem_name_tri));
    mf_rhs.set_finite_element(quad_among_cvx, 
                              getfem::fem_descriptor(data_fem_name_quad));
  }

  /* set boundary conditions
   * (Neuman on the upper face, Dirichlet elsewhere) */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  if (sol_ref == 0){
     for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
         mesh.region(SIMPLE_SUPPORT_BOUNDARY_NUM).add(i.cv(), i.f());
         mesh.region(CLAMPED_BOUNDARY_NUM).add(i.cv(), i.f()); 
     }
  }

  if (sol_ref == 1 || sol_ref == 2){
     for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
        base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
        un /= gmm::vect_norm2(un);
        //if ( (un[0] <= -0.9) || (un[0] >= 0.9) ) {
	if  (gmm::abs(un[0]) >= 0.999)  {
	   mesh.region(CLAMPED_BOUNDARY_NUM).add(i.cv(), i.f());
           //mesh.region(MOMENTUM_BOUNDARY_NUM).add(i.cv(), i.f());
	} else
           mesh.region(FORCE_BOUNDARY_NUM).add(i.cv(), i.f());
           mesh.region(MOMENTUM_BOUNDARY_NUM).add(i.cv(), i.f()); 
     }
  }

  // Updating things concerning the level-set
  scalar_type a = PARAM.real_value("CRACK_SEMI_LENGTH") ; 
  for (size_type d = 0; d < ls.get_mesh_fem().nb_dof(); ++d) {
    scalar_type x = ls.get_mesh_fem().point_of_dof(d)[0];
    scalar_type y = ls.get_mesh_fem().point_of_dof(d)[1];
    if (sol_ref == 0){
       ls.values(0)[d] = y  ; // + 1/4.*(x + .25);
       ls.values(1)[d] = x;}
    if (sol_ref == 1){
     ls.values(0)[d] = y  ;
     ls.values(1)[d] = x - a ; //x * x - a * a ;
    }
    if (sol_ref == 2){
     ls.values(0)[d] = y - x ;
     ls.values(1)[d] = gmm::abs(x + y) - a ; //x * x - a * a ;
    }
  }
  //ls.simplify(0.5);
  ls.touch(); 
  mls.add_level_set(ls); 
  mls.adapt();
  exact_sol.init(ls);
  mim.adapt();
  mfls_u.adapt();
  mfls_mortar.adapt();
  mfls_mortar_deriv.adapt();
  cout << "mfls_u.nb_dof()=" << mfls_u.nb_dof() << "\n";
}


