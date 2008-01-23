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
  
  size_type N_tri = 2, N_quad = 2, N = 2 ;
    
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
                                         "Name of simplex integration method");
  std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION");
  enrichment_option = PARAM.int_value("ENRICHMENT_OPTION",
                                      "Enrichment option");
    enr_area_radius = PARAM.real_value("RADIUS_ENR_AREA",
                                     "radius of the enrichment area");
  
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
    
    // Definition of the level-set
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
    if (PARAM.int_value("SOL_REF") == 2){
     ls.values(0)[d] = y - x ;
     ls.values(1)[d] = gmm::abs(x + y) - a ; //x * x - a * a ;
    }
  }
  //ls.simplify(0.5);
  ls.touch(); 
  mls.add_level_set(ls); 
  mls.adapt();
    
    if (PARAM.int_value("MOVE_NODES")){
       cout << "déplacement des noeuds \n" ;
       size_type nb_x_pos, nb_y_pos = 0 ;
       scalar_type seuil_select = PARAM.real_value("SEUIL_SELECT") ;
       scalar_type seuil_move = PARAM.real_value("SEUIL_MOVE") ;
    
    
// for(dal::bv_visitor i(mesh.convex_index()) ; !i.finished() ; ++i){
//    nb_x_pos = 0 ; // nombre de noeuds d'abscisse positive
//    nb_y_pos = 0 ; // nombre de noeuds d'ordonnée positive
//    for (int j=0; j<4 ; ++j){
//        if (mesh.points_of_convex(i)[j][0] > 0.) 
//            nb_x_pos += 1 ;
//        if (mesh.points_of_convex(i)[j][1] > 0.) 
//            nb_y_pos += 1 ;
//    }
//  
//    if (nb_x_pos < 4){ // juste pour éviter de traiter des éléments inutiles
//        if ( nb_y_pos == 1){
// 	  for (int j=0; j<4 ; ++j){
// 	      if ( (mesh.points_of_convex(i)[j][1] > 0.)
//                &&  (mesh.points_of_convex(i)[j][1] < seuil_select )  // 1./5./NX
// 	       ){
// 	         bgeot::base_node Q = mesh.points_of_convex(i)[j] ;
// 	         for (dal::bv_visitor ip(mesh.points().index()); !ip.finished(); ++ip) {
//                     bgeot::base_node& P = mesh.points()[ip];
// 	            if( gmm::vect_dist2(P, Q) < 1e-8){
// 		      cout << "déplacé de (" << P[0] << " ; " << P[1] << ") à : " ;
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
// 		      cout << "déplacé de (" << P[0] << " ; " << P[1] << ") à : " ;
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
		      cout << "déplacé de (" << P[0] << " ; " << P[1] << ") à : " ;
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
       cout << "quadrangles traversés par la fissure : \n" << quad_crossed << "\n" ;
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
    avg_area /= cpt ;
    avg_radius /= cpt ; 
    cout << "quality of mesh : " << quality << endl;
    cout << "average radius : " << avg_radius << endl;
    cout << "radius min : " << min_radius << " ; radius max : " << max_radius << endl;
    cout << "average area : " << avg_area << endl ;
    cout << "area min : " << min_area << " ; area max : " << max_area << endl;        
    
   /* read the parameters   */
  epsilon = PARAM.real_value("EPSILON", "thickness") ;
  nu = PARAM.real_value("NU", "nu") ;
  D = PARAM.real_value("D", "Flexure modulus") ;
  int dv = PARAM.int_value("DIRICHLET_VERSION", "Dirichlet version");
  int mv = PARAM.int_value("MORTAR_VERSION", "Mortar version");
  int cv = PARAM.int_value("CLOSING_VERSION");  
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
  if (PARAM.int_value("SOL_REF") == 0){
     for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
         mesh.region(SIMPLE_SUPPORT_BOUNDARY_NUM).add(i.cv(), i.f());
         mesh.region(CLAMPED_BOUNDARY_NUM).add(i.cv(), i.f()); 
     }
  }
  
  if (PARAM.int_value("SOL_REF") == 1 || PARAM.int_value("SOL_REF") == 2){
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
  exact_sol.init(ls);
  mim.adapt();
  mfls_u.adapt();
  mfls_mortar.adapt();
  mfls_mortar_deriv.adapt();
  cout << "mfls_u.nb_dof()=" << mfls_u.nb_dof() << "\n";
}


