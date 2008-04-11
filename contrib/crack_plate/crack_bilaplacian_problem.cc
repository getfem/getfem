#include "crack_bilaplacian.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_fourth_order.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_superlu.h"

/******** Exact Solution *******************************/

scalar_type D  = 1.  ;
scalar_type nu = 0.3 ;
scalar_type AAA = 0.1 ; // 1.0 ;
scalar_type BB = AAA * (3. * nu + 5.)/ (3. * (nu - 1.))   ;  // (-3.0+nu*nu-2.0*nu)/(nu*nu-2.0*nu+5.0);
scalar_type DD = 0.0 ; // 0.1 ;
scalar_type CC = DD * (nu + 7.)/ (3. * (nu - 1.))   ;   //  (-8.0*nu+3.0*AAA*nu*nu-6.0*nu*AAA+15.0*AAA)/(nu*nu-2.0*nu+5.0);
 

scalar_type sol_u(const base_node &x){
 scalar_type r = sqrt( x[0] * x[0] + x[1] * x[1] ) ;
 //scalar_type theta = 2. * atan( x[1] / ( x[0] + r ) ) ;
 scalar_type theta = atan2(x[1], x[0]);
 //return sqrt(r*r*r) * (sin(3.0/2.0*theta)+BB*sin(theta/2.0)+AAA*cos(3.0/2.0*theta)+CC*cos(theta/2.0));
 return sqrt(r*r*r)*(AAA*sin(theta/2.0)+BB*sin(3.0/2.0*theta)+CC*cos(3.0/2.0*theta)+DD*cos(theta/2.0));
}
 
scalar_type sol_lapl_u(const base_node &x) {
 scalar_type r = sqrt( x[0] * x[0] + x[1] * x[1] ) ;
 scalar_type theta = atan2(x[1], x[0]);
 return 2.0*(AAA*sin(theta/2.0)+DD*cos(theta/2.0))/sqrt(r);
 /* return 9.0/4.0/sqrt(r)*(sin(3.0/2.0*theta)+BB*sin(theta/2.0)+AAA*cos(3.0/2.0*theta)+CC*cos(theta/2.0))+1/sqrt(r)*(-9.0/4.0*sin(3.0/2.0*theta)-BB*sin(theta/
    2.0)/4.0-9.0/4.0*AAA*cos(3.0/2.0*theta)-CC*cos(theta/2.0)/4.0); */ }

scalar_type sol_f(const base_node &)
{ return 0. ; }

base_small_vector sol_du(const base_node &x) {
 base_small_vector res(x.size());
 scalar_type r = sqrt( x[0] * x[0] + x[1] * x[1] ) ;
 scalar_type theta = atan2(x[1], x[0]);
res[0] = 3.0/2.0*sqrt(r)*(BB*sin(3.0/2.0*theta)+AAA*sin(theta/2.0)+CC*cos(3.0/2.0*theta)
+DD*cos(theta/2.0))*cos(theta)-sqrt(r)*(3.0/2.0*BB*cos(3.0/2.0*theta)+
AAA*cos(theta/2.0)/2.0-3.0/2.0*CC*sin(3.0/2.0*theta)-DD*sin(theta/2.0)/2.0)*sin(theta );

res[1] = 3.0/2.0*sqrt(r)*(BB*sin(3.0/2.0*theta)+AAA*sin(theta/2.0)+CC*cos(3.0/2.0*theta)
+DD*cos(theta/2.0))*sin(theta)+sqrt(r)*(3.0/2.0*BB*cos(3.0/2.0*theta)+
AAA*cos(theta/2.0)/2.0-3.0/2.0*CC*sin(3.0/2.0*theta)-DD*sin(theta/2.0)/2.0)*cos(theta);

/*
res[0] =  3.0/2.0*sqrt(r)*(sin(3.0/2.0*theta)+BB*sin(theta/2.0)+AAA*cos(3.0/2.0*theta)+CC*cos(theta/2.0))*cos(theta)-sqrt(r)*(3.0/2.0*cos(3.0/2.0*theta)+BB*
cos(theta/2.0)/2.0-3.0/2.0*AAA*sin(3.0/2.0*theta)-CC*sin(theta/2.0)/2.0)*sin(
theta);
 
res[1] = 3.0/2.0*sqrt(r)*(sin(3.0/2.0*theta)+BB*sin(theta/2.0)+AAA*cos(3.0/2.0*theta)+CC*cos(theta/2.0))*sin(theta)+sqrt(r)*(3.0/2.0*cos(3.0/2.0*theta)+BB*
cos(theta/2.0)/2.0-3.0/2.0*AAA*sin(3.0/2.0*theta)-CC*sin(theta/2.0)/2.0)*cos(
theta);
*/
  return res;
}

base_small_vector neumann_val(const base_node &x)
{ base_small_vector res(x.size());
  res[0] = 0. ;
  res[1] = 0. ;
return res ; }

// base_matrix sol_hessian(const base_node &x) {
//   base_matrix m(x.size(), x.size());
//   // remplir la matrice hessienne de la solution exacte 
//   return m;
// }

// base_matrix sol_mtensor(const base_node &x) {
//   // moment de flexion de la solution exacte 
//   base_matrix m = sol_hessian(x), mm(x.size(), x.size());
//   scalar_type l = sol_lapl_u(x);
//   for (size_type i = 0; i < x.size(); ++i) mm(i,i) = l * nu;
//   gmm::scale(m, (1-nu));
//   gmm::add(mm, m);
//   gmm::scale(m, -D);
//   return m;
// }

// base_small_vector sol_bf(const base_node &x)
// { return -D * neumann_val(x); }


void exact_solution::init(getfem::level_set &ls) {
  std::vector<getfem::pglobal_function> cfun(4) ;
  for (unsigned j=0; j < 4; ++j)
    cfun[j] = bilaplacian_crack_singular(j, ls, nu) ;
  mf.set_functions(cfun);
  U.resize(4); assert(mf.nb_dof() == 4);
  // scalar_type A1 = 1., nu = 0.3 ;
  // scalar_type b1_ = 3. + (A2 / A1) * (24. * nu) / (3. * nu * nu - 6. * nu + 5. ) ; 
  U[0] = AAA ;
  U[1] = BB ;
  U[2] = CC ;
  U[3] = DD ;
}

scalar_type eval_fem_gradient_with_finite_differences(getfem::pfem pf, 
					       const base_vector &coeff,
					       size_type cv,
					       bgeot::pgeometric_trans pgt, 
					       bgeot::geotrans_inv_convex &gic,
					       const base_matrix &G, 
					       base_node X0, 
					       scalar_type h, unsigned dg) {
  X0[dg] -= h/2;
  base_node X0ref; gic.invert(X0, X0ref);
  getfem::fem_interpolation_context c(pgt, pf, X0ref, G, cv);

  base_vector val0(1);
  pf->interpolation(c, coeff, val0, 1);

  base_node X1(X0), X1ref; X1[dg] += h;
  gic.invert(X1, X1ref);
  c.set_xref(X1ref);

  base_vector val1(1);
  pf->interpolation(c, coeff, val1, 1);

  return (val1[0] - val0[0])/h;
}

scalar_type eval_fem_hessian_with_finite_differences(getfem::pfem pf, 
					      const base_vector &coeff,
					      size_type cv, 
					      bgeot::pgeometric_trans pgt, 
					      bgeot::geotrans_inv_convex &gic,
					      const base_matrix &G, 
					      base_node X0, 
					      scalar_type h, 
					      unsigned dg, unsigned dh) {
  X0[dh] -= h/2;
  scalar_type Gr0 = 
    eval_fem_gradient_with_finite_differences(pf, coeff, cv, pgt, gic, G, X0, h, dg);
  base_node X1(X0);
  X1[dh] += h;
  scalar_type Gr1 = 
    eval_fem_gradient_with_finite_differences(pf, coeff, cv, pgt, gic, G, X1, h, dg);
  return (Gr1 - Gr0)/h;
}

void validate_fem_derivatives(getfem::pfem pf, unsigned cv, 
			      bgeot::pgeometric_trans pgt, const base_matrix &G) {
  unsigned N = gmm::mat_nrows(G);
  scalar_type h = 1e-5;

  std::vector<base_node> pts(gmm::mat_ncols(G));
  for (unsigned j=0; j < pts.size(); ++j) {
    pts[j].resize(N); gmm::copy(gmm::mat_col(G, j), pts[j]);
  }
  cout << "validate_fem_derivatives: pf = " << &(*pf) << ", nbdof = "<< pf->nb_dof(cv) << ", cv = " << cv << " (~ at " << gmm::mean_value(pts) << ")\n";
  bgeot::geotrans_inv_convex gic(pts, pgt);

  //cout << "pts = " << pts << "\n";
  
  for (unsigned idof = 0; idof < pf->nb_dof(cv); ++idof) {    
    /* choose a random point in the convex */
    base_node X0(N), X0ref;
    base_node w(pgt->nb_points());
    do {
      for (unsigned i=0; i < w.size(); ++i) w[i] = 0.1 + 0.8*gmm::random(); 
      gmm::scale(w, 1/gmm::vect_norm1(w));
      gmm::mult(G, w, X0);

      //cout << "w = " << w << "\n";
      
      gic.invert(X0, X0ref);
      
      // avoid discontinuity lines in the HCT composite element..
      if (gmm::abs(X0ref[0] + X0ref[1] - 1) > 1e-2 &&
	  gmm::abs(X0ref[0] - X0ref[1]) > 1e-2 && 
	  gmm::abs(X0[0]) > 1e-3 && gmm::abs(X0[1])> 1e-3) break;
    } while (1);
    //cout << "testing X0 = " << X0 << " (X0ref=" << X0ref << ")\n";


    base_vector coeff(pf->nb_dof(cv)); coeff[idof] = 1;
    base_matrix grad(1,N), grad_fd(1,N);
    base_matrix hess(1,N*N), hess_fd(1,N*N);

    getfem::fem_interpolation_context c(pgt, pf, X0ref, G, cv);
    pf->interpolation_grad(c, coeff, grad, 1);
    pf->interpolation_hess(c, coeff, hess, 1);

    for (unsigned dg = 0; dg < N; ++dg) {
      grad_fd[dg] = 
	eval_fem_gradient_with_finite_differences(pf, coeff, cv, pgt, gic, G, X0, h, dg);
      for (unsigned dh = 0; dh < N; ++dh) {
	hess_fd(0,dg*N+dh) = 
	  eval_fem_hessian_with_finite_differences(pf, coeff, cv, pgt, gic, G, X0, h, dg, dh);
      }
    }
    
    scalar_type err_grad = 
      gmm::vect_dist2((base_vector&)grad, (base_vector&)grad_fd);
    scalar_type err_hess = 
      gmm::vect_dist2((base_vector&)hess, (base_vector&)hess_fd);
    
    if (err_grad > 1e-4 ||
	err_hess > 1e-4) {
      cout << "validate_fem_derivatives dof=" << idof << "/" << pf->nb_dof(cv) << " -- X0ref = " << X0ref << "\n";

      if (gmm::vect_dist2((base_vector&)grad, (base_vector&)grad_fd) > 1e-4)
	cout << "grad = " << (base_vector&)grad << "\ngrad_fd = " << (base_vector&)grad_fd << "\n";
      cout << "hess = " << (base_vector&)hess << "\nhess_fd = " << (base_vector&)hess_fd << "\n";
      if (err_grad + err_hess > 1.0) { cout << "---------> COMPLETEMENT FAUX!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"; abort(); }
    }
  }
}


void validate_fem_derivatives(const getfem::mesh_fem &mf) {
  bgeot::base_matrix G;
  for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
    //if (mf.nb_dof_of_element(cv) > 12) {
      vectors_to_base_matrix(G, mf.linked_mesh().points_of_convex(cv));
      validate_fem_derivatives(mf.fem_of_element(cv), cv, mf.linked_mesh().trans_of_convex(cv), G);
      //}
  }
}






/*                                                          */
/*****  Methods for class bilaplacian_crack_problem  ********/
/*                                                          */





void bilaplacian_crack_problem::init(void) {
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
    
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
					 "Name of simplex integration method");
  std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION");
  enrichment_option = PARAM.int_value("ENRICHMENT_OPTION",
				      "Enrichment option");
    enr_area_radius = PARAM.real_value("RADIUS_ENR_AREA",
				     "radius of the enrichment area");

sol_ref = PARAM.int_value("SOL_REF") ;
    
    /* First step : build the mesh */

 if (!MESH_FILE.empty()) {
    mesh.read_from_file(MESH_FILE);
    base_small_vector tt(N);
    tt[0] = PARAM.real_value("TRANSLAT_X") ;  
    tt[1] = PARAM.real_value("TRANSLAT_Y") ; 
    if (sol_ref == 1){
       tt[0] = - PARAM.real_value("CRACK_SEMI_LENGTH") ;
    }
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
    if (sol_ref == 1)  nsubdiv[0] = NX / 2 ; 
    getfem::regular_unit_mesh(mesh, nsubdiv, pgt, PARAM.int_value("MESH_NOISED") != 0);

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
    
 if (PARAM.int_value("MOVE_NODES")){
    cout << "d�placement des noeuds \n" ;
//    size_type nb_x_pos, nb_y_pos = 0 ;
    scalar_type seuil_select = PARAM.real_value("SEUIL_SELECT") ;
//    scalar_type seuil_move = PARAM.real_value("SEUIL_MOVE") ;
    
    
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

  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);
  getfem::pintegration_method sppi = 
    getfem::int_method_descriptor(SIMPLEX_INTEGRATION);
  getfem::pintegration_method sing_ppi = (SINGULAR_INTEGRATION.size() ?
		getfem::int_method_descriptor(SINGULAR_INTEGRATION) : 0);

  mim.set_integration_method(mesh.convex_index(), ppi);
  mls.add_level_set(ls);
  mim.set_simplex_im(sppi, sing_ppi);

  /* Setting the finite element on the mf_u */  
  mf_pre_u.set_finite_element(mesh.convex_index(), pf_u); 
  getfem::pfem pf_partition_of_unity = getfem::fem_descriptor(PARAM.string_value("PARTITION_OF_UNITY_FEM_TYPE")) ;
  mf_partition_of_unity.set_finite_element(mesh.convex_index(), pf_partition_of_unity);

  // set the mesh_fem of the multipliers (for the dirichlet condition)
  std::string dirichlet_fem_name = PARAM.string_value("DIRICHLET_FEM_TYPE");
  if (dirichlet_fem_name.size() == 0)
    mf_mult.set_finite_element(mesh.convex_index(), pf_u);
  else {
    cout << "DIRICHLET_FEM_TYPE="  << dirichlet_fem_name << "\n";
    mf_mult.set_finite_element(mesh.convex_index(), 
			       getfem::fem_descriptor(dirichlet_fem_name));
  }
  std::string dirichlet_der_fem_name
    = PARAM.string_value("DIRICHLET_DER_FEM_TYPE", "");  // for the dirichlet condition on the derivative
  if (dirichlet_der_fem_name.size() == 0)
    mf_mult_d.set_finite_element(mesh.convex_index(), pf_u);
  else {
    cout << "DIRICHLET_DER_FEM_TYPE="  << dirichlet_der_fem_name << "\n";
    mf_mult_d.set_finite_element(mesh.convex_index(), 
			     getfem::fem_descriptor(dirichlet_der_fem_name));
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

  // set the mesh_fem for the multipliers for the case the Integral Matching  
  mf_pre_mortar.set_finite_element(mesh.convex_index(),
             getfem::fem_descriptor(PARAM.string_value("MORTAR_FEM_TYPE")));
  mf_pre_mortar_deriv.set_finite_element(mesh.convex_index(),
             getfem::fem_descriptor(PARAM.string_value("MORTAR_DERIV_FEM_TYPE")));

  /* set boundary conditions: non-homogeneous Dirichlet */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  if (sol_ref == 0){
     for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
         mesh.region(SIMPLE_SUPPORT_BOUNDARY_NUM).add(i.cv(), i.f());
         mesh.region(CLAMPED_BOUNDARY_NUM).add(i.cv(), i.f()); 
     }
  }
  
  if (sol_ref == 1 ){
     for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
        base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
        un /= gmm::vect_norm2(un);
        //if ( (un[0] <= -0.9) || (un[0] >= 0.9) ) {
	if  ( -un[0] >= 0.999  )  {
	   mesh.region(CLAMPED_BOUNDARY_NUM).add(i.cv(), i.f());
        
	}
	if  (gmm::abs(un[1]) >= 0.999)  {
	   mesh.region(SIMPLE_SUPPORT_BOUNDARY_NUM).add(i.cv(), i.f());
	   mesh.region(MOMENTUM_BOUNDARY_NUM).add(i.cv(), i.f());
	}  
     }
  }

    if (sol_ref == 2 ){
     for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
        base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
        un /= gmm::vect_norm2(un);
        //if ( (un[0] <= -0.9) || (un[0] >= 0.9) ) {
	if  (gmm::abs(un[1]) >= 0.999)  {
	   mesh.region(SIMPLE_SUPPORT_BOUNDARY_NUM).add(i.cv(), i.f());
           mesh.region(MOMENTUM_BOUNDARY_NUM).add(i.cv(), i.f());
	} 
     }
    }
  

  exact_sol.init(ls);
  
  cout << "initialisation de la level-set : \n" ;
  
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
     ls.values(1)[d] = x  ; //x * x - a * a ;
    }
    if (sol_ref == 2){ // to modify if rotation is supported
     ls.values(0)[d] = y ;
     ls.values(1)[d] = gmm::abs(x) - a ; //x * x - a * a ;
    }
  }
  //ls.simplify(0.5);
  ls.touch();  
  mls.adapt();
  mim.adapt();
  mfls_u.adapt();
//   mfls_mult.adapt();
//   mfls_mult_d.adapt();
  mfls_mortar.adapt();
  mfls_mortar_deriv.adapt();
  cout << "mfls_u.nb_dof()=" << mfls_u.nb_dof() << "\n";
}


/* compute the error with respect to the exact solution */
void bilaplacian_crack_problem::compute_error(plain_vector &U) {

  if (PARAM.real_value("RADIUS_SPLIT_DOMAIN") == 0){
     plain_vector V(gmm::vect_size(U)) ;
     gmm::clear(V) ;

    cout << "\nL2 ERROR:"
       << getfem::asm_L2_dist(mim, mf_u(), U,
			      exact_sol.mf, exact_sol.U) << "\n";
    cout << "H1 ERROR:"
         << getfem::asm_H1_dist(mim, mf_u(), U,
  	  		      exact_sol.mf, exact_sol.U) << "\n";
    cout << "H2 ERROR:"
         << getfem::asm_H2_dist(mim, mf_u(), U, 
                              exact_sol.mf, exact_sol.U) << "\n"; 
    if ( PARAM.int_value("NORM_EXACT") ){
    cout << "L2 exact:"
         << getfem::asm_L2_dist(mim, mf_u(), V,
			      exact_sol.mf, exact_sol.U) << "\n";
    cout << "H1 exact:"
         << getfem::asm_H1_dist(mim, mf_u(), V,
			      exact_sol.mf, exact_sol.U) << "\n";
    cout << "H2 exact:"
         << getfem::asm_H2_dist(mim, mf_u(), V, 
                              exact_sol.mf, exact_sol.U) << "\n";
    } 
  }
  else {
  getfem::mesh_region r_center, r_ext ;
  scalar_type radius_split_domain = PARAM.real_value("RADIUS_SPLIT_DOMAIN") ;
  bool in_area ;
  for (dal::bv_visitor cv(mesh.convex_index()) ; !cv.finished() ; ++cv){
	in_area = true;
	/* For each element, we test all of its nodes. 
	   If all the nodes are inside the enrichment area,
	   then the element is completly inside the area too */ 
	for (unsigned j=0; j < mesh.nb_points_of_convex(cv); ++j) {
	  if (gmm::sqr(mesh.points_of_convex(cv)[j][0] ) + 
	      gmm::sqr(mesh.points_of_convex(cv)[j][1] ) > 
	      gmm::sqr(radius_split_domain)) 
	    in_area = false; break; 
	}
	  if (in_area) r_center.add(cv) ;
	  else r_ext.add(cv) ;
  }
  scalar_type L2_center, H1_center, H2_center;
  cout << "ERROR SPLITTED - RADIUS =  " << radius_split_domain << "\n";
  cout << "Error on the crack tip zone:\n" ;
        L2_center = getfem::asm_L2_dist(mim, mf_u(), U, exact_sol.mf, exact_sol.U, r_center) ;
  cout << "  L2 error:" << L2_center << "\n";
	H1_center = getfem::asm_H1_dist(mim, mf_u(), U, exact_sol.mf, exact_sol.U, r_center) ;
  cout << "  H1 error:" << H1_center << "\n";
	H2_center = getfem::asm_H2_dist(mim, mf_u(), U, exact_sol.mf, exact_sol.U, r_center) ;
  cout << "  H2 error:" << H2_center << "\n";
 	
  cout << "Error on the remaining part of the domain:\n"; 
  scalar_type L2_ext, H1_ext, H2_ext;
        L2_ext = getfem::asm_L2_dist(mim, mf_u(), U, exact_sol.mf, exact_sol.U, r_ext) ;
  cout << "  L2 error:" << L2_ext << "\n";
	H1_ext = getfem::asm_H1_dist(mim, mf_u(), U, exact_sol.mf, exact_sol.U, r_ext) ;
  cout << "  H1 error:" << H1_ext << "\n";
	H2_ext = getfem::asm_H2_dist(mim, mf_u(), U, exact_sol.mf, exact_sol.U, r_ext) ;
  cout << "  H2 error:" << H2_ext << "\n";

  cout << "Error on the hole domain:\n";
  cout << "L2 ERROR:"
       << gmm::sqrt( gmm::sqr(L2_center) + gmm::sqr(L2_ext) ) << "\n";

    cout << "H1 ERROR:"
         << gmm::sqrt( gmm::sqr(H1_center) + gmm::sqr(H1_ext) ) << "\n";
    cout << "H2 ERROR:"
         << gmm::sqrt( gmm::sqr(H2_center) + gmm::sqr(H2_ext) ) << "\n";
  }
}




/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool bilaplacian_crack_problem::solve(plain_vector &U) {  
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  
  // Setting the level-set
  ls.reinit();  
  for (size_type d = 0; d < ls.get_mesh_fem().nb_dof(); ++d) {
    scalar_type x = ls.get_mesh_fem().point_of_dof(d)[0];
    scalar_type y = ls.get_mesh_fem().point_of_dof(d)[1];
    ls.values(0)[d] = y  ; // + 1/4.*(x + .25);
    ls.values(1)[d] = x;
  }
  //ls.simplify(0.5);
  ls.touch();  
  mls.adapt();
  mim.adapt();
  mfls_u.adapt();
//   mfls_mult.adapt();
//   mfls_mult_d.adapt();
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
	for (size_type i = 0 ; i < ufunc.size()  ; ++i) {
        	ufunc[i] = bilaplacian_crack_singular(i + 4, ls, nu);
	}
  mf_sing_u.set_functions(ufunc);
  }

  
  
  // Setting the enrichment --------------------------------------------/
   
  switch(enrichment_option) {
  case 0 :  // No enrichment
    {
    mf_u_sum.set_mesh_fems(mfls_u);
      // an optional treatment : exporting a representation of the mesh     
      getfem::mesh_fem mf_enrich(mesh);
      getfem::pfem pf_mef = getfem::classical_fem(mesh.trans_of_convex(mesh.convex_index().first_true()), 1 );
      mf_enrich.set_finite_element(mesh.convex_index(), pf_mef) ;
      std::vector<scalar_type> UU(mf_enrich.nb_dof()) ;
      std::fill(UU.begin(), UU.end() ,0.) ;
      cout << "exporting mesh to " << "mesh_representation.vtk" << "..\n";
      getfem::vtk_export exp("mesh_representation.vtk", false);
      exp.exporting(mf_enrich); 
      exp.write_point_data(mf_enrich, UU, "mesh");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d mesh_representation.vtk -f "
	"WarpScalar -m BandedSurfaceMap -m Outline\n";
    }
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
      for (dal::bv_visitor i(mesh.convex_index()) ; !i.finished() ; ++i){
          getfem::pfem pf_mef = getfem::classical_fem(mesh.trans_of_convex(i), 1 );
          mf_enrich.set_finite_element(i, pf_mef) ;
      }
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
 
    }  // end of "enrichment_option = 3"
    break ;
  default : 
	GMM_ASSERT1(false, "Enrichment_option parameter is undefined");
	break ;  
	}  // end of switch
	
  mesh.write_to_file("toto.mesh");
  
  if (PARAM.int_value("SHOW_NAME_OF_DOF")==1) {  // printing the type of each dof
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

  //cout << "validate mf_sing_u():\n"; validate_fem_derivatives(mf_sing_u);

  //cout << "validate mf_u():\n"; validate_fem_derivatives(mf_u());
  
  cout << "Number of dof for u: " << mf_u().nb_dof() << endl;

  // Bilaplacian brick.
  getfem::mdbrick_bilaplacian<> BIL(mim, mf_u());
  BIL.D().set(D);
  if (KL) { BIL.set_to_KL(); BIL.nu().set(nu); }
  
  // Defining the normal derivative Dirichlet condition value.
  plain_vector F;


  /* WRONG !! 

    F.resize(nb_dof_rhs*N);
  getfem::interpolation_function(mf_rhs, F, sol_du, CLAMPED_BOUNDARY_NUM);  
   
  // Normal derivative Dirichlet condition brick. 
  getfem::mdbrick_normal_derivative_Dirichlet<>                   
    NDER_DIRICHLET(BIL, CLAMPED_BOUNDARY_NUM, mf_mult);       
 
  NDER_DIRICHLET.set_constraints_type(dirichlet_version);
  NDER_DIRICHLET.rhs().set(mf_rhs, F);
  */

  //gmm::resize(U, mf_u().nb_dof());  return true;


  // Normal derivative Dirichlet condition brick. 
  getfem::mdbrick_normal_derivative_Dirichlet<>                   
    NDER_DIRICHLET(BIL, CLAMPED_BOUNDARY_NUM, mf_mult_d);  // mfls_mult_d  
  NDER_DIRICHLET.set_constraints_type(dirichlet_version);
  NDER_DIRICHLET.R_must_be_derivated(); // hence we give the exact solution , and its gradient will be taken
  NDER_DIRICHLET.rhs().set(exact_sol.mf,exact_sol.U);
  if (dirichlet_version == getfem::PENALIZED_CONSTRAINTS)
    NDER_DIRICHLET.set_penalization_parameter(PARAM.real_value("EPS_DIRICHLET_PENAL")) ;
  
  // Defining the Dirichlet condition value.
  gmm::resize(F, nb_dof_rhs);
  getfem::interpolation_function(mf_rhs, F, sol_u,SIMPLE_SUPPORT_BOUNDARY_NUM);

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<>
    DIRICHLET(NDER_DIRICHLET, SIMPLE_SUPPORT_BOUNDARY_NUM, mf_mult); //mfls_mult
  DIRICHLET.rhs().set(exact_sol.mf,exact_sol.U);
  DIRICHLET.set_constraints_type(getfem::constraints_type(dirichlet_version));  
  if (dirichlet_version == getfem::PENALIZED_CONSTRAINTS)
    DIRICHLET.set_penalization_parameter(PARAM.real_value("EPS_DIRICHLET_PENAL")) ;
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
	unsigned n = gmm::mat_nrows(H);
	gmm::resize(H, n+1, gmm::mat_ncols(H));
	H(n, d) = 1;
      }
    }  
    gmm::resize(R, gmm::mat_nrows(H)); 
    mortar.set_constraints(H,R);
    final_model = &mortar;
    gmm::Harwell_Boeing_save("H.hb", H);

  }

  if ( PARAM.real_value("SEUIL_FINAL")!=0 ) { 
  // suppression of nodes with a very small term on the stiffness matrix diag
    getfem::mdbrick_constraint<> &extra = *(new getfem::mdbrick_constraint<>(DIRICHLET, 0));
    extra.set_constraints_type(getfem::constraints_type(dirichlet_version));  
    if (dirichlet_version == getfem::PENALIZED_CONSTRAINTS)
      extra.set_penalization_parameter(PARAM.real_value("EPS_DIRICHLET_PENAL"));

    sparse_matrix M2(mf_u().nb_dof(), mf_u().nb_dof());
    sparse_matrix H(0, mf_u().nb_dof());
    //getfem::asm_mass_matrix(M2, mim, mf_u(), mf_u());
    base_vector RR(mf_rhs.nb_dof(), 1.0);
    getfem::asm_stiffness_matrix_for_bilaplacian(M2, mim, mf_u(), 
                                                 mf_rhs, RR);
						 
    //cout << "stiffness_matrix_for_bilaplacian : " << M2 << "\n" ;
    cout << "termes diagonaux, de la matrice de rigidit�, inf�rieurs � 1e-10 : " ;
    for (size_type d = 0; d < mf_u().nb_dof(); ++d) {
        if (M2(d,d) < 1e-10) cout << M2(d,d) << " ; " ;
    }  
    cout << "\n" ;
    cout << "SEUIL_FINAL = " << PARAM.real_value("SEUIL_FINAL") << "\n" ;
    for (size_type d = 0; d < mf_u().nb_dof(); ++d) {
      if (M2(d,d) < PARAM.real_value("SEUIL_FINAL")) {
	cout << "OULALA " << d << " @ " << mf_u().point_of_dof(d) << " : " << M2(d,d) << "\n";	
        unsigned n = gmm::mat_nrows(H);
	gmm::resize(H, n+1, gmm::mat_ncols(H));
	H(n, d) = 1;
      }
    }
    base_vector R(gmm::mat_nrows(H)); 
    extra.set_constraints(H,R);
    final_model = &extra;
    gmm::Harwell_Boeing_save("M2.hb", M2);
  }

  
  cout << "Total number of variables : " << final_model->nb_dof() << endl;
  getfem::standard_model_state MS(*final_model);
  gmm::iteration iter(residual, 1, 40000);


  // getfem::useful_types<getfem::standard_model_state>::plsolver_type p;
  // p.reset(new getfem::linear_solver_cg_preconditioned_ildlt<sparse_matrix,plain_vector>);

  getfem::standard_solve(MS, *final_model, iter /* , p*/);

  // Solution extraction
  gmm::resize(U, mf_u().nb_dof());
  gmm::copy(BIL.get_solution(MS), U);
  return true;  
}

template<typename VEC1, typename VEC2>
void asm_H2_semi_dist_map(const getfem::mesh_im &mim, 
                          const getfem::mesh_fem &mf1, const VEC1 &U1,
                          const getfem::mesh_fem &mf2, const VEC2 &U2,
                          const getfem::mesh_fem &mf_P0, VEC1 &V,
                          getfem::mesh_region rg = getfem::mesh_region::all_convexes()) {
  mim.linked_mesh().intersect_with_mpi_region(rg);
  getfem::generic_assembly assem;    
  assem.set("u1=data$1(#1); u2=data$2(#2); "
            "V(#3)+=u1(i).u1(j).comp(Hess(#1).Hess(#1).Base(#3))(i,d,e,j,d,e,:)"
            "+ u2(i).u2(j).comp(Hess(#2).Hess(#2).Base(#3))(i,d,e,j,d,e,:)"
            "- 2*u1(i).u2(j).comp(Hess(#1).Hess(#2).Base(#3))(i,d,e,j,d,e,:)");
  
  assem.push_mi(mim);
  assem.push_mf(mf1);
  assem.push_mf(mf2);
  assem.push_mf(mf_P0);
  assem.push_data(U1);
  assem.push_data(U2);
  assem.push_vec(V);
  assem.assembly(rg);
}
  
void bilaplacian_crack_problem::compute_H2_error_field(const plain_vector &U) {

    getfem::mesh_fem mf_P0(mesh);
    mf_P0.set_finite_element(mesh.convex_index(), getfem::classical_fem(mesh.trans_of_convex(0), 0));
    plain_vector V(mf_P0.nb_dof());
    asm_H2_semi_dist_map(mim, mf_u(), U, exact_sol.mf, exact_sol.U, mf_P0, V);
    cout << "exporting H2 error map\n";
    getfem::vtk_export exp2(datafilename + "_H2.vtk");
    exp2.exporting(mf_P0);
    exp2.write_point_data(mf_P0, V, "H2 error map");
    
    mf_P0.write_to_file(datafilename + "_H2.meshfem", true);
    gmm::vecsave(datafilename + "_H2.V", V);
}

