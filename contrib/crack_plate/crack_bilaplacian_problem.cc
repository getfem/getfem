#include "crack_bilaplacian.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_fourth_order.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_superlu.h"

/******** Exact Solution *******************************/

scalar_type D  = 1.  ;
scalar_type nu = 0.3 ;
scalar_type BB = 0.1 ;
scalar_type AAA = BB * (3. * nu + 5.)/ (3. * (nu - 1.))   ;  // (-3.0+nu*nu-2.0*nu)/(nu*nu-2.0*nu+5.0);
scalar_type DD = 0.0 ;
scalar_type CC = DD * (nu + 7.)/ (3. * (nu - 1.))   ;   //  (-8.0*nu+3.0*BB*nu*nu-6.0*nu*BB+15.0*BB)/(nu*nu-2.0*nu+5.0);
 

scalar_type sol_u(const base_node &x){
 scalar_type r = sqrt( x[0] * x[0] + x[1] * x[1] ) ;
 //scalar_type theta = 2. * atan( x[1] / ( x[0] + r ) ) ;
 scalar_type theta = atan2(x[1], x[0]);
 //return sqrt(r*r*r) * (sin(3.0/2.0*theta)+AAA*sin(theta/2.0)+BB*cos(3.0/2.0*theta)+CC*cos(theta/2.0));
 return sqrt(r*r*r)*(AAA*sin(3.0/2.0*theta)+BB*sin(theta/2.0)+CC*cos(3.0/2.0*theta)+DD*cos(theta/2.0));
}
 
scalar_type sol_lapl_u(const base_node &x) {
 scalar_type r = sqrt( x[0] * x[0] + x[1] * x[1] ) ;
 scalar_type theta = atan2(x[1], x[0]);
 return 2.0*(BB*sin(theta/2.0)+DD*cos(theta/2.0))/sqrt(r);
 /* return 9.0/4.0/sqrt(r)*(sin(3.0/2.0*theta)+AAA*sin(theta/2.0)+BB*cos(3.0/2.0*theta)+CC*cos(theta/2.0))+1/sqrt(r)*(-9.0/4.0*sin(3.0/2.0*theta)-AAA*sin(theta/
    2.0)/4.0-9.0/4.0*BB*cos(3.0/2.0*theta)-CC*cos(theta/2.0)/4.0); */ }

scalar_type sol_f(const base_node &)
{ return 0. ; }

base_small_vector sol_du(const base_node &x) {
 base_small_vector res(x.size());
 scalar_type r = sqrt( x[0] * x[0] + x[1] * x[1] ) ;
 scalar_type theta = atan2(x[1], x[0]);
res[0] = 3.0/2.0*sqrt(r)*(AAA*sin(3.0/2.0*theta)+BB*sin(theta/2.0)+CC*cos(3.0/2.0*theta)
+DD*cos(theta/2.0))*cos(theta)-sqrt(r)*(3.0/2.0*AAA*cos(3.0/2.0*theta)+
BB*cos(theta/2.0)/2.0-3.0/2.0*CC*sin(3.0/2.0*theta)-DD*sin(theta/2.0)/2.0)*sin(theta );

res[1] = 3.0/2.0*sqrt(r)*(AAA*sin(3.0/2.0*theta)+BB*sin(theta/2.0)+CC*cos(3.0/2.0*theta)
+DD*cos(theta/2.0))*sin(theta)+sqrt(r)*(3.0/2.0*AAA*cos(3.0/2.0*theta)+
BB*cos(theta/2.0)/2.0-3.0/2.0*CC*sin(3.0/2.0*theta)-DD*sin(theta/2.0)/2.0)*cos(theta);

/*
res[0] =  3.0/2.0*sqrt(r)*(sin(3.0/2.0*theta)+AAA*sin(theta/2.0)+BB*cos(3.0/2.0*theta)+CC*cos(theta/2.0))*cos(theta)-sqrt(r)*(3.0/2.0*cos(3.0/2.0*theta)+AAA*
cos(theta/2.0)/2.0-3.0/2.0*BB*sin(3.0/2.0*theta)-CC*sin(theta/2.0)/2.0)*sin(
theta);
 
res[1] = 3.0/2.0*sqrt(r)*(sin(3.0/2.0*theta)+AAA*sin(theta/2.0)+BB*cos(3.0/2.0*theta)+CC*cos(theta/2.0))*sin(theta)+sqrt(r)*(3.0/2.0*cos(3.0/2.0*theta)+AAA*
cos(theta/2.0)/2.0-3.0/2.0*BB*sin(3.0/2.0*theta)-CC*sin(theta/2.0)/2.0)*cos(
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
    cfun[j] = bilaplacian_crack_singular(j, ls) ;
  mf.set_functions(cfun);
  U.resize(4); assert(mf.nb_dof() == 4);
  // scalar_type A1 = 1., nu = 0.3 ;
  // scalar_type b1_ = 3. + (A2 / A1) * (24. * nu) / (3. * nu * nu - 6. * nu + 5. ) ; 
  U[0] = BB ;
  U[1] = AAA ;
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




// functions for assembling the constraints of the integral matching 
namespace getfem {  
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

} // end namespace getfem

scalar_type crack_level(base_node P) {
return P[1];
}

scalar_type crack_tip_level(base_node P) {
return P[0];
}

/*                                                          */
/*****  Methods for class bilaplacian_crack_problem  ********/
/*                                                          */





void bilaplacian_crack_problem::init(void) {
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
  
  size_type N_tri = 2, N_quad = 2 ;
    
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
					 "Name of simplex integration method");
  std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION");
  enrichment_option = PARAM.int_value("ENRICHMENT_OPTION",
				      "Enrichment option");
    enr_area_radius = PARAM.real_value("RADIUS_ENR_AREA",
				     "radius of the enrichment area");
    
    /* First step : build the mesh */
    
    bgeot::pgeometric_trans pgt_tri = 
      bgeot::geometric_trans_descriptor(TRI_MESH_TYPE);
    bgeot::pgeometric_trans pgt_quad = 
      bgeot::geometric_trans_descriptor(QUAD_MESH_TYPE);
    N_tri = pgt_tri->dim();
    N_quad = pgt_quad->dim();
    GMM_ASSERT1(N_tri == 2, "For a plate problem, N should be 2");
    GMM_ASSERT1(N_quad == 2, "For a plate problem, N should be 2");
    std::vector<size_type> nsubdiv(N_tri);
    NX = PARAM.int_value("NX", "Number of space steps ") ;
    std::fill(nsubdiv.begin(),nsubdiv.end(), NX);
    if (!PARAM.int_value("QUAD"))
    getfem::regular_unit_mesh(mesh, nsubdiv, pgt_tri, PARAM.int_value("MESH_NOISED") != 0);
    else 
    getfem::regular_unit_mesh(mesh, nsubdiv, pgt_quad, PARAM.int_value("MESH_NOISED") != 0);
    
    bgeot::base_matrix M(N_tri,N_tri);
    for (size_type i=0; i < N_tri; ++i) {
      static const char *t[] = {"LX","LY","LZ"};
      M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
    }
    /* scale the unit mesh to [LX,LY,..] and incline it */
    mesh.transformation(M);
    
    base_small_vector tt(N_tri); 
    tt[0] = -0.5 ;
    tt[1] = -0.5;
    mesh.translation(tt); 
  
    dal::bit_vector quad_among_cvx, tri_among_cvx ;
    
    if (PARAM.real_value("QUAD") == 1 ) {
       cout << "transforming quad in two triangles\n";
       unsigned cpt ;
       base_small_vector values_ls(4), values_ls_tip(4) ;
       std::vector<bgeot::size_type> ind(3) ;
       std::vector<size_type> index_nodes_of_a_face(2), ind_sommets_cvx(4) ; 
       // For each quadrangle, we test wether or not it is crossed by the crack.
       // If it is so, we delete the quadrangle and create two triangles.
       for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i){
       
       
       
	   // we store the values of the level-set on the four nodes
	   //cout << "(mesh.points_of_convex(i)).size() = " << (mesh.points_of_convex(i)).size() << "\n";
	   for (unsigned node = 0 ; node < (mesh.points_of_convex(i)).size() ; ++node){
	       values_ls[node] = crack_level(mesh.points_of_convex(i)[node] ) ;
	       values_ls_tip[node] = crack_tip_level(mesh.points_of_convex(i)[node]) ;
	       }
	   // We iterate on the four faces. If 2 faces or more are
	   // crossed by the level-set, we conclude that the quadrangle is completly
	   // crossed by the level-set.
	   cpt = 0 ;
	   for (int j = 0 ; j < 4 ; j++) {
	       // we store the local indexes of the two nodes of the current face 
	       for (unsigned k = 0 ; k < (mesh.structure_of_convex(i)->ind_points_of_face(j)).size() ; k++ )
	         index_nodes_of_a_face[k] = (mesh.structure_of_convex(i)->ind_points_of_face(j))[k] ;
	       
	       //cout << "index_nodes_of_a_face = " << index_nodes_of_a_face << "\n" ;
		 
	       if ( ( values_ls[index_nodes_of_a_face[0]] * values_ls[index_nodes_of_a_face[1]] < 0 ) 
	           && (values_ls_tip[0] < 0) && (values_ls_tip[1] < 0) ) 
		  cpt++;
           }
	   
	   if (cpt >= 2) { // we add two triangles and delete the quadrangle
	      // we store the global indexes of the nodes of the i convex
	      //cout << "Convex " << i << "to split ; values of the level_set : " << values_ls << "\n" ;
	      for (unsigned j = 0 ; j < (mesh.ind_points_of_convex(i)).size() ; ++j)
	          ind_sommets_cvx[j] = mesh.ind_points_of_convex(i)[j] ;
              
	      mesh.sup_convex(i) ;   // delete pathologic quadrangle
              // Cut the quadrangle along the shortest diagonal
	      
	      std::vector<bgeot::base_node> P(4) ;
	      for (unsigned node = 0 ; node < (mesh.points_of_convex(i)).size() ; ++node)
	          P[node] = mesh.points_of_convex(i)[node] ;
	      scalar_type dist1 = 0., dist2 = 0. ;
	      for(int k = 0 ; k < 2 ; k++){ 
	          dist1 += (P[0][k] - P[3][k]) * (P[0][k] - P[3][k]) ;
		  dist2 += (P[2][k] - P[1][k]) * (P[2][k] - P[1][k]) ; 
	      }
	      
	      if ( dist1 > dist2){ //then cut along the [s1 s2] diagonal
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
	      else quad_among_cvx.add(i) ;
	   }
	   // save the information of what is triangle or quadrangle
	   for (dal::bv_visitor j(mesh.convex_index()); !j.finished(); ++j)
	       if (!quad_among_cvx[j]) tri_among_cvx.add(j) ;
// 	   cout << "convex indexes : " << mesh.convex_index() << "\n" ;
// 	   cout << "convex tri : " << tri_among_cvx << "\n" ;
// 	   cout << "convex quad : " << quad_among_cvx << "\n" ;
	   // The procedure used beyond was made possible to use thanks to the
	   // structure of the containers returned : the elements are ordered 
	   // in accordance with the numerotation of the vertex and faces in a quadrangle:
// 	   2____3        _2__
// 	   |    |     1 |    |  0
// 	   |____|       |____|
// 	   0    1         3
    }
    else {
    for(dal::bv_visitor i(mesh.convex_index()); !i.finished() ; ++i )
       tri_among_cvx.add(i) ;
    }
    
    
    scalar_type quality = 1.0;
    for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i)
    	quality = std::min(quality, mesh.convex_quality_estimate(i));
    cout << "quality of mesh : " << quality << endl;
    
   /* read the parameters   */
  epsilon = PARAM.real_value("EPSILON", "thickness") ;
  int dv = PARAM.int_value("DIRICHLET_VERSION", "Dirichlet version");
  int mv = PARAM.int_value("MORTAR_VERSION", "Mortar version");
  dirichlet_version = getfem::constraints_type(dv);
  mortar_version = getfem::constraints_type(mv);
  datafilename=PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual=PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;
  KL = (PARAM.int_value("KL", "Kirchhoff-Love model or not") != 0);
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
  mls.add_level_set(ls);
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
  if (dirichlet_fem_name_tri.size() == 0){
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
  if (dirichlet_der_fem_name_tri.size() == 0){
    mf_mult_d.set_finite_element(tri_among_cvx, pf_u_tri);
    mf_mult_d.set_finite_element(quad_among_cvx, pf_u_quad);
  }
  else {
    cout << "TRI_DIRICHLET_DER_FEM_TYPE="  << dirichlet_der_fem_name_tri << "\n";
    cout << "QUADDIRICHLET_DER_FEM_TYPE="  << dirichlet_der_fem_name_quad << "\n";
    mf_mult_d.set_finite_element(tri_among_cvx, 
			     getfem::fem_descriptor(dirichlet_der_fem_name_tri));
    mf_mult_d.set_finite_element(quad_among_cvx, 
			     getfem::fem_descriptor(dirichlet_der_fem_name_quad));
  }

  /* set the finite element on mf_rhs (same as mf_u if DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name_tri = PARAM.string_value("TRI_DATA_FEM_TYPE");
  std::string data_fem_name_quad = PARAM.string_value("QUAD_DATA_FEM_TYPE");
  if (data_fem_name_tri.size() == 0) {
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
  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
    mesh.region(SIMPLE_SUPPORT_BOUNDARY_NUM).add(i.cv(), i.f());
    mesh.region(CLAMPED_BOUNDARY_NUM).add(i.cv(), i.f()); 
  }
  
  exact_sol.init(ls);
}


/* compute the error with respect to the exact solution */
void bilaplacian_crack_problem::compute_error(plain_vector &U) {

  if (PARAM.real_value("RADIUS_SPLIT_DOMAIN") == 0){
// plain_vector V(gmm::vect_size(U)) ;
// gmm::clear(V) ;
//   cout << "L2 ERROR : "
//        << getfem::asm_L2_dist(mim, mf_u(), V,
// 			      exact_sol.mf, exact_sol.U) << "\n";
//   cout << "H1 ERROR : "
//        << getfem::asm_H1_dist(mim, mf_u(), V,
// 			      exact_sol.mf, exact_sol.U) << "\n";
//   cout << "H2 ERROR : "
//        << getfem::asm_H2_dist(mim, mf_u(), V, 
//                               exact_sol.mf, exact_sol.U) << "\n"; 

  cout << "L2 ERROR : "
       << getfem::asm_L2_dist(mim, mf_u(), U,
			      exact_sol.mf, exact_sol.U) << "\n";
  if (PARAM.int_value("ASM_H1_AND_H2_DIST") == 1){
    cout << "H1 ERROR : "
         << getfem::asm_H1_dist(mim, mf_u(), U,
  	  		      exact_sol.mf, exact_sol.U) << "\n";
    cout << "H2 ERROR : "
         << getfem::asm_H2_dist(mim, mf_u(), U, 
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
  cout << "Error on the crack tip zone : \n" ;
        L2_center = getfem::asm_L2_dist(mim, mf_u(), U, exact_sol.mf, exact_sol.U, r_center) ;
  cout << "  L2 error : " << L2_center << "\n";
	H1_center = getfem::asm_H1_dist(mim, mf_u(), U, exact_sol.mf, exact_sol.U, r_center) ;
  cout << "  H1 error : " << H1_center << "\n";
	H2_center = getfem::asm_H2_dist(mim, mf_u(), U, exact_sol.mf, exact_sol.U, r_center) ;
  cout << "  H2 error : " << H2_center << "\n";
 	
  cout << "Error on the remaining part of the domain : \n"; 
  scalar_type L2_ext, H1_ext, H2_ext;
        L2_ext = getfem::asm_L2_dist(mim, mf_u(), U, exact_sol.mf, exact_sol.U, r_ext) ;
  cout << "  L2 error : " << L2_ext << "\n";
	H1_ext = getfem::asm_H1_dist(mim, mf_u(), U, exact_sol.mf, exact_sol.U, r_ext) ;
  cout << "  H1 error : " << H1_ext << "\n";
	H2_ext = getfem::asm_H2_dist(mim, mf_u(), U, exact_sol.mf, exact_sol.U, r_ext) ;
  cout << "  H2 error : " << H2_ext << "\n";

  cout << "Error on the hole domain : \n";
  cout << "L2 ERROR : "
       << gmm::sqrt( gmm::sqr(L2_center) + gmm::sqr(L2_ext) ) << "\n";

    cout << "H1 ERROR : "
         << gmm::sqrt( gmm::sqr(H1_center) + gmm::sqr(H1_ext) ) << "\n";
    cout << "H2 ERROR : "
         << gmm::sqrt( gmm::sqr(H2_center) + gmm::sqr(H2_ext) ) << "\n";
  }
}




/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool bilaplacian_crack_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  //size_type N = mesh.dim();
  
  // Setting the level-set
  ls.reinit();  
  for (size_type d = 0; d < ls.get_mesh_fem().nb_dof(); ++d) {
    scalar_type x = ls.get_mesh_fem().point_of_dof(d)[0];
    scalar_type y = ls.get_mesh_fem().point_of_dof(d)[1];
    ls.values(0)[d] = y  ; //+ 1./4.*(x + .25);   // To modify in accordance with: scalar_type crack_level(bas_node P)
    ls.values(1)[d] = x;                          // idem, with: scalar_type crack_tip_level(bas_node P)
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
  std::vector<getfem::pglobal_function> ufunc(4);
  for (size_type i = 0 ; i < ufunc.size() ; ++i) {                              
    ufunc[i] = bilaplacian_crack_singular(i, ls);
  }
  mf_sing_u.set_functions(ufunc);
  
  
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
	     gmm::sqr(enr_area_radius)) {
	   pm_convexes.sup(i); break;
	 }
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
          cout << "enriched_dofs: " << enriched_dofs << "\n";
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
      if (PARAM.int_value("MORTAR_WITHOUT_SINGUL"))
         mf_u_sum.set_mesh_fems(mfls_u);
      else
         mf_u_sum.set_mesh_fems(mf_sing_u, mfls_u);
      
      cout << "cvlist_in_area: " << cvlist_in_area << "\n";
      cout << "mfls_u.nb_dof: " << mfls_u.nb_dof() << "\n";
      cout << "mf_u_sum.nb_dof: " << mf_u_sum.nb_dof() << "\n";
      cout << "MORTAR_BOUNDARY_IN: " << mesh.region(MORTAR_BOUNDARY_IN) << "\n";
      cout << "MORTAR_BOUNDARY_OUT: " << mesh.region(MORTAR_BOUNDARY_OUT) << "\n";
      
      // Searching the elements that are both crossed by the crack
      // and with one of his border which constitutes a part of the 
      // boundary between the enriched zone and the rest of the domain.
      getfem::mesh_region &boundary = mesh.region(MORTAR_BOUNDARY_IN);
      
      unsigned cpt = 0 ;
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
 
    }
    break ;
  default : 
	GMM_ASSERT1(false, "Enrichment_option parameter is undefined");
	break ;  
	}
  mesh.write_to_file("toto.mesh");
  
  if (0) {
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
    NDER_DIRICHLET(BIL, CLAMPED_BOUNDARY_NUM, mf_mult_d);    
  NDER_DIRICHLET.set_constraints_type(dirichlet_version);
  NDER_DIRICHLET.R_must_be_derivated(); // hence we give the exact solution , and its gradient will be taken
  NDER_DIRICHLET.rhs().set(exact_sol.mf,exact_sol.U);


  // Defining the Dirichlet condition value.
  gmm::resize(F, nb_dof_rhs);
  getfem::interpolation_function(mf_rhs, F, sol_u,SIMPLE_SUPPORT_BOUNDARY_NUM);

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<>
    DIRICHLET(NDER_DIRICHLET, SIMPLE_SUPPORT_BOUNDARY_NUM, mf_mult);
  DIRICHLET.rhs().set(exact_sol.mf,exact_sol.U);
  DIRICHLET.set_constraints_type(getfem::constraints_type(dirichlet_version));  
  if (dirichlet_version == getfem::PENALIZED_CONSTRAINTS)
    DIRICHLET.set_penalization_parameter(PARAM.real_value("EPS_DIRICHLET_PENAL")) ;
  
  getfem::mdbrick_abstract<> *final_model = &DIRICHLET ;
      
    getfem::mdbrick_constraint<> &mortar = 
      *(new getfem::mdbrick_constraint<>(DIRICHLET, 0));

          
//     getfem::mdbrick_constraint<> &extra = 
//       *(new getfem::mdbrick_constraint<>(mortar, 0));   

    getfem::mdbrick_constraint<> &extra = 
      *(new getfem::mdbrick_constraint<>(DIRICHLET, 0));   

  if (enrichment_option == 3 ) {
     /* add a constraint brick for the mortar junction between
       the enriched area and the rest of the mesh */
       int mult_with_H = PARAM.int_value("MULT_WITH_H") ;
       mortar_type = PARAM.int_value("MORTAR_TYPE") ;
       getfem::mesh_fem &mf_mortar = (mult_with_H == 1) ? mfls_mortar : mf_pre_mortar;
       getfem::mesh_fem &mf_mortar_deriv = (mult_with_H == 1) ? mfls_mortar_deriv : mf_pre_mortar_deriv;
       
//     getfem::mdbrick_constraint<> &mortar = 
//       *(new getfem::mdbrick_constraint<>(DIRICHLET, 0));
    mortar.set_constraints_type(getfem::constraints_type(mortar_version));  
    if (mortar_version == getfem::PENALIZED_CONSTRAINTS)
      mortar.set_penalization_parameter(PARAM.real_value("EPS_MORTAR_PENAL")) ;

    sparse_matrix H0(mf_mortar.nb_dof(), mf_u().nb_dof()), 
       H(mf_mortar_deriv.nb_dof(), mf_u().nb_dof());
    getfem::base_vector R(mf_mortar.nb_dof());
      
    
    if (mortar_type == 1) {
      // older version of integral matching (jan-feb 2007)      
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
      
      gmm::resize(H0, mf_mortar.nb_dof(), mf_u().nb_dof()) ;
      gmm::resize(H,  ind_union.size(),   mf_u().nb_dof()) ;
      
      cout << "bv_mortar = " << bv_mortar << "\n";
      cout << "bv_deriv = " << bv_deriv << "\n" ;
      cout << "bv_union = " << bv_union << "\n" ;
      
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
    else{
    
      // for the condition \int_\Gamma \mu \nabla \u = 0, with \mu vectorial
      mf_mortar_deriv.set_qdim(2) ;
      
      /* New version of the integral matching -----------------------------------------------*/
      
      // selecting nodes indices on the two meth. mult.
      
      // assembling matrices, selecting the lines corresponding to those holding information
      
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
      // 	  		mf_mortar_deriv, mf_pre_u, R, MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;
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
      //getfem::asm_normal_derivative_dirichlet_constraints
      asm_constraint_gradient_vectorial_mult
	(H0, mim, mf_u(), mf_mortar_deriv,
	 MORTAR_BOUNDARY_OUT, getfem::ASMDIR_BUILDH) ;
      gmm::add(gmm::sub_matrix(H0, sub_i1, sub_j),
	       gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
      
      cout << "first step \n" ;
      gmm::clear(H0);
      //getfem::asm_normal_derivative_dirichlet_constraints
      asm_constraint_gradient_vectorial_mult
	(H0, mim, mf_u(), mf_mortar_deriv, 
	 MORTAR_BOUNDARY_IN, getfem::ASMDIR_BUILDH) ;
      gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i1, sub_j), -1),
	       gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;    
	       
      // getfem::mdbrick_constraint<> &extra = *(new getfem::mdbrick_constraint<>(DIRICHLET, 0));
      
      
      // ------------------------------------------ end of new version
    }
    
    
    /* because of the discontinuous partition of mf_u(), some levelset 
       enriched functions do not contribute any more to the
       mass-matrix (the ones which are null on one side of the
       levelset, when split in two by the mortar partition, may create
       a "null" dof whose base function is all zero.. */
    sparse_matrix M2(mf_u().nb_dof(), mf_u().nb_dof());
    getfem::asm_mass_matrix(M2, mim, mf_u(), mf_u());
    //gmm::HarwellBoeing_IO::write("M2.hb", M2);
    for (size_type d = 0; d < mf_u().nb_dof(); ++d) {
      //       if (M2(d,d) < 1e-7) cout << "  weak mf_u() dof " << d << " @ " << 
      // 	  mf_u().point_of_dof(d) << " M2(d,d) = " << M2(d,d) << "\n";
      if (M2(d,d) < PARAM.real_value("SEUIL_MORTAR")) {
	cout << "removed : " << d << " @ " << mf_u().point_of_dof(d) << " : " << M2(d,d) << "\n";	
	unsigned n = gmm::mat_nrows(H);
	gmm::resize(H, n+1, gmm::mat_ncols(H));
	H(n, d) = 1;
      }
    }  
    gmm::resize(R, gmm::mat_nrows(H)); 
    mortar.set_constraints(H,R);
    final_model = &mortar;
    gmm::Harwell_Boeing_save("H.hb", H);        
    //------------------------------------------------------------------
    // Matching of the normal derivative
    //     getfem::mdbrick_constraint<> &mortar_derivative = 
    //       *(new getfem::mdbrick_constraint<>(mortar,0));
    //       
    //     gmm::clear(H0) ;
    //     gmm::clear(H) ;
    //     gmm::resize(H, ind_mortar.size(), mf_u().nb_dof());
    //     getfem::asm_normal_derivative_dirichlet_constraints(H0, R, mim, mf_u(),
    // 	     		mf_mortar_deriv, mf_pre_u, R, MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;
    //     gmm::copy(gmm::sub_matrix(H0, sub_i, sub_j), H) ;
    //     
    //     gmm::clear(H0) ;
    //     getfem::asm_normal_derivative_dirichlet_constraints(H0, R, mim, mf_u(),
    // 			mf_mortar_deriv, mf_pre_u, R, MORTAR_BOUNDARY_IN, 0, getfem::ASMDIR_BUILDH) ;
    //     gmm::add(gmm::sub_matrix(H0, sub_i, sub_j), H) ;
    //     
    //     //cout << "matrice de contraintes : \n" << H << "\n" ;
    //     cout << "matching of derivative values done.\n" ;  
    //     gmm::resize(R, ind_mortar.size());
    //     mortar_derivative.set_constraints(H,R);
    //     final_model = &mortar_derivative ;
  } 
  
  
  
  // suppression of nodes with a very small term on the mass matrix diag
  // (due to the fact that the elements are cut very close to the nodes by the level set).   
  //getfem::mdbrick_constraint<> &extra = *(new getfem::mdbrick_constraint<>(DIRICHLET, 0)); 
  
  if (PARAM.real_value("SEUIL") != 0. ) { 

    extra.set_constraints_type(getfem::constraints_type(dirichlet_version));  
    if (dirichlet_version == getfem::PENALIZED_CONSTRAINTS)
      extra.set_penalization_parameter(PARAM.real_value("EPS_DIRICHLET_PENAL"));

    sparse_matrix M2(mf_u().nb_dof(), mf_u().nb_dof());
    sparse_matrix H(0, mf_u().nb_dof());
    //getfem::asm_mass_matrix(M2, mim, mf_u(), mf_u());
    base_vector RR(mf_rhs.nb_dof(), 1.0);
    //getfem::asm_stiffness_matrix_for_bilaplacian(M2, mim, mf_u(), mf_rhs, RR);
    getfem::asm_mass_matrix(M2, mim, mf_u(), mf_u());
    
    for (size_type d = 0; d < mf_u().nb_dof(); ++d) {
      if (M2(d,d) < PARAM.real_value("SEUIL")) {
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
  getfem::standard_solve(MS, *final_model, iter);

  // Solution extraction
  gmm::resize(U, mf_u().nb_dof());
  gmm::copy(BIL.get_solution(MS), U);
  return (iter.converged());
}

// void bilaplacian_crack_problem::split_quadrangle(size_type cv_id, size_type dof) {
// // constitute the list of vertex of the 2 sub-triangles
// // itérer sur ce container mesh.ind_points_of_convex(i)
// //  
//  
// std::vector<bgeot::base_node> pts1(3), pts2(3);
// pts1[0] = mf_u().point_of_dof(dof) ;  
// pts2[0] = mf_u().point_of_dof(dof) ;  
// 
// pts[1] = mymesh.add_point(bgeot::base_node(0.0, 1.0, 0.0);
// pts[2] = mymesh.add_point(bgeot::base_node(0.0, 0.0, 1.0);      
// 
// 
// mesh.sup_convex(i) ; // delete the pathologic quadrangle
// bgeot::pgeometric_trans pgt = 
//       bgeot::geometric_trans_descriptor("GT_PK(2,1)");
// i = mesh.add_convex(pgt, it);
// 
// }

namespace getfem{
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

}
//function to save a vector for matlab
template<typename VEC> static void vecsave(std::string fname, const VEC& V) {
  std::ofstream f(fname.c_str()); f.precision(16);
  for (size_type i=0; i < V.size(); ++i) f << V[i] << "\n"; 
}
