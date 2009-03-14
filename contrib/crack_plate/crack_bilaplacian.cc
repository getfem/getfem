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

#include "crack_bilaplacian.h"
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_derivatives.h"

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

size_type is_global_dof_type(getfem::pdof_description dof){
size_type global_dof = 0 ;
   for (dim_type d = 0; d < 4 ; ++d){
       if (dof == getfem::global_dof(d)) {
          global_dof = 1;
	      }
   }
return global_dof ;
}
 
int main(int argc, char *argv[]) {

   try {
    bilaplacian_crack_problem p;
    p.PARAM.read_command_line(argc, argv);
    if (p.PARAM.int_value("MIXED_ELEMENTS"))
       p.init_mixed_elements();
    else
       p.init();
    plain_vector U;
    p.mesh.write_to_file("mesh.m") ;
    scalar_type ring_radius = p.PARAM.real_value("RING_RADIUS");
    if (p.PARAM.int_value("SOL_REF") == 0) {
       if (!p.solve(U)) GMM_ASSERT1(false, "Solve has failed");
       p.compute_sif(U, ring_radius);
       if (p.PARAM.int_value("COMPUTE_ERROR") == 1)     		   p.compute_error(U) ;
    }
    if (p.PARAM.int_value("SOL_REF") == 1) {
       if (!p.solve_moment(U)) GMM_ASSERT1(false, "Solve has failed");
        p.compute_sif(U, ring_radius);
/*       cout << "valeur des dofs sur les bords horizontaux : \n" ;
       for (size_type i = 0; i< p.mf_u().nb_dof() ; i++){
          if ( gmm::abs(p.mf_u().point_of_dof(i)[1]) > 0.499){
	     cout << U[i] << " ; \n" ;
	  }
       }*/ 
    }
    if (p.PARAM.int_value("SOL_REF") == 2) {
       if (!p.solve_moment(U)) GMM_ASSERT1(false, "Solve has failed");
       p.compute_error_beta(U) ;
    }
    
    //p.compute_H2_error_field(U) ;
    
    // -------------------------------------------------------------------
    // Print estimation of SIF (from the value of the coefficient 
    //  of the singularity, on integral matching)
    // -------------------------------------------------------------------
    
    unsigned q = p.mf_u().get_qdim();
    
    base_small_vector tab_fic(4) ;
    unsigned cpt = 0;
    if (p.PARAM.int_value("ENRICHMENT_OPTION") > 2){
    // affichage des coeffs devant les singularites, avec le raccord integral
	for (size_type d=0; d < p.mf_u().nb_dof(); d += q) {
		size_type cv = p.mf_u().first_convex_of_dof(d) ;
		getfem::pfem pf = p.mf_u().fem_of_element(cv);
		unsigned ld = unsigned(-1);
		for (unsigned dd = 0; dd < p.mf_u().nb_dof_of_element(cv); dd += q) {
		if (p.mf_u().ind_dof_of_element(cv)[dd] == d) {
			ld = dd/q; break;
		}
		}   
		if (ld == unsigned(-1)) {
		cout << "DOF " << d << "NOT FOUND in " << cv << " BUG BUG\n";
		} 
		else {
		if ( is_global_dof_type(pf->dof_types().at(ld)) ){
	// 	         printf("dof %4d @ %+6.2f:%+6.2f: ", d, p.mf_u().point_of_dof(d)[0], p.mf_u().point_of_dof(d)[1]);
	//                  printf(" %3d:%.16s", cv, name_of_dof(pf->dof_types().at(ld)).c_str());
	// 	         cout << " coeff donant le FIC;" << U[d] << "\n";
			cout << "coeff:" << U[d] << "\n" ;
			tab_fic[cpt] = U[d] ;
			cpt +=1 ;
			}
		}
	}
    // calcul des FIC k1 et k2 :
    scalar_type E, nu2 ;
    base_small_vector k(2) ;
    nu2 = p.nu * p.nu ;
    E = 3. * (1. - nu2) * p.D / (2. * p.epsilon * p.epsilon * p.epsilon) ;
    if (p.PARAM.int_value("SING_BASE_TYPE") == 0){
       k[0] = tab_fic[2] * (3.*(p.nu - 1.)) + tab_fic[3] * (3. * p.nu + 5.) ;
       k[0] *= - ( sqrt(2.) * E * p.epsilon ) / (4. * (1. - nu2) ) ;
       k[1] = (tab_fic[0] + 3. * tab_fic[1]) ;
       k[1] *= - E * p.epsilon * sqrt(2.) * (3. + p.nu) / (4. * (1. + p.nu) * (1. + p.nu) ) ;
    }
    if (p.PARAM.int_value("SING_BASE_TYPE") == 1){
       scalar_type A ;
       A =  sqrt(2) * p.epsilon * E * (p.nu + 3.) / (1. - nu2) ;
       k[0] = - tab_fic[0] * A ;
       k[1] =   tab_fic[1] * A ; 
    }
    cout << "value of the SIF k1:" << k[0] << "\n" ;
    cout << "value of the SIF k2:" << k[1] << "\n" ;
    
    // Compute relative error
    base_small_vector err_rel_coeff(4), k_exact(2), err_rel_fic(2) ;
    if (p.PARAM.int_value("SOL_REF") == 0.){
       if (p.PARAM.int_value("SING_BASE_TYPE") == 0){
          // 1.relative error on the coefs :
          for (int i=0 ; i< 4 ; ++i){
             if (p.exact_sol.U[i] != 0.)
	         err_rel_coeff[i] = gmm::abs(tab_fic[i] - p.exact_sol.U[i]) / gmm::abs(p.exact_sol.U[i]) ;
	     else
	         err_rel_coeff[i] = gmm::abs(tab_fic[i] - p.exact_sol.U[i]) ;
          }
          // 2. Compute exact fic
          k_exact[0] = p.exact_sol.U[2] * (3.*(p.nu - 1.)) + p.exact_sol.U[3] * (3. * p.nu + 5.) ;
          k_exact[0] *= - ( sqrt(2.) * E * p.epsilon ) / (4. * (1. - nu2) ) ;
          k_exact[1] = (p.exact_sol.U[0] + 3. * p.exact_sol.U[1]) ;
          k_exact[1] *= - E * p.epsilon * sqrt(2.) * (3. + p.nu) / (4. * (1. + p.nu) * (1. + p.nu) ) ;
	  // 3. Only in the 4 dofs case : measure the error of 
	  //    realisation of the singularity
	  base_small_vector err_to_singul(2) ;
	  cout << "FIC : Default of realisation of the singularity \n" ;  
	  if (p.exact_sol.U[2] != 0. && p.exact_sol.U[3] != 0.){
	     err_to_singul[0] = gmm::abs(tab_fic[2] - tab_fic[3] * (p.nu + 7.) / (3. * (p.nu - 1.))) ;
	     err_to_singul[0] /= gmm::abs(p.exact_sol.U[2]) ;
	     cout << "fic1 : err_to_singul:" << err_to_singul[0] << "\n" ;
	  } 
	  else cout << "coeffs 3 or 4 to compare with 0 => no comparison possible \n" ; 
	  if (p.exact_sol.U[0] != 0. && p.exact_sol.U[1] != 0.){
	     err_to_singul[1] = gmm::abs(tab_fic[1] - tab_fic[0] * (3. * p.nu + 5.) / (3. * (p.nu - 1.))) ;
	     err_to_singul[1] /= gmm::abs(p.exact_sol.U[1]) ;
	     cout << "fic 2 : err_to_singul:" << err_to_singul[1] << "\n" ;
	  } 
	  else cout << "coeffs 1 or 2 to compare with 0 => no comparison possible \n" ; 
       }
       else {
          // 1.relative error on the coefs :
          if (p.exact_sol.U[3] != 0.)
	      err_rel_coeff[0] = gmm::abs(tab_fic[0] - p.exact_sol.U[3]) / gmm::abs(p.exact_sol.U[3]) ;
	  else
	      err_rel_coeff[0] = gmm::abs(tab_fic[0] - p.exact_sol.U[3]) ;
          if (p.exact_sol.U[0] != 0.)
	      err_rel_coeff[1] = gmm::abs(tab_fic[1] - p.exact_sol.U[0]) / gmm::abs(p.exact_sol.U[0]) ;
	  else
	      err_rel_coeff[1] = gmm::abs(tab_fic[1] - p.exact_sol.U[0]) ;
          // 2. Compute exact fic
	  scalar_type A =  sqrt(2) * p.epsilon * E * (p.nu + 3.) / (1. - nu2) ;
          k_exact[0] = - p.exact_sol.U[3] * A ;
          k_exact[1] = p.exact_sol.U[0] * A ; 
       }
       cout << "Exact value of SIF : \n" ;
	  for (int i = 0; i<2 ; ++i){
	     cout << "k" << i+1 << "_exact = " << k_exact[i] << "\n" ; 
       } 
       // 2. Compute relative error on the fic
       for (int i=0 ; i< 2 ; ++i){
          err_rel_fic[i] = gmm::abs( k[i] - k_exact[i] ) ;
          if (k_exact[i] != 0.)
              err_rel_fic[i] /= gmm::abs(k_exact[i]) ;
       }
       
       cout << "FIC : Relative errors on the coeffs :\n" ;
       int nb_coeffs ;
       (p.PARAM.int_value("SING_BASE_TYPE") == 0) ? nb_coeffs = 4 : nb_coeffs = 2 ;
       for (int i=0 ; i< nb_coeffs ; ++i)
           cout << "err_rel_coeff:" << err_rel_coeff[i] << "\n" ;
       cout << "FIC : Relative errors on the fics :\n" ;
       for (int i=0 ; i < 2 ; ++i){
           cout << "err_rel_fic" << i + 1 << ":" << err_rel_fic[i] << "\n" ;
	   if (k_exact[i] == 0.){
	      cout << "Beware, absolute error is print, because exact value is 0.\n" ;
	   }
       }
    }



      if (p.PARAM.int_value("SOL_REF") == 1){
        scalar_type k1_exact = 3. * sqrt(p.PARAM.real_value("CRACK_SEMI_LENGTH")) / (2. * p.epsilon * p.epsilon) ;
        cout << "exact SIF k1:" << k1_exact << "\n" ;
        cout << "(To multipliy by Mo, which is usually set to 1).\n" ; 
        cout << "err_rel_fic1:" << gmm::abs(k[0] - k1_exact) / gmm::abs(k1_exact) << "\n" ;
        cout << "err_rel_fic2:" << gmm::abs(k[1])  << "(->0)\n" ;
      }

    } // End printing SIFs

    // --------------------------------------------------
    // Post-treatment for pretty printing of the solution
    // --------------------------------------------------

    int VTK_EXPORT = int(p.PARAM.int_value("VTK_EXPORT"));
    int MATLAB_EXPORT = int(p.PARAM.int_value("MATLAB_EXPORT"));
    int DX_EXPORT = int(p.PARAM.int_value("DX_EXPORT"));
	
    if (VTK_EXPORT || MATLAB_EXPORT || DX_EXPORT){
    // Post-traitement pour l'affichage, le plus simple :
    getfem::mesh mcut;
    p.mls.global_cut_mesh(mcut);
    unsigned Q = p.mf_u().get_qdim();
    assert( Q == 1 ) ;
    getfem::mesh_fem mf(mcut, dim_type(Q));
    mf.set_classical_discontinuous_finite_element(3, 0.001);
    plain_vector V(mf.nb_dof()) ;
    getfem::interpolation(p.mf_u(), mf, U, V);

    size_type N = mcut.dim() ;
    getfem::mesh_fem mf_grad(mcut, dim_type(Q));
    mf_grad.set_classical_discontinuous_finite_element(2, 0.001);
        // mf.set_finite_element
    //	(getfem::fem_descriptor("FEM_PK_DISCONTINUOUS(2, 2, 0.0001)"));
    plain_vector GRAD(N*mf_grad.nb_dof()) ;
/*    if (p.PARAM.int_value("TRY") == 1 ){ 
    getfem::compute_gradient(mf, mf_grad, U, GRAD) ;
    }
    if (p.PARAM.int_value("TRY") == 1 ){ 
    getfem::compute_gradient(mf, mf_grad, U, GRAD) ;
    }  */  


//     // Printing the normal rotation (= derivative of the vertical displacement)
//     if (p.PARAM.int_value("SHOW_DX") == 1 ){  // au lieu d'afficher la solution, on affiche sa derivee en y
//        getfem::compute_gradient(mf, mf_grad, V, GRAD) ;
//        gmm::resize(V, mf_grad.nb_dof()) ;
//        for (unsigned i=0 ; i< mf_grad.nb_dof(); i++){
//            V[i] = GRAD[2*i+1] ;
//        }
//       cout << "export to " << p.datafilename + ".vtk" << "..\n";
//       getfem::vtk_export exp(p.datafilename + ".vtk",
// 			     p.PARAM.int_value("VTK_EXPORT")==1);
//       exp.exporting(mf_grad); 
//       exp.write_point_data(mf_grad, V, "gradient_second_coordinate");
//       cout << "export done, you can view the data file with (for example)\n"
// 	   << "mayavi -d " << p.datafilename
// 	   << ".vtk -m BandedSurfaceMap -m Outline -f WarpScalar\n";
//     }
//     else{
//       cout << "export to " << p.datafilename + ".vtk" << "..\n";
//       getfem::vtk_export exp(p.datafilename + ".vtk",
// 			     p.PARAM.int_value("VTK_EXPORT")==1);
//       exp.exporting(p.mf_u()); 
//       exp.write_point_data(mf, V, "vertical_displacement");
//       cout << "export done, you can view the data file with (for example)\n"
// 	   << "mayavi -d " << p.datafilename
// 	   << ".vtk -m BandedSurfaceMap -m Outline -f WarpScalar\n";
//     }

    cout << "first part OK \n" ;

    getfem::stored_mesh_slice sl;
    getfem::mesh mcut_refined;

    unsigned NX = unsigned(p.PARAM.int_value("NX")), nn;
    if (NX < 6) nn = 12;
    else if (NX < 12) nn = 8;
    else if (NX < 30) nn = 3;
    else nn = 1;

    /* choose an adequate slice refinement based on the distance to the crack tip */
    std::vector<bgeot::short_type> nrefine(mcut.convex_index().last_true()+1);
    for (dal::bv_visitor cv(mcut.convex_index()); !cv.finished(); ++cv) {
      scalar_type dmin=0, d;
      base_node Pmin,P;
      for (unsigned i=0; i < mcut.nb_points_of_convex(cv); ++i) {
	P = mcut.points_of_convex(cv)[i];
	//d = gmm::vect_norm2(ls_function(P));
	d = gmm::vect_norm2(P) ; 
	if (d < dmin || i == 0) { dmin = d; Pmin = P; }
      }

      if (dmin < 1e-5)
	nrefine[cv] = short_type(nn*2);
      else if (dmin < .1) 
	nrefine[cv] = short_type(nn*2);
      else nrefine[cv] = short_type(nn);
      /*if (dmin < .01)
	cout << "cv: "<< cv << ", dmin = " << dmin << "Pmin=" << Pmin << " " << nrefine[cv] << "\n";*/
    }

    cout << "second part OK \n" ;

    getfem::mesh_slicer slicer(mcut);
    getfem::slicer_build_mesh bmesh(mcut_refined);
    slicer.push_back_action(bmesh);
    slicer.exec(nrefine, getfem::mesh_region::all_convexes());
     
    /*
      sl.build(mcut, 
      getfem::slicer_build_mesh(mcut_refined), nrefine);*/

    getfem::mesh_im mim_refined(mcut_refined); 
    mim_refined.set_integration_method(getfem::int_method_descriptor
				       ("IM_TRIANGLE(6)"));

    getfem::mesh_fem mf_refined(mcut_refined, dim_type(Q));
    mf_refined.set_classical_discontinuous_finite_element(2, 0.001);
    plain_vector W(mf_refined.nb_dof());

    cout << "OK before last interpolation \n" ;
    getfem::interpolation(mf, mf_refined, V, W);
    
//     getfem::mesh_fem mf_grad(mcut, 2) ;
//     plain_vector VGRAD(mf_grad.nb_dof() ) ;
//     getfem::compute_gradient(p.mls, p.mf_u(), mf_grad_u, U, GRAD) ;
//     getfem::mesh_fem mf_grad_refined(mcut_refined, 2) ;
//     plain_vector WGRAD(mf_grad_refined.nb_dof()) ;
//     getfem::interpolation(mf_grad, mf_grad_refined, GRAD, WGRAD) ;

      if (VTK_EXPORT){
      cout << "exporting solution to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk", false);
      exp.exporting(mf_refined); 
      exp.write_point_data(mf_refined, W, "vertical_displacement");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d " << p.datafilename  << ".vtk -f "
	"WarpScalar -m BandedSurfaceMap -m Outline\n";

      p.exact_sol.mf.set_qdim(dim_type(Q));
      assert(p.exact_sol.mf.nb_dof() == p.exact_sol.U.size());   
      plain_vector EXACT(mf_refined.nb_dof());
      getfem::interpolation(p.exact_sol.mf, mf_refined, 
			    p.exact_sol.U, EXACT);
	
      plain_vector DIFF(EXACT); gmm::add(gmm::scaled(W,-1),DIFF);
      if ((VTK_EXPORT & 2)) {
	cout << "exporting exact solution to VTK\n";
	getfem::vtk_export exp2(p.datafilename + "_exact.vtk");
	exp2.exporting(mf_refined);
	exp2.write_point_data(mf_refined, EXACT, "reference solution");
	cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d " << p.datafilename << "_exact.vtk -f "
	"WarpScalar -m BandedSurfaceMap -m Outline\n";
      }
      if ((VTK_EXPORT & 4)) {
	cout << "exporting difference with exact solution to VTK\n";
	getfem::vtk_export exp2(p.datafilename + "_diff.vtk");
	exp2.exporting(mf_refined);
	exp2.write_point_data(mf_refined, DIFF, "difference solution");
	cout << "export done, you can view the data file with (for example)\n"
	  "mayavi -d " << p.datafilename  << "_diff.vtk -f "
	"WarpScalar -m BandedSurfaceMap -m Outline\n";
      }
      if ((VTK_EXPORT == 5)) {
	cout << "exporting Von Mises \n";
	getfem::vtk_export exp3(p.datafilename + "_VM.vtk");
	exp3.exporting(mf_refined);
	plain_vector VM(mf_refined.nb_dof()) ;
	calcul_von_mises(p.mf_u(), U, mf_refined, VM, p.D) ;
	exp3.write_point_data(mf_refined, W, "difference solution");
	exp3.write_point_data(mf_refined, VM, "Von Mises");
	cout << "export done, you can view the data file with (for example)\n"
	  "mayavi -d " << p.datafilename  << "_diff.vtk -f "
	"WarpScalar -m BandedSurfaceMap -m Outline\n";
      }
      }
    
    if (MATLAB_EXPORT) {
      cout << "exporting solution to " << p.datafilename + ".mf" << " and " << p.datafilename << ".U\n";
      mf_refined.write_to_file(p.datafilename + ".mf", true);
      gmm::vecsave(p.datafilename + ".U", W);
      p.exact_sol.mf.set_qdim(dim_type(Q));
      assert(p.exact_sol.mf.nb_dof() == p.exact_sol.U.size());   // ??
      plain_vector EXACT(mf_refined.nb_dof());
      getfem::interpolation(p.exact_sol.mf, mf_refined, 
			    p.exact_sol.U, EXACT);
      gmm::vecsave(p.datafilename + ".EXACT", EXACT);

      cout << "exporting original mesh to " << p.datafilename + "_base.mesh\n";
      p.mesh.write_to_file(p.datafilename + "_base.mesh");

      cout << "matlab export done, you can view the results with\n";
      cout << "mf=gfMeshFem('load', 'bilaplacian.mf'); U=load('bilaplacian.U')'; "
	"EXACT=load('bilaplacian.EXACT')'; m0=gfMesh('load','bilaplacian_base.mesh');\n";
      cout << "gf_plot(mf, U-EXACT, 'refine',1); hold on; gf_plot_mesh(m0); hold off; colorbar;\n";
    }

    if (DX_EXPORT) {
    cout << "export solution to " << p.datafilename + ".dx" << "..\n";
    getfem::dx_export exp(p.datafilename + ".dx",
			   p.PARAM.int_value("DX_EXPORT")==1);
    exp.exporting(mf_refined); 
    exp.write_point_data(mf_refined, W, "vertical_displacement");
    cout << "export done for open dx \n" ;
    }

    
    
    
    
    /*************************************/
    // PRINTING IN 3D 
    /************************************/
    

    if (DX_EXPORT > 29){
    mf.clear() ;
    getfem::mesh_fem mf(mcut_refined, 1);
    mf.set_classical_discontinuous_finite_element(3, 0.001);
    
    
    plain_vector V3(mf.nb_dof()), VT(mf.nb_dof()*2), VTHETA(mf.nb_dof()*2); 
    if ( DX_EXPORT == 31){  // exporting the "finite element" solution on lagrange discont.
	// getfem::interpolation( "mesh_fem of exact sol", mf, "VECTOR of exact sol", VT); -> for membrane
	getfem::interpolation(mf_refined, mf, W, V3);
	plain_vector GRAD_U(mf_refined.nb_dof() * 2) ;
	getfem::compute_gradient(mf_refined, mf_refined, W, GRAD_U);
	getfem::interpolation(mf_refined, mf, GRAD_U, VTHETA);	
// 	getfem::interpolation(p.mf_u(), mf, U, V3);
// 	plain_vector GRAD_U(p.mf_rhs.nb_dof() * 2) ;
// 	getfem::compute_gradient(p.mf_u(), p.mf_rhs, U, GRAD_U);
// 	getfem::interpolation(p.mf_rhs, mf, GRAD_U, VTHETA);
	}
    else if (DX_EXPORT == 32) { // exporting the "exact" solution
        //p.exact_sol.mf_ut.set_qdim(2) ;  -> to modifify for taking account  of the membrane displacement
        //p.exact_sol.mf.set_qdim(2) ;
	//getfem::interpolation(p.exact_sol.mf_ut,    mf, p.exact_sol.UT, VT);
	getfem::interpolation(p.exact_sol.mf,    mf, p.exact_sol.U, V3);
	getfem::interpolation(p.exact_sol.mf,    mf_refined, p.exact_sol.U, W);
	plain_vector GRAD_U(mf_refined.nb_dof() * 2) ; // compute exact gradient on P3 mesh_fem
	getfem::compute_gradient(mf_refined, mf_refined, W, GRAD_U);
	getfem::interpolation(mf_refined, mf, GRAD_U, VTHETA);
    }
    
    getfem::mesh m3d; getfem::extrude(mcut_refined,m3d,1);
    getfem::base_matrix trans(3,3); 
    trans(0,0) = trans(1,1) = 1; trans(2,2) = p.epsilon;
    m3d.transformation(trans);
    getfem::mesh_fem mf3d(m3d);
    mf3d.set_classical_discontinuous_finite_element(3, 0.001);
    
    
    gmm::resize(V, mf3d.nb_dof()*3);
    bgeot::kdtree tree; tree.reserve(mf.nb_dof());
    for (unsigned i=0; i < mf.nb_dof(); ++i)
      tree.add_point_with_id(mf.point_of_dof(i),i);
    bgeot::kdtree_tab_type pts;
    for (unsigned i=0; i < mf3d.nb_dof(); ++i) {
      base_node P = mf3d.point_of_dof(i);
      base_node P2d0(2), P2d1(2); 
      scalar_type EPS = 1e-6;
      P2d0[0] = P[0]-EPS; P2d0[1] = P[1]-EPS;
      P2d1[0] = P[0]+EPS; P2d1[1] = P[1]+EPS;
      tree.points_in_box(pts, P2d0, P2d1);
      //cout << "P = " << P << ", P2d0=" << P2d0 << " " << P2d1 << ", pts.size=" << pts.size() << "\n";
      assert(pts.size() == 1);
      size_type j=pts[0].i;
      scalar_type x3 = P[2];
      assert(finite(VTHETA[2*j]));
      //assert(finite(VT[2*j]));
      V[3*i+0] = - x3 * VTHETA[2*j+0]; //V[3*i+0] = VT[2*j+0] + x3 * VTHETA.at(2*j+0);
      V[3*i+1] = - x3 * VTHETA[2*j+1]; //V[3*i+1] = VT[2*j+1] + x3 * VTHETA[2*j+1];
      V[3*i+2] = V3[j];
      assert(finite(V[3*i]));
      assert(finite(V[3*i+1]));
      assert(finite(V[3*i+2]));
    }
    
    cout << "export to " << p.datafilename + "_3d.dx" << "..\n";
	getfem::dx_export exp(p.datafilename + "_3d.dx",
			       p.PARAM.int_value("DX_EXPORT")>30);
    getfem::stored_mesh_slice sl_tetra;
    getfem::mesh mtetra;
    sl_tetra.build(m3d, getfem::slicer_build_mesh(mtetra), 1);
    getfem::mesh_fem mftetra(mtetra,3); 
    mftetra.set_classical_discontinuous_finite_element(3, 0.001);
    
    plain_vector Vtetra(mftetra.nb_dof());
    mf3d.set_qdim(3);
    getfem::interpolation(mf3d, mftetra, V, Vtetra);
    /* evaluating Von Mises*/
    getfem::mesh_fem mf_vm(mtetra,1) ; 
    mf_vm.set_classical_discontinuous_finite_element(3, 0.001);
    plain_vector VM(mf_vm.nb_dof()) ;
    calcul_von_mises(mftetra, Vtetra, mf_vm, VM) ;   // to do: add a mu value in additional last argument
        
    exp.exporting(mftetra);
    exp.write_point_data(mftetra, Vtetra, "plate_displacement");
    size_type deb=1;
    if (deb==1){ 
        exp.exporting(mf_vm);
        exp.write_point_data(mf_vm, VM, "von_mises_stress");
    }
    }
    }
   //getchar(); 
   } GMM_STANDARD_CATCH_ERROR;
  return 0; 
}
























