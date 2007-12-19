#include "crack_bilaplacian.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_fourth_order.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_superlu.h"

void bilaplacian_crack_problem::move_nodes_close_to_crack(void){
        
        cout << "passage par move_nodes_close_to_crack \n" ;
	sparse_matrix M(mfls_u.nb_dof(), mfls_u.nb_dof());
	getfem::asm_mass_matrix(M, mim, mfls_u, mfls_u);
	dal::bit_vector bad_dofs ;
	
	// Selection des dofs pathologiques
	for (int i=0; i < mfls_u.nb_dof(); i++){
	    if (M(i,i) < PARAM.real_value("MOVE_NODES_PARAMETER")) {
	       cout << "dof : " << mfls_u.point_of_dof(i)[0] << " ; " << mfls_u.point_of_dof(i)[1] << " sélectionné " ;
	       cout << "(M(" << i << ", " << i << ") = " << M(i,i) << "\n" ;
	       bad_dofs.add(i) ;
	    } 
	}
	
	cout << "dofs pathologiques : " << bad_dofs << "\n" ;
	size_type ind_new;
	std::vector<bgeot::base_node> pts(mesh.structure_of_convex(mesh.convex_index().first_true())->nb_points()) ;
	
	dal::bit_vector cv_to_update ;
	
	cout << "on passe en revue chaque dof pathologique \n" ;
	for( dal::bv_visitor i(bad_dofs) ; !i.finished(); ++i){
	   cout << "dof pathologique indice " << i << " \n" ;
	   cv_to_update.clear() ;
	   // Pour chaque dof pathologique, on constitue la liste des elements
	   // qui possèdent ce dof.
	   cout << "dof : " << mfls_u.point_of_dof(i)[0] << " ; " << mfls_u.point_of_dof(i)[1] << " : \n" ;
	   for( dal::bv_visitor j(mesh.convex_index()) ; !j.finished(); ++j){
	      // sur chaque element du maillage, on passe en revue tous les dofs,
	      // de manière à tester la présence ou l'absence du dof pathologique.
	      for( size_type k = 0 ; k < mfls_u.nb_dof_of_element(j) ; ++k){
		 if ( mfls_u.ind_dof_of_element(j)[k] == i ){ // ??
		    cv_to_update.add(j) ;
		    cout << "  un element ajouté. \n" ;
		    break;
	         }
	      } 
	}
	
	ind_new = mesh.add_point(bgeot::base_node(mfls_u.point_of_dof(i)[0], 0.0)) ;
	
	// On passe en revue tous les elements a mettre a jour :  
	for( dal::bv_visitor l(cv_to_update) ; !l.finished() ; ++l){
	// on crée les points du nouvel élément
	for( size_type m = 0 ; m < (mesh.structure_of_convex(l)->nb_points() ) ; ++m){
		if( mesh.points()[ mesh.ind_points_of_convex(i)[m] ][0] == mfls_u.point_of_dof(i)[0] 
		 && mesh.points()[ mesh.ind_points_of_convex(i)[m] ][1] == mfls_u.point_of_dof(i)[1] )
		pts[m] = mfls_u.point_of_dof(i);
		else
		pts[m] = mesh.points()[m] ;
	}
	// avant de créer l'élément, on récupère la transfo géométrique
	std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
	MESH_TYPE = bgeot::name_of_geometric_trans
		(mesh.trans_of_convex(mesh.convex_index().first_true()));
	bgeot::pgeometric_trans pgt = bgeot::geometric_trans_descriptor(MESH_TYPE);
	// on crée finalement l'élément
	mesh.add_convex_by_points(pgt, pts.begin()) ;
	
	// on détruit l'ancien élément
	mesh.sup_convex(l) ;
	}
	
	} 
}


