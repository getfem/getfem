#include <queue>
#include "getfem/dal_singleton.h"
#include "getfem/getfem_mesh_fem.h"

/* Fonction de renumérotation des degrés de libertés */

/// Parallel Enumeration of dofs
void mesh_fem::enumerate_dof_para(void) const {

    GMM_ASSERT1(linked_mesh_ != 0, "Uninitialized mesh_fem");
    context_check();
    if (fe_convex.card() == 0) {
      dof_enumeration_made = true;
      nb_total_dof = 0;
      return;
    }

    dof_structure.clear();
    MPI_STATUS mstatus;

/* Récupération du nombre total de régions (procs)!!!!*/	
/// ??????????
/// ??????????
	
	
/* Récupération du num de la région */
    size_type num_rg;
    num_rg = linked_mesh.mpi_region.id();
	
/* Création de la liste des ddl */
/// list_of_dof[1:nb_total_dof_mesh, 1:2]
    std::vector<size_type> list_of_dof_numeration;
    std::vector<size_type> list_of_dof_num_rg;
    dal::bit_vector enumeration_of_dof_made;

/* Création de la liste des ddl interface */
    std::vector<size_type> list_of_dof_linkable_index;
    std::vector<size_type> list_of_dof_linkable_to;

/* Création de la liste des ddl globaux */	
    std::vector<size_type> list_of_global_dof_index;
    std::vector<size_type> list_of_global_dof_local_index;

/* Initialisation de l'itérateur sur list_of_dof */
    size_type iter_dof = 0;

/* Construction de la liste des cv à la charge de chaque proc en fonction de la région */
/// list_of_cv[1:nb_cv_total_mesh, 1:2]
    std::vector<size_type> list_of_cv_num_rg = 0;
    std::vector<size_type> list_of_cv_first_index;
	
/// En parallèle sur les régions on récupère les indices des cv de chaque proc
    dal::bit_vector index_tab;
    index_tab = mesh_region.index();
	
    size_type nb_cv = 0;
    std::vector<size_type> neighboors;
    bgeot::pgeotrans_precomp pgp = 0;
    base_node P;
    bgeot::pgeotrans_precomp pgpj = 0;
    base_node Pj;
	
// Boucle i pour remplir la liste des cv:
    for(dal::bv_visitor icv(index_tab); !icv.finished(); ++icv)
    {
        list_of_cv_num_rg[icv] = num_rg;
	list_of_cv_fisrt_index[icv] = iter_dof;
		
	pfem pf = fem_of_element(icv);
	size_type nbd = pf->nb_dof(icv);
	iter_dof += nbd;
	nb_cv += 1;
    }
	
// Mise en commun par échange MPI_AllReduce du nombre totale de cv sur le mesh
    MPI_Allreduce (nb_cv, nb_cv_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);


// Mise en commun par échange MPI_AllReduce de la liste list_of_cv_num_rg
    MPI_Allreduce (list_of_cv_num_rg[0], list_of_cv_num_rg, nb_cv_tot, 
                   MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);


/* Construction de la liste des ddl à la charge de chaque proc */
    size_type nb_dof = 0;
    size_type nb_dof_tot;
    size_type nb_dof_inter = 0;
    size_type nb_global_dof = 0;
    size_type nb_global_dof_tot;
	
// Pour chaque cv dont ce proc a la charge : 
    for(size_type icv = 0; icv < list_of_cv_num_rg.size(); ++icv)
    {
	if (list_of_cv_num_rg[icv] == num_rg)
	{
	    pfem pf = fem_of_element(icv);
	    size_type nbd = pf->nb_dof(icv);
	    nb_dof += nbd;
	    pdof_description andof = global_dof(pf->dim());
			
// 	    pour chaque ddl associé à ce cv :
	    for (size_type i = list_of_cv_first_index[icv], i <= list_of_cv_fisrt_index[icv] + nbd; i++)
	    {
		fd.pnd = pf->dof_types()[i];
		fd.part = get_dof_partition(cv);
				
//	 	Test si le ddl est raccordable
		if (dof_linkable(fd.pnd))
		{
		     size_type bool_rg = 0;
		     size_type bool_inter = 0;
		     P.resize(linked_mesh().dim()); 
		     pgp->transform(linked_mesh().points_of_convex(icv), i, P);
					
//		     Récupération des voisins qui possèdent ce même ddl :
		     neighboors = linked_mesh.convex_to_point(i);
		     for (size_type jcv = neighboors[0]; jcv < neighboors.size(); ++jcv)
		     {
//			 Si le voisin appartient à la même région (ie pas ddl interface)
			 if (list_of_cv_num_rg[jcv] == num_rg)
			 {
			      bool_rg++;
			 }
//			 Sinon si c'est un dof interface "et" qui doit être à la charge de cette région
			 else if (list_of_cv_num_rg[jcv] > num_rg)
			 {						
			      bool_inter++;
			 }
		     }
//		     Test si pas ddl interface
		     if (bool_rg==neighboors.size()) // ie tout les voisins raccordable sont dans cette même region
		     {
			 list_of_dof_linkable_index[nb_dof_linkable] = list_of_cv_first_index[jcv]+j;
			 list_of_dof_linkable_to[nb_dof_linkable] = i;
			 nb_dof_linkable ++;
			 list_of_dof_num_rg[i] = num_rg;
		     }
		     else if ((bool_inter + bool_rg)==neighboors.size()) // ie tout les voisins raccordable doivent être associé à cette region
		     {
		         list_of_dof_num_rg[i] = num_rg;
			 for (size_type jcv = neighboors[0]; jcv < neighboors.size(); ++jcv)
			 {
///			     on associe le ddl correspondant au proc
			     pfem pfj = fem_of_element(jcv);
			     size_type nbdj = pfj->nb_dof(jcv);
			     for(size_type j = 0; j < nbdj; j++)
			     {
				 Pj.resize(linked_mesh().dim()); 
				 pgpj->transform(linked_mesh().points_of_convex(jcv), j, Pj);
				 if (P == Pj)
				 {
				     list_of_dof_linkable_index[nb_dof_linkable] = list_of_cv_first_index[jcv]+j;
				     list_of_dof_linkable_to[nb_dof_linkable] = i;
				     nb_dof_linkable++;
				     list_of_dof_num_rg[list_of_cv_first_index[jcv]+j] = num_rg;
				  }
			      }
			   }
		       }
					
		   }
//		   Si ddl global
		   else if(fd.pnd == andof)
		   {
///		       on recupère son numéro global : encountered_global_dof[1:3, num_region] = [num_global, icv, i]
		       size_type num = pf->index_of_global_dof(icv, i);
		       list_of_global_dof_index[nb_global_dof] = num;
		       list_of_global_dof_local_index [nb_global_dof] = i;
		       list_of_dof_num_rg[i] = num_rg;
		       nb_global_dof++;
		   }
//		   Si le ddl est non raccordable
		   else
		   {
///		       on associe ce ddl au proc => list_of_dof[i, 1] = num_region (ou qqch de remarquable!!!!)
		       list_of_dof_num_rg[i] = num_rg;
		   } //end if
	        } // end for
	    } // end if
	} // end for



// Mise en commun par échange MPI_AllReduce de la liste list_of _dof

// Mise en commun par échange MPI_AllReduce du nombre totale de cv sur le mesh
	MPI_Allreduce (nb_dof, nb_dof_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
	size_type nbd_p;
	for (size_type p = 0; p < nb_rg; p++)
	{
	  if (p < num_rg)
	  {
	    MPI_Recv(nbd_p, 1, MPI_INTEGER, p, 100, MPI_COMM_WORLD, mstatus);
	    numerot += nbd_p;
	  }
	  else if (p > num_rg)
	  {
	    MPI_Send(nb_dof, 1, MPI_INTEGER, p, 100, MPI_COMM_WORLD);
	  }
	}
	    


// Mise en commun par échange MPI_AllReduce de la liste list_of_cv_num_rg
	MPI_Allreduce (list_of_dof_num_rg[0], list_of_dof_num_rg, nb_dof_tot, 
                   MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);


// Mise en commun des information sur les global_dof
	MPI_Allreduce (nb_global_dof, nb_global_dof_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce (list_of_global_dof_index[0], list_of_global_dof_index[num_rg*1+nb_global_dof], nb_global_dof_tot,
		       MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce (list_of_global_dof_local_index[0], list_of_global_dof_local_index[num_rg*1+nb_global_dof], nb_global_dof_tot,
		       MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);



/* 	Numérotation des ddl en charge de chaque processeur */
// 	Pour chaque cv dont ce proc a la charge : 
	size_type numerot = 0;
	size_type ind_linkable = 0;
	
	for(size_type i = 0; i < list_of_dof_num_rg.size(); ++i)
	{
	    pfem pf = fem_of_element(icv);

	    if (list_of_dof_num_rg[i] == num_region && !enumeration_of_dof_made[i])
	    {
		 if (!dof_linkable(fd.pnd))
		 {
		      list_of_dof_numeration[i] = numerot;
		      numerot += Qdim / pf->target_dim();
		 }
//		 Test si ddl raccordable
		 else if (dof_linkable(fd.pnd))
		 {
///		      Recherche des cv ayant ce même point et dont le proc à la charge
		      while (list_of_dof_linkable_to[ind_linkable] == i)
		      {
///			   Récupération des indices correspondants dans list_of_dof et Numérotation dans list_of_dof(indices)
			   list_of_dof_numeration[list_of_dof_linkable_index[ind_linkable]] = numerot;
			   enumration_of_dof_made[list_of_dof_linkable_index[ind_linkable]] = true;
			   ind_linkable++;
		      }
		      list_of_dof_numeration[i] = numerot;
		      enumration_of_dof_made[i] = true;
		      numerot += Qdim / pf->target_dim();
		  } // Fin Si
	     } // Fin boucle
	} // Fin Boucle

// Traitement des ddl globaux
// Boucle sur list_of_global_dof_in_charge
	for (size_type i=0; i < list_of_global_dof_index.size(); i++)
	  {
	    pfem pf = fem_of_element(icv);

	    if(!enumeration_of_dof_made[i])
	    {
///	         Récupère les indices ayant le même num_global
	         for (size_type j = i; j < list_of_global_dof_index.size(); j++)
	         {
///	             Numérotation
	             if (list_of_global_dof_index[j] == list_of_global_dof_index[i] && !enumeration_of_dof_made[j])
		     {
	                 list_of_dof_numeration[j] = numerot;
	                 enumeration_of_dof_made[j] = true;
	             }
	         }
		 list_of_dof_numeration[i] = numerot;
		 enumeration_of_dof_made[i] = true;
		 numerot += Qdim / pf->target_dim();
	    }
	  }
// Fin boucle		

// Mise en commun de list_of_dof par échange avec MPI_AllReduce
	MPI_Allreduce(list_of_dof_numeration[0], list_of_dof_numeration, numerot
		      MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);

// Envoi de la structure numérotée
	size_type ind_tab = 0;
	for(size_type icv = 0; icv < list_of_cv_num_rg.size(); ++icv)
	{
	    if (list_of_cv_num_rg[icv] == num_rg)
	    {
	        pfem pf = fem_of_element(icv);
	        size_type nbd = pf->nb_dof(icv);
		tab.resize(nbd);
		for (size_type i = list_of_cv_first_index[icv]; i < list_of_cv_first_index[icv] + nbd; i++)
		{
		  tab[ind_tab] = list_of_dof_numeration[i];
		  ind_tab++;
		}
	    }
            dof_structure.add_convex_noverif(pf->structure(icv), tab.begin(), icv);
	}
	
	nb_total_dof = nb_dof_tot;
}
