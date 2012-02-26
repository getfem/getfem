/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2012-2012 Yves Renard, Julien Pommier.
 
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
#include <queue>
#include <queue>
#include "getfem/dal_singleton.h"
#include "getfem/getfem_mesh_fem.h"


namespace getfem {

#if GETFEM_PARA_LEVEL > 1

 struct fem_dof {
    size_type ind_node;
    pdof_description pnd;
    size_type part;
  };

/* Fonction de renumérotation des degrés de libertés */

/// Parallel Enumeration of dofs
void mesh_fem::enumerate_dof_para(void) const {

#if 0

  GMM_TRACE2("Enumeration dof para !!!!!!!!!!!!!!!!");
    GMM_ASSERT1(linked_mesh_ != 0, "Uninitialized mesh_fem");
    context_check();
    if (fe_convex.card() == 0) {
      dof_enumeration_made = true;
      nb_total_dof = 0;
      return;
    }

    dof_structure.clear();
    MPI_Status mstatus;
    fem_dof fd;

/* Récupération du nombre total de régions (procs)!!!!*/
    int num_rg, nb_rg;
    MPI_Comm_rank(MPI_COMM_WORLD, &num_rg);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_rg);
	
	
    cout<<"nb total region : "<<nb_rg<<" et nb_points = "<<linked_mesh().nb_points()<<endl;

/* Récupération du num de la région */
    //size_type num_rg;
    //num_rg = linked_mesh.mpi_region.id();
	
/* Création de la liste des ddl */
/// list_of_dof[1:nb_total_dof_mesh, 1:2]
    std::vector<size_type> list_of_dof_numeration;
    std::vector<size_type> list_of_dof_num_rg;
    dal::bit_vector enumeration_of_dof_made;

    list_of_dof_numeration.resize(100000);
    list_of_dof_num_rg.resize(100000);
    //    enumeration_of_dof_made(100000);

/* Création de la liste des ddl interface */
    std::vector<size_type> list_of_dof_linkable_index;
    std::vector<size_type> list_of_dof_linkable_to;

    list_of_dof_linkable_index.resize(10000);
    list_of_dof_linkable_to.resize(10000);

/* Création de la liste des ddl globaux */	
    std::vector<size_type> list_of_global_dof_index;
    std::vector<size_type> list_of_global_dof_local_index;

    list_of_global_dof_index.resize(10000);
    list_of_global_dof_local_index.resize(10000);

/* Initialisation de l'itérateur sur list_of_dof */
    size_type iter_dof = 0;

/* Construction de la liste des cv à la charge de chaque proc en fonction de la région */
/// list_of_cv[1:nb_cv_total_mesh, 1:2]
    std::vector<size_type> list_of_cv_num_rg;
    //list_of_cv_num_rg = size_type(0);
    std::vector<size_type> list_of_cv_first_index;

    std::vector<size_type> list_of_icv;
    std::vector<size_type> icv_in_list;
	
    list_of_cv_num_rg.resize(10000);
    list_of_cv_first_index.resize(10000);
    list_of_icv.resize(1000);
    icv_in_list.resize(1000);

/// En parallèle sur les régions on récupère les indices des cv de chaque proc
    //dal::bit_vector index_tab;
    const std::vector<size_type> &cmk = linked_mesh().cuthill_mckee_ordering();
    //index_tab = linked_mesh().region(num_rg).index();

    cout<<"cmk.size = "<<cmk.size()<<endl;
    cout<<"cmk[0] = "<<cmk[0]<<endl;
    cout<<"cml[10] = "<<cmk[10]<<endl;
    cout<<"nb_cv = "<<linked_mesh().convex_index().card()<<endl;

    size_type nb_cv = 0;
    std::vector<size_type> neighboors;
    bgeot::pgeotrans_precomp pgp = 0;
    base_node P;
    bgeot::pgeotrans_precomp pgpj = 0;
    base_node Pj;
	
// Boucle i pour remplir la liste des cv:
 GMM_TRACE2("Initialisation of lists cv");
 // for(dal::bv_visitor i(index_tab); !i.finished(); ++i)
 cout<<"me = "<<num_rg<<endl;

 bool entre = false;

 cout<<"bool avant : "<<entre<<endl;

 for(size_type i = cmk[0]; i<cmk.size(); i++)
      {
	size_type icv = cmk[i];
	if(linked_mesh().region(num_rg).is_in(icv))
	  {
	    GMM_TRACE2("ICI 0");
	    // cout<<"i = "<<i<<endl;
	    //cout<<"me "<<num_rg<<" et icv = "<<icv<<endl;
	    list_of_cv_num_rg[nb_cv] = num_rg;
	    GMM_TRACE2("ICI 1");
	    list_of_cv_first_index[nb_cv] = iter_dof;
	    GMM_TRACE2("ICI 2");
	    list_of_icv[nb_cv] = icv;
	    icv_in_list[icv] = nb_cv;
		
	    pfem pf = fem_of_element(icv);
	    size_type nbd = pf->nb_dof(icv);
	    iter_dof += nbd;
	    nb_cv += 1;
	    //cout<<"la!!!!!!"<<nb_cv<<endl;
	    entre = true;
	  }
	else
	  {
	    list_of_cv_num_rg[nb_cv] = 0;
	    list_of_cv_first_index[nb_cv] = 0;
	    list_of_icv[nb_cv]=0;
	    icv_in_list[icv]=0;
	    
	    pfem pf = fem_of_element(icv);
	    size_type nbd = pf->nb_dof(icv);
	    iter_dof += nbd;
	    nb_cv += 1;
	    //cout<<"ou la !!!!!!!!!"<<nb_cv<<endl;
	    
	  }
	//cout<<"me = "<<num_rg<<"iteration i : "<<i<<endl;
    }

 cout<<"bool : "<<entre<<endl;
	
 //int nb_cv_tot;
// Mise en commun par échange MPI_AllReduce du nombre totale de cv sur le mesh
    //GMM_TRACE2("Echange MPI 1");
    //   MPI_Allreduce(&nb_cv, &nb_cv_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);

    cout<<"nb_cv_tot = "<<nb_cv<<endl;




    std::vector<size_type> list_of_cv_num_rg_Recv;
    std::vector<size_type> list_of_cv_first_index_Recv;
    std::vector<size_type> list_of_icv_Recv;
    std::vector<size_type> icv_in_list_Recv;
    
    cout<<"size(cv_num_rg) = "<<list_of_cv_num_rg.size()<<endl;

    list_of_cv_num_rg.resize(nb_cv);
    list_of_cv_num_rg_Recv.resize(nb_cv);
    list_of_cv_first_index_Recv.resize(nb_cv);
    list_of_icv_Recv.resize(nb_cv);
    icv_in_list_Recv.resize(nb_cv);

    MPI_Barrier(MPI_COMM_WORLD);

// Mise en commun par échange MPI_AllReduce de la liste list_of_cv_num_rg
    GMM_TRACE2("Echange MPI 2");
    MPI_Allreduce (&list_of_cv_num_rg[0], &list_of_cv_num_rg_Recv, nb_cv, 
                   MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

    cout<<"num_rg[10] = "<<list_of_cv_num_rg[10]<<endl;
    cout<<"num_rg_Recv[10] = "<<list_of_cv_num_rg_Recv[10]<<endl;

    GMM_TRACE2("Echange 3");
// Mise en commun par échange MPI_AllReduce de la liste list_of_cv_num_rg
    MPI_Allreduce (&list_of_cv_first_index[0], &list_of_cv_first_index_Recv, nb_cv, 
                   MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

    GMM_TRACE2("Echange 4");
    MPI_Allreduce (&list_of_icv[0], &list_of_icv_Recv, nb_cv, 
		   MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

    GMM_TRACE2("Echange 5");

    MPI_Allreduce (&icv_in_list[0], &icv_in_list_Recv, nb_cv, 
		   MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

/* Construction de la liste des ddl à la charge de chaque proc */
    size_type nb_dof_rg = 0;
    size_type nb_dof_tot;
    //size_type nb_dof_inter = 0;
    size_type nb_global_dof = 0;
    size_type nb_global_dof_tot;
    size_type nb_dof_linkable = 0;
	
// Pour chaque cv dont ce proc a la charge : 
 GMM_TRACE2("Attributions des ddl aux rg");
    for(size_type cv = 0; cv < list_of_cv_num_rg_Recv.size(); ++cv)
    {
      size_type icv = list_of_icv_Recv[cv];
	if (list_of_cv_num_rg_Recv[cv] == num_rg)
	{
	    pfem pf = fem_of_element(icv);
	    size_type nbd = pf->nb_dof(icv);
	    nb_dof_rg += nbd;
	    pdof_description andof = global_dof(pf->dim());
			
// 	    pour chaque ddl associé à ce cv :
	    for (size_type i = list_of_cv_first_index_Recv[cv]; 
		 i <= list_of_cv_first_index_Recv[cv] + nbd; i++)
	    {
		fd.pnd = pf->dof_types()[i];
		fd.part = get_dof_partition(icv);
				
//	 	Test si le ddl est raccordable
		if (dof_linkable(fd.pnd))
		{
		     size_type bool_rg = 0;
		     size_type bool_inter = 0;
		     P.resize(linked_mesh().dim()); 
		     pgp->transform(linked_mesh().points_of_convex(icv), i, P);
					
//		     Récupération des voisins qui possèdent ce même ddl :
		     neighboors = linked_mesh().convex_to_point(i);
		     for (size_type jcv = neighboors[0]; jcv < neighboors.size(); ++jcv)
		     {
//			 Si le voisin appartient à la même région (ie pas ddl interface)
			 if (list_of_cv_num_rg_Recv[icv_in_list_Recv[jcv]] == num_rg)
			 {
			      bool_rg++;
			 }
//			 Sinon si c'est un dof interface "et" qui doit être à la charge de cette région
			 else if (list_of_cv_num_rg_Recv[icv_in_list_Recv[jcv]] > num_rg)
			 {						
			      bool_inter++;
			 }
		     }
//		     Test si pas ddl interface
		     if (bool_rg==neighboors.size() || bool_inter+bool_rg == neighboors.size()) 
		       // ie tout les voisins raccordable sont dans cette même region
		     {
		       /*   for(size_type jcv = neighboors[0]; jcv < neighboors.size(); ++jcv)
			 {
			   list_of_dof_linkable_index[nb_dof_linkable] = list_of_cv_first_index_Recv[jcv]+j;
			   list_of_dof_linkable_to[nb_dof_linkable] = i;
			   nb_dof_linkable ++;
			 }
			 list_of_dof_num_rg[i] = num_rg;
		     }
		     else if ((bool_inter + bool_rg)==neighboors.size()) 
		       // ie tout les voisins raccordable doivent être associé à cette region
		       {*/
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
				 if (&P == &Pj)
				 {
				     list_of_dof_linkable_index[nb_dof_linkable] = 
				       list_of_cv_first_index_Recv[icv_in_list_Recv[jcv]]+j;
				     list_of_dof_linkable_to[nb_dof_linkable] = i;
				     nb_dof_linkable++;
				     list_of_dof_num_rg[list_of_cv_first_index_Recv
							[icv_in_list_Recv[jcv]]+j] = num_rg;
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

// Mise en commun par échange MPI_AllReduce du nombre totale de dof sur le mesh

	MPI_Allreduce (&nb_dof_rg, &nb_dof_tot, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
	size_type nbd_p;
	size_type numerot = 0;
	for (int p = 0; p < nb_rg; p++)
	{
	  // Si il a des proc plus petit que num_rg, il recoit le nb de dof des autres plus petit pour mettre à jour son indice de début de numérotation
	  if (p < num_rg)
	  {
	    MPI_Recv(&nbd_p, 1, MPI_UNSIGNED, p, 100, MPI_COMM_WORLD, &mstatus);
	    numerot += nbd_p;
	  }
	  // Sinon il envoi le nombre de dof qu'il a à sa charge au autre qui lui sont supérieur
	  else if (p > num_rg)
	  {
	    MPI_Send(&nb_dof_rg, 1, MPI_UNSIGNED, p, 100, MPI_COMM_WORLD);
	  }
	}
	    

	std::vector<size_type> list_of_dof_num_rg_Recv;
	list_of_dof_num_rg_Recv.resize(nb_dof_tot);

// Mise en commun par échange MPI_AllReduce de la liste list_of_cv_num_rg
	MPI_Allreduce(&list_of_dof_num_rg[0], &list_of_dof_num_rg_Recv, nb_dof_tot, 
		       MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);


	std::vector<size_type> list_of_global_dof_index_Recv;
	std::vector<size_type> list_of_global_dof_local_index_Recv;
	size_type nb_global_dof_Recv;

// Mise en commun des information sur les global_dof
	MPI_Allreduce (&nb_global_dof, &nb_global_dof_tot, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

	list_of_global_dof_index_Recv.resize(nb_global_dof_tot);
	list_of_global_dof_local_index_Recv.resize(nb_global_dof_tot);

	for (int p = 0; p<nb_rg; p++)
	{
	    MPI_Send (&nb_global_dof, 1, MPI_UNSIGNED, p, 200, MPI_COMM_WORLD);
	    MPI_Recv (&nb_global_dof_Recv, 1, MPI_UNSIGNED, p, 200, MPI_COMM_WORLD, &mstatus);

	    MPI_Send(&list_of_global_dof_index[0], nb_global_dof, MPI_UNSIGNED, p, 300, MPI_COMM_WORLD);

	    MPI_Recv (&list_of_global_dof_index_Recv[0], nb_global_dof_Recv, 
		      MPI_UNSIGNED, p, 300, MPI_COMM_WORLD, &mstatus);

	    MPI_Send(&list_of_global_dof_local_index[0], nb_global_dof, MPI_UNSIGNED, p, 400, MPI_COMM_WORLD);

	    MPI_Recv (&list_of_global_dof_local_index_Recv[0], nb_global_dof_Recv, 
		      MPI_UNSIGNED, p, 400, MPI_COMM_WORLD, &mstatus);
	}



/* 	Numérotation des ddl en charge de chaque processeur */
// 	Pour chaque cv dont ce proc a la charge : 
	GMM_TRACE2("Numerotation des ddl");
	size_type ind_linkable = 0;
 for(size_type icv = 0; icv < list_of_cv_num_rg_Recv.size(); ++icv)
    {
       pfem pf = fem_of_element(icv);
       size_type nbd = pf->nb_dof(icv);
			
//     pour chaque ddl associé à ce cv :
       for (size_type i = list_of_cv_first_index_Recv[icv]; 
	    i <= list_of_cv_first_index_Recv[icv] + nbd; i++)
	 {
	    if (list_of_dof_num_rg_Recv[i] == num_rg && !enumeration_of_dof_made[i])
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
			   enumeration_of_dof_made[list_of_dof_linkable_index[ind_linkable]] = true;
			   ind_linkable++;
		      }
		      list_of_dof_numeration[i] = numerot;
		      enumeration_of_dof_made[i] = true;
		      numerot += Qdim / pf->target_dim();
		  } // Fin Si
	     } // Fin boucle
	} // Fin Boucle
    }// Fin Boucle

// Traitement des ddl globaux
// Boucle sur list_of_global_dof_in_charge
/*	if(num_rg == 0)// temporairement
	  {
	for (size_type i=0; i < list_of_global_dof_index_Recv.size(); i++)
	  {
	    pfem pf = fem_of_element(icv);

	    if(!enumeration_of_dof_made[i])
	    {
///	         Récupère les indices ayant le même num_global
	         for (size_type j = i; j < list_of_global_dof_index_Recv.size(); j++)
	         {
///	             Numérotation
	             if (list_of_global_dof_index_Recv[j] == list_of_global_dof_index_Recv[i] 
			 && !enumeration_of_dof_made[j])
		     {
	                 list_of_dof_numeration[list_of_global_dof_local_index_Recv[j]] = numerot;
	                 enumeration_of_dof_made[list_of_gloabl_dof_local_index_Recv[j]] = true;
	             }
	         }
		 list_of_dof_numeration[list_of_global_dof_local_index[i]] = numerot;
		 enumeration_of_dof_made[list_of_global_dof_local_index[i]] = true;
		 numerot += Qdim / pf->target_dim();
	    }
	  }
	  }*/
// Fin boucle		

// Mise en commun de list_of_dof par échange avec MPI_AllReduce
 std::vector<size_type> list_of_dof_numeration_Recv;
 list_of_dof_numeration_Recv.resize(nb_dof_tot);

	MPI_Allreduce(&list_of_dof_numeration[0], &list_of_dof_numeration_Recv, numerot,
		      MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

// Envoi de la structure numérotée
 GMM_TRACE2("Envoi de la numerotation");
	std::vector<size_type> tab;
	size_type ind_tab = 0;
	for(size_type icv = 0; icv < list_of_cv_num_rg.size(); ++icv)
	{
	   pfem pf = fem_of_element(icv);
	    if (list_of_cv_num_rg[icv] == num_rg)
	    {
	        size_type nbd = pf->nb_dof(icv);
		tab.resize(nbd);
		for (size_type i = list_of_cv_first_index_Recv[icv]; i < list_of_cv_first_index_Recv[icv] + nbd; i++)
		{
		  tab[ind_tab] = list_of_dof_numeration[i];
		  ind_tab++;
		}
	    }
            dof_structure.add_convex_noverif(pf->structure(icv), tab.begin(), icv);
	}
	
	nb_total_dof = nb_dof_tot;


#endif
}

#endif


} // end of getfem namespace
