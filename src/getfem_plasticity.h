/************************************************************************* */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  plasticity.h : perfect plasticity problem for isotropic      */
/*                           materials                                     */
/*                                                                         */
/* Date : June 10, 2004.                                                   */
/* Author : Marc ODUNLAMI                                                  */
/*                                                                         */
/* *********************************************************************** */

/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2004 Yves Renard, Julien Pommier.                    */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */

#ifndef GETFEM_PLASTICITY
#define GETFEM_PLASTICITY

#include <getfem_assembling.h> /* import assembly methods (and comp. of norms) */
#include <getfem_export.h>   /* export functions (save the solution in a file) */
#include <getfem_regular_meshes.h>
#include <gmm.h>

/* try to enable the SIGFPE if something evaluates to a Not-a-number of infinity
 * during computations
 */
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector;  /* special class for small (dim < 16) vectors */
using bgeot::base_node;   /* geometrical nodes (derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */
using bgeot::base_vector;
using bgeot::base_tensor;

using std::cout;
using std::cerr;
using std::cin;

/* definition of some matrix/vector types. These ones are built
   using the predefined types in Gmm++ */
typedef gmm::rsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef gmm::col_matrix<sparse_vector_type> col_sparse_matrix_type;

// calcul de tau_m*Id
template<typename MAT> MAT tau_m_Id(const MAT& tau){
  // calcul de la trace de tau
  scalar_type trace=gmm::mat_trace(tau);
  // on recupere sa dimension (lignes ou colonnes)
  size_type size_of_tau=gmm::mat_nrows(tau);
  MAT taumId(size_of_tau,size_of_tau);
  gmm::copy(gmm::identity_matrix(),taumId);
  gmm::scale(taumId,trace/size_of_tau);
  return taumId;
} 

// calcul de tau_d
template<typename MAT> MAT tau_d(const MAT& tau){
  // on recupere la dimension de tau (lignes ou colonnes)
  size_type size_of_tau=gmm::mat_nrows(tau);
  MAT taud(size_of_tau,size_of_tau);
  gmm::copy(gmm::scaled(tau_m_Id(tau),-1.0),taud); 
  gmm::add(tau,taud);
  return taud;
}

class type_proj {
  public :
    // if flag_proj=0 il output proj will be Proj(tau)
    // if flag_proj=1 il output proj will be gradProj(tau)
    // no others values allowed for flag_proj 
    virtual void compute_type_proj(const base_matrix& tau,const scalar_type stress_threshold, const scalar_type TOL, base_matrix& proj,const size_type flag_proj, const size_type flag_hyp)  const = 0;
};



class VM_projection : public type_proj  {
public :

  virtual void compute_type_proj(const base_matrix& tau ,const scalar_type stress_threshold, const scalar_type TOL, base_matrix& proj,const size_type flag_proj, const size_type flag_hyp)  const {

  

    
    // on verifie la valeur du flag_proj
    if(flag_proj!=0 && flag_proj!=1) DAL_THROW(dal::failure_error, "wrong value for the projection flag, must be 0 or 1 ");

    // on verifie que le seuil s est bien positif ou nul
    if(!(stress_threshold>=0.))
      DAL_THROW(dal::failure_error, "s is not a positive number "
		<< stress_threshold << ". You need to set "
		<< "s as a positive number");

    // on recupere la dimension de tau (lignes ou colonnes)
    size_type size_of_tau=gmm::mat_nrows(tau);

    // calcul du critere de determination du VM en fonction de flag_hyp : 
    // flag_hyp=0 : cas 3D, ou '2D plan' c'est normtaud directement
    // flag_hyp=1 : dans le cas 2D ca depend des valeurs propres de tau : - en contraintes planes : A IMPLEMENTER
    // flag_hyp= : dans le cas 2D ca depend des valeurs propres de tau :  - en deformations planes: A IMPLEMENTER 
 
    scalar_type normtaud;

 
    // on calcule tau_m*Id
    base_matrix taumId(size_of_tau,size_of_tau);
    gmm::copy(tau_m_Id(tau),taumId); 

    // calcul du deviateur de tau, taud
    base_matrix taud(size_of_tau,size_of_tau);
    gmm::copy(tau_d(tau),taud);


    // calcul de la norme de taud
    // depends on flag_hyp  : stress plane : 1
    //other for classical 3d : 
    //to be defined : plane strain, others ...
   

    // contraintes planes    
    if(flag_hyp==1){
      if(size_of_tau/=2) DAL_THROW(dal::failure_error, "wrong value for CALCULATION HYPOTHESIS, must be /=1 SINCE n/=2");
      // we form the 3D tau tensor considering that tau(3,j)=tau(i,3)=0
      base_matrix tau_aux(3,3); gmm::clear(tau_aux);
      gmm::copy(tau,gmm::sub_matrix(tau_aux,gmm::sub_interval(0,2)));
      // we calculate tau deviator and its norms
      base_matrix taud_aux(3,3);
      gmm::copy(tau_d(tau_aux),taud_aux);             
      normtaud=gmm::mat_euclidean_norm(taud_aux);  
    }

    else{
      normtaud=gmm::mat_euclidean_norm(taud);  
    }
    
    // on dimensionne et on initialise la matrice de projection ou sa derivee
    switch(flag_proj) {
      case 0:
	gmm::resize(proj,size_of_tau,size_of_tau);
        break;
      case 1:
	gmm::resize(proj,size_of_tau*size_of_tau, size_of_tau*size_of_tau);
	gmm::copy(gmm::identity_matrix(),proj);
	break;
      }
    
    // calcul de la projection ou de sa derivee

    
    // CAS OU normtaud=0 
    if (-TOL<=normtaud && normtaud<=TOL){
   
      switch(flag_proj) {
      case 0:
	gmm::copy(tau,proj);
        break;
      }

    }// FIN CAS OU normtaud=0    
 

  // CAS OU normtaud/=0 
  else{  
 
    //cas ou normtaud<s
    if (normtaud-stress_threshold< -TOL){
      switch(flag_proj) {
      case 0:
	gmm::copy(tau,proj);
        break;
      }
    }//fin cas ou normtaud<s
     

    //cas ou normtaud>s
    else if(normtaud-stress_threshold>TOL){
      // cout << "normtaud(>s) = " << normtaud <<"\n";
      // DAL_THROW(dal::failure_error, "norm(taud)>s ");
      switch(flag_proj) {
      case 0:
	gmm::copy(gmm::scaled(taud, stress_threshold/normtaud),proj);
	gmm::add(taumId,proj);
        break;
      case 1:

	// calcul de Igrad
         base_matrix Igrad(size_of_tau*size_of_tau,size_of_tau*size_of_tau);
	 gmm::copy(gmm::identity_matrix(),Igrad); 

 
	 // calcul de Igrad(*)Igrad
	 base_matrix Igrad2(size_of_tau*size_of_tau,size_of_tau*size_of_tau);
	 gmm::clear(Igrad2);
        
         // on construit le vecteur [1 0 0 1  0 0 1...] qui va etre copie dans certaines colonnes de Igrad(*)Igrad
         base_vector aux(size_of_tau*size_of_tau);
	 gmm::clear(aux);
         // boucle sur le vecteur aux pour placer les 1
	 for(size_type i=0;i<size_of_tau;++i)
	   for(size_type j=0;j<size_of_tau;++j)
	     if(i==j) aux[j*size_of_tau+i]=1.;

         // on copie aux dans les colonnes de Igrad(*)Igrad qu'il faut
	 for(size_type i=0;i<size_of_tau;++i)
	   for(size_type j=0;j<size_of_tau;++j)
	       if(i==j)  gmm::copy(aux, gmm::mat_col(Igrad2,j*size_of_tau+i)); 

	 //calcul de Id_grad
	 base_matrix Id_grad(size_of_tau*size_of_tau,size_of_tau*size_of_tau);
	 gmm::copy(gmm::scaled(Igrad2,-1./size_of_tau),Id_grad);
	 gmm::add(Igrad,Id_grad);         


         //calcul de ngrad(*)ngrad
         base_matrix ngrad2(size_of_tau*size_of_tau,size_of_tau*size_of_tau);
          // calcul de la normale n
	 base_matrix un(size_of_tau,size_of_tau);
	 gmm::copy(gmm::scaled(taud,1./normtaud),un);  

         // on copie la normale dasn un vecteur colonne range dans l'ordre du fortran
	 std::copy(un.begin(),un.end(),aux.begin());
	
	 // boucle sur les colonnes de ngrad(*)ngrad
	 for(size_type j=0;j<size_of_tau*size_of_tau;++j)
	      gmm::copy(gmm::scaled(aux,aux[j]), gmm::mat_col(ngrad2,j));
         
 
         //calcul final du gradient de la projection
         // NB : proj contient l'identite au depart (proj=Igrad)
         gmm::add(gmm::scaled(ngrad2,-1.),proj);
	 base_matrix aux2(size_of_tau*size_of_tau,size_of_tau*size_of_tau);
	 gmm::copy(gmm::scaled(proj,stress_threshold/normtaud),aux2);
	 gmm::mult(aux2,Id_grad,proj);
	 gmm::add(gmm::scaled(Igrad2,1./size_of_tau),proj);

         break;
      }

    
    }//fin cas ou normtaud>s


    //cas ou normtaud=s 
    else {

       switch(flag_proj) {
      
       case 0:
	 gmm::copy(tau,proj);
	 break;

       case 1:
	 // calcul de gradproj=Igrad
	 gmm::copy(gmm::identity_matrix(),proj); 
	 break;

       }

      
    }//fin cas ou normtaud=s 


  }// FIN CAS ou normtaud/=0 
  
  }


};


// sert a sauver les numéros de convexe dans leur ordre d'assemblage
// TODO : grouik
 std::ofstream sortie2("cv.txt");



// calcul de la projection de D*e + sigma_bar_ en un point de Gauss
class plasticity_projection_base  {
 protected:
  size_type N;
  const getfem::mesh_fem &mf;
  const std::vector<scalar_type> &U;
  const scalar_type stress_threshold, TOL, lambda, mu;  
  bgeot::multi_index sizes_;
  const type_proj *t_proj; 
  // on passe un reference sur sigma_bar pour des questions de rapidite, mais la reference est indispensable dans le constructeur car on veut modifier l'objet passe en argument a chaque fois ( possible vu que plasticity_projection_base sera utilise comme mutable) 
  std::vector<std::vector<scalar_type> > &sigma_bar_;

  // pour stocker la projection
  std::vector<std::vector<scalar_type> > &saved_proj_;
 
  const size_type flag_proj, flag_hyp;
  
  bool fill_sigma_bar, fill_saved_proj;
  
 public:  

  // methode pour retourner la variable sigma_bar en tous les points de Gauss de tous les convexes
   std::vector<std::vector<scalar_type> > &sigma_bar() {
    return sigma_bar_;
  }
    
  // methode pour retourner un element d'un tenseur sigma_bar en un point de Gauss d'un convexe donne
  // ne marche que pour les sigma_bar issus de la projection et non pas de sa derivee
  scalar_type &sigma_bar(size_type cv, size_type ii, int i, int j) {
    return sigma_bar_[cv][ii*N*N + j*N + i];
  }

   // methode pour retourner la projection en tous les points de Gauss de tous les convexes
   std::vector<std::vector<scalar_type> > &saved_proj() {
    return saved_proj_;
  }
    
  // methode pour retourner un element d'un tenseur projection en un point de Gauss d'un convexe donne
  // ne marche que pour les "proj" issus de la projection et non pas de sa derivee
  scalar_type &saved_proj(size_type cv, size_type ii, int i, int j) {
    return saved_proj_[cv][ii*N*N + j*N + i];
  }


  // constructeur
  plasticity_projection_base(const getfem::mesh_fem &mf_, const std::vector<scalar_type> &U_, 
			     const scalar_type stress_threshold_, const scalar_type TOL_, 
			     const scalar_type lambda_, const scalar_type mu_, 
			     const type_proj *t_proj_, std::vector<std::vector<scalar_type> > &sigma_bar__, 
                             std::vector<std::vector<scalar_type> > &saved_proj__,
			     const size_type flag_proj_, const size_type flag_hyp_, const bool fill_sigma) :
    mf(mf_), U(U_), stress_threshold(stress_threshold_), TOL(TOL_), lambda(lambda_), mu(mu_),t_proj(t_proj_),
    sigma_bar_( sigma_bar__),saved_proj_(saved_proj__), flag_proj(flag_proj_), flag_hyp(flag_hyp_)  {
    
    // mf(mf_) ; constructeur de recopie de l'objet getfem::mesh_fem &mf. NB : pour cet objet on aurait pau utiliser un tout autre constructeur(pas forcement un constructeur de recopie, ex mf(a,b,c)
    // lambda(lambda_) : lambda n'est pas un objet, simple initialisation en en-tete. 
    
    fill_sigma_bar = fill_sigma;   //en principe c false, sauf apres le newton;
    fill_saved_proj=false;  
    
    N = mf.linked_mesh().dim(); 
    
    if (mf.get_qdim() != N) DAL_THROW(dal::failure_error, "wrong qdim for the mesh_fem");      
    
    /* on fixe la taille du tenseur a renvoyer en fonction du cas traite */
    // calcul de la projection : tenseur N*N
    if(flag_proj==0){
      sizes_.resize(2); sizes_[0] = N; sizes_[1] = N; 
    }
    // calcul de la derivee de la projection : tenseur N*N*N*N
    else if(flag_proj==1){
      sizes_.resize(4); sizes_[0] = N; sizes_[1] = N; sizes_[2] = N; sizes_[3] = N;
      }
    // valeurs interdites pour flag_proj
    else
      DAL_THROW(dal::failure_error, "wrong value for the projection flag, must be 0 or 1 ");   
    
    // autant de sigma_bar_ que le nb max de convexes
    sigma_bar_.resize(mf.linked_mesh().convex_index().last_true()+1);
    
      // pareil pour la projection
    saved_proj_.resize(mf.linked_mesh().convex_index().last_true()+1);
    
  }

  // lorsqu'il faudra sauver dans sigma_bar_ le sigma_bar du dernier newton au temps n, utiliser set_fill_sigma_bar
  void set_fill_sigma_bar(bool b) { fill_sigma_bar = b; }

   // lorsqu'il faudra sauver dans saved_proj_ le saved_proj du dernier newton a certains temps n' inclus dans les instants n , utiliser set_fill_saved_proj
  void set_fill_saved_proj(bool b) { fill_saved_proj = b; }

  // methode sizes() de l'element non lineaire, renvoie la taille du tenseur en output
  const bgeot::multi_index &sizes() const { return sizes_; }

  // methode compute() de l'element non lineaire, calcule le tenseur en output
  void compute_proj(getfem::fem_interpolation_context& ctx, bgeot::base_tensor &t){ 

    size_type cv = ctx.convex_num();

    size_type ii = ctx.ii(); /* numero du pt d'integration dans le convexe */

    // on sauve les numéros de convexe dans leur ordre d'assemblage
    sortie2 << cv <<" "; 


    getfem::pfem pf = ctx.pf();
    base_vector coeff(mf.nb_dof_of_element(cv));
    base_matrix gradU(N, N), sigma(N,N);
    gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf.ind_dof_of_element(cv))), coeff);
    pf->interpolation_grad(ctx, coeff, gradU, mf.get_qdim());
    scalar_type ltrace_eps = lambda*gmm::mat_trace(gradU);
    
    // si besoin est on donne a sigma_bar[cv] une taille egale au nombre de points d'integration maximal par convexe. Apparement en pratique cette taille n'est que rarement atteinte
    
    if (sigma_bar_[cv].size() == 0){
      // cout << "\nctx..."<< mf.int_method_of_element(cv)->approx_method()->nb_points_on_convex() << endl;
      size_type nbgausspt = mf.int_method_of_element(cv)->approx_method()->nb_points_on_convex();
      sigma_bar_[cv].resize(N*N*nbgausspt);
      gmm::clear(sigma_bar_[cv]);

      // pour la projection egalement
      saved_proj_[cv].resize(N*N*nbgausspt);
      gmm::clear(saved_proj_[cv]);
    }


    t.adjust_sizes(sizes_);
   
    for (size_type i=0; i < N; ++i) {
      for (size_type j=0; j < N; ++j) {
	if(i==j)
	  sigma(i,j) = 2*mu*(gradU(i,j)+gradU(j,i))/2. + ltrace_eps + sigma_bar(cv,ii,i,j);
	else
	  sigma(i,j) = 2*mu*(gradU(i,j)+gradU(j,i))/2. + sigma_bar(cv,ii,i,j);
     
      }
   
    }
    

    base_matrix tau_star(N,N), gradproj(N,N), proj;
 
    t_proj->compute_type_proj(sigma, stress_threshold, TOL, proj, flag_proj, flag_hyp);

    // on ne remplit sigma_bar qu'a des moments precis (pas a tous les assemblages) et uniquement pour la projection, 
    // pas pour sa derivee 


 

    if (fill_sigma_bar && flag_proj==0) {

      for (size_type i=0; i < N; ++i)
	for (size_type j=0; j < N; ++j)
	    saved_proj(cv,ii,i,j) = proj(i,j);

      gmm::add(gmm::scaled(gradU, -mu), proj);
      gmm::add(gmm::scaled(gmm::transposed(gradU), -mu), proj);
      //cout << proj<<"\n";

      for (size_type i=0; i < N; ++i) proj(i,i) += ltrace_eps;

      //cout << proj<<"\n";  
    
      for (size_type i=0; i < N; ++i) {
	for (size_type j=0; j < N; ++j) {
          
	  //saved_proj ne peut etre sauve qu a certains temps n' forcement inclus dans l'ensemble des instants n de sauvegarde de sigma_bar qui lui est sauve a tous les temps
	    sigma_bar(cv,ii,i,j) = proj(i,j);
	}
      }
    }

   
    std::copy(proj.begin(),proj.end(),t.begin());

  }
};


class plasticity_projection : public getfem::nonlinear_elem_term {
  plasticity_projection_base *p;
public:  
  plasticity_projection(plasticity_projection_base *p_) : p(p_) {}
  virtual const bgeot::multi_index &sizes() const { return p->sizes(); }

  virtual void compute(getfem::fem_interpolation_context& ctx,
		       bgeot::base_tensor &t) { return p->compute_proj(ctx,t); }
};

#endif
