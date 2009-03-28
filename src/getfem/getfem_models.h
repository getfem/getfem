// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2009-2009 Yves Renard
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
// As a special exception, you  may use  this file  as it is a part of a free
// software  library  without  restriction.  Specifically,  if   other  files
// instantiate  templates  or  use macros or inline functions from this file,
// or  you compile this  file  and  link  it  with other files  to produce an
// executable, this file  does  not  by itself cause the resulting executable
// to be covered  by the GNU Lesser General Public License.  This   exception
// does not  however  invalidate  any  other  reasons why the executable file
// might be covered by the GNU Lesser General Public License.
//
//===========================================================================

/**
   @file getfem_models.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date March 21, 2009.
   @brief Model representation in Getfem.
*/

// TODO étape 1 :
// 1 - Refaire le partial mesh_fem en le faisant géré par la structure mesh
//     fem elle même dans le gout de se qui est fait pour la vectorisation
//     implicite. Faire qlq chose de plus carré et documenté qui permette
//     de séléctionner ou à la limite de faire des combinaisons linéaires de
//     dofs au niveau de la structure mesh_fem sans refaire des pfems ce qui
//     est très couteux. Devrait gérer aussi le dof_partition.
// 2 - Remplacer toute les utilisation de partial_mesh_fem avec la nouvelle
//     structure. Voir au niveau de Xfem ce que l'on peut gérer comme ca
//     (et fem sum ...).
// 3 - Compléter la gestion des variables dans les modèles, leur taille,
//     prise en compte des événements.
// 4 - documentation (avec avertissement -> pour version 4.0 de Getfem)
//     passer à 4.0 avant de committer.



#ifndef GETFEM_MODELS_H__
#define GETFEM_MODELS_H__

#include "getfem_assembling.h"
#include "getfem_partial_mesh_fem.h"

namespace getfem {


  typedef gmm::rsvector<scalar_type> model_real_sparse_vector;
  typedef gmm::rsvector<complex_type> model_complex_sparse_vector;
  typedef std::vector<scalar_type> model_real_plain_vector;
  typedef std::vector<complex_type> model_complex_plain_vector;
  
  typedef gmm::col_matrix<model_real_sparse_vector> model_real_sparse_matrix;
  typedef gmm::col_matrix<model_complex_sparse_vector>
    model_complex_sparse_matrix;

 
  class model {

    // State variables of the model
    bool complex_version;
    model_real_sparse_matrix rTM;      // tangent matrix, real version
    model_complex_sparse_matrix cTM;   // tangent matrix, complex version
    model_real_plain_vector rstate;
    model_complex_plain_vector cstate;
    model_real_plain_vector rrhs;
    model_complex_plain_vector crhs;

    // Variables and parameters of the model

    typedef enum  var_description_filter {
      VDESCRFILTER_NO,     // Variable being directly the dofs of a given fem
      VDESCRFILTER_REGION, // Variable being the dofs of a fem on a mesh region
      VDESCRFILTER_INFSUP  // Variable being the dofs of a fem on a mesh region
                           // with an additional filter on a mass matrix with
                           // respect to another fem.
    };

    struct var_description {
     
      // + format (scalaire, vectoriel, tensoriel ordre 2 ou 4)
      // + mim associée (avec possibilité de mim par défaut par fem ...).
      bool is_variable;  // This is a variable or a parameter.
      bool is_complex;   // The variable is complex numbers
      bool is_fem_dofs;  // The variable is the dofs of a fem
      var_description_filter filter; // A filter on the dofs is applied or not.
      int n_iter; //  number of version of the variable stored (for time
                  // integration schemes.

      // fem description of the variable
      const mesh_fem *mf;          // Principal fem of the variable.
      // la structure partial_mesh_fem n'est pas très efficace (elle définie
      // de nouveau fems), il faudrait
      // donner une certaine liberté entre fem et mesh_fem qui permette
      // de faire l'opération de séléction à l'interieur du mesh_fem
      // comme est faite l'opération qui consiste à vectoriser un fem non
      // vectoriel.
      partial_mesh_fem partial_mf; // Filter with repsect to mf.
      int m_region;                // Optional mesh_region for the filter.
      std::string filter_var;      // Optional variable name for the filter
                           // with the mass matrix of the correpsonding fem.

      size_type qdim;  // une donnée peut avoir un qdim != du fem.
                       // dim per dof for dof variables.

      gmm::sub_interval I; // For a variable : indices on the whole system
  
      model_real_plain_vector real_value;
      model_complex_plain_vector complex_value;



      const mesh_fem &associated_mf(void) const {
	GMM_ASSERT1(is_fem_dofs, "This variable is not linked to a fem");
	if (filter == VDESCRFILTER_NO) return *mf; else return partial_mf; 
      }

      size_type size(void) const // devrait contrôler que la variable
      // a bien été initialisée avec la bonne taille par actualize_sizes.
      { if (is_complex) complex_value.size(); else real_value.size(); }

      size_type set_size(size_type s) 
      { if (is_complex) complex_value.resize(s); else real_value.resize(s); }

      var_description(...) : partial_mf(*mf) {

      }
    };
    
    typedef typename std::map<std::string, var_description> VAR_SET;
    VAR_SET variables;

    actualize_sizes(void) { // corps à mettre dans le .cc;
      // mets à jour la taille des variables, les filtres ...
      // à refaire dés que : on ajoute une variable, un fem change, une
      //    variable non fem change de taille.

      // doit mettre aussi à jour les paramêtres en adaptant la taille
      // des vecteurs -> warning si la taille change.

      size_type tot_size = 0;
      
      for (VAR_SET::iterator it = variables.begin(); it != variables.end();
	   ++it) {
	if (it->is_fem_dofs) {
	  size_type s(0);
	  switch (it->filter) {
	  case VDESCRFILTER_NO:
	    s = it->mf->nb_dof();
	    if (!it->is_variable && it->mf->get_qdim() != it->qdim)
	      s *= it->qdim;
	    break;
	  case VDESCRFILTER_REGION: ...; break;
	  case VDESCRFILTER_INF_SUP: ...; break;
	  }
	  it->set_size(s);
	}
	if (it->is_variable) {
	  it->I = gmm::sub_interval(tot_size, it->size());
	  tot_size += it->size();
	}

	// + fixer la taille des vecteurs globaux du modèle ...

      }
      


    }

    // Faire la doc avec les différents types de variables gérées
    //   (fem, mult ...)

    // faire la fonction qui filtre toutes les variables à filtrer ... boucle
    //   + action globale sur les multiplicateurs d'une même variable.

    // gérer les dépendances entre l'objet et les mesh_fem, mesh_im
    //    : si un disparaît, l'objet est invalidé ...

    // faire la mise à jour des intervalles correspondant à chaque variable

    // centraliser les mesh_fems de base avec numérotation locale ?


    // il faudra gerer les changements de taille des variables si les fems
    // changent ... avec possibilité de ré-interpoler quand on raffine ...





    void init(void) { complex_version = false; }

  public :
    
    // interface for state variables

    bool is_complex_version(void) const { return complex_version; }

    const model_real_sparse_matrix &real_tangent_matrix() const {
      GMM_ASSERT1(!complex_version, "This model is a complex one");
      return rTM;
    }
    
    const model_complex_sparse_matrix &complex_tangent_matrix() const {
      GMM_ASSERT1(complex_version, "This model is a real one");
      return cTM;
    }
    
    const model_real_sparse_matrix &real_rhs() const {
      GMM_ASSERT1(!complex_version, "This model is a complex one");
      return rrhs;
    }
    
    const model_complex_sparse_matrix &complex_rhs() const {
      GMM_ASSERT1(complex_version, "This model is a real one");
      return crhs;
    }



    model(void) { init(); }

  };

  
  

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MODELS_H__  */
