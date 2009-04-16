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

#ifndef GETFEM_MODELS_H__
#define GETFEM_MODELS_H__

#include "getfem_assembling.h"
#include "getfem_partial_mesh_fem.h"

namespace getfem {


  typedef gmm::rsvector<scalar_type> model_real_sparse_vector;
  typedef gmm::rsvector<complex_type> model_complex_sparse_vector;
  typedef std::vector<scalar_type> model_real_plain_vector;
  typedef std::vector<complex_type> model_complex_plain_vector;

  // utiliser le même type que l'interface matlab/python pour représenter
  // les vecteurs/matrices ?
  // Cela faciliterait les échanges et réduirait les composantes de la 
  // classe model.
  
  typedef gmm::col_matrix<model_real_sparse_vector> model_real_sparse_matrix;
  typedef gmm::col_matrix<model_complex_sparse_vector>
    model_complex_sparse_matrix;

 
  class model : public context_dependencies {

    // State variables of the model
    bool complex_version;
    model_real_sparse_matrix rTM;      // tangent matrix, real version
    model_complex_sparse_matrix cTM;   // tangent matrix, complex version
    model_real_plain_vector rstate;
    model_complex_plain_vector cstate;
    model_real_plain_vector rrhs;
    model_complex_plain_vector crhs;
    mutable bool act_size_to_be_done;

    // Variables and parameters of the model

    enum  var_description_filter {
      VDESCRFILTER_NO,     // Variable being directly the dofs of a given fem
      VDESCRFILTER_REGION, /* Variable being the dofs of a fem on a mesh region
			    * (uses mf.dof_on_region). */
      VDESCRFILTER_INFSUP  /* Variable being the dofs of a fem on a mesh region
			    * with an additional filter on a mass matrix with
                            * respect to another fem. */
    };

    struct var_description {
     
      // + format (scalaire, vectoriel, tensoriel ordre 2 ou 4)
      // + mim associée (avec possibilité de mim par défaut par fem ...).
      bool is_variable;  // This is a variable or a parameter.
      bool is_complex;   // The variable is complex numbers
      bool is_fem_dofs;  // The variable is the dofs of a fem
      var_description_filter filter; // A filter on the dofs is applied or not.
      size_type n_iter; //  number of version of the variable stored (for time
                  // integration schemes.

      // fem description of the variable
      const mesh_fem *mf;          // Principal fem of the variable.
      const mesh_im *mim;          // Optional mesh_im for filter.
      partial_mesh_fem *partial_mf; // Filter with respect to mf.
      int m_region;                // Optional mesh_region for the filter.
      std::string filter_var;      // Optional variable name for the filter
                           // with the mass matrix of the corresponding fem.

      dim_type qdim;  // une donnée peut avoir un qdim != du fem.
                      // dim per dof for dof data.
      gmm::uint64_type v_num;

      gmm::sub_interval I; // For a variable : indices on the whole system
  
      std::vector<model_real_plain_vector> real_value;
      std::vector<model_complex_plain_vector> complex_value;

      var_description(bool is_var = false, bool is_com = false,
		      bool is_fem = false, size_type n_it = 1,
		      var_description_filter fil = VDESCRFILTER_NO,
		      const mesh_fem *mmf = 0, const mesh_im *im = 0,
		      int m_reg = 0, dim_type Q = 1,
		      const std::string &filter_v = std::string(""))
	: is_variable(is_var), is_complex(is_com), is_fem_dofs(is_fem),
	  filter(fil), n_iter(n_it), mf(mmf), mim(im), partial_mf(0),
	  m_region(m_reg),
	  filter_var(filter_v), qdim(Q), v_num(act_counter()) {
	if (filter != VDESCRFILTER_NO && mf != 0)
	  partial_mf = new partial_mesh_fem(*mf);
      }

      const mesh_fem &associated_mf(void) const {
	GMM_ASSERT1(is_fem_dofs, "This variable is not linked to a fem");
	if (filter == VDESCRFILTER_NO) return *mf; else return *partial_mf; 
      }

      ~var_description() { if (partial_mf) delete partial_mf; }

      size_type size(void) const // devrait contrôler que la variable
      // a bien été initialisée avec la bonne taille par actualize_sizes.
      { return is_complex ? complex_value[0].size() : real_value[0].size(); }

      void set_size(size_type s);
    };
    
    typedef std::map<std::string, var_description> VAR_SET;
    VAR_SET variables;

    void actualize_sizes(void);
    bool check_name_valitity(const std::string &name,
			     bool assert = true) const;


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



    void init(void) { complex_version = false; act_size_to_be_done = false; }

  public :

    // gerer aussi l'actualisation des mult qui n'est pas fait au début 
    // état du modèle -> actualize à faire ...
    void update_from_context(void) const {  act_size_to_be_done = true; };
    
    /** Boolean which says if the model deals with real or complex unknowns
	and data. */
    bool is_complex(void) const { return complex_version; }

    /** Add a fixed size variable to the model. niter is the number of version
	of the variable stored, for time integration schemes. */
    void add_fixed_size_variable(const std::string &name, size_type size,
				 size_type niter = 1);

    /** Add a fixed size data to the model. niter is the number of version
	of the data stored, for time integration schemes. */
    void add_fixed_size_data(const std::string &name, size_type size,
				 size_type niter = 1);

    /** Add a variable being the dofs of a finite element method to the model.
	niter is the number of version of the variable stored, for time
	integration schemes. */
    void add_fem_variable(const std::string &name, const mesh_fem &mf,
			  size_type niter = 1);

    /** Add a data being the dofs of a finite element method to the model.
	niter is the number of version of the data stored, for time
	integration schemes. */
     void add_fem_data(const std::string &name, const mesh_fem &mf,
			  dim_type qdim = 1, size_type niter = 1);
    
    /** Add a particular variable linked to a fem beeing a multiplier with
	respect to a primal variable. The dof will be filtered with a mass
	matrix to retain only linearly independant constraints on the primal
	variable. niter is the number of version of the data stored, for time
	integration schemes. */
    void add_mult_on_region(const std::string &name, const mesh_fem &mf,
			    const mesh_im &mim,
			    const std::string &primal_name, int region,
			    size_type niter = 1);

    /** Gives the access to the vector value of a variable. For the real
	version. */
    const model_real_plain_vector &
    real_variable(const std::string &name, size_type niter = 1) const;

    /** Gives the access to the vector value of a variable. For the complex
	version. */
    const model_complex_plain_vector &
    complex_variable(const std::string &name, size_type niter = 1) const;

    /** Gives the access to the vector value of a variable. For the real
	version. */
    model_real_plain_vector &
    real_variable(const std::string &name, size_type niter = 1);

    /** Gives the access to the vector value of a variable. For the complex
	version. */
    model_complex_plain_vector &
    complex_variable(const std::string &name, size_type niter = 1);

    /** Gives the access to the tangent matrix. For the real version. */
    const model_real_sparse_matrix &real_tangent_matrix() const {
      GMM_ASSERT1(!complex_version, "This model is a complex one");
      // context_check(); if (act_size_to_be_done) actualize_sizes();
      return rTM;
    }
    
    /** Gives the access to the tangent matrix. For the complex version. */
    const model_complex_sparse_matrix &complex_tangent_matrix() const {
      GMM_ASSERT1(complex_version, "This model is a real one");
      // context_check(); if (act_size_to_be_done) actualize_sizes();
      return cTM;
    }
    
    /** Gives the access to the right hand side of the tangent linear system.
	For the real version. */
    const model_real_plain_vector &real_rhs() const {
      GMM_ASSERT1(!complex_version, "This model is a complex one");
      // context_check(); if (act_size_to_be_done) actualize_sizes();
      return rrhs;
    }
    
    /** Gives the access to the right hand side of the tangent linear system.
	For the complex version. */
    const model_complex_plain_vector &complex_rhs() const {
      GMM_ASSERT1(complex_version, "This model is a real one");
      // context_check(); if (act_size_to_be_done) actualize_sizes();
      return crhs;
    }

    const model_real_plain_vector &real_state() const {
      GMM_ASSERT1(!complex_version, "This model is a complex one");
      // context_check(); if (act_size_to_be_done) actualize_sizes();
      return rstate;
    }
    
    const model_complex_plain_vector &complex_state() const {
      GMM_ASSERT1(complex_version, "This model is a real one");
      // context_check(); if (act_size_to_be_done) actualize_sizes();
      return cstate;
    }

    model_real_plain_vector &real_state() {
      GMM_ASSERT1(!complex_version, "This model is a complex one");
      // context_check(); if (act_size_to_be_done) actualize_sizes();
      return rstate;
    }
    
    model_complex_plain_vector &complex_state() {
      GMM_ASSERT1(complex_version, "This model is a real one");
      context_check(); if (act_size_to_be_done) actualize_sizes();
      return cstate;
    }

    /** List the model variables and constant */
    void varlist(std::ostream &ost) const;

    void clear(void) {
      variables.clear();
      rTM = model_real_sparse_matrix();
      cTM = model_complex_sparse_matrix();
      rstate = rrhs = model_real_plain_vector();
      cstate = crhs = model_complex_plain_vector();
    }

    model(bool comp_version = false)
    { init(); complex_version = comp_version; }

  };

  
  

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MODELS_H__  */
