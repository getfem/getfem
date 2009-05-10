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

#include "getfem_partial_mesh_fem.h"

namespace getfem {

  class virtual_brick;
  /** type of pointer on a brick */ 
  typedef boost::intrusive_ptr<const getfem::virtual_brick> pbrick;

  // Event management : The model has to react when something has changed in
  //    the context and ask for corresponding (linear) bricks to recompute
  //    some terms.
  //    For the moment two events are taken into account
  //      - Change in a mesh_fem
  //      - Change in the data of a variable
  //    For this, a brick has to declare on which variable it depends and
  //    on which data. When a linear brick depend on a variable, the
  //    recomputation is done when the eventual corresponding mesh_fem
  //    is changed (or the size of the variable for a fixed size variable).
  //    When a linear brick depend on a data, the recomputation is done
  //    when the corresponding vector value is changed. If a variable is used
  //    as a data, it has to be declared as a data by the brick.
  //    A nonlinear brick is recomputed at each assembly of the tangent system.
  //    Remember this behavior when some changed are done on the variable
  //    and/or data.
  //    The change on a mesh_im is not taken into account for the moment.
  //    The different versions of the variables is not taken into account
  //    separately.

  //=========================================================================
  //
  //  Model object.
  //
  //=========================================================================


  typedef gmm::wsvector<scalar_type> model_real_sparse_vector;
  typedef gmm::wsvector<complex_type> model_complex_sparse_vector;
  typedef std::vector<scalar_type> model_real_plain_vector;
  typedef std::vector<complex_type> model_complex_plain_vector;

  // utiliser le même type que l'interface matlab/python pour représenter
  // les vecteurs/matrices ?
  // Cela faciliterait les échanges et réduirait les composantes de la 
  // classe model.
  
  typedef gmm::col_matrix<model_real_sparse_vector> model_real_sparse_matrix;
  typedef gmm::col_matrix<model_complex_sparse_vector>
  model_complex_sparse_matrix;

  /** ``Model'' variables store the variables, the data and the
      description of a model. This includes the global tangent matrix, the
      right hand side and the constraints. There are two kinds of models, the
      ``real'' and the ``complex'' models.
  */
  class model : public context_dependencies {

    // State variables of the model
    bool complex_version;
    bool is_linear_;
    bool is_symmetric_;
    bool is_coercive_;
    mutable model_real_sparse_matrix rTM;    // tangent matrix, real version
    mutable model_complex_sparse_matrix cTM; // tangent matrix, complex version
    mutable model_real_plain_vector rrhs;
    mutable model_complex_plain_vector crhs;
    mutable bool act_size_to_be_done;
    dim_type leading_dim;

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
     
      bool is_variable;  // This is a variable or a parameter.
      bool is_complex;   // The variable is complex numbers
      bool is_fem_dofs;  // The variable is the dofs of a fem
      var_description_filter filter; // A filter on the dofs is applied or not.
      size_type n_iter; //  number of version of the variable stored (for time
                        // integration schemes.

      // fem description of the variable
      const mesh_fem *mf;           // Principal fem of the variable.
      const mesh_im *mim;           // Optional mesh_im for filter.
      ppartial_mesh_fem partial_mf; // Filter with respect to mf.
      size_type m_region;           // Optional mesh_region for the filter.
      std::string filter_var;       // Optional variable name for the filter
      // with the mass matrix of the corresponding fem.

      dim_type qdim;  // une donnée peut avoir un qdim != du fem.
      // dim per dof for dof data.
      gmm::uint64_type v_num, v_num_data;

      gmm::sub_interval I; // For a variable : indices on the whole system
  
      std::vector<model_real_plain_vector> real_value;
      std::vector<model_complex_plain_vector> complex_value;

      var_description(bool is_var = false, bool is_com = false,
		      bool is_fem = false, size_type n_it = 1,
		      var_description_filter fil = VDESCRFILTER_NO,
		      const mesh_fem *mmf = 0, const mesh_im *im = 0,
		      size_type m_reg = size_type(-1), dim_type Q = 1,
		      const std::string &filter_v = std::string(""))
	: is_variable(is_var), is_complex(is_com), is_fem_dofs(is_fem),
	  filter(fil), n_iter(n_it), mf(mmf), mim(im),
	  m_region(m_reg),
	  filter_var(filter_v), qdim(Q), v_num(0), v_num_data(act_counter()) {
	if (filter != VDESCRFILTER_NO && mf != 0)
	  partial_mf = new partial_mesh_fem(*mf);
	v_num_data = v_num;
      }

      const mesh_fem &associated_mf(void) const {
	GMM_ASSERT1(is_fem_dofs, "This variable is not linked to a fem");
	if (filter == VDESCRFILTER_NO) return *mf; else return *partial_mf; 
      }

      const mesh_fem *passociated_mf(void) const {
	if (!is_fem_dofs) return 0;
	if (filter == VDESCRFILTER_NO) return mf;
	else return partial_mf.get(); 
      }

      size_type size(void) const // devrait contrôler que la variable
      // a bien été initialisée avec la bonne taille par actualize_sizes.
      { return is_complex ? complex_value[0].size() : real_value[0].size(); }

      void set_size(size_type s);
    };

  public :

    typedef var_description *pvariable;
    typedef std::vector<std::string> varnamelist;
    typedef std::vector<const mesh_im *> mimlist;
    typedef std::vector<model_real_sparse_matrix> real_matlist;
    typedef std::vector<model_complex_sparse_matrix> complex_matlist;
    typedef std::vector<model_real_plain_vector> real_veclist;
    typedef std::vector<model_complex_plain_vector> complex_veclist;

    struct term_description {
      bool is_matrix_term; // tangent matrix term or rhs term.
      bool is_symmetric;   // Term have to be symmetrized.
      std::string var1, var2;
      term_description(const std::string &v)
	: is_matrix_term(false), var1(v) {}
      term_description(const std::string &v1, const std::string &v2,
		       bool issym)
	: is_matrix_term(true), is_symmetric(issym), var1(v1), var2(v2) {}
    };

    typedef std::vector<term_description> termlist;

  private :


    // rmatlist and cmatlist could be csc_matrix vectors to reduced the
    // amount of memory (but this should add a supplementary copy).
    struct brick_description {
      bool terms_to_be_computed;
      gmm::uint64_type v_num;
      pbrick pbr;               // brick pointer
      varnamelist vlist;        // List of variables used by the brick.
      varnamelist dlist;        // List of data used by the brick.
      termlist tlist;           // List of terms build by the brick
      mimlist mims;             // List of integration methods.
      size_type region;         // Optional region size_type(-1) for all.
      real_matlist rmatlist;    // List of matrices the brick have to fill in
      // (real version).
      real_veclist rveclist;    // List of rhs the brick have to fill in
      // (real version).
      complex_matlist cmatlist; // List of matrices the brick have to fill in
      // (complex version).
      complex_veclist cveclist; // List of rhs the brick have to fill in
      // (complex version).
      
      brick_description(pbrick p, const varnamelist &vl,
			const varnamelist &dl, const termlist &tl,
			const mimlist &mms, size_type reg)
	: terms_to_be_computed(true), v_num(0), pbr(p), vlist(vl), dlist(dl),
	  tlist(tl), mims(mms), region(reg) { }
    };
    
    typedef std::map<std::string, var_description> VAR_SET;
    mutable VAR_SET variables;
    std::vector<brick_description> bricks;

    void actualize_sizes(void) const;
    bool check_name_valitity(const std::string &name,
			     bool assert = true) const;

    void init(void) { complex_version = false; act_size_to_be_done = false; }

  public :

    // gerer aussi l'actualisation des mult qui n'est pas fait au début 
    // état du modèle -> actualize à faire ...
    void update_from_context(void) const {  act_size_to_be_done = true; }
    
    /** Boolean which says if the model deals with real or complex unknowns
	and data. */
    bool is_complex(void) const { return complex_version; }

    /** Return true if all the model terms do not affect the coercivity of
	the whole tangent system. */    
    bool is_coercive(void) const { return is_coercive_; }

    /** Return true if all the model terms are linear. */    
    bool is_linear(void) const { return is_linear_; }

    /** Total number of degrees of freedom in the model. */
    size_type nb_dof(void) const {
      context_check(); if (act_size_to_be_done) actualize_sizes();
      if (complex_version)
	return gmm::vect_size(crhs);
      else
	return gmm::vect_size(rrhs);
    }

    /** Leading dimension of the meshes used in the model. */
    dim_type leading_dimension(void) const { return leading_dim; }

    /** Gives a non already existing variable name bigining by `name`. */
    std::string new_name(const std::string &name);

    /** Gives the access to the vector value of a variable. For the real
	version. */
    const model_real_plain_vector &
    real_variable(const std::string &name, size_type niter = 0) const;

    /** Gives the access to the vector value of a variable. For the complex
	version. */
    const model_complex_plain_vector &
    complex_variable(const std::string &name, size_type niter = 0) const;

    /** Gives the write access to the vector value of a variable. Make a
	change flag of the variable set. For the real version. */
    model_real_plain_vector &
    set_real_variable(const std::string &name, size_type niter = 0);

    /** Gives the write access to the vector value of a variable. Make a
	change flag of the variable set. For the complex version. */
    model_complex_plain_vector &
    set_complex_variable(const std::string &name, size_type niter = 0);

    template<typename VECTOR, typename T>
    void from_variables(VECTOR &V, T) const {
      for (VAR_SET::iterator it = variables.begin();
	   it != variables.end(); ++it)
	if (it->second.is_variable)
	  gmm::copy(it->second.real_value[0],
		    gmm::sub_vector(V, it->second.I));
    }

    template<typename VECTOR, typename T>
    void from_variables(VECTOR &V, std::complex<T>) const {
      for (VAR_SET::iterator it = variables.begin();
	   it != variables.end(); ++it)
	if (it->second.is_variable)
	  gmm::copy(it->second.complex_value[0],
		    gmm::sub_vector(V, it->second.I));
    }

    template<typename VECTOR> void from_variables(VECTOR &V) const {
      typedef typename gmm::linalg_traits<VECTOR>::value_type T;
      context_check(); if (act_size_to_be_done) actualize_sizes();
      from_variables(V, T());
    }


    template<typename VECTOR, typename T>
    void to_variables(VECTOR &V, T) {
      for (VAR_SET::iterator it = variables.begin();
	   it != variables.end(); ++it)
	if (it->second.is_variable) {
	  gmm::copy(gmm::sub_vector(V, it->second.I),
		    it->second.real_value[0]);
	  it->second.v_num_data = act_counter();
	}
    }

    template<typename VECTOR, typename T>
    void to_variables(VECTOR &V, std::complex<T>) {
      for (VAR_SET::iterator it = variables.begin();
	   it != variables.end(); ++it)
	if (it->second.is_variable) {
	  gmm::copy(gmm::sub_vector(V, it->second.I),
		    it->second.complex_value[0]);
	  it->second.v_num_data = act_counter();
	}
    }

    template<typename VECTOR> void to_variables(VECTOR &V) {
      typedef typename gmm::linalg_traits<VECTOR>::value_type T;
      context_check(); if (act_size_to_be_done) actualize_sizes();
      to_variables(V, T());
    }

    /** Add a fixed size variable to the model. niter is the number of version
	of the variable stored, for time integration schemes. */
    void add_fixed_size_variable(const std::string &name, size_type size,
				 size_type niter = 1);

    /** Add a fixed size data to the model. niter is the number of version
	of the data stored, for time integration schemes. */
    void add_fixed_size_data(const std::string &name, size_type size,
			     size_type niter = 1);


    /** Add a fixed size data to the model initialized with V. */
    template <typename VECT>
    void add_initialized_fixed_size_data(const std::string &name, VECT &v) {
      this->add_fixed_size_data(name, gmm::vect_size(v), 1);
      if (this->is_complex()) // to be templated .. see later
	gmm::copy(v, this->set_complex_variable(name));
      else
	gmm::copy(gmm::real_part(v), this->set_real_variable(name));
    }

    /** Add a variable being the dofs of a finite element method to the model.
	niter is the number of version of the variable stored, for time
	integration schemes. */
    void add_fem_variable(const std::string &name, const mesh_fem &mf,
			  size_type niter = 1);

    /** Add a data being the dofs of a finite element method to the model.
	The data is initialized with V. */
    void add_fem_data(const std::string &name, const mesh_fem &mf,
		      dim_type qdim = 1, size_type niter = 1);

    /** Add a fixed size data to the model. niter is the number of version
	of the data stored, for time integration schemes. */
    template <typename VECT>
    void add_initialized_fem_data(const std::string &name, const mesh_fem &mf,
				  VECT &v) {
      this->add_fem_data(name, mf,
			 dim_type(gmm::vect_size(v) / mf.nb_dof()), 1);
      if (this->is_complex()) // to be templated .. see later
	gmm::copy(v, this->set_complex_variable(name));
      else
	gmm::copy(gmm::real_part(v), this->set_real_variable(name));
    }
    
    /** Add a particular variable linked to a fem being a multiplier with
	respect to a primal variable. The dof will be filtered with a mass
	matrix to retain only linearly independant constraints on the primal
	variable. niter is the number of version of the data stored, for time
	integration schemes. */
    void add_mult_on_region(const std::string &name, const mesh_fem &mf,
			    const mesh_im &mim,
			    const std::string &primal_name, size_type region,
			    size_type niter = 1);
    
    /** Gives the access to the mesh_fem of a variable if any. Throw an
	exception otherwise. */
    const mesh_fem &mesh_fem_of_variable(const std::string &name) const;

    /** Gives a pointer to the mesh_fem of a variable if any. 0 otherwise.*/
    const mesh_fem *pmesh_fem_of_variable(const std::string &name) const;


    /** Gives the access to the tangent matrix. For the real version. */
    const model_real_sparse_matrix &real_tangent_matrix(void) const {
      GMM_ASSERT1(!complex_version, "This model is a complex one");
      context_check(); if (act_size_to_be_done) actualize_sizes();
      return rTM;
    }
    
    /** Gives the access to the tangent matrix. For the complex version. */
    const model_complex_sparse_matrix &complex_tangent_matrix(void) const {
      GMM_ASSERT1(complex_version, "This model is a real one");
      context_check(); if (act_size_to_be_done) actualize_sizes();
      return cTM;
    }

//     template<typename MATRIX, typename T>
//     const MATRIX &tangent_matrix(T) const { return real_rhs(); }

//     template<typename MATRIX, typename T>
//     const MATRIX &tangent_matrix(std::complex<T>) const
//     { return complex_rhs(); }
    

//     template<typename MATRIX> const MATRIX &tangent_matrix(void) const {
//       typedef typename gmm::linalg_traits<MATRIX>::value_type T;
//       return tangent_matrix<MATRIX>(T());
//     }

    
    /** Gives the access to the right hand side of the tangent linear system.
	For the real version. */
    const model_real_plain_vector &real_rhs(void) const {
      GMM_ASSERT1(!complex_version, "This model is a complex one");
      context_check(); if (act_size_to_be_done) actualize_sizes();
      return rrhs;
    }
    
    /** Gives the access to the right hand side of the tangent linear system.
	For the complex version. */
    const model_complex_plain_vector &complex_rhs(void) const {
      GMM_ASSERT1(complex_version, "This model is a real one");
      context_check(); if (act_size_to_be_done) actualize_sizes();
      return crhs;
    }

//     template<typename VECTOR, typename T>
//     const VECTOR &rhs(T) const { return real_rhs(); }

//     template<typename VECTOR, typename T>
//     const VECTOR &rhs(std::complex<T>) const { return complex_rhs(); }
    

//     template<typename VECTOR> const VECTOR &rhs(void) const {
//       typedef typename gmm::linalg_traits<VECTOR>::value_type T;
//       return rhs<VECTOR>(T());
//     }

    /** List the model variables and constant. */
    void listvar(std::ostream &ost) const;

    /** List the model bricks. */
    void listbricks(std::ostream &ost) const;

    /** Force the re-computation of a brick for the next assembly. */ 
    void touch_brick(size_type ind_brick) {
      GMM_ASSERT1(ind_brick < bricks.size(), "Inexistent brick");
      bricks[ind_brick].terms_to_be_computed = true;
    }

    pbrick brick_pointer(size_type ind_brick) {
      GMM_ASSERT1(ind_brick < bricks.size(), "Inexistent brick");
      return bricks[ind_brick].pbr;
    }

    /** Add a brick to the model. varname is the list of variable used
        and datanames the data used. If a variable is used as a data, it
        should be declared in the datanames (it will depend on the value of
	the variable not only on the fem). */
    size_type add_brick(pbrick pbr, const varnamelist &varnames,
			const varnamelist &datanames,
			const termlist &terms, const mimlist &mims, 
			size_type region);

    /** Gives the name of the variable of index `ind_var` of the brick
	of index `ind_brick`. */
    const std::string &varname_of_brick(size_type ind_brick,
					size_type ind_var);

    /** Gives the name of the data of index `ind_data` of the brick
	of index `ind_brick`. */
    const std::string &dataname_of_brick(size_type ind_brick,
					 size_type ind_data);

    enum assembly_version { BUILD_RHS = 1, BUILD_MATRIX = 2, BUILD_ALL = 3 };

    /** Assembly of the tangent system taking into account the terms
	from all bricks. */
    void assembly(assembly_version version);

    void clear(void) {
      variables.clear();
      rTM = model_real_sparse_matrix();
      cTM = model_complex_sparse_matrix();
      rrhs = model_real_plain_vector();
      crhs = model_complex_plain_vector();
    }

    model(bool comp_version = false) {
      init(); complex_version = comp_version;
      is_linear_ = is_symmetric_ = is_coercive_ = true;
      leading_dim = 0;
    }

  };

  //=========================================================================
  //
  //  Brick object.
  //
  //=========================================================================

  // Si la brique a besoin d'un stockage "interne" elle doit le faire dans
  // le modèle en reservant une variable ? Si le besoin est une matrice
  //  creuse ? Prévoir qlq chose de spécifique dans le modèle ? (exemple de la
  // brique contraite ...)
  // A voir selon les besoins futurs.

  /** The virtual brick has to be derived to describe real model bricks.
      The set_flags method has to be called by the derived class.
      The virtual methods asm_real_tangent_terms and/or
      asm_complex_tangent_terms have to be defined.
      The brick should not store data. The data have to be stored in the
      model object. 
  **/
  class virtual_brick : virtual public dal::static_stored_object {
  private :
    bool islinear;    // The brick add a linear term or not.
    bool issymmetric; // The brick add a symmetric term or not.
    bool iscoercive;  // The brick add a potentialy coercive terms or not. 
    //   (in particular, not a term involving a multiplier)
    bool isreal;      // The brick admits a real version or not.
    bool iscomplex;   // The brick admits a complex version or not.
    bool isinit;      // internal flag.
    std::string name; // Name of the brick.
   
  public :

    typedef model::assembly_version nonlinear_version;
    
    virtual_brick(void) { isinit = false; }
    void set_flags(const std::string &bname, bool islin, bool issym,
		   bool iscoer, bool ire, bool isco) {
      name = bname;
      islinear = islin; issymmetric = issym; iscoercive = iscoer;
      isreal = ire; iscomplex = isco; isinit = true;
    }

#   define BRICK_NOT_INIT GMM_ASSERT1(isinit, "Set brick flags !")
    bool is_linear(void)    const { BRICK_NOT_INIT; return islinear;    }
    bool is_symmetric(void) const { BRICK_NOT_INIT; return issymmetric; }
    bool is_coercive(void)  const { BRICK_NOT_INIT; return iscoercive;  }
    bool is_real(void)      const { BRICK_NOT_INIT; return isreal;      }
    bool is_complex(void)   const { BRICK_NOT_INIT; return iscomplex;   }
    const std::string &brick_name(void) const { BRICK_NOT_INIT; return name; }

    virtual void asm_real_tangent_terms(const model &,
					const model::varnamelist &,
					const model::varnamelist &,
					const model::mimlist &,
					model::real_matlist &,
					model::real_veclist &,
					size_type, nonlinear_version) const
    { GMM_ASSERT1(false, "Brick has no real tangent terms !"); }

    virtual void asm_complex_tangent_terms(const model &,
					   const model::varnamelist &,
					   const model::varnamelist &,
					   const model::mimlist &,
					   model::complex_matlist &,
					   model::complex_veclist &,
					   size_type, nonlinear_version) const
    { GMM_ASSERT1(false, "Brick has no complex tangent terms !"); }
    
  };

  //=========================================================================
  //
  //  Functions adding standard bricks to the model.
  //
  //=========================================================================

  
  /** Add a Laplacian term on the variable `varname`. If it is a vector
      valued variable, the Laplacian term is componentwise. `region` is an
      optional mesh region on which the term is added. Return the brick index
      in the model.
  */
  size_type add_Laplacian_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   size_type region = size_type(-1));
  

  /** Add an elliptic term on the variable `varname`. The shape of the elliptic
      term depends both on the variable and the data. This corresponds to a
      term $-\text{div}(a\nabla u)$ where $a$ is the data and $u$ the variable.
      The data can be a scalar, a matrix or an order four tensor. The variable
      can be vector valued or not. If the data is a scalar or a matrix and
      the variable is vector valued then the term is added componentwise.
      An order four tensor data is allowed for vector valued variable only.
      The data can be constant or describbed on a fem. Of course, when
      the data is a tensor describe on a finite element method (a tensor
      field) the data can be a huge vector. The components of the
      matrix/tensor have to be stored with the fortran order (columnwise) in
      the data vector (compatibility with blas). The symmetry and coercivity
      of the given matrix/tensor is not verified (but assumed). `region` is an
      optional mesh region on which the term is added. Return the brick index
      in the model.
  */
  size_type add_generic_elliptic_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname, size_type region = size_type(-1));
  

  /** Add a source term on the variable `varname`. The source term is
      represented by the data `dataname` which could be constant or described
      on a fem.  `region` is an optional mesh region on which the term is
      added. An additional optional data `directdataname` can be provided. The
      corresponding data vector will be directly added to the right hand
      side without assembly. Return the brick index in the model.
  */
  size_type add_source_term_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname, size_type region = size_type(-1),
   const std::string &directdataname = std::string());

  /** Add a source term on the variable `varname` on a boundary `region`.
      The source term is
      represented by the data `dataname` which could be constant or described
      on a fem. A sclar product with the outward normal unit vector to
      the boundary is performed. The main aim of this brick is to represent
      a Neumann condition with a vector data without performing the
      scalar product with the normal as a pre-processing.
  */
  size_type add_normal_source_term_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname, size_type region);

  /** Add a Dirichlet condition on the variable `varname` and the mesh
      region `region`. This region should be a boundary. The Dirichlet
      condition is prescribed with a multiplier variable `multname` which
      should be first declared as a multiplier 
      variable on the mesh region in the model. `dataname` is the optional
      right hand side of  the Dirichlet condition. It could be constant or
      described on a fem; scalar or vector valued, depending on the variable
      on which the Dirichlet condition is prescribed. Return the brick index
      in the model.
  */
  size_type add_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &multname, size_type region,
   const std::string &dataname = std::string());

  /** Same function as the preceeding one but the multipliers variable will
      be declared to the brick by the function. `mf_mult` is the finite element
      method on which the multiplier will be build (it will be restricted to
      the mesh region `region` and eventually some conflicting dofs with some
      other multiplier variables will be suppressed).
  */
  size_type add_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   const mesh_fem &mf_mult, size_type region,
   const std::string &dataname = std::string());

  /** Same function as the preceeding one but the `mf_mult` parameter is
      replaced by `degree`. The multiplier will be described on a standard
      finite element method of the corresponding degree.
   */
  size_type add_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   dim_type degree, size_type region,
   const std::string &dataname = std::string());

  /** When `ind_brick` is the index of a Dirichlet brick with multiplier on
      the model `md`, the function return the name of the multiplier variable.
      Otherwise, it has an undefined behavior.
  */
  const std::string &mult_varname_Dirichlet(model &md, size_type ind_brick);

  /** Add a Dirichlet condition on the variable `varname` and the mesh
      region `region`. This region should be a boundary. The Dirichlet
      condition is prescribed with penalization. The penalization coefficient
      is intially `penalization_coeff` and will be added to the data of
      the model. `dataname` is the optional
      right hand side of  the Dirichlet condition. It could be constant or
      described on a fem; scalar or vector valued, depending on the variable
      on which the Dirichlet condition is prescribed. Return the brick index
      in the model.
  */
  size_type add_Dirichlet_condition_with_penalization
  (model &md, const mesh_im &mim, const std::string &varname,
   scalar_type penalization_coeff, size_type region,
   const std::string &dataname = std::string());

  /** Change the penalization coefficient of a Dirichlet condition with
      penalization brick. If the brick is not of this kind,
      this function has an undefined behavior.
  */
  void change_penalization_coeff(model &md, size_type ind_brick,
				 scalar_type penalisation_coeff); 

  /** Add a Helmoltz brick to the model. This corresponds to the scalar
      equation (@f$\Delta u + k^2u = 0@f$, with @f$K=k^2@f$).
      The weak formulation is (@f$\int k^2 u.v - \nabla u.\nabla v@f$)

      `dataname` should contain the wave number $k$. It can be real or
      complex, fem dependant or not.
  */
  size_type add_Helmholtz_brick(model &md, const mesh_im &mim,
				const std::string &varname,
				const std::string &dataname,
				size_type region = size_type(-1));


  /** Add a Fourier-Robin brick to the model. This correspond to the weak term
      (@f$\int (qu).v @f$) on a boundary. It is used to represent a
      Fourier-Robin boundary condition.

      `dataname` should contain the parameter $q$ which should be a
      (@f$N\times N@f$) matrix term, where $N$ is the dimension of the
      variable `varname`. Note that an additional right hand
      side can be added with a source term brick.
  */
  size_type add_Fourier_Robin_brick(model &md, const mesh_im &mim,
				    const std::string &varname,
				    const std::string &dataname,
				    size_type region);

  // Constraint brick.
  model_real_sparse_matrix &set_private_data_brick_real_matrix
  (model &md, size_type indbrick);
  model_real_plain_vector &set_private_data_brick_real_rhs
  (model &md, size_type indbrick);
  model_complex_sparse_matrix &set_private_data_brick_complex_matrix
  (model &md, size_type indbrick);
  model_complex_plain_vector &set_private_data_brick_complex_rhs
  (model &md, size_type indbrick);
  size_type add_constraint_with_penalization
  (model &md, const std::string &varname, scalar_type penalisation_coeff);
  size_type add_constraint_with_multipliers
  (model &md, const std::string &varname, const std::string &multname);

  template <typename VECT, typename T>
  void set_private_data_rhs(model &md, size_type ind,
				const VECT &L, T) {
    model_real_plain_vector &LL = set_private_data_brick_real_rhs(md, ind);
    gmm::resize(LL, gmm::vect_size(L));
    gmm::copy(L, LL);
  }

  template <typename VECT, typename T>
  void set_private_data_rhs(model &md, size_type ind, const VECT &L,
			   std::complex<T>) {
    model_complex_plain_vector &LL = set_private_data_brick_complex_rhs(md, ind);
    gmm::resize(LL, gmm::vect_size(L));
    gmm::copy(L, LL);
  }

  /** For some specific bricks having an internal right hand side vector
      (explicit bricks: 'constraint brick' and 'explicit rhs brick'),
      set this rhs. 
  */
  template <typename VECT>
  void set_private_data_rhs(model &md, size_type indbrick, const VECT &L) {
    typedef typename gmm::linalg_traits<VECT>::value_type T;
    set_private_data_rhs(md, indbrick, L, T());
  }

  template <typename MAT, typename T>
  void set_private_data_matrix(model &md, size_type ind,
				   const MAT &B, T) {
    model_real_sparse_matrix &BB = set_private_data_brick_real_matrix(md, ind);
    gmm::resize(BB, gmm::mat_nrows(B), gmm::mat_ncols(B));
    gmm::copy(B, BB);
  }

  template <typename MAT, typename T>
  void set_private_data_matrix(model &md, size_type ind, const MAT &B,
			      std::complex<T>) {
    model_complex_sparse_matrix &BB
      = set_private_data_brick_complex_matrix(md, ind);
    gmm::resize(BB, gmm::mat_nrows(B), gmm::mat_ncols(B));
    gmm::copy(B, BB);
  }

  /** For some specific bricks having an internal sparse matrix
      (explicit bricks: 'constraint brick' and 'explicit matrix brick'),
      set this matrix. @*/
  template <typename MAT>
  void set_private_data_matrix(model &md, size_type indbrick,
				   const MAT &B) {
    typedef typename gmm::linalg_traits<MAT>::value_type T;
    set_private_data_matrix(md, indbrick, B, T());
  }

  /** Add an additional explicit penalized constraint on the variable
      `varname`. The constraint is $BU=L$ with `B` being a rectangular
      sparse matrix.
      Be aware that `B` should not contain a palin row, otherwise the whole
      tangent matrix will be plain. It is possible to change the constraint
      at any time whith the methods set_private_matrix and set_private_rhs.
      The method change_penalization_coeff can also be used.
  */
  template <typename MAT, typename VECT>
  size_type add_constraint_with_penalization
  (model &md, const std::string &varname, scalar_type penalisation_coeff,
   const MAT &B, const VECT &L) {
    size_type ind
      = add_constraint_with_penalization(md, varname, penalisation_coeff);
    size_type n = gmm::mat_nrows(B), m = gmm::mat_ncols(B);
    set_private_data_rhs(md, ind, L);
    set_private_data_matrix(md, ind, B);
    return ind;
  }

  /** Add an additional explicit constraint on the variable `varname` thank to
    a multiplier `multname` peviously added to the model (should be a fixed
    size variable).
    The constraint is $BU=L$ with `B` being a rectangular sparse matrix.
    It is possible to change the constraint
    at any time whith the methods set_private_matrix
    and set_private_rhs.
  */
  template <typename MAT, typename VECT>
  size_type add_constraint_with_multipliers
  (model &md, const std::string &varname, const std::string &multname,
   const MAT &B, const VECT &L) {
    size_type ind = add_constraint_with_multipliers(md, varname, multname);
    set_private_data_rhs(md, ind, L);
    set_private_data_matrix(md, ind, B);
    return ind;
  }

  size_type add_explicit_matrix(model &md, const std::string &varname1,
				const std::string &varname2,
				bool issymmetric, bool iscoercive); 
  size_type add_explicit_rhs(model &md, const std::string &varname);
  
  /** Add a brick reprenting an explicit matrix to be added to the tangent
      linear system relatively to the variables 'varname1' and 'varname2'.
      The given matrix should have has many rows as the dimension of
      'varname1' and as many columns as the dimension of 'varname2'.
      If the two variables are different and if `issymmetric' is set to true
      then the transpose of the matrix is also added to the tangent system
      (default is false). set `iscoercive` to true if the term does not
      affect the coercivity of the tangent system (default is false).
      The matrix can be changed by the command set_private_matrix.
  */
  template <typename MAT>
  size_type add_explicit_matrix(model &md, const std::string &varname1,
				const std::string &varname2, const MAT &B,
				bool issymmetric = false,
				bool iscoercive = false) {
    size_type ind = add_explicit_matrix(md, varname1, varname2,
					issymmetric, iscoercive);
    set_private_data_matrix(md, ind, B);
    return ind;
  }

  /**  Add a brick reprenting an explicit right hand side to be added to
       the right hand side of the tangent
       linear system relatively to the variable 'varname'.
       The given rhs should have the same size than the dimension of
       'varname'. The rhs can be changed by the command set_private_rhs.
  */  
  template <typename VECT>
  size_type add_explicit_rhs(model &md, const std::string &varname,
			     const VECT &L) {
    size_type ind = add_explicit_rhs(md, varname);
    set_private_data_rhs(md, ind, L);
    return ind;
  }
  


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MODELS_H__  */
