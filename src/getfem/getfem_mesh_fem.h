// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 1999-2009 Yves Renard
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

/**@file getfem_mesh_fem.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date December 21, 1999.
   @brief Define the getfem::mesh_fem class
*/

#ifndef GETFEM_MESH_FEM_H__
#define GETFEM_MESH_FEM_H__

#include "getfem_mesh.h"
#include "getfem_fem.h"


namespace getfem {

  template <class CONT> struct tab_scal_to_vect_iterator {

    typedef typename CONT::const_iterator ITER;
    typedef typename std::iterator_traits<ITER>::value_type value_type;
    typedef typename std::iterator_traits<ITER>::pointer    pointer;
    typedef typename std::iterator_traits<ITER>::reference  reference;
    typedef typename std::iterator_traits<ITER>::difference_type
                                                            difference_type;
    typedef typename std::iterator_traits<ITER>::iterator_category
                                                            iterator_category;
    typedef size_t size_type;
    typedef tab_scal_to_vect_iterator<CONT> iterator;

    ITER it;
    dim_type N;
    dim_type ii;

    iterator &operator ++()
    { ++ii; if (ii == N) { ii = 0; ++it; } return *this; }
    iterator &operator --() 
    { if (ii == 0) { ii = N-1; --it; } else --ii; return *this; }
    iterator operator ++(int) { iterator tmp = *this; ++(*this); return tmp; }
    iterator operator --(int) { iterator tmp = *this; --(*this); return tmp; }
   
    iterator &operator +=(difference_type i)
    { it += (i+ii)/N; ii = dim_type((ii + i) % N); return *this; }
    iterator &operator -=(difference_type i)
    { it -= (i+N-ii-1)/N; ii = (ii - i + N * i) % N; return *this; }
    iterator operator +(difference_type i) const 
    { iterator itt = *this; return (itt += i); }
    iterator operator -(difference_type i) const
    { iterator itt = *this; return (itt -= i); }
    difference_type operator -(const iterator &i) const
    { return (it - i.it) * N + ii - i.ii; }

    value_type operator *() const { return (*it) + ii; }
    value_type operator [](int i) { return *(this + i); }

    bool operator ==(const iterator &i) const
    { return (it == i.it) && (ii == i.ii); }
    bool operator !=(const iterator &i) const { return !(i == *this); }
    bool operator < (const iterator &i) const
    { return (it < i.it) && (ii < i.ii); }

    tab_scal_to_vect_iterator(void) {}
    tab_scal_to_vect_iterator(const ITER &iter, dim_type n, dim_type i)
      : it(iter), N(n), ii(i) { }

  };

  /** @internal @brief structure for iteration over the dofs when Qdim
      != 1 and target_dim == 1
  */
  template <class CONT> class tab_scal_to_vect {
  public :
    typedef typename CONT::const_iterator ITER;
    typedef typename std::iterator_traits<ITER>::value_type value_type;
    typedef typename std::iterator_traits<ITER>::pointer    pointer;
    typedef typename std::iterator_traits<ITER>::pointer    const_pointer;
    typedef typename std::iterator_traits<ITER>::reference  reference;
    typedef typename std::iterator_traits<ITER>::reference  const_reference;
    typedef typename std::iterator_traits<ITER>::difference_type
            difference_type;
    typedef size_t size_type;
    typedef tab_scal_to_vect_iterator<CONT> iterator;
    typedef iterator                          const_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;


  protected :
    ITER it, ite;
    dim_type N;
    
  public :

    bool empty(void) const { return it == ite; }
    size_type size(void) const { return (ite - it) * N; }

    const_iterator begin(void) const { return iterator(it, N, 0); }
    const_iterator end(void) const { return iterator(ite, N, 0); }
    const_reverse_iterator rbegin(void) const
    { return const_reverse_iterator(end()); }
    const_reverse_iterator rend(void) const
    { return const_reverse_iterator(begin()); }
    
    value_type front(void) const { return *begin(); }
    value_type back(void) const { return *(--(end())); }

    tab_scal_to_vect(void) {}
    tab_scal_to_vect(const CONT &cc, dim_type n)
      : it(cc.begin()), ite(cc.end()), N(n) {}
    
    value_type operator [](size_type ii) const { return *(begin() + ii);}
  };
  
  /** Describe a finite element method linked to a mesh.
   *    
   *  @see mesh
   *  @see mesh_im
   */
  class mesh_fem : public context_dependencies {
  protected :
    
    typedef gmm::csc_matrix<scalar_type> REDUCTION_MATRIX;
    typedef gmm::csr_matrix<scalar_type> EXTENSION_MATRIX;
    
    dal::dynamic_array<pfem> f_elems;
    dal::bit_vector fe_convex;
    const mesh *linked_mesh_;
    REDUCTION_MATRIX R_;
    EXTENSION_MATRIX E_;
    mutable bgeot::mesh_structure dof_structure;
    mutable bool dof_enumeration_made;
    mutable size_type nb_total_dof;
    pfem auto_add_elt_pf; /* fem for automatic addition                   */
                       /* of element option. (0 = no automatic addition)  */
    dim_type auto_add_elt_K; /* Degree of the fem for automatic addition  */
                       /* of element option. (-1 = no automatic addition) */
    dim_type Qdim; /* this is the "global" target_dim */
    dim_type QdimM, QdimN; /* for matrix field with QdimM lines and QdimN */
                           /* columnsQdimM * QdimN = Qdim.                */
    std::vector<size_type> dof_partition;
    mutable gmm::uint64_type v_num_update, v_num;
    bool use_reduction;    /* A reduction matrix is applied or not.       */
    
  public :
    typedef base_node point_type;
    typedef tab_scal_to_vect<mesh::ind_cv_ct> ind_dof_ct;
    typedef tab_scal_to_vect<mesh::ind_pt_face_ct> ind_dof_face_ct;

    void update_from_context(void) const;

    gmm::uint64_type version_number(void) const
    { context_check(); return v_num; }

    /** Get the set of convexes where a finite element has been assigned.
     */
    inline const dal::bit_vector &convex_index(void) const
    { context_check(); return fe_convex; }
    
    /// Return true if a reduction matrix is applied to the dofs.
    bool is_reduced(void) const { return use_reduction; }
    
    /// Return the reduction matrix applied to the dofs.
    const REDUCTION_MATRIX &reduction_matrix(void) const { return R_; }

    /// Return the extension matrix corresponding to reduction applied (RE=I).
    const EXTENSION_MATRIX &extension_matrix(void) const { return E_; }
    
    /** Allows to set the reduction and the extension matrices.
     * Should satify (RR*EE=I). */
    template <typename MATR, typename MATE>
    void set_reduction_matrices(const MATR &RR, const MATE &EE) {
      context_check();
      GMM_ASSERT1(gmm::mat_ncols(RR) == nb_basic_dof() &&
		  gmm::mat_nrows(EE) == nb_basic_dof() &&
		  gmm::mat_nrows(RR) == gmm::mat_ncols(EE),
		  "Wrong dimension of reduction and/or extension matrices");
      R_ = REDUCTION_MATRIX(gmm::mat_nrows(RR), gmm::mat_ncols(RR));
      E_ = EXTENSION_MATRIX(gmm::mat_nrows(EE), gmm::mat_ncols(EE));
      gmm::copy(RR, R_);
      gmm::copy(EE, E_);
      use_reduction = true;
      touch(); v_num = act_counter();
    }

    /** Allows to set the reduction and the extension matrices in order
     *  to keep only a certain number of dof. */
    void reduce_to_basic_dof(const dal::bit_vector &kept_basic_dof);
    void reduce_to_basic_dof(const std::set<size_type> &kept_basic_dof);

    /// Validate or invalidate the reduction (keeping the reduction matrices).
    void set_reduction(bool r) {
      if (r != use_reduction) {
	use_reduction = r;
	if (use_reduction) {
	  context_check();
	  GMM_ASSERT1(gmm::mat_ncols(R_) == nb_basic_dof() &&
		      gmm::mat_nrows(E_) == nb_basic_dof() &&
		      gmm::mat_nrows(R_) == gmm::mat_ncols(E_),
		    "Wrong dimension of reduction and/or extension matrices");
	}
	touch(); v_num = act_counter();
      }
    }

    template<typename VECT1, typename VECT2>
    void reduce_vector(const VECT1 &V1, const VECT2 &V2) const {
      if (is_reduced()) {
	size_type qqdim = gmm::vect_size(V1) / nb_basic_dof();
	if (qqdim == 1)
	  gmm::mult(reduction_matrix(), V1, const_cast<VECT2 &>(V2));
	else
	  for (size_type k = 0; k < qqdim; ++k)
	    gmm::mult(reduction_matrix(),
		      gmm::sub_vector(V1, gmm::sub_slice(k, nb_basic_dof(),
							 qqdim)),
		      gmm::sub_vector(const_cast<VECT2 &>(V2),
				      gmm::sub_slice(k, nb_dof(),
						     qqdim)));
      }
      else gmm::copy(V1, const_cast<VECT2 &>(V2));
    }

    template<typename VECT1, typename VECT2>
    void extend_vector(const VECT1 &V1, const VECT2 &V2) const {
      if (is_reduced()) {
	size_type qqdim = gmm::vect_size(V1) / nb_dof();
	if (qqdim == 1)
	  gmm::mult(extension_matrix(), V1, const_cast<VECT2 &>(V2));
	else
	  for (size_type k = 0; k < qqdim; ++k)
	    gmm::mult(extension_matrix(),
		      gmm::sub_vector(V1, gmm::sub_slice(k, nb_dof(),
							 qqdim)),
		      gmm::sub_vector(const_cast<VECT2 &>(V2),
				      gmm::sub_slice(k, nb_basic_dof(),
						     qqdim)));
      }
      else gmm::copy(V1, const_cast<VECT2 &>(V2));
    }
    

    /// Return a reference to the underlying mesh.
    const mesh &linked_mesh(void) const { return *linked_mesh_; }

    /** Set the degree of the fem for automatic addition
     *  of element option. K=-1 disables the automatic addition.
     */
    void set_auto_add(dim_type K) { auto_add_elt_K = K; auto_add_elt_pf = 0; }

    /** Set the fem for automatic addition
     *  of element option. pf=0 disables the automatic addition.
     */
    void set_auto_add(pfem pf)
    { auto_add_elt_pf = pf; auto_add_elt_K = dim_type(-1);}

    /** Return the Q dimension. A mesh_fem used for scalar fields has
	Q=1, for vector fields, Q is typically equal to
	linked_mesh().dim().
     */
    dim_type get_qdim() const { return Qdim; }
    /** Change the Q dimension */
    void set_qdim(dim_type q) {
      if (q != Qdim || q != QdimM) {
	QdimM = Qdim = q; QdimN = 1;
	dof_enumeration_made = false; touch(); v_num = act_counter();
      }
    }

    /** Set the dimension for a matrix field. The Q dimension will
	always be the product of get_qdim_m and get_qdim_n */
    void set_qdim_mn(dim_type M, dim_type N) {
      if (M != QdimM || N != QdimN) {
	QdimM = M; QdimN = N; Qdim = dim_type(N*M);
	dof_enumeration_made = false; touch(); v_num = act_counter();
      }
    }

    /** for matrix fields, return the number of rows. */
    dim_type get_qdim_m() const { return QdimM; }
    /** for matrix fields, return the number of columns. */
    dim_type get_qdim_n() const { return QdimN; }
     

    /** Set the finite element method of a convex.
	@param cv the convex number.
	@param pf the finite element.
    */
    void set_finite_element(size_type cv, pfem pf);	
    /** Set the finite element on a set of convexes.
	@param cvs the set of convex indexes, as a dal::bit_vector.
	@param pf the finite element, typically obtained with
	@code getfem::fem_descriptor("FEM_SOMETHING(..)") 
	@endcode
    */
    void set_finite_element(const dal::bit_vector &cvs, pfem pf);
    /** shortcut for set_finite_element(linked_mesh().convex_index(),pf);
        and set_auto_add(pf). */
    void set_finite_element(pfem pf);
    /** Set a classical (i.e. lagrange polynomial) finite element on
	a convex.
	@param cv is the convex number.
	@param fem_degree the polynomial degree of the finite element.
    */
    void set_classical_finite_element(size_type cv, dim_type fem_degree);
    /** Set a classical (i.e. lagrange polynomial) finite element on
	a set of convexes.
	@param cvs the set of convexes, as a dal::bit_vector.
	@param fem_degree the polynomial degree of the finite element.
    */
    void set_classical_finite_element(const dal::bit_vector &cvs, 
				      dim_type fem_degree);
    /** Similar to set_classical_finite_element, but uses
	discontinuous lagrange elements. 

	@param cvs the set of convexes, as a dal::bit_vector.
	@param fem_degree the polynomial degree of the finite element.
	@param alpha the node inset, 0 <= alpha < 1, where 0 implies
	usual dof nodes, greater values move the nodes toward the
	center of gravity, and 1 means that all degrees of freedom
	collapse on the center of gravity.
     */
    void set_classical_discontinuous_finite_element(const dal::bit_vector &cvs,
						    dim_type fem_degree,
						    scalar_type alpha=0);
    /** Shortcut for
     * set_classical_finite_element(linked_mesh().convex_index(),...)
     */
    void set_classical_finite_element(dim_type fem_degree);
    /** Shortcut for
     *  set_classical_discontinuous_finite_element(linked_mesh().convex_index()
     *  ,...)
     */
    void set_classical_discontinuous_finite_element(dim_type fem_degree,
						    scalar_type alpha=0);   
    /** Return the basic fem associated with an element (in no fem is
     *	associated, the function will crash! use the convex_index() of
     *  the mesh_fem to check that a fem is associated to a given
     *  convex). This fem does not take into account the optional
     *  vectorization due to qdim nor the optional reduction.
     */
    virtual pfem fem_of_element(size_type cv) const
    { return  f_elems[cv]; }
    /** Give an array of the dof numbers a of convex.
     *  @param cv the convex number.
     *  @return a pseudo-container of the dof number.
     */
    virtual ind_dof_ct ind_basic_dof_of_element(size_type cv) const {
      context_check(); if (!dof_enumeration_made) enumerate_dof();
      return ind_dof_ct(dof_structure.ind_points_of_convex(cv),
			dim_type(Qdim/fem_of_element(cv)->target_dim()));
    }
    ind_dof_ct ind_dof_of_element(size_type cv) const IS_DEPRECATED
    { return ind_basic_dof_of_element(cv); }
    /** Give an array of the dof numbers lying of a convex face (all
      	degrees of freedom whose associated base function is non-zero
	on the convex face).
	@param cv the convex number.
	@param f the face number.
	@return a pseudo-container of the dof number.
    */
    virtual ind_dof_face_ct
    ind_basic_dof_of_face_of_element(size_type cv, short_type f) const {
      context_check(); if (!dof_enumeration_made) enumerate_dof();
      return ind_dof_face_ct
	(dof_structure.ind_points_of_face_of_convex(cv, f),
	 dim_type(Qdim/fem_of_element(cv)->target_dim()));
    }
    ind_dof_face_ct
    ind_dof_of_face_of_element(size_type cv,short_type f) const IS_DEPRECATED
    { return ind_basic_dof_of_face_of_element(cv, f); }
    /** Return the number of dof lying on the given convex face.
	@param cv the convex number.
	@param f the face number.
    */
    virtual size_type nb_basic_dof_of_face_of_element(size_type cv,
					      short_type f) const {
      pfem pf = f_elems[cv];
      return dof_structure.structure_of_convex(cv)->nb_points_of_face(f)
	* Qdim / pf->target_dim();
    }
    size_type nb_dof_of_face_of_element(size_type cv,
					short_type f) const IS_DEPRECATED
    { return nb_basic_dof_of_face_of_element(cv, f); }
    /** Return the number of  degrees of freedom attached to a given convex.
	@param cv the convex number.
    */
    virtual size_type nb_basic_dof_of_element(size_type cv) const
    { pfem pf = f_elems[cv]; return pf->nb_dof(cv) * Qdim / pf->target_dim(); }
    size_type nb_dof_of_element(size_type cv) const IS_DEPRECATED
    { return nb_basic_dof_of_element(cv); }
    
    /* Return the geometrical location of a degree of freedom in the
       reference convex.
	@param cv the convex number.
	@param i the local dof number.
    const base_node &reference_point_of_dof(size_type cv,size_type i) const {
      pfem pf = f_elems[cv];
      return pf->node_of_dof(cv, i * pf->target_dim() / Qdim);
    }
    */
    /** Return the geometrical location of a degree of freedom.
	@param cv the convex number.
	@param i the local dof number.
     */
    virtual base_node point_of_basic_dof(size_type cv, size_type i) const;
    base_node point_of_dof(size_type cv, size_type i) const IS_DEPRECATED
    { return point_of_basic_dof(cv, i); }
    /** Return the geometrical location of a degree of freedom.
	@param d the global dof number.
    */
    virtual base_node point_of_basic_dof(size_type d) const;
    base_node point_of_dof(size_type d) const IS_DEPRECATED
    { return point_of_basic_dof(d); }
    /** Return the dof component number (0<= x <Qdim) */
    virtual dim_type basic_dof_qdim(size_type d) const;
    dim_type dof_qdim(size_type d) const IS_DEPRECATED
    { return basic_dof_qdim(d); }
    /** Shortcut for convex_to_dof(d)[0] 
 	@param d the global dof number.
    */
    virtual size_type first_convex_of_basic_dof(size_type d) const;
    size_type first_convex_of_dof(size_type d) const IS_DEPRECATED
    { return first_convex_of_basic_dof(d); }
    /** Return the list of convexes attached to the specified dof 
	@param d the global dof number.
	@return an array of convex numbers.
     */
    virtual const mesh::ind_cv_ct &convex_to_basic_dof(size_type d) const;
    const mesh::ind_cv_ct &convex_to_dof(size_type d) const IS_DEPRECATED
    { return convex_to_basic_dof(d); }
    /** Renumber the degrees of freedom. You should not have
     * to call this function, as it is done automatically */
    void enumerate_dof(void) const;
    /** Return the total number of basic degrees of freedom (before the
     * optional redution). */
    virtual size_type nb_basic_dof(void) const {
      context_check(); if (!dof_enumeration_made) enumerate_dof();
      return nb_total_dof;
    }
    /// Return the total number of degrees of freedom.
    virtual size_type nb_dof(void) const {
      context_check(); if (!dof_enumeration_made) enumerate_dof();
      return use_reduction ? gmm::mat_nrows(R_) : nb_total_dof;
    }
    /** Get a list of basic dof lying on a given mesh_region.
	@param b the mesh_region.
	@return the list in a dal::bit_vector.
    */
    virtual dal::bit_vector basic_dof_on_region(const mesh_region &b) const;
    /** Get a list of dof lying on a given mesh_region. For a reduced mesh_fem
	a dof is lying on a region if its potential corresponding shape
	function is nonzero on this region. The extension matrix is used
	to make the correspondance between basic and reduced dofs.
	@param b the mesh_region.
	@return the list in a dal::bit_vector.
    */    
    dal::bit_vector dof_on_region(const mesh_region &b) const;
    dal::bit_vector dof_on_set(const mesh_region &b) const IS_DEPRECATED
    { return dof_on_region(b); }

    void set_dof_partition(size_type cv, unsigned partition_num) {
      if (dof_partition.empty() && partition_num == 0) return;      
      if (dof_partition.size() < linked_mesh().convex_index().last_true()+1) 
	dof_partition.resize(linked_mesh().convex_index().last_true()+1);
      if (dof_partition.at(cv) != partition_num) {
	dof_partition[cv] = partition_num;
	dof_enumeration_made = false;
      }
    }
    unsigned get_dof_partition(size_type cv) const {
      return (cv < dof_partition.size() ? unsigned(dof_partition[cv]) : 0); 
    }
    void clear_dof_partition() { dof_partition.clear(); }
    
    size_type memsize() const {
      return dof_structure.memsize() + 
	sizeof(mesh_fem) - sizeof(bgeot::mesh_structure) +
	f_elems.memsize() + fe_convex.memsize();
    }
    /** Build a new mesh_fem. A mesh object must be supplied. 
	@param me the linked mesh.
	@param Q the Q dimension (see mesh_fem::get_qdim).
    */
    explicit mesh_fem(const mesh &me, dim_type Q = 1);
    virtual ~mesh_fem();
    virtual void clear(void);
    /** Read the mesh_fem from a stream. 
	@param ist the stream.
     */
    virtual void read_from_file(std::istream &ist);
    /** Read the mesh_fem from a file.
        @param name the file name. */
    void read_from_file(const std::string &name);
    /* internal usage. */
    void write_basic_to_file(std::ostream &ost) const;
    /* internal usage. */
    void write_reduction_matrices_to_file(std::ostream &ost) const;
    /** Write the mesh_fem to a stream. */
    virtual void write_to_file(std::ostream &ost) const;
    /** Write the mesh_fem to a file. 

	@param name the file name

	@param with_mesh if set, then the linked_mesh() will also be
	saved to the file.
    */
    void write_to_file(const std::string &name, bool with_mesh=false) const;
  };

  /** Gives the descriptor of a classical finite element method of degree K
      on mesh.  
      
      The mesh_fem won't be destroyed until its linked_mesh is
      destroyed. All the mesh_fem built by this function are stored
      in a cache, which means that calling this function twice with
      the same arguments will return the same mesh_fem object. A
      consequence is that you should NEVER modify this mesh_fem!
   */
  const mesh_fem &classical_mesh_fem(const mesh &mesh, dim_type degree);

  /* Dummy mehs_fem for default parameter of functions. */
  const mesh_fem &dummy_mesh_fem(void);


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MESH_FEM_H__  */
