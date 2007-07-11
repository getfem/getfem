// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 1999-2007 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

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
      { it += (i+ii)/N; ii = (ii + i) % N; return *this; }
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
  class mesh_fem : public mesh_receiver, public context_dependencies {
  protected :
    
    dal::dynamic_array<pfem> f_elems;
    dal::bit_vector fe_convex;
    const mesh *linked_mesh_;
    mutable bgeot::mesh_structure dof_structure;
    mutable bool dof_enumeration_made;
    mutable size_type nb_total_dof;
    size_type auto_add_elt_K; /* Degree of the fem for automatic addition */
                       /* of element option. (-1 = no automatic addition) */
    dim_type Qdim; /* this is the "global" target_dim */
    dim_type QdimM, QdimN; /* for matrix field with QdimM lines and QdimN */
                           /* columnsQdimM * QdimN = Qdim.                */
    std::vector<size_type> dof_partition;
    
  public :
    typedef base_node point_type;
    typedef tab_scal_to_vect<mesh::ind_cv_ct> ind_dof_ct;
    typedef tab_scal_to_vect<mesh::ind_pt_face_ct> ind_dof_face_ct;

    void update_from_context(void) const {}

    /** Get the set of convexes where a finite element has been assigned.
     */
    inline const dal::bit_vector &convex_index(void) const
      { return fe_convex; }
    
    /// Return a reference to the underlying mesh.
    const mesh &linked_mesh(void) const { return *linked_mesh_; }

    /** Set the degree of the fem for automatic addition
     *  of element option. K=-1 disables the automatic addition.
     */
    void set_auto_add(size_type K) {
      auto_add_elt_K = K;
    }

    /** Return the Q dimension. A mesh_fem used for scalar fields has
	Q=1, for vector fields, Q is typically equal to
	linked_mesh().dim().
     */
    dim_type get_qdim() const { return Qdim; }
    /** Change the Q dimension */
    void set_qdim(dim_type q) {
      if (q != Qdim || q != QdimM) {
	QdimM = Qdim = q; QdimN = 1;
	dof_enumeration_made = false; touch();
      }
    }

    /** Set the dimension for a matrix field. The Q dimension will
	always be the product of get_qdim_m and get_qdim_n */
    void set_qdim_mn(dim_type M, dim_type N) {
      if (M != QdimM || N != QdimN) {
	QdimM = M; QdimN = N; Qdim = N*M;
	dof_enumeration_made = false; touch();
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
    /** shortcut for set_finite_element(linked_mesh().convex_index(),pf); */
    void set_finite_element(pfem pf);
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
    
    /** Return the fem associated with an element (in no fem is
	associated, the function will crash! use the convex_index() of
	the mesh_fem to check that a fem is associated to a given
	convex) */
    pfem fem_of_element(size_type cv) const { return  f_elems[cv]; }
    /** Give an array of the dof numbers a of convex.
     *  @param cv the convex number.
     *  @return a pseudo-container of the dof number.
     */
    ind_dof_ct
      ind_dof_of_element(size_type cv) const {
      if (!dof_enumeration_made) enumerate_dof();
      return ind_dof_ct(dof_structure.ind_points_of_convex(cv),
				 Qdim /fem_of_element(cv)->target_dim());
    }
    /** Give an array of the dof numbers lying of a convex face (all
      	degrees of freedom whose associated base function is non-zero
	on the convex face).
	@param cv the convex number.
	@param f the face number.
	@return a pseudo-container of the dof number.
    */
    
    ind_dof_face_ct
    ind_dof_of_face_of_element(size_type cv, short_type f) const {
      if (!dof_enumeration_made) enumerate_dof();
      return ind_dof_face_ct
	(dof_structure.ind_points_of_face_of_convex(cv, f),
	 Qdim /fem_of_element(cv)->target_dim());
    }
    /** Return the number of dof lying on the given convex face.
	@param cv the convex number.
	@param f the face number.
    */
    size_type nb_dof_of_face_of_element(size_type cv, short_type f) const {
      pfem pf = f_elems[cv];
      return dof_structure.structure_of_convex(cv)->nb_points_of_face(f)
	* Qdim / pf->target_dim();
    }

    /** Return the number of  degrees of freedom attached to a given convex.
	@param cv the convex number.
    */
    size_type nb_dof_of_element(size_type cv) const
    { pfem pf = f_elems[cv]; return pf->nb_dof(cv) * Qdim / pf->target_dim(); }
    
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
    base_node point_of_dof(size_type cv, size_type i) const;
    /** Return the geometrical location of a degree of freedom.
	@param d the global dof number.
    */
    base_node point_of_dof(size_type d) const;
    /** Return the dof component number (0<= x <Qdim) */
    dim_type dof_qdim(size_type d) const;
    /** Shortcut for convex_to_dof(d)[0] 
 	@param d the global dof number.
    */
    size_type first_convex_of_dof(size_type d) const;
    /** Return the list of convexes attached to the specified dof 
	@param d the global dof number.
	@return an array of convex numbers.
     */
    const mesh::ind_cv_ct &convex_to_dof(size_type d) const;
    /** Renumber the degrees of freedom. You should not have
     * to call this function, as it is done automatically */
    void enumerate_dof(void) const;
    /// Return the total number of degrees of freedom.
    size_type nb_dof(void) const
      { if (!dof_enumeration_made) enumerate_dof(); return nb_total_dof; }
    /** Get a list of dof lying on a given mesh_region.
	@param b the mesh_region.
	@return the list in a dal::bit_vector.
    */
    dal::bit_vector dof_on_region(const mesh_region &b) const;
    dal::bit_vector dof_on_set(const mesh_region &b) const 
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
      return (cv < dof_partition.size() ? dof_partition[cv] : 0); 
    }
    void clear_dof_partition() { dof_partition.clear(); }

    /* explicit calls to parent class 
       for HP aCC and mipspro CC who complain about hidden functions 
       (they're right)
    */
    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_DELETE &);
    void receipt(const MESH_ADD_CONVEX &m);
    void receipt(const MESH_SUP_CONVEX &m);
    void receipt(const MESH_SWAP_CONVEX &m);
    void receipt(const MESH_REFINE_CONVEX &m);
    
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
    void clear(void);
    /** Read the mesh_fem from a stream. 
	@param ist the stream.
     */
    void read_from_file(std::istream &ist);
    /** Read the mesh_fem from a file.
        @param name the file name. */
    void read_from_file(const std::string &name);
    /** Write the mesh_fem to a stream. */
    void write_to_file(std::ostream &ost) const;
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
