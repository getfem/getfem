/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_mesh_fem.h : Finite element methods on convex meshes. */
/*     									   */
/*                                                                         */
/* Date : December 21, 1999.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 1999-2002  Yves Renard.                                   */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
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


#ifndef GETFEM_MESH_FEM_H__
#define GETFEM_MESH_FEM_H__

/* *********************************************************************** */
/*									   */
/* Ameliorations :                                                         */
/*    - faire un vrai Cutill-Mc Kee pour la numerotation des ddl.          */
/*                                                                         */
/* *********************************************************************** */

#include <getfem_mesh.h>
#include <getfem_fem.h>

#undef IS_DEPRECATED
#define IS_DEPRECATED

namespace getfem
{

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

  /*
    structure for iteration over the dofs when Qdim != 1 and target_dim == 1
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
    CONT c;
    dim_type N;
    
  public :

    bool empty(void) const { return c.empty(); }
    size_type size(void) const { return c.size() * N; }

    const_iterator begin(void) const { return iterator(c.begin(), N, 0); }
    const_iterator end(void) const { return iterator(c.end(), N, 0); }
    const_reverse_iterator rbegin(void) const
      { return const_reverse_iterator(end()); }
    const_reverse_iterator rend(void) const
      { return const_reverse_iterator(begin()); }
    
    value_type front(void) const { return *begin(); }
    value_type back(void) const { return *(--(end())); }

    tab_scal_to_vect(void) {}
    tab_scal_to_vect(const CONT &cc, dim_type n) : c(cc), N(n) {}
    
    value_type operator [](size_type ii) const { return *(begin() + ii);}
  };



  struct fem_dof {
    base_node P;
    pdof_description pnd;
  };

  struct intfem  { // integrable fem
    pfem pf;
    pintegration_method pi;
    bool operator < (const intfem &l) const;
    intfem(pfem ppf, pintegration_method ppi) { pf = ppf; pi = ppi; }
    intfem(void) { }
  };
  
  typedef const intfem * pintfem;
  pintfem give_intfem(pfem ppf, pintegration_method ppi);
  
  typedef tab_scal_to_vect<bgeot::ref_mesh_point_ind_ct> ref_mesh_dof_ind_ct;
  typedef tab_scal_to_vect<bgeot::ind_ref_mesh_point_ind_ct> 
    ind_ref_mesh_dof_ind_ct;
  /// Describe a finite element method linked to a mesh.
  class mesh_fem : public getfem_mesh_receiver, public context_dependencies {
  protected :
    
    dal::dynamic_array<pintfem> f_elems;
    dal::bit_vector fe_convex;
    getfem_mesh *linked_mesh_;
    mutable bgeot::mesh_structure dof_structure;
    mutable bool dof_enumeration_made;
    mutable size_type nb_total_dof;
    bool is_valid;
    dim_type Qdim; /* this is the "global" target_dim */
    
  public :
    
    typedef base_node point_type;
    void update_from_context(void) const {}

    /** Gives in a structure dal::bit\_vector all convexes of the
     *          mesh where a finite element is defined.
     */
    inline const dal::bit_vector &convex_index(void) const
      { return fe_convex; }
    
    /// Gives a reference to the linked mesh of type getfem\_mesh.
    getfem_mesh &linked_mesh(void) const { return *linked_mesh_; }

    dim_type get_qdim() const { return Qdim; }
    void     set_qdim(dim_type q) { if (q != Qdim) 
      { Qdim = q; dof_enumeration_made = false; touch(); }}

    /** Set on the convex of index i the integrable finite element method
     *          with the description pif which is of type pintfem.
     */
    void set_finite_element(size_type cv, pintfem pif);
    /** Set on the convex of index i the finite element method
     *          with the description pf which is of type pfem and ppi of
     *          type pintegration_method.
     */
    void set_finite_element(size_type cv, pfem ppf,
			    pintegration_method ppi=0)
      { set_finite_element(cv, give_intfem(ppf, ppi)); }	
    /** Set on all the convexes of indexes in bv, which is of type
     *          dal::bit\_vector, the finite element method
     *          with the description pf which is of type pfem and ppi of
     *          type pintegration_method. 
     *  The argument ppi is optional. If omitted, the dummy integration
     * method IM_NONE() will be used.
     */
    void set_finite_element(const dal::bit_vector &cvs, pfem ppf,
			    pintegration_method ppi = 0);
    /** shortcut for set_finite_element(linked_mesh().convex_index(),pf,ppf); */
    void set_finite_element(pfem pf, pintegration_method ppi = 0);
    /** Set a classical (i.e. lagrange polynomial) finite element on
	the convexes listed in cvs (using getfem::classical_fem). If
	im_degree is not specified then IM_NONE will by used. If it is
	specified, the an appropriate approximated integration method
	will be selected (using getfem::classical_approx_im)
    */
    void set_classical_finite_element(const dal::bit_vector &cvs, 
				      dim_type fem_degree, dim_type im_degree=dim_type(-1));
    /** Similar to set_classical_finite_element, but uses discontinuous lagrange elements */
    void set_classical_discontinuous_finite_element(const dal::bit_vector &cvs, 
						    dim_type fem_degree, dim_type im_degree=dim_type(-1));
    /** shortcut for set_classical_finite_element(linked_mesh().convex_index(),...) */
    void set_classical_finite_element(dim_type fem_degree, dim_type im_degree=dim_type(-1));
    /** shortcut for set_classical_discontinuous_finite_element(linked_mesh().convex_index(),...) */
    void set_classical_discontinuous_finite_element(dim_type fem_degree, dim_type im_degree=dim_type(-1));
    
    /** return the fem associated with an element (in no fem is
	associated, the function will crash! use the convex_index() of
	the mesh_fem to check that a fem is associated to a given
	convex) */
    pfem fem_of_element(size_type cv) const
      { return  f_elems[cv]->pf; }
    pintegration_method int_method_of_element(size_type cv) const
      { return  f_elems[cv]->pi; }
    /** Gives an array of the degrees of freedom of the element
     *           of the convex of index i. 
     */
    ref_mesh_dof_ind_ct
      ind_dof_of_element(size_type ic) const {
      if (!dof_enumeration_made) enumerate_dof();
      return ref_mesh_dof_ind_ct(dof_structure.ind_points_of_convex(ic),
				 Qdim /fem_of_element(ic)->target_dim());
    }
    ind_ref_mesh_dof_ind_ct
    ind_dof_of_face_of_element(size_type cv, short_type f) const {
      if (!dof_enumeration_made) enumerate_dof();
      return ind_ref_mesh_dof_ind_ct
	(dof_structure.ind_points_of_face_of_convex(cv, f),
	 Qdim /fem_of_element(cv)->target_dim());
    }
    size_type nb_dof_of_face_of_element(size_type cv, short_type f) const {
      pfem pf = f_elems[cv]->pf;
      return dof_structure.structure_of_convex(cv)->nb_points_of_face(f)
	* Qdim / pf->target_dim();
    }

    /** Gives the number of  degrees of freedom of the element
     *           of the convex of index i. 
     */
    size_type nb_dof_of_element(size_type cv) const {
      pfem pf = f_elems[cv]->pf;
      return pf->nb_dof(cv) * Qdim / pf->target_dim();
    }
    /** Gives the point (base_node)  corresponding to the 
     *          degree of freedom i  of the element of index cv.
     */
    const base_node &reference_point_of_dof(size_type cv,size_type i) const {
      pfem pf = f_elems[cv]->pf;
      return pf->node_of_dof(cv, i * pf->target_dim() / Qdim);
    }
    /** Gives the point (base_node) corresponding to the degree of freedom
     *  i of the element of index cv in the element of reference.
     */
    base_node point_of_dof(size_type cv, size_type i) const;
    /** Gives the point (base_node)  corresponding to the 
     *          degree of freedom with global index i.
     */
    base_node point_of_dof(size_type d) const;
    /* Gives the dof component number (0<= x <Qdim) */
    dim_type dof_qdim(size_type d) const;
    /** Shortcut for convex_to_dof(d)[0] */
    size_type first_convex_of_dof(size_type d) const;
    size_type ind_in_first_convex_of_dof(size_type d) const;
    /** Return the list of convexes attached to the specified dof */
    bgeot::mesh_convex_ind_ct convex_to_dof(size_type ip) const;
    /** Renumbers the degrees of freedom. You should not have
     * to call this function */
    void enumerate_dof(void) const;
    /// Gives the total number of degrees of freedom.
    size_type nb_dof(void) const
      { if (!dof_enumeration_made) enumerate_dof(); return nb_total_dof; }
    dal::bit_vector dof_on_set(size_type b) const;
    dal::bit_vector dof_on_boundary(size_type b) const IS_DEPRECATED
    { return  dof_on_set(b); }
    void clear(void);
    
    /// Add to the boundary b the face f of the element i.
    void add_boundary_elt(size_type b, size_type c, short_type f) IS_DEPRECATED
    { linked_mesh().add_face_to_set(b, c, f); }
    /// Says whether or not element i is on the boundary b. 
    bool is_convex_on_boundary(size_type c, size_type b) const IS_DEPRECATED
    { return linked_mesh().is_convex_in_set(b, c); }
    bool is_face_on_boundary(size_type b, size_type c, short_type f)
      const IS_DEPRECATED { return linked_mesh().is_face_in_set(b,c,f); }
    /** returns the list of convexes on the boundary b */
    const dal::bit_vector &convex_on_boundary(size_type b) const IS_DEPRECATED
    { return linked_mesh().convexes_in_set(b); }
    const mesh_cvf_set::face_bitset
      &faces_of_convex_on_boundary(size_type c, size_type b) const 
      IS_DEPRECATED { return linked_mesh().faces_of_convex_in_set(c,b); }
    /** returns the list of boundary numbers */
    const dal::bit_vector &get_valid_boundaries() const IS_DEPRECATED
    { return linked_mesh().get_valid_sets(); }
    
    void sup_boundaries_of_convex(size_type c) IS_DEPRECATED 
    { linked_mesh().sup_convex_from_sets(c); }
    void sup_boundary_elt(size_type b, size_type c, short_type f)
      IS_DEPRECATED { linked_mesh().sup_face_from_set(b,c,f); }
    void sup_boundary(size_type b) IS_DEPRECATED
    { linked_mesh().sup_set(b); }
    void swap_boundaries_convex(size_type c1, size_type c2) IS_DEPRECATED
    { linked_mesh().swap_convex_in_sets(c1, c2); }

    /* explicit calls to parent class 
       for HP aCC and mipspro CC who complain about hidden functions 
       (they're right)
    */
    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_DELETE &);
    void receipt(const MESH_ADD_CONVEX &m) {getfem_mesh_receiver::receipt(m); }
    void receipt(const MESH_SUP_CONVEX &m);
    void receipt(const MESH_SWAP_CONVEX &m);
    
    size_type memsize() const {
      return dof_structure.memsize() + 
	sizeof(mesh_fem) - sizeof(bgeot::mesh_structure) +
	f_elems.memsize() + fe_convex.memsize();
    }
    
    mesh_fem(getfem_mesh &me, dim_type Q = 1);
    virtual ~mesh_fem();
    void read_from_file(std::istream &ist);
    void read_from_file(const std::string &name);
    void write_to_file(std::ostream &ost) const;
    void write_to_file(const std::string &name, bool with_mesh=false) const;
  };
  
}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MESH_FEM_H__  */
