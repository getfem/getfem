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


#ifndef __GETFEM_MESH_FEM_H
#define __GETFEM_MESH_FEM_H

/* *********************************************************************** */
/*									   */
/* Ameliorations :                                                         */
/*    - faire un vrai Cutill-Mc Kee pour la numerotation des ddl.          */
/*    - gerer le rafinement / derafinement et le suivi de bord.            */
/*      avec possibilite de bord courbes ... ?                             */
/*                                                                         */
/* *********************************************************************** */

#include <getfem_mesh.h>
#include <getfem_fem.h>
#include <getfem_mat_elem.h>

namespace getfem
{

  struct fem_dof {
    base_node P;
    pdof_description pnd;
  };
  
  ///  Describe a boundary as a list of faces of elements.
  struct boundary_description  {
    
    dal::bit_vector cvindex;
    dal::dynamic_tree_sorted<size_type> cv_in;
    dal::dynamic_array<dal::bit_vector> faces;
    
    /** Add a boudary element from the face f of the convex of index
     *          i of the mesh.
     */
    void add_elt(size_type c, short_type f);
    /** Delete a boudary element which is the face f of the convex of 
     *          index i of the mesh.
     */
    void sup_elt(size_type c, short_type f);
    /** Delete all boudary element linked with the convex of 
     *          index i of the mesh.
     */
    void sup_convex(size_type c);
    /** Gives in a structure dal::bit\_vector all the faces of convex
     *          of index i in the boundary.
     */
    const dal::bit_vector &faces_of_convex(size_type c) const;
    void swap_convex(size_type c1, size_type c2);
    size_type nb_convex(void) { return cvindex.card(); }
    void clear(void) { cv_in.clear(); faces.clear(); cvindex.clear(); }
    
  };

  struct intfem  { // integrable fem
    pfem pf;
    pintegration_method pi;
    bool operator < (const intfem &l) const;
    intfem(pfem ppf, pintegration_method ppi) { pf = ppf; pi = ppi; }
    intfem(void) { }
  };
  
  typedef const intfem * pintfem;
  pintfem give_intfem(pfem ppf, const pintegration_method ppi);
  
  typedef bgeot::ref_mesh_point_ind_ct ref_mesh_dof_ind_ct;
  
  /// Describe a finite element method linked to a mesh.
  class mesh_fem : public getfem_mesh_receiver {
  protected :
    
    dal::dynamic_array<boundary_description> boundaries;
    dal::bit_vector valid_boundaries;
    dal::dynamic_array<pintfem> f_elems;
    dal::bit_vector fe_convex;
    getfem_mesh *_linked_mesh;
    mutable bgeot::mesh_structure dof_structure;
    mutable bool dof_enumeration_made;
    mutable size_type nb_total_dof;
    bool is_valid;
    dim_type Qdim;
    
  public :
    
    typedef base_node point_type;
    
    /** Gives in a structure dal::bit\_vector all convexes of the
     *          mesh where a finite element is defined.
     */
    inline const dal::bit_vector &convex_index(void) const
      { return fe_convex; }
    
    /// Gives a pointer to the linked mesh of type getfem\_mesh.
    getfem_mesh &linked_mesh(void) const { return *_linked_mesh; }
    /** Set on the convex of index i the integrable finite element method
     *          with the description pif which is of type pintfem.
     */
    void set_finite_element(size_type cv, pintfem pif);
    /** Set on the convex of index i the finite element method
     *          with the description pf which is of type pfem and ppi of
     *          type pintegration_method.
     */
    void set_finite_element(size_type cv, pfem ppf,
			    const pintegration_method ppi)
      { set_finite_element(cv, give_intfem(ppf, ppi)); }	
    /** Set on all the convexes of indexes in bv, which is of type
     *          dal::bit\_vector, the finite element method
     *          with the description pf which is of type pfem and ppi of
     *          type pintegration_method.
     */
    void set_finite_element(const dal::bit_vector &cvs, pfem ppf,
			    const pintegration_method ppi);
    pfem fem_of_element(size_type cv) const
      { return  f_elems[cv]->pf; }
    const pintegration_method &int_method_of_element(size_type cv) const
      { return  f_elems[cv]->pi; }
    /** Gives an array of the degrees of freedom of the element
     *           of the convex of index i. 
     */
    ref_mesh_dof_ind_ct ind_dof_of_element(size_type ic) const {
      if (!dof_enumeration_made) enumerate_dof();
      return dof_structure.ind_points_of_convex(ic);
    }
    bgeot::ind_ref_mesh_point_ind_ct 
    ind_dof_of_face_of_element(size_type cv, short_type f)
      { return dof_structure.ind_points_of_face_of_convex(cv, f); }
    /** Gives the number of  degrees of freedom of the element
     *           of the convex of index i. 
     */
    size_type nb_dof_of_element(size_type cv) const {
      return f_elems[cv]->pf->nb_dof();
    }
    /** Gives the point (base_node)  corresponding to the 
     *          degree of freedom i  of the element of index cv.
     */
    const base_node &reference_point_of_dof(size_type cv,size_type i) const {
      return f_elems[cv]->pf->node_of_dof(i);
    }
    /** Gives the point (base_node) corresponding to the degree of freedom
     *  i of the element of index cv in the element of reference.
     */
    base_node point_of_dof(size_type cv, size_type i) const;
    /** Gives the point (base_node)  corresponding to the 
     *          degree of freedom with global index i.
     */
    base_node point_of_dof(size_type d) const;
    size_type first_convex_of_dof(size_type d) const
      { return dof_structure.first_convex_of_point(d); }
    size_type ind_in_first_convex_of_dof(size_type d) const {
      if (!dof_enumeration_made) enumerate_dof(); 
      return dof_structure.ind_in_first_convex_of_point(d);
    }
    void enumerate_dof(void) const;
    /// Gives the total number of degrees of freedom.
    size_type nb_dof(void) const
      { if (!dof_enumeration_made) enumerate_dof(); return nb_total_dof; }
    dal::bit_vector dof_on_boundary(size_type b) const;
    void clear(void);
    
    /// Add to the boundary b the face f of the element i.
    void add_boundary_elt(size_type b, size_type c, short_type f)
      { valid_boundaries.add(b); boundaries[b].add_elt(c, f); }
    /// Says whether or not element i is on the boundary b. 
    bool is_convex_on_boundary(size_type c, size_type b) const
      { return (valid_boundaries[b] && boundaries[b].cvindex[c]); }
    const dal::bit_vector &convex_on_boundary(size_type b) const;
    const dal::bit_vector &faces_of_convex_on_boundary(size_type c,
						       size_type b) const;
    /* returns the list of boundary numbers  [JP] */
    const dal::bit_vector &get_valid_boundaries() const 
      { return valid_boundaries; }
    
    void sup_boundaries_of_convex(size_type c);
    void sup_boundary_elt(size_type b, size_type c, short_type f)
      { if (valid_boundaries[b]) boundaries[b].sup_elt(c,f); }
    void sup_boundary(size_type b)
      { valid_boundaries.sup(b); boundaries[b].clear(); }
    void swap_boundaries_convex(size_type c1, size_type c2);
    
    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_DELETE &);
    void receipt(const MESH_SUP_CONVEX &m);
    void receipt(const MESH_SWAP_CONVEX &m);
    void receipt(const MESH_REFINE_CONVEX &m);
    void receipt(const MESH_UNREFINE_CONVEX &m);
    void receipt(const MESH_FEM_TOUCH &m);
    
    size_type memsize() const {
      return dof_structure.memsize() + 
	sizeof(mesh_fem) - sizeof(bgeot::mesh_structure) +
	boundaries.memsize() + valid_boundaries.memsize() +
	f_elems.memsize() + fe_convex.memsize();
    }
    
    mesh_fem(getfem_mesh &me, dim_type Q = 1); 
    virtual ~mesh_fem();
    void read_from_file(std::istream &ist);
    void read_from_file(const std::string &name);
    void write_to_file(std::ostream &ost) const;
    void write_to_file(const std::string &name) const;
  };
  
}  /* end of namespace getfem.                                             */


#endif /* __GETFEM_MESH_FEM_H  */
