/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_mesh.h : meshes for computations.                     */
/*     									   */
/*                                                                         */
/* Date : November 05, 1999.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; version 2 of the License.                 */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software Foundation, */
/* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.         */
/*                                                                         */
/* *********************************************************************** */


#ifndef __GETFEM_MESH_H
#define __GETFEM_MESH_H

#include <bgeot_mesh.h>
#include <bgeot_geometric_trans.h>
#include <linkmsg.h>
#include <getfem_config.h>

namespace getfem
{

  /* ********************************************************************* */
  /*								   	   */
  /*	I. Classes of message                                		   */
  /*									   */
  /* ********************************************************************* */

  struct MESH_CLEAR  /* clear message for the structure.                   */
  { operator int(void) const { return 0; } };
  struct MESH_ADD_POINT /* point suppression message.                      */
  { 
    size_t ipt;
    operator int(void) const { return 1; }
    MESH_ADD_POINT(size_t i) { ipt = i; }
    MESH_ADD_POINT(void) {}
  };
  struct MESH_SUP_POINT
  { 
    size_t ipt;
    operator int(void) const { return 2; }
    MESH_SUP_POINT(size_t i) { ipt = i; }
    MESH_SUP_POINT(void) {}
  };
  struct MESH_SWAP_POINT
  { 
    size_t ipt1, ipt2;
    operator int(void) const { return 3; }
    MESH_SWAP_POINT(size_t i, size_t j) { ipt1 = i; ipt2 = j; }
    MESH_SWAP_POINT(void) {}
  };
  struct MESH_ADD_CONVEX
  { 
    size_t icv;
    operator int(void) const { return 4; }
    MESH_ADD_CONVEX(size_t i) { icv = i; }
    MESH_ADD_CONVEX(void) {}
  };
  struct MESH_SUP_CONVEX
  { 
    size_t icv;
    operator int(void) const { return 5; }
    MESH_SUP_CONVEX(size_t i) { icv = i; }
    MESH_SUP_CONVEX(void) {}
  };
  struct MESH_SWAP_CONVEX
  { 
    size_t icv1, icv2;
    operator int(void) const { return 6; }
    MESH_SWAP_CONVEX(size_t i, size_t j) { icv1 = i; icv2 = j; }
    MESH_SWAP_CONVEX(void) {}
  };
  struct MESH_REFINE_CONVEX
  { 
    size_t icv, nb;
    size_t *alist;
    int mtype;
    operator int(void) const { return 7; }
    MESH_REFINE_CONVEX(size_t i, size_t n, size_t *l, int m)
    { icv = i; nb = n; alist = l; mtype = m; }
    MESH_REFINE_CONVEX(void) {}
  };
  struct MESH_UNREFINE_CONVEX
  { 
    size_t icv, nb;
    size_t *alist;
    int mtype;
    operator int(void) const { return 8; }
    MESH_UNREFINE_CONVEX(size_t i, size_t n, size_t *l, int m)
    { icv = i; nb = n; alist = l; mtype = m; }
    MESH_UNREFINE_CONVEX(void) {}
  };
  struct MESH_WRITE_TO_FILE
  { 
    std::ostream *ost;
    operator int(void) const { return 9; }
    MESH_WRITE_TO_FILE(std::ostream &o) { ost = &o; }
    MESH_WRITE_TO_FILE(void) {}
  };
  struct MESH_READ_FROM_FILE
  { 
    std::istream *ist;
    operator int(void)  const { return 10; }
    MESH_READ_FROM_FILE(std::istream &i) { ist = &i; }
    MESH_READ_FROM_FILE(void) {}
  };
  struct MESH_FEM_CHANGE
  { 
    void *ptr;
    operator int(void)  const { return 11; }
    MESH_FEM_CHANGE(void *p) : ptr(p) {}
    MESH_FEM_CHANGE(void) {}
  };
  struct MESH_FEM_DELETE
  { 
    void *ptr;
    operator int(void)  const { return 12; }
    MESH_FEM_DELETE(void *p) : ptr(p) {}
    MESH_FEM_DELETE(void) {}
  };
  struct MESH_FEM_TOUCH
  { 
    void *ptr;
    operator int(void)  const { return 13; }
    MESH_FEM_TOUCH(void *p) : ptr(p) {}
    MESH_FEM_TOUCH(void) {}
  };


  class getfem_mesh_receiver : public lmsg::virtual_linkmsg_receiver
  {
    public :

      virtual void receipt(const MESH_CLEAR           &)
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_ADD_POINT       &) 
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_SUP_POINT       &) 
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_SWAP_POINT      &) 
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_ADD_CONVEX      &) 
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_SUP_CONVEX      &) 
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_SWAP_CONVEX     &)
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_REFINE_CONVEX   &) 
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_UNREFINE_CONVEX &) 
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_WRITE_TO_FILE   &) 
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_READ_FROM_FILE  &) 
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_FEM_CHANGE      &) 
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_FEM_DELETE      &) 
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_FEM_TOUCH       &) 
      { DAL_THROW(internal_error, "internal error");}
  };

  /* refinement methods are :                                              */
  /* mtype = 0 : simplexification.                                         */
  /* mtype = 1 : bank total.                                               */
  /* mtype = 2 : bank green.                                               */


  /* ********************************************************************* */
  /*								   	   */
  /*	II. Class getfem_mesh                                 		   */
  /*									   */
  /* ********************************************************************* */

  /** Describe a mesh for the computation of pde problems. This
   *      mesh is able to be link with classes which define computation
   *      methods. 
   */
  class getfem_mesh : public bgeot::mesh<base_node>
  {
    public :

      typedef lmsg::linkmsg_sender<getfem_mesh_receiver> msg_sender;

    protected :

      double eps_p;  /* infinity distance under wich two points are equal. */
      msg_sender lkmsg; /* gestionnaire de msg.                            */
      dal::dynamic_array<bgeot::pgeometric_trans> gtab;
      dal::bit_vector trans_exists;

    public :

      /// Constructor.
      getfem_mesh(dim_type NN = dim_type(-1)); 
      double eps(void) const { return eps_p; }
      const msg_sender &lmsg_sender(void) const { return lkmsg; }
      msg_sender &lmsg_sender(void) { return lkmsg; }

      /** Add the point pt to the mesh and return the index of the
       *          point. If the point is to close to an existing point, the
       *          function do not add the point and return the index of the
       *          already existing point. pt should be of type base\_node.
       */
      size_type add_point(const base_node &pt);
      /// Gives the number of points in the mesh.
      size_type nb_points(void) const { return pts.card(); }
      /// Delete the point of index i from the mesh.
      void sup_point(size_type i);
      /// Swap the indexes of points of index i and j in the whole structure.
      void swap_points(size_type i, size_type j);
      /** Search if the point pt is in (or approximatively in)
       *          the mesh, and return the index of the point, or
       *          size\_type(-1) if not found.
       */
      size_type search_point(const base_node &pt) const
      { return pts.search(pt); }

      bgeot::pgeometric_trans trans_of_convex(size_type ic) const
      { 
	if (!(trans_exists[ic]))
	  DAL_THROW(internal_error, "internal error");
	return gtab[ic]; 
      }

      /** Add a convex to the mesh. cvs is of type 
       *          bgeot::convex\_structure and "it" is an iterator on a list
       *          of indexes of points. Return the index
       *          of the convex in the mesh.
       */
      template<class ITER>
	size_type add_convex(bgeot::pgeometric_trans pgt, ITER ipts)
      { 
	bool present;
	size_type i = bgeot::mesh<base_node>::add_convex(pgt->structure(),
							 ipts, &present);
	gtab[i] = pgt; trans_exists[i] = true;
	if (!present) lmsg_sender().send(MESH_ADD_CONVEX(i));
	return i;
      }

      /** Add a convex to the mesh. cvs is of type 
       *          bgeot::convex\_structure and "it" is an iterator on a list
       *          of points of type base\_node. Return the index
       *          of the convex in the mesh.
       */
      template<class ITER>
	size_type add_convex_by_points(bgeot::pgeometric_trans pgt, ITER ipts);
	
      /** Add a simplex of dimension dim to the mesh. 
       *          "it" is an iterator on a list of indexes of the points.
       *          Return the index of the convex in the mesh.
       */
      template<class ITER>
	size_type add_simplex(dim_type dim, ITER ipts)
      { return add_convex(bgeot::simplex_trans(dim, 1), ipts); }
      /** Add a simplex of dimension dim to the mesh. 
       *          "it" is an iterator on a list of points of type base\_node.
       *          Return the index of the convex in the mesh.
       */
      template<class ITER>
	size_type add_simplex_by_points(dim_type dim, ITER ipts);
      size_type add_segment(size_type a, size_type b);
      size_type add_segment_by_points(const base_node &pt1,
				      const base_node &pt2)
      { return add_segment(add_point(pt1), add_point(pt2)); }
      size_type add_triangle(size_type a,size_type b, size_type c);
      size_type add_triangle_by_points(const base_node &p1,
				       const base_node &p2,
				       const base_node &p3);
      size_type add_tetrahedron(size_type a,
				size_type b, size_type c, size_type d);
      size_type add_tetrahedron_by_points(const base_node &p1,
					  const base_node &p2,
					  const base_node &p3,
					  const base_node &p4);
      /** Add a parallelepiped of dimension dim to the mesh. 
       *          "it" is an iterator on a list of indexes of the points.
       *          Return the index of the convex in the mesh.
       */
      template<class ITER>
	size_type add_parallelepiped(dim_type dim, const ITER &ipts);
      /** Add a parallelepiped of dimension dim to the mesh. 
       *          "it" is an iterator on a list of points of type base\_node.
       *          Return the index of the convex in the mesh.
       */
      template<class ITER>
	size_type add_parallelepiped_by_points(dim_type dim, const ITER &pts);
      /** Add a parallelepiped of dimension dim to the
       *          mesh. org is the point of type base\_node representing
       *          the origine and "it" is an iterator on a list of
       *          vectors of type base\_vector.
       *          Return the index of the convex in the mesh.
       */
      template<class ITER>
	size_type add_parallelepiped_by_vectors(dim_type dim,
				    const base_node &org, const ITER &vects);

      /** Add a prism of dimension dim to the mesh. 
       *          "it" is an iterator on a list of indexes of the points.
       *          Return the index of the convex in the mesh.
       */
      template<class ITER>
	size_type add_prism(dim_type dim, const ITER &ipts);

       /** Add a prism of dimension dim to the mesh. 
       *          "it" is an iterator on a list of points of type base\_node.
       *          Return the index of the convex in the mesh.
       */
      template<class ITER>
	size_type add_prism_by_points(dim_type dim, const ITER &pts);

      /// Delete the convex of index i from the mesh.
      void sup_convex(size_type ic);
      /** Swap the indexes of the convex of indexes i and j 
       *          in the whole structure.
       */
      void swap_convex(size_type i, size_type j);

      /* returns the normal of face 'f' evaluated at the point 'pt'       */
      /* (pt is a position in the reference convex)                       */
      base_vector normal_of_face_of_convex(size_type ic, short_type f,
					   const base_node &pt) const;
      base_vector normal_of_face_of_convex(size_type ic, short_type f,
					   size_type n) const;

      void translation(base_vector);
      void transformation(base_matrix);
  
      void optimize_structure(void);
      void clear(void);
      
      int write_to_file(const std::string &name) const;
      int write_to_file(std::ostream &ost) const;
      int read_from_file(const std::string &name);
      int read_from_file(std::istream &ist);
     
  };

  template<class ITER>
    size_type getfem_mesh::add_convex_by_points(bgeot::pgeometric_trans pgt,
					                           ITER ipts)
  {
    short_type nb = pgt->nb_points();
    static std::vector<size_type> *ind;
    static bool isinit = false;
    if (!isinit) { ind = new std::vector<size_type>(); isinit = true; }
    if (ind->size() != nb) ind->resize(nb);
    for (short_type i = 0; i < nb; ++ipts, ++i) (*ind)[i] = add_point(*ipts);
    return add_convex(pgt, ind->begin());
  }

  template<class ITER>
   size_type getfem_mesh::add_simplex_by_points(dim_type dim, ITER ipts)
  {
    return add_convex_by_points(bgeot::simplex_trans(dim, 1), ipts);
  }

  template<class ITER>
    size_type getfem_mesh::add_parallelepiped(dim_type dim, const ITER &ipts)
  { return add_convex(bgeot::parallelepiped_trans(dim, 1), ipts); }

  template<class ITER>
    size_type getfem_mesh::add_parallelepiped_by_points
    (dim_type dim, const ITER &pts)
  { return add_convex_by_points(bgeot::parallelepiped_trans(dim, 1), pts); }

  template<class ITER>
    size_type getfem_mesh::add_parallelepiped_by_vectors
    (dim_type dim, const base_node &org, const ITER &vects)
  {
    size_type nbp = (size_type(1) << size_type(dim)), i, j;
    std::vector<size_type> ipt;
    ipt.resize(nbp);
    base_node a; ITER b;

    for (i = 0; i < nbp; i++)
    {
      for (a = org, b = vects, j = 0; j < dim; ++j, ++b)
	if (i & (1 << j)) a += *b;
      ipt[i] = add_point(a);
    }

    return add_parallelepiped(dim, ipt.begin());
  }

  template<class ITER>
    size_type getfem_mesh::add_prism(dim_type dim, const ITER &ipts)
  { return add_convex(bgeot::prism_trans(dim, 1), ipts); }

  template<class ITER>
    size_type getfem_mesh::add_prism_by_points
    (dim_type dim, const ITER &pts)
  { return add_convex_by_points(bgeot::prism_trans(dim, 1), pts); }

}  /* end of namespace getfem.                                             */


#endif /* __GETFEM_MESH_H  */
