/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_mesh.h : meshes with points.                           */
/*     									   */
/*                                                                         */
/* Date : December 20, 1999.                                               */
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

#ifndef BGEOT_MESH_H__
#define BGEOT_MESH_H__

#include <bgeot_mesh_structure.h>
#include <bgeot_convex.h>

namespace bgeot
{

  template<class PT> struct mesh_point_ct
          : public dal::dynamic_tree_sorted<PT,
                   dal::lexicographical_less<PT,
                   dal::approx_less<typename PT::value_type> > > { };

  template<class PT, class PT_TAB = mesh_point_ct<PT> > class mesh
    : public mesh_structure
  {
    public :

      typedef PT point_type;
      typedef typename PT::vector_type vector_type;

      typedef dal::tab_ref_index_ref< typename PT_TAB::const_iterator,
                     ref_mesh_point_ind_ct::const_iterator> ref_mesh_pt_ct;
      typedef dal::tab_ref_index_ref< typename PT_TAB::const_iterator, 
                     ind_ref_mesh_point_ind_ct::const_iterator> ref_mesh_face_pt_ct;
      typedef bgeot::convex<PT, ref_mesh_pt_ct> ref_convex;


    protected :

      dim_type dimension;
      PT_TAB pts;
      
    public :

      const PT_TAB &points(void) const { return pts; }
      PT_TAB &points(void) { return pts; }
      void swap_points(int i, int j)
      { pts.swap(i,j); mesh_structure::swap_points(i,j); }
     
      /* acces aux points directeurs non fait : le faire d'abord dans      */
      /* bgeot_mesh_structure avec container adapte.                       */

      /**
	  returns a (pseudo)container of the points of a given convex 
      */
      ref_mesh_pt_ct points_of_convex(size_type ic) const {
        ref_mesh_point_ind_ct rct = ind_points_of_convex(ic);
        return ref_mesh_pt_ct(pts.begin(), rct.begin(), rct.end());
      } 
    
      /** 
	  returns a (pseudo)container of points of face of a given convex 
      */
      ref_mesh_face_pt_ct points_of_face_of_convex(size_type ic, size_type f) const {
        ind_ref_mesh_point_ind_ct rct = ind_points_of_face_of_convex(ic,f);
        return ref_mesh_face_pt_ct(pts.begin(), rct.begin(), rct.end());
      }

      PT dir_point_of_convex(size_type ic, size_type j) const
      {	return pts[ind_dir_point_of_convex(ic, j)]; }

      ref_convex convex(size_type ic) const
      { return ref_convex(structure_of_convex(ic), points_of_convex(ic)); }

      dim_type dim(void) const { return dimension; }
      void clear(void)
      { dimension = dim_type(-1); mesh_structure::clear(); pts.clear(); }

    mesh(void) { dimension = dim_type(-1); }
    size_type memsize(void) const { return mesh_structure::memsize() + 
				    sizeof(mesh) - sizeof(mesh_structure) + 
				    pts.memsize(); }
  };

}  /* end of namespace bgeot.                                              */


#endif /* BGEOT_MESH_H__                                                   */
