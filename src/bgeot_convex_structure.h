/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_convex_structure.h :  structure of a convex.           */
/*     									   */
/*                                                                         */
/* Date : December 20, 1999.                                               */
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


#ifndef __BGEOT_CONVEX_STRUCTURE_H
#define __BGEOT_CONVEX_STRUCTURE_H

#include <dal_ref.h>
#include <dal_alloc.h>
#include <ftool.h>
#include <bgeot_config.h>
#include <bgeot_tensor.h>
#include <bgeot_poly.h>
#include <dal_fonc_tables.h>

namespace bgeot
{
  class convex_structure;

  /// Pointer on a convex structure description. 
  typedef const convex_structure * pconvex_structure;
  ///
  typedef std::vector<pconvex_structure>        convex_structure_faces_ct;
  typedef std::vector<short_type>               convex_ind_ct;
  typedef dal::tab_ref_index_ref< convex_ind_ct::const_iterator,
	         convex_ind_ct::const_iterator> ref_convex_ind_ct;

  /**  Structure of a convex.  This class 
   *       is not to be manipulate by itself. Use pconvex\_structure and
   *       the functions written to produce the convex structures from
   *       classicals convexes (simplexes, polygonals ...). The reason is
   *       that there is no need for having more than one convex structure
   *       for the same type of convex. 
   */
  class convex_structure
  {
    protected :

      dim_type Nc;
      short_type nbpt, nbf;
      convex_structure_faces_ct  faces_struct;
      std::vector<convex_ind_ct> faces;
      convex_ind_ct              _dir_points;
      pconvex_structure basic_pcvs;

    public :

      /// Number of faces.
      inline short_type nb_faces(void)  const { return nbf;  }
      /// Dimension of the convex.
      inline dim_type  dim(void)        const { return Nc;   }
      /// Number of vertexes.
      inline short_type nb_points(void) const { return nbpt; }
      /// Original structure (if concerned).
      pconvex_structure basic_structure(void) const 
      { return basic_pcvs; }
      /// Number of vertexes of face i.
      inline short_type nb_points_of_face(short_type i) const
      { return faces[i].size(); }
      /** Gives an array of the indexes of the vertexes of face i into
       *  the array of the vertexes of the convex. ind\_points\_of\_face(i)[j]
       *  is the vertex j of the face i.
       */
      inline const convex_ind_ct &ind_points_of_face(short_type i) const
      { return faces[i]; }

      inline const convex_ind_ct &ind_dir_points() const
      { return _dir_points; }
      /** Gives a pointer array on the structures of the faces.
       *   faces\_structure()[i] is a pointer on the structure of the face i.
       */
      inline const convex_structure_faces_ct &faces_structure(void) const
      { return faces_struct; }

      inline ref_convex_ind_ct ind_dir_points_of_face(short_type i) const
      {
	return ref_convex_ind_ct(faces[i].begin(),
				 faces_struct[i]->ind_dir_points().begin(),
				 faces_struct[i]->ind_dir_points().end());
      }

      void init_for_adaptative(pconvex_structure cvs);
      void add_point_adaptative(short_type i, short_type f);

  };

  /** @name functions on convex structures
   */
  //@{

  /** Print the details of the convex structure cvs to the output stream o.
   *   For debuging purpose.
   */
  std::ostream &operator << (std::ostream &o,
				   const convex_structure &cv);

  /// Gives a pointer on the structures of a simplex of dimension d.
  pconvex_structure simplex_structure(dim_type d);
  /// Gives a pointer on the structures of a parallelepiped of dimension d.
  pconvex_structure parallelepiped_structure(dim_type d);
  /// Gives a pointer on the structures of a polygon with n vertex.
  pconvex_structure polygon_structure(short_type);
  /** Gives a pointer on the structures of a convex which is the direct
   *   product of the convexes represented by *pcvs1 and *pcvs2.
   */
  pconvex_structure convex_product_structure(pconvex_structure,
					     pconvex_structure);
  /** Gives a pointer on the structures of a prism of dimension d.
   *   i.e. the direct product of a simplex of dimension d-1 and a segment.
   */
  inline pconvex_structure prism_structure(dim_type nc)
  { 
    return convex_product_structure(simplex_structure(nc-1),
				    simplex_structure(1));
  }

  /// Simplex structure with the Lagrange grid of degree k.
  pconvex_structure simplex_structure(dim_type n, short_type k);

  //@}

}  /* end of namespace bgeot.                                             */


#endif /* __BGEOT_CONVEX_STRUCTURE_H                                      */
