/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot) version 1.0                    */
/* File    :  bgeot_convex_ref.h :  convexes of reference                  */
/*     									   */
/* Date : Septembre 28, 2001.                                              */
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


#ifndef __BGEOT_CONVEX_REF_H
#define __BGEOT_CONVEX_REF_H

#include <bgeot_convex.h>

namespace bgeot
{

  /* structures de reference.                                             */

  class convex_of_reference : public convex<base_node>
  {
    protected : 

      std::vector<base_vector> _normals;

    public :

      virtual scalar_type is_in(const base_node &) const = 0;
      virtual scalar_type is_in_face(short_type, const base_node &) const =0;
      const std::vector<base_vector> &normals(void) const { return _normals; }
  };

  typedef const convex_of_reference * pconvex_ref;

  pconvex_ref simplex_of_reference(dim_type nc, short_type k);
  inline pconvex_ref simplex_of_reference(dim_type nc)
  { return  simplex_of_reference(nc, 1); }
  pconvex_ref parallelepiped_of_reference(dim_type nc);
  pconvex_ref convex_ref_product(pconvex_ref a, pconvex_ref b);

  /* fonctions en sursis ... */
  // pconvex_ref multiply_convex_of_reference(pconvex_ref a, dim_type n);

  // pconvex_ref nonconforming_triangle_ref(void);

}  /* end of namespace bgeot.                                             */


#endif /* __BGEOT_CONVEX_REF_H                                            */
