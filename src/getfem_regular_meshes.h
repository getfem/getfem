/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_regular_meshes.h : basic refinement of meshes.        */
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


#ifndef __GETFEM_REGULAR_MESHES_H
#define __GETFEM_REGULAR_MESHES_H

#include <getfem_mesh.h>
#include <bgeot_simplexify.h>

namespace getfem
{
  /* ******************************************************************** */
  /*	    Generation de maillages reguliers de simplexes.               */
  /* ******************************************************************** */

  /* ******************************************************************** */
  /* ajoute au maillage "me" un maillage regulier de dimension "N"        */
  /* a partir du point "org" defini par les parallelepipedes engendres    */
  /* par la liste des vecteurs "vect". La liste iref est une liste        */
  /* d'entiers correspondant au nombre de parallelepidedes sur chaque     */
  /* dimension. Chaque parallelepidede est simplexifie.                   */
  /* ******************************************************************** */

  void _parallelepiped_regular_simplex_mesh(getfem_mesh &me, dim_type N,
    const base_node &org, const base_vector *ivect, const size_type *iref);

  template<class ITER1, class ITER2>
    void parallelepiped_regular_simplex_mesh(getfem_mesh &me,
					     dim_type N,
	     const base_node &org, ITER1 ivect, ITER2 iref)
  { 
    std::vector<base_vector> vect(N);
    std::copy(ivect, ivect+N, vect.begin());
    std::vector<size_type> ref(N);
    std::copy(iref, iref+N, ref.begin());
    _parallelepiped_regular_simplex_mesh(me, N, org, &(vect[0]),
					 &(ref[0]));
  } 

  void _parallelepiped_regular_prism_mesh(getfem_mesh &me, dim_type N,
      const base_node &org, const base_vector *ivect, const size_type *iref);

  template<class ITER1, class ITER2>
    void parallelepiped_regular_prism_mesh(getfem_mesh &me,
					     dim_type N,
	     const base_node &org, ITER1 ivect, ITER2 iref)
  { 
    std::vector<base_vector> vect(N);
    std::copy(ivect, ivect+N, vect.begin());
    std::vector<size_type> ref(N);
    std::copy(iref, iref+N, ref.begin());
    _parallelepiped_regular_prism_mesh(me, N, org, &(vect[0]),
					 &(ref[0]));
  } 

    void _parallelepiped_regular_mesh(getfem_mesh &me, dim_type N,
    const base_node &org, const base_vector *ivect, const size_type *iref);

  template<class ITER1, class ITER2>
    void parallelepiped_regular_mesh(getfem_mesh &me,
					     dim_type N,
	     const base_node &org, ITER1 ivect, ITER2 iref)
  { 
    std::vector<base_vector> vect(N);
    std::copy(ivect, ivect+N, vect.begin());
    std::vector<size_type> ref(N);
    std::copy(iref, iref+N, ref.begin());
    _parallelepiped_regular_mesh(me, N, org, &(vect[0]),
					 &(ref[0]));
  } 

}  /* end of namespace getfem.                                             */


#endif /* __GETFEM_REGULAR_MESHES_H  */
