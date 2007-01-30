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

/**@file getfem_regular_meshes.h
   @author  Yves Renard <Yves.Renard@insa-toulouse.fr>
   @date December 20, 1999.
   @brief Build regular meshes.
*/
#ifndef GETFEM_REGULAR_MESHES_H__
#define GETFEM_REGULAR_MESHES_H__

#include "getfem_mesh.h"

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

  void parallelepiped_regular_simplex_mesh_(mesh &me, dim_type N,
    const base_node &org, const base_small_vector *ivect, const size_type *iref);

  template<class ITER1, class ITER2>
    void parallelepiped_regular_simplex_mesh(mesh &me,
					     dim_type N,
	     const base_node &org, ITER1 ivect, ITER2 iref)
  { 
    std::vector<base_small_vector> vect(N);
    std::copy(ivect, ivect+N, vect.begin());
    std::vector<size_type> ref(N);
    std::copy(iref, iref+N, ref.begin());
    parallelepiped_regular_simplex_mesh_(me, N, org, &(vect[0]),
					 &(ref[0]));
  } 

  void parallelepiped_regular_prism_mesh_(mesh &me, dim_type N,
      const base_node &org, const base_small_vector *ivect, const size_type *iref);

  template<class ITER1, class ITER2>
    void parallelepiped_regular_prism_mesh(mesh &me,
					     dim_type N,
	     const base_node &org, ITER1 ivect, ITER2 iref)
  { 
    std::vector<base_small_vector> vect(N);
    std::copy(ivect, ivect+N, vect.begin());
    std::vector<size_type> ref(N);
    std::copy(iref, iref+N, ref.begin());
    parallelepiped_regular_prism_mesh_(me, N, org, &(vect[0]),
					 &(ref[0]));
  } 

  void parallelepiped_regular_mesh_(mesh &me, dim_type N,
    const base_node &org, const base_small_vector *ivect, const size_type *iref, bool linear_gt);

  template<class ITER1, class ITER2>
    void parallelepiped_regular_mesh(mesh &me,
					     dim_type N,
                                     const base_node &org, ITER1 ivect, ITER2 iref, bool linear_gt=false)
  { 
    std::vector<base_small_vector> vect(N);
    std::copy(ivect, ivect+N, vect.begin());
    std::vector<size_type> ref(N);
    std::copy(iref, iref+N, ref.begin());
    parallelepiped_regular_mesh_(me, N, org, &(vect[0]), &(ref[0]), linear_gt);
  } 

  /**
     Build a regular mesh of the unit square/cube/, etc.
     @param m the output mesh.

     @param pgt the geometric transformation to use. For example, use 
     @code
     pgt = geometric_trans_descriptor("GT_PK(2,1"); // to build a mesh of triangles
     pgt = geometric_trans_descriptor("QK(3,2)"); // to build a mesh of order 2 parallelepipeded
     @endcode

     @param nsubdiv is the number of cells in each direction.  

     @param noised if set, will cause the interior nodes to be randomly "shaken".
   */
  void regular_unit_mesh(mesh& m, std::vector<size_type> nsubdiv, 
			 bgeot::pgeometric_trans pgt, bool noised = false);

  /**
     Build a regular mesh parametrized by the string st.
     The format of st is the following
     std::string st("GT='GT_PK(2,1)'; NSUBDIV=[5,5]; ORG=[0,0]; SIZES=[1,1]; NOISED=0");
     where GT is the geometric transformation, NSUBDIV a vector of the number
     of subdivisions in each coordinate (default value 2), ORG is the origin
     of the mesh (default value [0,0,...]), SIZES is a vector of the sizes
     in each direction (default value [1, 1, ...] and if NOISED=1 the nodes
     of the interior of the mesh are randomly "shaken"(default value NOISED=0).
     All the parameters are optional but GT.

     @param m the output mesh.    
  */
  void regular_mesh(mesh& m, const std::string &st);

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_REGULAR_MESHES_H__  */
