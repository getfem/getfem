/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2000-2015 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/

/**@file getfem_Navier_Stokes.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>, Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date April 15, 2005.
   @brief Navier-Stokes dynamic brick.
*/
#ifndef GETFEM_NAVIER_STOKES_H__
#define GETFEM_NAVIER_STOKES_H__


#include "getfem_assembling_tensors.h"

namespace getfem {

  /**
     assembly of Tangent matrix for Navier-Stokes.
     @ingroup asm
  */
  template<typename MAT, typename VECT>
  void asm_navier_stokes_tgm(const MAT &M, 
			     const mesh_im &mim, 
			     const mesh_fem &mf,
			     const VECT &U,
			     const mesh_region &rg = mesh_region::all_convexes()) {
    GMM_ASSERT1(mf.get_qdim() == mf.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

    getfem::generic_assembly assem;
    assem.set("u=data(#1);"
	      "t=comp(vGrad(#1).vBase(#1).vBase(#1));"
	      "M(#1, #1) += u(i).t(i,k,j,:,k,:,j);"
	      "M(#1, #1) += u(i).t(:,j,k,:,k,i,j);"
	      "M(#1, #1) += u(i).t(i,j,j,:,k,:,k)/2;"
	      "M(#1, #1) += u(i).t(:,k,k,:,j,i,j)/2;");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mat(const_cast<MAT&>(M));
    assem.push_data(U);
    assem.assembly(rg);
  }

  /**
     assembly of right hand side for Navier-Stokes.
     @ingroup asm
  */
  template<typename VECT1, typename VECT2> 
  void asm_navier_stokes_rhs(const VECT1 &V, 
			     const mesh_im &mim, 
			     const mesh_fem &mf,
			     const VECT2 &U,
			     const mesh_region &rg = mesh_region::all_convexes()) {
    GMM_ASSERT1(mf.get_qdim() == mf.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

    getfem::generic_assembly assem;
    assem.set("u=data(#1);"
	      "t=comp(vBase(#1).vGrad(#1).vBase(#1));"
	      "V(#1) += u(i).u(j).t(i,k,j,k,l,:,l);"
	      "V(#1) += u(i).u(j).t(i,k,j,l,l,:,k)/2;");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_vec(const_cast<VECT1&>(V));
    assem.push_data(U);
    assem.assembly(rg);
  }

 

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_NAVIER_STOKES_H__ */
