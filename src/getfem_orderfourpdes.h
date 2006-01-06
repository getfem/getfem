// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_orderfourspdes.h : assembly procedures and bricks 
//                                     for order four pdes.
// Date    : January 6, 2006.
// Authors : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2006-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================

/**@file getfem_orderfourpdes.h
   @brief assembly procedures and bricks for order four pdes.
*/
#ifndef GETFEM_ORDERFOURPDES_H_
#define GETFEM_ORDERFOURPDES_H__

#include <getfem_modeling.h>
#include <getfem_assembling_tensors.h>

namespace getfem {

  /* real version. */
  template<typename VECT1, typename VECT2, typename T>
  void asm_normal_derivative_source_term_(const VECT1 &B, const mesh_im &mim,
					  const mesh_fem &mf,
					  const mesh_fem &mf_data,
					  const VECT2 &F,
					  const mesh_region &rg, T) {
    generic_assembly assem;
    if (mf.get_qdim() == 1)
      assem.set("F=data(#2);"
		"V(#1)+=comp(Grad(#1).Normal().Base(#2))(:,i,i,j).F(j);");
    else
      assem.set("F=data(qdim(#1),#2);"
		"V(#1)+=comp(vGrad(#1).Normal.Base(#2))(:,i,k,k,j).F(i,j);");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mf_data);
    assem.push_data(F);
    assem.push_vec(const_cast<VECT1&>(B));
    assem.assembly(rg);
  }

  /* complex version. */
  template<typename VECT1, typename VECT2, typename T>
  void asm_normal_derivative_source_term_(VECT1 &B, const mesh_im &mim,
					  const mesh_fem &mf,
					  const mesh_fem &mf_data,
					  const VECT2 &F,
					  const mesh_region &rg,
					  std::complex<T>) {
    asm_normal_derivative_source_term_(gmm::real_part(B), mim, mf, mf_data,
				       gmm::real_part(F), rg, T());
    asm_normal_derivative_source_term_(gmm::imag_part(B), mim, mf, mf_data,
				       gmm::imag_part(F), rg, T());
  }

  /**
     assembly of @f$\int_\Gamma{\partial_n u f}@f$.
     @ingroup asm
  */
  template<typename VECT1, typename VECT2>
  void asm_normal_derivative_source_term(VECT1 &B, const mesh_im &mim,
					 const mesh_fem &mf,
					 const mesh_fem &mf_data,
					 const VECT2 &F,
					 const mesh_region &rg) {
    if (mf_data.get_qdim() != 1)
      DAL_THROW(invalid_argument, "invalid data mesh fem (Qdim=1 required)");
    asm_normal_derivative_source_term_(B, mim, mf, mf_data, F, rg,
		     typename gmm::linalg_traits<VECT2>::value_type());
  }

  


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_ORDERFOURPDES_H__  */
