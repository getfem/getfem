/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2012-2012 Liang Jin Lim
 
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
/**
@file getfem_torus.h
@brief Provides mesh and mesh fem of torus.
@date May 2014
@author Liang Jin Lim
*/

#pragma once

#ifndef GETFEM_TORUS_H__
#define GETFEM_TORUS_H__

#include "getfem/getfem_mesh_fem.h"


namespace bgeot{

/**An adaptor that adapts a two dimensional geometric_trans to include radial dimension.*/
struct torus_geom_trans : public geometric_trans{

  virtual void poly_vector_val(const base_node &, bgeot::base_vector &) const;
  virtual void poly_vector_val(const base_node &, const bgeot::convex_ind_ct &,
    bgeot::base_vector &) const;
  virtual void poly_vector_grad(const base_node &, bgeot::base_matrix &) const;
  inline virtual void poly_vector_grad(const base_node &,
    const bgeot::convex_ind_ct &, bgeot::base_matrix &) const;
  inline virtual void compute_K_matrix
    (const bgeot::base_matrix &, const bgeot::base_matrix &, bgeot::base_matrix &) const;

  virtual void poly_vector_hess(const base_node &, bgeot::base_matrix &) const;

  torus_geom_trans(bgeot::pgeometric_trans poriginal_trans);

private:
  pgeometric_trans poriginal_trans_;
};

pconvex_structure torus_structure_descriptor(pconvex_structure ori_structure);

bool is_torus_structure(pconvex_structure cvs);

pgeometric_trans torus_geom_trans_descriptor(pgeometric_trans poriginal_trans);

bool is_torus_geom_trans(pgeometric_trans pgt);

}

namespace getfem
{
  /**Torus fem, the real grad base value is modified to compute radial grad of F/R.
     It stores a reference to the original fem object. By default, torus_fem is vectorial.
     There is an option to change it to a scalar form through set_to_scalar(bool is_scalar).
     torus_mesh_fem will automatically check qdim of itself and set the form accordingly.
  */
  class torus_fem : public virtual_fem{

  public :
    virtual size_type index_of_global_dof(size_type cv, size_type i) const;
    void base_value(const base_node &, base_tensor &) const;
    void grad_base_value(const base_node &, base_tensor &) const;
    void hess_base_value(const base_node &, base_tensor &) const;
    void real_base_value(const fem_interpolation_context& c,
      base_tensor &t, bool = true) const;
    void real_grad_base_value(const fem_interpolation_context& c,
      base_tensor &t, bool = true) const;
    void real_hess_base_value(const fem_interpolation_context&,
      base_tensor &, bool = true) const;

    torus_fem(pfem pf) : virtual_fem(*pf), poriginal_fem_(pf), is_scalar_(false){
      init();
    }

    void set_to_scalar(bool is_scalar);

  protected :
    void init();

  private :
    pfem poriginal_fem_;
    bool is_scalar_;
  };
  
  /**Copy an original 2D mesh to become a torus mesh with radial dimension.*/
  class torus_mesh : public mesh
  {
  private:
    bool is_adapted_;
    
  public:    
    void adapt();
    void adapt(const getfem::mesh &original_mesh);
  };

  /**Mesh fem object that adapts */
  class torus_mesh_fem : public mesh_fem{
  public:

    torus_mesh_fem(torus_mesh &mesh, bgeot::dim_type dim) : mesh_fem(mesh, dim){}    
    void enumerate_dof(void) const;

  private:
    void adapt_to_torus_();
    void del_torus_fem_();
  };
  
  /**Perform reduction of a matrix from full 3d of torus structure to axisymmetric 2d.*/
  template <typename MAT1, typename MAT2>
  void torus_matrix_reduction(const MAT1 &M1, MAT2 &M2){
    GMM_ASSERT1(gmm::mat_ncols(M1)%3 == 0 && gmm::mat_nrows(M1)%3 == 0, 
      "Origin matrix must be 3 dimension.");

    bgeot::size_type nb_pts = gmm::mat_ncols(M1)/3;

    GMM_ASSERT1(gmm::mat_ncols(M2) == nb_pts*2 && gmm::mat_nrows(M2) == nb_pts*2, 
      "Invalid target matrix size.");

    for(int i=0; i <M2.nrows()/2; ++i){
      for(int row = 0; row < 2; ++row){
        for(int j=0; j < M2.ncols()/2; ++j){
          for(int col = 0; col < 2; ++col){
            M2(i*2 + row, j * 2 + col) = M1(i*3 + row, j * 3 + col);
            if(row == 0) M2(i*2, j*2 + col) += M1(i*3 + 2, j*3 + col);
            if(col == 0) M2(i*2 + row, j*2) += M1(i*3 + row, j*3 + 2);
            if(row == 0 && col == 0) M2(i*2, j*2) += M1(i*3 + 2, j*3 + 2);
          }
        }
      }
    }
  }

  /**Perform reduction of a vector from full 3d of torus structure to axisymmetric 2d*/
  template <typename VECT1, typename VECT2>
  void torus_vector_reduction(const VECT1 &V1, VECT2 &V2){
    GMM_ASSERT1(gmm::vect_size(V1)%3 == 0, "Origin vector must be 3 dimension.");
    bgeot::size_type nb_pts = gmm::vect_size(V1)/3;
    GMM_ASSERT1(gmm::vect_size(V2) == nb_pts*2, "Invalid target vector size.");

    for(int i=0; i <V2.size()/2; ++i){
      for(int row = 0; row < 2; ++row){
        V2[i*2 + row] = V1[i*3 + row];
        if(row == 0) V2[i*2] += V1[i*3 + 2];
      }
    }
  }
}

#endif /* GETFEM_TORUS_H__  */
