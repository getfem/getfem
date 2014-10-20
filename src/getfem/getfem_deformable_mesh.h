/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

Copyright (C) 2012-2012 Andriy Andreykiv

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

/** @file getfem_deformable_mesh.h
@author "Andriy Andreykiv" <andriy.andreykiv@gmail.com>
@date August 7, 2012.
@brief A class adaptor to deform a mesh.
*/

#pragma once
#ifndef GETFEM_DEFORMABLE_MESH_H__
#define GETFEM_DEFORMABLE_MESH_H__

#include <getfem/getfem_mesh.h>
#include <getfem/getfem_mesh_fem.h>
#include <getfem/getfem_models.h>

namespace getfem {

  /** An object function that first deforms and then remembers 
  to restore a mesh if it has to be restored 
  for other bricks. By default the mesh is deformed on 
  construct and undeformed in the destructor (by RAII principle)
  but it's also possible to specify deform_on_construct = false 
  and then call explicitely deform() and undeform() methods.
  Optional to_be_restored flag will control whether the mesh will be restored
  when the deformator destructs.
  */
  template<class VECTOR = model_real_plain_vector> 
  class temporary_mesh_deformator
  {
  public:
    temporary_mesh_deformator(const mesh& m, const mesh_fem &mf, 
      const VECTOR &dU, bool deform_on_construct = true, 
      bool to_be_restored = true) : 
      dU_(mf.nb_basic_dof()),
      mf_(mf), 
      m_(const_cast<getfem::mesh &>(mf.linked_mesh())),
      deform_on_construct_(deform_on_construct),
      is_deformed_(false),
      to_be_restored_(to_be_restored){
      mf.extend_vector(dU, dU_);
      if (deform_on_construct_) deform();
    }

    void deform(){
      if (is_deformed_) return;
      deforming_mesh_(dU_);
      is_deformed_ = true;
    }

    void undeform(){
      if (!is_deformed_) return;
      VECTOR dU_inverted(dU_);
      gmm::scale(dU_inverted, scalar_type(-1.0));
      deforming_mesh_(dU_inverted);
      is_deformed_ = false;
    }

    ~temporary_mesh_deformator(){
      if (to_be_restored_ && deform_on_construct_){
        undeform();
      }
    }

  private:
    void deforming_mesh_(VECTOR &dU){
      auto &ppts = m_.points();
      size_type ddim = mf_.get_qdim();
      auto init_nb_points = ppts.card();

      dal::bit_vector conv_indices = mf_.convex_index();
      //this vector will track if a point can be deformed
      std::vector<bool> deform_pt_flag(ppts.size(), true);
      size_type cv;
      for (cv << conv_indices; cv != bgeot::size_type(-1); cv << conv_indices)
      {
        getfem::mesh::ind_cv_ct pt_index = m_.ind_points_of_convex(cv);
        getfem::mesh_fem::ind_dof_ct dof = mf_.ind_basic_dof_of_element(cv);
        bgeot::size_type num_points = m_.structure_of_convex(cv)->nb_points();

        GMM_ASSERT2(dof.size() % num_points == 0,
          "mesh_fem should be isoparametric to the mesh, "
          "with nb_points() of convex == size of ind_basic_dof_of_element / qdim()");


        for (size_type pt = 0; pt < num_points; ++pt)
        {
          /** iterate through each components of point [pt]and deform the component*/
          if (deform_pt_flag[pt_index[pt]])
          for (size_type comp = 0; comp < ddim; ++comp)
            //move pts by dU;
            ppts[pt_index[pt]][comp] += dU[dof[pt*ddim + comp]];

          //flag current [pt] to deformed
          deform_pt_flag[pt_index[pt]] = false;
        }
        ppts.resort();
      }
      GMM_ASSERT1(ppts.card() == init_nb_points, 
                  "Error, after deforming the mesh, number of nodes are different.");
    }

    VECTOR dU_;
    const mesh_fem &mf_;
    mesh &m_;
    bool deform_on_construct_;
    bool is_deformed_;
    bool to_be_restored_;
  };

}//end of getfem namespace

#endif //GETFEM_DEFORMABLE_MESH_H__
