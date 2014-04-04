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
@brief This is a normal mesh, whith one extra method, allowing to displace points.
*/

#pragma once
#include <getfem/getfem_mesh.h>
#include <getfem/getfem_mesh_fem.h>
#include <getfem/getfem_models.h>

namespace getfem {

  template<class VECTOR> class temporary_mesh_deformator;


  /** This is a normal mesh, whith one extra method, allowing to displace points.       
  The mesh can only be deformed by instance of class temporary_mesh_deformator, 
  that restores the mesh on it's (deformator) destruction
  */
  class deformable_mesh : public mesh {
  public:
    mutable bool must_be_restored;
  private:
    template <class VECTOR> friend class temporary_mesh_deformator;

    /** says that if the mesh was deformed, it should be deformed back to the 
    underformed state, as other bricks don't know they are dealing with a deformed mesh
    This mesh is used in Updated Lagrane based formulations, but the restore feature
    allows to use it with Total Lagrange as well*/
    inline bool to_be_restored() const {return must_be_restored;}

    /**displace the points by a given displacement vector
    @param U displacement vector as described by mf using dof index (NOT pts index)
    @param &mf mesh_fem object that corresponds to &U, should be compatible with the mesh
    */
    template<typename VEC>
    void deform_mesh(const VEC &dU, const mesh_fem& mf)
    {   
      PT_TAB& ppts = points();
      size_type ddim = ppts.dim();

      GMM_ASSERT1((&mf.linked_mesh())==this,"in deform_mesh mf should be defined on the same mesh");

      GMM_ASSERT1(mf.get_qdim() == ddim, "input mesh_fem and the mesh dim are not compatible");

      dal::bit_vector conv_indices = mf.convex_index(); 
      //this vector will track if a point can be deformed
      std::vector<bool> deform_pt_flag(ppts.size(), true);
      size_type cv;
      for(cv << conv_indices; 
        cv!=bgeot::size_type(-1); cv << conv_indices) 
      {
        getfem::mesh::ind_cv_ct pt_index
          =  mf.linked_mesh().ind_points_of_convex(cv);
        getfem::mesh_fem::ind_dof_ct dof=mf.ind_basic_dof_of_element(cv);
        bgeot::size_type num_points = 
          mf.linked_mesh().structure_of_convex(cv)->nb_points(); 

        GMM_ASSERT2(dof.size() == num_points*ddim, 
          "mesh_fem should be isoparametric to the mesh, "
          "with nb_points() of convex * dim == size of ind_basic_dof_of_element");

        for(size_type pt = 0; pt < num_points; ++pt) 
        { 
          /** iterate through each components of point [pt]and deform the component*/
          if(deform_pt_flag[pt_index[pt]])
            for (size_type comp = 0; comp < ddim; ++comp)
              //move pts by dU;
                ppts[pt_index[pt]][comp] += dU[dof[pt*ddim + comp]];

          //flag current [pt] to deformed
          deform_pt_flag[pt_index[pt]] = false;
        }
        ppts.resort();
      }
    }

  public:

    deformable_mesh(bool _must_be_restored = true, const std::string &name = std::string());
    deformable_mesh(const deformable_mesh&);
  };

  /**cast a conventional mesh into deformable one and remove the const*/
  deformable_mesh& make_deformable_mesh(const mesh&);


  /** An object function that first deformes and then remembers 
  to restore a deformable mesh if it has to be restored 
  for other bricks. By default the mesh is deformed on 
  construct and undeformed in the destructor (by RAII principle)
  but it's also possible to specify deform_on_construct = false 
  and then call explicitely deform() and undeform() methods
  */
  template<class VECTOR = model_real_plain_vector> 
  class temporary_mesh_deformator
  {
    VECTOR dU;
    const mesh_fem& mf;
    deformable_mesh& m;
    bool deform_on_construct_;
    bool is_deformed_;
  public:
    temporary_mesh_deformator(const mesh& _m, const mesh_fem &_mf, 
      const VECTOR &_dU, bool deform_on_construct = true) : 
      dU(_dU),
      mf(_mf), 
      m(make_deformable_mesh(_m)),
      deform_on_construct_(deform_on_construct),
      is_deformed_(false)
    {
      if (deform_on_construct_)
      { 
        m.deform_mesh(dU,mf);
        is_deformed_ = true;
      }
    }

    void deform() 
    {
      GMM_ASSERT1(!is_deformed_, "trying to deformed an already deformed mesh");
      m.deform_mesh(dU,mf);
      is_deformed_ = true;
    }

    void undeform() 
    {
      GMM_ASSERT1(is_deformed_, "trying to restored yet undeformed mesh");
      GMM_ASSERT1(!deform_on_construct_, 
        "undeformed() should not be called if it was asked to deformed on construct"); 
      m.deform_mesh(gmm::scaled(dU,scalar_type(-1.0)),mf);
      is_deformed_ = false;
    }

    ~temporary_mesh_deformator()
    {
      if (m.to_be_restored() && deform_on_construct_)
      {
        m.deform_mesh(gmm::scaled(dU,scalar_type(-1.0)),mf); 
      }
    }
  };

}//end of getfem namespace
