/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2012-2015 Liang Jin Lim

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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
@file getfem_im_data.h
@brief Provides indexing of integration points for mesh_im.
@date Feb 2014
@author Liang Jin Lim
*/

#pragma once

#ifndef GETFEM_IM_DATA_H__
#define GETFEM_IM_DATA_H__

#include <getfem/getfem_mesh_im.h>

namespace getfem{
  using bgeot::size_type;
  using bgeot::scalar_type;

  /**check if a given tensor size is equivalent to a vector*/
  bool is_equivalent_with_vector(const bgeot::multi_index &sizes, size_type vector_size);
  /**check if a given tensor size is equivalent to a matrix*/
  bool is_equivalent_with_matrix(const bgeot::multi_index &sizes, size_type nrows, size_type ncols);

  /** im_data provides indexing to the integration points of a mesh
  im object. The im_data data contains a reference of mesh_im object
  . The index can be filtered by region, and each im_data has its 
  own tensorial size.

  Filtered methods will provide filtered index on the region.
  This class also provides reading and writing tensor( including
  matrix, vector and scalar) from a vector data (generally a
  fixed-size variable from the model.)
  
  im_data can be used to provide integration point index on convex or
  on faces of convex, but not both. To create an im_data that represents
  integration points on a face, the filter region provided has to contain
  only faces.
  */
  class im_data : public context_dependencies, virtual public dal::static_stored_object {
  public:
    /**
    * Constructor
    * @param mim Reference mesh_im object
    * @param tensor_size tensor dimension of each integration points
    * @param filtered_region index not in the region will be filtered
    *        out. If filtered_region can contain only convexes or only
    *        faces or both convexes and faces.
    */
    im_data(const mesh_im& mim_, bgeot::multi_index tensor_size,
            size_type filtered_region_ = size_type(-1));

    /**
    * Constructor. The tensor size by default is a scalar value.
    * @param mim Reference mesh_im object
    * @param filtered_region index not in the region will be filtered
    *        out. If filtered_region can contain only convexes or only
    *        faces or both convexes and faces.
    */
    im_data(const mesh_im& mim_, size_type filtered_region_ = size_type(-1));

    /**set filtered region id*/
    void set_region(size_type region);

    /**return filtered region id*/
    inline size_type filtered_region() const {return region_;}

    /**Returns the index of an integration point with no filtering*/
    size_type index_of_point(size_type cv, size_type i,
                             bool use_filter = false) const;

    /**Returns the index of an integration point with filtering*/
    size_type filtered_index_of_point(size_type cv, size_type i) const;

    /**Returns the index of the first integration point with no filtering*/
    size_type index_of_first_point(size_type cv, short_type f = short_type(-1),
                                   bool use_filter = false) const;

    /**Returns the index of the first integration point with filtering*/
    size_type filtered_index_of_first_point(size_type cv, short_type f = short_type(-1)) const;

    /**Total numbers of index (integration points)*/
    size_type nb_index(bool use_filter=false) const;

    /**Total numbers of filtered index (integration points)*/
    size_type nb_filtered_index() const
    { return nb_index(true); }

    /**Total number of points in element cv*/
    size_type nb_points_of_element(size_type cv, bool use_filter=false) const;

    /**Number of points in element cv, on face f (or in the interior)*/
    size_type nb_points_of_element(size_type cv, short_type f, bool use_filter=false) const;

    /**Total number of points in element cv, which lie in filtered_region()*/
    size_type nb_filtered_points_of_element(size_type cv) const
    { return nb_points_of_element(cv, true); }

    /**Number of points in element cv, on face f (or in the interior),
    which lie in filtered_region()*/
    size_type nb_filtered_points_of_element(size_type cv, short_type f) const
    { return nb_points_of_element(cv, f, true); }

    /**Number of (active) faces in element cv*/
    short_type nb_faces_of_element(size_type cv) const;

    /**sum of tensor elements, M(3,3) will have 3*3=9 elements*/
    size_type nb_tensor_elem() const;

    /**List of convexes*/
    dal::bit_vector convex_index(bool use_filter=false) const;

    /**List of convex in filtered region*/
    dal::bit_vector filtered_convex_index() const
    { return convex_index(true); }

    /**called automatically when there is a change in dependencies*/
    void update_from_context () const;

    /**linked mesh im*/
    inline const mesh_im &linked_mesh_im() const { return im_; }

    /**implicit conversion to mesh im*/
    inline operator const mesh_im &() const { return im_; }

    /**linked mesh*/
    inline const mesh &linked_mesh () const { return im_.linked_mesh(); }

    getfem::papprox_integration approx_int_method_of_element(size_type cv) const
    { return im_.int_method_of_element(cv)->approx_method(); }

    inline const bgeot::multi_index& tensor_size () const { return tensor_size_;}

    void set_tensor_size (const bgeot::multi_index& tensor_size);

    inline gmm::uint64_type version_number() const { context_check(); return v_num_; }

    /**Extend a vector from filtered size to full size and copy the data to correct index*/
    template <typename VECT>
    void extend_vector(const VECT &V1, VECT &V2) const {
      if (V1.size() == 0 && V2.size() == 0)
        return;
      size_type nb_data = V1.size()/nb_filtered_index();
      GMM_ASSERT1(V1.size() == nb_data*nb_filtered_index(), "Invalid size of vector V1");
      GMM_ASSERT1(V2.size() == nb_data*nb_index(), "Invalid size of vector V2");
      if (nb_filtered_index() == nb_index()) {
        gmm::copy(V1, V2);
        return;
      }

      const getfem::mesh_region &rg = im_.linked_mesh().region(filtered_region());
      for (getfem::mr_visitor v(rg); !v.finished(); ++v) {
        size_type nb_pts = nb_points_of_element(v.cv(), v.f());
        size_type first_id = (v.f() == short_type(-1))
                           ? convexes[v.cv()].first_int_pt_id
                           : convexes[v.cv()].first_int_pt_onface_id[v.f()];
        size_type first_fid = (v.f() == short_type(-1))
                            ? convexes[v.cv()].first_int_pt_fid
                            : convexes[v.cv()].first_int_pt_onface_fid[v.f()];
        if (first_fid != size_type(-1))
          gmm::copy(
            gmm::sub_vector(V1, gmm::sub_interval(first_fid*nb_data, nb_pts*nb_data)),
            gmm::sub_vector(V2, gmm::sub_interval(first_id*nb_data, nb_pts*nb_data)));
      }
    }

    /**Filter a vector from full size to filtered size and copy the data to correct index*/    
    template <typename VECT>
    void reduce_vector(const VECT &V1, VECT &V2) const {
      if (V1.size() == 0 && V2.size() == 0)
        return;
      size_type nb_data = V1.size()/nb_index();
      GMM_ASSERT1(V1.size() == nb_data*nb_index(), "Invalid size of vector V1");
      GMM_ASSERT1(V2.size() == nb_data*nb_filtered_index(), 
                               "Invalid size of vector V2");
      if (nb_filtered_index() == nb_index()) {
        gmm::copy(V1, V2);
        return;
      }

      const getfem::mesh_region &rg = im_.linked_mesh().region(filtered_region());
      for (getfem::mr_visitor v(rg); !v.finished(); ++v) {
        size_type nb_pts = nb_points_of_element(v.cv(), v.f());
        size_type first_id = (v.f() == short_type(-1))
                           ? convexes[v.cv()].first_int_pt_id
                           : convexes[v.cv()].first_int_pt_onface_id[v.f()];
        size_type first_fid = (v.f() == short_type(-1))
                            ? convexes[v.cv()].first_int_pt_fid
                            : convexes[v.cv()].first_int_pt_onface_fid[v.f()];
        if (first_fid != size_type(-1))
          gmm::copy(
            gmm::sub_vector(V1, gmm::sub_interval(first_id*nb_data, nb_pts*nb_data)),
            gmm::sub_vector(V2, gmm::sub_interval(first_fid*nb_data, nb_pts*nb_data)));
      }
    }

    /**get a scalar value of an integration point 
    from a raw vector data, described by the tensor size.*/
    template <typename VECT>
    typename VECT::value_type get_value(const VECT &V1, size_type cv,
                                        size_type i, bool use_filter = true) const {
      GMM_ASSERT1(nb_tensor_elem_*nb_index(use_filter) == V1.size(),
                  "Invalid tensorial size for vector V1");
      GMM_ASSERT1(nb_tensor_elem_ == 1, "im_data is not of scalar type");
      size_type ptid = index_of_point(cv,i,use_filter);
      GMM_ASSERT2(ptid != size_type(-1), "Point index of gauss point not found");
      return V1[ptid];
    }

    /**get a vector of an integration point 
    from a raw vector data, described by the tensor size.*/
    template <typename VECT1, typename VECT2>
    void get_vector(const VECT1 &V1, size_type cv, size_type i,
                    VECT2& V2, bool use_filter = true) const {
      if (V1.size() == 0 && V2.size() == 0)
        return;
      GMM_ASSERT1(nb_tensor_elem_*nb_index(use_filter) == V1.size(),
                  "Invalid tensorial size for vector V1");
      GMM_ASSERT1(is_equivalent_with_vector(tensor_size_, V2.size()),
                  "V2 is incompatible with im_data tensor size");
      size_type ptid = index_of_point(cv,i,use_filter);
      GMM_ASSERT2(ptid != size_type(-1), "Point index of gauss point not found");
      gmm::copy(gmm::sub_vector(V1, gmm::sub_interval(ptid*nb_tensor_elem_,
                                                      nb_tensor_elem_)),
                V2);
    }

    /**get a matrix of an integration point 
    from a raw vector data, described by the tensor size.*/
    template <typename VECT, typename MAT>
    void get_matrix(const VECT &V1, size_type cv, size_type i, 
                    MAT& M, bool use_filter = true) const {
      if (V1.size() == 0 && M.size() == 0)
        return;
      GMM_ASSERT1(nb_tensor_elem_*nb_index(use_filter) == V1.size(),
                  "Invalid tensorial size for vector V1");
      GMM_ASSERT1(is_equivalent_with_matrix(tensor_size_, M.nrows(), M.ncols()),
                  "M is incompatible with im_data tensor size");
      size_type ptid = index_of_point(cv,i,use_filter);
      GMM_ASSERT2(ptid != size_type(-1), "Point index of gauss point not found");
      gmm::copy(gmm::sub_vector(V1, gmm::sub_interval(ptid*nb_tensor_elem_,
                                                      nb_tensor_elem_)),
                M.as_vector());
    }

    /**get a tensor of an integration point 
    from a raw vector data, described by the tensor size.*/
    template <typename VECT, typename TENSOR>
    void get_tensor(const VECT &V1, size_type cv, size_type i, 
                    TENSOR& T, bool use_filter = true) const {
      if (V1.size() == 0 && T.size() == 0)
        return;
      GMM_ASSERT1(nb_tensor_elem_*nb_index(use_filter) == V1.size(),
                  "Invalid tensorial size for vector V1");
      GMM_ASSERT1(tensor_size_ == T.sizes(),
                  "T is incompatible with im_data tensor size");
      size_type ptid = index_of_point(cv,i,use_filter);
      GMM_ASSERT2(ptid != size_type(-1), "Point index of gauss point not found");
      gmm::copy(gmm::sub_vector(V1, gmm::sub_interval(ptid*nb_tensor_elem_,
                                                      nb_tensor_elem_)),
                T.as_vector());
    }

    /**set a value of an integration point 
    from a raw vector data, described by the tensor size.*/
    template <typename VECT>
    typename VECT::value_type &set_value(VECT &V1, size_type cv, size_type i,
                                         bool use_filter = true) const {
      GMM_ASSERT1(nb_tensor_elem_*nb_index(use_filter) == V1.size(),
                  "Invalid tensorial size for vector V1");
      GMM_ASSERT1(nb_tensor_elem_ == 1, "im_data is not of scalar type");
      size_type ptid = index_of_point(cv,i,use_filter);
      GMM_ASSERT2(ptid != size_type(-1), "Point index of gauss point not found");
      return V1[ptid];
    }

    /**set a vector of an integration point 
    from a raw vector data, described by the tensor size.*/
    template <typename VECT1, typename VECT2>
    void set_vector(VECT1 &V1, size_type cv, size_type i, 
                    const VECT2& V2, bool use_filter = true) const {
      if (V1.size() == 0 && V2.size() == 0)
        return;
      GMM_ASSERT1(nb_tensor_elem_*nb_index(use_filter) == V1.size(),
                  "Invalid tensorial size for vector V1");
      GMM_ASSERT1(is_equivalent_with_vector(tensor_size_, V2.size()),
                  "V2 is incompatible with im_data tensor size");
      size_type ptid = index_of_point(cv,i,use_filter);
      GMM_ASSERT2(ptid != size_type(-1), "Point index of gauss point not found");
      gmm::copy(V2,
                gmm::sub_vector(V1, gmm::sub_interval(ptid*nb_tensor_elem_,
                                                      nb_tensor_elem_)));
    }

    /**set a matrix of an integration point 
    from a raw vector data, described by the tensor size.*/
    template <typename VECT, typename MAT>
    void set_matrix(VECT &V1, size_type cv, size_type i, 
                    const MAT& M, bool use_filter = true) const {
      if (V1.size() == 0 && M.size() == 0)
        return;
      GMM_ASSERT1(nb_tensor_elem_*nb_index(use_filter) == V1.size(),
                  "Invalid tensorial size for vector V1");
      GMM_ASSERT1(is_equivalent_with_matrix(tensor_size_, M.nrows(), M.ncols()),
                  "M is incompatible with im_data tensor size");
      size_type ptid = index_of_point(cv,i,use_filter);
      GMM_ASSERT2(ptid != size_type(-1), "Point index of gauss point not found");
      gmm::copy(M.as_vector(),
                gmm::sub_vector(V1, gmm::sub_interval(ptid*nb_tensor_elem_,
                                                      nb_tensor_elem_)));
    }

    /**set a tensor of an integration point 
    from a raw vector data, described by the tensor size.*/
    template <typename VECT, typename TENSOR>
    void set_tensor(VECT &V1, size_type cv, size_type i,
                    const TENSOR& T, bool use_filter = true) const {
      if (V1.size() == 0 && T.size() == 0)
        return;
      GMM_ASSERT1(nb_tensor_elem_*nb_index(use_filter) == V1.size(),
                  "Invalid tensorial size for vector V1");
      GMM_ASSERT1(tensor_size_ == T.sizes(),
                  "T is incompatible with im_data tensor size");
      size_type ptid = index_of_point(cv,i,use_filter);
      GMM_ASSERT2(ptid != size_type(-1), "Point index of gauss point not found");
      gmm::copy(T.as_vector(),
                gmm::sub_vector(V1, gmm::sub_interval(ptid*nb_tensor_elem_,
                                                      nb_tensor_elem_)));
    }


    template <typename VECT1, typename VECT2>
    void set_vector(VECT1 &V1, size_type ptid, const VECT2& V2) const {
      GMM_ASSERT1(V1.size() != 0, "V1 of zero size");
      GMM_ASSERT1(V2.size() != 0, "V2 of zero size");
      GMM_ASSERT1(is_equivalent_with_vector(tensor_size_, V2.size()),
                  "V2 is incompatible with im_data tensor size");
      gmm::copy(V2,
                gmm::sub_vector(V1, gmm::sub_interval(ptid*nb_tensor_elem_,
                                                      nb_tensor_elem_)));
    }

    template <typename VECT1, typename TENSOR>
    void set_tensor(VECT1 &V1, size_type ptid, const TENSOR& T) const {
      GMM_ASSERT1(V1.size() != 0, "V1 of zero size");
      GMM_ASSERT1(T.size() != 0, "V2 of zero size");
      GMM_ASSERT1(tensor_size_ == T.sizes(),
                  "T is incompatible with im_data tensor size");
      gmm::copy(T.as_vector(),
                gmm::sub_vector(V1, gmm::sub_interval(ptid*nb_tensor_elem_,
                                                      nb_tensor_elem_)));
    }

  private:
    const mesh_im &im_;

    size_type              region_;

    mutable size_type      nb_int_pts_intern;
    mutable size_type      nb_int_pts_onfaces;
    mutable size_type      nb_filtered_int_pts_intern;
    mutable size_type      nb_filtered_int_pts_onfaces;

    struct convex_data {
      size_type first_int_pt_id;   // index
      size_type first_int_pt_fid;  // filtered index
      size_type nb_int_pts;        // number of internal integration points
      std::vector<size_type> first_int_pt_onface_id;
      std::vector<size_type> first_int_pt_onface_fid;
      std::vector<size_type> nb_int_pts_onface;

      convex_data()
        : first_int_pt_id(-1), first_int_pt_fid(-1), nb_int_pts(0),
          first_int_pt_onface_id(0), first_int_pt_onface_fid(0), nb_int_pts_onface(0)
      {}
    };

    mutable std::vector<convex_data> convexes;

    mutable gmm::uint64_type v_num_;

    bgeot::multi_index     tensor_size_;
    size_type              nb_tensor_elem_;
    lock_factory           locks_;
  };
}
#endif /* GETFEM_IM_DATA_H__  */
