/*===========================================================================

 Copyright (C) 2012-2015 Liang Jin Lim.

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

===========================================================================*/

#include "getfem/getfem_im_data.h"
#include "getfem/getfem_omp.h"

namespace getfem
{

  im_data::im_data(const getfem::mesh_im& mim,
                   bgeot::multi_index tsize,
                   size_type filtered_region_)
    : im_(mim), region_(filtered_region_),
      nb_int_pts_intern(0), nb_int_pts_onfaces(0),
      nb_filtered_int_pts_intern(0), nb_filtered_int_pts_onfaces(0),
      locks_()
  {
    set_tensor_size(tsize);
    add_dependency(im_);
    update_from_context();
  }

  im_data::im_data(const getfem::mesh_im& mim, size_type filtered_region_)
    : im_(mim), region_(filtered_region_),
      nb_int_pts_intern(0), nb_int_pts_onfaces(0),
      nb_filtered_int_pts_intern(0), nb_filtered_int_pts_onfaces(0),
      locks_()
  {
    tensor_size_.resize(1);
    tensor_size_[0] = 1;
    nb_tensor_elem_ = 1;
    add_dependency(im_);
    update_from_context();
  }

  size_type im_data::nb_index(bool use_filter) const {
    context_check();
    if (use_filter)
      return nb_filtered_int_pts_intern + nb_filtered_int_pts_onfaces;
    else
      return nb_int_pts_intern + nb_int_pts_onfaces;
  }

  size_type im_data::nb_points_of_element(size_type cv, bool use_filter) const {
    context_check();
    if (cv < convexes.size()) {
      size_type nb_int_pts(0);
      if (use_filter) {
        for (short_type f=0, nb_faces=nb_faces_of_element(cv);
             f < nb_faces; ++f)
          if (convexes[cv].first_int_pt_onface_fid[f] != size_type(-1))
            nb_int_pts += convexes[cv].nb_int_pts_onface[f];
        if (convexes[cv].first_int_pt_fid != size_type(-1)) // convex in filtered_region()
          nb_int_pts += convexes[cv].nb_int_pts;
      }
      else {
        for (auto nb_pts : convexes[cv].nb_int_pts_onface)
          nb_int_pts += nb_pts;
        if (nb_int_pts_intern > 0) // has_convexes
          nb_int_pts += convexes[cv].nb_int_pts;
      }
      return nb_int_pts;
    }
    return 0;
  }

  size_type im_data::nb_points_of_element(size_type cv, short_type f,
                                          bool use_filter) const {
    context_check();
    if (cv < convexes.size()) {
      if (f == short_type(-1)) {
        if (!use_filter || convexes[cv].first_int_pt_fid != size_type(-1))
          return convexes[cv].nb_int_pts;
      }
      else if (f < convexes[cv].nb_int_pts_onface.size()) {
        if (!use_filter || convexes[cv].first_int_pt_onface_fid[f] != size_type(-1))
          return convexes[cv].nb_int_pts_onface[f];
      }
    }
    return 0;
  }

  short_type im_data::nb_faces_of_element(size_type cv) const {
    context_check();
    if (cv < convexes.size())
      return short_type(convexes[cv].first_int_pt_onface_id.size());
    return 0;
  }

  size_type im_data::index_of_point(size_type cv, size_type i,
                                    bool use_filter) const {
    context_check();
    if (cv < convexes.size()) {
      if (i < convexes[cv].nb_int_pts) { // internal point
        size_type int_pt_id = use_filter ? convexes[cv].first_int_pt_fid
                                         : convexes[cv].first_int_pt_id;
        if (int_pt_id != size_type(-1))
          return int_pt_id + i;
      }
      else if (nb_faces_of_element(cv) != 0) {
        const getfem::papprox_integration pim = approx_int_method_of_element(cv);
        for (short_type f=0, nb_faces=pim->nb_convex_faces();
             f < nb_faces; ++f)
          if (i < pim->repart()[f+1]) {
            size_type int_pt_id = use_filter ? convexes[cv].first_int_pt_onface_fid[f]
                                             : convexes[cv].first_int_pt_onface_id[f];
            if (int_pt_id != size_type(-1))
              return int_pt_id + i - pim->ind_first_point_on_face(f);
            break;
          }
      }
    }
    return size_type(-1);
  }

  size_type im_data::filtered_index_of_point(size_type cv, size_type i) const
  { return index_of_point(cv, i, true); }

  size_type im_data::index_of_first_point(size_type cv, short_type f,
                                          bool use_filter) const {
    context_check();
    if (cv < convexes.size()) {
      if (f == short_type(-1)) {
        return use_filter ? convexes[cv].first_int_pt_fid
                          : convexes[cv].first_int_pt_id;
      }
      else if (f < nb_faces_of_element(cv)) {
        return use_filter ? convexes[cv].first_int_pt_onface_fid[f]
                          : convexes[cv].first_int_pt_onface_id[f];
      }
    }
    return size_type(-1);
  }

  size_type im_data::filtered_index_of_first_point(size_type cv, short_type f) const
  { return index_of_first_point(cv, f, true); }

  dal::bit_vector im_data::convex_index(bool use_filter) const {
    context_check();
    dal::bit_vector ind = im_.convex_index();
    if (use_filter && filtered_region() != size_type(-1))
      ind &= linked_mesh().region(filtered_region()).index();
    return ind;
  }

  void im_data::set_region(size_type rg) {
    region_ = rg;
    update_from_context();
  }

  void im_data::update_from_context() const {
    local_guard lock = locks_.get_lock();
    bool no_region = (filtered_region() == size_type(-1));
    const getfem::mesh_region &rg = no_region
                                  ? getfem::mesh_region::all_convexes()
                                  : im_.linked_mesh().region(filtered_region());
    bool has_convexes = (no_region || !rg.is_only_faces());
    bool has_faces = (!no_region && !rg.is_only_convexes());

    size_type nb_cv = im_.convex_index().last_true() + 1;
    convexes.clear();
    convexes.resize(nb_cv);

    for (dal::bv_visitor cv(im_.convex_index()); !cv.finished(); ++cv)
      convexes[cv].nb_int_pts
        = im_.int_method_of_element(cv)->approx_method()->nb_points_on_convex();

    nb_int_pts_intern = 0;
    nb_filtered_int_pts_intern = 0;
    if (has_convexes)
      // The global indexing of integration points follows the indexing of convexes
      for (dal::bv_visitor cv(im_.convex_index()); !cv.finished(); ++cv) {
        convexes[cv].first_int_pt_id = nb_int_pts_intern;
        nb_int_pts_intern += convexes[cv].nb_int_pts;
        if (no_region || rg.is_in(cv)) {
          convexes[cv].first_int_pt_fid = nb_filtered_int_pts_intern;
          nb_filtered_int_pts_intern += convexes[cv].nb_int_pts;
        }
      }

    nb_int_pts_onfaces = 0;
    nb_filtered_int_pts_onfaces = 0;
    if (has_faces) {
      for (dal::bv_visitor cv(im_.convex_index()); !cv.finished(); ++cv) {
        short_type nb_faces = im_.linked_mesh().nb_faces_of_convex(cv);
        convexes[cv].first_int_pt_onface_id.assign(nb_faces, size_type(-1));
        convexes[cv].first_int_pt_onface_fid.assign(nb_faces, size_type(-1));
        convexes[cv].nb_int_pts_onface.assign(nb_faces, size_type(-1));
        const getfem::papprox_integration pim(approx_int_method_of_element(cv));
        for (short_type f = 0; f < nb_faces; ++f) {
          convexes[cv].first_int_pt_onface_id[f] = nb_int_pts_intern +
                                                   nb_int_pts_onfaces;
          size_type nb_pts = pim->nb_points_on_face(f);
          nb_int_pts_onfaces += nb_pts;
          if (rg.is_in(cv, f)) {
            convexes[cv].first_int_pt_onface_fid[f] = nb_filtered_int_pts_intern +
                                                      nb_filtered_int_pts_onfaces;
            nb_filtered_int_pts_onfaces += nb_pts;
          }
          convexes[cv].nb_int_pts_onface[f] = nb_pts;
        }
      }
    }
    v_num_ = act_counter();
    touch();
  }

  size_type im_data::nb_tensor_elem() const {
    return nb_tensor_elem_;
  }

  void im_data::set_tensor_size(const bgeot::multi_index& tsize) {
    tensor_size_ = tsize;    
    nb_tensor_elem_ = tensor_size_.total_size();
  }

  bool is_equivalent_with_vector(const bgeot::multi_index &sizes, size_type vector_size) {
    bool checked = false;
    size_type size = 1;
    for (size_type i = 0; i < sizes.size(); ++i) {
      if (sizes[i] > 1 && checked) return false;
      if (sizes[i] > 1) { 
        checked = true;
        size = sizes[i];
        if (size != vector_size) return false;
      }
    }
    return (vector_size == size);
  }

  bool is_equivalent_with_matrix(const bgeot::multi_index &sizes, size_type nrows, size_type ncols) {
    if (nrows == 1 || ncols == 1) {
      return is_equivalent_with_vector(sizes, nrows + ncols - 1);
    }    
    size_type tensor_row = 1;
    size_type tensor_col = 1;
    bool first_checked = false;
    bool second_checked = false;
    for (size_type i = 0; i < sizes.size(); ++i) {
      if (sizes[i] > 1 && !first_checked) {
        first_checked = true;
        tensor_row = sizes[i];
        if (tensor_row != nrows) return false;
      } else if (sizes[i] > 1 && !second_checked) {
        second_checked = true;
        tensor_col = sizes[i];
        if (tensor_col != ncols) return false;
      }
      else if (sizes[i] > 1 && first_checked && second_checked) return false;
    }
    return (nrows == tensor_row && ncols == tensor_col);
  }
}
