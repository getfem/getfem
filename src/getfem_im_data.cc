#include "getfem/getfem_im_data.h"
#include "getfem/getfem_omp.h"

namespace getfem
{

  im_data::im_data(const getfem::mesh_im& meshIm,
                   bgeot::multi_index tensorSize,
                   size_type filteredRegion)
    :im_(meshIm), nb_filtered_index_(0), nb_index_(0),
    filtered_region_(filteredRegion), is_on_face_(false),
    locks_()
  {
    set_tensor_size(tensorSize);
    add_dependency(im_);
    if (filtered_region_!=size_type(-1) && 
        im_.linked_mesh().region(filtered_region_).is_only_faces()){
        update_index_for_face_();
        is_on_face_ = true;
    }
    else update_index_();
  }

  im_data::im_data(const getfem::mesh_im& meshIm, size_type filteredRegion)
    :im_(meshIm), nb_filtered_index_(0), nb_index_(0),
    filtered_region_(filteredRegion), is_on_face_(false) {
    tensor_size_.resize(1);
    tensor_size_[0] = 1;
    nb_tensor_elem_ = 1;
    add_dependency(im_);
    if (filtered_region_ != size_type(-1) &&
      im_.linked_mesh().region(filtered_region_).is_only_faces()){
      update_index_for_face_();
      is_on_face_ = true;
    }
    else update_index_();
  }
    
  void im_data::update_index_() const {
    local_guard lock = locks_.get_lock();
    nb_index_         = 0;
    face_container_.clear();
    filtered_face_container_.clear();

    size_type nElement = im_.convex_index().last_true() + 1;
    int_point_index_.clear();
    int_point_index_.resize(nElement, size_type(-1));
    
    nb_filtered_index_ = 0;
    filtered_int_point_index_.clear();
    filtered_int_point_index_.resize(nElement, size_type(-1));
    filtered_convex_index_.clear();

    for(dal::bv_visitor cv(im_.convex_index()); !cv.finished(); ++cv)
    {
      size_type nPoint = nb_points_of_element(cv);


      int_point_index_[cv] = nb_index_;
      nb_index_          += nPoint;

      if(filtered_region_ == size_type(-1) 
        || im_.linked_mesh().region(filtered_region_).is_in(cv))
      {
        filtered_convex_index_.add(cv);
        filtered_int_point_index_[cv]  = nb_filtered_index_;
        nb_filtered_index_            += nPoint;
      }
    }
    v_num_ = act_counter();
  }

  void im_data::update_index_for_face_() const {
    local_guard lock = locks_.get_lock();
    nb_index_ = 0;
    int_point_index_.clear();
    filtered_int_point_index_.clear();

    size_type nElement = im_.convex_index().last_true() + 1;
    face_container_.clear();
    face_container_.resize(nElement);

    nb_filtered_index_ = 0;
    filtered_face_container_.clear();
    filtered_face_container_.resize(nElement);
    filtered_convex_index_.clear();

    const getfem::mesh &linked_mesh = im_.linked_mesh();
    const getfem::mesh_region &rg = linked_mesh.region(filtered_region_);
    
    for (dal::bv_visitor cv(im_.convex_index()); !cv.finished(); ++cv)
    {
      short_type nb_faces = linked_mesh.nb_faces_of_convex(cv);      
      const getfem::papprox_integration pim = im_.int_method_of_element(cv)->approx_method();
      size_type nb_points = pim->nb_points_on_convex();
      size_type nb_points_on_faces = pim->nb_points() - nb_points;
      face_container_[cv].index_of_face_.resize(nb_faces);
      face_container_[cv].inactive_ = false;
      face_container_[cv].start_index_ = nb_points;
      face_container_[cv].point_index_.resize(nb_points_on_faces);
      size_type nb_index_convex = 0;

      for (short_type i_face = 0; i_face < nb_faces; ++i_face){
        size_type nb_pt_face = pim->nb_points_on_face(i_face);
        for (size_type k = 0; k < nb_pt_face; ++k){
          face_container_[cv].point_index_[nb_index_convex] = nb_index_;
          face_container_[cv].index_of_face_[i_face].add(nb_index_);
          ++nb_index_;
          ++nb_index_convex;
        }
      }
      if (rg.index().is_in(cv)){
        filtered_convex_index_.add(cv);
        filtered_face_container_[cv].index_of_face_.resize(nb_faces);
        filtered_face_container_[cv].inactive_ = false;
        filtered_face_container_[cv].start_index_ = nb_points - nb_points;
        filtered_face_container_[cv].point_index_.resize(nb_points_on_faces);
        nb_index_convex = 0;

        for (short_type i_face = 0; i_face < nb_faces; ++i_face){
          size_type nb_pt_face = pim->nb_points_on_face(i_face);
          for (size_type k = 0; k < nb_pt_face; ++k){
            if (!rg.is_in(cv, i_face)){
              filtered_face_container_[cv].point_index_[nb_index_convex] = size_type(-1);
            }
            else{
              filtered_face_container_[cv].point_index_[nb_index_convex] = nb_filtered_index_;
              filtered_face_container_[cv].index_of_face_[i_face].add(nb_filtered_index_);
              ++nb_filtered_index_;
            }
            ++nb_index_convex;
          }
        }
      }
    }
    v_num_ = act_counter();
  }

  size_type im_data::nb_index() const{
    context_check();
    return nb_index_;
  }

  size_type im_data::nb_filtered_index() const{
    context_check();
    return nb_filtered_index_;
  }

  size_type im_data::nb_points_of_element(size_type cv, short_type f) const{
    context_check();
    if (!im_.convex_index().is_in(cv)) return 0;
    if (f == short_type(-1)){
      size_type nb_points = im_.int_method_of_element(cv)
                              ->approx_method()->nb_points_on_convex();
      if (!is_on_face_) return nb_points;
      else return im_.int_method_of_element(cv)->approx_method()->nb_points() - nb_points;

    }
    GMM_ASSERT1(is_on_face_, "im_data is not defined on face.");
    return im_.int_method_of_element(cv)->approx_method()->nb_points_on_face(f);
  }

  size_type im_data::index_of_point(size_type cv, size_type i) const{
    context_check();
    if (!is_on_face_){
      if (cv < int_point_index_.size()) return int_point_index_[cv] + i;
      else return size_type(-1);
    }
    if (cv < face_container_.size()){
      size_type return_index = i - face_container_[cv].start_index_;
      if (face_container_[cv].inactive_ ||
          face_container_[cv].start_index_ > i ||
          face_container_[cv].point_index_.size() < return_index){
        return size_type(-1);
      }
      return face_container_[cv].point_index_[return_index];
    }
    return size_type(-1);
  }

  size_type im_data::nb_tensor_elem() const{
    return nb_tensor_elem_;
  }

  void im_data::set_tensor_size(const bgeot::multi_index& tensorSize){
    tensor_size_  = tensorSize;    
    nb_tensor_elem_ = 1;
    bgeot::multi_index::iterator it = tensor_size_.begin();
    bgeot::multi_index::iterator itEnd = tensor_size_.end();
    for(; it != itEnd; ++it) nb_tensor_elem_ *= *it;
  }

  size_type im_data::filtered_index_of_point(size_type cv, size_type i) const{
    context_check();
    if (!is_on_face_){
      if(cv < filtered_int_point_index_.size())
      {
        if (filtered_int_point_index_[cv] == size_type(-1)) return size_type(-1);
        return filtered_int_point_index_[cv] + i;
      }
      else return size_type(-1);
    }
    if (cv < filtered_face_container_.size()){
      size_type return_index = i - face_container_[cv].start_index_;
      if (filtered_face_container_[cv].inactive_ ||
          filtered_face_container_[cv].start_index_ > i ||
          filtered_face_container_[cv].point_index_.size() < return_index){
        return size_type(-1);
      }
      return filtered_face_container_[cv].point_index_[return_index];
    }
    return size_type(-1);
  }

  dal::bit_vector im_data::filtered_convex_index() const{
    context_check();
    return filtered_convex_index_;
  }

  std::vector<size_type> im_data::filtered_index_of_first_point() const{
    context_check();
    return filtered_int_point_index_;
  }

  void im_data::set_region(size_type rg) const{
    filtered_region_ = rg;
    touch();
  }

  void im_data::update_from_context() const{
    is_on_face_ = false;
    if (filtered_region_ != size_type(-1) && 
        im_.linked_mesh().region(filtered_region_).is_only_faces()){
      update_index_for_face_();
      is_on_face_ = true;
    }
    else update_index_();
    touch();
  }

  bool is_equivalent_with_vector(const bgeot::multi_index &sizes, size_type vector_size){
    bool checked = false;
    size_type size = 1;
    for (size_type i = 0; i < sizes.size(); ++i){
      if (sizes[i] > 1 && checked) return false;
      if (sizes[i] > 1){ 
        checked = true;
        size = sizes[i];
        if (size != vector_size) return false;
      }
    }
    return (vector_size == size);
  }

  bool is_equivalent_with_matrix(const bgeot::multi_index &sizes, size_type nrows, size_type ncols){
    if (nrows == 1 || ncols == 1){
      return is_equivalent_with_vector(sizes, nrows + ncols - 1);
    }    
    size_type tensor_row = 1;
    size_type tensor_col = 1;
    bool first_checked = false;
    bool second_checked = false;
    for (size_type i = 0; i < sizes.size(); ++i){
      if (sizes[i] > 1 && !first_checked){
        first_checked = true;
        tensor_row = sizes[i];
        if (tensor_row != nrows) return false;
      }else if (sizes[i] > 1 && !second_checked){
        second_checked = true;
        tensor_col = sizes[i];
        if (tensor_col != ncols) return false;
      }
      else if (sizes[i] > 1 && first_checked && second_checked) return false;
    }
    return (nrows == tensor_row && ncols == tensor_col);
  }
}
