#include "getfem/getfem_im_data.h"

namespace getfem
{

  im_data::im_data (const getfem::mesh_im& meshIm, bgeot::multi_index tensorSize,
    size_type filteredRegion)
    :im_(meshIm), filtered_region_(filteredRegion), 
    nb_filtered_index_(0), nb_index_(0), require_update_(true){
    set_tensor_size(tensorSize);
    add_dependency(im_);
    update_index_();
  }

  im_data::im_data (const getfem::mesh_im& meshIm, size_type filteredRegion)
    :im_(meshIm), filtered_region_(filteredRegion), 
    nb_filtered_index_(0), nb_index_(0), require_update_(true){
    tensor_size_.resize(1);
    tensor_size_[0] = 1;
    nb_tensor_elem_ = 1;
    add_dependency(im_);
    update_index_();
  }

  void im_data::update_index_() const{
    nb_index_         = 0;      
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

    require_update_ = false;
    touch();
    v_num_ = act_counter();
  }

  size_type im_data::nb_index() const{
    context_check();
    if (require_update_) update_index_();
    return nb_index_;
  }

  size_type im_data::nb_filtered_index() const{
    context_check();
    if (require_update_) update_index_();
    return nb_filtered_index_;
  }

  size_type im_data::nb_points_of_element(size_type cv) const{
    return im_.int_method_of_element(cv)->approx_method()->nb_points_on_convex();
  }

  size_type im_data::index_of_point(size_type cv, size_type i) const{
    context_check();
    if (require_update_) update_index_();
    if(cv < int_point_index_.size()) return int_point_index_[cv] + i;
    else return size_type(-1);
  }

  size_type im_data::nb_tensor_elem() const{
    return nb_tensor_elem_;
  }

  void im_data::set_tensor_size (const bgeot::multi_index& tensorSize){
    tensor_size_  = tensorSize;    
    nb_tensor_elem_ = 1;
    auto it      = tensor_size_.begin();
    auto itEnd   = tensor_size_.end();
    for(; it != itEnd; ++it) nb_tensor_elem_ *= *it;
  }

  size_type im_data::filtered_index_of_point(size_type cv, size_type i) const{
    context_check();
    if (require_update_) update_index_();
    if(cv < filtered_int_point_index_.size())
    {
      if (filtered_int_point_index_[cv] == size_type(-1)) return size_type(-1);
      return filtered_int_point_index_[cv] + i;
    }
    else return size_type(-1);
  }

  dal::bit_vector im_data::filtered_convex_index() const{
    context_check();
    if (require_update_) update_index_();
    return filtered_convex_index_;
  }

  std::vector<size_type> im_data::filtered_index_of_first_point() const{
    context_check();
    if (require_update_) update_index_();
    return filtered_int_point_index_;
  }

  void im_data::set_region(size_type rg){
    filtered_region_ = rg;
    update_index_();
  }

  void im_data::update_from_context() const{
    require_update_ = true;
  }
}
