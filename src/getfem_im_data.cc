#include "getfem/getfem_im_data.h"
#include "getfem/getfem_omp.h"

namespace getfem
{

  im_data::im_data(const getfem::mesh_im& meshIm,
                   bgeot::multi_index tensorSize,
                   size_type filteredRegion)
    :im_(meshIm), nb_filtered_index_(0), nb_index_(0),
    filtered_region_(filteredRegion),
    locks_()
  {
    set_tensor_size(tensorSize);
    {local_guard lock = locks_.get_lock(); add_dependency(im_);}
    update_index_();
  }

  im_data::im_data(const getfem::mesh_im& meshIm, size_type filteredRegion)
    :im_(meshIm), nb_filtered_index_(0), nb_index_(0),
    filtered_region_(filteredRegion) {
    tensor_size_.resize(1);
    tensor_size_[0] = 1;
    nb_tensor_elem_ = 1;
    {local_guard lock = locks_.get_lock(); add_dependency(im_);}
    update_index_();
  }
    
  void im_data::update_index_() const {
    local_guard lock = locks_.get_lock();
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

  size_type im_data::nb_points_of_element(size_type cv) const{
    context_check();
    if (!im_.convex_index().is_in(cv)) return 0;
    return im_.int_method_of_element(cv)->approx_method()->nb_points_on_convex();
  }

  size_type im_data::index_of_point(size_type cv, size_type i) const{
    context_check();
    if (cv < int_point_index_.size()) return int_point_index_[cv] + i;
    else return size_type(-1);
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
    if(cv < filtered_int_point_index_.size())
    {
      if (filtered_int_point_index_[cv] == size_type(-1)) return size_type(-1);
      return filtered_int_point_index_[cv] + i;
    }
    else return size_type(-1);
  }

  dal::bit_vector im_data::filtered_convex_index() const{
    context_check();
    return filtered_convex_index_;
  }

  std::vector<size_type> im_data::filtered_index_of_first_point() const{
    context_check();
    return filtered_int_point_index_;
  }

  void im_data::set_region(size_type rg){
    filtered_region_ = rg;
    touch();
  }

  void im_data::update_from_context() const{
    update_index_();
    touch();
  }
}
