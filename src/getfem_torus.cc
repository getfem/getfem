/*===========================================================================

 Copyright (C) 2014-2020 Liang Jin Lim

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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

#include "getfem/bgeot_torus.h"
#include "getfem/getfem_torus.h"

namespace getfem
{

  void torus_fem::init(){
    cvr = poriginal_fem_->ref_convex(0);
    dim_ = cvr->structure()->dim();
    is_standard_fem = false;
    is_equiv = real_element_defined = true;
    is_pol = poriginal_fem_->is_polynomial();
    is_polycomp = poriginal_fem_->is_polynomialcomp();
    is_lag = poriginal_fem_->is_lagrange();
    es_degree = poriginal_fem_->estimated_degree();
    ntarget_dim = 3;

    std::stringstream nm;
    nm << "FEM_TORUS(" << poriginal_fem_->debug_name() << ")";
    debug_name_ = nm.str();

    init_cvs_node();
    GMM_ASSERT1(poriginal_fem_->target_dim() == 1, "Vectorial fems not supported");
    size_type nb_dof_origin = poriginal_fem_->nb_dof(0);
    for (size_type k = 0; k < nb_dof_origin; ++k)
    {
      for(size_type j = 0; j < 2; ++j){
        add_node(xfem_dof(poriginal_fem_->dof_types()[k], j), 
          poriginal_fem_->node_of_dof(0, k));
      }
    }
  }

  size_type torus_fem::index_of_global_dof
    (size_type /*cv*/, size_type i) const
  { return i; }

  void torus_fem::base_value(const base_node &, base_tensor &) const
  { GMM_ASSERT1(false, "No base values, real only element."); }

  void torus_fem::grad_base_value(const base_node &, base_tensor &) const
  { GMM_ASSERT1(false, "No grad values, real only element."); }

  void torus_fem::hess_base_value(const base_node &, base_tensor &) const{ 
    GMM_ASSERT1(false, "No hess values, real only element.");
  }

  void torus_fem::real_base_value(const fem_interpolation_context& c,
    base_tensor &t, bool) const{

    GMM_ASSERT1(!(poriginal_fem_->is_on_real_element()), "Original FEM must not be real.");

    base_tensor u_orig;
    poriginal_fem_->base_value(c.xref(), u_orig);
    if (!(poriginal_fem_->is_equivalent())){ 
      base_tensor u_temp = u_orig; 
      u_orig.mat_transp_reduction(u_temp, c.M(), 0); 
    }

    if(is_scalar_){
      t = u_orig;
      return;
    }

    //expand original base of [nb_base, 1] 
    //to vectorial form [nb_base * dim_, dim_ + 1]
    bgeot::multi_index tensor_size(u_orig.sizes());
    tensor_size[0] *= dim_;
    tensor_size[1] = ntarget_dim;
    t.adjust_sizes(tensor_size);
    for (size_type i = 0; i < u_orig.sizes()[0]; ++i) {
      for (dim_type j = 0; j < dim_; ++j) {
        t(i*dim_ + j, j) = u_orig(i, 0);
      }
    }
  }

  void torus_fem::real_grad_base_value
    (const getfem::fem_interpolation_context& c, base_tensor &t, bool) const 
  {
    GMM_ASSERT1(!(poriginal_fem_->is_on_real_element()), "Original FEM must not be real.");

    bgeot::scalar_type radius = std::abs(c.xreal()[0]);
    base_tensor u;
    bgeot::pstored_point_tab ppt = c.pgp()->get_ppoint_tab();
    getfem::pfem_precomp pfp = getfem::fem_precomp(poriginal_fem_, ppt, 0);
    base_tensor u_origin = pfp->grad(c.ii());
    //poriginal_fem_->grad_base_value(c.xref(), u_origin);
    GMM_ASSERT1(!u_origin.empty(), "Original FEM is unable to provide grad base value!");

    base_tensor n_origin = pfp->val(c.ii());
    //poriginal_fem_->base_value(c.xref(), n_origin);
    GMM_ASSERT1(!n_origin.empty(), "Original FEM is unable to provide base value!");

    //expand original grad of [nb_base, 1, dim_] 
    //to vectorial form [nb_base * dim_, dim_ + 1, dim_ + 1]
    const bgeot::multi_index &origin_size = u_origin.sizes();
    bgeot::multi_index tensor_size(origin_size);
    dim_type dim_size = is_scalar_ ? 1 : dim_;
    tensor_size[0] *= dim_size;
    tensor_size[1] = ntarget_dim;
    tensor_size[2] = dim_ + 1;
    u.adjust_sizes(tensor_size);
    for (size_type i = 0; i < origin_size[0]; ++i) { //dof
      for (dim_type j = 0; j < dim_size; ++j) {
        for (dim_type k = 0; k < dim_; ++k) {
          u(i*dim_size+j, j, k) = u_origin(i, 0, k);
        }
      }
    }
    t = u;
    t.mat_transp_reduction(u, c.B(), 2);

    if(is_scalar_) return;

    for (size_type i = 0; i < origin_size[0]; ++i) {
      t(i*dim_size, dim_, dim_) = n_origin[i] / radius;
    }
  }

  void torus_fem::real_hess_base_value
    (const getfem::fem_interpolation_context & /*c*/, base_tensor & /*t*/, bool) const 
  {
    GMM_ASSERT1(false, "Hessian not yet implemented in torus fem.");
  }

  void torus_fem::set_to_scalar(bool is_scalar){
    if(is_scalar_ == is_scalar) return;
    
    is_scalar_ = is_scalar;

    if(is_scalar_){
      ntarget_dim = 1;
      dof_types_.clear();
      init_cvs_node();
      size_type nb_dof_origin = poriginal_fem_->nb_dof(0);
      for (size_type k = 0; k < nb_dof_origin; ++k){
          add_node(poriginal_fem_->dof_types()[k], poriginal_fem_->node_of_dof(0, k));
      }
    }else{
      ntarget_dim = 3;
      dof_types_.clear();
      init_cvs_node();
      size_type nb_dof_origin = poriginal_fem_->nb_dof(0);
      for (size_type k = 0; k < nb_dof_origin; ++k)
      {
        for(size_type j = 0; j < 2; ++j){
          add_node(xfem_dof(poriginal_fem_->dof_types()[k], j), 
            poriginal_fem_->node_of_dof(0, k));
        }
      }
    }
  }

  pfem torus_fem::get_original_pfem() const{return poriginal_fem_;}

  DAL_SIMPLE_KEY(torus_fem_key, bgeot::size_type);

  getfem::pfem new_torus_fem(getfem::pfem pf) {
    static bgeot::size_type key_count = 0;
    ++key_count;
    getfem::pfem pfem_torus = std::make_shared<torus_fem>(pf);
    dal::pstatic_stored_object_key
      pk = std::make_shared<torus_fem_key>(key_count);
    dal::add_stored_object(pk, pfem_torus, pfem_torus->node_tab(0));
    return pfem_torus;
  }

  void del_torus_fem(getfem::pfem pf){
    const torus_fem *ptorus_fem = dynamic_cast<const torus_fem*>(pf.get());
    if (ptorus_fem != 0) dal::del_stored_object(pf);
  }

  void torus_mesh_fem::adapt_to_torus_(){

    for (dal::bv_visitor cv(linked_mesh().convex_index()); !cv.finished(); ++cv){
      pfem poriginal_fem = fem_of_element(cv);
      if(poriginal_fem == 0) continue;

      del_torus_fem(poriginal_fem);

      pfem pf = new_torus_fem(poriginal_fem);
      torus_fem *pf_torus = dynamic_cast<torus_fem*>(const_cast<virtual_fem*>(pf.get()));
      pf_torus->set_to_scalar((Qdim != 3));
      set_finite_element(cv, pf);
    }
    touch();
  }

  void torus_mesh_fem::enumerate_dof(void) const
  {
    const_cast<torus_mesh_fem*>(this)->adapt_to_torus_();

    for (dal::bv_visitor cv(linked_mesh().convex_index()); !cv.finished(); ++cv){
      pfem pf = fem_of_element(cv);
      if(pf == 0) continue;
      torus_fem *pf_torus = dynamic_cast<torus_fem*>(const_cast<virtual_fem*>(pf.get()));
      if(pf_torus == 0) continue;
      pf_torus->set_to_scalar((Qdim != 3));
    }

    mesh_fem::enumerate_dof();
  }

  torus_mesh::torus_mesh(std::string name) : mesh(std::move(name)){}

  scalar_type torus_mesh::convex_radius_estimate(size_type ic) const
  {
    base_matrix G;
    bgeot::vectors_to_base_matrix(G, points_of_convex(ic));
    G.resize(2, G.ncols());
    auto pgt_torus = std::dynamic_pointer_cast<const bgeot::torus_geom_trans>(trans_of_convex(ic));
    GMM_ASSERT2(pgt_torus, "Internal error, convex is not a torus transformation.");
    return getfem::convex_radius_estimate(pgt_torus->get_original_transformation(), G);
  }

  void torus_mesh::adapt(const getfem::mesh &original_mesh){
    clear();
    GMM_ASSERT1(original_mesh.dim() == 2, "Adapting torus feature must be a 2d mesh");
    mesh::copy_from(original_mesh);
    adapt();
  }

  void torus_mesh::adapt(){

    bgeot::node_tab node_tab_copy(pts);
    this->pts.clear();
    for(size_type pt = 0; pt < node_tab_copy.size(); ++pt){
      node_tab_copy[pt].resize(3);
      this->pts.add_node(node_tab_copy[pt]);
    }

    for(size_type i = 0; i < this->convex_tab.size(); ++i){
      bgeot::pconvex_structure pstructure 
        = bgeot::torus_structure_descriptor(convex_tab[i].cstruct);
      convex_tab[i].cstruct = pstructure;
    }

    for(size_type i = 0; i < this->gtab.size(); ++i){
      bgeot::pgeometric_trans pgeom_trans = bgeot::torus_geom_trans_descriptor(gtab[i]);
      gtab[i] = pgeom_trans;
    }
  }
}
