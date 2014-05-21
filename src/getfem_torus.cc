#include "getfem/getfem_torus.h"

namespace bgeot{

  /**torus_structure which extends a 2 dimensional structure with a radial dimension*/
  class torus_structure : public convex_structure{

    friend pconvex_structure torus_structure_descriptor(pconvex_structure);
  };

  class torus_reference : public convex_of_reference{

  public :
    scalar_type is_in(const base_node& point) const{
      GMM_ASSERT1(point.size() >= 2, "Invalid dimension of pt " << point);
      base_node point_2d = point;
      point_2d.resize(2);
      return ori_ref_convex_->is_in(point_2d);
    }

    scalar_type is_in_face(bgeot::short_type f, const base_node& point) const{
      GMM_ASSERT1(point.size() >= 2, "Invalid dimension of pt " << point);
      base_node point2D = point;
      point2D.resize(2);
      return ori_ref_convex_->is_in_face(f, point2D);
    }  
    torus_reference(bgeot::pconvex_ref ori_ref_convex){
      ori_ref_convex_ = ori_ref_convex;
      cvs = torus_structure_descriptor(ori_ref_convex->structure());
      convex<base_node>::points().resize(cvs->nb_points());
      normals_.resize(ori_ref_convex->normals().size());

      const std::vector<base_small_vector> &ori_normals = ori_ref_convex->normals();
      const stored_point_tab &ori_points = ori_ref_convex->points();
      for(size_type n = 0; n < ori_normals.size(); ++n){
       normals_[n] = ori_normals[n];
       normals_[n].resize(3);
      }            

      std::copy(ori_points.begin(), ori_points.end(), convex<base_node>::points().begin());
      for(size_type pt = 0; pt < convex<base_node>::points().size(); ++pt){
        convex<base_node>::points()[pt].resize(3);
      }
      ppoints = store_point_tab(convex<base_node>::points());
    }

  private:
    bgeot::pconvex_ref ori_ref_convex_;
  };

  DAL_SIMPLE_KEY(torus_structure_key, pconvex_structure);

  pconvex_structure torus_structure_descriptor(pconvex_structure ori_structure){

    dal::pstatic_stored_object o = dal::search_stored_object(torus_structure_key(ori_structure));
    if (o) return dal::stored_cast<convex_structure>(o);

    torus_structure *p = new torus_structure;
    p->Nc = ori_structure->dim() + 1;
    p->nbpt = ori_structure->nb_points();
    p->nbf = ori_structure->nb_faces();

    p->faces_struct.resize(p->nbf);
    p->faces.resize(p->nbf);

    for (size_type j = 0; j < p->nbf; ++j){
      p->faces_struct[j] = ori_structure->faces_structure()[j];
      short_type nbIndex = ori_structure->nb_points_of_face(j);
      p->faces[j].resize(nbIndex);
      p->faces[j] = ori_structure->ind_points_of_face(j);
    }

    p->dir_points_.resize(ori_structure->ind_dir_points().size());
    p->dir_points_ = ori_structure->ind_dir_points();

    p->basic_pcvs = ori_structure->basic_structure().get();

    torus_structure_key *key = new torus_structure_key(ori_structure);
    dal::add_stored_object(key, p, dal::PERMANENT_STATIC_OBJECT);
    return p;
  }

  DAL_SIMPLE_KEY(torus_reference_key, pconvex_ref);

  pconvex_ref ptorus_reference(pconvex_ref ori_convex_reference)
  {
    dal::pstatic_stored_object o 
      = dal::search_stored_object(torus_reference_key(ori_convex_reference));

    if (o) return dal::stored_cast<bgeot::convex_of_reference>(o);
    pconvex_ref p = new torus_reference(ori_convex_reference);
    dal::add_stored_object(new torus_reference_key(ori_convex_reference), p,
      p->structure(), &(p->points()), dal::PERMANENT_STATIC_OBJECT);
    return p;
  }

  void torus_geom_trans::poly_vector_val(const base_node &pt, bgeot::base_vector &val) const{
    base_node pt_2d(pt);
    pt_2d.resize(2);
    poriginal_trans_->poly_vector_val(pt_2d, val);
  }

  void torus_geom_trans::poly_vector_val(const base_node &pt, const bgeot::convex_ind_ct &ind_ct,
    bgeot::base_vector &val) const{
      base_node pt_2d(pt);
      pt_2d.resize(2);
      poriginal_trans_->poly_vector_val(pt_2d, ind_ct, val);     
  }

  void torus_geom_trans::poly_vector_grad(const base_node &pt, bgeot::base_matrix &pc) const{
    base_node pt2D(pt);
    pt2D.resize(2);
    bgeot::base_matrix pc2D(nb_points(), 2);
    poriginal_trans_->poly_vector_grad(pt2D, pc2D);

    bgeot::base_vector base_value;
    poriginal_trans_->poly_vector_val(pt2D, base_value);

    pc.resize(nb_points(), 3);
    for (size_type i = 0; i < nb_points(); ++i)
    {
      for (bgeot::dim_type n = 0; n < 2; ++n)
      {
        pc(i, n) = pc2D(i, n);
      }
      pc(i, 2) = base_value[i]; // radial direction, pc = base_x;
    }
  }

  void torus_geom_trans::poly_vector_grad(const base_node &pt,
    const bgeot::convex_ind_ct &ind_ct, bgeot::base_matrix &pc) const{
      base_node pt2D(pt);
      pt2D.resize(2);
      bgeot::base_matrix pc2D(ind_ct.size(), 2);
      poriginal_trans_->poly_vector_grad(pt2D, pc2D);
      pc.resize(ind_ct.size(), dim());
      for (size_type i = 0; i < ind_ct.size(); ++i)
      {
        for (bgeot::dim_type n = 0; n < 2; ++n)
        {
          pc(i, n) = pc2D(i, n);
        }
      }
  }

  void torus_geom_trans::compute_K_matrix
    (const bgeot::base_matrix &G, const bgeot::base_matrix &pc, bgeot::base_matrix &K) const{
      bgeot::geometric_trans::compute_K_matrix(G, pc, K);
      K(2, 2) = 0.0;
      for(int j = 0; j < nb_points(); ++ j) K(2, 2) += G(0, j) * pc(j, 2);
      for(int i = 0; i < 2; ++i) K(2, i) = K(i, 2) = 0;
  }

  void torus_geom_trans::poly_vector_hess(const base_node &pt, bgeot::base_matrix &pc) const{
    GMM_ASSERT1(false, "Sorry, Hessian is not supported in sxisymmetric transformation.");
  }

  torus_geom_trans::torus_geom_trans(pgeometric_trans poriginal_trans) 
    : poriginal_trans_(poriginal_trans){
      geometric_trans::is_lin = poriginal_trans_->is_linear();
      geometric_trans::cvr = ptorus_reference(poriginal_trans_->convex_ref());
      complexity_ = poriginal_trans_->complexity();
      fill_standard_vertices();
  }

  DAL_SIMPLE_KEY(torus_geom_trans_key, pgeometric_trans);

  pgeometric_trans torus_geom_trans_descriptor(pgeometric_trans poriginal_trans){
    dal::pstatic_stored_object o 
      = dal::search_stored_object(torus_geom_trans_key(poriginal_trans));

    if (o) return dal::stored_cast<torus_geom_trans>(o);

    bgeot::pgeometric_trans p = new torus_geom_trans(poriginal_trans);
    dal::add_stored_object(new torus_geom_trans_key(poriginal_trans), 
      p, dal::PERMANENT_STATIC_OBJECT);
    return p;
  }
}

namespace getfem
{

  void torus_fem::init()
  {
    cvr = poriginal_fem_->ref_convex(0);
    dim_ = cvr->structure()->dim();
    is_equiv = real_element_defined = true;
    is_polycomp = is_pol = is_lag = false;
    es_degree = 5;
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
    //expand original base of [nb_base, 1] 
    //to vectorial form [nb_base * dim_, dim_ + 1]
    bgeot::multi_index tensor_size(u_orig.sizes());
    tensor_size[0] *= dim_;
    tensor_size[1] = dim_ + 1;
    t.adjust_sizes(tensor_size);
    for(int i = 0; i < u_orig.sizes()[0]; ++i){
      for(int j = 0; j < dim_; ++j){
        t(i*dim_ + j, j) = u_orig(i, 0);
      }
    }
  }

  void torus_fem::real_grad_base_value
    (const getfem::fem_interpolation_context& c, base_tensor &t, bool) const 
  {
    GMM_ASSERT1(!(poriginal_fem_->is_on_real_element()), "Original FEM must not be real.");

    base_tensor u_origin;
    base_tensor u;
    poriginal_fem_->grad_base_value(c.xref(), u_origin);
    GMM_ASSERT1(!u_origin.empty(), "Original FEM is unable to provide grad base value!");

    base_tensor n_origin;
    poriginal_fem_->base_value(c.xref(), n_origin);
    GMM_ASSERT1(!n_origin.empty(), "Original FEM is unable to provide base value!");

    //expand original grad of [nb_base, 1, dim_] 
    //to vectorial form [nb_base * dim_, dim_ + 1, dim_ + 1]
    const bgeot::multi_index &origin_size = u_origin.sizes();
    bgeot::multi_index tensor_size(origin_size);
    tensor_size[0] *= dim_;
    tensor_size[1] = dim_ + 1;
    tensor_size[2] += 1;
    u.adjust_sizes(tensor_size);
    for(int i = 0; i < origin_size[0]; ++i){ //dof
      for(int j = 0; j < dim_; ++j){
        for(int k = 0; k < dim_; ++k){
          u(i*dim_+j, j, k) = u_origin(i, 0, k);
        }
      }
    }
    t = u;
    t.mat_transp_reduction(u, c.B(), 2);

    for(int i = 0; i < origin_size[0]; ++i){
      t(i*dim_, dim_, dim_) = n_origin[i] /c.xreal()[0];
    }
  }

  void torus_fem::real_hess_base_value
    (const getfem::fem_interpolation_context &c, base_tensor &t, bool) const 
  {
    GMM_ASSERT1(false, "Hessian not yet implemented in torus fem.");
  }
    
  DAL_SIMPLE_KEY(torus_fem_key, bgeot::size_type);

  getfem::pfem new_torus_fem(getfem::pfem pf){
    static bgeot::size_type key_count = 0;
    ++key_count;
    getfem::pfem pfem_torus= new torus_fem(pf);
    dal::add_stored_object(new torus_fem_key(key_count), pfem_torus);
    return pfem_torus;
  }

  void del_torus_fem(getfem::pfem pfem){    
    const torus_fem *ptorus_fem = dynamic_cast<const torus_fem*>(pfem.get());
    if (ptorus_fem != 0) dal::del_stored_object(pfem);
  }

  void torus_mesh_fem::adapt_to_torus(){

    for (dal::bv_visitor cv(linked_mesh().convex_index()); !cv.finished(); ++cv){
      pfem poriginal_fem = fem_of_element(cv);
      if(poriginal_fem == 0) continue;

      del_torus_fem(poriginal_fem);

      pfem pfem_torus = new_torus_fem(poriginal_fem);
      set_finite_element(cv, pfem_torus);
    }
    touch();
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
