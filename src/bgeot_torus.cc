/*===========================================================================

 Copyright (C) 2014-2017 Liang Jin Lim

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
#include "getfem/bgeot_torus.h"

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
    torus_reference(bgeot::pconvex_ref ori_ref_convex) :
      convex_of_reference(torus_structure_descriptor(ori_ref_convex->structure()), false)
    {
      ori_ref_convex_ = ori_ref_convex;
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

    dal::pstatic_stored_object_key
      pk = std::make_shared<torus_structure_key>(ori_structure);
    dal::pstatic_stored_object o = dal::search_stored_object(pk);
    if (o) return std::dynamic_pointer_cast<const convex_structure>(o);

    auto p = std::make_shared<torus_structure>();
    pconvex_structure pcvs(p);
    p->Nc = dim_type(ori_structure->dim() + 1);
    p->nbpt = ori_structure->nb_points();
    p->nbf = ori_structure->nb_faces();

    p->faces_struct.resize(p->nbf);
    p->faces.resize(p->nbf);

    for (short_type j = 0; j < p->nbf; ++j){
      p->faces_struct[j] = ori_structure->faces_structure()[j];
      short_type nbIndex = ori_structure->nb_points_of_face(j);
      p->faces[j].resize(nbIndex);
      p->faces[j] = ori_structure->ind_points_of_face(j);
    }

    p->dir_points_.resize(ori_structure->ind_dir_points().size());
    p->dir_points_ = ori_structure->ind_dir_points();

    p->basic_pcvs = basic_structure(ori_structure);

    dal::add_stored_object(pk, pcvs, dal::PERMANENT_STATIC_OBJECT);
    return pcvs;
  }

  DAL_SIMPLE_KEY(torus_reference_key, pconvex_ref);

  pconvex_ref ptorus_reference(pconvex_ref ori_convex_reference)
  {
    dal::pstatic_stored_object_key
      pk = std::make_shared<torus_reference_key>(ori_convex_reference);
    dal::pstatic_stored_object o = dal::search_stored_object(pk);

    if (o) return std::dynamic_pointer_cast<const bgeot::convex_of_reference>(o);
    pconvex_ref p = std::make_shared<torus_reference>(ori_convex_reference);
    dal::add_stored_object(pk, p, p->structure(), p->pspt(),
			   dal::PERMANENT_STATIC_OBJECT);
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

    for (size_type i = 0; i < nb_points(); ++i){
      for (bgeot::dim_type n = 0; n < 2; ++n) pc(i, n) = pc2D(i, n);
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
      for (size_type i = 0; i < ind_ct.size(); ++i){
        for (bgeot::dim_type n = 0; n < 2; ++n) pc(i, n) = pc2D(i, n);
      }
  }

  void torus_geom_trans::compute_K_matrix
    (const bgeot::base_matrix &G, const bgeot::base_matrix &pc, bgeot::base_matrix &K) const{
      bgeot::geometric_trans::compute_K_matrix(G, pc, K);
      K(2, 2) = 0.0;
      for (short_type j = 0; j < nb_points(); ++ j) K(2, 2) += G(0, j) * pc(j, 2);
      for (short_type i = 0; i < 2; ++i) K(2, i) = K(i, 2) = 0;
  }

  void torus_geom_trans::poly_vector_hess(const base_node & /*pt*/,
                                          bgeot::base_matrix & /*pc*/) const{
    GMM_ASSERT1(false, "Sorry, Hessian is not supported in axisymmetric transformation.");
  }

  torus_geom_trans::torus_geom_trans(pgeometric_trans poriginal_trans) 
    : poriginal_trans_(poriginal_trans){
      geometric_trans::is_lin = poriginal_trans_->is_linear();
      geometric_trans::cvr = ptorus_reference(poriginal_trans_->convex_ref());
      complexity_ = poriginal_trans_->complexity();
      fill_standard_vertices();
  }

  pgeometric_trans torus_geom_trans::get_original_transformation() const{
    return poriginal_trans_;
  }

  bool is_torus_structure(pconvex_structure cvs){
    const torus_structure *cvs_torus = dynamic_cast<const torus_structure *>(cvs.get());
    return cvs_torus != NULL;
  }

  DAL_SIMPLE_KEY(torus_geom_trans_key, pgeometric_trans);

  pgeometric_trans torus_geom_trans_descriptor(pgeometric_trans poriginal_trans){
    dal::pstatic_stored_object_key
      pk = std::make_shared<torus_geom_trans_key>(poriginal_trans);
    dal::pstatic_stored_object o = dal::search_stored_object(pk);

    if (o) return std::dynamic_pointer_cast<const torus_geom_trans>(o);

    bgeot::pgeometric_trans p = std::make_shared<torus_geom_trans>(poriginal_trans);
    dal::add_stored_object(pk, p, dal::PERMANENT_STATIC_OBJECT);
    return p;
  }

  bool is_torus_geom_trans(pgeometric_trans pgt){
    const torus_geom_trans *pgt_torus = dynamic_cast<const torus_geom_trans *>(pgt.get());
    return pgt_torus != NULL;
  }

}