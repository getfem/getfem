/*===========================================================================

 Copyright (C) 2013-2017 Yves Renard, Konstantinos Poulios and Andriy Andreykiv.

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

#include "getfem/getfem_generic_assembly.h"
#include "getfem/getfem_models.h"

namespace getfem {

// Structure describing a contact boundary (or contact body)
struct contact_boundary {
  size_type region;            // boundary region for the slave (source)
                               // and volume region for the master (target)
  const getfem::mesh_fem *mfu; // F.e.m. for the displacement.
  std::string dispname;        // Variable name for the displacement
  mutable const model_real_plain_vector *U;// Displacement
  mutable model_real_plain_vector U_unred; // Unreduced displacement

  contact_boundary(size_type r, const mesh_fem *mf, const std::string &dn)
    : region(r), mfu(mf), dispname(dn)
  {}
};

//extract element displacements from a contact boundary object
base_small_vector element_U(const contact_boundary &cb, size_type cv)
{
  auto U_elm = base_small_vector{};
  slice_vector_on_basic_dof_of_element(*(cb.mfu), *cb.U, cv, U_elm);
  return U_elm;
}

//Returns an iterator of a box which centre is closest to the given point
auto most_central_box(const bgeot::rtree::pbox_set &bset,
                      const bgeot::base_node       &pt) -> decltype(begin(bset))
{
  using namespace std;

  auto itmax = begin(bset);

  auto it = itmax;
  if (bset.size() > 1) {
    auto rate_max = scalar_type{-1};
    for (; it != end(bset); ++it) {
      auto rate_box = scalar_type{1};
      for (size_type i = 0; i < pt.size(); ++i) {
        auto h = (*it)->max[i] - (*it)->min[i];
        if (h > 0.) {
          auto rate = min((*it)->max[i] - pt[i], pt[i] - (*it)->min[i]) / h;
          rate_box = min(rate, rate_box);
        }
      }
      if (rate_box > rate_max) {
        itmax = it;
        rate_max = rate_box;
      }
    }
  }

  return itmax;
}

//Transformation that creates identity mapping between two contact boundaries,
//deformed with provided displacement fields
class  interpolate_transformation_on_deformed_domains
  : public virtual_interpolate_transformation {

  contact_boundary master;//also marked with a target or Y prefix/suffix
  contact_boundary slave; //also marked with a source or X prefix/suffix

  mutable bgeot::rtree element_boxes;
  mutable std::vector<size_type> box_to_convex; //index to obtain
                                                //a convex number from a box number
  mutable bgeot::geotrans_inv_convex gic;
  mutable fem_precomp_pool fppool;

  //Create a box tree based on the deformed elements of the master (target)
  void compute_element_boxes() const { // called by init
    base_matrix G;
    model_real_plain_vector Uelm; //element displacement
    element_boxes.clear();

    auto bnum = master.region;
    auto &mfu = *(master.mfu);
    auto &U   = *(master.U);
    auto &m   = mfu.linked_mesh();
    auto N    = m.dim();

    base_node Xdeformed(N), bmin(N), bmax(N);
    auto region = m.region(bnum);

    //the box tree creation and subsequent transformation inversion
    //should be done for all elements of the master, while integration
    //will be performed only on a thread partition of the slave
    region.prohibit_partitioning();

    GMM_ASSERT1(mfu.get_qdim() == N, "Wrong mesh_fem qdim");

    dal::bit_vector points_already_interpolated;
    std::vector<base_node> transformed_points(m.nb_max_points());
    box_to_convex.clear();
    box_to_convex.reserve(region.size());

    for (getfem::mr_visitor v(region, m); !v.finished(); ++v) {
      auto cv   = v.cv();
      auto pgt  = m.trans_of_convex(cv);
      auto pf_s = mfu.fem_of_element(cv);
      auto pfp  = fppool(pf_s, pgt->pgeometric_nodes());

      slice_vector_on_basic_dof_of_element(mfu, U, cv, Uelm);
      mfu.linked_mesh().points_of_convex(cv, G);

      auto ctx   = fem_interpolation_context{pgt, pfp, size_type(-1), G, cv};
      auto nb_pt = pgt->structure()->nb_points();

      for (size_type k = 0; k < nb_pt; ++k) {
        auto ind = m.ind_points_of_convex(cv)[k];

        // computation of a transformed vertex
        ctx.set_ii(k);
        if (points_already_interpolated.is_in(ind)) {
          Xdeformed = transformed_points[ind];
        } else {
          pf_s->interpolation(ctx, Uelm, Xdeformed, dim_type{N});
          Xdeformed += ctx.xreal(); //Xdeformed = U + Xo
          transformed_points[ind] = Xdeformed;
          points_already_interpolated.add(ind);
        }

        if (k == 0) // computation of bounding box
          bmin = bmax = Xdeformed;
        else {
          for (size_type l = 0; l < N; ++l) {
            bmin[l] = std::min(bmin[l], Xdeformed[l]);
            bmax[l] = std::max(bmax[l], Xdeformed[l]);
          }
        }
      }

      // Safety coefficient of 1.3 (for nonlinear transformations)
      // Expanding the box by 15% of the largest edge in every direction
      auto h = bmax[0] - bmin[0];
      for (size_type k = 1; k < N; ++k) h = std::max(h, bmax[k] - bmin[k]);
      for (size_type k = 0; k < N; ++k) {bmin[k] -= h * 0.15; bmax[k] += h * 0.15;}

      // Store the bounding box and additional information.
      element_boxes.add_box(bmin, bmax, box_to_convex.size());
      box_to_convex.push_back(cv);
    }
  }

  fem_interpolation_context deformed_master_context(size_type cv) const
  {
    auto &mfu  = *(master.mfu);
    auto G     = base_matrix{};
    auto pfu   = mfu.fem_of_element(cv);
    auto pgt   = master.mfu->linked_mesh().trans_of_convex(cv);
    auto pfp   = fppool(pfu, pgt->pgeometric_nodes());
    master.mfu->linked_mesh().points_of_convex(cv, G);
    return {pgt, pfp, size_type(-1), G, cv};
  }

  std::vector<bgeot::base_node> deformed_master_nodes(size_type cv) const {
    using namespace bgeot;
    using namespace std;

    auto nodes = vector<base_node>{};

    auto U_elm = element_U(master, cv);
    auto &mfu  = *(master.mfu);
    auto G     = base_matrix{};
    auto pfu   = mfu.fem_of_element(cv);
    auto pgt   = master.mfu->linked_mesh().trans_of_convex(cv);
    auto pfp   = fppool(pfu, pgt->pgeometric_nodes());
    auto N     = mfu.linked_mesh().dim();
    auto pt    = base_node(N);
    auto U     = base_small_vector(N);
    master.mfu->linked_mesh().points_of_convex(cv, G);
    auto ctx = fem_interpolation_context{pgt, pfp, size_type(-1), G, cv};
    auto nb_pt = pgt->structure()->nb_points();
    nodes.reserve(nb_pt);
    for (size_type k = 0; k < nb_pt; ++k) {
      ctx.set_ii(k);
      pfu->interpolation(ctx, U_elm, U, dim_type{N});
      gmm::add(ctx.xreal(), U, pt);
      nodes.push_back(pt);
    }

    return nodes;
  }

public:

  interpolate_transformation_on_deformed_domains(
    size_type              source_region,
    const getfem::mesh_fem &mf_source,
    const std::string      &source_displacements,
    size_type              target_region,
    const getfem::mesh_fem &mf_target,
    const std::string      &target_displacements)
    :
      slave{source_region, &mf_source, source_displacements},
      master{target_region, &mf_target, target_displacements}
{}


  void extract_variables(const ga_workspace           &workspace,
                         std::set<var_trans_pair>     &vars,
                         bool                         ignore_data,
                         const mesh                   &m_x,
                         const std::string            &interpolate_name) const override {
    if (!ignore_data || !(workspace.is_constant(master.dispname))){
      vars.emplace(master.dispname, interpolate_name);
      vars.emplace(slave.dispname, "");
    }
  }

  void init(const ga_workspace &workspace) const override {

    for (auto pcb : std::list<const contact_boundary*>{&master, &slave}) {
      auto &mfu = *(pcb->mfu);
      if (mfu.is_reduced()) {
        gmm::resize(pcb->U_unred, mfu.nb_basic_dof());
        mfu.extend_vector(workspace.value(pcb->dispname), pcb->U_unred);
        pcb->U = &(pcb->U_unred);
      } else {
        pcb->U = &(workspace.value(pcb->dispname));
      }
    }
    compute_element_boxes();
  };

  void finalize() const override {
    element_boxes.clear();
    box_to_convex.clear();
    master.U_unred.clear();
    slave.U_unred.clear();
    fppool.clear();
  }

  int transform(const ga_workspace                    &workspace,
                const mesh                            &m_x,
                fem_interpolation_context             &ctx_x,
                const base_small_vector               &/*Normal*/,
                const mesh                            **m_t,
                size_type                             &cv,
                short_type                            &face_num,
                base_node                             &P_ref,
                base_small_vector                     &N_y,
                std::map<var_trans_pair, base_tensor> &derivatives,
                bool                                  compute_derivatives) const override {

    auto &target_mesh = master.mfu->linked_mesh();
    *m_t = &target_mesh;
    auto transformation_success = false;

    using namespace gmm;
    using namespace bgeot;
    using namespace std;

    //compute a deformed point of the slave
    auto cv_x    = ctx_x.convex_num();
    auto U_elm_x = element_U(slave, cv_x);
    auto &mfu_x  = *(slave.mfu);
    auto pfu_x   = mfu_x.fem_of_element(cv_x);
    auto N       = mfu_x.linked_mesh().dim();
    auto U_x     = base_small_vector(N);
    auto G_x     = base_matrix{}; //coordinates of the source element nodes
    m_x.points_of_convex(cv_x, G_x);
    ctx_x.set_pf(pfu_x);
    pfu_x->interpolation(ctx_x, U_elm_x, U_x, dim_type{N});
    auto pt_x = base_small_vector(N); //deformed point of the slave
    add(ctx_x.xreal(), U_x, pt_x);

    //Find the best box from the master (target) that
    //corresponds to this point (The box which centre is the closest to the point).
    //Obtain the corresponding element number using the box id and box_to_convex
    //indices. Compute deformed nodes of the target element. Invert the geometric
    //transformation of the target element with deformed nodes, obtaining this way
    //reference coordinates of the target element
    auto bset = rtree::pbox_set{};
    element_boxes.find_boxes_at_point(pt_x, bset);
    while (!bset.empty())
    {
      auto itmax = most_central_box(bset, pt_x);
      cv = box_to_convex[(*itmax)->id];
      auto deformed_nodes_y = deformed_master_nodes(cv);
      gic.init(deformed_nodes_y, target_mesh.trans_of_convex(cv));
      auto converged = true;
      auto is_in = gic.invert(pt_x, P_ref, converged);
      if (is_in && converged) {
        face_num = static_cast<short_type>(-1);
        transformation_success = true;
        break;
      }
      if (bset.size() == 1) break;
      bset.erase(itmax);
    }

    //Since this transformation can be seen as Xsource + Usource - Utarget,
    //the corresponding stiffnesses are identity matrix for Usource and
    //minus identity for Utarget. The required answer in this function is
    //stiffness X shape function. Hence, returning shape function for Usource
    //and min shape function for Utarget
    if (compute_derivatives && transformation_success) {
      GMM_ASSERT2(derivatives.size() == 2,
                  "Expecting to return derivatives only for Umaster and Uslave");

      for (auto &pair : derivatives)
      {
        if (pair.first.varname == slave.dispname)
        {
          auto base_ux = base_tensor{};
          auto vbase_ux = base_matrix{} ;
          ctx_x.base_value(base_ux);
          auto qdim_ux = pfu_x->target_dim();
          auto ndof_ux = pfu_x->nb_dof(cv_x) * N / qdim_ux;
          vectorize_base_tensor(base_ux, vbase_ux, ndof_ux, qdim_ux, N);
          pair.second.adjust_sizes(ndof_ux, N);
          copy(vbase_ux.as_vector(), pair.second.as_vector());
        }
        else
        if (pair.first.varname == master.dispname)
        {
          auto ctx_y = deformed_master_context(cv);
          ctx_y.set_xref(P_ref);
          auto base_uy = base_tensor{};
          auto vbase_uy = base_matrix{} ;
          ctx_y.base_value(base_uy);
          auto pfu_y   = master.mfu->fem_of_element(cv);
          auto dim_y = master.mfu->linked_mesh().dim();
          auto qdim_uy = pfu_y->target_dim();
          auto ndof_uy = pfu_y->nb_dof(cv) * dim_y / qdim_uy;
          vectorize_base_tensor(base_uy, vbase_uy, ndof_uy, qdim_uy, dim_y);
          pair.second.adjust_sizes(ndof_uy, dim_y);
          copy(vbase_uy.as_vector(), pair.second.as_vector());
          scale(pair.second.as_vector(), -1.);
        }
        else GMM_ASSERT2(false, "unexpected derivative variable");
      }
    }

    return transformation_success ? 1 : 0;
  }

};

  void add_interpolate_transformation_on_deformed_domains
  (ga_workspace &workspace, const std::string &transname,
   const mesh &source_mesh, const std::string &source_displacements,
   const mesh_region &source_region, const mesh &target_mesh,
   const std::string &target_displacements, const mesh_region &target_region)
  {
    auto pmf_source = workspace.associated_mf(source_displacements);
    auto pmf_target = workspace.associated_mf(target_displacements);
    auto p_transformation
      = std::make_shared<interpolate_transformation_on_deformed_domains>(source_region.id(),
                                                                         *pmf_source,
                                                                         source_displacements,
                                                                         target_region.id(),
                                                                         *pmf_target,
                                                                         target_displacements);
    workspace.add_interpolate_transformation(transname, p_transformation);
  }

  void add_interpolate_transformation_on_deformed_domains
  (model &md, const std::string &transname,
   const mesh &source_mesh, const std::string &source_displacements,
   const mesh_region &source_region, const mesh &target_mesh,
   const std::string &target_displacements, const mesh_region &target_region)
  {
    auto &mf_source = md.mesh_fem_of_variable(source_displacements);
    auto mf_target = md.mesh_fem_of_variable(target_displacements);
    auto p_transformation
      = std::make_shared<interpolate_transformation_on_deformed_domains>(source_region.id(),
                                                                         mf_source,
                                                                         source_displacements,
                                                                         target_region.id(),
                                                                         mf_target,
                                                                         target_displacements);
    md.add_interpolate_transformation(transname, p_transformation);
  }

}  /* end of namespace getfem.                                             */