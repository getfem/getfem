/*===========================================================================

 Copyright (C) 2013-2018 Yves Renard

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

#include "getfem/getfem_generic_assembly_tree.h"
#include "getfem/getfem_generic_assembly_semantic.h"
#include "getfem/getfem_generic_assembly_compile_and_exec.h"
#include "getfem/getfem_generic_assembly_functions_and_operators.h"

namespace getfem {

  //=========================================================================
  // Interpolation functions
  //=========================================================================

  // general Interpolation
  void ga_interpolation(ga_workspace &workspace,
                        ga_interpolation_context &gic) {
    ga_instruction_set gis;
    ga_compile_interpolation(workspace, gis);
    ga_interpolation_exec(gis, workspace, gic);
  }

  // Interpolation on a Lagrange fem on the same mesh
  struct ga_interpolation_context_fem_same_mesh
    : public ga_interpolation_context {
    base_vector &result;
    std::vector<int> dof_count;
    const mesh_fem &mf;
    bool initialized;
    bool is_torus;
    size_type s;

    virtual bgeot::pstored_point_tab
    ppoints_for_element(size_type cv, short_type f,
                       std::vector<size_type> &ind) const {
      pfem pf = mf.fem_of_element(cv);
      GMM_ASSERT1(pf->is_lagrange(),
                  "Only Lagrange fems can be used in interpolation");

      if (f != short_type(-1)) {

        for (size_type i = 0;
             i < pf->node_convex(cv).structure()->nb_points_of_face(f); ++i)
          ind.push_back
            (pf->node_convex(cv).structure()->ind_points_of_face(f)[i]);
      } else {
        for (size_type i = 0; i < pf->node_convex(cv).nb_points(); ++i)
          ind.push_back(i);
      }

      return pf->node_tab(cv);
    }

    virtual bool use_pgp(size_type) const { return true; }
    virtual bool use_mim() const { return false; }

    void init_(size_type si, size_type q, size_type qmult) {
      s = si;
      gmm::resize(result, qmult * mf.nb_basic_dof());
      gmm::clear(result);
      gmm::resize(dof_count, mf.nb_basic_dof()/q);
      gmm::clear(dof_count);
      initialized = true;
    }

    void store_result_for_torus(size_type cv, size_type i, base_tensor &t) {
      size_type target_dim = mf.fem_of_element(cv)->target_dim();
      GMM_ASSERT2(target_dim == 3, "Invalid torus fem.");
      size_type qdim = 1;
      size_type result_dim = 2;
      if (!initialized) {init_(qdim, qdim, qdim);}
      size_type idof = mf.ind_basic_dof_of_element(cv)[i];
      result[idof] = t[idof%result_dim];
      ++dof_count[idof];
    }

    virtual void store_result(size_type cv, size_type i, base_tensor &t) {
      if (is_torus){
        store_result_for_torus(cv, i, t);
        return;
      }
      size_type si = t.size();
      size_type q = mf.get_qdim();
      size_type qmult = si / q;
      GMM_ASSERT1( (si % q) == 0, "Incompatibility between the mesh_fem and "
                   "the size of the expression to be interpolated");
      if (!initialized) { init_(si, q, qmult); }
      GMM_ASSERT1(s == si, "Internal error");
      size_type idof = mf.ind_basic_dof_of_element(cv)[i*q];
      gmm::add(t.as_vector(),
               gmm::sub_vector(result, gmm::sub_interval(qmult*idof, s)));
      (dof_count[idof/q])++;
    }

    virtual void finalize() {
      std::vector<size_type> data(3);
      data[0] = initialized ? result.size() : 0;
      data[1] = initialized ? dof_count.size() : 0;
      data[2] = initialized ? s : 0;
      MPI_MAX_VECTOR(data);
      if (!initialized) {
        if (data[0]) {
          gmm::resize(result, data[0]);
          gmm::resize(dof_count, data[1]);
          gmm::clear(dof_count);
          s = data[2];
        }
        gmm::clear(result);
      }
      if (initialized)
        GMM_ASSERT1(gmm::vect_size(result) == data[0] &&  s == data[2] &&
                  gmm::vect_size(dof_count) == data[1], "Incompatible sizes");
      MPI_SUM_VECTOR(result);
      MPI_SUM_VECTOR(dof_count);
      for (size_type i = 0; i < dof_count.size(); ++i)
        if (dof_count[i])
          gmm::scale(gmm::sub_vector(result, gmm::sub_interval(s*i, s)),
                     scalar_type(1) / scalar_type(dof_count[i]));
    }

    virtual const mesh &linked_mesh() { return mf.linked_mesh(); }

    ga_interpolation_context_fem_same_mesh(const mesh_fem &mf_, base_vector &r)
      : result(r), mf(mf_), initialized(false), is_torus(false) {
      GMM_ASSERT1(!(mf.is_reduced()),
                  "Interpolation on reduced fem is not allowed");
      if (dynamic_cast<const torus_mesh_fem*>(&mf)){
        auto first_cv = mf.first_convex_of_basic_dof(0);
        auto target_dim = mf.fem_of_element(first_cv)->target_dim();
        if (target_dim == 3) is_torus = true;
      }
    }
  };

  void ga_interpolation_Lagrange_fem
  (ga_workspace &workspace, const mesh_fem &mf, base_vector &result) {
    ga_interpolation_context_fem_same_mesh gic(mf, result);
    ga_interpolation(workspace, gic);
  }

  void ga_interpolation_Lagrange_fem
  (const getfem::model &md, const std::string &expr, const mesh_fem &mf,
   base_vector &result, const mesh_region &rg) {
    ga_workspace workspace(md);
    workspace.add_interpolation_expression(expr, mf.linked_mesh(), rg);
    ga_interpolation_Lagrange_fem(workspace, mf, result);
  }

  // Interpolation on a cloud of points
  struct ga_interpolation_context_mti
    : public ga_interpolation_context {
    base_vector &result;
    const mesh_trans_inv &mti;
    bool initialized;
    size_type s, nbdof;


    virtual bgeot::pstored_point_tab
    ppoints_for_element(size_type cv, short_type,
                       std::vector<size_type> &ind) const {
      std::vector<size_type> itab;
      mti.points_on_convex(cv, itab);
      std::vector<base_node> pt_tab(itab.size());
      for (size_type i = 0; i < itab.size(); ++i) {
        pt_tab[i] = mti.reference_coords()[itab[i]];
        ind.push_back(i);
      }
      return store_point_tab(pt_tab);
    }

    virtual bool use_pgp(size_type) const { return false; }
    virtual bool use_mim() const { return false; }

    virtual void store_result(size_type cv, size_type i, base_tensor &t) {
      size_type si = t.size();
      if (!initialized) {
        s = si;
        gmm::resize(result, s * nbdof);
        gmm::clear(result);
        initialized = true;
      }
      GMM_ASSERT1(s == si, "Internal error");
      size_type ipt = mti.point_on_convex(cv, i);
      size_type dof_t = mti.id_of_point(ipt);
      gmm::add(t.as_vector(),
               gmm::sub_vector(result, gmm::sub_interval(s*dof_t, s)));
    }

    virtual void finalize() {
      std::vector<size_type> data(2);
      data[0] = initialized ? result.size() : 0;
      data[1] = initialized ? s : 0;
      MPI_MAX_VECTOR(data);
      if (!initialized) {
        if (data[0]) {
          gmm::resize(result, data[0]);
          s = data[1];
        }
        gmm::clear(result);
      }
      if (initialized)
        GMM_ASSERT1(gmm::vect_size(result) == data[0] &&  s == data[1],
                    "Incompatible sizes");
      MPI_SUM_VECTOR(result);
    }

    virtual const mesh &linked_mesh() { return mti.linked_mesh(); }

    ga_interpolation_context_mti(const mesh_trans_inv &mti_, base_vector &r,
                                 size_type nbdof_ = size_type(-1))
      : result(r), mti(mti_), initialized(false), nbdof(nbdof_) {
      if (nbdof == size_type(-1)) nbdof = mti.nb_points();
    }
  };

  // Distribute to be parallelized
  void ga_interpolation_mti
  (const getfem::model &md, const std::string &expr, mesh_trans_inv &mti,
   base_vector &result, const mesh_region &rg, int extrapolation,
   const mesh_region &rg_source, size_type nbdof) {

    ga_workspace workspace(md);
    workspace.add_interpolation_expression(expr, mti.linked_mesh(), rg);

    mti.distribute(extrapolation, rg_source);
    ga_interpolation_context_mti gic(mti, result, nbdof);
    ga_interpolation(workspace, gic);
  }


  // Interpolation on a im_data
  struct ga_interpolation_context_im_data
    : public ga_interpolation_context {
    base_vector &result;
    const im_data &imd;
    bool initialized;
    size_type s;

    virtual bgeot::pstored_point_tab
    ppoints_for_element(size_type cv, short_type f,
                        std::vector<size_type> &ind) const {
      pintegration_method pim =imd.linked_mesh_im().int_method_of_element(cv);
      if (pim->type() == IM_NONE) return bgeot::pstored_point_tab();
      GMM_ASSERT1(pim->type() == IM_APPROX, "Sorry, exact methods cannot "
                  "be used in high level generic assembly");
      size_type i_start(0), i_end(0);
      if (f == short_type(-1))
        i_end = pim->approx_method()->nb_points_on_convex();
      else {
        i_start = pim->approx_method()->ind_first_point_on_face(f);
        i_end = i_start + pim->approx_method()->nb_points_on_face(f);
      }
      for (size_type i = i_start; i < i_end; ++i) ind.push_back(i);
      return pim->approx_method()->pintegration_points();
    }

    virtual bool use_pgp(size_type cv) const {
      pintegration_method pim =imd.linked_mesh_im().int_method_of_element(cv);
      if (pim->type() == IM_NONE) return false;
      GMM_ASSERT1(pim->type() == IM_APPROX, "Sorry, exact methods cannot "
                  "be used in high level generic assembly");
      return !(pim->approx_method()->is_built_on_the_fly());
    }
    virtual bool use_mim() const { return true; }

    virtual void store_result(size_type cv, size_type i, base_tensor &t) {
      size_type si = t.size();
      if (!initialized) {
        s = si;
        GMM_ASSERT1(imd.tensor_size() == t.sizes() ||
                    (imd.tensor_size().size() == size_type(1) &&
                     imd.tensor_size()[0] == size_type(1) &&
                     si == size_type(1)),
                    "Im_data tensor size " << imd.tensor_size() <<
                    " does not match the size of the interpolated "
                    "expression " << t.sizes() << ".");
        gmm::resize(result, s * imd.nb_filtered_index());
        gmm::clear(result);
        initialized = true;
      }
      GMM_ASSERT1(s == si, "Internal error");
      size_type ipt = imd.filtered_index_of_point(cv, i);
      GMM_ASSERT1(ipt != size_type(-1),
                  "Im data with no data on the current integration point.");
      gmm::add(t.as_vector(),
               gmm::sub_vector(result, gmm::sub_interval(s*ipt, s)));
    }

    virtual void finalize() {
      std::vector<size_type> data(2);
      data[0] = initialized ? result.size() : 0;
      data[1] = initialized ? s : 0;
      MPI_MAX_VECTOR(data);
      if (initialized) {
        GMM_ASSERT1(gmm::vect_size(result) == data[0] &&  s == data[1],
                    "Incompatible sizes");
      } else {
        if (data[0]) {
          gmm::resize(result, data[0]);
          s = data[1];
        }
        gmm::clear(result);
      }
      MPI_SUM_VECTOR(result);
    }

    virtual const mesh &linked_mesh() { return imd.linked_mesh(); }

    ga_interpolation_context_im_data(const im_data &imd_, base_vector &r)
      : result(r), imd(imd_), initialized(false) { }
  };

  void ga_interpolation_im_data
  (ga_workspace &workspace, const im_data &imd, base_vector &result) {
    ga_interpolation_context_im_data gic(imd, result);
    ga_interpolation(workspace, gic);
  }

  void ga_interpolation_im_data
  (const getfem::model &md, const std::string &expr, const im_data &imd,
   base_vector &result, const mesh_region &rg) {
    ga_workspace workspace(md);
    workspace.add_interpolation_expression
      (expr, imd.linked_mesh_im(), rg);

    ga_interpolation_im_data(workspace, imd, result);
  }


  // Interpolation on a stored_mesh_slice
  struct ga_interpolation_context_mesh_slice
    : public ga_interpolation_context {
    base_vector &result;
    const stored_mesh_slice &sl;
    bool initialized;
    size_type s;
    std::vector<size_type> first_node;

    virtual bgeot::pstored_point_tab
    ppoints_for_element(size_type cv, short_type f,
                        std::vector<size_type> &ind) const {
      GMM_ASSERT1(f == short_type(-1), "No support for interpolation on faces"
                                       " for a stored_mesh_slice yet.");
      size_type ic = sl.convex_pos(cv);
      const mesh_slicer::cs_nodes_ct &nodes = sl.nodes(ic);
      std::vector<base_node> pt_tab(nodes.size());
      for (size_type i=0; i < nodes.size(); ++i) {
        pt_tab[i] = nodes[i].pt_ref;
        ind.push_back(i);
      }
      return store_point_tab(pt_tab);
    }

    virtual bool use_pgp(size_type /* cv */) const { return false; } // why not?
    virtual bool use_mim() const { return false; }

    virtual void store_result(size_type cv, size_type i, base_tensor &t) {
      size_type si = t.size();
      if (!initialized) {
        s = si;
        gmm::resize(result, s * sl.nb_points());
        gmm::clear(result);
        initialized = true;
        first_node.resize(sl.nb_convex());
        for (size_type ic=0; ic < sl.nb_convex()-1; ++ic)
          first_node[ic+1] = first_node[ic] + sl.nodes(ic).size();
      }
      GMM_ASSERT1(s == si && result.size() == s * sl.nb_points(), "Internal error");
      size_type ic = sl.convex_pos(cv);
      size_type ipt = first_node[ic] + i;
      gmm::add(t.as_vector(),
               gmm::sub_vector(result, gmm::sub_interval(s*ipt, s)));
    }

    virtual void finalize() {
      std::vector<size_type> data(2);
      data[0] = initialized ? result.size() : 0;
      data[1] = initialized ? s : 0;
      MPI_MAX_VECTOR(data);
      if (initialized) {
        GMM_ASSERT1(gmm::vect_size(result) == data[0] &&  s == data[1],
                    "Incompatible sizes");
      } else {
        if (data[0]) {
          gmm::resize(result, data[0]);
          s = data[1];
        }
        gmm::clear(result);
      }
      MPI_SUM_VECTOR(result);
    }

    virtual const mesh &linked_mesh() { return sl.linked_mesh(); }

    ga_interpolation_context_mesh_slice(const stored_mesh_slice &sl_, base_vector &r)
      : result(r), sl(sl_), initialized(false) { }
  };

  void ga_interpolation_mesh_slice
  (ga_workspace &workspace, const stored_mesh_slice &sl, base_vector &result) {
    ga_interpolation_context_mesh_slice gic(sl, result);
    ga_interpolation(workspace, gic);
  }

  void ga_interpolation_mesh_slice
  (const getfem::model &md, const std::string &expr, const stored_mesh_slice &sl,
   base_vector &result, const mesh_region &rg) {
    ga_workspace workspace(md);
    workspace.add_interpolation_expression(expr, sl.linked_mesh(), rg);
    ga_interpolation_mesh_slice(workspace, sl, result);
  }


  //=========================================================================
  // Local projection functions
  //=========================================================================

  void ga_local_projection(const getfem::model &md, const mesh_im &mim,
                           const std::string &expr, const mesh_fem &mf,
                           base_vector &result, const mesh_region &region) {

    // The case where the expression is a vector one and mf a scalar fem is
    // not taken into account for the moment.

    // Could be improved by not performing the assembly of the global mass matrix
    // working locally. This means a specific assembly.
    model_real_sparse_matrix M(mf.nb_dof(), mf.nb_dof());
    asm_mass_matrix(M, mim, mf, region);

    ga_workspace workspace(md);
    size_type nbdof = md.nb_dof();
    gmm::sub_interval I(nbdof, mf.nb_dof());
    workspace.add_fem_variable("c__dummy_var_95_", mf, I, base_vector(nbdof));
    if (mf.get_qdims().size() > 1)
      workspace.add_expression("("+expr+"):Test_c__dummy_var_95_",mim,region,2);
    else
      workspace.add_expression("("+expr+").Test_c__dummy_var_95_",mim,region,2);
    base_vector residual(nbdof+mf.nb_dof());
    workspace.set_assembled_vector(residual);
    workspace.assembly(1);
    getfem::base_vector F(mf.nb_dof());
    gmm::resize(result, mf.nb_dof());
    gmm::copy(gmm::sub_vector(residual, I), F);

    getfem::base_matrix loc_M;
    getfem::base_vector loc_U;
    for (mr_visitor v(region, mf.linked_mesh(), true); !v.finished(); ++v) {
      if (mf.convex_index().is_in(v.cv())) {
        size_type nd = mf.nb_basic_dof_of_element(v.cv());
        loc_M.base_resize(nd, nd); gmm::resize(loc_U, nd);
        gmm::sub_index J(mf.ind_basic_dof_of_element(v.cv()));
        gmm::copy(gmm::sub_matrix(M, J, J), loc_M);
        gmm::lu_solve(loc_M, loc_U, gmm::sub_vector(F, J));
        gmm::copy(loc_U, gmm::sub_vector(result, J));
      }
    }
    MPI_SUM_VECTOR(result);
  }

  //=========================================================================
  // Interpolate transformation with an expression
  //=========================================================================

  class interpolate_transformation_expression
    : public virtual_interpolate_transformation, public context_dependencies {

    struct workspace_gis_pair : public std::pair<ga_workspace, ga_instruction_set> {
      inline ga_workspace &workspace() { return this->first; }
      inline ga_instruction_set &gis() { return this->second; }
    };

    const mesh &source_mesh;
    const mesh &target_mesh;
    const size_type target_region;
    std::string expr;
    mutable bgeot::rtree element_boxes;
    mutable bool recompute_elt_boxes;
    mutable ga_workspace local_workspace;
    mutable ga_instruction_set local_gis;
    mutable bgeot::geotrans_inv_convex gic;
    mutable base_node P;
    mutable std::set<var_trans_pair> used_vars;
    mutable std::set<var_trans_pair> used_data;
    mutable std::map<var_trans_pair,
                     workspace_gis_pair> compiled_derivatives;
    mutable bool extract_variable_done;
    mutable bool extract_data_done;

  public:
    void update_from_context() const {
      recompute_elt_boxes = true;
    }

    void extract_variables(const ga_workspace &workspace,
                           std::set<var_trans_pair> &vars,
                           bool ignore_data, const mesh &/* m */,
                           const std::string &/* interpolate_name */) const {
      if ((ignore_data && !extract_variable_done) ||
          (!ignore_data && !extract_data_done)) {
        if (ignore_data)
          used_vars.clear();
        else
          used_data.clear();
        ga_workspace aux_workspace;
        aux_workspace = ga_workspace(true, workspace);
        aux_workspace.clear_expressions();
        aux_workspace.add_interpolation_expression(expr, source_mesh);
        for (size_type i = 0; i < aux_workspace.nb_trees(); ++i)
          ga_extract_variables(aux_workspace.tree_info(i).ptree->root,
                               aux_workspace, source_mesh,
                               ignore_data ? used_vars : used_data,
                               ignore_data);
        if (ignore_data)
          extract_variable_done = true;
        else
          extract_data_done = true;
      }
      if (ignore_data)
        vars.insert(used_vars.begin(), used_vars.end());
      else
        vars.insert(used_data.begin(), used_data.end());
    }

    void init(const ga_workspace &workspace) const {
      size_type N = target_mesh.dim();

      // Expression compilation
      local_workspace = ga_workspace(true, workspace);
      local_workspace.clear_expressions();

      local_workspace.add_interpolation_expression(expr, source_mesh);
      local_gis = ga_instruction_set();
      ga_compile_interpolation(local_workspace, local_gis);

      // In fact, transformations are not allowed  ... for future compatibility
      for (const std::string &transname : local_gis.transformations)
        local_workspace.interpolate_transformation(transname)
          ->init(local_workspace);

      if (!extract_variable_done) {
        std::set<var_trans_pair> vars;
        extract_variables(workspace, vars, true, source_mesh, "");
      }

      for (const var_trans_pair &var : used_vars) {
        workspace_gis_pair &pwi = compiled_derivatives[var];
        pwi.workspace() = local_workspace;
        pwi.gis() = ga_instruction_set();
        if (pwi.workspace().nb_trees()) {
          ga_tree tree = *(pwi.workspace().tree_info(0).ptree);
          ga_derivative(tree, pwi.workspace(), source_mesh,
                        var.varname, var.transname, 1);
          if (tree.root)
            ga_semantic_analysis(tree, local_workspace, dummy_mesh(),
				 1, false, true);
          ga_compile_interpolation(pwi.workspace(), pwi.gis());
        }
      }

      // Element_boxes update (if necessary)
      if (recompute_elt_boxes) {

        element_boxes.clear();
        base_node bmin(N), bmax(N);
        const dal::bit_vector&
          convex_index = (target_region == mesh_region::all_convexes().id())
                       ? target_mesh.convex_index()
                       : target_mesh.region(target_region).index();
        for (dal::bv_visitor cv(convex_index); !cv.finished(); ++cv) {

          bgeot::pgeometric_trans pgt = target_mesh.trans_of_convex(cv);

          size_type nbd_t = pgt->nb_points();
          if (nbd_t) {
            gmm::copy(target_mesh.points_of_convex(cv)[0], bmin);
            gmm::copy(bmin, bmax);
          } else {
            gmm::clear(bmin);
            gmm::clear(bmax);
          }
          for (short_type ip = 1; ip < nbd_t; ++ip) {
            // size_type ind = target_mesh.ind_points_of_convex(cv)[ip];
            const base_node &pt = target_mesh.points_of_convex(cv)[ip];

            for (size_type k = 0; k < N; ++k) {
              bmin[k] = std::min(bmin[k], pt[k]);
              bmax[k] = std::max(bmax[k], pt[k]);
            }
          }

          scalar_type h = bmax[0] - bmin[0];
          for (size_type k = 1; k < N; ++k) h = std::max(h, bmax[k]-bmin[k]);
          if (pgt->is_linear()) h *= 1E-4;
          for (auto &&val : bmin) val -= h*0.2;
          for (auto &&val : bmax) val += h*0.2;

          element_boxes.add_box(bmin, bmax, cv);
        }
        recompute_elt_boxes = false;
      }
    }

    void finalize() const {
      for (const std::string &transname : local_gis.transformations)
        local_workspace.interpolate_transformation(transname)->finalize();
      local_gis = ga_instruction_set();
    }

    std::string expression() const { return expr; }

    int transform(const ga_workspace &/*workspace*/, const mesh &m,
                  fem_interpolation_context &ctx_x,
                  const base_small_vector &Normal,
                  const mesh **m_t,
                  size_type &cv, short_type &face_num,
                  base_node &P_ref,
                  base_small_vector &/*N_y*/,
                  std::map<var_trans_pair, base_tensor> &derivatives,
                  bool compute_derivatives) const {
      int ret_type = 0;

      ga_interpolation_single_point_exec(local_gis, local_workspace, ctx_x,
                                         Normal, m);

      GMM_ASSERT1(local_workspace.assembled_tensor().size() == m.dim(),
                  "Wrong dimension of the transformation expression");
      P.resize(m.dim());
      gmm::copy(local_workspace.assembled_tensor().as_vector(), P);

      bgeot::rtree::pbox_set bset;
      element_boxes.find_boxes_at_point(P, bset);
      *m_t = &target_mesh;

      while (bset.size()) {
        bgeot::rtree::pbox_set::iterator it = bset.begin(), itmax = it;

        if (bset.size() > 1) {
          // Searching the box for which the point is the most in the interior
          scalar_type rate_max = scalar_type(-1);
          for (; it != bset.end(); ++it) {

            scalar_type rate_box = scalar_type(1);
            for (size_type i = 0; i < m.dim(); ++i) {
              scalar_type h = (*it)->max[i] - (*it)->min[i];
              if (h > scalar_type(0)) {
                scalar_type rate
                  = std::min((*it)->max[i] - P[i], P[i] - (*it)->min[i]) / h;
                rate_box = std::min(rate, rate_box);
              }
            }
            if (rate_box > rate_max) {
              itmax = it;
              rate_max = rate_box;
            }
          }
        }

        cv = (*itmax)->id;
        gic.init(target_mesh.points_of_convex(cv),
                 target_mesh.trans_of_convex(cv));

        bool converged = true;
        bool is_in = gic.invert(P, P_ref, converged, 1E-4);
        // cout << "cv = " << cv << " P = " << P << " P_ref = " << P_ref << endl;
        // cout << " is_in = " << int(is_in) << endl;
        // for (size_type iii = 0;
        //     iii < target_mesh.points_of_convex(cv).size(); ++iii)
        //  cout << target_mesh.points_of_convex(cv)[iii] << endl;

        if (is_in && converged) {
          face_num = short_type(-1); // Should detect potential faces ?
          ret_type = 1;
          break;
        }

        if (bset.size() == 1) break;
        bset.erase(itmax);
      }

      // Note on derivatives of the transformation : for efficiency and
      // simplicity reasons, the derivative should be computed with
      // the value of corresponding test functions. This means that
      // for a transformation F(u) the computed derivative is F'(u).Test_u
      // including the Test_u.
      if (compute_derivatives) { // Could be optimized ?
        for (auto &&d : derivatives) {
          workspace_gis_pair &pwi = compiled_derivatives[d.first];

          gmm::clear(pwi.workspace().assembled_tensor().as_vector());
          ga_function_exec(pwi.gis());
          d.second = pwi.workspace().assembled_tensor();
        }
      }
      return ret_type;
    }

    interpolate_transformation_expression
    (const mesh &sm, const mesh &tm, size_type trg, const std::string &expr_)
      : source_mesh(sm), target_mesh(tm), target_region(trg), expr(expr_),
        recompute_elt_boxes(true), extract_variable_done(false),
        extract_data_done(false)
    { this->add_dependency(tm); }

  };


  void add_interpolate_transformation_from_expression
  (ga_workspace &workspace, const std::string &name, const mesh &sm,
   const mesh &tm, const std::string &expr) {
    add_interpolate_transformation_from_expression
    (workspace, name, sm, tm, size_type(-1), expr);
  }

  void add_interpolate_transformation_from_expression
  (ga_workspace &workspace, const std::string &name, const mesh &sm,
   const mesh &tm, size_type trg, const std::string &expr) {
    pinterpolate_transformation
      p = std::make_shared<interpolate_transformation_expression>
          (sm, tm, trg, expr);
    workspace.add_interpolate_transformation(name, p);
  }

  void add_interpolate_transformation_from_expression
  (model &md, const std::string &name, const mesh &sm, const mesh &tm,
   const std::string &expr) {
    add_interpolate_transformation_from_expression
    (md, name, sm, tm, size_type(-1), expr);
  }

  void add_interpolate_transformation_from_expression
  (model &md, const std::string &name, const mesh &sm, const mesh &tm,
   size_type trg, const std::string &expr) {
    pinterpolate_transformation
      p = std::make_shared<interpolate_transformation_expression>
          (sm, tm, trg, expr);
    md.add_interpolate_transformation(name, p);
  }

  //=========================================================================
  // Interpolate transformation on neighbour element (for internal faces)
  //=========================================================================

  class interpolate_transformation_neighbour
    : public virtual_interpolate_transformation, public context_dependencies {

  public:
    void update_from_context() const {}
    void extract_variables(const ga_workspace &/* workspace */,
                           std::set<var_trans_pair> &/* vars */,
                           bool /* ignore_data */, const mesh &/* m */,
                           const std::string &/* interpolate_name */) const {}
    void init(const ga_workspace &/* workspace */) const {}
    void finalize() const {}
    
    std::string expression(void) const { return "X"; }

    int transform(const ga_workspace &/*workspace*/, const mesh &m_x,
                  fem_interpolation_context &ctx_x,
                  const base_small_vector &/*Normal*/, const mesh **m_t,
                  size_type &cv, short_type &face_num,
                  base_node &P_ref,
                  base_small_vector &/*N_y*/,
                  std::map<var_trans_pair, base_tensor> &/*derivatives*/,
                  bool compute_derivatives) const {

      int ret_type = 0;
      *m_t = &m_x;
      size_type cv_x = ctx_x.convex_num();
      short_type face_x = ctx_x.face_num();
      GMM_ASSERT1(face_x != short_type(-1), "Neighbour transformation can "
                  "only be applied to internal faces");

      auto adj_face = m_x.adjacent_face(cv_x, face_x);

      if (adj_face.cv != size_type(-1)) {
        bgeot::geotrans_inv_convex gic;
        gic.init(m_x.points_of_convex(adj_face.cv),
                 m_x.trans_of_convex(adj_face.cv));
        bool converged = true;
        gic.invert(ctx_x.xreal(), P_ref, converged);
	bool is_in = (ctx_x.pgt()->convex_ref()->is_in(P_ref) < 1E-4);
	GMM_ASSERT1(is_in && converged, "Geometric transformation inversion "
                    "has failed in neighbour transformation");
        face_num = adj_face.f;
        cv = adj_face.cv;
        ret_type = 1;
      }
      GMM_ASSERT1(!compute_derivatives,
                  "No derivative for this transformation");
      return ret_type;
    }

    interpolate_transformation_neighbour() { }

  };


  pinterpolate_transformation interpolate_transformation_neighbour_instance() {
    return (std::make_shared<interpolate_transformation_neighbour>());
  }

  //=========================================================================
  // Interpolate transformation on neighbour element (for extrapolation)
  //=========================================================================

  class interpolate_transformation_element_extrapolation
    : public virtual_interpolate_transformation, public context_dependencies {

    const mesh &sm;
    std::map<size_type, size_type> elt_corr;

  public:
    void update_from_context() const {}
    void extract_variables(const ga_workspace &/* workspace */,
                           std::set<var_trans_pair> &/* vars */,
                           bool /* ignore_data */, const mesh &/* m */,
                           const std::string &/* interpolate_name */) const {}
    void init(const ga_workspace &/* workspace */) const {}
    void finalize() const {}
    std::string expression(void) const { return "X"; }

    int transform(const ga_workspace &/*workspace*/, const mesh &m_x,
                  fem_interpolation_context &ctx_x,
                  const base_small_vector &/*Normal*/, const mesh **m_t,
                  size_type &cv, short_type &face_num,
                  base_node &P_ref,
                  base_small_vector &/*N_y*/,
                  std::map<var_trans_pair, base_tensor> &/*derivatives*/,
                  bool compute_derivatives) const {
      int ret_type = 0;
      *m_t = &m_x;
      GMM_ASSERT1(&sm == &m_x, "Bad mesh");
      size_type cv_x = ctx_x.convex_num(), cv_y = cv_x;
      auto it = elt_corr.find(cv_x);
      if (it != elt_corr.end()) cv_y = it->second;

      if (cv_x != cv_y) {
        bgeot::geotrans_inv_convex gic;
        gic.init(m_x.points_of_convex(cv_y),
                 m_x.trans_of_convex(cv_y));
        bool converged = true;
        gic.invert(ctx_x.xreal(), P_ref, converged, 1E-4);
        GMM_ASSERT1(converged, "Geometric transformation inversion "
                    "has failed in element extrapolation transformation");
        face_num = short_type(-1);
        cv = cv_y;
        ret_type = 1;
      } else {
        cv = cv_x;
        face_num = short_type(-1);
        P_ref = ctx_x.xref();
        ret_type = 1;
      }
      GMM_ASSERT1(!compute_derivatives,
                  "No derivative for this transformation");
      return ret_type;
    }

    void set_correspondance(const std::map<size_type, size_type> &ec) {
      elt_corr = ec;
    }

    interpolate_transformation_element_extrapolation
    (const mesh &sm_, const std::map<size_type, size_type> &ec)
      : sm(sm_), elt_corr(ec) { }
  };


  void add_element_extrapolation_transformation
  (model &md, const std::string &name, const mesh &sm,
   std::map<size_type, size_type> &elt_corr) {
    pinterpolate_transformation
      p = std::make_shared<interpolate_transformation_element_extrapolation>
      (sm, elt_corr);
    md.add_interpolate_transformation(name, p);
  }

  void add_element_extrapolation_transformation
  (ga_workspace &workspace, const std::string &name, const mesh &sm,
   std::map<size_type, size_type> &elt_corr) {
    pinterpolate_transformation
      p = std::make_shared<interpolate_transformation_element_extrapolation>
      (sm, elt_corr);
    workspace.add_interpolate_transformation(name, p);
  }

  void set_element_extrapolation_correspondance
  (ga_workspace &workspace, const std::string &name,
   std::map<size_type, size_type> &elt_corr) {
    GMM_ASSERT1(workspace.interpolate_transformation_exists(name),
                "Unknown transformation");
    auto pit=workspace.interpolate_transformation(name).get();
    auto cpext
      = dynamic_cast<const interpolate_transformation_element_extrapolation *>
      (pit);
    GMM_ASSERT1(cpext,
                "The transformation is not of element extrapolation type");
    const_cast<interpolate_transformation_element_extrapolation *>(cpext)
      ->set_correspondance(elt_corr);
  }

  void set_element_extrapolation_correspondance
  (model &md, const std::string &name,
   std::map<size_type, size_type> &elt_corr) {
    GMM_ASSERT1(md.interpolate_transformation_exists(name),
                "Unknown transformation");
    auto pit=md.interpolate_transformation(name).get();
    auto cpext
      = dynamic_cast<const interpolate_transformation_element_extrapolation *>
      (pit);
    GMM_ASSERT1(cpext,
                "The transformation is not of element extrapolation type");
    const_cast<interpolate_transformation_element_extrapolation *>(cpext)
      ->set_correspondance(elt_corr);
  }


  //=========================================================================
  // Secondary domains
  //=========================================================================


  class standard_secondary_domain : public virtual_secondary_domain {
    
  public:

    virtual const mesh_region &give_region(const mesh &,
					   size_type, short_type) const
    { return region; }
    // virtual void init(const ga_workspace &workspace) const = 0;
    // virtual void finalize() const = 0;

    standard_secondary_domain(const mesh_im &mim__, const mesh_region &region_)
      : virtual_secondary_domain(mim__, region_) {}
  };

  void add_standard_secondary_domain
  (model &md, const std::string &name, const mesh_im &mim,
   const mesh_region &rg) { 
    psecondary_domain p = std::make_shared<standard_secondary_domain>(mim, rg);
    md.add_secondary_domain(name, p);
  }
  
  void add_standard_secondary_domain
  (ga_workspace &workspace, const std::string &name, const mesh_im &mim,
   const mesh_region &rg) { 
    psecondary_domain p = std::make_shared<standard_secondary_domain>(mim, rg);
    workspace.add_secondary_domain(name, p);
  }
  

} /* end of namespace */
