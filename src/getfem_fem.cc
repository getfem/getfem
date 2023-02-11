/*===========================================================================

 Copyright (C) 1999-2020 Yves Renard

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

/** @file getfem_fem.cc
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @date December 21, 1999.
    @brief implementation of some finite elements.
 */

#include "getfem/bgeot_torus.h"
#include "getfem/dal_singleton.h"
#include "getfem/dal_tree_sorted.h"
#include "gmm/gmm_algobase.h"
#include "getfem/dal_naming_system.h"
#include "getfem/getfem_fem.h"
#include "getfem/getfem_integration.h"
#include "getfem/getfem_omp.h"
#include "getfem/getfem_torus.h"

namespace getfem {

  using bgeot::read_base_poly;
  using bgeot::base_rational_fraction;

  const base_matrix& fem_interpolation_context::M() const {
    if (gmm::mat_nrows(M_) == 0) {
      GMM_ASSERT2(have_pgt() && have_G() && have_pf(), "cannot compute M");
      M_.resize(pf_->nb_dof(convex_num()), pf_->nb_base(convex_num()));
      pf_->mat_trans(M_,G(),pgt());
    }
    return M_;
  }

  size_type fem_interpolation_context::convex_num() const {
    GMM_ASSERT3(convex_num_ != size_type(-1), "");
    return convex_num_;
  }

  bool fem_interpolation_context::is_convex_num_valid() const {
    return convex_num_ != size_type(-1);
  }

  short_type fem_interpolation_context::face_num() const {
    GMM_ASSERT3(face_num_ != short_type(-1),
                "Face number is asked but not defined");
    return face_num_;
  }

  bool fem_interpolation_context::is_on_face() const {
    return (face_num_ != short_type(-1));
  }

  // Comment regarding fem defined on the real element:
  // Precomputed values, gradients and hessians of fem defined on the real
  // element are not dealt with within fem_interpolation_context.
  // In that case, any precomputations can be performed within the fem
  // itself, which has to store the corresponding values internally.
  // The methods real_base_value, real_grad_base_value or real_hess_base_value
  // of the fem have the possibility to check if the passed ctx has a pfp
  // and extract the corresponding internally stored results based on
  // ctx.convex_num(), ctx.pfp()->get_ppoint_tab() and ctx.ii().
  // In that case, the storage available in ctx.pfp()->c, ctx.pfp()->pc
  // and ctx.pfp()->hpc is not used.


  // Specific multiplication for fem_interpolation_context use.
  static inline void spec_mat_tmult_(const base_tensor &g, const base_matrix &B,
                                     base_tensor &t) {
    size_type P = B.nrows(), N = B.ncols();
    size_type M = t.adjust_sizes_changing_last(g, P);
    bgeot::mat_tmult(&(*(g.begin())), &(*(B.begin())), &(*(t.begin())),M,N,P);
  }

  void fem_interpolation_context::pfp_base_value(base_tensor& t,
                                                 const pfem_precomp &pfp__) {
    const pfem &pf__ = pfp__->get_pfem();
    GMM_ASSERT1(ii_ != size_type(-1), "Internal error");

    if (pf__->is_standard())
      t = pfp__->val(ii());
    else {
      if (pf__->is_on_real_element())
        pf__->real_base_value(*this, t);
      else {
        switch(pf__->vectorial_type()) {
        case virtual_fem::VECTORIAL_NOTRANSFORM_TYPE:
          t = pfp__->val(ii()); break;
        case virtual_fem::VECTORIAL_PRIMAL_TYPE:
          t.mat_transp_reduction(pfp__->val(ii()), K(), 1); break;
        case virtual_fem::VECTORIAL_DUAL_TYPE:
          t.mat_transp_reduction(pfp__->val(ii()), B(), 1); break;
        }
        if (!(pf__->is_equivalent())) {
          set_pfp(pfp__);
          { base_tensor u = t; t.mat_transp_reduction(u, M(), 0); }
        }
      }
    }
  }


  void fem_interpolation_context::base_value(base_tensor& t,
                                             bool withM) const {
    if (pfp_ && ii_ != size_type(-1) && pf_->is_standard())
      t = pfp_->val(ii());
    else {
      if (pf_->is_on_real_element())
        pf_->real_base_value(*this, t);
      else {
        if (pfp_ && ii_ != size_type(-1)) {
          switch(pf_->vectorial_type()) {
          case virtual_fem::VECTORIAL_NOTRANSFORM_TYPE:
            t = pfp_->val(ii()); break;
          case virtual_fem::VECTORIAL_PRIMAL_TYPE:
            t.mat_transp_reduction(pfp_->val(ii()), K(), 1); break;
          case virtual_fem::VECTORIAL_DUAL_TYPE:
            t.mat_transp_reduction(pfp_->val(ii()), B(), 1); break;
          }
        }
        else {
          switch(pf_->vectorial_type()) {
          case virtual_fem::VECTORIAL_NOTRANSFORM_TYPE:
            pf_->base_value(xref(), t); break;
          case virtual_fem::VECTORIAL_PRIMAL_TYPE:
            {
              base_tensor u; pf_->base_value(xref(), u);
              t.mat_transp_reduction(u,K(),1);
            } break;
          case virtual_fem::VECTORIAL_DUAL_TYPE:
            {
              base_tensor u; pf_->base_value(xref(), u);
              t.mat_transp_reduction(u,B(),1);
            } break;
          }
        }
        if (withM && !(pf_->is_equivalent()))
          { base_tensor u = t; t.mat_transp_reduction(u, M(), 0); }
      }
    }
  }

  void fem_interpolation_context::pfp_grad_base_value
  (base_tensor& t, const pfem_precomp &pfp__) {
    const pfem &pf__ = pfp__->get_pfem();
    GMM_ASSERT1(ii_ != size_type(-1), "Internal error");

    if (pf__->is_standard()) {
      // t.mat_transp_reduction(pfp__->grad(ii()), B(), 2);
      spec_mat_tmult_(pfp__->grad(ii()), B(), t);
    } else {
      if (pf__->is_on_real_element())
        pf__->real_grad_base_value(*this, t);
      else {
        switch(pf__->vectorial_type()) {
        case virtual_fem::VECTORIAL_PRIMAL_TYPE:
          {
            base_tensor u;
            // u.mat_transp_reduction(pfp__->grad(ii()), B(), 2);
            spec_mat_tmult_(pfp__->grad(ii()), B(), u);
            t.mat_transp_reduction(u, K(), 1);
          }
          break;
        case virtual_fem::VECTORIAL_DUAL_TYPE:
          {
            base_tensor u;
            // u.mat_transp_reduction(pfp__->grad(ii()), B(), 2);
            spec_mat_tmult_(pfp__->grad(ii()), B(), u);
            t.mat_transp_reduction(u, B(), 1);
          }
          break;
        default:
          // t.mat_transp_reduction(pfp__->grad(ii()), B(), 2);
          spec_mat_tmult_(pfp__->grad(ii()), B(), t);
        }
        if (!(pf__->is_equivalent())) {
          set_pfp(pfp__);
          base_tensor u = t; t.mat_transp_reduction(u, M(), 0);
        }
      }
    }
  }


  void fem_interpolation_context::grad_base_value(base_tensor& t,
                                                  bool withM) const {
    if (pfp_ && ii_ != size_type(-1) && pf_->is_standard()) {
      // t.mat_transp_reduction(pfp_->grad(ii()), B(), 2);
      spec_mat_tmult_(pfp_->grad(ii()), B(), t);
    } else {
      if (pf()->is_on_real_element())
        pf()->real_grad_base_value(*this, t);
      else {
        if (have_pfp() && ii() != size_type(-1)) {
          switch(pf()->vectorial_type()) {
          case virtual_fem::VECTORIAL_PRIMAL_TYPE:
            {
              base_tensor u;
              // u.mat_transp_reduction(pfp_->grad(ii()), B(), 2);
              spec_mat_tmult_(pfp_->grad(ii()), B(), u);
              t.mat_transp_reduction(u, K(), 1);
            }
            break;
          case virtual_fem::VECTORIAL_DUAL_TYPE:
            {
              base_tensor u;
              // u.mat_transp_reduction(pfp_->grad(ii()), B(), 2);
              spec_mat_tmult_(pfp_->grad(ii()), B(), u);
              t.mat_transp_reduction(u, B(), 1);
            }
            break;
          default:
            // t.mat_transp_reduction(pfp_->grad(ii()), B(), 2);
            spec_mat_tmult_(pfp_->grad(ii()), B(), t);
          }

        } else {
          base_tensor u;
          pf()->grad_base_value(xref(), u);
          if (u.size()) { /* only if the FEM can provide grad_base_value */
            // t.mat_transp_reduction(u, B(), 2);
            spec_mat_tmult_(u, B(), t);
            switch(pf()->vectorial_type()) {
            case virtual_fem::VECTORIAL_PRIMAL_TYPE:
              u = t; t.mat_transp_reduction(u, K(), 1); break;
            case virtual_fem::VECTORIAL_DUAL_TYPE:
              u = t; t.mat_transp_reduction(u, B(), 1);  break;
            default: break;
            }
          }
        }
        if (withM && !(pf()->is_equivalent()))
          { base_tensor u = t; t.mat_transp_reduction(u, M(), 0); }
      }
    }
  }

  void fem_interpolation_context::hess_base_value(base_tensor& t,
                                                  bool withM) const {
    if (pf()->is_on_real_element())
      pf()->real_hess_base_value(*this, t);
    else {
      base_tensor tt;
      if (have_pfp() && ii() != size_type(-1))
        tt = pfp()->hess(ii());
      else
        pf()->hess_base_value(xref(), tt);

      switch(pf()->vectorial_type()) {
      case virtual_fem::VECTORIAL_PRIMAL_TYPE:
        { base_tensor u = tt; tt.mat_transp_reduction(u, K(), 1); } break;
      case virtual_fem::VECTORIAL_DUAL_TYPE:
        { base_tensor u = tt; tt.mat_transp_reduction(u, B(), 1); } break;
      default: break;
      }

      if (tt.size()) { /* only if the FEM can provide hess_base_value */
        tt.adjust_sizes(tt.sizes()[0], tt.sizes()[1], gmm::sqr(tt.sizes()[2]));
        t.mat_transp_reduction(tt, B3(), 2);
        if (!pgt()->is_linear()) {
          if (have_pfp()) {
            tt.mat_transp_reduction(pfp()->grad(ii()), B32(), 2);
          } else {
            base_tensor u;
            pf()->grad_base_value(xref(), u);
            tt.mat_transp_reduction(u, B32(), 2);
          }
          t -= tt;
        }
        if (!(pf()->is_equivalent()) && withM)
          { tt = t; t.mat_transp_reduction(tt, M(), 0); }
      }
    }
  }

  void fem_interpolation_context::set_pfp(pfem_precomp newpfp) {
    if (pfp_ != newpfp) {
      pfp_ = newpfp;
      if (pfp_) { pf_ = pfp()->get_pfem(); }
      else pf_ = 0;
      M_.resize(0,0);
    }
  }

  void fem_interpolation_context::set_pf(pfem newpf) {
    if (pf_ != newpf || have_pfp()) {
      set_pfp(0);
      pf_ = newpf;
    }
  }

  void virtual_fem::real_base_value(const fem_interpolation_context &c,
                                    base_tensor &t, bool withM) const
  { c.base_value(t, withM); }

  void virtual_fem::real_grad_base_value(const fem_interpolation_context &c,
                                    base_tensor &t, bool withM) const
  { c.grad_base_value(t, withM);}

  void virtual_fem::real_hess_base_value(const fem_interpolation_context &c,
                                         base_tensor &t, bool withM) const
  { c.hess_base_value(t, withM); }

  /* ******************************************************************** */
  /*        Class for description of an interpolation dof.                */
  /* ******************************************************************** */

  enum ddl_type { LAGRANGE, NORMAL_DERIVATIVE, DERIVATIVE, MEAN_VALUE,
                  BUBBLE1, LAGRANGE_NONCONFORMING, GLOBAL_DOF,
                  SECOND_DERIVATIVE, NORMAL_COMPONENT, EDGE_COMPONENT,
                  IPK_CENTER};

  struct ddl_elem {
    ddl_type t;
    gmm::int16_type hier_degree;
    short_type hier_raff;
    size_type spec;
    bool operator < (const ddl_elem &l) const {
      if (t < l.t) return true;
      if (t > l.t) return false;
      if (hier_degree < l.hier_degree) return true;
      if (hier_degree > l.hier_degree) return false;
      if (hier_raff < l.hier_raff) return true;
      if (hier_raff > l.hier_raff) return false;
      if (spec < l.spec) return true;
      return false;
    }
    ddl_elem(ddl_type s = LAGRANGE, gmm::int16_type k = -1, short_type l = 0,
             size_type spec_ = 0)
      : t(s), hier_degree(k), hier_raff(l), spec(spec_) {}
  };

  struct dof_description {
    std::vector<ddl_elem> ddl_desc;
    bool linkable;
    dim_type coord_index;
    size_type xfem_index;
    bool all_faces;

    dof_description()
    { linkable = true; all_faces = false; coord_index = 0; xfem_index = 0; }
  };

  struct dof_description_comp__ {
    int operator()(const dof_description &m, const dof_description &n) const;
  };

  // ATTENTION : en cas de modif, changer aussi dof_description_compare,
  //             product_dof, et dof_hierarchical_compatibility.
  int dof_description_comp__::operator()(const dof_description &m,
                                         const dof_description &n) const {
    int nn = gmm::lexicographical_less<std::vector<ddl_elem> >()
      (m.ddl_desc, n.ddl_desc);
    if (nn < 0) return -1;
    if (nn > 0) return 1;
    nn = int(m.linkable) - int(n.linkable);
    if (nn < 0) return -1;
    if (nn > 0) return 1;
    nn = int(m.coord_index) - int(n.coord_index);
    if (nn < 0) return -1;
    if (nn > 0) return 1;
    nn = int(m.xfem_index) - int(n.xfem_index);
    if (nn < 0) return -1;
    if (nn > 0) return 1;
    nn = int(m.all_faces) - int(n.all_faces);
    if (nn < 0) return -1;
    if (nn > 0) return 1;
    return 0;
  }

  typedef dal::dynamic_tree_sorted<dof_description, dof_description_comp__> dof_d_tab;

  pdof_description lagrange_dof(dim_type n) {
    THREAD_SAFE_STATIC dim_type n_old = dim_type(-2);
    THREAD_SAFE_STATIC pdof_description p_old = nullptr;
    if (n != n_old) {
      dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
      dof_description l;
      l.ddl_desc.resize(n);
      std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(LAGRANGE));
      p_old = &(tab[tab.add_norepeat(l)]);
      n_old = n;
    }
    return p_old;
  }

  pdof_description lagrange_0_dof(dim_type n) {
    THREAD_SAFE_STATIC dim_type n_old = dim_type(-2);
    THREAD_SAFE_STATIC pdof_description p_old = nullptr;
    if (n != n_old) {
      dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
      dof_description l;
      l.all_faces = true;
      l.ddl_desc.resize(n);
      l.linkable = false;
      std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(LAGRANGE));
      p_old = &(tab[tab.add_norepeat(l)]);
      n_old = n;
    }
    return p_old;
  }

  pdof_description deg_hierarchical_dof(pdof_description p, int deg) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l = *p;
    for (size_type i = 0; i < l.ddl_desc.size(); ++i)
      l.ddl_desc[i].hier_degree = gmm::int16_type(deg);
    return &(tab[tab.add_norepeat(l)]);
  }

  size_type reserve_xfem_index() {
    static size_type ind = 100;
    return ind += 1000;
  }

  pdof_description xfem_dof(pdof_description p, size_type ind) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l = *p; l.xfem_index = ind;
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description ipk_center_dof(dim_type n, size_type k_ind) {
    THREAD_SAFE_STATIC dim_type n_old = dim_type(-2);
    THREAD_SAFE_STATIC pdof_description p_old = nullptr;
    if (n != n_old) {
      dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
      dof_description l;
      l.ddl_desc.resize(n);
      std::fill(l.ddl_desc.begin(), l.ddl_desc.end(),
                ddl_elem(IPK_CENTER, -1, 0, k_ind));
      p_old = &(tab[tab.add_norepeat(l)]);
      n_old = n;
    }
    return p_old;
  }
  


  pdof_description to_coord_dof(pdof_description p, dim_type ct) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l = *p;
    l.coord_index = ct;
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description raff_hierarchical_dof(pdof_description p, short_type deg) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l = *p;
    for (size_type i = 0; i < l.ddl_desc.size(); ++i)
      l.ddl_desc[i].hier_raff = deg;
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description lagrange_nonconforming_dof(dim_type n) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l; l.linkable = false;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(LAGRANGE));
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description bubble1_dof(dim_type n) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(BUBBLE1));
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description derivative_dof(dim_type n, dim_type num_der) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(LAGRANGE));
    l.ddl_desc.at(num_der) = ddl_elem(DERIVATIVE);
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description second_derivative_dof(dim_type n, dim_type num_der1,
                                         dim_type num_der2) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(LAGRANGE));
    l.ddl_desc[num_der1] = ddl_elem(SECOND_DERIVATIVE);
    l.ddl_desc[num_der2] = ddl_elem(SECOND_DERIVATIVE);
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description normal_derivative_dof(dim_type n) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(),
              ddl_elem(NORMAL_DERIVATIVE));
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description normal_component_dof(dim_type n) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(),
              ddl_elem(NORMAL_COMPONENT));
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description edge_component_dof(dim_type n) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(),
              ddl_elem(EDGE_COMPONENT));
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description mean_value_dof(dim_type n) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l;
    l.ddl_desc.resize(n);
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(), ddl_elem(MEAN_VALUE));
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description global_dof(dim_type n) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    dof_description l;
    l.all_faces = true;
    l.ddl_desc.resize(n);
    l.linkable = false;
    std::fill(l.ddl_desc.begin(), l.ddl_desc.end(),ddl_elem(GLOBAL_DOF));
    return &(tab[tab.add_norepeat(l)]);
  }

  pdof_description product_dof(pdof_description a, pdof_description b) {
    dof_d_tab& tab = dal::singleton<dof_d_tab>::instance();
    size_type nb1 = a->ddl_desc.size(), nb2 = b->ddl_desc.size();
    dof_description l;
    l.linkable = a->linkable && b->linkable;
    l.coord_index = std::max(a->coord_index, b->coord_index); // logique ?
    l.xfem_index = a->xfem_index;
    l.all_faces = a->all_faces || b->all_faces;
    GMM_ASSERT2(a->xfem_index == b->xfem_index, "Invalid product of dof");
    l.ddl_desc.resize(nb1+nb2);
    std::copy(a->ddl_desc.begin(), a->ddl_desc.end(), l.ddl_desc.begin());
    std::copy(b->ddl_desc.begin(), b->ddl_desc.end(), l.ddl_desc.begin()+nb1);

    {
      gmm::int16_type deg = -1;
      for (size_type i = 0; i < l.ddl_desc.size(); ++i)
        deg = std::max(deg, l.ddl_desc[i].hier_degree);
      for (size_type i = 0; i < l.ddl_desc.size(); ++i)
        l.ddl_desc[i].hier_degree = deg;
    }
    {
      short_type deg = 0;
      for (size_type i = 0; i < l.ddl_desc.size(); ++i)
        deg = std::max(deg, l.ddl_desc[i].hier_raff);
      for (size_type i = 0; i < l.ddl_desc.size(); ++i)
        l.ddl_desc[i].hier_raff = deg;
    }
    return &(tab[tab.add_norepeat(l)]);
  }

  // ATTENTION : en cas de modif, changer aussi
  //             dof_description_comp__::operator,
  //             product_dof, et dof_hierarchical_compatibility.
  //             et dof_description_compare
  int dof_weak_compatibility(pdof_description a, pdof_description b) {
    int nn;
    std::vector<ddl_elem>::const_iterator
      ita = a->ddl_desc.begin(), itae = a->ddl_desc.end(),
      itb = b->ddl_desc.begin(), itbe = b->ddl_desc.end();
    for (; ita != itae && itb != itbe; ++ita, ++itb)
    {
      if ((nn = int(ita->t) - int (itb->t)) != 0) return nn;
      if ((nn = int(ita->hier_degree) - int (itb->hier_degree)) != 0)
        return nn;
    }
    for (; ita != itae; ++ita) if (ita->t != LAGRANGE) return 1;
    for (; itb != itbe; ++itb) if (itb->t != LAGRANGE) return -1;
    return 0;
  }


  // ATTENTION : en cas de modif, changer aussi
  //             dof_description_comp__::operator,
  //             product_dof, et dof_hierarchical_compatibility.
  int dof_description_compare(pdof_description a, pdof_description b) {
    if (a == b) return 0;
    int nn;
    if ((nn = int(a->coord_index) - int(b->coord_index)) != 0) return nn;
    if ((nn = int(a->linkable) - int(b->linkable)) != 0) return nn;
    if ((nn = int(a->xfem_index) - int(b->xfem_index)) != 0) return nn;
    return dof_weak_compatibility(a,b);
  }

  bool dof_linkable(pdof_description a)
  { return a->linkable; }

  bool dof_compatibility(pdof_description a, pdof_description b)
  { return (dof_linkable(a) && dof_description_compare(a, b) == 0); }

  size_type dof_xfem_index(pdof_description a)
  { return a->xfem_index; }

  dim_type coord_index_of_dof(pdof_description a)
  { return a->coord_index; }

  bool dof_hierarchical_compatibility(pdof_description a, pdof_description b)
  {
    if (a->coord_index != b->coord_index) return false;
    if (a->linkable != b->linkable) return false;
    if (a->xfem_index != b->xfem_index) return false;
    std::vector<ddl_elem>::const_iterator
      ita = a->ddl_desc.begin(), itae = a->ddl_desc.end(),
      itb = b->ddl_desc.begin(), itbe = b->ddl_desc.end();
    for (; ita != itae && itb != itbe; ++ita, ++itb)
    { if (ita->t != itb->t) return false; }
    for (; ita != itae; ++ita) if (ita->t != LAGRANGE) return false;
    for (; itb != itbe; ++itb) if (itb->t != LAGRANGE) return false;
    return true;
  }


  /* ******************************************************************** */
  /*        Members methods of virtual_fem           .                    */
  /* ******************************************************************** */


  void virtual_fem::add_node(const pdof_description &d, const base_node &pt,
                             const dal::bit_vector &faces) {
    short_type nb = cv_node.nb_points();
    cv_node.points().resize(nb+1);
    cv_node.points()[nb] = pt;
    dof_types_.resize(nb+1);
    face_tab.resize(nb+1);
    dof_types_[nb] = d;
    cvs_node->add_point_adaptative(nb, short_type(-1));
    for (dal::bv_visitor f(faces); !f.finished(); ++f) {
      cvs_node->add_point_adaptative(nb, short_type(f));
      face_tab[nb].push_back(short_type(f));
    }
    pspt_valid = false;
  }


  void virtual_fem::add_node(const pdof_description &d, const base_node &pt) {
    dal::bit_vector faces;
    for (short_type f = 0; f < cvs_node->nb_faces(); ++f)
      if (d->all_faces || gmm::abs(cvr->is_in_face(f, pt)) < 1.0E-7)
        faces.add(f);
     add_node(d, pt, faces);
  }

  void virtual_fem::init_cvs_node() {
    cvs_node->init_for_adaptative(cvr->structure());
    cv_node = bgeot::convex<base_node>(cvs_node);
    face_tab.resize(0);
    pspt_valid = false;
  }

  void virtual_fem::unfreeze_cvs_node() {
    cv_node.structure() = cvs_node;
    pspt_valid = false;
  }

  const std::vector<short_type> &
  virtual_fem::faces_of_dof(size_type /*cv*/, size_type i) const {
    static const std::vector<short_type> no_faces;
    return (i < face_tab.size()) ? face_tab[i] : no_faces;
  }

  void virtual_fem::copy(const virtual_fem &f) {
    dof_types_ = f.dof_types_;
    cvs_node = bgeot::new_convex_structure();
    *cvs_node = *f.cvs_node; // deep copy
    cv_node = f.cv_node;
    cv_node.structure() = cvs_node;
    pspt = 0;
    pspt_valid = false;
    cvr = f.cvr;
    dim_ = f.dim_;
    ntarget_dim = f.ntarget_dim;
    vtype = f.vtype;
    is_equiv = f.is_equiv;
    is_lag = f.is_lag;
    is_pol = f.is_pol;
    is_polycomp = f.is_polycomp;
    real_element_defined = f.real_element_defined;
    is_standard_fem = f.is_standard_fem;
    es_degree = f.es_degree;
    hier_raff = f.hier_raff;
    debug_name_ = f.debug_name_;
    face_tab = f.face_tab;
  }

  /* ******************************************************************** */
  /*        PK class.                                                         */
  /* ******************************************************************** */

  class PK_fem_ : public fem<base_poly> {
  public :
    void calc_base_func(base_poly &p, size_type i, short_type K) const;
    PK_fem_(dim_type nc, short_type k);
    ~PK_fem_() {}
  };

  void PK_fem_::calc_base_func(base_poly &p, size_type i, short_type K) const {
    dim_type N = dim();
    base_poly l0(N, 0), l1(N, 0);
    bgeot::power_index w(short_type(N+1));
    l0.one(); l1.one(); p = l0;

    if (K != 0) {
      for (short_type nn = 0; nn < N; ++nn) l0 -= base_poly(N, 1, nn);

      w[0] = K;
      for (short_type nn = 1; nn <= N; ++nn) {
        w[nn]=short_type(floor(0.5+bgeot::to_scalar((cv_node.points()[i])[nn-1]*opt_long_scalar_type(K))));
        w[0]=short_type(w[0] - w[nn]);
      }

      for (int nn = 0; nn <= N; ++nn)
        for (int j = 0; j < w[nn]; ++j) {
          if (nn == 0)
            p *= (l0 * (opt_long_scalar_type(K) / opt_long_scalar_type(j+1)))
              - (l1 * (opt_long_scalar_type(j) / opt_long_scalar_type(j+1)));
          else
            p *= (base_poly(N,1,short_type(nn-1)) * (opt_long_scalar_type(K)/opt_long_scalar_type(j+1)))
              - (l1 * (opt_long_scalar_type(j) / opt_long_scalar_type(j+1)));
        }
    }
  }

  PK_fem_::PK_fem_(dim_type nc, short_type k) {
    cvr = bgeot::simplex_of_reference(nc);
    dim_ = cvr->structure()->dim();
    is_standard_fem = is_equiv = is_pol = is_lag = true;
    es_degree = k;

    init_cvs_node();
    bgeot::pconvex_ref cvn = bgeot::simplex_of_reference(nc, k);
    size_type R = cvn->nb_points();
    for (size_type i = 0; i < R; ++i)
      add_node(k==0 ? lagrange_0_dof(nc) : lagrange_dof(nc), cvn->points()[i]);

    base_.resize(R);
    for (size_type r = 0; r < R; r++) calc_base_func(base_[r], r, k);
  }

  static pfem PK_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0,
                "Bad type of parameters");
    int n = dim_type(::floor(params[0].num() + 0.01));
    int k = short_type(::floor(params[1].num() + 0.01));
    GMM_ASSERT1(n > 0 && n < 100 && k >= 0 && k <= 150 &&
                double(n) == params[0].num() && double(k) == params[1].num(),
                "Bad parameters");
    pfem p = std::make_shared<PK_fem_>(dim_type(n), short_type(k));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }


  /* ******************************************************************** */
  /*        Tensorial product of fem (for polynomial fem).                    */
  /* ******************************************************************** */

  struct tproduct_femi : public fem<base_poly> {
    tproduct_femi(ppolyfem fi1, ppolyfem fi2);
  };

  tproduct_femi::tproduct_femi(ppolyfem fi1, ppolyfem fi2) {
    if (fi2->target_dim() != 1) std::swap(fi1, fi2);
    GMM_ASSERT1(fi2->target_dim() == 1, "dimensions mismatch");

    is_pol = true;
    is_equiv = fi1->is_equivalent() && fi2->is_equivalent();
    GMM_ASSERT1(is_equiv,
                "Product of non equivalent elements not available, sorry.");
    is_lag = fi1->is_lagrange() && fi2->is_lagrange();
    is_standard_fem = fi1->is_standard() && fi2->is_standard();
    es_degree = short_type(fi1->estimated_degree() + fi2->estimated_degree());

    bgeot::convex<base_node> cv
      = bgeot::convex_direct_product(fi1->node_convex(0), fi2->node_convex(0));
    cvr = bgeot::convex_ref_product(fi1->ref_convex(0), fi2->ref_convex(0));
    dim_ = cvr->structure()->dim();
    init_cvs_node();

    ntarget_dim = fi2->target_dim();
    base_.resize(cv.nb_points() * ntarget_dim);
    size_type i, j, r;
    for (j = 0, r = 0; j < fi2->nb_dof(0); ++j)
      for (i = 0; i < fi1->nb_dof(0); ++i, ++r)
        add_node(product_dof(fi1->dof_types()[i], fi2->dof_types()[j]),
                 cv.points()[r]);

    for (j = 0, r = 0; j < fi2->nb_base_components(0); j++)
      for (i = 0; i < fi1->nb_base_components(0); i++, ++r) {
        base_[r] = fi1->base()[i];
        base_[r].direct_product(fi2->base()[j]);
      }
  }

  static pfem product_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 1 && params[1].type() == 1,
                "Bad type of parameters");
    pfem pf1 = params[0].method();
    pfem pf2 = params[1].method();
    GMM_ASSERT1(pf1->is_polynomial() && pf2->is_polynomial(),
                "Both arguments to FEM_PRODUCT must be polynomial FEM");
    pfem p = std::make_shared<tproduct_femi>(ppolyfem(pf1.get()),
                                             ppolyfem(pf2.get()));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*    Generic Hierarchical fem (for polynomial fem). To be interfaced.  */
  /* ******************************************************************** */

  struct thierach_femi : public fem<base_poly> {
    thierach_femi(ppolyfem fi1, ppolyfem fi2);
  };

  thierach_femi::thierach_femi(ppolyfem fi1, ppolyfem fi2)
    : fem<base_poly>(*fi1) {
    grad_computed_ = false;
    hess_computed_ = false;
    GMM_ASSERT1(fi2->target_dim()==fi1->target_dim(), "dimensions mismatch.");
    GMM_ASSERT1(fi2->basic_structure(0) == fi1->basic_structure(0),
                "Incompatible elements.");
    GMM_ASSERT1(fi1->is_equivalent() &&  fi2->is_equivalent(), "Sorry, "
                "no hierachical construction for non tau-equivalent fems.");
    es_degree = fi2->estimated_degree();
    is_lag = false;
    unfreeze_cvs_node();
    for (size_type i = 0; i < fi2->nb_dof(0); ++i) {
      bool found = false;
      for (size_type j = 0; j < fi1->nb_dof(0); ++j) {
        if ( gmm::vect_dist2(fi2->node_of_dof(0,i),
                             fi1->node_of_dof(0,j)) < 1e-13
             && dof_hierarchical_compatibility(fi2->dof_types()[i],
                                               fi1->dof_types()[j]))
            { found = true; break; }
      }
      if (!found) {
        add_node(deg_hierarchical_dof(fi2->dof_types()[i],
                                      fi1->estimated_degree()),
                 fi2->node_of_dof(0,i));
        base_.resize(nb_dof(0));
        base_[nb_dof(0)-1] = (fi2->base())[i];
      }
    }
  }

  struct thierach_femi_comp : public fem<bgeot::polynomial_composite> {
    thierach_femi_comp(ppolycompfem fi1, ppolycompfem fi2);
  };

  thierach_femi_comp::thierach_femi_comp(ppolycompfem fi1, ppolycompfem fi2)
    : fem<bgeot::polynomial_composite>(*fi1) {
    GMM_ASSERT1(fi2->target_dim()==fi1->target_dim(), "dimensions mismatch.");
    GMM_ASSERT1(fi2->basic_structure(0) == fi1->basic_structure(0),
                "Incompatible elements.");
    GMM_ASSERT1(fi1->is_equivalent() &&  fi2->is_equivalent(), "Sorry, "
                "no hierachical construction for non tau-equivalent fems.");

    es_degree = std::max(fi2->estimated_degree(), fi1->estimated_degree());

    is_lag = false;
    hier_raff = short_type(fi1->hierarchical_raff() + 1);
    unfreeze_cvs_node();
    for (size_type i = 0; i < fi2->nb_dof(0); ++i) {
      bool found = false;
      for (size_type j = 0; j < fi1->nb_dof(0); ++j) {
        if ( gmm::vect_dist2(fi2->node_of_dof(0,i),
                             fi1->node_of_dof(0,j)) < 1e-13
             && dof_hierarchical_compatibility(fi2->dof_types()[i],
                                               fi1->dof_types()[j]))
          { found = true; break; }
      }
      if (!found) {
        add_node(raff_hierarchical_dof(fi2->dof_types()[i], hier_raff),
                 fi2->node_of_dof(0,i));
        base_.resize(nb_dof(0));
        base_[nb_dof(0)-1] = (fi2->base())[i];
      }
    }
  }

  static pfem gen_hierarchical_fem
  (fem_param_list &params, std::vector<dal::pstatic_stored_object> &deps) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 1 && params[1].type() == 1,
                "Bad type of parameters");
    pfem pf1 = params[0].method();
    pfem pf2 = params[1].method();
    if (pf1->is_polynomial() && pf2->is_polynomial())
      return std::make_shared<thierach_femi>(ppolyfem(pf1.get()),
                                             ppolyfem(pf2.get()));
    GMM_ASSERT1(pf1->is_polynomialcomp() && pf2->is_polynomialcomp(),
                "Bad parameters");
    pfem p  = std::make_shared<thierach_femi_comp>(ppolycompfem(pf1.get()),
                                                   ppolycompfem(pf2.get()));
    deps.push_back(p->ref_convex(0));
    deps.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /* PK hierarchical fem.                                                 */
  /* ******************************************************************** */

  static pfem PK_hierarch_fem(fem_param_list &params,
                              std::vector<dal::pstatic_stored_object> &) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0,
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01)), s;
    GMM_ASSERT1(n > 0 && n < 100 && k > 0 && k <= 150 &&
                double(n) == params[0].num() && double(k) == params[1].num(),
                "Bad parameters");
    std::stringstream name;
    if (k == 1)
      name << "FEM_PK(" << n << ",1)";
    else {
      for (s = 2; s <= k; ++s) if ((k % s) == 0) break;
      name << "FEM_GEN_HIERARCHICAL(FEM_PK_HIERARCHICAL(" << n << ","
           << k/s << "), FEM_PK(" << n << "," << k << "))";
    }
    return fem_descriptor(name.str());
  }

  static pfem QK_hierarch_fem(fem_param_list &params,
                              std::vector<dal::pstatic_stored_object> &) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0,
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    GMM_ASSERT1(n > 0 && n < 100 && k > 0 && k <= 150 &&
                double(n) == params[0].num() && double(k) == params[1].num(),
                "Bad parameters");
    std::stringstream name;
    if (n == 1)
      name << "FEM_PK_HIERARCHICAL(1," << k << ")";
    else
      name << "FEM_PRODUCT(FEM_PK_HIERARCHICAL(" << n-1 << "," << k
           << "),FEM_PK_HIERARCHICAL(1," << k << "))";
    return fem_descriptor(name.str());
  }

  static pfem prism_PK_hierarch_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0,
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    GMM_ASSERT1(n > 1 && n < 100 && k >= 0 && k <= 150 &&
                double(n) == params[0].num() && double(k) == params[1].num(),
                "Bad parameters");
    std::stringstream name;
    if (n == 2)
      name << "FEM_QK_HIERARCHICAL(1," << k << ")";
    else
      name << "FEM_PRODUCT(FEM_PK_HIERARCHICAL(" << n-1 << "," << k
           << "),FEM_PK_HIERARCHICAL(1," << k << "))";
    return fem_descriptor(name.str());
  }


  /* ******************************************************************** */
  /* parallelepiped fems.                                                 */
  /* ******************************************************************** */

  static pfem QK_fem_(fem_param_list &params, bool discontinuous) {
    const char *fempk = discontinuous ? "FEM_PK_DISCONTINUOUS" : "FEM_PK";
    const char *femqk = discontinuous ? "FEM_QK_DISCONTINUOUS" : "FEM_QK";
    GMM_ASSERT1(params.size() == 2 || (discontinuous && params.size() == 3),
                "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0 &&
                (params.size() != 3 || params[2].type() == 0),
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    char alpha[128]; alpha[0] = 0;
    if (discontinuous && params.size() == 3) {
      scalar_type v = params[2].num();
      GMM_ASSERT1(v >= 0 && v <= 1, "Bad value for alpha: " << v);
      sprintf(alpha, ",%g", v);
    }
    GMM_ASSERT1(n > 0 && n < 100 && k >= 0 && k <= 150 &&
                double(n) == params[0].num() && double(k) == params[1].num(),
                "Bad parameters");
    std::stringstream name;
    if (n == 1)
      name << fempk << "(1," << k << alpha << ")";
    else
      name << "FEM_PRODUCT(" << femqk << "(" << n-1 << ","
           << k << alpha << ")," << fempk << "(1," << k << alpha << "))";
    return fem_descriptor(name.str());
  }

  static pfem QK_fem(fem_param_list &params,
                     std::vector<dal::pstatic_stored_object> &) {
    return QK_fem_(params, false);
  }

  static pfem QK_discontinuous_fem(fem_param_list &params,
                                   std::vector<dal::pstatic_stored_object> &) {
    return QK_fem_(params, true);
  }


  /* ******************************************************************** */
  /* prims fems.                                                          */
  /* ******************************************************************** */

  static pfem prism_PK_fem(fem_param_list &params,
                           std::vector<dal::pstatic_stored_object> &) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0,
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    GMM_ASSERT1(n > 1 && n < 100 && k >= 0 && k <= 150 &&
                double(n) == params[0].num() && double(k) == params[1].num(),
                "Bad parameters");
    std::stringstream name;
    if (n == 2)
      name << "FEM_QK(1," << k << ")";
    else
      name << "FEM_PRODUCT(FEM_PK(" << n-1 << "," << k << "),FEM_PK(1,"
           << k << "))";
    return fem_descriptor(name.str());
  }

  static pfem
  prism_PK_discontinuous_fem(fem_param_list &params,
                             std::vector<dal::pstatic_stored_object> &) {
    GMM_ASSERT1(params.size() == 2 || params.size() == 3,
                "Bad number of parameters : "
                << params.size() << " should be 2 or 3.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0 &&
                (params.size() != 3 || params[2].type() == 0),
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    char alpha[128]; alpha[0] = 0;
    if (params.size() == 3) {
      scalar_type v = params[2].num();
      GMM_ASSERT1(v >= 0 && v <= 1, "Bad value for alpha: " << v);
      sprintf(alpha, ",%g", v);
    }
    GMM_ASSERT1(n > 1 && n < 100 && k >= 0 && k <= 150 &&
                double(n) == params[0].num() && double(k) == params[1].num(),
                "Bad parameters");
    std::stringstream name;
    if (n == 2)
      name << "FEM_QK_DISCONTINUOUS(1," << k << alpha << ")";
    else
      name << "FEM_PRODUCT(FEM_PK_DISCONTINUOUS(" << n-1 << "," << k << alpha
           << "),FEM_PK_DISCONTINUOUS(1,"
           << k << alpha << "))";
    return fem_descriptor(name.str());
  }

  /* ******************************************************************** */
  /*        P1 NON CONFORMING (dim 2)                                         */
  /* ******************************************************************** */

   static pfem P1_nonconforming_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 0, "Bad number of parameters ");
    auto p = std::make_shared<fem<base_poly>>();
    p->mref_convex() = bgeot::simplex_of_reference(2);
    p->dim() = 2;
    p->is_standard() = p->is_equivalent() = true;
    p->is_polynomial() = p->is_lagrange() = true;
    p->estimated_degree() = 1;
    p->init_cvs_node();
    p->base().resize(3);

    p->add_node(lagrange_dof(2), base_small_vector(0.5, 0.5));
    read_poly(p->base()[0], 2, "2*x + 2*y - 1");
    p->add_node(lagrange_dof(2), base_small_vector(0.0, 0.5));
    read_poly(p->base()[1], 2, "1 - 2*x");
    p->add_node(lagrange_dof(2), base_small_vector(0.5, 0.0));
    read_poly(p->base()[2], 2, "1 - 2*y");
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));

    return pfem(p);
  }


  /* ******************************************************************** */
  /*     Quad8/Hexa20 SERENDIPITY ELEMENT (dim 2 or 3) (incomplete Q2)    */
  /* ******************************************************************** */

  // local dof numeration for 2D:
  // 5--6--7
  // |     |
  // 3     4
  // |     |
  // 0--1--2
  //
  // local dof numeration for 3D:
  //
  //       17---18---19
  //      /|        /|
  //     / 10      / 11
  //   15  |      16 |
  //   /   5----6/---7
  //  /   /     /   /
  // 12---13---14  /
  // |  3      |  4
  // 8 /       9 /
  // |/        |/
  // 0----1----2

   static pfem
   build_Q2_incomplete_fem(fem_param_list &params,
                           std::vector<dal::pstatic_stored_object> &deps,
                           bool discontinuous) {
    GMM_ASSERT1(params.size() <= 1, "Bad number of parameters");
    dim_type n = 2;
    if (params.size() > 0) {
      GMM_ASSERT1(params[0].type() == 0, "Bad type of parameters");
      n = dim_type(::floor(params[0].num() + 0.01));
       GMM_ASSERT1(n == 2 || n == 3, "Bad parameter, expected value 2 or 3");
    }
    auto p = std::make_shared<fem<base_poly>>();
    p->mref_convex() = bgeot::parallelepiped_of_reference(n);
    p->dim() = n;
    p->is_standard() = p->is_equivalent() = true;
    p->is_polynomial() = p->is_lagrange() = true;
    p->estimated_degree() = 2;
    p->init_cvs_node();
    p->base().resize(n == 2 ? 8 : 20);
    auto lag_dof = discontinuous ? lagrange_nonconforming_dof(n) 
                                 : lagrange_dof(n);

    if (n == 2) {
      std::stringstream s
        ( "1 - 2*x^2*y - 2*x*y^2 + 2*x^2 + 5*x*y + 2*y^2 - 3*x - 3*y;"
          "4*(x^2*y - x^2 - x*y + x);"
          "2*x*y*y - 2*x*x*y + 2*x*x - x*y - x;"
          "4*(x*y*y - x*y - y*y + y);"
          "4*(x*y - x*y*y);"
          "2*x*x*y - 2*x*y*y - x*y + 2*y*y - y;"
          "4*(x*y - x*x*y);"
          "2*x*x*y + 2*x*y*y - 3*x*y;");

      for (int i = 0; i < 8; ++i) p->base()[i] = read_base_poly(2, s);

      p->add_node(lag_dof, base_small_vector(0.0, 0.0));
      p->add_node(lag_dof, base_small_vector(0.5, 0.0));
      p->add_node(lag_dof, base_small_vector(1.0, 0.0));
      p->add_node(lag_dof, base_small_vector(0.0, 0.5));
      p->add_node(lag_dof, base_small_vector(1.0, 0.5));
      p->add_node(lag_dof, base_small_vector(0.0, 1.0));
      p->add_node(lag_dof, base_small_vector(0.5, 1.0));
      p->add_node(lag_dof, base_small_vector(1.0, 1.0));
    } else {
      std::stringstream s
        ("1 + 2*x^2*y*z + 2*x*y^2*z + 2*x*y*z^2"
           " - 2*x^2*y - 2*x^2*z - 2*x*y^2 - 2*y^2*z - 2*y*z^2 - 2*x*z^2 - 7*x*y*z"
           " + 2*x^2 + 2*y^2 + 2*z^2 + 5*y*z + 5*x*z + 5*x*y - 3*x - 3*y - 3*z;"
         "4*( - x^2*y*z + x*y*z + x^2*z - x*z + x^2*y - x*y - x^2 + x);"
         "2*x^2*y*z - 2*x*y^2*z - 2*x*y*z^2"
           " - 2*x^2*y - 2*x^2*z + 2*x*y^2 + 2*x*z^2 + 3*x*y*z + 2*x^2 - x*y - x*z - x;"
         "4*( - x*y^2*z + x*y^2 + y^2*z + x*y*z - x*y - y^2 - y*z + y);"
         "4*(x*y^2*z - x*y^2 - x*y*z + x*y);"
         " - 2*x^2*y*z + 2*x*y^2*z - 2*x*y*z^2"
           " + 2*x^2*y - 2*x*y^2 - 2*y^2*z + 2*y*z^2 + 3*x*y*z - x*y + 2*y^2 - y*z - y;"
         "4*(x^2*y*z - x^2*y - x*y*z + x*y);"
         " - 2*x^2*y*z - 2*x*y^2*z + 2*x*y*z^2 + 2*x^2*y + 2*x*y^2 + x*y*z - 3*x*y;"
         "4*( - x*y*z^2 + x*z^2 + y*z^2 + x*y*z - x*z - y*z - z^2 + z);"
         "4*(x*y*z^2 - x*y*z - x*z^2 + x*z);"
         "4*(x*y*z^2 - x*y*z - y*z^2 + y*z);"
         "4*( - x*y*z^2 + x*y*z);"
         " - 2*x^2*y*z - 2*x*y^2*z + 2*x*y*z^2"
           " + 2*x^2*z + 2*y^2*z - 2*x*z^2 - 2*y*z^2 + 3*x*y*z - x*z - y*z + 2*z^2 - z;"
         "4*(x^2*y*z - x^2*z - x*y*z + x*z);"
         " - 2*x^2*y*z + 2*x*y^2*z - 2*x*y*z^2 + 2*x^2*z + 2*x*z^2 + x*y*z - 3*x*z;"
         "4*(x*y^2*z - y^2*z - x*y*z + y*z);"
         "4*( - x*y^2*z + x*y*z);"
         "2*x^2*y*z - 2*x*y^2*z - 2*x*y*z^2 + 2*y^2*z + 2*y*z^2 + x*y*z - 3*y*z;"
         "4*( - x^2*y*z + x*y*z);"
         "2*x^2*y*z + 2*x*y^2*z + 2*x*y*z^2 - 5*x*y*z;");

      for (int i = 0; i < 20; ++i) p->base()[i] = read_base_poly(3, s);

      p->add_node(lag_dof, base_small_vector(0.0, 0.0, 0.0));
      p->add_node(lag_dof, base_small_vector(0.5, 0.0, 0.0));
      p->add_node(lag_dof, base_small_vector(1.0, 0.0, 0.0));
      p->add_node(lag_dof, base_small_vector(0.0, 0.5, 0.0));
      p->add_node(lag_dof, base_small_vector(1.0, 0.5, 0.0));
      p->add_node(lag_dof, base_small_vector(0.0, 1.0, 0.0));
      p->add_node(lag_dof, base_small_vector(0.5, 1.0, 0.0));
      p->add_node(lag_dof, base_small_vector(1.0, 1.0, 0.0));

      p->add_node(lag_dof, base_small_vector(0.0, 0.0, 0.5));
      p->add_node(lag_dof, base_small_vector(1.0, 0.0, 0.5));
      p->add_node(lag_dof, base_small_vector(0.0, 1.0, 0.5));
      p->add_node(lag_dof, base_small_vector(1.0, 1.0, 0.5));

      p->add_node(lag_dof, base_small_vector(0.0, 0.0, 1.0));
      p->add_node(lag_dof, base_small_vector(0.5, 0.0, 1.0));
      p->add_node(lag_dof, base_small_vector(1.0, 0.0, 1.0));
      p->add_node(lag_dof, base_small_vector(0.0, 0.5, 1.0));
      p->add_node(lag_dof, base_small_vector(1.0, 0.5, 1.0));
      p->add_node(lag_dof, base_small_vector(0.0, 1.0, 1.0));
      p->add_node(lag_dof, base_small_vector(0.5, 1.0, 1.0));
      p->add_node(lag_dof, base_small_vector(1.0, 1.0, 1.0));
    }
    deps.push_back(p->ref_convex(0));
    deps.push_back(p->node_tab(0));

    return pfem(p);
  }

  static pfem Q2_incomplete_fem
  (fem_param_list &params, std::vector<dal::pstatic_stored_object> &deps)
  { return build_Q2_incomplete_fem(params, deps, false); }

  static pfem Q2_incomplete_discontinuous_fem
  (fem_param_list &params, std::vector<dal::pstatic_stored_object> &deps)
  { return build_Q2_incomplete_fem(params, deps, true); }

  /* ******************************************************************** */
  /*    Lagrange Pyramidal element of degree 0, 1 and 2                   */
  /* ******************************************************************** */

  // local dof numeration for K=1:
  //    4
  //   /|||
  //  / || |
  // 2-|--|-3
  // | |  | |
  // ||    ||
  // ||    ||
  // 0------1
  //
  // local dof numeration for K=2:
  //
  //    13
  //   /  |
  //  11---12
  //  |    |
  //  9----10
  //  /     |
  // 6---7---8
  // |       |
  // 3   4   5
  // |       |
  // 0---1---2

  static pfem
  build_pyramid_QK_fem(short_type k, bool disc, scalar_type alpha=0) {
    auto p = std::make_shared<fem<base_rational_fraction>>();
    p->mref_convex() = bgeot::pyramid_QK_of_reference(1);
    p->dim() = 3;
    p->is_standard() = p->is_equivalent() = true;
    p->is_polynomial() = false;
    p->is_lagrange() = true;
    p->estimated_degree() = k;
    p->init_cvs_node();
    auto lag_dof = disc ? lagrange_nonconforming_dof(3) : lagrange_dof(3);

    if (k == 0) {
      p->base().resize(1);
      p->base()[0] = read_base_poly(3, "1");
      p->add_node(lagrange_0_dof(3), base_small_vector(0.0, 0.0, 0.5));
    } else if (k == 1) {
      p->base().resize(5);
      base_rational_fraction // Q = xy/(1-z)
        Q(read_base_poly(3, "x*y"), read_base_poly(3, "1-z"));
      p->base()[0] = (read_base_poly(3, "1-x-y-z") + Q)*0.25;
      p->base()[1] = (read_base_poly(3, "1+x-y-z") - Q)*0.25;
      p->base()[2] = (read_base_poly(3, "1-x+y-z") - Q)*0.25;
      p->base()[3] = (read_base_poly(3, "1+x+y-z") + Q)*0.25;
      p->base()[4] = read_base_poly(3, "z");

      p->add_node(lag_dof, base_small_vector(-1.0, -1.0, 0.0));
      p->add_node(lag_dof, base_small_vector( 1.0, -1.0, 0.0));
      p->add_node(lag_dof, base_small_vector(-1.0,  1.0, 0.0));
      p->add_node(lag_dof, base_small_vector( 1.0,  1.0, 0.0));
      p->add_node(lag_dof, base_small_vector( 0.0,  0.0, 1.0));

    } else if (k == 2) {
      p->base().resize(14);

      base_poly xi0  = read_base_poly(3, "(1-z-x)*0.5");
      base_poly xi1  = read_base_poly(3, "(1-z-y)*0.5");
      base_poly xi2  = read_base_poly(3, "(1-z+x)*0.5");
      base_poly xi3  = read_base_poly(3, "(1-z+y)*0.5");
      base_poly x    = read_base_poly(3, "x");
      base_poly y    = read_base_poly(3, "y");
      base_poly z    = read_base_poly(3, "z");
      base_poly ones = read_base_poly(3, "1");
      base_poly un_z = read_base_poly(3, "1-z");

      std::vector<base_node> points = { base_node(-1.0, -1.0, 0.0),
                                        base_node( 0.0, -1.0, 0.0),
                                        base_node( 1.0, -1.0, 0.0),
                                        base_node(-1.0,  0.0, 0.0),
                                        base_node( 0.0,  0.0, 0.0),
                                        base_node( 1.0,  0.0, 0.0),
                                        base_node(-1.0,  1.0, 0.0),
                                        base_node( 0.0,  1.0, 0.0),
                                        base_node( 1.0,  1.0, 0.0),
                                        base_node(-0.5, -0.5, 0.5),
                                        base_node( 0.5, -0.5, 0.5),
                                        base_node(-0.5,  0.5, 0.5),
                                        base_node( 0.5,  0.5, 0.5),
                                        base_node( 0.0,  0.0, 1.0) };

      if (disc && alpha != scalar_type(0)) {
        base_node G =
          gmm::mean_value(points.begin(), points.end());
        for (auto && pt : points)
          pt = (1-alpha)*pt + alpha*G;
        for (size_type d = 0; d < 3; ++d) {
          base_poly S(1,2);
          S[0] = -alpha * G[d] / (1-alpha);
          S[1] = 1. / (1-alpha);
          xi0 = bgeot::poly_substitute_var(xi0, S, d);
          xi1 = bgeot::poly_substitute_var(xi1, S, d);
          xi2 = bgeot::poly_substitute_var(xi2, S, d);
          xi3 = bgeot::poly_substitute_var(xi3, S, d);
          x = bgeot::poly_substitute_var(x, S, d);
          y = bgeot::poly_substitute_var(y, S, d);
          z = bgeot::poly_substitute_var(z, S, d);
          un_z = bgeot::poly_substitute_var(un_z, S, d);
        }
      }
      base_rational_fraction Q(read_base_poly(3, "1"), un_z);

      p->base()[ 0] = Q*Q*xi0*xi1*(x*y-z*un_z);
      p->base()[ 1] = -Q*Q*xi0*xi1*xi2*y*4.;
      p->base()[ 2] = Q*Q*xi1*xi2*(-x*y-z*un_z);
      p->base()[ 3] = -Q*Q*xi3*xi0*xi1*x*4.;
      p->base()[ 4] = Q*Q*xi0*xi1*xi2*xi3*16.;
      p->base()[ 5] = Q*Q*xi1*xi2*xi3*x*4.;
      p->base()[ 6] = Q*Q*xi3*xi0*(-x*y-z*un_z);
      p->base()[ 7] = Q*Q*xi2*xi3*xi0*y*4.;
      p->base()[ 8] = Q*Q*xi2*xi3*(x*y-z*un_z);
      p->base()[ 9] = Q*z*xi0*xi1*4.;
      p->base()[10] = Q*z*xi1*xi2*4.;
      p->base()[11] = Q*z*xi3*xi0*4.;
      p->base()[12] = Q*z*xi2*xi3*4.;
      p->base()[13] = z*(z*2-ones);

      for (const auto &pt : points)
        p->add_node(lag_dof, pt);

    } else GMM_ASSERT1(false, "Sorry, pyramidal Lagrange fem "
                       "implemented only for degree 0, 1 or 2");

    return pfem(p);
  }


  static pfem pyramid_QK_fem
  (fem_param_list &params, std::vector<dal::pstatic_stored_object> &deps) {
    GMM_ASSERT1(params.size() <= 1, "Bad number of parameters");
    short_type k = 2;
    if (params.size() > 0) {
      GMM_ASSERT1(params[0].type() == 0, "Bad type of parameters");
      k = dim_type(::floor(params[0].num() + 0.01));
    }
    pfem p = build_pyramid_QK_fem(k, false);
    deps.push_back(p->ref_convex(0));
    deps.push_back(p->node_tab(0));
    return p;
  }

  static pfem pyramid_QK_disc_fem
  (fem_param_list &params, std::vector<dal::pstatic_stored_object> &deps) {
    GMM_ASSERT1(params.size() <= 2, "Bad number of parameters");
    short_type k = 2;
    if (params.size() > 0) {
      GMM_ASSERT1(params[0].type() == 0, "Bad type of parameters");
      k = dim_type(::floor(params[0].num() + 0.01));
    }
    scalar_type alpha(0);
    if (params.size() > 1) {
      GMM_ASSERT1(params[1].type() == 0, "Bad type of parameters");
      alpha = params[1].num();
    }
    pfem p = build_pyramid_QK_fem(k, true, alpha);
    deps.push_back(p->ref_convex(0));
    deps.push_back(p->node_tab(0));
    return p;
  }


  /* ******************************************************************** */
  /*    Incomplete quadratic Lagrange pyramidal element.                  */
  /* ******************************************************************** */

  // local dof numeration:
  //
  //    12
  //   /  |
  //  10---11
  //  |    |
  //  8----9
  //  /     |
  // 5---6---7
  // |       |
  // 3       4
  // |       |
  // 0---1---2

  static pfem build_pyramid_Q2_incomplete_fem(bool disc) {
    auto p = std::make_shared<fem<base_rational_fraction>>();
    p->mref_convex() = bgeot::pyramid_QK_of_reference(1);
    p->dim() = 3;
    p->is_standard() = p->is_equivalent() = true;
    p->is_polynomial() = false;
    p->is_lagrange() = true;
    p->estimated_degree() = 2;
    p->init_cvs_node();
    auto lag_dof = disc ? lagrange_nonconforming_dof(3) : lagrange_dof(3);

    p->base().resize(13);

    base_poly xy = read_base_poly(3, "x*y");
    base_poly z = read_base_poly(3, "z");

    base_poly p00m = read_base_poly(3, "1-z");

    base_poly pmmm = read_base_poly(3, "1-x-y-z");
    base_poly ppmm = read_base_poly(3, "1+x-y-z");
    base_poly pmpm = read_base_poly(3, "1-x+y-z");
    base_poly pppm = read_base_poly(3, "1+x+y-z");

    base_poly mmm0_ = read_base_poly(3, "-1-x-y")*0.25;
    base_poly mpm0_ = read_base_poly(3, "-1+x-y")*0.25;
    base_poly mmp0_ = read_base_poly(3, "-1-x+y")*0.25;
    base_poly mpp0_ = read_base_poly(3, "-1+x+y")*0.25;

    base_poly pp0m = read_base_poly(3, "1+x-z");
    base_poly pm0m = read_base_poly(3, "1-x-z");
    base_poly p0mm = read_base_poly(3, "1-y-z");
    base_poly p0pm = read_base_poly(3, "1+y-z");

    // (Bedrosian, 1992)
    p->base()[ 0] = mmm0_ * pmmm + base_rational_fraction(mmm0_ * xy, p00m); // N5
    p->base()[ 1] = base_rational_fraction(pp0m * pm0m * p0mm * 0.5, p00m);  // N6
    p->base()[ 2] = mpm0_ * ppmm - base_rational_fraction(mpm0_ * xy, p00m); // N7
    p->base()[ 3] = base_rational_fraction(p0pm * p0mm * pm0m * 0.5, p00m);  // N4
    p->base()[ 4] = base_rational_fraction(p0pm * p0mm * pp0m * 0.5, p00m);  // N8
    p->base()[ 5] = mmp0_ * pmpm - base_rational_fraction(mmp0_ * xy, p00m); // N3
    p->base()[ 6] = base_rational_fraction(pp0m * pm0m * p0pm * 0.5, p00m);  // N2
    p->base()[ 7] = mpp0_ * pppm + base_rational_fraction(mpp0_ * xy, p00m); // N1
    p->base()[ 8] = base_rational_fraction(pm0m * p0mm * z, p00m);           // N11
    p->base()[ 9] = base_rational_fraction(pp0m * p0mm * z, p00m);           // N12
    p->base()[10] = base_rational_fraction(pm0m * p0pm * z, p00m);           // N10
    p->base()[11] = base_rational_fraction(pp0m * p0pm * z, p00m);           // N9
    p->base()[12] = read_base_poly(3, "2*z^2-z");                            // N13

    p->add_node(lag_dof, base_small_vector(-1.0, -1.0, 0.0));
    p->add_node(lag_dof, base_small_vector( 0.0, -1.0, 0.0));
    p->add_node(lag_dof, base_small_vector( 1.0, -1.0, 0.0));
    p->add_node(lag_dof, base_small_vector(-1.0,  0.0, 0.0));
    p->add_node(lag_dof, base_small_vector( 1.0,  0.0, 0.0));
    p->add_node(lag_dof, base_small_vector(-1.0,  1.0, 0.0));
    p->add_node(lag_dof, base_small_vector( 0.0,  1.0, 0.0));
    p->add_node(lag_dof, base_small_vector( 1.0,  1.0, 0.0));
    p->add_node(lag_dof, base_small_vector(-0.5, -0.5, 0.5));
    p->add_node(lag_dof, base_small_vector( 0.5, -0.5, 0.5));
    p->add_node(lag_dof, base_small_vector(-0.5,  0.5, 0.5));
    p->add_node(lag_dof, base_small_vector( 0.5,  0.5, 0.5));
    p->add_node(lag_dof, base_small_vector( 0.0,  0.0, 1.0));

    return pfem(p);
  }


  static pfem pyramid_Q2_incomplete_fem
  (fem_param_list &params, std::vector<dal::pstatic_stored_object> &deps) {
    GMM_ASSERT1(params.size() == 0, "Bad number of parameters");
    pfem p = build_pyramid_Q2_incomplete_fem(false);
    deps.push_back(p->ref_convex(0));
    deps.push_back(p->node_tab(0));
    return p;
  }

  static pfem pyramid_Q2_incomplete_disc_fem
  (fem_param_list &params, std::vector<dal::pstatic_stored_object> &deps) {
    GMM_ASSERT1(params.size() == 0, "Bad number of parameters");
    pfem p = build_pyramid_Q2_incomplete_fem(true);
    deps.push_back(p->ref_convex(0));
    deps.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*    Incomplete quadratic Lagrange prism element.                      */
  /* ******************************************************************** */

  // local dof numeration:
  //
  //    14
  //    /|`
  //  12 | 13
  //  /  8  `
  // 9--10--11
  // |   |   |
  // |   5   |
  // 6  / `  7
  // | 3   4 |
  // |/     `|
  // 0---1---2

  static pfem build_prism_incomplete_P2_fem(bool disc) {
    auto p = std::make_shared<fem<base_rational_fraction>>();
    p->mref_convex() = bgeot::prism_of_reference(3);
    p->dim() = 3;
    p->is_standard() = p->is_equivalent() = true;
    p->is_polynomial() = false;
    p->is_lagrange() = true;
    p->estimated_degree() = 2;
    p->init_cvs_node();
    auto lag_dof = disc ? lagrange_nonconforming_dof(3) : lagrange_dof(3);

    p->base().resize(15);

    std::stringstream s
      ( "-2*y*z^2-2*x*z^2+2*z^2-2*y^2*z-4*x*y*z+5*y*z-2*x^2*z+5*x*z"
          "-3*z+2*y^2+4*x*y-3*y+2*x^2-3*x+1;"
        "4*(x*y*z+x^2*z-x*z-x*y-x^2+x);"
        "2*x*z^2-2*x^2*z-x*z+2*x^2-x;"
        "4*(y^2*z+x*y*z-y*z-y^2-x*y+y);"
        "4*(x*y-x*y*z);"
        "2*y*z^2-2*y^2*z-y*z+2*y^2-y;"
        "4*(y*z^2+x*z^2-z^2-y*z-x*z+z);"
        "4*(x*z-x*z^2);"
        "4*(y*z-y*z^2);"
        "-2*y*z^2-2*x*z^2+2*z^2+2*y^2*z+4*x*y*z-y*z+2*x^2*z-x*z-z;"
        "4*(-x*y*z-x^2*z+x*z);"
        "2*x*z^2+2*x^2*z-3*x*z;"
        "4*(-y^2*z-x*y*z+y*z);"
        "4*x*y*z;"
        "2*y*z^2+2*y^2*z-3*y*z;");

    for (int i = 0; i < 15; ++i)
      p->base()[i] = read_base_poly(3, s);

    p->add_node(lag_dof, base_small_vector(0.0, 0.0, 0.0));
    p->add_node(lag_dof, base_small_vector(0.5, 0.0, 0.0));
    p->add_node(lag_dof, base_small_vector(1.0, 0.0, 0.0));
    p->add_node(lag_dof, base_small_vector(0.0, 0.5, 0.0));
    p->add_node(lag_dof, base_small_vector(0.5, 0.5, 0.0));
    p->add_node(lag_dof, base_small_vector(0.0, 1.0, 0.0));
    p->add_node(lag_dof, base_small_vector(0.0, 0.0, 0.5));
    p->add_node(lag_dof, base_small_vector(1.0, 0.0, 0.5));
    p->add_node(lag_dof, base_small_vector(0.0, 1.0, 0.5));
    p->add_node(lag_dof, base_small_vector(0.0, 0.0, 1.0));
    p->add_node(lag_dof, base_small_vector(0.5, 0.0, 1.0));
    p->add_node(lag_dof, base_small_vector(1.0, 0.0, 1.0));
    p->add_node(lag_dof, base_small_vector(0.0, 0.5, 1.0));
    p->add_node(lag_dof, base_small_vector(0.5, 0.5, 1.0));
    p->add_node(lag_dof, base_small_vector(0.0, 1.0, 1.0));

    return pfem(p);
  }


  static pfem prism_incomplete_P2_fem
  (fem_param_list &params, std::vector<dal::pstatic_stored_object> &deps) {
    GMM_ASSERT1(params.size() == 0, "Bad number of parameters");
    pfem p = build_prism_incomplete_P2_fem(false);
    deps.push_back(p->ref_convex(0));
    deps.push_back(p->node_tab(0));
    return p;
  }

  static pfem prism_incomplete_P2_disc_fem
  (fem_param_list &params, std::vector<dal::pstatic_stored_object> &deps) {
    GMM_ASSERT1(params.size() == 0, "Bad number of parameters");
    pfem p = build_prism_incomplete_P2_fem(true);
    deps.push_back(p->ref_convex(0));
    deps.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*        P1 element with a bubble base function on a face              */
  /* ******************************************************************** */

  struct P1_wabbfoaf_ : public PK_fem_ {
    P1_wabbfoaf_(dim_type nc);
  };

  P1_wabbfoaf_::P1_wabbfoaf_(dim_type nc) : PK_fem_(nc, 1) {
    is_lag = false; es_degree = 2;
    base_node pt(nc); pt.fill(0.5);
    unfreeze_cvs_node();
    add_node(bubble1_dof(nc), pt);
    base_.resize(nb_dof(0));
    base_[nc+1] = base_[1]; base_[nc+1] *= scalar_type(1 << nc);
    for (int i = 2; i <= nc; ++i) base_[nc+1] *= base_[i];
    // Le raccord assure la continuite
    // des possibilités de raccord avec du P2 existent mais il faudrait
    // modifier qlq chose (transformer les fct de base P1)
  }

  static pfem P1_with_bubble_on_a_face(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 1, "Bad number of parameters : "
                << params.size() << " should be 1.");
    GMM_ASSERT1(params[0].type() == 0, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    GMM_ASSERT1(n > 1 && n < 100 && double(n) == params[0].num(),
                "Bad parameter");
    pfem p  = std::make_shared<P1_wabbfoaf_>(dim_type(n));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*    Element Pk on a square with internal nodes (with possible node    */
  /*  correspondance for surface in 3D (build for HHO elements).          */
  /* ******************************************************************** */
  // Not tested. To be tested

  struct CIPK_SQUARE_ : public fem<base_poly> {
    dim_type k;   // 
    mutable base_matrix K;
    base_small_vector norient;
    mutable bgeot::pgeotrans_precomp pgp;
    mutable bgeot::pgeometric_trans pgt_stored;
    // mutable pfem_precomp pfp;

    bgeot::pstored_point_tab pC;

    virtual void mat_trans(base_matrix &M, const base_matrix &G,
                           bgeot::pgeometric_trans pgt) const;
    CIPK_SQUARE_(dim_type nc_);
  };

  void CIPK_SQUARE_::mat_trans(base_matrix &M,
                              const base_matrix &G,
                              bgeot::pgeometric_trans pgt) const {
    // All the dof of this element are concentrated at the center of the element
    // This is a PK element on a square. The base functions are the monomials
    // The trans_mat is here to eliminate a potential rotation of angle k PI/2
    // The idea is to compute the gradient at the center and to sort and
    // orientate the derivative in the two directions, and exchange the base
    // functions accordingly.
    dim_type N = dim_type(G.nrows());
    gmm::copy(gmm::identity_matrix(), M);
    if (pgt != pgt_stored) {
      pgt_stored = pgt;
      pgp = bgeot::geotrans_precomp(pgt, pC, 0);
    }
    gmm::resize(K, N, 2);
    gmm::mult(G, pgp->grad(0), K);
    scalar_type s0(0), m0(M_PI);
    for (size_type i = 0; i < N; ++i, m0*=M_PI) s0 += m0 * K(i, 0);
    scalar_type s1(0), m1(M_PI);
    for (size_type i = 0; i < N; ++i, m1*=M_PI) s1 += m1 * K(i, 1);
    scalar_type a0 = gmm::sgn(s0), a1 = gmm::sgn(s1);

    bool inv = false;
    for (size_type i = 0; i < N; ++i) {
      if (K(i, 0) * a0 < K(i, 1) * a1 - 1e-6) break;
      if (K(i, 0) * a0 > K(i, 1) * a1 + 1e-6) { inv = true; break; }
    }

    if (a0 < 0.) {
      for (size_type i = 1, l = 1; i < k; ++i)
        for (size_type j = 0; j <= i; ++j, ++l)
          { if (((i-j) % 2) == 1) M(l, l) = -1.0; } // X^(i-j) Y^j
    }

    if (a1 < 0.) {
      for (size_type i = 1, l = 1; i < k; ++i)
        for (size_type j = 0; j <= i; ++j, ++l)
          { if ((j % 2) == 1) M(l, l) = -1.0; } // X^(i-j) Y^j
    }

    if (inv) {
      // std::swap(a0, a1);
      for (size_type i = 1, l = 1; i < k; ++i)
        for (size_type j = 0; j <= i; ++j, ++l)
          { M(l, l) = 0.0; M(l, l+i-2*j) = 1.0; }
    }
  }

  CIPK_SQUARE_::CIPK_SQUARE_(dim_type k_) {
    k = k_;
    pgt_stored = 0;

    cvr = bgeot::parallelepiped_of_reference(2);
    dim_ = cvr->structure()->dim();
    init_cvs_node();
    es_degree = k;
    is_pol = true;
    is_standard_fem = is_lag = is_equiv = false;
    base_.resize((k+1)*(k+2)/2);

    std::vector<base_node> C(1);
    C[0] = base_node(0.5, 0.5);
    pC = bgeot::store_point_tab(C);

    base_poly X, Y;
    read_poly(X, 2, "x-0.5"); read_poly(Y, 2, "y-0.5");

    base_[0] = bgeot::one_poly(2);
    add_node(ipk_center_dof(2,0), C[0]);

    for (size_type i = 1, l = 1; i < k; ++i)
      for (size_type j = 0; j <= i; ++j, ++l) { // X^(i-j) Y^j
        base_[l] = base_[0];
        for (size_type n = 0; n < i-j; ++n) base_[l] *= X;
        for (size_type n = 0; n < j; ++n) base_[l] *= Y;
        
        add_node(ipk_center_dof(2,l), C[0]);
      }
  }

  static pfem CIPK_SQUARE(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 1, "Bad number of parameters : "
                << params.size() << " should be 1.");
    GMM_ASSERT1(params[0].type() == 0, "Bad type of parameters");
    int k = int(::floor(params[0].num() + 0.01));
    GMM_ASSERT1(k >= 0 && k < 50 && double(k) == params[0].num(),
                "Bad parameter");
    pfem p = std::make_shared<CIPK_SQUARE_>(dim_type(k));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }



  

  /* ******************************************************************** */
  /*    Element RT0 on the simplexes.                                     */
  /* ******************************************************************** */

  struct P1_RT0_ : public fem<base_poly> {
    dim_type nc;
    mutable base_matrix K;
    base_small_vector norient;
    mutable bgeot::pgeotrans_precomp pgp;
    mutable bgeot::pgeometric_trans pgt_stored;
    // mutable pfem_precomp pfp;

    virtual void mat_trans(base_matrix &M, const base_matrix &G,
                           bgeot::pgeometric_trans pgt) const;
    P1_RT0_(dim_type nc_);
  };

  void P1_RT0_::mat_trans(base_matrix &M,
                          const base_matrix &G,
                          bgeot::pgeometric_trans pgt) const {
    dim_type N = dim_type(G.nrows());
    gmm::copy(gmm::identity_matrix(), M);
    if (pgt != pgt_stored) {
      pgt_stored = pgt;
      pgp = bgeot::geotrans_precomp(pgt, node_tab(0), 0);
      // pfp = fem_precomp(this, node_tab(0), 0);
    }
    GMM_ASSERT1(N == nc, "Sorry, this element works only in dimension " << nc);

    gmm::mult(G, pgp->grad(0), K); gmm::lu_inverse(K);
    for (unsigned i = 0; i <= nc; ++i) {
      if (!(pgt->is_linear()))
        { gmm::mult(G, pgp->grad(i), K); gmm::lu_inverse(K); }
      bgeot::base_small_vector n(nc);
      gmm::mult(gmm::transposed(K), cvr->normals()[i], n);

      M(i,i) = gmm::vect_norm2(n);
      n /= M(i,i);
      scalar_type ps = gmm::vect_sp(n, norient);
      if (ps < 0) M(i, i) *= scalar_type(-1);
      if (gmm::abs(ps) < 1E-8)
        GMM_WARNING2("RT0 : The normal orientation may be incorrect");
    }
  }

  // The dof of this RT0 are not the integral on the edges or faces of the
  // normal component but the normal component on the midpoint of each edge/face
  // The reason : easier to deal in case of nonlinear transformation (otherwise
  // an integral with several Gauss points would be necessary to compute on
  // edges / faces)
  // Shape fonctions on the reference element for nc = 2
  // sqrt(2)(x, y) ; (x-1, y) ; (x, y-1)
  // The shape functions on the real element are K \phi ||K^{-T}n_i||, where
  //   K is the gradient of the transformation. 
  P1_RT0_::P1_RT0_(dim_type nc_) {
    nc = nc_;
    pgt_stored = 0;
    gmm::resize(K, nc, nc);
    gmm::resize(norient, nc);
    norient[0] = M_PI;
    for (unsigned i = 1; i < nc; ++i) norient[i] = norient[i-1]*M_PI;

    cvr = bgeot::simplex_of_reference(nc);
    dim_ = cvr->structure()->dim();
    init_cvs_node();
    es_degree = 1;
    is_pol = true;
    is_standard_fem = is_lag = is_equiv = false;
    ntarget_dim = nc;
    vtype = VECTORIAL_PRIMAL_TYPE;
    base_.resize(nc*(nc+1));


    for (size_type j = 0; j < nc; ++j)
      for (size_type i = 0; i <= nc; ++i) {
        base_[i+j*(nc+1)] = base_poly(nc, 1, short_type(j));
        if (i-1 == j) base_[i+j*(nc+1)] -= bgeot::one_poly(nc);
        if (i == 0) base_[i+j*(nc+1)] *= sqrt(opt_long_scalar_type(nc));
      }

    base_node pt(nc);
    pt.fill(scalar_type(1)/scalar_type(nc));
    add_node(normal_component_dof(nc), pt);

    for (size_type i = 0; i < nc; ++i) {
      pt[i] = scalar_type(0);
      add_node(normal_component_dof(nc), pt);
      pt[i] = scalar_type(1)/scalar_type(nc);
    }
  }

  static pfem P1_RT0(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 1, "Bad number of parameters : "
                << params.size() << " should be 1.");
    GMM_ASSERT1(params[0].type() == 0, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    GMM_ASSERT1(n > 1 && n < 100 && double(n) == params[0].num(),
                "Bad parameter");
    pfem p = std::make_shared<P1_RT0_>(dim_type(n));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }


  /* ******************************************************************** */
  /*    Element RT0 on parallelepideds.                                   */
  /* ******************************************************************** */

  struct P1_RT0Q_ : public fem<base_poly> {
    dim_type nc;
    mutable base_matrix K;
    base_small_vector norient;
    mutable bgeot::pgeotrans_precomp pgp;
    mutable bgeot::pgeometric_trans pgt_stored;
    // mutable pfem_precomp pfp;

    virtual void mat_trans(base_matrix &M, const base_matrix &G,
                           bgeot::pgeometric_trans pgt) const;
    P1_RT0Q_(dim_type nc_);
  };

  void P1_RT0Q_::mat_trans(base_matrix &M,
                          const base_matrix &G,
                          bgeot::pgeometric_trans pgt) const {
    dim_type N = dim_type(G.nrows());
    gmm::copy(gmm::identity_matrix(), M);
    if (pgt != pgt_stored) {
      pgt_stored = pgt;
      pgp = bgeot::geotrans_precomp(pgt, node_tab(0),0);
      // pfp = fem_precomp(this, node_tab(0), 0);
    }
    GMM_ASSERT1(N == nc, "Sorry, this element works only in dimension " << nc);

    gmm::mult(G, pgp->grad(0), K); gmm::lu_inverse(K);
    for (unsigned i = 0; i < unsigned(2*nc); ++i) {
      if (!(pgt->is_linear()))
        { gmm::mult(G, pgp->grad(i), K); gmm::lu_inverse(K); }
      bgeot::base_small_vector n(nc);
      gmm::mult(gmm::transposed(K), cvr->normals()[i], n);

      M(i,i) = gmm::vect_norm2(n);
      n /= M(i,i);
      scalar_type ps = gmm::vect_sp(n, norient);
      if (ps < 0) M(i, i) *= scalar_type(-1);
      if (gmm::abs(ps) < 1E-8)
        GMM_WARNING2("RT0Q : The normal orientation may be incorrect");
    }
  }

  P1_RT0Q_::P1_RT0Q_(dim_type nc_) {
    nc = nc_;
    pgt_stored = 0;
    gmm::resize(K, nc, nc);
    gmm::resize(norient, nc);
    norient[0] = M_PI;
    for (unsigned i = 1; i < nc; ++i) norient[i] = norient[i-1]*M_PI;

    cvr = bgeot::parallelepiped_of_reference(nc);
    dim_ = cvr->structure()->dim();
    init_cvs_node();
    es_degree = 1;
    is_pol = true;
    is_standard_fem = is_lag = is_equiv = false;
    ntarget_dim = nc;
    vtype = VECTORIAL_PRIMAL_TYPE;
    base_.resize(nc*2*nc);

    for (size_type j = 0; j < size_type(nc*2*nc); ++j)
      base_[j] = bgeot::null_poly(nc);

    for (short_type i = 0; i < nc; ++i) {
      base_[2*i+i*2*nc] = base_poly(nc, 1, i);
      base_[2*i+1+i*2*nc] = base_poly(nc, 1, i) - bgeot::one_poly(nc);
    }

    base_node pt(nc); pt.fill(0.5);

    for (size_type i = 0; i < nc; ++i) {
      pt[i] = 1.0;
      add_node(normal_component_dof(nc), pt);
      pt[i] = 0.0;
      add_node(normal_component_dof(nc), pt);
      pt[i] = 0.5;
    }
  }

  static pfem P1_RT0Q(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 1, "Bad number of parameters : "
                << params.size() << " should be 1.");
    GMM_ASSERT1(params[0].type() == 0, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    GMM_ASSERT1(n > 1 && n < 100 && double(n) == params[0].num(),
                "Bad parameter");
    pfem p = std::make_shared<P1_RT0Q_>(dim_type(n));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }


  /* ******************************************************************** */
  /*    Nedelec Element.                                                  */
  /* ******************************************************************** */

  struct P1_nedelec_ : public fem<base_poly> {
    dim_type nc;
    base_small_vector norient;
    std::vector<base_small_vector> tangents;
    mutable bgeot::pgeotrans_precomp pgp;
    mutable bgeot::pgeometric_trans pgt_stored;
    mutable pfem_precomp pfp;

    virtual void mat_trans(base_matrix &M, const base_matrix &G,
                           bgeot::pgeometric_trans pgt) const;
    P1_nedelec_(dim_type nc_);
  };

  void P1_nedelec_::mat_trans(base_matrix &M, const base_matrix &G,
                              bgeot::pgeometric_trans pgt) const {
    bgeot::base_small_vector t(nc), v(nc);
    GMM_ASSERT1(G.nrows() == nc,
                "Sorry, this element works only in dimension " << nc);

    if (pgt != pgt_stored) {
      pgt_stored = pgt;
      pgp = bgeot::geotrans_precomp(pgt, node_tab(0), 0);
      pfp = fem_precomp(std::make_shared<P1_nedelec_>(nc), node_tab(0), 0);
    }
    fem_interpolation_context ctx(pgp,pfp,0,G,0);

    for (unsigned i = 0; i < nb_dof(0); ++i) {
      ctx.set_ii(i);
      gmm::mult(ctx.K(), tangents[i], t);
      t /= gmm::vect_norm2(t);
      gmm::mult(gmm::transposed(ctx.B()), t, v);
      scalar_type ps = gmm::vect_sp(t, norient);
      if (ps < 0) v *= scalar_type(-1);
      if (gmm::abs(ps) < 1E-8)
        GMM_WARNING2("Nedelec element: "
                     "The normal orientation may be incorrect");

      const bgeot::base_tensor &tt = pfp->val(i);
      for (size_type j = 0; j < nb_dof(0); ++j) {
        scalar_type a = scalar_type(0);
        for (size_type k = 0; k < nc; ++k) a += tt(j, k) * v[k];
        M(j, i) = a;
      }
    }
    // In fact matrix M is diagonal (at least for linear transformations).
    // The computation can be simplified.
    gmm::lu_inverse(M);
  }

  P1_nedelec_::P1_nedelec_(dim_type nc_) {
    nc = nc_;
    pgt_stored = 0;
    gmm::resize(norient, nc);
    norient[0] = M_PI;
    for (unsigned i = 1; i < nc; ++i) norient[i] = norient[i-1]*M_PI;

    cvr = bgeot::simplex_of_reference(nc);
    dim_ = cvr->structure()->dim();
    init_cvs_node();
    es_degree = 1;
    is_pol = true;
    is_standard_fem = is_lag = is_equiv = false;
    ntarget_dim = nc;
    vtype = VECTORIAL_DUAL_TYPE;
    base_.resize(nc*(nc+1)*nc/2);
    tangents.resize(nc*(nc+1)*nc/2);

    std::vector<base_poly> lambda(nc+1);
    std::vector<base_vector> grad_lambda(nc+1);
    lambda[0] = bgeot::one_poly(nc);
    gmm::resize(grad_lambda[0], nc);
    gmm::fill(grad_lambda[0], scalar_type(-1));
    for (size_type i = 1; i <= nc; ++i) {
      lambda[i] = base_poly(nc, 1, short_type(i-1));
      lambda[0] -= lambda[i];
      gmm::resize(grad_lambda[i], nc);
      grad_lambda[i][i-1] = 1;
    }

    size_type j = 0;
    for (size_type k = 0; k <= nc; ++k)
      for (size_type l = k+1; l <= nc; ++l, ++j) {
        for (size_type i = 0; i < nc; ++i) {
          base_[j+i*(nc*(nc+1)/2)] = lambda[k] * grad_lambda[l][i]
            - lambda[l] * grad_lambda[k][i];
          // cout << "base(" << j << "," << i << ") = " <<  base_[j+i*(nc*(nc+1)/2)] << endl;
        }

        base_node pt = (cvr->points()[k] + cvr->points()[l]) / scalar_type(2);
        add_node(edge_component_dof(nc), pt);
        tangents[j] = cvr->points()[l] - cvr->points()[k];
        tangents[j] /= gmm::vect_norm2(tangents[j]);
        // cout << "tangent(" << j << ") = " << tangents[j] << endl;
      }
  }

  static pfem P1_nedelec(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 1, "Bad number of parameters : "
                << params.size() << " should be 1.");
    GMM_ASSERT1(params[0].type() == 0, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    GMM_ASSERT1(n > 1 && n < 100 && double(n) == params[0].num(),
                "Bad parameter");
    pfem p = std::make_shared<P1_nedelec_>(dim_type(n));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }


  /* ******************************************************************** */
  /*    P1 element with a bubble base function on a face : type lagrange  */
  /* ******************************************************************** */

  struct P1_wabbfoafla_ : public PK_fem_
  { // idem elt prec mais avec raccord lagrange. A faire en dim. quelconque ..
    P1_wabbfoafla_();
  };

  P1_wabbfoafla_::P1_wabbfoafla_() : PK_fem_(2, 1) {
    unfreeze_cvs_node();
    es_degree = 2;
    base_node pt(2); pt.fill(0.5);
    add_node(lagrange_dof(2), pt);
    base_.resize(nb_dof(0));

    read_poly(base_[0], 2, "1 - y - x");
    read_poly(base_[1], 2, "x*(1 - 2*y)");
    read_poly(base_[2], 2, "y*(1 - 2*x)");
    read_poly(base_[3], 2, "4*x*y");
  }

  static pfem P1_with_bubble_on_a_face_lagrange(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 0, "Bad number of parameters");
    pfem p = std::make_shared<P1_wabbfoafla_>();
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }


  /* ******************************************************************** */
  /*    PK Gauss-Lobatto element on the segment                           */
  /* ******************************************************************** */

  static const double fem_coef_gausslob_1[4] =
    { 1.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e-01,
      1.000000000000000e+00 };

  static const double fem_coef_gausslob_2[9] =
    { 1.000000000000000e+00, -3.000000000000000e+00, 2.000000000000000e+00,
      0.000000000000000e-01, 4.000000000000000e+00, -4.000000000000000e+00,
      0.000000000000000e-01, -1.000000000000000e+00, 2.000000000000000e+00 };

  static const double fem_coef_gausslob_3[16] =
    { 1.000000000000000e+00, -6.000000000000000e+00, 1.000000000000000e+01,
      -5.000000000000000e+00, 0.000000000000000e-01, 8.090169943749474e+00,
      -1.927050983124842e+01, 1.118033988749895e+01, 0.000000000000000e-01,
      -3.090169943749474e+00, 1.427050983124842e+01, -1.118033988749895e+01,
      0.000000000000000e-01, 1.000000000000000e+00, -5.000000000000000e+00,
      5.000000000000000e+00 };

  static const double fem_coef_gausslob_4[25] =
    { 1.000000000000000e+00, -1.000000000000000e+01, 3.000000000000000e+01,
      -3.500000000000000e+01, 1.400000000000000e+01, 0.000000000000000e-01,
      1.351300497744848e+01, -5.687234826567877e+01, 7.602600995489696e+01,
      -3.266666666666667e+01, 0.000000000000000e-01, -5.333333333333333e+0,
      4.266666666666667e+01, -7.466666666666667e+01, 3.733333333333333e+01,
      0.000000000000000e-01, 2.820328355884853e+00, -2.479431840098789e+01,
      5.464065671176971e+01, -3.266666666666667e+01, 0.000000000000000e-01,
      -1.000000000000000e+00, 9.000000000000000e+00, -2.100000000000000e+01,
      1.400000000000000e+01 };

  static const double fem_coef_gausslob_5[36] =
    { 1.000000000000000e+00, -1.500000000000000e+01, 7.000000000000000e+01,
      -1.400000000000000e+02, 1.260000000000000e+02, -4.200000000000000e+01,
      0.000000000000000e-01, 2.028283187263934e+01, -1.315819844629968e+02,
      2.996878573260126e+02, -2.884609153819731e+02, 1.000722106463179e+02,
      0.000000000000000e-01, -8.072374540610696e+00, 9.849823568607377e+01,
      -2.894574355386068e+02, 3.201990314777874e+02, -1.211674570846437e+02,
      0.000000000000000e-01, 4.489369296352334e+00, -6.035445290945900e+01,
      2.203358804738940e+02, -2.856382539454310e+02, 1.211674570846437e+02,
      0.000000000000000e-01, -2.699826628380976e+00, 3.743820168638206e+01,
      -1.465663022612998e+02, 2.119001378496167e+02, -1.000722106463179e+02,
      0.000000000000000e-01, 1.000000000000000e+00, -1.400000000000000e+01,
      5.600000000000000e+01, -8.400000000000000e+01, 4.200000000000000e+01 };

  static const double fem_coef_gausslob_6[49] =
    { 1.000000000000000e+00, -2.100000000000000e+01, 1.400000000000000e+02,
      -4.200000000000000e+02, 6.300000000000000e+02, -4.620000000000000e+02,
      1.320000000000000e+02, 0.000000000000000e-01, 2.840315320583963e+01,
      -2.618707992013254e+02, 8.915455952644434e+02, -1.426720320066430e+03,
      1.086906031987037e+03, -3.182636611895653e+02, 0.000000000000000e-01,
      -1.133797045109102e+01, 1.954053162272040e+02, -8.515354909154440e+02,
      1.555570646517052e+03, -1.285566162567286e+03, 3.974636611895653e+02,
      0.000000000000000e-01, 6.400000000000000e+00, -1.216000000000000e+02,
      6.528000000000000e+02, -1.382400000000000e+03, 1.267200000000000e+03,
      -4.224000000000000e+02, 0.000000000000000e-01, -4.099929626153486e+00,
      8.051601475380118e+01, -4.643586932712080e+02, 1.089694751524101e+03,
      -1.099215804570106e+03, 3.974636611895653e+02, 0.000000000000000e-01,
      2.634746871404869e+00, -5.245053177967978e+01, 3.115485889222086e+02,
      -7.661450779747230e+02, 8.226759351503546e+02, -3.182636611895653e+02,
      0.000000000000000e-01, -1.000000000000000e+00, 2.000000000000000e+01,
      -1.200000000000000e+02, 3.000000000000000e+02, -3.300000000000000e+02,
      1.320000000000000e+02 };

  static const double fem_coef_gausslob_7[64] =
    { 1.000000000000000e+00, -2.800000000000000e+01, 2.520000000000000e+02,
      -1.050000000000000e+03, 2.310000000000000e+03, -2.772000000000000e+03,
      1.716000000000000e+03, -4.290000000000000e+02, 0.000000000000000e-01,
      3.787519721423474e+01, -4.699045385400572e+02, 2.217166528901653e+03,
      -5.195916436930171e+03, 6.469992588404291e+03, -4.101225846986446e+03,
      1.042012507936495e+03, 0.000000000000000e-01, -1.513857963869697e+01,
      3.497259983590744e+02, -2.101837798737552e+03, 5.599947843186200e+03,
      -7.539551580657127e+03, 5.032695631593761e+03, -1.325841514105659e+03,
      0.000000000000000e-01, 8.595816328530350e+00, -2.189405837038675e+02,
      1.612357003153179e+03, -4.947308401212060e+03, 7.342604827965170e+03,
      -5.255204822337236e+03, 1.457896159806284e+03, 0.000000000000000e-01,
      -5.620377978515898e+00, 1.480753190084385e+02, -1.171440824431867e+03,
      3.964008996775197e+03, -6.427195249873722e+03, 4.950068296306753e+03,
      -1.457896159806284e+03, 0.000000000000000e-01, 3.883318851088245e+00,
      -1.038534676200790e+02, 8.481025943868683e+02, -3.011828579891082e+03,
      5.186049587313396e+03, -4.248194967145850e+03, 1.325841514105659e+03,
      0.000000000000000e-01, -2.595374776640464e+00, 6.989727249649080e+01,
      -5.793475032722815e+02, 2.106096578071916e+03, -3.744900173152008e+03,
      3.192861708569018e+03, -1.042012507936495e+03, 0.000000000000000e-01,
      1.000000000000000e+00, -2.700000000000000e+01, 2.250000000000000e+02,
      -8.250000000000000e+02, 1.485000000000000e+03, -1.287000000000000e+03,
      4.290000000000000e+02 };

  static const double fem_coef_gausslob_8[81] =
    { 1.000000000000000e+00, -3.600000000000000e+01, 4.200000000000000e+02,
      -2.310000000000000e+03, 6.930000000000000e+03, -1.201200000000000e+04,
      1.201200000000000e+04, -6.435000000000000e+03, 1.430000000000000e+03,
      0.000000000000000e-01, 4.869949034318613e+01, -7.815432549965997e+02,
      4.860656931912704e+03, -1.551737634456316e+04, 2.788918324236022e+04,
      -2.854121637661796e+04, 1.553203650368708e+04, -3.490440192125469e+03,
      0.000000000000000e-01, -1.947740331442309e+01, 5.805138088476765e+02,
      -4.583922431826645e+03, 1.659300238583299e+04, -3.217606831061953e+04,
      3.461498040096808e+04, -1.950464316590829e+04, 4.495614716020147e+03,
      0.000000000000000e-01, 1.108992781389876e+01, -3.644117397881921e+02,
      3.513408770397022e+03, -1.458458798126453e+04, 3.105326914181443e+04,
      -3.569574046476449e+04, 2.111700401254369e+04, -5.050031666751820e+03,
      0.000000000000000e-01, -7.314285714285714e+00, 2.486857142857143e+02,
      -2.574628571428571e+03, 1.174674285714286e+04, -2.719451428571429e+04,
      3.347017142857143e+04,  -2.091885714285714e+04, 5.229714285714286e+03,
      0.000000000000000e-01, 5.181491353118710e+00, -1.789312751409961e+02,
      1.913693930879634e+03, -9.161425477258226e+03, 2.246586272145708e+04,
      -2.927759904600966e+04, 1.928324932147088e+04, -5.050031666751820e+03,
      0.000000000000000e-01, -3.748881746893966e+00, 1.304893011815078e+02,
      -1.418925315009559e+03, 6.967886161876575e+03, -1.767073170824305e+04,
      2.395969028817416e+04, -1.646027456225289e+04, 4.495614716020147e+03,
      0.000000000000000e-01, 2.569661265399177e+00, -8.980255438911062e+01,
      9.847166850754155e+02, -4.899241601766511e+03, 1.264999919894514e+04,
      -1.754928623032154e+04, 1.239148503331668e+04, -3.490440192125469e+03,
      0.000000000000000e-01, -1.000000000000000e+00, 3.500000000000000e+01,
      -3.850000000000000e+02, 1.925000000000000e+03, -5.005000000000000e+03,
      7.007000000000000e+03, -5.005000000000000e+03, 1.430000000000000e+03 };

  static const double fem_coef_gausslob_9[100] =
    { 1.000000000000000e+00, -4.500000000000000e+01, 6.600000000000000e+02,
      -4.620000000000000e+03, 1.801800000000000e+04, -4.204200000000000e+04,
      6.006000000000000e+04, -5.148000000000000e+04, 2.431000000000000e+04,
      -4.862000000000000e+03, 0.000000000000000e-01, 6.087629005856384e+01,
      -1.226341297551814e+03, 9.697405499616728e+03, -4.021760400071670e+04,
      9.725280603247945e+04, -1.421240159771644e+05, 1.237105621496494e+05,
      -5.906188663911807e+04, 1.190819794274692e+04, 0.000000000000000e-01,
      -2.435589341485964e+01, 9.095415690555422e+02, -9.111255870205338e+03,
      4.276661412893713e+04, -1.114146601604227e+05, 1.709573510663859e+05,
      -1.539309822784481e+05, 7.531473186776967e+04, -1.546698442965720e+04,
      0.000000000000000e-01, 1.388757697026791e+01, -5.717395054991898e+02,
      6.975542885191629e+03, -3.743822964286359e+04, 1.068054892026842e+05,
      -1.747039036097183e+05, 1.648204809285228e+05, -8.352714678107793e+04,
      1.762561894579012e+04, 0.000000000000000e-01, -9.198709522206266e+00,
      3.919017281073292e+02, -5.132147808492991e+03, 3.020136030457121e+04,
      -9.337958558844686e+04, 1.629937209113228e+05, -1.619399017586618e+05,
      8.553993533073839e+04, -1.866608440961585e+04, 0.000000000000000e-01,
      6.589286067498368e+00, -2.852085019651226e+02, 3.859417687766574e+03,
      -2.381847798089535e+04, 7.784545414265617e+04, -1.434184925463668e+05,
      1.495994578589254e+05, -8.245482435580428e+04, 1.866608440961585e+04,
      0.000000000000000e-01, -4.905768350885374e+00, 2.141208512030290e+02,
      -2.948048350518040e+03, 1.867520721718195e+04, -6.311993447254502e+04,
      1.208313444661297e+05, -1.311255887283438e+05, 7.510342373103317e+04,
      -1.762561894579012e+04, 0.000000000000000e-01, 3.659127863806492e+00,
      -1.604518938956052e+02, 2.230466872754075e+03, -1.433960781600324e+04,
      4.943623515122416e+04, -9.697372467640532e+04, 1.082245668039501e+05,
      -6.388812799914516e+04, 1.546698442965720e+04, 0.000000000000000e-01,
      -2.551909672185327e+00, 1.121770505458310e+02, -1.567380916112636e+03,
      1.015673778978859e+04, -3.539780430762936e+04, 7.040572036581639e+04,
      -7.991059497559391e+04, 4.811189484560421e+04, -1.190819794274692e+04,
      0.000000000000000e-01, 1.000000000000000e+00, -4.400000000000000e+01,
      6.160000000000000e+02, -4.004000000000000e+03, 1.401400000000000e+04,
      -2.802800000000000e+04, 3.203200000000000e+04, -1.944800000000000e+04,
      4.862000000000000e+03 };

  static const double fem_coef_gausslob_10[121] =
    { 1.000000000000000e+00, -5.500000000000000e+01, 9.900000000000000e+02,
      -8.580000000000000e+03, 4.204200000000000e+04, -1.261260000000000e+05,
      2.402400000000000e+05, -2.917200000000000e+05, 2.187900000000000e+05,
      -9.237800000000000e+04, 1.679600000000000e+04, 0.000000000000000e-01,
      7.440573476392379e+01, -1.837547309495916e+03, 1.797721877249297e+04,
      -9.362519219895049e+04, 2.909773746534557e+05, -5.668103968059107e+05,
      6.987888266998443e+05, -5.297630128145374e+05, 2.254581472606030e+05,
      -4.123982399226536e+04, 0.000000000000000e-01, -2.977479259019735e+01,
      1.361302600480045e+03, -1.684411461834327e+04, 9.915381814787320e+04,
      -3.316413447923031e+05, 6.777335995676564e+05, -8.637075009663581e+05,
      6.706702925312933e+05, -2.905859066780586e+05, 5.388962900035043e+04,
      0.000000000000000e-01, 1.699123898926703e+01, -8.563552206749944e+02,
      1.288193007592445e+04, -8.652550760766087e+04, 3.163119148667456e+05,
      -6.879421746683590e+05, 9.173107022760963e+05, -7.368809317678005e+05,
      3.277210563158369e+05, -6.203762550909711e+04, 0.000000000000000e-01,
      -1.128077599537177e+01, 5.884060270969327e+02, -9.496934299975062e+03,
      6.981839702203256e+04, -2.759867860090975e+05, 6.390150586777585e+05,
      -8.953334100127128e+05, 7.481407085083284e+05, -3.434511859876541e+05,
      6.671702685021839e+04, 0.000000000000000e-01, 8.126984126984127e+00,
      -4.307301587301587e+02, 7.184253968253968e+03, -5.536101587301587e+04,
      2.309526349206349e+05, -5.631187301587302e+05, 8.261892063492063e+05,
      -7.184253968253968e+05, 3.412520634920635e+05, -6.825041269841270e+04,
      0.000000000000000e-01, -6.131078401763724e+00, 3.277460052754944e+02,
      -5.563892337052540e+03, 4.401679638240306e+04, -1.898229640674875e+05,
      4.802970424048780e+05, -7.315927845245721e+05, 6.593462428792687e+05,
      -3.237190825145298e+05, 6.671702685021839e+04, 0.000000000000000e-01,
      4.719539826177156e+00, -2.535442364309725e+02, 4.348374949258725e+03,
      -3.493745849693315e+04, 1.537770968392435e+05, -3.987659746141970e+05,
      6.242937855878344e+05, -5.790845728346388e+05, 2.926551987751342e+05,
      -6.203762550909711e+04, 0.000000000000000e-01, -3.595976071589556e+00,
      1.937484128537626e+02, -3.342868418261306e+03, 2.710687970739982e+04,
      -1.208013807254569e+05, 3.181552127960248e+05, -5.073176789159280e+05,
      4.804304374445346e+05, -2.483103833254456e+05, 5.388962900035043e+04,
      0.000000000000000e-01, 2.539125352570295e+00, -1.370261203741934e+02,
      2.372031907702074e+03, -1.933271708314826e+04, 8.675745431426523e+04,
      -2.305316371991209e+05, 3.716008535065897e+05, -3.564317671210514e+05,
      1.869400926620506e+05, -4.123982399226536e+04, 0.000000000000000e-01,
      -1.000000000000000e+00, 5.400000000000000e+01, -9.360000000000000e+02,
      7.644000000000000e+03, -3.439800000000000e+04, 9.172800000000000e+04,
      -1.485120000000000e+05, 1.432080000000000e+05, -7.558200000000000e+04,
      1.679600000000000e+04 };

  static const double fem_coef_gausslob_11[144] =
    { 1.000000000000000e+00, -6.600000000000000e+01, 1.430000000000000e+03,
      -1.501500000000000e+04, 9.009000000000000e+04, -3.363360000000000e+05,
      8.168160000000000e+05, -1.312740000000000e+06, 1.385670000000000e+06,
      -9.237800000000000e+05, 3.527160000000000e+05, -5.878600000000000e+04,
      0.000000000000000e-01, 8.928790437897459e+01, -2.652104227906408e+03,
      3.141784850363868e+04, -2.002791716233576e+05, 7.743819221375186e+05,
      -1.922871126014229e+06, 3.137025618338122e+06, -3.346678943928588e+06,
      2.248625248986712e+06, -8.636670995581690e+05, 1.446085194818801e+05,
      0.000000000000000e-01, -3.573451504477809e+01, 1.963011183403142e+03,
      -2.937610003978132e+04, 2.114542549822465e+05, -8.792000471014727e+05,
      2.288870840565664e+06, -3.858042797257275e+06, 4.213930780912515e+06,
      -2.881507044264936e+06, 1.121761824282757e+06, -1.898189887480762e+05,
      0.000000000000000e-01, 2.040222049137883e+01,  -1.235400296692680e+03,
      2.244501976036562e+04, -1.840644193090062e+05, 8.352985709834925e+05,
      -2.311501023251837e+06, 4.072373742578464e+06, -4.597525083798586e+06,
      3.224564284996090e+06, -1.280533828091591e+06, 2.201577342088092e+05,
      0.000000000000000e-01, -1.356395972416294e+01, 8.500434612026409e+02,
      -1.656519758544980e+04, 1.484886632288921e+05, -7.274015626850215e+05,
      2.139270125221522e+06, -3.953929242730888e+06, 4.636483776329532e+06,
      -3.352298848558998e+06, 1.364514095203156e+06, -2.393982879242230e+05,
      0.000000000000000e-01, 9.803369220897211e+00, -6.243148521745071e+02,
      1.257271925920887e+04, -1.180754347844014e+05, 6.096877374094766e+05,
      -1.885008008173690e+06, 3.641310114569300e+06, -4.434918589444567e+06,
      3.311645240511804e+06, -1.385401912847929e+06, 2.488026449837513e+05,
      0.000000000000000e-01, -7.447686911212630e+00, 4.784415903433910e+02,
      -9.808275331323619e+03, 9.456732952354300e+04, -5.045513347291155e+05,
      1.617062776783016e+06, -3.237833360324115e+06, 4.079238919323811e+06,
      -3.141771586138831e+06, 1.351427181973335e+06, -2.488026449837513e+05,
      0.000000000000000e-01, 5.819619912607334e+00, -3.757783851118444e+02,
      7.785050290867036e+03, -7.625636515182549e+04, 4.153153824802320e+05,
      -1.363841143951918e+06, 2.804561170833446e+06, -3.631789084056234e+06,
      2.874063732359706e+06, -1.268867071963298e+06, 2.393982879242230e+05,
      0.000000000000000e-01, -4.587084978600364e+00, 2.971291972095372e+02,
      -6.195597982053585e+03, 6.128651035627590e+04, -3.381847677955126e+05,
      1.128582073344179e+06, -2.364480249965060e+06, 3.125557361498120e+06,
      -2.527901385564680e+06, 1.141201248205309e+06, -2.201577342088092e+05,
      0.000000000000000e-01, 3.549738472143437e+00, -2.303803824096343e+02,
      4.822860472410032e+03, -4.799737703690675e+04, 2.670306747478879e+05,
      -9.003482951717774e+05, 1.909697516429202e+06, -2.560483668180433e+06,
      2.103933182581560e+06, -9.662470519460815e+05, 1.898189887480762e+05,
      0.000000000000000e-01, -2.529605817247381e+00, 1.643527121363626e+02,
      -3.448327347881922e+03, 3.443600981453996e+04, -1.924805754474854e+05,
      6.528637806490704e+05, -1.394862512471197e+06, 1.886334531344429e+06,
      -1.565422824908427e+06, 7.270266147425120e+05, -1.446085194818801e+05,
      0.000000000000000e-01, 1.000000000000000e+00, -6.500000000000000e+01,
      1.365000000000000e+03, -1.365000000000000e+04, 7.644000000000000e+04,
      -2.598960000000000e+05, 5.569200000000000e+05, -7.558200000000000e+05,
      6.298500000000000e+05, -2.939300000000000e+05, 5.878600000000000e+04 };

  static const double fem_coef_gausslob_12[169] =
    { 1.000000000000000e+00, -7.800000000000000e+01, 2.002000000000000e+03,
      -2.502500000000000e+04, 1.801800000000000e+05, -8.168160000000000e+05,
      2.450448000000000e+06, -4.988412000000000e+06, 6.928350000000000e+06,
      -6.466460000000000e+06, 3.879876000000000e+06, -1.352078000000000e+06,
      2.080120000000000e+05, 0.000000000000000e-01, 1.055228477237708e+02,
      -3.710649283521791e+03, 5.230891096204477e+04, -4.000264992878070e+05,
      1.877737871587605e+06, -5.758751500257082e+06, 1.189876582421287e+07,
      -1.670085336595664e+07, 1.570840983751716e+07, -9.480360124877962e+06,
      3.318799039873028e+06, -5.124248673374161e+05, 0.000000000000000e-01,
      -4.223530748650643e+01, 2.744602756239140e+03, -4.883026612748412e+04,
      4.213447894212839e+05, -2.125569662252139e+06, 6.831231541213443e+06,
      -1.457745185889376e+07, 2.094131826427298e+07, -2.004063010887189e+07,
      1.225627209853367e+07, -4.335340118801754e+06, 6.749529540569050e+05,
      0.000000000000000e-01, 2.412126940135046e+01, -1.727728082317750e+03,
      3.727953470399318e+04, -3.660428937429080e+05, 2.013286677665669e+06,
      -6.871455230919223e+06, 1.531439984708661e+07, -2.272429867457359e+07,
      2.229291249432404e+07, -1.390087225188044e+07, 4.993770899110666e+06,
      -7.872767949619026e+05, 0.000000000000000e-01, -1.605014152757175e+01,
      1.189832345018574e+03, -2.753035300171167e+04, 2.951729666043859e+05,
      -1.750245295939125e+06, 6.340418364700963e+06, -1.480658414760223e+07,
      2.279585234765233e+07, -2.303126179585568e+07, 1.470734589654657e+07,
      -5.387526099148482e+06, 8.631843338394749e+05, 0.000000000000000e-01,
      1.162284069412542e+01, -8.756167723813171e+02, 2.093616690467957e+04,
      -2.350848402511507e+05, 1.467905987404880e+06, -5.583024412516426e+06,
      1.360724316266146e+07, -2.172800271110029e+07, 2.264080373496021e+07,
      -1.484050586374976e+07, 5.558088637639386e+06, -9.074958680213037e+05,
      0.000000000000000e-01, -8.865800865800866e+00, 6.738008658008658e+02,
      -1.640173160173160e+04, 1.890632034632035e+05, -1.219313593073593e+06,
      4.803100813852814e+06, -1.211898237229437e+07, 1.998830268398268e+07,
      -2.144876606060606e+07, 1.443281454545455e+07, -5.532578909090909e+06,
      9.220964848484848e+05, 0.000000000000000e-01, 6.984318980773296e+00,
      -5.335955916987455e+02, 1.312836634638491e+04, -1.537652068533838e+05,
      1.012269836848560e+06, -4.084347297774684e+06, 1.057602476941974e+07,
      -1.790936242524428e+07, 1.971847079705798e+07, -1.359625813912256e+07,
      5.331861778616258e+06, -9.074958680213037e+05, 0.000000000000000e-01,
      -5.596679201904262e+00, 4.289927383042951e+02, -1.062596939367885e+04,
      1.257256555151274e+05, -8.388435073376655e+05, 3.440109149730765e+06,
      -9.074697250266149e+06, 1.567950042058767e+07, -1.762881516112804e+07,
      1.241472483931862e+07, -4.970685906925217e+06, 8.631843338394749e+05,
      0.000000000000000e-01, 4.489137849846207e+00, -3.448281543178960e+02,
      8.578250876817155e+03, -1.021659510074640e+05, 6.876731048919487e+05,
      -2.851145716716003e+06, 7.618634882796075e+06, -1.335715271315875e+07,
      1.525930546501226e+07, -1.092966082914868e+07, 4.453550640432165e+06,
      -7.872767949619026e+05, 0.000000000000000e-01, -3.514808358523940e+00,
      2.703477424140389e+02, -6.743800341396192e+03, 8.065306199056715e+04,
      -5.459331868312813e+05, 2.279586137602263e+06, -6.143562568432354e+06,
      1.087848437431968e+07, -1.256803423488744e+07, 9.114425759470104e+06,
      -3.764095329881106e+06, 6.749529540569050e+05, 0.000000000000000e-01,
      2.522322790441051e+00, -1.941585635394145e+02, 4.850890672082830e+03,
      -5.815428585185421e+04, 3.949277670351431e+05, -1.655905848916830e+06,
      4.485333711312109e+06, -7.989838200781784e+06, 9.294715032477446e+06,
      -6.793611930544113e+06, 2.830299368175965e+06, -5.124248673374161e+05,
      0.000000000000000e-01, -1.000000000000000e+00, 7.700000000000000e+01,
      -1.925000000000000e+03, 2.310000000000000e+04, -1.570800000000000e+05,
      6.597360000000000e+05, -1.790712000000000e+06, 3.197700000000000e+06,
      -3.730650000000000e+06, 2.735810000000000e+06, -1.144066000000000e+06,
      2.080120000000000e+05 };

  static const double fem_coef_gausslob_13[196] =
    { 1.000000000000000e+00, -9.100000000000000e+01, 2.730000000000000e+03,
      -4.004000000000000e+04, 3.403400000000000e+05, -1.837836000000000e+06,
      6.651216000000000e+06, -1.662804000000000e+07, 2.909907000000000e+07,
      -3.556553000000000e+07, 2.974571600000000e+07, -1.622493600000000e+07,
      5.200300000000000e+06, -7.429000000000000e+05, 0.000000000000000e-01,
      1.231105960095249e+02, -5.057514000718615e+03, 8.362619816879477e+04,
      -7.548172444491492e+05, 4.219784854542762e+06, -1.560990588170226e+07,
      3.960523865471669e+07, -7.003645336672149e+07, 8.625846242015506e+07,
      -7.256273938600024e+07, 3.975792117062384e+07, -1.278833059440045e+07,
      1.832147578471145e+06, 0.000000000000000e-01, -4.927732466948325e+01,
      3.738734011094227e+03, -7.796486012083982e+04, 7.935560370804071e+05,
      -4.765562636392240e+06, 1.846680880578437e+07, -4.837506802694851e+07,
      8.753277424234758e+07, -1.096660444159086e+08, 9.346798694406962e+07,
      -5.173889595224656e+07, 1.677849814552107e+07, -2.419777739872722e+06,
      0.000000000000000e-01, 2.814884111282567e+01, -2.353904695653627e+03,
      5.948276501192433e+04, -6.883051529064705e+05, 4.502895601024915e+06,
      -1.851735708867110e+07, 5.063079101609929e+07, -9.458217148027056e+07,
      1.214199817163488e+08, -1.054743067017811e+08, 5.927651305189011e+07,
      -1.946011726944643e+07, 2.834919298555227e+06, 0.000000000000000e-01,
      -1.874042032746388e+01, 1.621968976936385e+03, -4.394233902947216e+04,
      5.547892506224127e+05, -3.908876121817967e+06, 1.704431633766850e+07,
      -4.878627691196220e+07, 9.448003328782282e+07, -1.248200449226441e+08,
      1.109678546438889e+08, -6.355498649510845e+07, 2.119358718443647e+07,
      -3.128057142433586e+06, 0.000000000000000e-01, 1.358783086373605e+01,
      -1.195146716951513e+03, 3.345811195278515e+04, -4.422483366082136e+05,
      3.278781781106079e+06, -1.499532485039955e+07, 4.474689669452104e+07,
      -8.978037194074775e+07, 1.222039761714010e+08, -1.114085938050346e+08,
      6.517880226738924e+07, -2.213159774622623e+07, 3.317403211532315e+06,
      0.000000000000000e-01, -1.039065134180050e+01, 9.220321824274881e+02,
      -2.627964906979335e+04, 3.565631308351814e+05, -2.729347397416285e+06,
      1.291899991743804e+07, -3.987098284679090e+07, 8.253644504125924e+07,
      -1.155541226713629e+08, 1.080161713852185e+08, -6.460510650895925e+07,
      2.236736005147009e+07, -3.410612094153076e+06, 0.000000000000000e-01,
      8.225051804237637e+00, -7.337438600306013e+02, 2.113982918743374e+04,
      -2.914573383452468e+05, 2.277144431399071e+06, -1.103660550091489e+07,
      3.493361321847434e+07, -7.418006034176242e+07, 1.064417028079657e+08,
      -1.018292957440870e+08, 6.222452923525811e+07, -2.197059717251990e+07,
      3.410612094153076e+06, 0.000000000000000e-01,  -6.651370533890156e+00,
      5.953674398464600e+02, -1.727143617692029e+04, 2.405949101477657e+05,
      -1.905359074321061e+06, 9.386078020968711e+06, -3.025905095154839e+07,
      6.552811535463406e+07, -9.594395490329804e+07, 9.365009838355809e+07,
      -5.835707981219508e+07, 2.099464400369387e+07, -3.317403211532315e+06,
      0.000000000000000e-01, 5.430796129527214e+00, -4.871978584294892e+02,
      1.419769025501944e+04, -1.991370307462581e+05, 1.591472100521582e+06,
      -7.928247009292726e+06, 2.589562028603176e+07, -5.690356974983853e+07,
      8.463743197871011e+07, -8.398458536550263e+07, 5.322039739169053e+07,
      -1.947115566720015e+07, 3.128057142433586e+06, 0.000000000000000e-01,
      -4.414467779036191e+00, 3.966097971621681e+02, -1.159268854505275e+04,
      1.633445673328560e+05, -1.313458735263122e+06, 6.593624701202045e+06,
      -2.173390279168494e+07, 4.826160481317836e+07, -7.262663174126566e+07,
      7.298651647234048e+07, -4.687881110584063e+07, 1.739383361177152e+07,
      -2.834919298555227e+06, 0.000000000000000e-01, 3.487743182287134e+00,
      -3.136500314357888e+02, 9.185689380767778e+03, -1.298134042763866e+05,
      1.048017196494400e+06, -5.287706376855868e+06, 1.753577438278919e+07,
      -3.921741432164137e+07, 5.949694434313328e+07, -6.033542452985016e+07,
      3.913958191606598e+07, -1.467861247282431e+07, 2.419777739872722e+06,
      0.000000000000000e-01, -2.516624450464659e+00, 2.264447557529062e+02,
      -6.639311014646825e+03, 9.399061131310191e+04, -7.605959998781342e+05,
      3.848998924774717e+06, -1.281093272369737e+07, 2.877371846174006e+07,
      -4.386952078323465e+07, 4.473878170318014e+07, -2.920546515856783e+07,
      1.102958792572444e+07, -1.832147578471145e+06, 0.000000000000000e-01,
      1.000000000000000e+00, -9.000000000000000e+01, 2.640000000000000e+03,
      -3.740000000000000e+04, 3.029400000000000e+05, -1.534896000000000e+06,
      5.116320000000000e+06, -1.151172000000000e+07, 1.758735000000000e+07,
      -1.797818000000000e+07, 1.176753600000000e+07, -4.457400000000000e+06,
      7.429000000000000e+05 };

  static const double fem_coef_gausslob_14[225] =
    { 1.000000000000000e+00, -1.050000000000000e+02, 3.640000000000000e+03,
      -6.188000000000000e+04, 6.126120000000000e+05, -3.879876000000000e+06,
      1.662804000000000e+07, -4.988412000000000e+07, 1.066965900000000e+08,
      -1.636014380000000e+08, 1.784742960000000e+08, -1.352078000000000e+08,
      6.760390000000000e+07, -2.005830000000000e+07, 2.674440000000000e+06,
      0.000000000000000e-01, 1.420511699552723e+02, -6.740724197514907e+03,
      1.291563810650317e+05, -1.357536885890637e+06, 8.899789739175488e+06,
      -3.898284725729828e+07, 1.186784002318284e+08, -2.564866709310747e+08,
      3.962828733570607e+08, -4.348015522628213e+08, 3.308658899677545e+08,
      -1.660170000478420e+08, 4.939776003229131e+07, -6.601663651220939e+06,
      0.000000000000000e-01, -5.686066782897496e+01, 4.980782943180571e+03,
      -1.202886759366057e+05, 1.425067611910860e+06, -1.003204881513823e+07,
      4.601736532886305e+07, -1.446080878593359e+08, 3.197256124417587e+08,
      -5.024241189222256e+08, 5.584380089699472e+08, -4.292685524711295e+08,
      2.171358326952990e+08, -6.503152653250142e+07, 8.737812306213077e+06,
      0.000000000000000e-01, 3.248522677540712e+01, -3.136209432322820e+03,
      9.172216152088957e+04, -1.234458149638161e+06, 9.460578149039340e+06,
      -4.602710110035357e+07, 1.508977142152744e+08, -3.443000936718750e+08,
      5.541917935770469e+08, -6.276282328622294e+08, 4.896977188070307e+08,
      -2.507043130148833e+08, 7.583045543918066e+07, -1.027267982590785e+07,
      0.000000000000000e-01, -2.163547691954172e+01, 2.161829695966707e+03,
      -6.777232432848026e+04, 9.945601317987818e+05, -8.202377657972944e+06,
      4.227975782999525e+07, -1.449995018490026e+08, 3.427553248466447e+08,
      -5.674380093336526e+08, 6.573468625138306e+08, -5.224446315267705e+08,
      2.715766154508912e+08, -8.319461131269656e+07, 1.139164303704407e+07,
      0.000000000000000e-01, 1.569978564006007e+01, -1.594280678736617e+03,
      5.164364569735884e+04, -7.932250762767995e+05, 6.879605840139093e+06,
      -3.716431671531196e+07, 1.327627119052969e+08, -3.248633554644047e+08,
      5.536614212284358e+08, -6.572276042331931e+08, 5.332102141003966e+08,
      -2.820526436178851e+08, 8.770029066945359e+07, -1.216316370145453e+07,
      0.000000000000000e-01, -1.202518395631915e+01, 1.231993083463710e+03,
      -4.063141778456666e+04, 6.405521495953326e+05, -5.734055805513453e+06,
      3.204057307537486e+07, -1.182863821502541e+08, 2.983631937198568e+08,
      -5.225422090403253e+08, 6.354190974001179e+08, -5.265537919663995e+08,
      2.837551732890861e+08, -8.968009595208465e+07, 1.261735673043107e+07,
      0.000000000000000e-01, 9.547785547785548e+00, -9.834219114219114e+02,
      3.278709557109557e+04, -5.252427785547786e+05, 4.798602442890443e+06,
      -2.744701911421911e+07, 1.038669217715618e+08, -2.685490364568765e+08,
      4.816180870862471e+08, -5.987952711608392e+08, 5.064437616783217e+08,
      -2.780475554312354e+08, 8.937242853146853e+07, -1.276748979020979e+07,
      0.000000000000000e-01, -7.763592643695499e+00, 8.024013730981778e+02,
      -2.693903659000494e+04, 4.360799342555947e+05, -4.038452021764897e+06,
      2.347605516322762e+07, -9.046087127754393e+07, 2.384165701554799e+08,
      -4.360078989902932e+08, 5.526354677146948e+08, -4.761786531169400e+08,
      2.660933883812130e+08, -8.696289827395033e+07, 1.261735673043107e+07,
      0.000000000000000e-01, 6.402657115106499e+00, -6.632652203626004e+02,
      2.237191510124190e+04, -3.647008347515229e+05, 3.408912135873189e+06,
      -2.004238739169268e+07, 7.824760241311718e+07, -2.092325230710293e+08,
      3.885803431690498e+08, -5.004334616015021e+08, 4.381904244262927e+08,
      -2.487967617473505e+08, 8.258400115090981e+07, -1.216316370145453e+07,
      0.000000000000000e-01, -5.303584550798654e+00, 5.502727059685037e+02,
      -1.861988467344103e+04, 3.050015660700474e+05, -2.869271809771707e+06,
      1.700462338220501e+07, -6.701518639186630e+07, 1.811217810430339e+08,
      -3.403535526125257e+08, 4.438883801280693e+08, -3.938531369776322e+08,
      2.266861847568459e+08, -7.628839120592035e+07, 1.139164303704407e+07,
      0.000000000000000e-01, 4.356127432886088e+00, -4.524531153029460e+02,
      1.534317874813675e+04, -2.521565333504139e+05, 2.382646312619400e+06,
      -1.419908547341188e+07, 5.633074020272954e+07, -1.534171364617007e+08,
      2.907942363862238e+08, -3.828802350952767e+08, 3.432339697459343e+08,
      -1.997222564631490e+08, 6.798706212352923e+07, -1.027267982590785e+07,
      0.000000000000000e-01, -3.466327490120249e+00, 3.602867454806401e+02,
      -1.223518159744950e+04, 2.015152850494290e+05, -1.909713800970267e+06,
      1.142278748458884e+07, -4.551909122318173e+07, 1.246206804545239e+08,
      -2.376275441309645e+08, 3.149824199011395e+08, -2.844660497989077e+08,
      1.668669076381706e+08, -5.729784575448167e+07, 8.737812306213077e+06,
      0.000000000000000e-01, 2.512080922932578e+00, -2.612119914965073e+02,
      8.878143206793750e+03, -1.464124202177336e+05, 1.389929291394544e+06,
      -8.332053211967151e+06, 3.329158201137647e+07, -9.143262460433713e+07,
      1.749809182259230e+08, -2.329047114119376e+08, 2.113183971320489e+08,
      -1.245975118891604e+08, 4.302553108480184e+07, -6.601663651220939e+06,
      0.000000000000000e-01, -1.000000000000000e+00, 1.040000000000000e+02,
      -3.536000000000000e+03, 5.834400000000000e+04, -5.542680000000000e+05,
      3.325608000000000e+06, -1.330243200000000e+07, 3.658168800000000e+07,
      -7.011490200000000e+07, 9.348653600000000e+07, -8.498776000000000e+07,
      5.022004000000000e+07, -1.738386000000000e+07, 2.674440000000000e+06 };

  static const double fem_coef_gausslob_16[289] =
    { 1.000000000000000e+00, -1.360000000000000e+02, 6.120000000000000e+03,
      -1.356600000000000e+05, 1.763580000000000e+06, -1.481407200000000e+07,
      8.535727200000000e+07, -3.505745100000000e+08, 1.051723530000000e+09,
      -2.337163400000000e+09, 3.866943080000000e+09, -4.745793780000000e+09,
      4.259045700000000e+09, -2.714556600000000e+09, 1.163381400000000e+09,
      -3.005401950000000e+08, 3.535767000000000e+07, 0.000000000000000e-01,
      1.839908474110340e+02, -1.132675577023645e+04, 2.828774748172428e+05,
      -3.903228397113182e+06, 3.393217989215092e+07, -1.997936813387011e+08,
      8.326183533203730e+08, -2.523654441487500e+09, 5.650496513849965e+09,
      -9.402288564814163e+09, 1.159004125992500e+10, -1.043753556167033e+10,
      6.671119113102151e+09, -2.865585732582308e+09, 7.416762022559170e+08,
      -8.739414676532781e+07, 0.000000000000000e-01, -7.365158560983221e+01,
      8.363752486870456e+03, -2.630512922725150e+05, 4.088268892503375e+06,
      -3.814296099037553e+07, 2.350888858038512e+08, -1.010915772862718e+09,
      3.133749897384450e+09, -7.134585101325640e+09, 1.202392366962845e+10,
      -1.496979454364896e+10, 1.358833701849443e+10, -8.740756338653663e+09,
      3.774397412355719e+09, -9.811765345365573e+08, 1.160408606498937e+08,
      0.000000000000000e-01, 4.208515410469884e+01, -5.266887544418802e+03,
      2.004067164817591e+05, -3.534528216742270e+06, 3.586506864599599e+07,
      -2.342572871648561e+08, 1.050195570375994e+09, -3.357626237614668e+09,
      7.826157982970905e+09, -1.343314504492630e+10, 1.696909064108448e+10,
      -1.558480944857237e+10, 1.012167818105129e+10, -4.405611470443590e+09,
      1.152926420144654e+09, -1.371250292488943e+08, 0.000000000000000e-01,
      -2.804154671449467e+01, 3.632134644860107e+03, -1.481030987143817e+05,
      2.845430205439061e+06, -3.103475950537538e+07, 2.145184209516461e+08,
      -1.004950951655481e+09, 3.325504487572721e+09, -7.965634686855131e+09,
      1.397533464514550e+10, -1.797134014690632e+10, 1.674913384367476e+10,
      -1.101142723314074e+10, 4.842306972881947e+09, -1.278281396603255e+09,
      1.531698732399162e+08, 0.000000000000000e-01, 2.036791285538472e+01,
      -2.681212491930484e+03, 1.129589658306899e+05, -2.270501508829428e+06,
      2.601887714173907e+07, -1.882644410377163e+08, 9.175357517098464e+08,
      -3.139134215049928e+09, 7.731773974122787e+09, -1.388518223176127e+10,
      1.820883321323428e+10, -1.725391802961138e+10, 1.150422262078233e+10,
      -5.120395806896533e+09, 1.365808881323657e+09, -1.651383905702324e+08,
      0.000000000000000e-01, -1.562975981216744e+01, 2.075857199589987e+03,
      -8.904128307330524e+04, 1.836683464233000e+06, -2.171339625537650e+07,
      1.623702309776589e+08, -8.168673425840398e+08, 2.877184244039678e+09,
      -7.272633187920932e+09, 1.336161532561319e+10, -1.787465369706492e+10,
      1.723415209556249e+10, -1.166677727073840e+10, 5.262202876034756e+09,
      -1.420107794637259e+09, 1.734782145645626e+08, 0.000000000000000e-01,
      1.245158740600359e+01, -1.662689739514946e+03, 7.210078018619273e+04,
      -1.511262926531350e+06, 1.823010400853032e+07, -1.394732137018548e+08,
      7.186625912313578e+08, -2.591802100670653e+09, 6.699969496687534e+09,
      -1.256822106464210e+10, 1.713562111741701e+10, -1.680796707246793e+10,
      1.155571599892056e+10, -5.285087289024682e+09, 1.444204616664480e+09,
      -1.784123720377501e+08, 0.000000000000000e-01, -1.018430458430458e+01,
      1.364696814296814e+03, -5.959855042735043e+04, 1.262405659052059e+06,
      -1.543602456068376e+07, 1.199989722604507e+08, -6.293065120124320e+08,
      2.311744565308469e+09, -6.087583637383061e+09, 1.162721665412277e+10,
      -1.612769282864336e+10, 1.607722369253147e+10, -1.122097126220979e+10,
      5.203928701314685e+09, -1.440373122685315e+09, 1.800466403356643e+08,
      0.000000000000000e-01, 8.484036081962663e+00, -1.139564172918365e+03,
      5.000628114583172e+04, -1.066865683702600e+06, 1.316848914076596e+07,
      -1.035421266853741e+08, 5.500823983897466e+08, -2.049399257475671e+09,
      5.477078740096760e+09, -1.061962761323504e+10, 1.495184835525608e+10,
      -1.512401891411349e+10, 1.070494963879464e+10, -5.031502683587495e+09,
      1.410393335939522e+09, -1.784123720377501e+08, 0.000000000000000e-01,
      -7.151250283691663e+00, 9.621468006776697e+02, -4.236328367402523e+04,
      9.083924059514159e+05, -1.128778312795541e+07, 8.948673396532559e+07,
      -4.799806642461592e+08, 1.807454740270899e+09, -4.886699456126079e+09,
      9.591077381959552e+09, -1.367409274689175e+10, 1.400781324267719e+10,
      -1.004054471299105e+10, 4.777971704223384e+09, -1.355543638395742e+09,
      1.734782145645626e+08, 0.000000000000000e-01, 6.060147075943005e+00,
      -8.163167553056976e+02, 3.602890134025652e+04, -7.753708269797396e+05,
      9.681484150831325e+06, -7.721340002337801e+07, 4.170906086685726e+08,
      -1.583343835764609e+09, 4.319156776154927e+09, -8.559301071696637e+09,
      1.232825943540154e+10, -1.276387222258454e+10, 9.248884856115279e+09,
      -4.449869455469566e+09, 1.276405367800062e+09, -1.651383905702324e+08,
      0.000000000000000e-01, -5.123522643555085e+00, 6.907394286218810e+02,
      -3.053901287681582e+04, 6.589382296622791e+05, -8.256408016568484e+06,
      6.613528180437879e+07, -3.591109349416936e+08, 1.371451685008549e+09,
      -3.766497172573159e+09, 7.519828793959344e+09, -1.091857986975193e+10,
      1.140164818726860e+10, -8.336452758217789e+09, 4.048470812623064e+09,
      -1.172436575235404e+09, 1.531698732399162e+08, 0.000000000000000e-01,
      4.271888353674923e+00, -5.762713061877649e+02, 2.550919052304026e+04,
      -5.514258520673562e+05, 6.926418073801134e+06, -5.565457178175802e+07,
      3.033329010654617e+08, -1.163492208450654e+09, 3.211252052317023e+09,
      -6.446888262886903e+09, 9.417864122967573e+09, -9.899668972442366e+09,
      7.289624669351123e+09, -3.566718678141100e+09, 1.041074047837655e+09,
      -1.371250292488943e+08, 0.000000000000000e-01, -3.434977405385288e+00,
      4.635617485626016e+02, -2.053688028669370e+04, 4.444943513906060e+05,
      -5.592632684586483e+06, 4.503253959356561e+07, -2.460675245081598e+08,
      9.466718497394927e+08, -2.621823641133609e+09, 5.284002694315866e+09,
      -7.752423037115038e+09, 8.187712309040204e+09, -6.060153271928369e+09,
      2.981652672294607e+09, -8.754772358617425e+08, 1.160408606498937e+08,
      0.000000000000000e-01, 2.505373764729175e+00, -3.381913429670035e+02,
      1.499009100007397e+04, -3.246847962658711e+05, 4.089321087106837e+06,
      -3.296978262323839e+07, 1.804331430493320e+08, -6.954301078105766e+08,
      1.930060872117711e+09, -3.899125665782255e+09, 5.735918309736332e+09,
      -6.075963842786729e+09, 4.511802094762441e+09, -2.227740310582889e+09,
      6.566301459893279e+08, -8.739414676532781e+07, 0.000000000000000e-01,
      -1.000000000000000e+00, 1.350000000000000e+02, -5.985000000000000e+03,
      1.296750000000000e+05, -1.633905000000000e+06, 1.318016700000000e+07,
      -7.217710500000000e+07, 2.783974050000000e+08, -7.733261250000000e+08,
      1.563837275000000e+09, -2.303105805000000e+09, 2.442687975000000e+09,
      -1.816357725000000e+09, 8.981988750000000e+08, -2.651825250000000e+08,
      3.535767000000000e+07 };

  static const double fem_coef_gausslob_24[625] =
    { 1.000000000000000e+00, -3.000000000000000e+02, 2.990000000000000e+04,
      -1.480050000000000e+06, 4.351347000000000e+07, -8.412604200000000e+08,
      1.141710570000000e+10, -1.137633032250000e+11, 8.595449577000000e+11,
      -5.042663751840000e+12, 2.337962284944000e+13, -8.678799391080000e+13,
      2.603639817324000e+14, -6.351736697208000e+14, 1.264298066396640e+15,
      -2.054484357894540e+15, 2.719170473683950e+15, -2.914666390092600e+15,
      2.505590405518200e+15, -1.701164012167620e+15, 8.910859111354200e+14,
      -3.471763290138000e+14, 9.468445336740000e+13, -1.612380184155000e+13,
      1.289904147324000e+12, 0.000000000000000e-01, 4.058641271792565e+02,
      -5.527892747506228e+04, 3.080680865091051e+06, -9.608544697986574e+07,
      1.921815551511292e+09, -2.664514371195282e+10, 2.693344195607988e+11,
      -2.055620699316044e+12, 1.214897536784061e+13, -5.664114215003022e+13,
      2.111637034506248e+14, -6.356401824720572e+14, 1.554903769748071e+15,
      -3.101863582054997e+15, 5.049759901096915e+15, -6.693732300831693e+15,
      7.184248666000200e+15, -6.182729556103940e+15, 4.201716500300169e+15,
      -2.202700032171432e+15, 8.588067464630930e+14, -2.343664562972425e+14,
      3.993223187351215e+13, -3.196139551478179e+12, 0.000000000000000e-01,
      -1.624756704002590e+02, 4.076566744849383e+04, -2.856559180759908e+06,
      1.002242332396244e+08, -2.149192199655287e+09, 3.116591524589708e+10,
      -3.248555446646729e+11, 2.534404897909630e+12, -1.522400128496926e+13,
      7.186059820369571e+13, -2.704952845330370e+14, 8.204876425573068e+14,
      -2.019503996546278e+15, 4.049106649016715e+15, -6.619542027745681e+15,
      8.805466617020362e+15, -9.478911036302426e+15, 8.178282706527972e+15,
      -5.570062256804815e+15, 2.925593222099543e+15, -1.142548313369483e+15,
      3.122535907953383e+14, -5.327144417235540e+13, 4.268671053543734e+12,
      0.000000000000000e-01, 9.285790471348058e+01, -2.567292232023452e+04,
      2.172505008713478e+06, -8.632693629032278e+07, 2.009759260312978e+09,
      -3.083881003217180e+10, 3.346965047874316e+11, -2.690206053619648e+12,
      1.652940573703274e+13, -7.940281485029880e+13, 3.030598877175188e+14,
      -9.295754861643471e+14, 2.308921995608790e+15, -4.664329550291535e+15,
      7.673380080847903e+15, -1.026158278560169e+16, 1.109639715400819e+16,
      -9.611028022799504e+15, 6.567868332535721e+15, -3.459765425094874e+15,
      1.354622808613475e+15, -3.710479852977045e+14, 6.342841599973944e+13,
      -5.091588188794019e+12, 0.000000000000000e-01, -6.190263201810896e+01,
      1.771295817306943e+04, -1.605426884400296e+06, 6.937138005497782e+07,
      -1.732266817595608e+09, 2.807090718452186e+10, -3.177491729312126e+11,
      2.638958096874092e+12, -1.663806369766778e+13, 8.158796865135113e+13,
      -3.166342096198069e+14, 9.845663378570931e+14, -2.473337869447428e+15,
      5.044014753850598e+15, -8.364663482728635e+15, 1.126254226203764e+16,
      -1.225027170970453e+16, 1.066427200479037e+16, -7.319785257851210e+15,
      3.870748078365627e+15, -1.520689469004945e+15, 4.177883147633562e+14,
      -7.160927970988626e+13, 5.762006100161049e+12, 0.000000000000000e-01,
      4.501042175430940e+01, -1.308958048352837e+04, 1.225547369430789e+06,
      -5.535760993480313e+07, 1.449945861634959e+09, -2.454369660090688e+10,
      2.883865544213490e+11, -2.470900824125910e+12, 1.598637745987517e+13,
      -8.009303566237079e+13, 3.164491556576764e+14, -9.988976501079246e+14,
      2.541436559580714e+15, -5.239264041161704e+15, 8.769361705836964e+15,
      -1.190220432935039e+16, 1.303612471279278e+16, -1.141726255498977e+16,
      7.878336148801056e+15, -4.185657107040718e+15, 1.651238357905402e+15,
      -4.553312625304459e+14, 7.830179054235796e+13, -6.319165567954503e+12,
      0.000000000000000e-01, -3.460823180120079e+01, 1.015469610544880e+04,
      -9.679531870145017e+05, 4.485134755229398e+07, -1.210735945225935e+09,
      2.114609948086286e+10, -2.559531795223062e+11, 2.252595670603948e+12,
      -1.492191456512238e+13, 7.630937241452699e+13, -3.068986895462340e+14,
      9.837309848042438e+14, -2.536329824251836e+15, 5.289429958022946e+15,
      -8.942835505684019e+15, 1.224496508118832e+16, -1.351567177210176e+16,
      1.191830727462528e+16, -8.273935552638878e+15, 4.419525275428065e+15,
      -1.751884239411906e+15, 4.851652818736986e+14, -8.375521716296401e+13,
      6.782865257514837e+12, 0.000000000000000e-01, 2.766582877253528e+01,
      -8.161941725332768e+03, 7.865526415735044e+05, -3.702889412503837e+07,
      1.019390720188427e+09, -1.819645578676150e+10, 2.252249000861092e+11,
      -2.025483052842698e+12, 1.369084263890056e+13, -7.131369560391163e+13,
      2.915943262851778e+14, -9.485944639479267e+14, 2.478119277999363e+15,
      -5.228788122037598e+15, 8.932615077176731e+15, -1.234455468758280e+16,
      1.373835480205531e+16, -1.220422101884528e+16, 8.528504524807064e+15,
      -4.582579019928659e+15, 1.826238758482371e+15, -5.082003015125067e+14,
      8.811547249928783e+13, -7.164301017225998e+12, 0.000000000000000e-01,
      -2.275670591599455e+01, 6.737590226239552e+03, -6.539504206517175e+05,
      3.111139106156313e+07, -8.679722967903971e+08, 1.573365420642779e+10,
      -1.979909637705803e+11, 1.810880611538906e+12, -1.244463025448231e+13,
      6.585374856883214e+13, -2.732735801573156e+14, 9.011914571652850e+14,
      -2.383831456114384e+15, 5.087291841992724e+15, -8.780953301699755e+15,
      1.224889816733013e+16, -1.374781660040675e+16, 1.230672077620068e+16,
      -8.660226214361099e+15, 4.682883135541219e+15, -1.876981572450838e+15,
      5.250670186065034e+14, -9.147680701891190e+13, 7.470231264323625e+12,
      0.000000000000000e-01, 1.912965144613660e+01, -5.677631395079384e+03,
      5.537935689545739e+05, -2.653927821530689e+07, 7.474036306614312e+08,
      -1.369940651896189e+10, 1.745319507247928e+11, -1.617301631376405e+12,
      1.126327452432594e+13, -6.039298017983792e+13, 2.538313181524926e+14,
      -8.473116290006231e+14, 2.267097817730376e+15, -4.890112526330036e+15,
      8.524655979372967e+15, -1.200076621013731e+16, 1.358349523272153e+16,
      -1.225446423734989e+16, 8.685297402148473e+15, -4.727405822382261e+15,
      1.906317065479628e+15, -5.362492238636365e+14, 9.390516767804315e+13,
      -7.704880889559378e+12, 0.000000000000000e-01, -1.635497697368700e+01,
      4.862656953497025e+03, -4.759804636482565e+05, 2.293041632316653e+07,
      -6.502015538842764e+08, 1.201606379408697e+10, -1.545199230484324e+11,
      1.446437458063616e+12, -1.018096106218994e+13, 5.518468570532575e+13,
      -2.344620385089091e+14, 7.909885622725649e+14, -2.138165417503776e+15,
      4.657339854616682e+15, -8.194527986057532e+15, 1.163730420250179e+16,
      -1.328057957790062e+16, 1.207345046837354e+16, -8.618482475976474e+15,
      4.722435616707614e+15, -1.916176031917556e+15, 5.421466879431477e+14,
      -9.544980641239802e+13, 7.870911362255844e+12, 0.000000000000000e-01,
      1.417066207624127e+01, -4.218698704765500e+03, 4.140273580736256e+05,
      -2.002373113338682e+07, 5.706909517514257e+08, -1.061235738178459e+10,
      1.374488760853391e+11, -1.296867134188179e+12, 9.206001697332345e+12,
      -5.034424091354864e+13, 2.158419782218546e+14, -7.348173021384713e+14,
      2.004252206567683e+15, -4.404149060382797e+15, 7.815178748758451e+15,
      -1.118956431220418e+16, 1.286957131473114e+16, -1.178684055128811e+16,
      8.473168167228462e+15, -4.673705358375533e+15, 1.908297467260950e+15,
      -5.431049066701974e+14, 9.614927445703371e+13, -7.969947411626695e+12,
      0.000000000000000e-01, -1.240846755882427e+01, 3.697723332529632e+03,
      -3.636177333437864e+05, 1.763791694375029e+07, -5.046596469793725e+08,
      9.429433336134134e+09, -1.228099190218494e+11, 1.166008419408402e+12,
      -8.333618884154624e+12, 4.590449180645647e+13, -1.982963080519099e+14,
      6.803133908337801e+14, -1.870091239145240e+15, 4.141349396659428e+15,
      -7.405302748248103e+15, 1.068239700855010e+16, -1.237594459251991e+16,
      1.141465416121965e+16, -8.261228940134621e+15, 4.586380576952005e+15,
      -1.884249466545215e+15, 5.394272826690079e+14, -9.603440259637641e+13,
      8.002866883031367e+12, 0.000000000000000e-01, 1.095558179466470e+01,
      -3.267249008495488e+03, 3.217786804294041e+05, -1.564425754031200e+07,
      4.489762785029423e+08, -8.420409805207488e+09, 1.101506623685581e+11,
      -1.051033145830932e+12, 7.553210200363392e+12, -4.185258869687206e+13,
      1.819278279112502e+14, -6.282336154227903e+14, 1.738507104235632e+15,
      -3.876120340749998e+15, 6.978305450391071e+15, -1.013471839778258e+16,
      1.182005159615111e+16, -1.097352991130819e+16, 7.992846614334175e+15,
      -4.464988119249731e+15, 1.845417602986293e+15, -5.313770797673897e+14,
      9.512946342200696e+13, -7.969947411626695e+12, 0.000000000000000e-01,
      -9.733404845062563e+00, 2.904495367336991e+03, -2.863957452026712e+05,
      1.394908621565543e+07, -4.012835563442479e+08, 7.548227155070059e+09,
      -9.908687724892865e+10, 9.492474279583278e+11, -6.852122094227041e+12,
      3.815223404048325e+13, -1.667054040084672e+14, 5.788252023371522e+14,
      -1.610924210936436e+15, 3.612762299121459e+15, -6.543084510709907e+15,
      9.560030643675167e+15, -1.121725598566102e+16, 1.047660019645038e+16,
      -7.676343347474552e+15, 4.313320840279778e+15, -1.792974677700824e+15,
      5.191726764406062e+14, -9.345206628174223e+13, 7.870911362255844e+12,
      0.000000000000000e-01, 8.685152146887849e+00, -2.592917301003156e+03,
      2.559159080755008e+05, -1.248235381982623e+07, 3.597715754514247e+08,
      -6.783361410773893e+09, 8.929618961811727e+10, -8.582135820103382e+11,
      6.217423158689400e+12, -3.475607425638306e+13, 1.525197220141970e+14,
      -5.320009405626489e+14, 1.487763351500956e+15, -3.353349463718517e+15,
      6.104800041187937e+15, -8.967037300250809e+15, 1.057820092203265e+16,
      -9.933456336414407e+15, 7.318037585534789e+15, -4.134330534453634e+15,
      1.727837357443638e+15, -5.029774927870322e+14, 9.101197367138191e+13,
      -7.704880889559378e+12, 0.000000000000000e-01, -7.768227503689160e+00,
      2.320048262018063e+03, -2.291579825895192e+05, 1.118998178170061e+07,
      -3.230127412726743e+08, 6.101825984880789e+09, -8.050592964465223e+10,
      7.757517929934399e+11, -5.636578411502579e+12, 3.161187883403309e+13,
      -1.392153217132889e+14, 4.874510273126649e+14, -1.368719368735070e+15,
      3.098228256202772e+15, -5.665515780120550e+15, 8.360206053329256e+15,
      -9.909089467205730e+15, 9.350136076665120e+15, -6.922098442148973e+15,
      3.930003596385763e+15, -1.650608740098542e+15, 4.828842861248500e+14,
      -8.780874332485509e+13, 7.470231264323625e+12, 0.000000000000000e-01,
      6.949251795654158e+00, -2.076080736426886e+03, 2.051850668187582e+05,
      -1.002851556374879e+07, 2.898385274463058e+08, -5.483488755332865e+09,
      7.247948130849744e+10, -6.998845716368053e+11, 5.097508994937879e+12,
      -2.866481138828522e+13, 1.266058930533223e+14, -4.447041797385597e+14,
      1.252927498139101e+15, -2.846337193255717e+15, 5.224630116918579e+15,
      -7.740148279482860e+15, 9.211839907560781e+15, -8.729031202403334e+15,
      6.490342376708594e+15, -3.701195553992629e+15, 1.561498591338377e+15,
      -4.588915147832623e+14, 8.382775191413614e+13, -7.164301017225998e+12,
      0.000000000000000e-01, -6.200545901231066e+00, 1.852852310460863e+03,
      -1.832115058838774e+05, 8.961081566680860e+06, -2.592406841374813e+08,
      4.910586592056635e+09, -6.500190149790611e+10, 6.287466876199138e+11,
      -4.588252565931313e+12, 2.585696625207004e+13, -1.144768233740891e+14,
      4.031460020145847e+14, -1.139023606360627e+15, 2.595327972551887e+15,
      -4.779021130665880e+15, 7.103675353349177e+15, -8.483943776806492e+15,
      8.068562477979859e+15, -6.021870692282251e+15, 3.447372991345801e+15,
      -1.460201300789598e+15, 4.308660981996213e+14, -7.903354901739207e+13,
      6.782865257514837e+12, 0.000000000000000e-01, 5.497265260133587e+00,
      -1.643010914375697e+03, 1.625245542407611e+05, -7.953853255844527e+06,
      2.302798044521235e+08, -4.366227075820252e+09, 5.786337032883471e+10,
      -5.604566478207931e+11, 4.096239643506256e+12, -2.312433390636899e+13,
      1.025754064938261e+14, -3.619933583977279e+14, 1.025085118853682e+15,
      -2.341436140046740e+15, 4.322778692729082e+15, -6.443312152718818e+15,
      7.717747123310279e+15, -7.362354286134902e+15, 5.512353176524103e+15,
      -3.166155510128892e+15, 1.345687520087759e+15, -3.984797768116558e+14,
      7.335818308855012e+13, -6.319165567954503e+12, 0.000000000000000e-01,
      -4.814420396783881e+00, 1.439137261510244e+03, -1.424001050225189e+05,
      6.972107765750786e+06, -2.019777804539574e+08, 3.832494908082456e+09,
      -5.083618287186299e+10, 4.929144471099037e+11, -3.606960360420349e+12,
      2.039001482865880e+13, -9.058350359554217e+13, 3.202053356214614e+14,
      -9.083926524953095e+14, 2.078951013598256e+15, -3.846222877659622e+15,
      5.745792065253601e+15, -6.898564020656367e+15, 6.597335811773795e+15,
      -4.952528004193797e+15, 2.852412393699805e+15, -1.215806035913631e+15,
      3.610885650804218e+14, -6.667886669397892e+13, 5.762006100161049e+12,
      0.000000000000000e-01, 4.122502768575858e+00, -1.232445305915745e+03,
      1.219756720552666e+05, -5.974119340666288e+06, 1.731450552812862e+08,
      -3.287266438655767e+09, 4.363384253154989e+10, -4.234185297825935e+11,
      3.101259925467429e+12, -1.754945237689139e+13, 7.805398522969187e+13,
      -2.762644893141324e+14, 7.848217578377621e+14, -1.798840649081661e+15,
      3.333370625337169e+15, -4.988259107765304e+15, 6.000070847701146e+15,
      -5.749271353120764e+15, 4.324788417805135e+15, -2.496262406568326e+15,
      1.066418114121039e+15, -3.174727574108465e+14, 5.876970053131701e+13,
      -5.091588188794019e+12, 0.000000000000000e-01, -3.378098091794216e+00,
      1.009981094027902e+03, -9.997415292076045e+04, 4.897701324351262e+06,
      -1.419932385393241e+08, 2.696914741487987e+09, -3.581511558046563e+10,
      3.477438352489987e+11, -2.548653261914288e+12, 1.443296944374584e+13,
      -6.424560812216257e+13, 2.275969917931459e+14, -6.472060159591627e+14,
      1.485016620793322e+15, -2.755030677045359e+15, 4.127938042773006e+15,
      -4.971860898078912e+15, 4.770796079822418e+15, -3.594142516031414e+15,
      2.077829100777859e+15, -8.891455208945626e+14, 2.651635856092349e+14,
      -4.917666111269423e+13, 4.268671053543734e+12, 0.000000000000000e-01,
      2.493031698759675e+00, -7.454011644125376e+02, 7.379166798110049e+04,
      -3.615566630393531e+06, 1.048426846839652e+08, -1.991802210178859e+09,
      2.645916950749144e+10, -2.569938254028299e+11, 1.884300408926143e+12,
      -1.067564580288449e+13, 4.754491961430585e+13, -1.685282542848975e+14,
      4.795322158961934e+14, -1.101030336950954e+15, 2.044141709763629e+15,
      -3.065196721720791e+15, 3.694953407319221e+15, -3.548705940334544e+15,
      2.676012339710790e+15, -1.548605987126603e+15, 6.633870802695019e+14,
      -1.980596394144405e+14, 3.677511736196415e+13, -3.196139551478179e+12,
      0.000000000000000e-01, -1.000000000000000e+00, 2.990000000000000e+02,
      -2.960100000000000e+04, 1.450449000000000e+06, -4.206302100000000e+07,
      7.991973990000000e+08, -1.061790830100000e+10, 1.031453949240000e+11,
      -7.563995627760000e+11, 4.286264189064000e+12, -1.909335866037600e+13,
      6.769463525042400e+13, -1.926693464819760e+14, 4.425043232388240e+14,
      -8.217937431578160e+14, 1.232690614736724e+15, -1.486479858947226e+15,
      1.428186531145374e+15, -1.077403874372826e+15, 6.237601377947940e+14,
      -2.673257733406260e+14, 7.985055567317400e+13, -1.483389769422600e+13,
      1.289904147324000e+12 };

  static const double fem_coef_gausslob_32[1089] =
    { 1.000000000000000e+00, -5.280000000000000e+02, 9.275200000000000e+04,
      -8.115800000000000e+06, 4.236447600000000e+08, -1.462986571200000e+10,
      3.573867195360000e+11, -6.471252385884000e+12, 8.987850535950000e+13,
      -9.826716585972000e+14, 8.629643838226320e+15, -6.184578084062196e+16,
      3.663173172867608e+17, -1.811459261308158e+18, 7.539120925634905e+18,
      -2.657540126286304e+19, 7.972620378858912e+19, -2.042658293145551e+20,
      4.479513800757788e+20, -8.416770667739633e+20, 1.354699278902855e+21,
      -1.864910695632502e+21, 2.189242990525111e+21, -2.181310950704368e+21,
      1.832301198591669e+21, -1.285429763935079e+21, 7.434251911077520e+20,
      -3.481117958361696e+20, 1.286127324517868e+20, -3.607069737728273e+19,
      7.214139475456547e+18, -9.163120704712953e+17, 5.553406487704820e+16,
      0.000000000000000e-01, 7.143214682034540e+02, -1.714133648848133e+05,
      1.688198752795928e+07, -9.347155797218437e+08, 3.338933321339879e+10,
      -8.331872976904955e+11, 1.530331866800279e+13, -2.146891160964000e+14,
      2.364531150945150e+15,  -2.087974965454055e+16, 1.502766434956253e+17,
      -8.930913114593062e+17, 4.428285547618542e+18, -1.847053591535166e+19,
      6.522650112096549e+19, -1.959752299775458e+20, 5.027465769775919e+20,
      -1.103711061615762e+21, 2.075747280986809e+21, -3.343657099618051e+21,
      4.606186697816854e+21, -5.410587449851726e+21, 5.393904740095235e+21,
      -4.533054367123241e+21, 3.181470251931304e+21, -1.840698151594781e+21,
      8.622096114618880e+20, -3.186487205259235e+20, 8.939310241964093e+19,
      -1.788315002192243e+19, 2.271972224754086e+18, -1.377242653770284e+17,
      0.000000000000000e-01, -2.859599139140576e+02, 1.263498115660723e+05,
      -1.563762114667972e+07, 9.735262028286150e+08, -3.727077001705718e+10,
      9.724728419064618e+11, -1.841437932011095e+13, 2.640184847320257e+14,
      -2.955001754223333e+15, 2.641501188215526e+16, -1.919333187312935e+17,
      1.149300695265413e+18, -5.733477927220639e+18, 2.403400433666191e+19,
      -8.522440623029794e+19, 2.569470483161848e+20, -6.610937005531425e+20,
      1.454972285683492e+21, -2.742254919882459e+21, 4.425523119288710e+21,
      -6.106476088763322e+21, 7.183124069502222e+21, -7.170000472640595e+21,
      6.032425263760434e+21, -4.238001893187230e+21, 2.454158919196046e+21,
      -1.150483843370953e+21, 4.254938488755735e+20, -1.194452209653931e+20,
      2.390927098916194e+19, -3.039204212011252e+18, 1.843238572103207e+17,
      0.000000000000000e-01, 1.634368826450262e+02, -7.956978254181237e+04,
      1.188506225092293e+07, -8.373897346830345e+08, 3.478333876388285e+10,
      -9.598393772867426e+11, 1.891593013667364e+13, -2.793128536050818e+14,
      3.196655214719607e+15, -2.907291848273606e+16, 2.141468780998689e+17,
      -1.296440193900301e+18, 6.525499955570070e+18, -2.755634128784674e+19,
      9.831738858873329e+19, -2.979627701541287e+20, 7.700118175018154e+20,
      -1.701110443974214e+21, 3.216663602425906e+21, -5.205922417479511e+21,
      7.201199888459469e+21, -8.489441777094031e+21, 8.490373918346112e+21,
      -7.155631502801035e+21, 5.034836798969603e+21, -2.919618180375494e+21,
      1.370386674687711e+21, -5.073912120857659e+20, 1.425804015450344e+20,
      -2.856660788223791e+19, 3.634277715764878e+18, -2.205841595819251e+17,
      0.000000000000000e-01, -1.089624682152292e+02, 5.490286557349962e+04,
      -8.781653856365252e+06, 6.724119974793967e+08, -2.993574836604601e+10,
      8.717420300323828e+11, -1.790617757376325e+13, 2.730388108651176e+14,
      -3.204823752640508e+15, 2.974036495407882e+16, -2.226577421064578e+17,
      1.366028725355318e+18, -6.951897867287077e+18, 2.962837694696566e+19,
      -1.065339992825695e+20, 3.250038142137898e+20, -8.446629720680499e+20,
      1.875178351120271e+21, -3.560917379168056e+21, 5.784532874816462e+21,
      -8.027770029015426e+21, 9.491255762425382e+21, -9.516676697643140e+21,
      8.038962462344214e+21, -5.667965587919967e+21, 3.292818343972966e+21,
      -1.548126823651629e+21, 5.740631512413311e+20, -1.615361972833789e+20,
      3.240477969405117e+19, -4.127262310667621e+18, 2.507669351857205e+17,
      0.000000000000000e-01, 7.924220268285228e+01, -4.057928795465258e+04,
      6.704332297502558e+06, -5.364604956918865e+08, 2.503646198522573e+10,
      -7.610195551213778e+11, 1.621371476972958e+13, -2.548664517300070e+14,
      3.067722734616486e+15, -2.906734241270257e+16, 2.214249975700081e+17,
      -1.378338812356282e+18, 7.101001568235638e+18, -3.058038376939241e+19,
      1.109399051188740e+20, -3.410471802537065e+20, 8.922579803599067e+20,
      -1.992320553154499e+21, 3.802564715982721e+21, -6.204660474711801e+21,
      8.644825045503473e+21, -1.025665248235139e+22, 1.031630159337718e+22,
      -8.738843624083011e+21, 6.176944575139345e+21, -3.596664164242111e+21,
      1.694457654020145e+21, -6.294971162057861e+20, 1.774358410868178e+20,
      -3.564953335849152e+19, 4.546981308698057e+18, -2.766285114923478e+17,
      0.000000000000000e-01, -6.094810716099949e+01, 3.149084615481936e+04,
      -5.296674502372112e+06, 4.346997745238327e+08, -2.090081537881765e+10,
      6.551264846374850e+11, -1.436792716105805e+13, 2.318076448212311e+14,
      -2.854539793592080e+15, 2.758693070051896e+16, -2.137570364709653e+17,
      1.350278421513405e+18, -7.045142142783022e+18, 3.067459655658656e+19,
      -1.123483766643411e+20, 3.482651662900422e+20, -9.178174833360337e+20,
      2.062604456625362e+21, -3.959135109240898e+21, 6.492786887393485e+21,
      -9.086987198956428e+21, 1.082464243206250e+22, -1.092690369776145e+22,
      9.286162259274268e+21, -6.583075803335110e+21, 3.843334757170004e+21,
      -1.815042549975247e+21, 6.757786079688690e+20, -1.908640091400089e+20,
      3.841802724253007e+19, -4.908370532055551e+18, 2.990786503802649e+17,
      0.000000000000000e-01, 4.874762271398016e+01, -2.532459924379460e+04,
      4.306289111489640e+06, -3.590409832532830e+08, 1.760136771205027e+10,
      -5.636351004936871e+11, 1.263327395525212e+13, -2.081295681322519e+14,
      2.613155513828582e+15, -2.570230208401814e+16, 2.023153786109967e+17,
      -1.296022489490302e+18, 6.846470251956364e+18, -3.013872322115386e+19,
      1.114644422978128e+20, -3.485183237291328e+20, 9.255531281005032e+20,
      -2.094244828319674e+21, 4.044473148145946e+21, -6.669094801072411e+21,
      9.379689895649970e+21, -1.122286017421572e+22, 1.137425276024759e+22,
      -9.701397530693913e+21, 6.900090712824613e+21, -4.040491989417283e+21,
      1.913372195656766e+21, -7.141717651945982e+20, 2.021704584417841e+20,
      -4.077964478229480e+19, 5.220210907153068e+18, -3.186495777814819e+17,
      0.000000000000000e-01, -4.013085745911482e+01, 2.092265556924549e+04,
      -3.583307397906160e+06, 3.019036845311046e+08, -1.499682567334736e+10,
      4.875419883971985e+11, -1.110534210962378e+13, 1.859662144990967e+14,
      -2.372232833495577e+15, 2.368570555021916e+16, -1.890606460337203e+17,
      1.226710982466906e+18, -6.556236864370271e+18, 2.916718348444145e+19,
      -1.089043462231537e+20, 3.434548724554294e+20, -9.192120758302752e+20,
      2.094521354663323e+21, -4.070706954543361e+21, 6.750946202167051e+21,
      -9.544297656724117e+21, 1.147387384551836e+22, -1.167874642179135e+22,
      1.000023491239336e+22, -7.138163819152158e+21, 4.193633752898508e+21,
      -1.991877420085289e+21, 7.455334138647318e+20, -2.115867266165027e+20,
      4.277941431548355e+19, -5.488110031116298e+18, 3.356769581362759e+17,
      0.000000000000000e-01, 3.377658303742900e+01, -1.765321539233928e+04,
      3.038340440514934e+06, -2.578584606807322e+08, 1.292884612296167e+10,
      -4.249332478575872e+11, 9.796453340445514e+12, -1.661322016030880e+14,
      2.146412285165976e+15, -2.170063067679739e+16, 1.753071524663780e+17,
      -1.150445213494342e+18, 6.214123504683469e+18, -2.791805131931260e+19,
      1.051885109311114e+20, -3.345072958360405e+20, 9.021185531229185e+20,
      -2.069976464960887e+21, 4.048800325514122e+21, -6.754020775002269e+21,
      9.599959543792231e+21, -1.159763408249418e+22, 1.185806715330365e+22,
      -1.019594067425241e+22, 7.305655578610319e+21, -4.307135697275206e+21,
      2.052428738060058e+21, -7.705007825196310e+20, 2.192794049399991e+20,
      -4.444874922276286e+19, 5.715873383290923e+18, -3.503832522669480e+17,
      0.000000000000000e-01, -2.892964123166297e+01, 1.514678464159201e+04,
      -2.616230196259168e+06, 2.232056370284014e+08, -1.126780271399590e+10,
      3.733563821201946e+11, -8.686293054983950e+12, 1.487584698214981e+14,
      -1.941627928264148e+15, 1.983312707085937e+16, -1.618550810194144e+17,
      1.072675099399285e+18, -5.848903282733471e+18, 2.651290247911275e+19,
      -1.007365700808825e+20, 3.228755039109242e+20, -8.771430385090572e+20,
      2.026394505724555e+21, -3.988616003559215e+21, 6.692583651440911e+21,
      -9.564190176741721e+21, 1.161237649908247e+22, -1.192826855978758e+22,
      1.030040441263462e+22, -7.409917359938946e+21, 4.384750316004612e+21,
      -2.096582122064279e+21, 7.895856427194460e+20, -2.253772189916786e+20,
      4.581095810605600e+19, -5.906211978906061e+18, 3.629208802827826e+17,
      0.000000000000000e-01, 2.513021588285159e+01, -1.317482888832123e+04,
      2.281636380391524e+06, -1.954241069277852e+08, 9.915879553852877e+09,
      -3.305907224268159e+11, 7.745610540950591e+12, -1.336744677758923e+14,
      1.759053047347007e+15, -1.812022604295268e+16, 1.491398073535739e+17,
      -9.967823580408796e+17, 5.480122870501913e+18, -2.504020332812559e+19,
      9.587106380784488e+19, -3.095239755572209e+20, 8.466795723536765e+20,
      -1.968748626171725e+21, 3.898845137794306e+21, -6.579450581556874e+21,
      9.452948974562387e+21, -1.153486704390072e+22, 1.190416296509360e+22,
      -1.032457212621363e+22, 7.457660044656362e+21, -4.429851231444799e+21,
      2.125705105699995e+21, -8.032242691112874e+20, 2.299857495593084e+20,
      -4.688429048988357e+19, 6.061137852550284e+18, -3.733965028540281e+17,
      0.000000000000000e-01, -2.208386882313948e+01,  1.158936145601140e+04,
      -2.011104324155879e+06, 1.727696977814699e+08, -8.800873732349352e+09,
      2.948204525687799e+11, -6.945679729271358e+12, 1.206045724118162e+14,
      -1.597549338531331e+15, 1.657073913780072e+16, -1.373597907391489e+17,
      9.246697496818755e+17, -5.120171153308929e+18, 2.356084529994870e+19,
      -9.082843603687474e+19, 2.951965214544448e+20, -8.126535950009281e+20,
      1.901181843699431e+21, -3.786945462658442e+21, 6.425892713139728e+21,
      -9.280555773871743e+21, 1.138038464790677e+22, -1.179939706531010e+22,
      1.027858871223671e+22, -7.455106011504378e+21, 4.445547883990988e+21,
      -2.141041178542259e+21, 8.118044630932303e+20, -2.331957228005324e+20,
      4.768370096691469e+19, -6.182197530650295e+18, 3.818855272205113e+17,
      0.000000000000000e-01, 1.959414243804128e+01, -1.029081802559516e+04,
      1.788568172854703e+06, -1.540118150914350e+08, 7.869521604579143e+09,
      -2.646147373681910e+11, 6.261419532268355e+12, -1.092584911616051e+14,
      1.455025808801897e+15, -1.517863644621452e+16, 1.265704706962783e+17,
      -8.572524657708660e+17, 4.776247526570586e+18, -2.211426166036566e+19,
      8.577380339701724e+19, -2.804435586124510e+20, 7.765584710933440e+20,
      -1.827036131219099e+21, 3.659136530638906e+21, -6.241580538370836e+21,
      9.059595877078572e+21, -1.116262879266932e+22, 1.162640401193756e+22,
      -1.017180635719499e+22, 7.408032173818105e+21, -4.434731590925861e+21,
      2.143740476291448e+21, -8.156798285468939e+20, 2.350876961959084e+20,
      -4.822189338581046e+19, 6.270614615621034e+18, -3.884411477495515e+17,
      0.000000000000000e-01, -1.752536768073206e+01, 9.210004044516183e+03,
      -1.602710362254713e+06, 1.382643168640165e+08, -7.082209191374808e+09,
      2.388593247336584e+11, -5.671955040601066e+12, 9.936819687224059e+13,
      -1.329133616843959e+15, 1.393095353694366e+16, -1.167468019716584e+17,
      7.948230830253194e+17, -4.451987054929648e+18, 2.072406047410437e+19,
      -8.081630829587058e+19, 2.656550395781670e+20, -7.395103797274482e+20,
      1.748920815324121e+21, -3.520456085792643e+21, 6.034597075619834e+21,
      -8.800877008017673e+21, 1.089363982096666e+22, -1.139632655550776e+22,
      1.001274028996965e+22, -7.321760958016250e+21, 4.400085249649023e+21,
      -2.134871514399776e+21, 8.151766211905410e+20, -2.357346041272461e+20,
      4.850993390791870e+19, -6.327377892472496e+18, 3.931001229299567e+17,
      0.000000000000000e-01, 1.578106132320405e+01, -8.297465341460668e+03,
      1.445356637015484e+06, -1.248763055461196e+08, 6.409121288934440e+09,
      -2.166867323836831e+11, 5.160255428324136e+12, -9.069980923376013e+13,
      1.217593153323901e+15, -1.281217702947176e+16, 1.078222149705347e+17,
      -7.373025946691738e+17, 4.148685852384210e+18, -1.940267191328230e+19,
      7.602301795140850e+19, -2.510934651183994e+20, 7.023105241061891e+20,
      -1.668804442314510e+21, 3.374862747204858e+21, -5.811516441898685e+21,
      8.513453677248631e+21, -1.058376696971405e+22, 1.111895644585412e+22,
      -9.809014479986894e+21, 7.201130357372178e+21, -4.344074707924148e+21,
      2.115422234515678e+21, -8.105961395642809e+20, 2.352029716821803e+20,
      -4.855759087275033e+19, 6.353294711027987e+18, -3.958864781209368e+17,
      0.000000000000000e-01, -1.429082487951404e+01, 7.516973886624383e+03,
      -1.310468641451437e+06, 1.133605392742571e+08, -5.827511997735239e+09,
      1.974178249055285e+11, -4.712515373341917e+12, 8.305450385112180e+13,
      -1.118328972822835e+15, 1.180653064142852e+16, -9.971166758181266e+16,
      6.844038883665075e+17, -3.866168854945464e+18, 1.815490936983781e+19,
      -7.143043815405257e+19, 2.369235292422868e+20, -6.655061581666021e+20,
      1.588114879269808e+21, -3.225364968659972e+21, 5.577529629049808e+21,
      -8.204710901105034e+21, 1.024169036500672e+22, -1.080270746628453e+22,
      9.567317871713345e+21, -7.050459812170524e+21, 4.268932026970228e+21,
      -2.086295163199683e+21, 8.022143863884770e+20, -2.335532639673231e+20,
      4.837349156604751e+19, -6.349020768043736e+18, 3.968137980027335e+17,
      0.000000000000000e-01, 1.300209931372656e+01, -6.841393840416191e+03,
      1.193492661387020e+06, -1.033456200764548e+08, 5.319778623982719e+09,
      -1.805161947596514e+11, 4.317533185036128e+12, -7.626509475155619e+13,
      1.029508945956098e+15, -1.089906772387182e+16, 9.232461935811958e+16,
      -6.357336264452389e+17, 3.603376239255576e+18, -1.698055636889742e+19,
      6.705347306221760e+19, -2.232368248073979e+20, 6.294452022548759e+20,
      -1.507836188614350e+21, 3.074157987189278e+21, -5.336595912526426e+21,
      7.880489249366001e+21, -9.874488341898865e+21, 1.045462363245188e+22,
      -9.293378659359744e+21, 6.873520006757566e+21, -4.176636208226636e+21,
      2.048299449271481e+21, -7.902800175855315e+20, 2.308396453521622e+20,
      -4.796514797886738e+19, 6.315072588841990e+18, -3.958864781209368e+17,
      0.000000000000000e-01, -1.187480374787823e+01, 6.249975469245804e+03,
      -1.090926975816432e+06, 9.454341715249161e+07, -4.872094426264295e+09,
      1.655534657183789e+11, -3.966168303111664e+12, 7.019129535830873e+13,
      -9.495382387492220e+14, 1.007610862599221e+16, -8.557186874129430e+16,
      5.908530248892806e+17, -3.358744210862016e+18, 1.587616779782261e+19,
      -6.289207145157216e+19, 2.100713184820985e+20, -5.943219992524428e+20,
      1.428595119033459e+21, -2.922754970902040e+21, 5.091600520853282e+21,
      -7.545231007708337e+21, 9.487734903251288e+21, -1.008041543577854e+22,
      8.991956191593799e+21, -6.673509171379833e+21, 4.068893946683280e+21,
      -2.002141237919232e+21, 7.750111453424084e+20, -2.271093028431890e+20,
      4.733888021452981e+19, -6.251826041286116e+18, 3.931001229299567e+17,
      0.000000000000000e-01, 1.087774446529347e+01, -5.726532522119743e+03,
      1.000026921141720e+06, -8.672640379209894e+07, 4.473426627332310e+09,
      -1.521830785057648e+11, 3.650893458666403e+12, -6.471493237339299e+13,
      8.770338012270925e+14, -9.325329559951301e+15, 7.936874744339234e+16,
      -5.493120646406107e+17, 3.130441935246294e+18, -1.483627518494573e+19,
      5.893595494590624e+19, -1.974260043524307e+20, 5.602136541583943e+20,
      -1.350733620606971e+21, 2.772103374812838e+21, -4.844503527518292e+21,
      7.202129004924270e+21, -9.085610385118606e+21, 9.684512731454793e+21,
      -8.666845324741530e+21, 6.453034816508831e+21, -3.947120866745059e+21,
      1.948412902571092e+21, -7.565912375504262e+20, 2.224014019524002e+20,
      -4.649964958533596e+19, 6.159502112364614e+18, -3.884411477495515e+17,
      0.000000000000000e-01, -9.986140223631332e+00, 5.258180249016786e+03,
      -9.185985927746997e+05, 7.971153565650261e+07, -4.114819554974130e+09,
      1.401203840137855e+11, -3.365432249133715e+12, 5.973558126785002e+13,
      -8.107918476058353e+14, 8.635671857055418e+15, -7.363618322103649e+16,
      5.106667921634938e+17,  -2.916510065095704e+18, 1.385415474244305e+19,
      -5.516783142970931e+19, 1.852714215024819e+20, -5.271074464383260e+20,
      1.274366202537370e+21, -2.622681385943975e+21, 4.596469304424025e+21,
      -6.853262795168416e+21, 8.671008970964122e+21, -9.270120999117532e+21,
      8.320885075782802e+21, -6.214097196642535e+21, 3.812422050430118e+21,
      -1.887580884371838e+21, 7.351640810621922e+20, -2.167456694682573e+20,
      4.545079901812915e+19, -6.038139340406067e+18, 3.818855272205113e+17,
      0.000000000000000e-01, 9.179865311132163e+00, -4.834435688158695e+03,
      8.448504512588955e+05, -7.334848338052211e+07, 3.788859743400152e+09,
      -1.291272970868570e+11, 3.104465491388026e+12, -5.516672358352230e+13,
      7.497538905042167e+14, -7.997160527339879e+15, 6.830050930752028e+16,
      -4.744858033268795e+17, 2.714931990926294e+18, -1.292227727069609e+19,
      5.156543365822560e+19, -1.735567332289304e+20, 4.949202035461436e+20,
      -1.199422235955881e+21, 2.474571823687724e+21, -4.347969142753134e+21,
      6.499709691677232e+21, -8.245627749154246e+21, 8.839266680637617e+21,
      -7.955961175047700e+21, 5.958068545674167e+21, -3.665569136784129e+21,
      1.819971125502232e+21, -7.108274904080239e+20, 2.101605178572964e+20,
      -4.419368247642274e+19, 5.887550238778617e+18, -3.733965028540281e+17,
      0.000000000000000e-01, -8.442157387023658e+00, 4.446553379007233e+03,
      -7.772828492878541e+05, 6.751075382325961e+07, -3.489264209336539e+09,
      1.190001385182876e+11, -2.863388547750731e+12, 5.093235749821151e+13,
      -6.929732092529343e+14, 7.400674318693824e+15, -6.329249553366740e+16,
      4.403495020338399e+17, -2.523657530207839e+18, 1.203252080967825e+19,
      -4.810263202989546e+19, 1.622139284826088e+20, -4.635104707077193e+20,
      1.125673633860567e+21, -2.327511860731306e+21, 4.098851177315484e+21,
      -6.141619421494262e+21, 7.810022020551529e+21, -8.392815674436741e+21,
      7.572989471395764e+21, -5.685659534933809e+21, 3.506969466398297e+21,
      -1.745750145592470e+21, 6.836250778812437e+20, -2.026505202012846e+20,
      4.272714338022826e+19, -5.707256190142981e+18, 3.629208802827826e+17,
      0.000000000000000e-01, 7.758620432176431e+00, -4.087010780643976e+03,
      7.146017483117621e+05, -6.208866298624787e+07, 3.210548205967185e+09,
      -1.095595506626254e+11, 2.638102073062174e+12, -4.696390598303079e+13,
      6.395814992099562e+14, -6.837680401038084e+15, 5.854580721511378e+16,
      -4.078439182746133e+17, 2.340589674143266e+18, -1.117619206099512e+19,
      4.474976818167516e+19, -1.511594767319364e+20, 4.326839217130804e+20,
      -1.052747833818875e+21, 2.180916376866737e+21, -3.848371023940319e+21,
      5.778239784350781e+21, -7.363608932672221e+21, 7.930445122237345e+21,
      -7.171862635842173e+21, 5.396860454297140e+21, -3.336620071079666e+21,
      1.664898414527208e+21, -6.535348447882521e+20, 1.942028797566696e+20,
      -4.104676746515045e+19, 5.496390689251412e+18, -3.503832522669480e+17,
      0.000000000000000e-01,  -7.116397718865410e+00, 3.749079648317452e+03,
      -6.556462179536947e+05, 5.698334876317732e+07, -2.947736407946278e+09,
      1.006414850064029e+11, -2.424817814258759e+12, 4.319719539696606e+13,
      -5.887538442457942e+14, 6.299924892323114e+15, -5.399488828717867e+16,
      3.765493816021260e+17, -2.163536905559044e+18, 1.034386803998833e+19,
      -4.147323866260055e+19, 1.402934443266815e+20, -4.021917192563760e+20,
      9.801245774469403e+20, -2.033870287254597e+21, 3.595172619362371e+21,
      -5.407874926278053e+21, 6.904593808205439e+21, -7.450539640286401e+21,
      6.751333782133479e+21, -5.090837485270753e+21, 3.154034661016916e+21,
      -1.577170275266693e+21, 6.204523939342194e+20, -1.847822507348537e+20,
      3.914377458647115e+19, -5.253552629244530e+18, 3.356769581362759e+17,
      0.000000000000000e-01, 6.503405049947615e+00, -3.426426844095358e+03,
      5.993202798185401e+05, -5.210105929407407e+07, 2.696081626301935e+09,
      -9.208817752140222e+10, 2.219856964226088e+12, -3.956916806058024e+13,
      5.396682472814272e+14, -5.779046608150560e+15, 4.957204191401071e+16,
      -3.460227309661112e+17, 1.990124376198793e+18, -9.525027063011665e+18,
      3.823419916527245e+19, -1.294955870987350e+20, 3.717202435921960e+20,
      -9.071121012468153e+20, 1.885079625280803e+21, -3.337199379905083e+21,
      5.027744029565647e+21, -6.429775852763725e+21, 6.949963688436526e+21,
      -6.308792489007701e+21, 4.765750352591050e+21, -2.958122806424166e+21,
      1.482030100296912e+21, -5.841647400501416e+20, 1.743227189970331e+20,
      -3.700329724016468e+19, 4.976575581854351e+18, -3.186495777814819e+17,
      0.000000000000000e-01, -5.907497514106532e+00, 3.112678595829013e+03,
      -5.445178292388210e+05, 4.734677219183480e+07, -2.450744429063559e+09,
      8.373760838189771e+10, -2.019407140707212e+12, 3.601376581781814e+13,
      -4.914525865076907e+14, 5.266042927583444e+15, -4.520313655093684e+16,
      3.157692701941383e+17, -1.817642904902012e+18, 8.707370090912989e+18,
      -3.498599158499494e+19, 1.186170493957938e+20, -3.408681448268105e+20,
      8.327925173944575e+20, -1.732759335262471e+21, 3.071495239478500e+21,
      -4.633677555133575e+21, 5.934150774303194e+21, -6.423619232920034e+21,
      5.839849717449368e+21, -4.418427829351906e+21, 2.746981777317360e+21,
      -1.378544874829898e+21, 5.443069194938188e+20, -1.627146166161762e+20,
      3.460155133741941e+19, -4.662146280112927e+18, 2.990786503802649e+17,
      0.000000000000000e-01, 5.315369270978286e+00, -2.800843064095075e+03,
      4.900224139921877e+05, -4.261558203433840e+07, 2.206354210150574e+09,
      -7.540878769771313e+10, 1.819175365822707e+12, -3.245589496559513e+13,
      4.431044900769581e+14, -4.750435903707279e+15, 4.080066038394767e+16,
      -2.851956961225284e+17, 1.642785899790599e+18, -7.875595013484737e+18,
      3.166934141851351e+19, -1.074644294009643e+20, 3.091013384752239e+20,
      -7.559132271023600e+20, 1.574408999942342e+21, -2.793807988826481e+21,
      4.219517259473626e+21, -5.410137053845813e+21, 5.863599367756335e+21,
      -5.337558212358710e+21, 4.043769184497504e+21, -2.517518733791446e+21,
      1.265191806068099e+21, -5.002850262989436e+20, 1.497812681253764e+20,
      -3.190085448905626e+19, 4.305131059057072e+18, -2.766285114923478e+17,
      0.000000000000000e-01, -4.710772351834030e+00, 2.482373368683492e+03,
      -4.343438634071739e+05, 3.777856440947868e+07, -1.956282178028853e+09,
      6.687710881826573e+10, -1.613799071688983e+12, 2.880102840436904e+13,
      -3.933509953446171e+14, 4.218785748675981e+15, -3.625111116519231e+16,
      2.535230396849270e+17, -1.461153893751222e+18, 7.009048271913450e+18,
      -2.820301208876073e+19, 9.576835312471221e+19, -2.756632919785344e+20,
      6.746687832309763e+20, -1.406360250100074e+21, 2.497787658823857e+21,
      -3.775905442111079e+21, 4.846020652107188e+21, -5.257496778824818e+21,
      4.790865278530114e+21, -3.633565158982738e+21, 2.264711978580421e+21,
      -1.139484606704540e+21, 4.511274997631573e+20, -1.352342175988863e+20,
      2.884004791547230e+19, -3.897279615275436e+18, 2.507669351857205e+17,
      0.000000000000000e-01, 4.070989770987400e+00, -2.145310206509185e+03,
      3.753936962817004e+05, -3.265459454227024e+07, 1.691185508550510e+09,
      -5.782472303217436e+10, 1.395652621524724e+12, -2.491398584416296e+13,
      3.403599167085817e+14, -3.651608452505086e+15, 3.138862675548854e+16,
      -2.196030665561983e+17, 1.266200972268635e+18, -6.076691812326445e+18,
      2.446363025550984e+19, -8.311520079137181e+19, 2.393790730533568e+20,
      -5.862224225217512e+20, 1.222781062900224e+21, -2.173219858365351e+21,
      3.287615107047680e+21, -4.222527251629046e+21, 4.584681175928391e+21,
      -4.181215238485294e+21, 3.173915831158933e+21, -1.979997675938663e+21,
      9.971596317435381e+20, -3.951620422561584e+20, 1.185761286177830e+20,
      -2.531374184616153e+19, 3.424415390856724e+18, -2.205841595819251e+17,
      0.000000000000000e-01, -3.358090451409268e+00, 1.769674233094550e+03,
      -3.096791496404924e+05, 2.694027470522958e+07, -1.395380783042074e+09,
      4.771664530498290e+10, -1.151859937920620e+12, 2.056566436202355e+13,
      -2.810129791321012e+14, 3.015587336970416e+15, -2.592812452685370e+16,
      1.814511218788989e+17, -1.046544742872117e+18, 5.024209499496014e+18,
      -2.023384009065594e+19, 6.877115067771995e+19, -1.981490581745830e+20,
      4.854671646250063e+20, -1.013093139325803e+21, 1.801437605201838e+21,
      -2.726610427990654e+21, 3.503909183048419e+21, -3.806619619529808e+21,
      3.473717880327268e+21, -2.638522642188920e+21, 1.647082019674793e+21,
      -8.300649678444784e+20, 3.291782934571713e+20, -9.884928188742345e+19,
      2.111857359313218e+19, -2.859159218719010e+18, 1.843238572103207e+17,
      0.000000000000000e-01, 2.488636218117664e+00, -1.311502616747683e+03,
      2.295099147207366e+05, -1.996696431094197e+07, 1.034261165808075e+09,
      -3.537054923192207e+10, 8.539117568460436e+11, -1.524770635014292e+13,
      2.083740755848528e+14, -2.236412246676083e+15, 1.923184048547382e+16,
      -1.346128075294067e+17, 7.765487558375990e+17, -3.728808938617731e+18,
      1.502032959138465e+19, -5.106384693104742e+19, 1.471677691807719e+20,
      -3.606628515829967e+20, 7.528686576360846e+20, -1.339136443296313e+21,
      2.027551807534738e+21, -2.606468672347013e+21, 2.832680005398635e+21,
      -2.585940609411422e+21, 1.964981315222245e+21, -1.227139920687813e+21,
      6.186996825719858e+20, -2.454684425809199e+20, 7.374666999744300e+19,
      -1.576324668155186e+19, 2.135204267310823e+18, -1.377242653770284e+17,
      0.000000000000000e-01, -1.000000000000000e+00, 5.270000000000000e+02,
      -9.222500000000000e+04, 8.023575000000000e+06, -4.156211850000000e+08,
      1.421424452700000e+10, -3.431724750090000e+11, 6.128079910875000e+12,
      -8.375042544862500e+13, 8.989212331485750e+14, -7.730722605077745e+15,
      5.411505823554421e+16, -3.122022590512166e+17, 1.499257002256941e+18,
      -6.039863923377964e+18, 2.053553733948508e+19, -5.919066644910405e+19,
      1.450751628654511e+20, -3.028762172103277e+20, 5.388008495636356e+20,
      -8.158984293392197e+20, 1.049012266293282e+21, -1.140230724231829e+21,
      1.041080226472539e+21, -7.912209721191298e+20, 4.942087918159488e+20,
      -2.492163992918032e+20, 9.889539654436636e+19, -2.971733590742043e+19,
      6.353361469862300e+18, -8.607780055942471e+17, 5.553406487704820e+16 };

  static const double *fem_coeff_gausslob[] =
    { 0, fem_coef_gausslob_1, fem_coef_gausslob_2, fem_coef_gausslob_3,
      fem_coef_gausslob_4, fem_coef_gausslob_5, fem_coef_gausslob_6, fem_coef_gausslob_7, fem_coef_gausslob_8,
      fem_coef_gausslob_9, fem_coef_gausslob_10, fem_coef_gausslob_11, fem_coef_gausslob_12, fem_coef_gausslob_13,
      fem_coef_gausslob_14, 0, fem_coef_gausslob_16, 0, 0, 0, 0, 0, 0, 0, fem_coef_gausslob_24, 0, 0, 0, 0, 0, 0, 0,
      fem_coef_gausslob_32, };

  const unsigned fem_coeff_gausslob_max_k = 33;


  class PK_GL_fem_ : public fem<base_poly> {
  public :
    PK_GL_fem_(unsigned k);
  };

  PK_GL_fem_::PK_GL_fem_(unsigned k) {
    cvr = bgeot::simplex_of_reference(1);
    dim_ = cvr->structure()->dim();
    is_standard_fem = is_equiv = is_pol = is_lag = true;
    es_degree = short_type(k);
    GMM_ASSERT1(k < fem_coeff_gausslob_max_k && fem_coeff_gausslob[k],
                "try another degree");
    init_cvs_node();
    std::stringstream sstr; sstr << "IM_GAUSSLOBATTO1D(" << k*2-1 << ")";
    pintegration_method gl_im = int_method_descriptor(sstr.str());
    std::vector<base_node> points(k+1);
    for (size_type i = 0; i < k+1; ++i) {
      points[i] = gl_im->approx_method()->point(i);
    }
    std::sort(points.begin(),points.end());
    for (size_type i = 0; i < k+1; ++i) {
      // cout << points[i][0] << endl;
      add_node(lagrange_dof(1), points[i]);
    }
    base_.resize(k+1);
    const double *coefs = fem_coeff_gausslob[k];
    for (size_type r = 0; r < k+1; r++) {
      base_[r] = base_poly(1,short_type(k));
      std::copy(coefs, coefs+k+1, base_[r].begin());
      coefs += k+1;
    }
  }

  static pfem PK_GL_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 1, "Bad number of parameters : "
                << params.size() << " should be 1.");
    GMM_ASSERT1(params[0].type() == 0, "Bad type of parameters");
    int k = int(::floor(params[0].num() + 0.01));
    pfem p = std::make_shared<PK_GL_fem_>(k);
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*        Hermite element on the segment                                    */
  /* ******************************************************************** */

  struct hermite_segment__ : public fem<base_poly> {
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
                           bgeot::pgeometric_trans pgt) const;
    hermite_segment__();
  };

  void hermite_segment__::mat_trans(base_matrix &M,
                                    const base_matrix &G,
                                    bgeot::pgeometric_trans pgt) const {
    THREAD_SAFE_STATIC bgeot::pgeotrans_precomp pgp;
    THREAD_SAFE_STATIC bgeot::pgeometric_trans pgt_stored = nullptr;
    THREAD_SAFE_STATIC base_matrix K(1, 1);
    THREAD_SAFE_STATIC base_vector r(1);
    dim_type N = dim_type(G.nrows());

    if (pgt != pgt_stored) {
      gmm::resize(r, N);
      for (size_type i = 0; i < N; ++i) r[i] = ::exp(double(i));
      pgt_stored = pgt;
      pgp = bgeot::geotrans_precomp(pgt, node_tab(0), 0);
      gmm::resize(K, N, 1);
    }
    gmm::copy(gmm::identity_matrix(), M);
    // gradient at point 0
    gmm::mult(G, pgp->grad(1), K);
    if (N == 1) M(1, 1) = K(0,0);
    else  M(1, 1) = gmm::mat_euclidean_norm(K)
      * gmm::sgn(gmm::vect_sp(gmm::mat_col(K, 0), r));
    // gradient at point 1
    if (!(pgt->is_linear())) gmm::mult(G, pgp->grad(3), K);
    if (N == 1) M(3, 3) = K(0,0);
    else M(3, 3) = gmm::mat_euclidean_norm(K)
      * gmm::sgn(gmm::vect_sp(gmm::mat_col(K, 0), r));
  }

  // Hermite element on the segment. when the real element lies in
  // a 2 or 3 dimensional domain, the element should still work if
  // the tangent coincides.
  hermite_segment__::hermite_segment__() {
    base_node pt(1);
    cvr = bgeot::simplex_of_reference(1);
    dim_ = cvr->structure()->dim();
    init_cvs_node();
    es_degree = 3;
    is_pol = true;
    is_standard_fem = is_lag = is_equiv = false;
    base_.resize(4);

    pt[0] = 0.0; add_node(lagrange_dof(1), pt);
    read_poly(base_[0], 1, "(1 - x)^2*(2*x + 1)");

    pt[0] = 0.0; add_node(derivative_dof(1, 0), pt);
    read_poly(base_[1], 1, "x*(x - 1)*(x - 1)");

    pt[0] = 1.0; add_node(lagrange_dof(1), pt);
    read_poly(base_[2], 1, "x*x*(3  - 2*x)");

    pt[0] = 1.0; add_node(derivative_dof(1, 0), pt);
    read_poly(base_[3], 1, "x*x*(x - 1)");
  }

  /* ******************************************************************** */
  /*    Hermite element on the triangle                                   */
  /* ******************************************************************** */

  struct hermite_triangle__ : public fem<base_poly> {
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
                           bgeot::pgeometric_trans pgt) const;
    hermite_triangle__();
  };

  void hermite_triangle__::mat_trans(base_matrix &M,
                                    const base_matrix &G,
                                    bgeot::pgeometric_trans pgt) const {

    THREAD_SAFE_STATIC bgeot::pgeotrans_precomp pgp;
    THREAD_SAFE_STATIC bgeot::pgeometric_trans pgt_stored = nullptr;
    THREAD_SAFE_STATIC base_matrix K(2, 2);
    dim_type N = dim_type(G.nrows());

    GMM_ASSERT1(N == 2, "Sorry, this version of hermite "
                "element works only in dimension two.")
    if (pgt != pgt_stored)
      { pgt_stored = pgt; pgp = bgeot::geotrans_precomp(pgt, node_tab(0), 0); }
    gmm::copy(gmm::identity_matrix(), M);

    gmm::mult(G, pgp->grad(0), K);
    for (size_type i = 0; i < 3; ++i) {
      if (i && !(pgt->is_linear())) gmm::mult(G, pgp->grad(i*3), K);
      gmm::copy(K, gmm::sub_matrix(M, gmm::sub_interval(1+3*i, 2)));
    }
  }

  hermite_triangle__::hermite_triangle__() {
    cvr = bgeot::simplex_of_reference(2);
    dim_ = cvr->structure()->dim();
    init_cvs_node();
    es_degree = 3;
    is_pol = true;
    is_standard_fem = is_lag = is_equiv = false;
    base_.resize(10);

    add_node(lagrange_dof(2), base_node(0.0, 0.0));
    read_poly(base_[0], 2, "(1 - x - y)*(1 + x + y - 2*x*x - 11*x*y - 2*y*y)");

    add_node(derivative_dof(2, 0), base_node(0.0, 0.0));
    read_poly(base_[1], 2, "x*(1 - x - y)*(1 - x - 2*y)");

    add_node(derivative_dof(2, 1), base_node(0.0, 0.0));
    read_poly(base_[2], 2, "y*(1 - x - y)*(1 - 2*x - y)");

    add_node(lagrange_dof(2), base_node(1.0, 0.0));
    read_poly(base_[3], 2, "-2*x*x*x + 7*x*x*y + 7*x*y*y + 3*x*x - 7*x*y");

    add_node(derivative_dof(2, 0), base_node(1.0, 0.0));
    read_poly(base_[4], 2, "x*x*x - 2*x*x*y - 2*x*y*y - x*x + 2*x*y");

    add_node(derivative_dof(2, 1), base_node(1.0, 0.0));
    read_poly(base_[5], 2, "x*y*(2*x + y - 1)");

    add_node(lagrange_dof(2), base_node(0.0, 1.0));
    read_poly(base_[6], 2, "7*x*x*y + 7*x*y*y - 2*y*y*y + 3*y*y - 7*x*y");

    add_node(derivative_dof(2, 0), base_node(0.0, 1.0));
    read_poly(base_[7], 2, "x*y*(x + 2*y - 1)");

    add_node(derivative_dof(2, 1), base_node(0.0, 1.0));
    read_poly(base_[8], 2, "y*y*y - 2*y*y*x - 2*y*x*x - y*y + 2*x*y");

    add_node(lagrange_dof(2), base_node(1.0/3.0, 1.0/3.0));
    read_poly(base_[9], 2, "27*x*y*(1 - x - y)");

 }

  /* ******************************************************************** */
  /*    Hermite element on the tetrahedron                                */
  /* ******************************************************************** */

  struct hermite_tetrahedron__ : public fem<base_poly> {
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
                           bgeot::pgeometric_trans pgt) const;
    hermite_tetrahedron__();
  };

  void hermite_tetrahedron__::mat_trans(base_matrix &M,
                                    const base_matrix &G,
                                    bgeot::pgeometric_trans pgt) const {
    THREAD_SAFE_STATIC bgeot::pgeotrans_precomp pgp;
    THREAD_SAFE_STATIC bgeot::pgeometric_trans pgt_stored = nullptr;
    THREAD_SAFE_STATIC base_matrix K(3, 3);
    dim_type N = dim_type(G.nrows());
    GMM_ASSERT1(N == 3, "Sorry, this version of hermite "
                "element works only on dimension three.")
    if (pgt != pgt_stored)
      { pgt_stored = pgt; pgp = bgeot::geotrans_precomp(pgt, node_tab(0), 0); }
    gmm::copy(gmm::identity_matrix(), M);

    gmm::mult(G, pgp->grad(0), K);
    for (size_type k = 0; k < 4; ++k) {
      if (k && !(pgt->is_linear())) gmm::mult(G, pgp->grad(k*4), K);
      gmm::copy(K, gmm::sub_matrix(M, gmm::sub_interval(1+4*k, 3)));
    }
  }

  hermite_tetrahedron__::hermite_tetrahedron__() {
    cvr = bgeot::simplex_of_reference(3);
    dim_ = cvr->structure()->dim();
    init_cvs_node();
    es_degree = 3;
    is_pol = true;
    is_standard_fem = is_lag = is_equiv = false;
    base_.resize(20);
    std::stringstream s
      ( "1 - 3*x*x - 13*x*y - 13*x*z - 3*y*y - 13*y*z - 3*z*z + 2*x*x*x"
        "+ 13*x*x*y + 13*x*x*z + 13*x*y*y + 33*x*y*z + 13*x*z*z + 2*y*y*y"
        "+ 13*y*y*z + 13*y*z*z + 2*z*z*z;"
        "x - 2*x*x - 3*x*y - 3*x*z + x*x*x + 3*x*x*y + 3*x*x*z + 2*x*y*y"
        "+ 4*x*y*z + 2*x*z*z;"
        "y - 3*x*y - 2*y*y - 3*y*z + 2*x*x*y + 3*x*y*y + 4*x*y*z"
        "+ y*y*y + 3*y*y*z + 2*y*z*z;"
        "z - 3*x*z - 3*y*z - 2*z*z + 2*x*x*z + 4*x*y*z + 3*x*z*z"
        "+ 2*y*y*z + 3*y*z*z + z*z*z;"
        "3*x*x - 7*x*y - 7*x*z - 2*x*x*x + 7*x*x*y + 7*x*x*z + 7*x*y*y"
        "+ 7*x*y*z + 7*x*z*z;"
        "-x*x + 2*x*y + 2*x*z + x*x*x - 2*x*x*y - 2*x*x*z - 2*x*y*y"
        "- 2*x*y*z - 2*x*z*z;"
        "-x*y + 2*x*x*y + x*y*y;"
        "-x*z + 2*x*x*z + x*z*z;"
        "-7*x*y + 3*y*y - 7*y*z + 7*x*x*y + 7*x*y*y + 7*x*y*z - 2*y*y*y"
        "+ 7*y*y*z + 7*y*z*z;"
        "-x*y + x*x*y + 2*x*y*y;"
        "2*x*y - y*y + 2*y*z - 2*x*x*y - 2*x*y*y - 2*x*y*z + y*y*y"
        "- 2*y*y*z - 2*y*z*z;"
        "-y*z + 2*y*y*z + y*z*z;"
        "-7*x*z - 7*y*z + 3*z*z + 7*x*x*z + 7*x*y*z + 7*x*z*z + 7*y*y*z"
        "+ 7*y*z*z - 2*z*z*z;"
        "-x*z + x*x*z + 2*x*z*z;"
        "-y*z + y*y*z + 2*y*z*z;"
        "2*x*z + 2*y*z - z*z - 2*x*x*z - 2*x*y*z - 2*x*z*z - 2*y*y*z"
        "- 2*y*z*z + z*z*z;"
        "27*x*y*z;"
        "27*y*z - 27*x*y*z - 27*y*y*z - 27*y*z*z;"
        "27*x*z - 27*x*x*z - 27*x*y*z - 27*x*z*z;"
        "27*x*y - 27*x*x*y - 27*x*y*y - 27*x*y*z;");

    base_node pt(3);
    for (unsigned k = 0; k < 5; ++k) {
      for (unsigned i = 0; i < 4; ++i) {
        base_[k*4+i] = read_base_poly(3, s);
        pt[0] = pt[1] = pt[2] = ((k == 4) ? 1.0/3.0 : 0.0);
        if (k == 4 && i) pt[i-1] = 0.0;
        if (k < 4 && k) pt[k-1] = 1.0;
        if (k == 4 || i == 0)  add_node(lagrange_dof(3), pt);
        else add_node(derivative_dof(3, dim_type(i-1)), pt);
      }
    }
  }

  static pfem Hermite_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 1, "Bad number of parameters : "
                << params.size() << " should be 1.");
    GMM_ASSERT1(params[0].type() == 0, "Bad type of parameters");
    int d = int(::floor(params[0].num() + 0.01));
    pfem p;
    switch(d) {
    case 1 : p = std::make_shared<hermite_segment__>(); break;
    case 2 : p = std::make_shared<hermite_triangle__>(); break;
    case 3 : p = std::make_shared<hermite_tetrahedron__>(); break;
    default : GMM_ASSERT1(false, "Sorry, Hermite element in dimension "
                          << d << " not available");
    }
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return pfem(p);
  }

  /* ******************************************************************** */
  /*    Argyris element on the triangle                                   */
  /* ******************************************************************** */

  struct argyris_triangle__ : public fem<base_poly> {
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
                           bgeot::pgeometric_trans pgt) const;
    argyris_triangle__();
  };

  void argyris_triangle__::mat_trans(base_matrix &M,
                                    const base_matrix &G,
                                    bgeot::pgeometric_trans pgt) const {

    THREAD_SAFE_STATIC bgeot::pgeotrans_precomp pgp;
    THREAD_SAFE_STATIC pfem_precomp pfp;
    THREAD_SAFE_STATIC bgeot::pgeometric_trans pgt_stored = nullptr;
    THREAD_SAFE_STATIC base_matrix K(2, 2);
    dim_type N = dim_type(G.nrows());
    GMM_ASSERT1(N == 2, "Sorry, this version of argyris "
                "element works only on dimension two.")
    if (pgt != pgt_stored) {
      pgt_stored = pgt;
      pgp = bgeot::geotrans_precomp(pgt, node_tab(0), 0);
      pfp = fem_precomp(std::make_shared<argyris_triangle__>(), node_tab(0), 0);
    }
    gmm::copy(gmm::identity_matrix(), M);

    gmm::mult(G, pgp->grad(0), K); // gradient at point (0, 0)
    for (unsigned k = 0; k < 3; ++k) {
      if (k && !(pgt->is_linear())) gmm::mult(G, pgp->grad(6*k), K);
      M(1+6*k, 1+6*k) = K(0,0); M(1+6*k, 2+6*k) = K(0,1);
      M(2+6*k, 1+6*k) = K(1,0); M(2+6*k, 2+6*k) = K(1,1);
      if (!(pgt->is_linear())) {
        base_matrix XX[2], H(2,4), B(2,2), X(2,2);
        XX[0] = XX[1] = base_matrix(2,2);
        gmm::copy(gmm::transposed(K), B); gmm::lu_inverse(B);
        gmm::mult(G, pgp->hessian(6*k), H);
        for (unsigned j = 0; j < 2; ++j) {
          XX[j](0,0) = B(0, j)*H(0, 0) + B(1, j)*H(1, 0);
          XX[j](0,1) = XX[j](1,0) = B(0, j)*H(0, 1) + B(1, j)*H(1, 1);
          XX[j](1,1) = B(0, j)*H(0, 3) + B(1, j)*H(1, 3);
        }
        for (unsigned j = 0; j < 2; ++j) {
          gmm::copy(gmm::scaled(XX[0], K(j,0)), X);
          gmm::add(gmm::scaled(XX[1], K(j,1)), X);
          M(1+j+6*k, 3+6*k) = X(0,0); M(1+j+6*k, 4+6*k) = X(1, 0);
          M(1+j+6*k, 5+6*k) = X(1, 1);
        }
      }
      scalar_type a = K(0,0), b = K(0,1), c = K(1,0), d = K(1,1);
      M(3+6*k, 3+6*k) = a*a;     M(3+6*k, 4+6*k) = a*b;       M(3+6*k, 5+6*k) = b*b;
      M(4+6*k, 3+6*k) = 2.0*a*c; M(4+6*k, 4+6*k) = b*c + a*d; M(4+6*k, 5+6*k) = 2.0*b*d;
      M(5+6*k, 3+6*k) = c*c;     M(5+6*k, 4+6*k) = c*d;       M(5+6*k, 5+6*k) = d*d;
    }

    THREAD_SAFE_STATIC base_matrix W(3, 21);
    base_small_vector norient(M_PI, M_PI * M_PI);
    if (pgt->is_linear()) gmm::lu_inverse(K);
    for (unsigned i = 18; i < 21; ++i) {
      if (!(pgt->is_linear()))
        { gmm::mult(G, pgp->grad(i), K); gmm::lu_inverse(K); }
      bgeot::base_small_vector n(2), v(2);
      gmm::mult(gmm::transposed(K), cvr->normals()[i-18], n);
      n /= gmm::vect_norm2(n);

      scalar_type ps = gmm::vect_sp(n, norient);
      if (ps < 0) n *= scalar_type(-1);
      if (gmm::abs(ps) < 1E-8)
        GMM_WARNING2("Argyris : The normal orientation may be incorrect");
      gmm::mult(K, n, v);
      const bgeot::base_tensor &t = pfp->grad(i);
      for (unsigned j = 0; j < 21; ++j)
        W(i-18, j) = t(j, 0, 0) * v[0] + t(j, 0, 1) * v[1];
    }
    THREAD_SAFE_STATIC base_matrix A(3,3);
    THREAD_SAFE_STATIC bgeot::base_vector w(3);
    THREAD_SAFE_STATIC bgeot::base_vector coeff(3);
    THREAD_SAFE_STATIC gmm::sub_interval SUBI(18, 3);
    THREAD_SAFE_STATIC gmm::sub_interval SUBJ(0, 3);
    gmm::copy(gmm::sub_matrix(W, SUBJ, SUBI), A);
    gmm::lu_inverse(A);
    gmm::copy(gmm::transposed(A), gmm::sub_matrix(M, SUBI));

    for (unsigned j = 0; j < 18; ++j) {
      gmm::mult(W, gmm::mat_row(M, j), w);
      gmm::mult(A, gmm::scaled(w, -1.0), coeff);
      gmm::copy(coeff, gmm::sub_vector(gmm::mat_row(M, j), SUBI));
    }
  }

  argyris_triangle__::argyris_triangle__() {
    cvr = bgeot::simplex_of_reference(2);
    dim_ = cvr->structure()->dim();
    init_cvs_node();
    es_degree = 5;
    is_pol = true;
    is_lag = false;
    is_standard_fem = is_equiv = false;
    base_.resize(21);

    std::stringstream s
      ("1 - 10*x^3 - 10*y^3 + 15*x^4 - 30*x*x*y*y"
       "+ 15*y*y*y*y - 6*x^5 + 30*x*x*x*y*y + 30*x*x*y*y*y - 6*y^5;"
       "x - 6*x*x*x - 11*x*y*y + 8*x*x*x*x + 10*x*x*y*y"
       "+ 18*x*y*y*y - 3*x*x*x*x*x + x*x*x*y*y - 10*x*x*y*y*y - 8*x*y*y*y*y;"
       "y - 11*x*x*y - 6*y*y*y + 18*x*x*x*y + 10*x*x*y*y"
       "+ 8*y*y*y*y - 8*x*x*x*x*y - 10*x*x*x*y*y + x*x*y*y*y - 3*y*y*y*y*y;"
       "0.5*x*x - 1.5*x*x*x + 1.5*x*x*x*x - 1.5*x*x*y*y"
       "- 0.5*x*x*x*x*x + 1.5*x*x*x*y*y + x*x*y*y*y;"
       "x*y - 4*x*x*y - 4*x*y*y + 5*x*x*x*y + 10*x*x*y*y"
       "+ 5*x*y*y*y - 2*x*x*x*x*y - 6*x*x*x*y*y - 6*x*x*y*y*y - 2*x*y*y*y*y;"
       "0.5*y*y - 1.5*y*y*y - 1.5*x*x*y*y + 1.5*y*y*y*y + x*x*x*y*y"
       "+ 1.5*x*x*y*y*y - 0.5*y*y*y*y*y;"
       "10*x^3 - 15*x^4 + 15*x*x*y*y + 6*x^5 - 15*x*x*x*y*y - 15*x*x*y*y*y;"
       "-4*x*x*x + 7*x*x*x*x - 3.5*x*x*y*y - 3*x*x*x*x*x + 3.5*x*x*x*y*y"
       "+ 3.5*x*x*y*y*y;"
       "-5*x*x*y + 14*x*x*x*y + 18.5*x*x*y*y - 8*x*x*x*x*y"
       "- 18.5*x*x*x*y*y - 13.5*x*x*y*y*y;"
       "0.5*x*x*x - x*x*x*x + 0.25*x*x*y*y + 0.5*x*x*x*x*x"
       "- 0.25*x*x*x*y*y - 0.25*x*x*y*y*y;"
       "x*x*y - 3*x*x*x*y - 3.5*x*x*y*y + 2*x*x*x*x*y + 3.5*x*x*x*y*y"
       "+ 2.5*x*x*y*y*y;"
       "1.25*x*x*y*y - 0.75*x*x*x*y*y - 1.25*x*x*y*y*y;"
       "10*y*y*y + 15*x*x*y*y - 15*y^4 - 15*x*x*x*y*y - 15*x*x*y*y*y + 6*y^5;"
       "-5*x*y*y + 18.5*x*x*y*y + 14*x*y*y*y - 13.5*x*x*x*y*y"
       "- 18.5*x*x*y*y*y - 8*x*y*y*y*y;"
       "-4*y*y*y - 3.5*x*x*y*y + 7*y*y*y*y + 3.5*x*x*x*y*y"
       "+ 3.5*x*x*y*y*y - 3*y*y*y*y*y;"
       "1.25*x*x*y*y - 1.25*x*x*x*y*y - 0.75*x*x*y*y*y;"
       "x*y*y - 3.5*x*x*y*y - 3*x*y*y*y + 2.5*x*x*x*y*y + 3.5*x*x*y*y*y"
       "+ 2*x*y*y*y*y;"
       "0.5*y*y*y + 0.25*x*x*y*y - y*y*y*y - 0.25*x*x*x*y*y"
       "- 0.25*x*x*y*y*y + 0.5*y*y*y*y*y;"
       "sqrt(2) * (-8*x*x*y*y + 8*x*x*x*y*y + 8*x*x*y*y*y);"
       "-16*x*y*y + 32*x*x*y*y + 32*x*y*y*y - 16*x*x*x*y*y"
       "- 32*x*x*y*y*y - 16*x*y*y*y*y;"
       "-16*x*x*y + 32*x*x*x*y + 32*x*x*y*y - 16*x*x*x*x*y"
       "- 32*x*x*x*y*y - 16*x*x*y*y*y;");

    base_node pt(2);
    for (unsigned k = 0; k < 7; ++k) {
      for (unsigned i = 0; i < 3; ++i) {
        base_[k*3+i] = read_base_poly(2, s);
        if (k == 6) {
          pt[0] = pt[1] = 0.5; if (i) pt[i-1] = 0.0;
          add_node(normal_derivative_dof(2), pt);
        }
        else {
          pt[0] = pt[1] = 0.0; if (k/2) pt[k/2-1] = 1.0;
          if (k & 1)
            add_node(second_derivative_dof(2, dim_type((i) ? 1:0),
                                           dim_type((i == 2) ? 1:0)), pt);
          else {
            if (i) add_node(derivative_dof(2, dim_type(i-1)), pt);
            else add_node(lagrange_dof(2), pt);
          }
        }
      }
    }
  }

  static pfem triangle_Argyris_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 0, "Bad number of parameters");
    pfem p = std::make_shared<argyris_triangle__>();
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*    Morley element on the triangle                                    */
  /* ******************************************************************** */

  struct morley_triangle__ : public fem<base_poly> {
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
                           bgeot::pgeometric_trans pgt) const;
    morley_triangle__();
  };

  void morley_triangle__::mat_trans(base_matrix &M,
                                    const base_matrix &G,
                                    bgeot::pgeometric_trans pgt) const {

    THREAD_SAFE_STATIC bgeot::pgeotrans_precomp pgp;
    THREAD_SAFE_STATIC pfem_precomp pfp;
    THREAD_SAFE_STATIC bgeot::pgeometric_trans pgt_stored = nullptr;
    THREAD_SAFE_STATIC base_matrix K(2, 2);
    dim_type N = dim_type(G.nrows());
    GMM_ASSERT1(N == 2, "Sorry, this version of morley "
                "element works only on dimension two.")

    if (pgt != pgt_stored) {
      pgt_stored = pgt;
      pgp = bgeot::geotrans_precomp(pgt, node_tab(0), 0);
      pfp = fem_precomp(std::make_shared<morley_triangle__>(), node_tab(0), 0);
    }
    gmm::copy(gmm::identity_matrix(), M);
    THREAD_SAFE_STATIC base_matrix W(3, 6);
    base_small_vector norient(M_PI, M_PI * M_PI);
    if (pgt->is_linear())
      { gmm::mult(G, pgp->grad(0), K); gmm::lu_inverse(K); }
    for (unsigned i = 3; i < 6; ++i) {
      if (!(pgt->is_linear()))
        { gmm::mult(G, pgp->grad(i), K); gmm::lu_inverse(K); }
      bgeot::base_small_vector n(2), v(2);
      gmm::mult(gmm::transposed(K), cvr->normals()[i-3], n);
      n /= gmm::vect_norm2(n);

      scalar_type ps = gmm::vect_sp(n, norient);
      if (ps < 0) n *= scalar_type(-1);
      if (gmm::abs(ps) < 1E-8)
        GMM_WARNING2("Morley : The normal orientation may be incorrect");
      gmm::mult(K, n, v);
      const bgeot::base_tensor &t = pfp->grad(i);
      for (unsigned j = 0; j < 6; ++j)
        W(i-3, j) = t(j, 0, 0) * v[0] + t(j, 0, 1) * v[1];
    }
    //    cout << "W = " << W << endl; getchar();

    THREAD_SAFE_STATIC base_matrix A(3, 3);
    THREAD_SAFE_STATIC base_vector w(3);
    THREAD_SAFE_STATIC base_vector coeff(3);
    THREAD_SAFE_STATIC gmm::sub_interval SUBI(3, 3);
    THREAD_SAFE_STATIC gmm::sub_interval SUBJ(0, 3);
    gmm::copy(gmm::sub_matrix(W, SUBJ, SUBI), A);
    gmm::lu_inverse(A);
    gmm::copy(gmm::transposed(A), gmm::sub_matrix(M, SUBI));

    for (unsigned j = 0; j < 3; ++j) {
      gmm::mult(W, gmm::mat_row(M, j), w);
      gmm::mult(A, gmm::scaled(w, -1.0), coeff);
      gmm::copy(coeff, gmm::sub_vector(gmm::mat_row(M, j), SUBI));
    }
  }

  morley_triangle__::morley_triangle__() {
    cvr = bgeot::simplex_of_reference(2);
    dim_ = cvr->structure()->dim();
    init_cvs_node();
    es_degree = 2;
    is_pol = true;
    is_standard_fem = is_lag = is_equiv = false;
    base_.resize(6);

    std::stringstream s("1 - x - y + 2*x*y;  (x + y + x^2 - 2*x*y - y^2)/2;"
                        "(x + y - x^2 - 2*x*y + y^2)/2;"
                        "((x+y)^2 - x - y)*sqrt(2)/2;  x*(x-1);  y*(y-1);");

    for (unsigned k = 0; k < 6; ++k)
      base_[k] = read_base_poly(2, s);

    add_node(lagrange_dof(2), base_node(0.0, 0.0));
    add_node(lagrange_dof(2), base_node(1.0, 0.0));
    add_node(lagrange_dof(2), base_node(0.0, 1.0));
    add_node(normal_derivative_dof(2), base_node(0.5, 0.5));
    add_node(normal_derivative_dof(2), base_node(0.0, 0.5));
    add_node(normal_derivative_dof(2), base_node(0.5, 0.0));
  }

  static pfem triangle_Morley_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 0, "Bad number of parameters");
    pfem p = std::make_shared<morley_triangle__>();
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*    DISCONTINUOUS PK                                                  */
  /* ******************************************************************** */

  struct PK_discont_ : public PK_fem_ {
  public :

    PK_discont_(dim_type nc, short_type k, scalar_type alpha=scalar_type(0))
      : PK_fem_(nc, k) {

      std::fill(dof_types_.begin(),
                dof_types_.end(), lagrange_nonconforming_dof(nc));
     
      if (alpha != scalar_type(0)) {
        base_node G =
          gmm::mean_value(cv_node.points().begin(), cv_node.points().end());
        for (size_type i=0; i < cv_node.nb_points(); ++i)
          cv_node.points()[i] = (1-alpha)*cv_node.points()[i] + alpha*G;
        for (size_type d = 0; d < nc; ++d) {
          base_poly S(1,2);
          S[0] = -alpha * G[d] / (1-alpha);
          S[1] = 1. / (1-alpha);
          for (size_type j=0; j < nb_base(0); ++j)
            base_[j] = bgeot::poly_substitute_var(base_[j],S,d);
        }
      }
    }
  };

  static pfem PK_discontinuous_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 2 || params.size() == 3,
                "Bad number of parameters : "
                << params.size() << " should be 2 or 3.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0 &&
                (params.size() != 3 || params[2].type() == 0),
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    scalar_type alpha = 0.;
    if (params.size() == 3) alpha = params[2].num();
    GMM_ASSERT1(n > 0 && n < 100 && k >= 0 && k <= 150 && alpha >= 0 &&
                alpha < 1 && double(n) == params[0].num()
                && double(k) == params[1].num(), "Bad parameters");
    pfem p = std::make_shared<PK_discont_>(dim_type(n), short_type(k), alpha);
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*    Connectable interior PK                                           */
  /* ******************************************************************** */

  struct PK_conn_int_ : public PK_fem_ {
  public :

    PK_conn_int_(dim_type nc, short_type k)
      : PK_fem_(nc, k) {
      scalar_type alpha = 0.1;

      std::fill(dof_types_.begin(), dof_types_.end(), lagrange_dof(nc));
      if (alpha != scalar_type(0)) {
        base_node G =
          gmm::mean_value(cv_node.points().begin(), cv_node.points().end());
        for (size_type i=0; i < cv_node.nb_points(); ++i)
          cv_node.points()[i] = (1-alpha)*cv_node.points()[i] + alpha*G;
        for (size_type d = 0; d < nc; ++d) {
          base_poly S(1,2);
          S[0] = -alpha * G[d] / (1-alpha);
          S[1] = 1. / (1-alpha);
          for (size_type j=0; j < nb_base(0); ++j)
            base_[j] = bgeot::poly_substitute_var(base_[j],S,d);
        }
      }
    }
  };

  static pfem conn_int_PK_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0,
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    GMM_ASSERT1(n > 0 && n < 100 && k >= 0 && k <= 150
                && double(n) == params[0].num()
                && double(k) == params[1].num(), "Bad parameters");
    pfem p = std::make_shared<PK_conn_int_>(dim_type(n), short_type(k));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*    Interior PK                                                       */
  /* ******************************************************************** */
  /* Interior PK element for arbitrary shape                              */
  /* (equivalent to DISCONTINUOUS_PK(n,k,0.1) for simplicies              */

  
  struct PK_int_ : public PK_fem_ {
  public :

    PK_int_(dim_type nc, short_type k, bgeot::pconvex_ref cvr_)
      : PK_fem_(nc, k) {
      cvr = cvr_;
      
      scalar_type alpha = 0.1;
      std::fill(dof_types_.begin(),
                  dof_types_.end(), lagrange_nonconforming_dof(nc));

      if (alpha != scalar_type(0)) {
        base_node G =
          gmm::mean_value(cv_node.points().begin(), cv_node.points().end());
        for (size_type i=0; i < cv_node.nb_points(); ++i)
          cv_node.points()[i] = (1-alpha)*cv_node.points()[i] + alpha*G;
        for (size_type d = 0; d < nc; ++d) {
          base_poly S(1,2);
          S[0] = -alpha * G[d] / (1-alpha);
          S[1] = 1. / (1-alpha);
          for (size_type j=0; j < nb_base(0); ++j)
            base_[j] = bgeot::poly_substitute_var(base_[j],S,d);
        }
      }
    }
  };

  static pfem simplex_IPK_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0,
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    GMM_ASSERT1(n > 0 && n < 100 && k >= 0 && k <= 150
                && double(n) == params[0].num()
                && double(k) == params[1].num(), "Bad parameters");
    pfem p = std::make_shared<PK_int_>(dim_type(n), short_type(k),
                                     bgeot::simplex_of_reference(dim_type(n)));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  static pfem parallelepiped_IPK_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0,
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    GMM_ASSERT1(n > 0 && n < 100 && k >= 0 && k <= 150
                && double(n) == params[0].num()
                && double(k) == params[1].num(), "Bad parameters");
    pfem p = std::make_shared<PK_int_>(dim_type(n), short_type(k),
                              bgeot::parallelepiped_of_reference(dim_type(n)));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  static pfem prism_IPK_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0,
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    GMM_ASSERT1(n > 0 && n < 100 && k >= 0 && k <= 150
                && double(n) == params[0].num()
                && double(k) == params[1].num(), "Bad parameters");
    pfem p = std::make_shared<PK_int_>(dim_type(n), short_type(k),
                                       bgeot::prism_of_reference(dim_type(n)));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }


  /* ******************************************************************** */
  /*    PK element with a bubble base function                            */
  /* ******************************************************************** */

  struct PK_with_cubic_bubble_ : public PK_fem_ {
    PK_with_cubic_bubble_(dim_type nc, short_type k);
  };

  PK_with_cubic_bubble_::PK_with_cubic_bubble_(dim_type nc, short_type k)
    : PK_fem_(nc, k) {
    unfreeze_cvs_node();
    is_lag = false; es_degree = short_type(nc+1);
    base_node pt(nc);
    size_type j;
    PK_fem_ P1(nc, 1);

    pt.fill(1./(nc+1)); /* barycenter of the convex */

    add_node(bubble1_dof(nc), pt);
    base_.resize(nb_dof(0));

    j = nb_dof(0) - 1;
    base_[j] = base_poly(nc, 0);
    base_[j].one();
    for (size_type i = 0; i < P1.nb_dof(0); i++) base_[j] *= P1.base()[i];
    // cout << "bubble = " << base_[j] << endl;
  }

  static pfem PK_with_cubic_bubble(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0,
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    GMM_ASSERT1(k < n+1, "dimensions mismatch");
    GMM_ASSERT1(n > 0 && n < 100 && k >= 0 && k <= 150 &&
                double(n) == params[0].num() && double(k) == params[1].num(),
                "Bad parameters");
    pfem p = std::make_shared<PK_with_cubic_bubble_>(dim_type(n),
                                                     short_type(k));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*    classical fem                                                     */
  /* ******************************************************************** */

  static pfem classical_fem_(bgeot::pgeometric_trans pgt,
                             short_type k, bool complete=false,
                             bool discont=false, scalar_type alpha=0) {
    THREAD_SAFE_STATIC bgeot::pgeometric_trans pgt_last = nullptr;
    THREAD_SAFE_STATIC short_type k_last = short_type(-1);
    THREAD_SAFE_STATIC pfem fem_last = nullptr;
    THREAD_SAFE_STATIC char complete_last = 0;
    THREAD_SAFE_STATIC char discont_last = 0;
    THREAD_SAFE_STATIC scalar_type alpha_last = 0;

    bool found = false;
    if (pgt_last == pgt && k_last == k && complete_last == complete &&
        discont_last == discont && alpha_last == alpha)
      return fem_last;

    dim_type n = pgt->structure()->dim();
    dim_type nbp = dim_type(pgt->basic_structure()->nb_points());
    std::stringstream name;
    std::string suffix(discont ? "_DISCONTINUOUS" : "");
    GMM_ASSERT2(discont || alpha == scalar_type(0),
                "Cannot use an alpha parameter in continuous elements.");
    std::stringstream arg_;
    if (discont && alpha != scalar_type(0))
      arg_ << "," << alpha;
    std::string arg(arg_.str());

    // Identifying if it is a torus structure
    if (bgeot::is_torus_geom_trans(pgt) && n == 3) n = 2;

    /* Identifying P1-simplexes.                                          */
    if (nbp == n+1)
      if (pgt->basic_structure() == bgeot::simplex_structure(n)) {
        name << "FEM_PK" << suffix << "(" << n << "," << k << arg << ")";
        found = true;
      }

    /* Identifying Q1-parallelepiped.                                     */
    if (!found && nbp == (size_type(1) << n))
      if (pgt->basic_structure() == bgeot::parallelepiped_structure(n)) {
        if (!complete && k == 2 && (n == 2 || n == 3) &&
            pgt->structure() == bgeot::Q2_incomplete_structure(n)) {
          GMM_ASSERT2(alpha == scalar_type(0),
                      "Cannot use an alpha parameter in incomplete Q2"
                      " elements.");
          name << "FEM_Q2_INCOMPLETE" << suffix << "(" << n << ")";
        } else
          name << "FEM_QK" << suffix << "(" << n << "," << k << arg << ")";
        found = true;
      }

    /* Identifying Q1-prisms.                                             */
    if (!found && nbp == 2 * n)
      if (pgt->basic_structure() == bgeot::prism_P1_structure(n)) {
        if (!complete && k == 2 && n == 3 &&
            pgt->structure() == bgeot::prism_incomplete_P2_structure()) {
          GMM_ASSERT2(alpha == scalar_type(0),
                      "Cannot use an alpha parameter in incomplete prism"
                      " elements.");
          name << "FEM_PRISM_INCOMPLETE_P2" << suffix;
        } else
          name << "FEM_PRISM_PK" << suffix << "(" << n << "," << k << arg << ")";
        found = true;
      }

    /* Identifying pyramids.                                              */
    if (!found && nbp == 5)
      if (pgt->basic_structure() == bgeot::pyramid_QK_structure(1)) {
        if (!complete && k == 2 &&
            pgt->structure() == bgeot::pyramid_Q2_incomplete_structure()) {
          GMM_ASSERT2(alpha == scalar_type(0),
                      "Cannot use an alpha parameter in incomplete pyramid"
                      " elements.");
          name << "FEM_PYRAMID_Q2_INCOMPLETE" << suffix;
        } else
          name << "FEM_PYRAMID_QK" << suffix << "(" << k << arg << ")";
        found = true;;
      }

    // To be completed

    GMM_ASSERT1(found, "This element is not taken into account. Contact us");
    fem_last = fem_descriptor(name.str());
    pgt_last = pgt;
    k_last = k;
    complete_last = complete;
    discont_last = discont;
    alpha_last = alpha;
    return fem_last;
  }

  pfem classical_fem(bgeot::pgeometric_trans pgt, short_type k,
                     bool complete) {
    return classical_fem_(pgt, k, complete);
  }

  pfem classical_discontinuous_fem(bgeot::pgeometric_trans pgt, short_type k,
                                   scalar_type alpha, bool complete) {
    return classical_fem_(pgt, k, complete, true, alpha);
  }

  /* ******************************************************************** */
  /*    Naming system                                                     */
  /* ******************************************************************** */

  pfem structured_composite_fem_method(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies);
  pfem PK_composite_hierarch_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies);
  pfem PK_composite_full_hierarch_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies);
  pfem HCT_triangle_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies);
  pfem reduced_HCT_triangle_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies);
  pfem reduced_quadc1p3_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies);
  pfem quadc1p3_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies);
  pfem P1bubbletriangle_fem(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies);
  pfem hho_method(fem_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies);
  
  struct fem_naming_system : public dal::naming_system<virtual_fem> {
    fem_naming_system() : dal::naming_system<virtual_fem>("FEM") {
      add_suffix("HERMITE", Hermite_fem);
      add_suffix("ARGYRIS", triangle_Argyris_fem);
      add_suffix("MORLEY", triangle_Morley_fem);
      add_suffix("PK", PK_fem);
      add_suffix("QK", QK_fem);
      add_suffix("QK_DISCONTINUOUS", QK_discontinuous_fem);
      add_suffix("PRISM_PK", prism_PK_fem);
      add_suffix("PK_PRISM", prism_PK_fem); // for backwards compatibility
      add_suffix("PK_DISCONTINUOUS", PK_discontinuous_fem);
      add_suffix("PRISM_PK_DISCONTINUOUS", prism_PK_discontinuous_fem);
      add_suffix("PK_PRISM_DISCONTINUOUS", prism_PK_discontinuous_fem); // for backwards compatibility
      add_suffix("SIMPLEX_IPK", simplex_IPK_fem);
      add_suffix("PRISM_IPK", prism_IPK_fem);
      add_suffix("QUAD_IPK", parallelepiped_IPK_fem);
      add_suffix("SIMPLEX_CIPK", conn_int_PK_fem);
      add_suffix("PK_WITH_CUBIC_BUBBLE", PK_with_cubic_bubble);
      add_suffix("PRODUCT", product_fem);
      add_suffix("P1_NONCONFORMING", P1_nonconforming_fem);
      add_suffix("P1_BUBBLE_FACE", P1_with_bubble_on_a_face);
      add_suffix("P1_BUBBLE_FACE_LAG", P1_with_bubble_on_a_face_lagrange);
      add_suffix("P1_PIECEWISE_LINEAR_BUBBLE", P1bubbletriangle_fem);
      add_suffix("GEN_HIERARCHICAL", gen_hierarchical_fem);
      add_suffix("PK_HIERARCHICAL", PK_hierarch_fem);
      add_suffix("QK_HIERARCHICAL", QK_hierarch_fem);
      add_suffix("PRISM_PK_HIERARCHICAL", prism_PK_hierarch_fem);
      add_suffix("PK_PRISM_HIERARCHICAL", prism_PK_hierarch_fem); // for backwards compatibility
      add_suffix("STRUCTURED_COMPOSITE", structured_composite_fem_method);
      add_suffix("PK_HIERARCHICAL_COMPOSITE", PK_composite_hierarch_fem);
      add_suffix("PK_FULL_HIERARCHICAL_COMPOSITE",
                 PK_composite_full_hierarch_fem);
      add_suffix("PK_GAUSSLOBATTO1D", PK_GL_fem);
      add_suffix("QUAD_CIPK", CIPK_SQUARE);
      add_suffix("Q2_INCOMPLETE", Q2_incomplete_fem);
      add_suffix("Q2_INCOMPLETE_DISCONTINUOUS", Q2_incomplete_discontinuous_fem);
      add_suffix("HCT_TRIANGLE", HCT_triangle_fem);
      add_suffix("REDUCED_HCT_TRIANGLE", reduced_HCT_triangle_fem);
      add_suffix("QUADC1_COMPOSITE", quadc1p3_fem);
      add_suffix("REDUCED_QUADC1_COMPOSITE", reduced_quadc1p3_fem);
      add_suffix("HHO", hho_method);
      add_suffix("RT0", P1_RT0);
      add_suffix("RT0Q", P1_RT0Q);
      add_suffix("NEDELEC", P1_nedelec);
      add_suffix("PYRAMID_QK", pyramid_QK_fem);
      add_suffix("PYRAMID_QK_DISCONTINUOUS", pyramid_QK_disc_fem);
      add_suffix("PYRAMID_LAGRANGE", pyramid_QK_fem); // for backwards compatibility
      add_suffix("PYRAMID_DISCONTINUOUS_LAGRANGE", pyramid_QK_disc_fem); // for backwards compatibility
      add_suffix("PYRAMID_Q2_INCOMPLETE", pyramid_Q2_incomplete_fem);
      add_suffix("PYRAMID_Q2_INCOMPLETE_DISCONTINUOUS",
                 pyramid_Q2_incomplete_disc_fem);
      add_suffix("PRISM_INCOMPLETE_P2", prism_incomplete_P2_fem);
      add_suffix("PRISM_INCOMPLETE_P2_DISCONTINUOUS",
                 prism_incomplete_P2_disc_fem);
    }
  };

  // get a fem descriptor from a string name of a fem.
  pfem fem_descriptor(const std::string &name) {
    size_type i = 0;
    pfem  pf = dal::singleton<fem_naming_system>::instance().method(name, i);
    const_cast<virtual_fem &>(*pf).debug_name()
      = dal::singleton<fem_naming_system>::instance().shorter_name_of_method(pf);
    return pf;
  }

  // get the string name of a fem descriptor.
  std::string name_of_fem(pfem p) {

    auto &instance = dal::singleton<fem_naming_system>::instance();
    auto *p_torus = dynamic_cast<const torus_fem *>(p.get());
    if (p_torus) return instance.shorter_name_of_method(p_torus->get_original_pfem());
    return instance.shorter_name_of_method(p);

    return dal::singleton<fem_naming_system>::instance().
      shorter_name_of_method(p);
  }

  // allows the add of a fem.
  void add_fem_name(std::string name,
                    dal::naming_system<virtual_fem>::pfunction f) {
    dal::singleton<fem_naming_system>::instance().add_suffix(name, f);
  }

  /* ******************************************************************** */
  /*    Aliases functions                                                 */
  /* ******************************************************************** */

  pfem PK_fem(size_type n, short_type k) {
    THREAD_SAFE_STATIC pfem pf = nullptr;
    THREAD_SAFE_STATIC size_type d = size_type(-2);
    THREAD_SAFE_STATIC short_type r = short_type(-2);
    if (d != n || r != k) {
      std::stringstream name;
      name << "FEM_PK(" << n << "," << k << ")";
      pf = fem_descriptor(name.str());
      d = n; r = k;
    }
    return pf;
  }

  pfem QK_fem(size_type n, short_type k) {
    THREAD_SAFE_STATIC pfem pf = nullptr;
    THREAD_SAFE_STATIC size_type d = size_type(-2);
    THREAD_SAFE_STATIC short_type r = short_type(-2);
    if (d != n || r != k) {
      std::stringstream name;
      name << "FEM_QK(" << n << "," << k << ")";
      pf = fem_descriptor(name.str());
      d = n; r = k;
    }
    return pf;
  }

  pfem prism_PK_fem(size_type n, short_type k) {
    THREAD_SAFE_STATIC pfem pf = nullptr;
    THREAD_SAFE_STATIC size_type d = size_type(-2);
    THREAD_SAFE_STATIC short_type r = short_type(-2);
    if (d != n || r != k) {
      std::stringstream name;
      name << "FEM_PRISM_PK(" << n << "," << k << ")";
      pf = fem_descriptor(name.str());
      d = n; r = k;
    }
    return pf;
  }


  /* ********************************************************************* */
  /*       Precomputation on fem.                                          */
  /* ********************************************************************* */

  DAL_DOUBLE_KEY(pre_fem_key_, pfem, bgeot::pstored_point_tab);

  fem_precomp_::fem_precomp_(const pfem pff, const bgeot::pstored_point_tab ps) :
    pf(pff), pspt(ps) {
    DAL_STORED_OBJECT_DEBUG_CREATED(this, "Fem_precomp");
    for (const auto &p : *pspt)
      GMM_ASSERT1(p.size() == pf->dim(), "dimensions mismatch");
  }

  void fem_precomp_::init_val() const {
    c.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i)
      pf->base_value((*pspt)[i], c[i]);
  }

  void fem_precomp_::init_grad() const {
    pc.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i)
      pf->grad_base_value((*pspt)[i], pc[i]);
  }

  void fem_precomp_::init_hess() const {
    hpc.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i)
      pf->hess_base_value((*pspt)[i], hpc[i]);
  }

  pfem_precomp fem_precomp(pfem pf, bgeot::pstored_point_tab pspt,
                           dal::pstatic_stored_object dep) {
    dal::pstatic_stored_object_key pk = std::make_shared<pre_fem_key_>(pf,pspt);
    dal::pstatic_stored_object o = dal::search_stored_object(pk);
    if (o) return std::dynamic_pointer_cast<const fem_precomp_>(o);
    pfem_precomp p = std::make_shared<fem_precomp_>(pf, pspt);
    dal::add_stored_object(pk, p, pspt, dal::AUTODELETE_STATIC_OBJECT);
    if (dal::exists_stored_object(pf)) dal::add_dependency(p, pf);
    if (dep) dal::add_dependency(p, dep);
    return p;
  }

  void fem_precomp_pool::clear() {
    for (const pfem_precomp &p : precomps)
      dal::del_stored_object(p, true);
    precomps.clear();
  }


}  /* end of namespace getfem.                                            */
