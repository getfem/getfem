/*===========================================================================

 Copyright (C) 2013-2020 Yves Renard

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

#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_generic_assembly_tree.h"
#include "getfem/getfem_generic_assembly_semantic.h"
#include "getfem/getfem_generic_assembly_compile_and_exec.h"
#include "getfem/getfem_generic_assembly_functions_and_operators.h"

// #define GA_USES_BLAS // not so interesting, at least for debian blas

// #define GA_DEBUG_INFO(a) { cout << a << endl; }
#define GA_DEBUG_INFO(a)



namespace getfem {


  template <class VEC1, class VEC2>
  inline void copy_scaled_4(const VEC1 &v1, const scalar_type a, VEC2 &v2) {
    auto it1 = v1.begin();
    auto it2 = v2.begin(), it2e = v2.end();
    size_type nd = v1.size() >> 2;
    for (size_type i = 0; i < nd; ++i) {
      *it2++ = (*it1++) * a;
      *it2++ = (*it1++) * a;
      *it2++ = (*it1++) * a;
      *it2++ = (*it1++) * a;
    }
    for (; it2 != it2e;)
      *it2++ = (*it1++) * a;
  }

  template <class VEC1, class VEC2>
  inline void add_scaled_4(const VEC1 &v1, const scalar_type a, VEC2 &v2) {
    auto it1 = v1.begin();
    auto it2 = v2.begin(), it2e = v2.end();
    size_type nd = v1.size() >> 2;
    for (size_type i = 0; i < nd; ++i) {
      *it2++ += (*it1++) * a;
      *it2++ += (*it1++) * a;
      *it2++ += (*it1++) * a;
      *it2++ += (*it1++) * a;
    }
    for (; it2 != it2e;)
      *it2++ += (*it1++) * a;
  }

  template <class VEC1, class VEC2>
  inline void copy_scaled_8(const VEC1 &v1, const scalar_type a, VEC2 &v2) {
    auto it1 = v1.begin();
    auto it2 = v2.begin(), it2e = v2.end();
    size_type nd = v1.size() >> 3;
    for (size_type i = 0; i < nd; ++i) {
      *it2++ = (*it1++) * a;
      *it2++ = (*it1++) * a;
      *it2++ = (*it1++) * a;
      *it2++ = (*it1++) * a;
      *it2++ = (*it1++) * a;
      *it2++ = (*it1++) * a;
      *it2++ = (*it1++) * a;
      *it2++ = (*it1++) * a;
    }
    for (; it2 != it2e;)
      *it2++ = (*it1++) * a;
  }

  template <class VEC1, class VEC2>
  inline void add_scaled_8(const VEC1 &v1, const scalar_type a, VEC2 &v2) {
    auto it1 = v1.begin();
    auto it2 = v2.begin(), it2e = v2.end();
    size_type nd = v1.size() >> 3;
    for (size_type i = 0; i < nd; ++i) {
      *it2++ += (*it1++) * a;
      *it2++ += (*it1++) * a;
      *it2++ += (*it1++) * a;
      *it2++ += (*it1++) * a;
      *it2++ += (*it1++) * a;
      *it2++ += (*it1++) * a;
      *it2++ += (*it1++) * a;
      *it2++ += (*it1++) * a;
    }
    for (; it2 != it2e;)
      *it2++ += (*it1++) * a;
  }

  bool operator <(const gauss_pt_corresp &gpc1,
                  const gauss_pt_corresp &gpc2) {
    if (gpc1.pai != gpc2.pai)
      return (gpc1.pai  <  gpc2.pai );
    if (gpc1.nodes.size() !=  gpc2.nodes.size())
      return (gpc1.nodes.size() < gpc2.nodes.size());
    for (size_type i = 0; i < gpc1.nodes.size(); ++i)
      if (gpc1.nodes[i] != gpc2.nodes[i])
        return (gpc1.nodes[i] < gpc2.nodes[i]);
    if (gpc1.pgt1 != gpc2.pgt1)
      return (gpc1.pgt1 <  gpc2.pgt1);
    if (gpc1.pgt2 !=  gpc2.pgt2)
      return (gpc1.pgt2 <  gpc2.pgt2);
    return false;
  }

  bool operator <(const ga_instruction_set::region_mim &rm1,
                  const ga_instruction_set::region_mim &rm2) {
    if (rm1.mim() != rm2.mim()) return (rm1.mim() < rm2.mim());
    if (rm1.region() != rm2.region()) return (rm1.region() < rm2.region());
    return (rm1.psd() < rm2.psd());
  }

  //=========================================================================
  // Instructions for compilation: basic optimized operations on tensors
  //=========================================================================

  struct ga_instruction_extract_local_im_data : public ga_instruction {
    base_tensor &t;
    const im_data &imd;
    papprox_integration &pai;
    const base_vector &U;
    const fem_interpolation_context &ctx;
    size_type qdim, cv_old;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: extract local im data");
      size_type cv = ctx.convex_num();
      if (cv != cv_old) {
        cv_old = cv;
        GMM_ASSERT1(imd.linked_mesh_im().int_method_of_element(cv)
                    ->approx_method() == pai, "Im data have to be used only "
                    "on their original integration method.");
      }
      size_type ipt = imd.filtered_index_of_point(cv, ctx.ii());
      GMM_ASSERT1(ipt != size_type(-1),
                  "Im data with no data on the current integration point.");
      auto it = U.begin()+ipt*qdim;
      std::copy(it, it+qdim, t.begin());
      return 0;
    }
    ga_instruction_extract_local_im_data
    (base_tensor &t_, const im_data &imd_, const base_vector &U_,
     papprox_integration &pai_, const fem_interpolation_context &ctx_,
     size_type qdim_)
      : t(t_), imd(imd_), pai(pai_), U(U_), ctx(ctx_), qdim(qdim_),
        cv_old(-1)
    {}
  };

  struct ga_instruction_slice_local_dofs : public ga_instruction {
    const mesh_fem &mf;
    const base_vector &U;
    const fem_interpolation_context &ctx;
    base_vector &coeff;
    size_type qmult1, qmult2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Slice local dofs");
      GMM_ASSERT1(qmult1 != 0 && qmult2 != 0, "Internal error");
      slice_vector_on_basic_dof_of_element(mf, U, ctx.convex_num(),
                                           coeff, qmult1, qmult2);
      return 0;
    }
    ga_instruction_slice_local_dofs(const mesh_fem &mf_, const base_vector &U_,
                                    const fem_interpolation_context &ctx_,
                                    base_vector &coeff_,
                                    size_type qmult1_, size_type qmult2_)
      : mf(mf_), U(U_), ctx(ctx_), coeff(coeff_),
        qmult1(qmult1_), qmult2(qmult2_) {}
  };

  struct ga_instruction_update_pfp : public ga_instruction {
    const mesh_fem &mf;
    const fem_interpolation_context &ctx;
    fem_precomp_pool &fp_pool;
    pfem_precomp &pfp;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Pfp update");
      if (ctx.have_pgp()) {
        size_type cv = ctx.is_convex_num_valid()
                     ? ctx.convex_num() : mf.convex_index().first_true();
        pfem pf = mf.fem_of_element(cv);
        if (!pfp || pf != pfp->get_pfem() ||
            ctx.pgp()->get_ppoint_tab() != pfp->get_ppoint_tab()) {
          pfp = fp_pool(pf, ctx.pgp()->get_ppoint_tab());
        }
      } else {
        pfp = 0;
      }
      return 0;
    }

    ga_instruction_update_pfp(const mesh_fem &mf_, pfem_precomp &pfp_,
                              const fem_interpolation_context &ctx_,
                              fem_precomp_pool &fp_pool_)
      : mf(mf_), ctx(ctx_), fp_pool(fp_pool_), pfp(pfp_) {}
  };

  struct ga_instruction_first_ind_tensor : public ga_instruction {
    base_tensor &t;
    const fem_interpolation_context &ctx;
    size_type qdim;
    const mesh_fem *mfn, **mfg;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: adapt first index of tensor");
      const mesh_fem &mf = *(mfg ? *mfg : mfn);
      GA_DEBUG_ASSERT(mfg ? *mfg : mfn, "Internal error");
      size_type cv_1 = ctx.is_convex_num_valid()
                     ? ctx.convex_num() : mf.convex_index().first_true();
      pfem pf = mf.fem_of_element(cv_1);
      GMM_ASSERT1(pf, "An element without finite element method defined");
      size_type Qmult = qdim / pf->target_dim();
      size_type s = pf->nb_dof(cv_1) * Qmult;
      if (t.sizes()[0] != s)
        { bgeot::multi_index mi = t.sizes(); mi[0] = s; t.adjust_sizes(mi); }
      return 0;
    }

    ga_instruction_first_ind_tensor(base_tensor &t_,
                                    const fem_interpolation_context &ctx_,
                                    size_type qdim_, const mesh_fem *mfn_,
                                    const mesh_fem **mfg_)
      : t(t_),  ctx(ctx_), qdim(qdim_), mfn(mfn_), mfg(mfg_) {}
  };

  struct ga_instruction_second_ind_tensor
    : public ga_instruction_first_ind_tensor {

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: adapt second index of tensor");
      const mesh_fem &mf = *(mfg ? *mfg : mfn);
      size_type cv_1 = ctx.is_convex_num_valid()
                     ? ctx.convex_num() : mf.convex_index().first_true();
      pfem pf = mf.fem_of_element(cv_1);
      GMM_ASSERT1(pf, "An element without finite element methode defined");
      size_type Qmult = qdim / pf->target_dim();
      size_type s = pf->nb_dof(cv_1) * Qmult;
      if (t.sizes()[1] != s)
        { bgeot::multi_index mi = t.sizes(); mi[1] = s; t.adjust_sizes(mi); }
      return 0;
    }

    ga_instruction_second_ind_tensor(base_tensor &t_,
                                     fem_interpolation_context &ctx_,
                                     size_type qdim_, const mesh_fem *mfn_,
                                     const mesh_fem **mfg_)
      : ga_instruction_first_ind_tensor(t_, ctx_, qdim_, mfn_, mfg_) {}

  };

  struct ga_instruction_two_first_ind_tensor : public ga_instruction {
    base_tensor &t;
    const fem_interpolation_context &ctx1, &ctx2;
    size_type qdim1;
    const mesh_fem *mfn1, **mfg1;
    size_type qdim2;
    const mesh_fem *mfn2, **mfg2;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: adapt two first indices of tensor");
      const mesh_fem &mf1 = *(mfg1 ? *mfg1 : mfn1);
      const mesh_fem &mf2 = *(mfg2 ? *mfg2 : mfn2);
      size_type cv_1 = ctx1.is_convex_num_valid()
                     ? ctx1.convex_num() : mf1.convex_index().first_true();
      size_type cv_2 = ctx2.is_convex_num_valid()
                     ? ctx2.convex_num() : mf2.convex_index().first_true();
      pfem pf1 = mf1.fem_of_element(cv_1);
      GMM_ASSERT1(pf1, "An element without finite element method defined");
      pfem pf2 = mf2.fem_of_element(cv_2);
      GMM_ASSERT1(pf2, "An element without finite element method defined");
      size_type Qmult1 = qdim1 / pf1->target_dim();
      size_type s1 = pf1->nb_dof(cv_1) * Qmult1;
      size_type Qmult2 = qdim2 / pf2->target_dim();
      size_type s2 = pf2->nb_dof(cv_2) * Qmult2;
      GMM_ASSERT1(s1 > 0 && s2 >0, "Element without degrees of freedom");
      if (t.sizes()[0] != s1 || t.sizes()[1] != s2) {
        bgeot::multi_index mi = t.sizes();
        mi[0] = s1; mi[1] = s2;
        t.adjust_sizes(mi);
      }
      return 0;
    }

    ga_instruction_two_first_ind_tensor
    (base_tensor &t_, const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     size_type qdim1_, const mesh_fem *mfn1_, const mesh_fem **mfg1_,
     size_type qdim2_, const mesh_fem *mfn2_, const mesh_fem **mfg2_)
      : t(t_),  ctx1(ctx1_), ctx2(ctx2_), qdim1(qdim1_), mfn1(mfn1_),
        mfg1(mfg1_), qdim2(qdim2_), mfn2(mfn2_), mfg2(mfg2_) {}
  };


  struct ga_instruction_X_component : public ga_instruction {
    scalar_type &t;
    const fem_interpolation_context &ctx;
    size_type n;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: X component");
      t = ctx.xreal()[n];
      return 0;
    }

    ga_instruction_X_component
    (scalar_type &t_, const fem_interpolation_context &ctx_, size_type n_)
      : t(t_),  ctx(ctx_), n(n_) {}
  };

  struct ga_instruction_X : public ga_instruction {
    base_tensor &t;
    const fem_interpolation_context &ctx;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: X");
      GA_DEBUG_ASSERT(t.size() == ctx.xreal().size(), "dimensions mismatch");
      gmm::copy(ctx.xreal(), t.as_vector());
      return 0;
    }

    ga_instruction_X(base_tensor &t_, const fem_interpolation_context &ctx_)
      : t(t_),  ctx(ctx_) {}
  };

  struct ga_instruction_copy_small_vect : public ga_instruction {
    base_tensor &t;
    const base_small_vector &vec;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: copy small vector");
      GMM_ASSERT1(t.size() == vec.size(), "Invalid vector size.");
      gmm::copy(vec, t.as_vector());
      return 0;
    }
    ga_instruction_copy_small_vect(base_tensor &t_,
                                   const base_small_vector &vec_)
      : t(t_), vec(vec_)  {}
  };

  struct ga_instruction_copy_Normal : public ga_instruction_copy_small_vect {

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: unit normal vector");
      GMM_ASSERT1(t.size() == vec.size(), "Invalid outward unit normal "
                  "vector. Possible reasons: not on boundary or "
                  "transformation failed.");
      gmm::copy(vec, t.as_vector());
      return 0;
    }
    ga_instruction_copy_Normal(base_tensor &t_,
                               const base_small_vector &Normal_)
      : ga_instruction_copy_small_vect(t_, Normal_)  {}
  };

  struct ga_instruction_level_set_normal_vector : public ga_instruction {
    base_tensor &t;
    const mesh_im_level_set *mimls;
    const fem_interpolation_context &ctx;
    base_small_vector vec;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: unit normal vector to a level-set");
      mimls->compute_normal_vector(ctx, vec);
      GMM_ASSERT1(t.size() == vec.size(), "Invalid outward unit normal "
                  "vector. Possible reasons: not on boundary or "
                  "transformation failed.");
      gmm::copy(vec, t.as_vector());
      return 0;
    }
    ga_instruction_level_set_normal_vector
    (base_tensor &t_, const mesh_im_level_set *mimls_,
     const fem_interpolation_context &ctx_)
      : t(t_), mimls(mimls_), ctx(ctx_), vec(t.size())  {}
  };

  struct ga_instruction_element_size : public ga_instruction {
    base_tensor &t;
    scalar_type &es;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: element_size");
      GMM_ASSERT1(t.size() == 1, "Invalid element size.");
      t[0] = es;
      return 0;
    }
    ga_instruction_element_size(base_tensor &t_, scalar_type &es_)
      : t(t_), es(es_)  {}
  };

  struct ga_instruction_element_K : public ga_instruction {
    base_tensor &t;
    const fem_interpolation_context &ctx;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: element_K");
      GMM_ASSERT1(t.size() == (ctx.K()).size(), "Invalid tensor size.");
      gmm::copy(ctx.K().as_vector(), t.as_vector());
      return 0;
    }
    ga_instruction_element_K(base_tensor &t_,
                             const fem_interpolation_context &ct)
      : t(t_), ctx(ct)  {}
  };

  struct ga_instruction_element_B : public ga_instruction {
    base_tensor &t;
    const fem_interpolation_context &ctx;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: element_B");
      GMM_ASSERT1(t.size() == (ctx.B()).size(), "Invalid tensor size.");
      gmm::copy(ctx.B().as_vector(), t.as_vector());
      return 0;
    }
    ga_instruction_element_B(base_tensor &t_,
                             const fem_interpolation_context &ct)
      : t(t_), ctx(ct)  {}
  };

  struct ga_instruction_val_base : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;
    const mesh_fem &mf;
    const pfem_precomp &pfp;

    virtual int exec() { // --> t(ndof,target_dim)
      GA_DEBUG_INFO("Instruction: compute value of base functions");
      // if (ctx.have_pgp()) ctx.set_pfp(pfp);
      // else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      // GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      // ctx.base_value(t);
      if (ctx.have_pgp()) ctx.pfp_base_value(t, pfp);
      else {
        ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
        GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
        ctx.base_value(t);
      }
      return 0;
    }

    ga_instruction_val_base(base_tensor &tt, fem_interpolation_context &ct,
                            const mesh_fem &mf_, const pfem_precomp &pfp_)
      : t(tt), ctx(ct), mf(mf_), pfp(pfp_) {}
  };

  struct ga_instruction_xfem_plus_val_base : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;
    const mesh_fem &mf;
    pfem_precomp &pfp;

    virtual int exec() { // --> t(ndof,target_dim)
      GA_DEBUG_INFO("Instruction: compute value of base functions");
      if (ctx.have_pgp()) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      int old_xfem_side = ctx.xfem_side();
      ctx.set_xfem_side(1);
      ctx.base_value(t);
      ctx.set_xfem_side(old_xfem_side);
      return 0;
    }

    ga_instruction_xfem_plus_val_base(base_tensor &tt,
                                      fem_interpolation_context &ct,
                                      const mesh_fem &mf_, pfem_precomp &pfp_)
      : t(tt), ctx(ct), mf(mf_), pfp(pfp_) {}
  };

  struct ga_instruction_xfem_minus_val_base : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;
    const mesh_fem &mf;
    pfem_precomp &pfp;

    virtual int exec() { // --> t(ndof,target_dim)
      GA_DEBUG_INFO("Instruction: compute value of base functions");
      if (ctx.have_pgp()) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      int old_xfem_side = ctx.xfem_side();
      ctx.set_xfem_side(-1);
      ctx.base_value(t);
      ctx.set_xfem_side(old_xfem_side);
      return 0;
    }

    ga_instruction_xfem_minus_val_base
    (base_tensor &tt, fem_interpolation_context &ct,
     const mesh_fem &mf_, pfem_precomp &pfp_)
      : t(tt), ctx(ct), mf(mf_), pfp(pfp_) {}
  };

  struct ga_instruction_grad_base : public ga_instruction_val_base {

    virtual int exec() { // --> t(ndof,target_dim,N)
      GA_DEBUG_INFO("Instruction: compute gradient of base functions");
      // if (ctx.have_pgp()) ctx.set_pfp(pfp);
      // else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      // GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      // ctx.grad_base_value(t);
      if (ctx.have_pgp()) ctx.pfp_grad_base_value(t, pfp);
      else {
        ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
        GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
        ctx.grad_base_value(t);
      }
      return 0;
    }

    ga_instruction_grad_base(base_tensor &tt, fem_interpolation_context &ct,
                             const mesh_fem &mf_, pfem_precomp &pfp_)
    : ga_instruction_val_base(tt, ct, mf_, pfp_)
    {}
  };

  struct ga_instruction_xfem_plus_grad_base : public ga_instruction_val_base {

    virtual int exec() { // --> t(ndof,target_dim,N)
      GA_DEBUG_INFO("Instruction: compute gradient of base functions");
      if (ctx.have_pgp()) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      int old_xfem_side = ctx.xfem_side();
      ctx.set_xfem_side(1);
      ctx.grad_base_value(t);
      ctx.set_xfem_side(old_xfem_side);
      return 0;
    }

    ga_instruction_xfem_plus_grad_base
    (base_tensor &tt, fem_interpolation_context &ct,
     const mesh_fem &mf_, pfem_precomp &pfp_)
    : ga_instruction_val_base(tt, ct, mf_, pfp_)
    {}
  };

  struct ga_instruction_xfem_minus_grad_base : public ga_instruction_val_base {

    virtual int exec() { // --> t(ndof,target_dim,N)
      GA_DEBUG_INFO("Instruction: compute gradient of base functions");
      if (ctx.have_pgp()) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      int old_xfem_side = ctx.xfem_side();
      ctx.set_xfem_side(-1);
      ctx.grad_base_value(t);
      ctx.set_xfem_side(old_xfem_side);
      return 0;
    }

    ga_instruction_xfem_minus_grad_base
    (base_tensor &tt, fem_interpolation_context &ct,
     const mesh_fem &mf_, pfem_precomp &pfp_)
    : ga_instruction_val_base(tt, ct, mf_, pfp_)
    {}
  };


  struct ga_instruction_hess_base : public ga_instruction_val_base {

    virtual int exec() { // --> t(ndof,target_dim,N*N)
      GA_DEBUG_INFO("Instruction: compute Hessian of base functions");
      if (ctx.have_pgp()) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      ctx.hess_base_value(t);
      return 0;
    }

    ga_instruction_hess_base(base_tensor &tt, fem_interpolation_context &ct,
                             const mesh_fem &mf_, pfem_precomp &pfp_)
    : ga_instruction_val_base(tt, ct, mf_, pfp_)
    {}
  };

  struct ga_instruction_xfem_plus_hess_base : public ga_instruction_val_base {

    virtual int exec() { // --> t(ndof,target_dim,N*N)
      GA_DEBUG_INFO("Instruction: compute Hessian of base functions");
      if (ctx.have_pgp()) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      int old_xfem_side = ctx.xfem_side();
      ctx.set_xfem_side(1);
      ctx.hess_base_value(t);
      ctx.set_xfem_side(old_xfem_side);
      return 0;
    }

    ga_instruction_xfem_plus_hess_base
    (base_tensor &tt, fem_interpolation_context &ct,
     const mesh_fem &mf_, pfem_precomp &pfp_)
    : ga_instruction_val_base(tt, ct, mf_, pfp_)
    {}
  };

  struct ga_instruction_xfem_minus_hess_base : public ga_instruction_val_base {

    virtual int exec() { // --> t(ndof,target_dim,N*N)
      GA_DEBUG_INFO("Instruction: compute Hessian of base functions");
      if (ctx.have_pgp()) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      int old_xfem_side = ctx.xfem_side();
      ctx.set_xfem_side(-1);
      ctx.hess_base_value(t);
      ctx.set_xfem_side(old_xfem_side);
      return 0;
    }

    ga_instruction_xfem_minus_hess_base
    (base_tensor &tt, fem_interpolation_context &ct,
     const mesh_fem &mf_, pfem_precomp &pfp_)
    : ga_instruction_val_base(tt, ct, mf_, pfp_)
    {}
  };

  struct ga_instruction_val : public ga_instruction {
    scalar_type &a;
    base_tensor &t;
    const base_tensor &Z;
    const base_vector &coeff;
    size_type qdim;
    // Z(ndof,target_dim), coeff(Qmult,ndof) --> t(target_dim*Qmult)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: variable value");
      size_type ndof = Z.sizes()[0];
      if (!ndof) { gmm::clear(t.as_vector()); return 0; }
      GA_DEBUG_ASSERT(t.size() == qdim, "dimensions mismatch");

      if (qdim == 1) {
        GA_DEBUG_ASSERT(gmm::vect_size(coeff) == ndof,
                        "Wrong size for coeff vector");
        auto itc = coeff.begin(); auto itZ = Z.begin();
        a = (*itc++) * (*itZ++);
        while (itc != coeff.end()) a += (*itc++) * (*itZ++);
      } else {
        size_type target_dim = Z.sizes()[1];
        if (target_dim == 1) {
          GA_DEBUG_ASSERT(gmm::vect_size(coeff) == ndof*qdim,
                          "Wrong size for coeff vector");
          auto itc = coeff.begin(); auto itZ = Z.begin();
          for (auto it = t.begin(); it != t.end(); ++it)
            *it = (*itc++) * (*itZ);
          ++itZ;
          for (size_type j = 1; j < ndof; ++j, ++itZ) {
            for (auto it = t.begin(); it != t.end(); ++it)
              *it += (*itc++) * (*itZ);
          }
        } else {
          size_type Qmult = qdim / target_dim;
          GA_DEBUG_ASSERT(gmm::vect_size(coeff) == ndof*Qmult,
                          "Wrong size for coeff vector");

          gmm::clear(t.as_vector());
          auto itc = coeff.begin();
          for (size_type j = 0; j < ndof; ++j) {
            auto it = t.begin();
            for (size_type q = 0; q < Qmult; ++q, ++itc) {
              for (size_type r = 0; r < target_dim; ++r)
                *it++ += (*itc) * Z[j + r*ndof];
            }
          }
        }
      }
      return 0;
    }

    ga_instruction_val(base_tensor &tt, const base_tensor &Z_,
                       const base_vector &co, size_type q)
      : a(tt[0]), t(tt), Z(Z_), coeff(co), qdim(q) {}
  };

  struct ga_instruction_grad : public ga_instruction_val {
    // Z(ndof,target_dim,N), coeff(Qmult,ndof) --> t(target_dim*Qmult,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: gradient");
      size_type ndof = Z.sizes()[0];
      if (!ndof) { gmm::clear(t.as_vector()); return 0; }
      size_type N = Z.sizes()[2];
      if (qdim == 1) {
        GA_DEBUG_ASSERT(t.size() == N, "dimensions mismatch");
        GA_DEBUG_ASSERT(coeff.size() == ndof, "Wrong size for coeff vector");
        auto itZ = Z.begin();
        for (auto it = t.begin(); it != t.end(); ++it) {
          auto itc = coeff.begin();
          *it =  (*itc++) * (*itZ++);
          while (itc != coeff.end()) *it += (*itc++) * (*itZ++);
        }
      } else {
        size_type target_dim = Z.sizes()[1];
        if (target_dim == 1) {
          GA_DEBUG_ASSERT(t.size() == N*qdim, "dimensions mismatch");
          GA_DEBUG_ASSERT(coeff.size() == ndof*qdim,
                          "Wrong size for coeff vector");
          for (size_type q = 0; q < qdim; ++q) {
            auto itZ = Z.begin(); auto it = t.begin() + q;
            for (size_type k = 0; k < N; ++k) {
              if (k)  it += qdim;
              auto itc = coeff.begin() + q;
              *it = (*itc) * (*itZ++);
              for (size_type j = 1; j < ndof; ++j)
                { itc += qdim; *it += (*itc) * (*itZ++); }
            }
          }
        } else {
          size_type Qmult = qdim / target_dim;
          GA_DEBUG_ASSERT(t.size() == N*qdim, "dimensions mismatch");
          GA_DEBUG_ASSERT(coeff.size() == ndof*Qmult,
                          "Wrong size for coeff vector");
          gmm::clear(t.as_vector());
          for (size_type q = 0; q < Qmult; ++q) {
            auto itZ = Z.begin();
            for (size_type k = 0; k < N; ++k)
              for (size_type r = 0; r < target_dim; ++r)
                for (size_type j = 0; j < ndof; ++j)
                  t[r + q*target_dim + k*qdim] += coeff[j*Qmult+q] * (*itZ++);
          }
        }
      }
      return 0;
    }

    ga_instruction_grad(base_tensor &tt, const base_tensor &Z_,
                        const base_vector &co, size_type q)
    : ga_instruction_val(tt, Z_, co, q)
    {}

  };

  struct ga_instruction_hess : public ga_instruction_val {
    // Z(ndof,target_dim,N*N), coeff(Qmult,ndof) --> t(target_dim*Qmult,N,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Hessian");
      size_type ndof = Z.sizes()[0];
      if (!ndof) { gmm::clear(t.as_vector()); return 0; }
      size_type NN = gmm::sqr(t.sizes().back());
      GA_DEBUG_ASSERT(NN == Z.sizes()[2], "Internal error");
      if (qdim == 1) {
        GA_DEBUG_ASSERT(gmm::vect_size(coeff) == ndof,
                        "Wrong size for coeff vector");
        auto it = Z.begin(); auto itt = t.begin();
        for (size_type kl = 0; kl < NN; ++kl, ++itt) {
          *itt = scalar_type(0);
          for (auto itc = coeff.begin(); itc != coeff.end(); ++itc, ++it)
            *itt += (*itc) * (*it);
        }
        GMM_ASSERT1(itt == t.end(),  "dimensions mismatch");
      } else {
        size_type target_dim = Z.sizes()[1];
        if (target_dim == 1) {
          GA_DEBUG_ASSERT(t.size() == NN*qdim, "dimensions mismatch");
          GA_DEBUG_ASSERT(gmm::vect_size(coeff) == ndof*qdim,
                          "Wrong size for coeff vector");
          gmm::clear(t.as_vector());
          for (size_type q = 0; q < qdim; ++q) {
            base_tensor::const_iterator it = Z.begin();
            for (size_type kl = 0; kl < NN; ++kl)
              for (size_type j = 0; j < ndof; ++j, ++it)
                t[q + kl*qdim] += coeff[j*qdim+q] * (*it);
          }
        } else {
          size_type Qmult = qdim / target_dim;
          GA_DEBUG_ASSERT(t.size() == NN*qdim, "dimensions mismatch");
          GA_DEBUG_ASSERT(gmm::vect_size(coeff) == ndof*Qmult,
                          "Wrong size for coeff vector");
          gmm::clear(t.as_vector());
          for (size_type q = 0; q < Qmult; ++q) {
            base_tensor::const_iterator it = Z.begin();
            for (size_type kl = 0; kl < NN; ++kl)
              for (size_type r = 0; r < target_dim; ++r)
                for (size_type j = 0; j < ndof; ++j, ++it)
                  t[r + q*target_dim + kl*qdim] += coeff[j*Qmult+q] * (*it);
          }
        }
      }
      return 0;
    }

    ga_instruction_hess(base_tensor &tt, const base_tensor &Z_,
                        const base_vector &co, size_type q)
    : ga_instruction_val(tt, Z_, co, q)
    {}
  };

  struct ga_instruction_diverg : public ga_instruction_val {
    // Z(ndof,target_dim,N), coeff(Qmult,ndof) --> t(1)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: divergence");
      size_type ndof = Z.sizes()[0];
      if (!ndof) { gmm::clear(t.as_vector()); return 0; }
      size_type target_dim = Z.sizes()[1];
      size_type N = Z.sizes()[2];
      size_type Qmult = qdim / target_dim;
      GA_DEBUG_ASSERT(Qmult*target_dim == N && (Qmult == 1 || target_dim == 1),
                      "Dimensions mismatch for divergence operator");
      GA_DEBUG_ASSERT(gmm::vect_size(coeff) == ndof*Qmult,
                      "Wrong size for coeff vector");

      t[0] = scalar_type(0);
      base_tensor::const_iterator it = Z.begin();
      if (Qmult == 1)
        for (size_type k = 0; k < N; ++k) {
          if (k) it += (N*ndof + 1);
          for (size_type j = 0; j < ndof; ++j) {
            if (j) ++it;
            t[0] += coeff[j] * (*it);
          }
        }
      else // if (target_dim() == 1)
        for (size_type k = 0; k < N; ++k) {
          if (k) ++it;
          for (size_type j = 0; j < ndof; ++j) {
            if (j) ++it;
            t[0] += coeff[j*N+k] * (*it);
          }
        }
      return 0;
    }

    ga_instruction_diverg(base_tensor &tt, const base_tensor &Z_,
                          const base_vector &co, size_type q)
    : ga_instruction_val(tt, Z_, co, q)
    {}
  };

  struct ga_instruction_copy_val_base : public ga_instruction {
    base_tensor &t;
    const base_tensor &Z;
    size_type qdim;
    // Z(ndof,target_dim) --> t(Qmult*ndof,Qmult*target_dim)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: value of test functions");
      if (qdim == 1) {
        GA_DEBUG_ASSERT(t.size() == Z.size(), "Wrong size for base vector");
        std::copy(Z.begin(), Z.end(), t.begin());
      } else {
        size_type target_dim = Z.sizes()[1];
        size_type Qmult = qdim / target_dim;
        if (Qmult == 1) {
          std::copy(Z.begin(), Z.end(), t.begin());
        } else {
          if (target_dim == 1) {
            size_type ndof = Z.sizes()[0];
            GA_DEBUG_ASSERT(t.size() == Z.size() * Qmult * Qmult,
                            "Wrong size for base vector");
            std::fill(t.begin(), t.end(), scalar_type(0));
            auto itZ = Z.begin();
            size_type s = t.sizes()[0], sss = s+1;

            // Performs t(i*Qmult+j, k*Qmult + j) = Z(i,k);
            auto it = t.begin();
            for (size_type i = 0; i < ndof; ++i, ++itZ) {
              if (i) it += Qmult;
              auto it2 = it;
              *it2 = *itZ;
              for (size_type j = 1; j < Qmult; ++j) { it2 += sss; *it2 = *itZ; }
            }
          } else {
            size_type ndof = Z.sizes()[0];
            GA_DEBUG_ASSERT(t.size() == Z.size() * Qmult * Qmult,
                            "Wrong size for base vector");
            std::fill(t.begin(), t.end(), scalar_type(0));
            auto itZ = Z.begin();
            size_type s = t.sizes()[0], ss = s * Qmult, sss = s+1;

            // Performs t(i*Qmult+j, k*Qmult + j) = Z(i,k);
            for (size_type k = 0; k < target_dim; ++k) {
              auto it = t.begin() + (ss * k);
              for (size_type i = 0; i < ndof; ++i, ++itZ) {
                if (i) it += Qmult;
                auto it2 = it;
                *it2 = *itZ;
                for (size_type j = 1; j < Qmult; ++j)
                  { it2 += sss; *it2 = *itZ; }
              }
            }
          }
        }
      }
      return 0;
    }

    ga_instruction_copy_val_base(base_tensor &tt, const base_tensor &Z_,
                                 size_type q) : t(tt), Z(Z_), qdim(q) {}
  };

  struct ga_instruction_copy_grad_base : public ga_instruction_copy_val_base {
    // Z(ndof,target_dim,N) --> t(Qmult*ndof,Qmult*target_dim,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: gradient of test functions");
      if (qdim == 1) {
        std::copy(Z.begin(), Z.end(), t.begin());
      } else {
        size_type target_dim = Z.sizes()[1];
        size_type Qmult = qdim / target_dim;
        if (Qmult == 1) {
          std::copy(Z.begin(), Z.end(), t.begin());
        } else {
          if (target_dim == 1) {
            size_type ndof = Z.sizes()[0];
            size_type N = Z.sizes()[2];
            GA_DEBUG_ASSERT(t.size() == Z.size() * Qmult * Qmult,
                            "Wrong size for gradient vector");
            std::fill(t.begin(), t.end(), scalar_type(0));
            base_tensor::const_iterator itZ = Z.begin();
            size_type s = t.sizes()[0], sss = s+1, ssss = s*target_dim*Qmult;

            // Performs t(i*Qmult+j, k*Qmult + j, l) = Z(i,k,l);
            for (size_type l = 0; l < N; ++l) {
              base_tensor::iterator it = t.begin() + (ssss*l);
              for (size_type i = 0; i < ndof; ++i, ++itZ) {
                if (i) it += Qmult;
                base_tensor::iterator it2 = it;
                *it2 = *itZ;
                for (size_type j = 1; j < Qmult; ++j) { it2+=sss; *it2=*itZ; }
              }
            }
          } else {
            size_type ndof = Z.sizes()[0];
            size_type N = Z.sizes()[2];
            GA_DEBUG_ASSERT(t.size() == Z.size() * Qmult * Qmult,
                            "Wrong size for gradient vector");
            std::fill(t.begin(), t.end(), scalar_type(0));
            base_tensor::const_iterator itZ = Z.begin();
            size_type s = t.sizes()[0], ss = s * Qmult, sss = s+1;
            size_type ssss = ss*target_dim;

            // Performs t(i*Qmult+j, k*Qmult + j, l) = Z(i,k,l);
            for (size_type l = 0; l < N; ++l)
              for (size_type k = 0; k < target_dim; ++k) {
                base_tensor::iterator it = t.begin() + (ss * k + ssss*l);
                for (size_type i = 0; i < ndof; ++i, ++itZ) {
                  if (i) it += Qmult;
                  base_tensor::iterator it2 = it;
                  *it2 = *itZ;
                  for (size_type j = 1; j < Qmult; ++j) { it2+=sss; *it2=*itZ; }
                }
              }
          }
        }
      }
      return 0;
    }

    ga_instruction_copy_grad_base(base_tensor &tt, const base_tensor &Z_,
                                  size_type q)
      : ga_instruction_copy_val_base(tt,Z_,q) {}
  };

  struct ga_instruction_copy_vect_val_base : public ga_instruction {
    base_tensor &t;
    const base_tensor &Z;
    size_type qdim;
    // Z(ndof) --> t(qdim*ndof,qdim*target_dim)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: vectorized value of test functions");

      size_type ndof = Z.sizes()[0];
      GA_DEBUG_ASSERT(t.size() == Z.size() * qdim * qdim,
                      "Wrong size for base vector");
      // std::fill(t.begin(), t.end(), scalar_type(0)); // Factorized
      auto itZ = Z.begin();
      size_type s = t.sizes()[0], sss = s+1;

      // Performs t(i*qdim+j, k*qdim + j) = Z(i,k);
      auto it = t.begin();
      for (size_type i = 0; i < ndof; ++i, ++itZ) {
        if (i) it += qdim;
        auto it2 = it;
        *it2 = *itZ;
        for (size_type j = 1; j < qdim; ++j) { it2 += sss; *it2 = *itZ; }
      }
      return 0;
    }

    ga_instruction_copy_vect_val_base(base_tensor &tt, const base_tensor &Z_,
                                      size_type q) : t(tt), Z(Z_), qdim(q) {}
  };

  struct ga_instruction_copy_vect_grad_base
    : public ga_instruction_copy_vect_val_base {
    // Z(ndof,N) --> t(qdim*ndof,qdim,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: vectorized gradient of test functions");
      size_type ndof = Z.sizes()[0];
      size_type N = Z.sizes()[2];
      GA_DEBUG_ASSERT(t.size() == Z.size() * qdim * qdim,
                      "Wrong size for gradient vector");
      // std::fill(t.begin(), t.end(), scalar_type(0)); // Factorized
      base_tensor::const_iterator itZ = Z.begin();
      size_type s = t.sizes()[0], sss = s+1, ssss = s*qdim;

      // Performs t(i*qdim+j, k*qdim + j, l) = Z(i,k,l);
      for (size_type l = 0; l < N; ++l) {
        base_tensor::iterator it = t.begin() + (ssss*l);
        for (size_type i = 0; i < ndof; ++i, ++itZ) {
          if (i) it += qdim;
          base_tensor::iterator it2 = it;
          *it2 = *itZ;
          for (size_type j = 1; j < qdim; ++j) { it2+=sss; *it2=*itZ; }
        }
      }
      return 0;
    }

    ga_instruction_copy_vect_grad_base(base_tensor &tt, const base_tensor &Z_,
                                       size_type q)
      : ga_instruction_copy_vect_val_base(tt,Z_,q) {}
  };

  struct ga_instruction_copy_hess_base : public ga_instruction_copy_val_base {
    // Z(ndof,target_dim,N*N) --> t(Qmult*ndof,Qmult*target_dim,N,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Hessian of test functions");
      size_type target_dim = Z.sizes()[1];
      size_type Qmult = qdim / target_dim;
      if (Qmult == 1) {
        gmm::copy(Z.as_vector(), t.as_vector());
      } else {
        size_type ndof = Z.sizes()[0];
        GA_DEBUG_ASSERT(t.size() == Z.size() * Qmult * Qmult,
                      "Wrong size for Hessian vector");
        gmm::clear(t.as_vector());
        base_tensor::const_iterator itZ = Z.begin();
        size_type s = t.sizes()[0], ss = s * Qmult, sss = s+1;

        // Performs t(i*Qmult+j, k*Qmult + j, l, m) = Z(i,k,l*N+m)
        size_type NNdim = Z.sizes()[2]*target_dim;
        for (size_type klm = 0; klm < NNdim; ++klm) {
          base_tensor::iterator it = t.begin() + (ss * klm);
          for (size_type i = 0; i < ndof; ++i, ++itZ) {
            if (i) it += Qmult;
            base_tensor::iterator it2 = it;
            *it2 = *itZ;
            for (size_type j = 1; j < Qmult; ++j) { it2 += sss; *it2 = *itZ; }
          }
        }
      }
      return 0;
    }

    ga_instruction_copy_hess_base(base_tensor &tt, const base_tensor &Z_,
                                  size_type q)
    : ga_instruction_copy_val_base(tt, Z_, q) {}
  };

  struct ga_instruction_copy_diverg_base : public ga_instruction_copy_val_base {
    // Z(ndof,target_dim,N) --> t(Qmult*ndof)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: divergence of test functions");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type N = Z.sizes()[2];
      size_type Qmult = qdim / target_dim;
      GA_DEBUG_ASSERT(Qmult*target_dim == N && (Qmult == 1 || target_dim == 1),
                      "Dimensions mismatch for divergence operator");
      GA_DEBUG_ASSERT(t.size() == ndof * Qmult,
                      "Wrong size for divergence vector");
      gmm::clear(t.as_vector());
      base_tensor::const_iterator itZ = Z.begin();
      if (Qmult == 1) { // target_dim == N
        // Performs t(i) = Trace(Z(i,:,:))
        for (size_type l = 0; l < N; ++l) {
          base_tensor::iterator it = t.begin();
          if (l) itZ += target_dim*ndof+1;
          for (size_type i = 0; i < ndof; ++i) {
            if (i) { ++it; ++itZ; }
            *it += *itZ;
          }
        }
      } else { // Qmult == N
        // Performs t(i*Qmult+j) = Z(i,1,j)
        for (size_type j = 0; j < N; ++j) {
          base_tensor::iterator it = t.begin() + j;
          if (j) ++itZ;
          for (size_type i = 0; i < ndof; ++i) {
            if (i) { it += Qmult; ++itZ; }
            *it += *itZ;
          }
        }
      }
      return 0;
    }

    ga_instruction_copy_diverg_base(base_tensor &tt, const base_tensor &Z_,
                                    size_type q)
      : ga_instruction_copy_val_base(tt, Z_, q) {}
  };

  struct ga_instruction_elementary_trans {
    const base_vector &coeff_in;
    base_vector coeff_out;
    pelementary_transformation elemtrans;
    const mesh_fem &mf1, &mf2;
    const fem_interpolation_context &ctx;
    base_matrix &M;
    size_type &icv;

    void do_transformation(size_type n, size_type m) {
      if (icv != ctx.convex_num() || M.size() == 0) {
        M.base_resize(m, n);
        icv = ctx.convex_num();
        elemtrans->give_transformation(mf1, mf2, icv, M);
      }
      coeff_out.resize(gmm::mat_nrows(M));
      gmm::mult(M, coeff_in, coeff_out); // remember: coeff == coeff_out
    }

    ga_instruction_elementary_trans
    (const base_vector &co, pelementary_transformation e,
     const mesh_fem &mf1_, const mesh_fem &mf2_,
     const fem_interpolation_context &ctx_, base_matrix &M_,
     size_type &icv_)
      : coeff_in(co), elemtrans(e), mf1(mf1_), mf2(mf2_), ctx(ctx_),
        M(M_), icv(icv_) {}
    ~ga_instruction_elementary_trans() {};
  };

  struct ga_instruction_elementary_trans_val
    : public ga_instruction_val, ga_instruction_elementary_trans {
    // Z(ndof,target_dim), coeff_in(Qmult,ndof) --> t(target_dim*Qmult)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: variable value with elementary "
                    "transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      do_transformation(coeff_in.size(), ndof*Qmult);
      return ga_instruction_val::exec();
    }

    ga_instruction_elementary_trans_val
    (base_tensor &tt, const base_tensor &Z_, const base_vector &co, size_type q,
     pelementary_transformation e, const mesh_fem &mf1_, const mesh_fem &mf2_,
     fem_interpolation_context &ctx_, base_matrix &M_, size_type &icv_)
      : ga_instruction_val(tt, Z_, coeff_out, q),
        ga_instruction_elementary_trans(co, e, mf1_, mf2_, ctx_, M_, icv_) {}
  };

  struct ga_instruction_elementary_trans_grad
    : public ga_instruction_grad, ga_instruction_elementary_trans {
    // Z(ndof,target_dim,N), coeff_in(Qmult,ndof) --> t(target_dim*Qmult,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: gradient with elementary transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      do_transformation(coeff_in.size(), ndof*Qmult);
      return ga_instruction_grad::exec();
    }

    ga_instruction_elementary_trans_grad
    (base_tensor &tt, const base_tensor &Z_, const base_vector &co, size_type q,
     pelementary_transformation e, const mesh_fem &mf1_, const mesh_fem &mf2_,
     fem_interpolation_context &ctx_, base_matrix &M_, size_type &icv_)
      : ga_instruction_grad(tt, Z_, coeff_out, q),
        ga_instruction_elementary_trans(co, e, mf1_, mf2_, ctx_, M_, icv_) {}
  };

  struct ga_instruction_elementary_trans_hess
    : public ga_instruction_hess, ga_instruction_elementary_trans {
    // Z(ndof,target_dim,N,N), coeff_in(Qmult,ndof) --> t(target_dim*Qmult,N,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Hessian with elementary transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      do_transformation(coeff_in.size(), ndof*Qmult);
      return ga_instruction_hess::exec();
    }

    ga_instruction_elementary_trans_hess
    (base_tensor &tt, const base_tensor &Z_, const base_vector &co, size_type q,
     pelementary_transformation e, const mesh_fem &mf1_, const mesh_fem &mf2_,
     fem_interpolation_context &ctx_, base_matrix &M_, size_type &icv_)
      : ga_instruction_hess(tt, Z_, coeff_out, q),
        ga_instruction_elementary_trans(co, e, mf1_, mf2_, ctx_, M_, icv_) {}
  };

  struct ga_instruction_elementary_trans_diverg
    : public ga_instruction_diverg, ga_instruction_elementary_trans {
    // Z(ndof,target_dim,N), coeff_in(Qmult,ndof) --> t(1)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: divergence with elementary transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      do_transformation(coeff_in.size(), ndof*Qmult);
      return ga_instruction_diverg::exec();
    }

    ga_instruction_elementary_trans_diverg
    (base_tensor &tt, const base_tensor &Z_, const base_vector &co, size_type q,
     pelementary_transformation e, const mesh_fem &mf1_, const mesh_fem &mf2_,
     fem_interpolation_context &ctx_, base_matrix &M_, size_type &icv_)
      : ga_instruction_diverg(tt, Z_, coeff_out, q),
        ga_instruction_elementary_trans(co, e, mf1_, mf2_, ctx_, M_, icv_) {}
  };

  struct ga_instruction_update_group_info : public ga_instruction {
    const ga_workspace &workspace;
    const ga_instruction_set &gis;
    const ga_instruction_set::interpolate_info &inin;
    const std::string gname;
    ga_instruction_set::variable_group_info &vgi;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Update group info for "+gname);
      if (vgi.cached_mesh && vgi.cached_mesh == inin.m)
        return 0;

      vgi.cached_mesh = inin.m;
      const std::string &varname
        = inin.m ? workspace.variable_in_group(gname, *(inin.m))
                 : workspace.first_variable_of_group(gname);
      vgi.varname = &varname;
      vgi.mf = workspace.associated_mf(varname);
      GA_DEBUG_ASSERT(vgi.mf, "Group variable should always have a mesh_fem");
      vgi.reduced_mf = vgi.mf->is_reduced();
      if (vgi.reduced_mf) {
        const auto it = gis.really_extended_vars.find(varname);
        GA_DEBUG_ASSERT(it != gis.really_extended_vars.end(),
                        "Variable " << varname << " not in extended variables");
        vgi.U = &(it->second);
        vgi.I = &(workspace.temporary_interval_of_variable(varname));
      } else {
        vgi.U = &(workspace.value(varname));
        vgi.I = &(workspace.interval_of_variable(varname));
      }
      vgi.alpha = workspace.factor_of_variable(varname);
      return 0;
    }

    ga_instruction_update_group_info
    (const ga_workspace &workspace_, const ga_instruction_set &gis_,
     const ga_instruction_set::interpolate_info &inin_,
     const std::string &gname_, ga_instruction_set::variable_group_info &vgi_)
      : workspace(workspace_), gis(gis_), inin(inin_), gname(gname_), vgi(vgi_)
    {}
  };

  struct ga_instruction_interpolate_filter : public ga_instruction {
    base_tensor &t;
    const ga_instruction_set::interpolate_info &inin;
    const size_type pt_type;
    const int nb;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: interpolated filter");
      if ((pt_type == size_type(-1) && inin.pt_type) ||
          (pt_type != size_type(-1) && inin.pt_type == pt_type)) {
        GA_DEBUG_INFO("Instruction: interpolated filter: pass");
        return 0;
      }
      else {
        GA_DEBUG_INFO("Instruction: interpolated filter: filtered");
        gmm::clear(t.as_vector());
        return nb;
      }
      return 0;
    }

    ga_instruction_interpolate_filter
    (base_tensor &t_, const ga_instruction_set::interpolate_info &inin_,
     size_type ind_, int nb_)
      : t(t_), inin(inin_), pt_type(ind_), nb(nb_) {}
  };

  struct ga_instruction_copy_interpolated_small_vect : public ga_instruction {
    base_tensor &t;
    const base_small_vector &vec;
    const ga_instruction_set::interpolate_info &inin;
    
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: copy small vector");
      GMM_ASSERT1(!(inin.has_ctx) || inin.ctx.is_convex_num_valid(),
                  "Invalid element, probably transformation failed");
      GMM_ASSERT1(t.size() == vec.size(), "Invalid vector size.");
      gmm::copy(vec, t.as_vector());
      return 0;
    }
    ga_instruction_copy_interpolated_small_vect
    (base_tensor &t_, const base_small_vector &vec_,
     const ga_instruction_set::interpolate_info &inin_)
      : t(t_), vec(vec_), inin(inin_)  {}
  };
  
  struct ga_instruction_interpolate : public ga_instruction {
    base_tensor &t;
    const mesh **m;
    const mesh_fem *mfn, **mfg;
    const base_vector *Un, **Ug;
    fem_interpolation_context &ctx;
    base_vector coeff;
    size_type qdim;
    const size_type &ipt;
    fem_precomp_pool &fp_pool;
    ga_instruction_set::interpolate_info &inin;

    virtual int exec() {
      GMM_ASSERT1(ctx.is_convex_num_valid(), "No valid element for the "
                  "transformation. Probably transformation failed");
      const mesh_fem &mf = *(mfg ? *mfg : mfn);
      const base_vector &U = *(Ug ? *Ug : Un);
      GMM_ASSERT1(&(mf.linked_mesh()) == *m, "Interpolation of a variable "
                  "on another mesh than the one it is defined on");
      slice_vector_on_basic_dof_of_element(mf, U, ctx.convex_num(), coeff);
      pfem pf = mf.fem_of_element(ctx.convex_num());
      GMM_ASSERT1(pf, "Undefined finite element method");
      if (ctx.have_pgp()) {
        if (ipt == 0)
          inin.pfps[&mf] = fp_pool(pf, ctx.pgp()->get_ppoint_tab());
        ctx.set_pfp(inin.pfps[&mf]);
      } else {
        ctx.set_pf(pf);
      }
      return 0;
    }

    ga_instruction_interpolate
    (base_tensor &tt, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_, const base_vector *Un_, const base_vector **Ug_,
     fem_interpolation_context &ctx_, size_type q, const size_type &ipt_,
     fem_precomp_pool &fp_pool_, ga_instruction_set::interpolate_info &inin_)
      : t(tt), m(m_), mfn(mfn_), mfg(mfg_), Un(Un_), Ug(Ug_),
        ctx(ctx_), qdim(q), ipt(ipt_), fp_pool(fp_pool_), inin(inin_) {}
  };

  struct ga_instruction_interpolate_val : public ga_instruction_interpolate {
    // --> t(target_dim*Qmult)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: interpolated variable value");
      ga_instruction_interpolate::exec();
      ctx.pf()->interpolation(ctx, coeff, t.as_vector(), dim_type(qdim));
      // cout << "interpolate " << &U << " result : " << t.as_vector() << endl;
      return 0;
    }

    ga_instruction_interpolate_val
    (base_tensor &tt, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_, const base_vector *Un_, const base_vector **Ug_,
     fem_interpolation_context &ctx_, size_type q, size_type &ipt_,
     fem_precomp_pool &fp_pool_, ga_instruction_set::interpolate_info &inin_)
    : ga_instruction_interpolate(tt, m_, mfn_, mfg_, Un_, Ug_,ctx_, q, ipt_,
                                 fp_pool_, inin_)
    {}
  };

  struct ga_instruction_interpolate_grad : public ga_instruction_interpolate {
    // --> t(target_dim*Qmult,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: interpolated variable grad");
      ga_instruction_interpolate::exec();
      base_matrix v(qdim, ctx.N());
      ctx.pf()->interpolation_grad(ctx, coeff, v, dim_type(qdim));
      gmm::copy(v.as_vector(), t.as_vector());
      return 0;
    }

    ga_instruction_interpolate_grad
    (base_tensor &tt, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_, const base_vector *Un_, const base_vector **Ug_,
     fem_interpolation_context &ctx_, size_type q, size_type &ipt_,
     fem_precomp_pool &fp_pool_, ga_instruction_set::interpolate_info &inin_)
    : ga_instruction_interpolate(tt, m_, mfn_, mfg_, Un_, Ug_, ctx_, q, ipt_,
                                 fp_pool_, inin_)
    {}
  };

  struct ga_instruction_interpolate_hess : public ga_instruction_interpolate {
    // --> t(target_dim*Qmult,N,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: interpolated variable hessian");
      ga_instruction_interpolate::exec();
      base_matrix v(qdim, ctx.N()*ctx.N()); // To be optimized
      ctx.pf()->interpolation_hess(ctx, coeff, v, dim_type(qdim));
      gmm::copy(v.as_vector(), t.as_vector());
      return 0;
    }

    ga_instruction_interpolate_hess
    (base_tensor &tt, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_, const base_vector *Un_, const base_vector **Ug_,
     fem_interpolation_context &ctx_, size_type q, size_type &ipt_,
     fem_precomp_pool &fp_pool_, ga_instruction_set::interpolate_info &inin_)
    : ga_instruction_interpolate(tt, m_, mfn_, mfg_, Un_, Ug_, ctx_, q, ipt_,
                                 fp_pool_, inin_)
    {}
  };

  struct ga_instruction_interpolate_diverg : public ga_instruction_interpolate {
    // --> t(1)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: interpolated variable divergence");
      ga_instruction_interpolate::exec();
      ctx.pf()->interpolation_diverg(ctx, coeff, t[0]);
      return 0;
    }

    ga_instruction_interpolate_diverg
    (base_tensor &tt, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_, const base_vector *Un_, const base_vector **Ug_,
     fem_interpolation_context &ctx_, size_type q, size_type &ipt_,
     fem_precomp_pool &fp_pool_, ga_instruction_set::interpolate_info &inin_)
    : ga_instruction_interpolate(tt, m_, mfn_, mfg_, Un_, Ug_, ctx_, q, ipt_,
                                 fp_pool_, inin_)
    {}
  };

  struct ga_instruction_interpolate_base {
    base_tensor ZZ;
    const mesh **m;
    const mesh_fem *mfn, **mfg;
    const size_type &ipt;
    ga_instruction_set::interpolate_info &inin;
    fem_precomp_pool &fp_pool;

    virtual int exec() {
      GMM_ASSERT1(inin.ctx.is_convex_num_valid(), "No valid element for "
                  "the transformation. Probably transformation failed");
      const mesh_fem &mf = *(mfg ? *mfg : mfn);
      GMM_ASSERT1(&(mf.linked_mesh()) == *m, "Interpolation of a variable "
                  "on another mesh than the one it is defined on");

      pfem pf = mf.fem_of_element(inin.ctx.convex_num());
      GMM_ASSERT1(pf, "Undefined finite element method");

      if (inin.ctx.have_pgp()) {
        if (ipt == 0)
          inin.pfps[&mf] = fp_pool(pf, inin.ctx.pgp()->get_ppoint_tab());
        inin.ctx.set_pfp(inin.pfps[&mf]);
      } else {
        inin.ctx.set_pf(pf);
      }
      return 0;
    }

    ga_instruction_interpolate_base
    (const mesh **m_, const mesh_fem *mfn_, const mesh_fem **mfg_,
     const size_type &ipt_, ga_instruction_set::interpolate_info &inin_,
     fem_precomp_pool &fp_pool_)
      : m(m_), mfn(mfn_), mfg(mfg_), ipt(ipt_), inin(inin_),
        fp_pool(fp_pool_) {}
  };

  struct ga_instruction_interpolate_val_base
    : public ga_instruction_copy_val_base, ga_instruction_interpolate_base {
    // ctx --> Z(ndof,target_dim) --> t(Qmult*ndof,Qmult*target_dim)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: interpolated base value");
      ga_instruction_interpolate_base::exec();
      inin.ctx.pf()->real_base_value(inin.ctx, ZZ); // remember Z == ZZ
      return ga_instruction_copy_val_base::exec();
    }

    ga_instruction_interpolate_val_base
    (base_tensor &t_, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_, const size_type &ipt_, size_type q,
     ga_instruction_set::interpolate_info &inin_, fem_precomp_pool &fp_pool_)
      : ga_instruction_copy_val_base(t_, ZZ, q),
        ga_instruction_interpolate_base(m_, mfn_, mfg_, ipt_,
                                        inin_, fp_pool_) {}
  };

  struct ga_instruction_interpolate_grad_base
    : public ga_instruction_copy_grad_base, ga_instruction_interpolate_base {
    // ctx --> Z(ndof,target_dim,N) --> t(Qmult*ndof,Qmult*target_dim,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: interpolated base grad");
      ga_instruction_interpolate_base::exec();
      inin.ctx.pf()->real_grad_base_value(inin.ctx, ZZ); // remember Z == ZZ
      return ga_instruction_copy_grad_base::exec();
    }

    ga_instruction_interpolate_grad_base
    (base_tensor &t_, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_, const size_type &ipt_, size_type q,
     ga_instruction_set::interpolate_info &inin_, fem_precomp_pool &fp_pool_)
      : ga_instruction_copy_grad_base(t_, ZZ, q),
        ga_instruction_interpolate_base(m_, mfn_, mfg_, ipt_,
                                        inin_, fp_pool_) {}
  };

  struct ga_instruction_interpolate_hess_base
    : public ga_instruction_copy_hess_base, ga_instruction_interpolate_base {
    // ctx --> Z(ndof,target_dim,N*N) --> t(Qmult*ndof,Qmult*target_dim,N,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: interpolated base hessian");
      ga_instruction_interpolate_base::exec();
      inin.ctx.pf()->real_hess_base_value(inin.ctx, ZZ); // remember Z == ZZ
      return ga_instruction_copy_hess_base::exec();
    }

    ga_instruction_interpolate_hess_base
    (base_tensor &t_, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_, const size_type &ipt_, size_type q,
     ga_instruction_set::interpolate_info &inin_, fem_precomp_pool &fp_pool_)
      : ga_instruction_copy_hess_base(t_, ZZ, q),
        ga_instruction_interpolate_base(m_, mfn_, mfg_, ipt_,
                                        inin_, fp_pool_) {}
  };

  struct ga_instruction_interpolate_diverg_base
    : public ga_instruction_copy_diverg_base, ga_instruction_interpolate_base {
    // ctx --> Z(ndof,target_dim,N*N) --> t(Qmult*ndof)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: interpolated base divergence");
      ga_instruction_interpolate_base::exec();
      inin.ctx.pf()->real_grad_base_value(inin.ctx, ZZ); // remember Z == ZZ
      return ga_instruction_copy_diverg_base::exec();
    }

    ga_instruction_interpolate_diverg_base
    (base_tensor &t_, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_, const size_type &ipt_, size_type q,
     ga_instruction_set::interpolate_info &inin_, fem_precomp_pool &fp_pool_)
      : ga_instruction_copy_diverg_base(t_, ZZ, q),
        ga_instruction_interpolate_base(m_, mfn_, mfg_, ipt_,
                                        inin_, fp_pool_) {}
  };


  struct ga_instruction_elementary_trans_base {
    base_tensor t_in;
    base_tensor &t_out;
    pelementary_transformation elemtrans;
    const mesh_fem &mf1, &mf2;
    const fem_interpolation_context &ctx;
    base_matrix &M;
    size_type &icv;

    void do_transformation(size_type n, size_type m) {
      if (icv != ctx.convex_num() || M.size() == 0) {
        M.base_resize(m, n);
        icv = ctx.convex_num();
        elemtrans->give_transformation(mf1, mf2, icv, M);
      }
      t_out.mat_reduction(t_in, M, 0);
    }

    ga_instruction_elementary_trans_base
    (base_tensor &t_, pelementary_transformation e, const mesh_fem &mf1_,
     const mesh_fem &mf2_,
     const fem_interpolation_context &ctx_, base_matrix &M_, size_type &icv_)
      : t_out(t_), elemtrans(e), mf1(mf1_), mf2(mf2_), ctx(ctx_),
        M(M_), icv(icv_) {}
  };

  struct ga_instruction_elementary_trans_val_base
    : public ga_instruction_copy_val_base,
             ga_instruction_elementary_trans_base {
    // Z(ndof,target_dim) --> t_in --> t_out(Qmult*ndof,Qmult*target_dim)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: value of test functions with elementary "
                    "transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      t_in.adjust_sizes(Qmult*ndof, Qmult*Z.sizes()[1]);
      ga_instruction_copy_val_base::exec();
      do_transformation(t_out.sizes()[0], ndof*Qmult);
      return 0;
    }

    ga_instruction_elementary_trans_val_base
    (base_tensor &t_, const base_tensor &Z_, size_type q,
     pelementary_transformation e, const mesh_fem &mf1_, const mesh_fem &mf2_,
     fem_interpolation_context &ctx_, base_matrix &M_, size_type &icv_)
      : ga_instruction_copy_val_base(t_in, Z_, q),
        ga_instruction_elementary_trans_base(t_, e, mf1_, mf2_, ctx_,
                                             M_, icv_) {}
  };

  struct ga_instruction_elementary_trans_grad_base
    : public ga_instruction_copy_grad_base,
             ga_instruction_elementary_trans_base {
    // Z(ndof,target_dim,N) --> t_in --> t_out(Qmult*ndof,Qmult*target_dim,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: gradient of test functions with elementary "
                    "transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      t_in.adjust_sizes(Qmult*ndof, Qmult*Z.sizes()[1], Z.sizes()[2]);
      ga_instruction_copy_grad_base::exec();
      do_transformation(t_out.sizes()[0], ndof*Qmult);
      return 0;
    }

    ga_instruction_elementary_trans_grad_base
    (base_tensor &t_, const base_tensor &Z_, size_type q,
     pelementary_transformation e, const mesh_fem &mf1_, const mesh_fem &mf2_,
     fem_interpolation_context &ctx_, base_matrix &M_, size_type &icv_)
      : ga_instruction_copy_grad_base(t_in, Z_, q),
        ga_instruction_elementary_trans_base(t_, e, mf1_, mf2_, ctx_,
                                             M_,  icv_) {}
  };

  struct ga_instruction_elementary_trans_hess_base
    : public ga_instruction_copy_hess_base,
             ga_instruction_elementary_trans_base {
    // Z(ndof,target_dim,N*N) --> t_out(Qmult*ndof,Qmult*target_dim,N,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Hessian of test functions with elementary "
                    "transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      t_in.adjust_sizes(Qmult*ndof, Qmult*Z.sizes()[1], Z.sizes()[2]);
      ga_instruction_copy_hess_base::exec();
      do_transformation(t_out.sizes()[0], ndof*Qmult);
      return 0;
    }

    ga_instruction_elementary_trans_hess_base
    (base_tensor &t_, const base_tensor &Z_, size_type q,
     pelementary_transformation e, const mesh_fem &mf1_, const mesh_fem &mf2_,
     fem_interpolation_context &ctx_, base_matrix &M_, size_type &icv_)
      : ga_instruction_copy_hess_base(t_in, Z_, q),
        ga_instruction_elementary_trans_base(t_, e, mf1_, mf2_, ctx_,
                                             M_, icv_) {}
  };

  struct ga_instruction_elementary_trans_diverg_base
    : public ga_instruction_copy_diverg_base,
             ga_instruction_elementary_trans_base {
    // Z(ndof,target_dim,N) --> t_out(Qmult*ndof)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: divergence of test functions with elementary "
                    "transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      t_in.adjust_sizes(Qmult*ndof);
      ga_instruction_copy_diverg_base::exec();
      do_transformation(t_out.sizes()[0], ndof*Qmult);
      return 0;
    }

    ga_instruction_elementary_trans_diverg_base
    (base_tensor &t_, const base_tensor &Z_, size_type q,
     pelementary_transformation e, const mesh_fem &mf1_, const mesh_fem &mf2_,
     fem_interpolation_context &ctx_, base_matrix &M_, size_type &icv_)
      : ga_instruction_copy_diverg_base(t_in, Z_, q),
        ga_instruction_elementary_trans_base(t_, e, mf1_, mf2_, ctx_,
                                             M_, icv_) {}
  };


  struct ga_instruction_add : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: addition");
      GA_DEBUG_ASSERT(t.size() == tc1.size(),
                      "internal error " << t.size() << " != " << tc1.size());
      GA_DEBUG_ASSERT(t.size() == tc2.size(),
                      "internal error " << t.size() << " != " << tc2.size());
      gmm::add(tc1.as_vector(), tc2.as_vector(), t.as_vector());
      return 0;
    }
    ga_instruction_add(base_tensor &t_,
                       const base_tensor &tc1_, const base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  struct ga_instruction_add_to : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: addition");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "internal error " << t.size()
                      << " incompatible with " << tc1.size());
      gmm::add(tc1.as_vector(), t.as_vector());
      return 0;
    }
    ga_instruction_add_to(base_tensor &t_, const base_tensor &tc1_)
      : t(t_), tc1(tc1_) {}
  };

  struct ga_instruction_add_to_coeff : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1;
    scalar_type &coeff;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: addition with scale");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "internal error " << t.size()
                      << " incompatible with " << tc1.size());
      gmm::add(gmm::scaled(tc1.as_vector(), coeff), t.as_vector());
      return 0;
    }
    ga_instruction_add_to_coeff(base_tensor &t_, const base_tensor &tc1_,
                                scalar_type &coeff_)
      : t(t_), tc1(tc1_), coeff(coeff_) {}
  };

  struct ga_instruction_sub : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: subtraction");
      GA_DEBUG_ASSERT(t.size() == tc1.size() && t.size() == tc2.size(),
                      "internal error");
      gmm::add(tc1.as_vector(), gmm::scaled(tc2.as_vector(), scalar_type(-1)),
               t.as_vector());
      return 0;
    }
    ga_instruction_sub(base_tensor &t_,
                       const base_tensor &tc1_, const base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  struct ga_instruction_opposite : public ga_instruction {
    base_tensor &t;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: multiplication with -1");
      gmm::scale(t.as_vector(), scalar_type(-1));
      return 0;
    }
    ga_instruction_opposite(base_tensor &t_) : t(t_) {}
  };

  struct ga_instruction_print_tensor : public ga_instruction {
    base_tensor &t;
    pga_tree_node pnode;
    const fem_interpolation_context &ctx;
    size_type &nbpt, &ipt;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: tensor print");
      cout << "Print term "; ga_print_node(pnode, cout);
      cout << " on Gauss point " << ipt << "/" << nbpt << " of element "
           << ctx.convex_num() << ": " << t << endl;
      return 0;
    }
    ga_instruction_print_tensor(base_tensor &t_, pga_tree_node pnode_,
                                const fem_interpolation_context &ctx_,
                                size_type &nbpt_, size_type &ipt_)
      : t(t_), pnode(pnode_), ctx(ctx_), nbpt(nbpt_), ipt(ipt_) {}
  };

  struct ga_instruction_copy_tensor : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: tensor copy");
      std::copy(tc1.begin(), tc1.end(), t.begin());
      // gmm::copy(tc1.as_vector(), t.as_vector());
      return 0;
    }
    ga_instruction_copy_tensor(base_tensor &t_, const base_tensor &tc1_)
      : t(t_), tc1(tc1_) {}
  };

  struct ga_instruction_clear_tensor : public ga_instruction {
    base_tensor &t;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: clear tensor");
      std::fill(t.begin(), t.end(), scalar_type(0));
      return 0;
    }
    ga_instruction_clear_tensor(base_tensor &t_) : t(t_) {}
  };

  struct ga_instruction_copy_tensor_possibly_void : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: tensor copy possibly void");
      if (tc1.size())
        gmm::copy(tc1.as_vector(), t.as_vector());
      else
        gmm::clear(t.as_vector());
      return 0;
    }
    ga_instruction_copy_tensor_possibly_void(base_tensor &t_,
                                             const base_tensor &tc1_)
      : t(t_), tc1(tc1_) {}
  };

  struct ga_instruction_copy_scalar : public ga_instruction {
    scalar_type &t; const scalar_type &t1;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: scalar copy");
      t = t1;
      return 0;
    }
    ga_instruction_copy_scalar(scalar_type &t_, const scalar_type &t1_)
      : t(t_), t1(t1_) {}
  };

  struct ga_instruction_copy_vect : public ga_instruction {
    base_vector &t;
    const base_vector &t1;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: fixed size tensor copy");
      gmm::copy(t1, t);
      return 0;
    }
    ga_instruction_copy_vect(base_vector &t_, const base_vector &t1_)
      : t(t_), t1(t1_) {}
  };

  struct ga_instruction_trace : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1;
    size_type n;
    // tc1(:,:,...,n,n) --> t(:,:,...)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Trace");
      GA_DEBUG_ASSERT(t.size()*n*n == tc1.size(), "Wrong sizes");
      size_type s = t.size() * (n+1);
      auto it = t.begin();
      auto it1 = tc1.begin();
      for (; it != t.end(); ++it, ++it1) {
        auto it2 = it1;
        *it = *it2;
        for (size_type i = 1; i < n; ++i) { it2 += s; *it += *it2; }
      }
      return 0;
    }

    ga_instruction_trace(base_tensor &t_, const base_tensor &tc1_, size_type n_)
      : t(t_), tc1(tc1_), n(n_) {}
  };

  struct ga_instruction_deviator : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1;
    size_type n;
    // tc1(:,:,...,n,n) --> t(:,:,...,n,n)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Deviator");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");

      gmm::copy(tc1.as_vector(), t.as_vector());

      size_type nb = t.size()/(n*n);
      size_type s = nb * (n+1), j = 0;
      base_tensor::iterator it = t.begin();
      base_tensor::const_iterator it1 = tc1.begin();
      for (; j < nb; ++it, ++it1, ++j) {
        scalar_type tr(0);
        base_tensor::const_iterator it2 = it1;
        tr += *it2;
        for (size_type i = 1; i < n; ++i) { it2 += s; tr += *it2; }
        tr /= scalar_type(n);

        base_tensor::iterator it3 = it;
        *it3 -= tr;
        for (size_type i = 1; i < n; ++i) { it3 += s; *it3 -= tr; }
      }
      return 0;
    }

    ga_instruction_deviator(base_tensor &t_, const base_tensor &tc1_,
                            size_type n_)
      : t(t_), tc1(tc1_), n(n_) {}
  };

  struct ga_instruction_transpose : public ga_instruction { // To be optimized
    base_tensor &t;
    const base_tensor &tc1;
    size_type n1, n2, nn;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: transpose");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");

      size_type n0 = tc1.size() / (n1*n2*nn);
      auto it = t.begin();
      for (size_type i = 0; i < nn; ++i) {
        size_type s1 = i*n1*n2*n0;
        for (size_type j = 0; j < n1; ++j) {
          size_type s2 = s1 + j*n0;
          for (size_type k = 0; k < n2; ++k) {
            size_type s3 = s2 + k*n1*n0;
            for (size_type l = 0; l < n0; ++l, ++it)
              *it = tc1[s3+l];
          }
        }
      }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_transpose(base_tensor &t_, const base_tensor &tc1_,
                             size_type n1_, size_type n2_, size_type nn_)
      : t(t_), tc1(tc1_), n1(n1_), n2(n2_), nn(nn_) {}
  };

  struct ga_instruction_swap_indices : public ga_instruction {// To be optimized
    base_tensor &t;
    const base_tensor &tc1;
    size_type nn1, nn2, ii2, ii3;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: swap indices");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      size_type ii1 = t.size() / (nn1*nn2*ii2*ii3);

      auto it = t.begin();
      for (size_type i = 0; i < ii3; ++i)
        for (size_type j = 0; j < nn1; ++j)
          for (size_type k = 0; k < ii2; ++k)
            for (size_type l = 0; l < nn2; ++l) {
              size_type ind = j*ii1+k*ii1*nn1+l*ii1*nn1*ii2+i*ii1*nn1*ii2*nn2;
              for (size_type m = 0; m < ii1; ++m, ++it)
                *it = tc1[m+ind];
            }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_swap_indices(base_tensor &t_, const base_tensor &tc1_,
                                size_type n1_, size_type n2_,
                                size_type i2_, size_type i3_)
      : t(t_), tc1(tc1_), nn1(n1_), nn2(n2_), ii2(i2_), ii3(i3_) {}
  };

  struct ga_instruction_index_move_last : public ga_instruction {// To be optimized
    base_tensor &t;
    const base_tensor &tc1;
    size_type nn, ii2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: swap indices");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      size_type ii1 = t.size() / (nn*ii2);

      auto it = t.begin();
      for (size_type i = 0; i < nn; ++i)
        for (size_type j = 0; j < ii2; ++j) {
          size_type ind = i*ii1+j*ii1*nn;
          for (size_type k = 0; k < ii1; ++k, ++it)
            *it = tc1[k+ind];
        }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_index_move_last(base_tensor &t_, const base_tensor &tc1_,
                                   size_type n_, size_type i2_)
      : t(t_), tc1(tc1_), nn(n_), ii2(i2_) {}
  };

  struct ga_instruction_transpose_no_test : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1;
    size_type n1, n2, nn;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: transpose");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");

      auto it = t.begin();
      for (size_type i = 0; i < nn; ++i) {
        size_type s1 = i*n1*n2;
        for (size_type j = 0; j < n1; ++j) {
          size_type s2 = s1 + j;
          for (size_type k = 0; k < n2; ++k, ++it)
            *it = tc1[s2 + k*n1];
        }
      }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_transpose_no_test(base_tensor &t_, const base_tensor &tc1_,
                                     size_type n1_, size_type n2_,
                                     size_type nn_)
      : t(t_), tc1(tc1_), n1(n1_), n2(n2_), nn(nn_) {}
  };

  struct ga_instruction_transpose_test : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: copy tensor and transpose test functions");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      GA_DEBUG_ASSERT(t.sizes().size() >= 2, "Wrong sizes");

      size_type s1 = t.sizes()[0], s2 = t.sizes()[1], s3 = s1*s2;
      size_type s = t.size() / s3;
      base_tensor::iterator it = t.begin();
      for (size_type k = 0; k < s; ++k)
        for (size_type j = 0; j < s2;  ++j)
          for (size_type i = 0; i < s1; ++i, ++it)
            *it = tc1[j+s2*i+k*s3];
      return 0;
    }
    ga_instruction_transpose_test(base_tensor &t_, const base_tensor &tc1_)
      : t(t_), tc1(tc1_) {}
  };

  struct ga_instruction_sym : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: symmetric part");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      size_type order = t.sizes().size();
      size_type s1 = t.sizes()[order-2], s2 = t.sizes()[order-1];
      size_type s = t.size() / (s1*s2);
      for (size_type i = 0; i < s1;  ++i)
        for (size_type j = 0; j < s2;  ++j) {
          base_tensor::iterator it = t.begin() + s*(i + s1*j);
          base_tensor::const_iterator it1 = tc1.begin() + s*(i + s1*j),
                                      it1T = tc1.begin() + s*(j + s2*i);
          for (size_type k = 0; k < s; ++k) *it++ = 0.5*(*it1++ + *it1T++);
        }
      return 0;
    }
    ga_instruction_sym(base_tensor &t_, const base_tensor &tc1_)
      : t(t_), tc1(tc1_) {}
  };

  struct ga_instruction_skew : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: skew-symmetric part");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      size_type order = t.sizes().size();
      size_type s1 = t.sizes()[order-2], s2 = t.sizes()[order-1];
      size_type s = t.size() / (s1*s2);
      for (size_type i = 0; i < s1;  ++i)
        for (size_type j = 0; j < s2;  ++j) {
          base_tensor::iterator it = t.begin() + s*(i + s1*j);
          base_tensor::const_iterator it1 = tc1.begin() + s*(i + s1*j),
                                      it1T = tc1.begin() + s*(j + s2*i);
          for (size_type k = 0; k < s; ++k) *it++ = 0.5*(*it1++ - *it1T++);
        }
      return 0;
    }
    ga_instruction_skew(base_tensor &t_, const base_tensor &tc1_)
      : t(t_), tc1(tc1_) {}
  };

  struct ga_instruction_scalar_add : public ga_instruction {
    scalar_type &t;
    const scalar_type &c, &d;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: scalar addition");
      t = c + d;
      return 0;
    }
    ga_instruction_scalar_add(scalar_type &t_, const scalar_type &c_,
                              const  scalar_type &d_)
      : t(t_), c(c_), d(d_) {}
  };

  struct ga_instruction_scalar_sub : public ga_instruction {
    scalar_type &t;
    const scalar_type &c, &d;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: scalar subtraction");
      t = c - d;
      return 0;
    }
    ga_instruction_scalar_sub(scalar_type &t_, const scalar_type &c_,
                              const  scalar_type &d_)
      : t(t_), c(c_), d(d_) {}
  };

  struct ga_instruction_scalar_scalar_mult : public ga_instruction {
    scalar_type &t;
    const scalar_type &c, &d;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: scalar multiplication");
      t = c * d;
      return 0;
    }
    ga_instruction_scalar_scalar_mult(scalar_type &t_, const scalar_type &c_,
                                      const  scalar_type &d_)
      : t(t_), c(c_), d(d_) {}
  };

  struct ga_instruction_scalar_scalar_div : public ga_instruction {
    scalar_type &t;
    const scalar_type &c, &d;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: scalar division");
      t = c / d;
      return 0;
    }
    ga_instruction_scalar_scalar_div(scalar_type &t_, const scalar_type &c_,
                                     const  scalar_type &d_)
      : t(t_), c(c_), d(d_) {}
  };

  struct ga_instruction_scalar_mult : public ga_instruction {
    base_tensor &t, &tc1;
    const scalar_type &c;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: multiplication of a tensor by a scalar " << c);
      gmm::copy(gmm::scaled(tc1.as_vector(), c), t.as_vector());
      return 0;
    }
    ga_instruction_scalar_mult(base_tensor &t_, base_tensor &tc1_,
                               const scalar_type &c_)
      : t(t_), tc1(tc1_), c(c_) {}
  };

  struct ga_instruction_scalar_div : public ga_instruction {
    base_tensor &t, &tc1;
    const scalar_type &c;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: division of a tensor by a scalar");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");

      base_tensor::iterator it = t.begin(), it1 = tc1.begin();
      for (; it != t.end(); ++it, ++it1) *it = *it1/c;
      return 0;
    }
    ga_instruction_scalar_div(base_tensor &t_, base_tensor &tc1_,
                               const scalar_type &c_)
      : t(t_), tc1(tc1_), c(c_) {}
  };

  // Performs Cross product in the presence of test functions
  struct ga_instruction_cross_product_tf : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    bool inv;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Cross product with test functions");

      size_type n1 = tc1.size() / 3, n2 =  tc2.size() / 3, nn=n1*n2;
      GA_DEBUG_ASSERT(t.size() == nn*3, "Bad tensor size for cross product");
      size_type mm=2*nn, n1_2 = 2*n1, n2_2 = 2*n2;
      base_tensor::iterator it = t.begin(), it2 = tc2.begin();

      if (inv) {
        for (size_type i = 0; i < n2; ++i, ++it2) {
          base_tensor::iterator it1 = tc1.begin();
          for (size_type j = 0; j < n1; ++j, ++it, ++it1) {
            *it    = - it1[n1]  *it2[n2_2] + it1[n1_2]*it2[n2];
            it[nn] = - it1[n1_2]*it2[0]    + it1[0]   *it2[n2_2];
            it[mm] = - it1[0]   *it2[n2]   + it1[n1]  *it2[0];
          }
        }
      } else {
        for (size_type i = 0; i < n2; ++i, ++it2) {
          base_tensor::iterator it1 = tc1.begin();
          for (size_type j = 0; j < n1; ++j, ++it, ++it1) {
            *it    = it1[n1]  *it2[n2_2] - it1[n1_2]*it2[n2];
            it[nn] = it1[n1_2]*it2[0]    - it1[0]   *it2[n2_2];
            it[mm] = it1[0]   *it2[n2]   - it1[n1]  *it2[0];
          }
        }
      } 
      return 0;
    }
    ga_instruction_cross_product_tf(base_tensor &t_, base_tensor &tc1_,
                                    base_tensor &tc2_, bool inv_)
      : t(t_), tc1(tc1_), tc2(tc2_), inv(inv_) {}
  };

  // Performs Cross product in the absence of test functions
  struct ga_instruction_cross_product : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Cross product with test functions");
      GA_DEBUG_ASSERT(t.size() == 3 && tc1.size() == 3 && tc2.size() == 3,
                       "Bad tensor size for cross product");
      t[0] = tc1[1]*tc2[2] - tc1[2]*tc2[1];
      t[1] = tc1[2]*tc2[0] - tc1[0]*tc2[2];
      t[2] = tc1[0]*tc2[1] - tc1[1]*tc2[0];
      return 0;
    }
    ga_instruction_cross_product(base_tensor &t_, base_tensor &tc1_,
                                 base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };



  
  struct ga_instruction_dotmult : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: componentwise multiplication");
      size_type s2 = tc2.size(), s1_1 = tc1.size() / s2;
      GA_DEBUG_ASSERT(t.size() == s1_1*s2, "Wrong sizes");

      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < s2; ++i)
        for (size_type m = 0; m < s1_1; ++m, ++it)
          *it = tc1[m+s1_1*i] * tc2[i];
      return 0;
    }
    ga_instruction_dotmult(base_tensor &t_, base_tensor &tc1_,
                           base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  struct ga_instruction_dotdiv : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: componentwise division");
      size_type s2 = tc2.size(), s1_1 = tc1.size() / s2;
      GA_DEBUG_ASSERT(t.size() == s1_1*s2, "Wrong sizes");

      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < s2; ++i)
        for (size_type m = 0; m < s1_1; ++m, ++it)
          *it = tc1[m+s1_1*i] / tc2[i];
      return 0;
    }
    ga_instruction_dotdiv(base_tensor &t_, base_tensor &tc1_,
                          base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  // Performs Ami Bni -> Cmni
  struct ga_instruction_dotmult_spec : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: specific componentwise multiplication");
      size_type s2_1 = tc2.sizes()[0], s2_2 = tc2.size() / s2_1;
      size_type s1_1 = tc1.size() / s2_2;

      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < s2_2; ++i)
        for (size_type n = 0; n < s2_1; ++n)
          for (size_type m = 0; m < s1_1; ++m, ++it)
            *it = tc1[m+s1_1*i] * tc2[n+s2_1*i];
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_dotmult_spec(base_tensor &t_, base_tensor &tc1_,
                           base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  // Performs Amijik -> Cmjk. To be optimized
  struct ga_instruction_contract_1_1 : public ga_instruction {
    base_tensor &t, &tc1;
    size_type nn, ii2, ii3;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: single contraction on a single tensor");

      size_type ii1 = tc1.size() / (nn*nn*ii2*ii3);

      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < ii3; ++i)
        for (size_type j = 0; j < ii2; ++j)
          for (size_type k = 0; k < ii1; ++k, ++it) {
            *it = scalar_type(0);
            size_type pre_ind = k+j*ii1*nn+i*ii1*nn*ii2*nn;
            for (size_type n = 0; n < nn; ++n)
              *it += tc1[pre_ind+n*ii1+n*ii1*nn*ii2];
          }

      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_contract_1_1(base_tensor &t_, base_tensor &tc1_,
                                size_type n_, size_type i2_, size_type i3_)
      : t(t_), tc1(tc1_), nn(n_), ii2(i2_), ii3(i3_)  {}
  };

  // Performs Amijk Bnljp -> Cmniklp. To be optimized
  struct ga_instruction_contract_2_1 : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn, ii1, ii2, ii3, ii4;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: single contraction on two tensors");

      size_type ift1 = tc1.size() / (nn*ii1*ii2);
      size_type ift2 = tc2.size() / (nn*ii3*ii4);

      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < ii4; ++i)
        for (size_type j = 0; j < ii3; ++j)
          for (size_type k = 0; k < ii2; ++k)
            for (size_type l = 0; l < ii1; ++l)
              for (size_type p = 0; p < ift2; ++p)
                for (size_type q = 0; q < ift1; ++q, ++it) {
                  *it = scalar_type(0);
                  size_type ind1 = q+l*ift1+k*ift1*ii1*nn;
                  size_type ind2 = p+j*ift2+i*ift2*ii3*nn;
                  for (size_type n = 0; n < nn; ++n)
                    *it += tc1[ind1+n*ift1*ii1] * tc2[ind2+n*ift2*ii3];
                }

      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_contract_2_1(base_tensor &t_, base_tensor &tc1_,
                                base_tensor &tc2_,
                                size_type n_, size_type i1_, size_type i2_,
                                size_type i3_, size_type i4_)
      : t(t_), tc1(tc1_), tc2(tc2_), nn(n_),
        ii1(i1_), ii2(i2_), ii3(i3_), ii4(i4_)  {}
  };

  // Performs Amijk Bnljp -> Cnmiklp. To be optimized
  struct ga_instruction_contract_2_1_rev : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn, ii1, ii2, ii3, ii4;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: single contraction on two tensors");

      size_type ift1 = tc1.size() / (nn*ii1*ii2);
      size_type ift2 = tc2.size() / (nn*ii3*ii4);

      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < ii4; ++i)
        for (size_type j = 0; j < ii3; ++j)
          for (size_type k = 0; k < ii2; ++k)
            for (size_type l = 0; l < ii1; ++l)
              for (size_type q = 0; q < ift1; ++q)
                for (size_type p = 0; p < ift2; ++p, ++it) {
                  *it = scalar_type(0);
                  size_type ind1 = q+l*ift1+k*ift1*ii1*nn;
                  size_type ind2 = p+j*ift2+i*ift2*ii3*nn;
                  for (size_type n = 0; n < nn; ++n)
                    *it += tc1[ind1+n*ift1*ii1] * tc2[ind2+n*ift2*ii3];
                }

      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_contract_2_1_rev(base_tensor &t_, base_tensor &tc1_,
                                    base_tensor &tc2_,
                                    size_type n_, size_type i1_, size_type i2_,
                                    size_type i3_, size_type i4_)
      : t(t_), tc1(tc1_), tc2(tc2_), nn(n_),
        ii1(i1_), ii2(i2_), ii3(i3_), ii4(i4_)  {}
  };

  // Performs Amijklp Bnqjrls -> Cmnikpqrs. To be optimized
  struct ga_instruction_contract_2_2 : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn1, nn2, ii1, ii2, ii3, ii4, ii5, ii6;
    bool inv_tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: single contraction on two tensors");

      size_type ift1 = tc1.size() / (nn1*nn2*ii1*ii2*ii3);
      size_type ift2 = tc2.size() / (nn1*nn2*ii3*ii4*ii5);

      size_type sn1 = ift2*ii4, sn2 = ift2*ii4*nn1*ii5;
      if (inv_tc2) std::swap(sn1, sn2);

      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < ii6; ++i)
        for (size_type j = 0; j < ii5; ++j)
          for (size_type k = 0; k < ii4; ++k)
            for (size_type l = 0; l < ii3; ++l)
              for (size_type p = 0; p < ii2; ++p)
                for (size_type q = 0; q < ii1; ++q)
                  for (size_type r = 0; r < ift2; ++r)
                    for (size_type s = 0; s < ift1; ++s, ++it) {
                      *it = scalar_type(0);
                      size_type ind1
                        = s+q*ift1+p*ift1*ii1*nn1+l*ift1*ii1*nn1*ii2*nn2;
                      size_type ind2
                        = r+k*ift2+j*ift2*ii4*nn1+i*ift2*ii4*nn1*ii5*nn2;
                      for (size_type n1 = 0; n1 < nn1; ++n1)
                        for (size_type n2 = 0; n2 < nn2; ++n2)
                          *it += tc1[ind1+n1*ift1*ii1+n2*ift1*ii1*nn1*ii2]
                            * tc2[ind2+n1*sn1+n2*sn2];
                    }

      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_contract_2_2(base_tensor &t_, base_tensor &tc1_,
                                base_tensor &tc2_,
                                size_type n1_, size_type n2_,
                                size_type i1_, size_type i2_, size_type i3_,
                                size_type i4_, size_type i5_, size_type i6_,
                                bool intc2)
      : t(t_), tc1(tc1_), tc2(tc2_), nn1(n1_), nn2(n2_),
        ii1(i1_), ii2(i2_), ii3(i3_), ii4(i4_), ii5(i5_), ii6(i6_),
        inv_tc2(intc2)  {}
  };

  // Performs Amijklp Bnqjrls -> Cnmikpqrs. To be optimized
  struct ga_instruction_contract_2_2_rev : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn1, nn2, ii1, ii2, ii3, ii4, ii5, ii6;
    bool inv_tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: single contraction on two tensors");

      size_type ift1 = tc1.size() / (nn1*nn2*ii1*ii2*ii3);
      size_type ift2 = tc2.size() / (nn1*nn2*ii3*ii4*ii5);

      size_type sn1 = ift2*ii4, sn2 = ift2*ii4*nn1*ii5;
      if (inv_tc2) std::swap(sn1, sn2);

      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < ii6; ++i)
        for (size_type j = 0; j < ii5; ++j)
          for (size_type k = 0; k < ii4; ++k)
            for (size_type l = 0; l < ii3; ++l)
              for (size_type p = 0; p < ii2; ++p)
                for (size_type q = 0; q < ii1; ++q)
                  for (size_type s = 0; s < ift1; ++s)
                    for (size_type r = 0; r < ift2; ++r, ++it) {
                      *it = scalar_type(0);
                      size_type ind1
                        = s+q*ift1+p*ift1*ii1*nn1+l*ift1*ii1*nn1*ii2*nn2;
                      size_type ind2
                        = r+k*ift2+j*ift2*ii4*nn1+i*ift2*ii4*nn1*ii5*nn2;
                      for (size_type n1 = 0; n1 < nn1; ++n1)
                        for (size_type n2 = 0; n2 < nn2; ++n2)
                          *it += tc1[ind1+n1*ift1*ii1+n2*ift1*ii1*nn1*ii2]
                            * tc2[ind2+n1*sn1+n2*sn2];
                    }

      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_contract_2_2_rev(base_tensor &t_, base_tensor &tc1_,
                                    base_tensor &tc2_,
                                    size_type n1_, size_type n2_,
                                    size_type i1_, size_type i2_, size_type i3_,
                                    size_type i4_, size_type i5_, size_type i6_,
                                    bool intc2)
      : t(t_), tc1(tc1_), tc2(tc2_), nn1(n1_), nn2(n2_),
        ii1(i1_), ii2(i2_), ii3(i3_), ii4(i4_), ii5(i5_), ii6(i6_),
        inv_tc2(intc2)  {}
  };


  // Performs Amj Bjk -> Cmk. To be optimized
  struct ga_instruction_matrix_mult : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type n;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: order one contraction "
                    "(dot product or matrix multiplication)");

      size_type s1 = tc1.size() / n;
      size_type s2 = tc2.size() / n;

      base_tensor::iterator it = t.begin();
      for (size_type k = 0; k < s2; ++k)
        for (size_type i = 0; i < s1; ++i, ++it) {
            *it = scalar_type(0);
            for (size_type j = 0; j < n; ++j)
              *it += tc1[i+j*s1] * tc2[j+k*n];
          }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_matrix_mult(base_tensor &t_, base_tensor &tc1_,
                               base_tensor &tc2_, size_type n_)
      : t(t_), tc1(tc1_), tc2(tc2_), n(n_) {}
  };

  // Performs Amij Bnjk -> Cmnik. To be optimized
  struct ga_instruction_matrix_mult_spec : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type n, m, p; // tc1 of size q*m*n, tc2 of size l*n*p
                       // t of size q*l*m*p
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: specific order one contraction "
                    "(dot product or matrix multiplication)");
      size_type q = tc1.size() / (m * n);
      size_type l = tc2.size() / (p * n);

      base_tensor::iterator it = t.begin();
      for (size_type r = 0; r < p; ++r)
        for (size_type k = 0; k < m; ++k)
          for (size_type j = 0; j < l; ++j)
            for (size_type i = 0; i < q; ++i, ++it) {
              *it = scalar_type(0);
              for (size_type s = 0; s < n; ++s)
                *it += tc1[i+k*q+s*q*m] * tc2[j+s*l+r*l*n];
            }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_matrix_mult_spec(base_tensor &t_, base_tensor &tc1_,
                                    base_tensor &tc2_, size_type n_,
                                    size_type m_, size_type p_)
      : t(t_), tc1(tc1_), tc2(tc2_), n(n_), m(m_), p(p_) {}
  };

  // Performs Amij Bnjk -> Cnmik. To be optimized
  struct ga_instruction_matrix_mult_spec2 : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type n, m, p; // tc1 of size q*m*n, tc2 of size l*n*p
                       // t of size l*q*m*p
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: specific order one contraction "
                    "(dot product or matrix multiplication)");
      size_type q = tc1.size() / (m * n);
      size_type l = tc2.size() / (p * n);

      base_tensor::iterator it = t.begin();
      for (size_type r = 0; r < p; ++r)
        for (size_type k = 0; k < m; ++k)
          for (size_type i = 0; i < q; ++i)
            for (size_type j = 0; j < l; ++j, ++it) {
              *it = scalar_type(0);
              for (size_type s = 0; s < n; ++s)
                *it += tc1[i+k*q+s*q*m] * tc2[j+s*l+r*l*n];
            }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_matrix_mult_spec2(base_tensor &t_, base_tensor &tc1_,
                                     base_tensor &tc2_, size_type n_,
                                     size_type m_, size_type p_)
      : t(t_), tc1(tc1_), tc2(tc2_), n(n_), m(m_), p(p_) {}
  };

  // Performs Ani Bmi -> Cmn
  struct ga_instruction_contraction : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: contraction operation of size " << nn);
#if GA_USES_BLAS
      long m = int(tc1.size()/nn), k = int(nn), n = int(tc2.size()/nn);
      long lda = m, ldb = n, ldc = m;
      char T = 'T', N = 'N';
      scalar_type alpha(1), beta(0);
      gmm::dgemm_(&N, &T, &m, &n, &k, &alpha, &(tc1[0]), &lda, &(tc2[0]), &ldb,
                  &beta, &(t[0]), &ldc);
#else
      size_type s1 = tc1.size()/nn, s2 = tc2.size()/nn;
      GA_DEBUG_ASSERT(t.size() == s1*s2, "Internal error");

      auto it1=tc1.begin(), it2=tc2.begin(), it2end=it2 + s2;
      for (auto it = t.begin(); it != t.end(); ++it) {
        auto it11 = it1, it22 = it2;
        scalar_type a = (*it11) * (*it22);
        for (size_type i = 1; i < nn; ++i)
          { it11 += s1; it22 += s2; a += (*it11) * (*it22); }
        *it = a;
        ++it2; if (it2 == it2end) { it2 = tc2.begin(), ++it1; }
      }
      // auto it = t.begin(); // Unoptimized version.
      // for (size_type i = 0; i < s1; ++i)
      //   for (size_type j = 0; j < s2; ++j, ++it) {
      //     *it = scalar_type(0);
      //     for (size_type k = 0; k < nn; ++k)
      //       *it += tc1[i+k*s1] * tc2[j+k*s2];
      //   }
#endif
      return 0;
    }
    ga_instruction_contraction(base_tensor &t_, base_tensor &tc1_,
                             base_tensor &tc2_, size_type n_)
      : t(t_), tc1(tc1_), tc2(tc2_), nn(n_) {}
  };

  // Performs Ani Bmi -> Cmn
  struct ga_instruction_contraction_opt0_2 : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type n, q;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: contraction operation of size " << n*q <<
                    " optimized for vectorized second tensor of type 2");
      size_type nn = n*q, s1 = tc1.size()/nn, s2 = tc2.size()/nn, s2_q = s2/q;
      size_type s1_qq = s1*q, s2_qq = s2*q;
      GA_DEBUG_ASSERT(t.size() == s1*s2, "Internal error");

      auto it = t.begin(), it1 = tc1.begin();
      for (size_type i = 0; i < s1; ++i, ++it1) {
        auto it2 = tc2.begin();
        for (size_type j = 0; j < s2_q; ++j) {
          if (j) it2+=q;
          auto itt1 = it1;
          for (size_type l = 0; l < q; ++l, ++it) {
            if (l) itt1 += s1;
            auto ittt1 = itt1, ittt2 = it2;
            *it = *ittt1 * (*ittt2);
            for (size_type m = 1; m < n; ++m) {
              ittt1 += s1_qq, ittt2 += s2_qq; *it += *ittt1 * (*ittt2);
            }
          }
        }
      }
      // base_tensor u = t;
      // ga_instruction_contraction toto(t, tc1, tc2, n*q);
      // toto.exec();
      // GMM_ASSERT1(gmm::vect_dist2(t.as_vector(), u.as_vector()) < 1E-9, "Erroneous");
      return 0;
    }
    ga_instruction_contraction_opt0_2(base_tensor &t_, base_tensor &tc1_,
                                    base_tensor &tc2_, size_type n_,
                                    size_type q_)
      : t(t_), tc1(tc1_), tc2(tc2_), n(n_), q(q_) {}
  };

  // Performs Ani Bmi -> Cmn
  template <int N>
  struct ga_instruction_contraction_opt0_2_unrolled : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type q;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: unrolled contraction of size " << N*q <<
                    " optimized for vectorized second tensor of type 2");
      size_type nn = N*q, s1 = tc1.size()/nn, s2 = tc2.size()/nn, s2_q = s2/q;
      size_type s1_qq = s1*q, s2_qq = s2*q;
      GA_DEBUG_ASSERT(t.size() == s1*s2, "Internal error");

      auto it = t.begin(), it1 = tc1.begin();
      for (size_type i = 0; i < s1; ++i, ++it1) {
        auto it2 = tc2.begin();
        for (size_type j = 0; j < s2_q; ++j) {
          if (j) it2+=q;
          auto itt1 = it1;
          for (size_type l = 0; l < q; ++l, ++it) {
            if (l) itt1 += s1;
            auto ittt1 = itt1, ittt2 = it2;
            *it = *ittt1 * (*ittt2);
            for (size_type m = 1; m < N; ++m) {
              ittt1 += s1_qq, ittt2 += s2_qq; *it += *ittt1 * (*ittt2);
            }
          }
        }
      }
      return 0;
    }
    ga_instruction_contraction_opt0_2_unrolled(base_tensor &t_, base_tensor &tc1_,
                                             base_tensor &tc2_, size_type q_)
      : t(t_), tc1(tc1_), tc2(tc2_), q(q_) {}
  };

  // Performs Ani Bmi -> Cmn
  template <int N, int Q>
  struct ga_instruction_contraction_opt0_2_dunrolled : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: unrolled contraction of size " << N*Q
                    << " optimized for vectorized second tensor of type 2");
      size_type s1 = tc1.size()/(N*Q), s2 = tc2.size()/(N*Q), s2_q = s2/Q;
      size_type s1_qq = s1*Q, s2_qq = s2*Q;
      GA_DEBUG_ASSERT(t.size() == s1*s2, "Internal error");

      auto it = t.begin(), it1 = tc1.begin();
      for (size_type i = 0; i < s1; ++i, ++it1) {
        auto it2 = tc2.begin();
        for (size_type j = 0; j < s2_q; ++j) {
          if (j) it2+=Q;
          auto itt1 = it1;
          for (size_type l = 0; l < Q; ++l, ++it) {
            if (l) itt1 += s1;
            auto ittt1 = itt1, ittt2 = it2;
            *it = *ittt1 * (*ittt2);
            for (size_type m = 1; m < N; ++m) {
              ittt1 += s1_qq, ittt2 += s2_qq; *it += *ittt1 * (*ittt2);
            }
          }
        }
      }
      return 0;
    }
    ga_instruction_contraction_opt0_2_dunrolled
    (base_tensor &t_, base_tensor &tc1_, base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  // Performs Ani Bmi -> Cmn
  struct ga_instruction_contraction_opt2_0 : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type n, q;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: contraction operation of size " << n*q <<
                    " optimized for vectorized second tensor of type 2");
      size_type nn = n*q, s1 = tc1.size()/nn, s2 = tc2.size()/nn;
      size_type s1_q = s1/q, s1_qq = s1*q, s2_qq = s2*q;
      GA_DEBUG_ASSERT(t.size() == s1*s2, "Internal error");

      auto it = t.begin();
      for (size_type i = 0; i < s1_q; ++i)  {
        auto it1 = tc1.begin() + i*q;
        for (size_type l = 0; l < q; ++l) {
          auto it2 = tc2.begin() + l*s2;
          for (size_type j = 0; j < s2; ++j, ++it, ++it2) {
            auto itt1 = it1, itt2 = it2;
            *it = *itt1 * (*itt2);
            for (size_type m = 1; m < n; ++m) {
              itt1 += s1_qq, itt2 += s2_qq; *it += *itt1 * (*itt2);
            }
          }
        }
      }
      return 0;
    }
    ga_instruction_contraction_opt2_0(base_tensor &t_, base_tensor &tc1_,
                                    base_tensor &tc2_, size_type n_,
                                    size_type q_)
      : t(t_), tc1(tc1_), tc2(tc2_), n(n_), q(q_) { }
  };

  // Performs Ani Bmi -> Cmn
  template <int N>
  struct ga_instruction_contraction_opt2_0_unrolled : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type q;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: unrolled contraction of size " << N*q
                    << " optimized for vectorized second tensor of type 2");
      size_type nn = N*q, s1 = tc1.size()/nn, s2 = tc2.size()/nn;
      size_type s1_q = s1/q, s1_qq = s1*q, s2_qq = s2*q;
      GA_DEBUG_ASSERT(t.size() == s1*s2, "Internal error");

      auto it = t.begin(), it1 = tc1.begin();
      for (size_type i = 0; i < s1_q; ++i, it1 += q)  {
        for (size_type l = 0; l < q; ++l) {
          auto it2 = tc2.begin() + l*s2;
          for (size_type j = 0; j < s2; ++j, ++it, ++it2) {
            auto itt1 = it1, itt2 = it2;
            *it = *itt1 * (*itt2);
            for (size_type m = 1; m < N; ++m) {
              itt1 += s1_qq, itt2 += s2_qq; *it += *itt1 * (*itt2);
            }
          }
        }
      }
      return 0;
    }
    ga_instruction_contraction_opt2_0_unrolled(base_tensor &t_, base_tensor &tc1_,
                                             base_tensor &tc2_, size_type q_)
      : t(t_), tc1(tc1_), tc2(tc2_), q(q_) {}
  };

  // Performs Ani Bmi -> Cmn
  template <int N, int Q>
  struct ga_instruction_contraction_opt2_0_dunrolled : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: unrolled contraction of size " << N*Q
                    << " optimized for vectorized second tensor of type 2");
      size_type s1 = tc1.size()/(N*Q), s2 = tc2.size()/(N*Q);
      size_type s1_q = s1/Q, s1_qq = s1*Q, s2_qq = s2*Q;
      GA_DEBUG_ASSERT(t.size() == s1*s2, "Internal error");

      auto it = t.begin(), it1 = tc1.begin();
      for (size_type i = 0; i < s1_q; ++i, it1 += Q)  {
        for (size_type l = 0; l < Q; ++l) {
          auto it2 = tc2.begin() + l*s2;
          for (size_type j = 0; j < s2; ++j, ++it, ++it2) {
            auto itt1 = it1, itt2 = it2;
            *it = *itt1 * (*itt2);
            for (size_type m = 1; m < N; ++m) {
              itt1 += s1_qq, itt2 += s2_qq; *it += *itt1 * (*itt2);
            }
          }
        }
      }
      return 0;
    }
    ga_instruction_contraction_opt2_0_dunrolled
    (base_tensor &t_, base_tensor &tc1_, base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  // Performs Ani Bmi -> Cmn
  struct ga_instruction_contraction_opt0_1 : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: contraction operation of size " << nn <<
                    " optimized for vectorized second tensor of type 1");
      size_type ss1=tc1.size(), s1 = ss1/nn, s2=tc2.size()/nn, s2_n=s2/nn;

      auto it = t.begin(), it1 = tc1.begin();
      for (size_type i = 0; i < s1; ++i, ++it1) {
        auto it2 = tc2.begin();
        for (size_type j = 0; j < s2_n; ++j) {
          if (j) it2 += nn;
          auto itt1 = it1;
          *it++ = (*itt1) * (*it2);
          for (size_type k = 1; k < nn; ++k)
            { itt1 += s1; *it++ = (*itt1) * (*it2); }
        }
      }
      return 0;
    }
    ga_instruction_contraction_opt0_1(base_tensor &t_, base_tensor &tc1_,
                                    base_tensor &tc2_, size_type n_)
      : t(t_), tc1(tc1_), tc2(tc2_), nn(n_) {}
  };

  template<int N> inline void reduc_elem_unrolled_opt1_
  (const base_vector::iterator &it, const base_vector::iterator &it1,
   scalar_type a, size_type s1) {
    it[N-1] = it1[(N-1)*s1] * a;
    reduc_elem_unrolled_opt1_<N-1>(it, it1, a, s1);
  }
  template<> inline void reduc_elem_unrolled_opt1_<1>
  (const base_vector::iterator &it, const base_vector::iterator &it1,
   scalar_type a, size_type /* s1 */)
  { *it = (*it1) * a; }

  // Performs Ani Bmi -> Cmn
  template <int N>
  struct ga_instruction_contraction_opt0_1_unrolled : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: unrolled contraction operation of size " << N
                    << " optimized for vectorized second tensor of type 1");
      size_type s1 = tc1.size()/N, s2 = tc2.size()/N;
      auto it = t.begin(), it1 = tc1.begin();
      for (size_type i = 0; i < s1; ++i, ++it1) {
        auto it2 = tc2.begin(), it2e = it2 + s2;
        for (; it2 != it2e; it2 += N, it += N)
          reduc_elem_unrolled_opt1_<N>(it, it1, *it2, s1);
      }
      return 0;
    }
    ga_instruction_contraction_opt0_1_unrolled(base_tensor &t_, base_tensor &tc1_,
                                             base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  // Performs Ani Bmi -> Cmn
  struct ga_instruction_contraction_opt1_1 : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: contraction operation of size " << nn <<
                    " optimized for both vectorized tensor of type 1");
      size_type s1 = tc1.size()/nn, s2 = tc2.size()/nn, s2_1 = s2+1;
      GA_DEBUG_ASSERT(t.size() == s2*s1, "Internal error");
      size_type ss1 = s1/nn, ss2 = s2/nn;

      // std::fill(t.begin(), t.end(), scalar_type(0)); // Factorized
      auto it2 = tc2.begin();
      for (size_type j = 0; j < ss2; ++j) {
        if (j) it2 += nn;
        auto it1 = tc1.begin(), it = t.begin() + j*nn;
        for (size_type i = 0; i < ss1; ++i) {
          if (i) { it1 += nn, it += s2*nn; }
          scalar_type a = (*it1) * (*it2);
          auto itt = it;
          *itt = a; itt += s2_1; *itt = a;
          for (size_type k = 2; k < nn; ++k) { itt += s2_1; *itt = a; }
        }
      }
      return 0;
    }
    ga_instruction_contraction_opt1_1(base_tensor &t_, base_tensor &tc1_,
                                    base_tensor &tc2_, size_type n_)
      : t(t_), tc1(tc1_), tc2(tc2_), nn(n_) {}
  };



  template<int N> inline scalar_type reduc_elem_unrolled__
  (base_tensor::iterator &it1, base_tensor::iterator &it2,
   size_type s1, size_type s2) {
    return (it1[(N-1)*s1])*(it2[(N-1)*s2])
      + reduc_elem_unrolled__<N-1>(it1, it2, s1, s2);
  }
  template<> inline scalar_type reduc_elem_unrolled__<1>
  (base_tensor::iterator &it1, base_tensor::iterator &it2,
   size_type /*s1*/, size_type /*s2*/)
  { return (*it1)*(*it2); }

  // Performs Ani Bmi -> Cmn. Unrolled operation.
  template<int N> struct ga_instruction_contraction_unrolled
    : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: unrolled contraction operation of size " << N);
      size_type s1 = tc1.size()/N, s2 = tc2.size()/N;
      GA_DEBUG_ASSERT(t.size() == s1*s2, "Internal error, " << t.size()
                      << " != " << s1 << "*" << s2);
      base_tensor::iterator it1=tc1.begin(), it2=tc2.begin(), it2end=it2 + s2;
      for (base_tensor::iterator it = t.begin(); it != t.end(); ++it) {
        *it = reduc_elem_unrolled__<N>(it1, it2, s1, s2);
        ++it2; if (it2 == it2end) { it2 = tc2.begin(), ++it1; }
      }
      return 0;
    }
    ga_instruction_contraction_unrolled(base_tensor &t_, base_tensor &tc1_,
                                      base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  template<int N, int S2> inline void reduc_elem_d_unrolled__
  (base_tensor::iterator &it, base_tensor::iterator &it1,
   base_tensor::iterator &it2, size_type s1, size_type s2)  {
    *it++ = reduc_elem_unrolled__<N>(it1, it2, s1, s2);
    reduc_elem_d_unrolled__<N, S2-1>(it, it1, ++it2, s1, s2);
  }
  // A Repeated definition is following because partial specialization
  // of functions is not allowed in C++ for the moment.
  // The gain in assembly time is small compared to the simply unrolled version
  template<> inline void reduc_elem_d_unrolled__<1, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<2, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<3, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<4, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<5, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<6, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<7, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<8, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<9, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<10, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<11, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<12, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<13, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<14, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<15, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }
  template<> inline void reduc_elem_d_unrolled__<16, 0>
  (base_tensor::iterator &/* it */, base_tensor::iterator &/* it1 */,
   base_tensor::iterator &/* it2 */, size_type /* s1 */, size_type /* s2 */) { }

  // Performs Ani Bmi -> Cmn. Automatically doubly unrolled operation
  // (for uniform meshes).
  template<int N, int S2> struct ga_ins_red_d_unrolled
    : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: doubly unrolled contraction operation of size "
                    << S2 << "x" << N);
      size_type s1 = tc1.size()/N, s2 = tc2.size()/N;
      GA_DEBUG_ASSERT(s2 == S2, "Internal error");
      GA_DEBUG_ASSERT(t.size() == s1*s2, "Internal error, " << t.size()
                      << " != " << s1 << "*" << s2);
      base_tensor::iterator it = t.begin(), it1 = tc1.begin();
      for (size_type ii = 0; ii < s1; ++ii, ++it1) {
        base_tensor::iterator it2 = tc2.begin();
        reduc_elem_d_unrolled__<N, S2>(it, it1, it2, s1, s2);
      }
      GA_DEBUG_ASSERT(it == t.end(), "Internal error");
      return 0;
    }
    ga_ins_red_d_unrolled(base_tensor &t_, base_tensor &tc1_, base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };


  pga_instruction ga_instruction_contraction_switch
  (assembly_tensor &t_, assembly_tensor &tc1_, assembly_tensor &tc2_,
   size_type n, bool &to_clear) {
    base_tensor &t = t_.tensor(), &tc1 = tc1_.tensor(), &tc2 = tc2_.tensor();

    if (tc1_.sparsity() == 1 && tc2_.sparsity() == 1 &&
        tc1_.qdim() == n && tc2_.qdim() == n) {
      to_clear = true;
      t_.set_sparsity(10, tc1_.qdim());
      return std::make_shared<ga_instruction_contraction_opt1_1>(t, tc1, tc2, n);
    }

    if (tc2_.sparsity() == 1) {
      switch(n) {
      case 2:
        return std::make_shared<ga_instruction_contraction_opt0_1_unrolled<2>>
          (t, tc1, tc2);
      case 3:
        return std::make_shared<ga_instruction_contraction_opt0_1_unrolled<3>>
          (t, tc1, tc2);
      case 4:
        return std::make_shared<ga_instruction_contraction_opt0_1_unrolled<4>>
          (t, tc1, tc2);
      case 5:
        return std::make_shared<ga_instruction_contraction_opt0_1_unrolled<5>>
          (t, tc1, tc2);
      default:
        return std::make_shared<ga_instruction_contraction_opt0_1>(t,tc1,tc2, n);
      }
    }
    if (tc2_.sparsity() == 2) {
      size_type q2 = tc2.sizes()[1];
      size_type n2 = (tc2.sizes().size() > 2) ? tc2.sizes()[2] : 1;
      if (n2*q2 == n) {
        switch (n2) {
        case 1:
          switch (q2) {
          case 2:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<1,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<1,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<1,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_contraction_opt0_2_unrolled<1>>
              (t, tc1, tc2, q2);
          }
        case 2:
          switch (q2) {
          case 2:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<2,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<2,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<2,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_contraction_opt0_2_unrolled<2>>
              (t, tc1, tc2, q2);
          }
        case 3:
          switch (q2) {
          case 2:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<3,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<3,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<3,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_contraction_opt0_2_unrolled<3>>
              (t, tc1, tc2, q2);
          }
        case 4:
          return std::make_shared<ga_instruction_contraction_opt0_2_unrolled<4>>
            (t, tc1, tc2, q2);
        case 5:
          return std::make_shared<ga_instruction_contraction_opt0_2_unrolled<5>>
            (t, tc1, tc2, q2);
        default:
          return std::make_shared<ga_instruction_contraction_opt0_2>
            (t,tc1,tc2,n2,q2);
        }
      }
    }
    if (tc1_.sparsity() == 2) {
      size_type q1 = tc1.sizes()[1];
      size_type n1 = (tc1.sizes().size() > 2) ? tc1.sizes()[2] : 1;
      if (n1*q1 == n) {
        switch (n1) {
        case 1:
          switch (q1) {
          case 2:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<1,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<1,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<1,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_contraction_opt2_0_unrolled<1>>
              (t, tc1, tc2, q1);
          }
        case 2:
          switch (q1) {
          case 2:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<2,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<2,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<2,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_contraction_opt2_0_unrolled<2>>
              (t, tc1, tc2, q1);
          }
        case 3:
          switch (q1) {
          case 2:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<3,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<3,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<3,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_contraction_opt2_0_unrolled<3>>
              (t, tc1, tc2, q1);
          }
          return std::make_shared<ga_instruction_contraction_opt2_0_unrolled<3>>
            (t, tc1, tc2, q1);
        case 4:
          return std::make_shared<ga_instruction_contraction_opt2_0_unrolled<4>>
            (t, tc1, tc2, q1);
        case 5:
          return std::make_shared<ga_instruction_contraction_opt2_0_unrolled<5>>
            (t, tc1, tc2, q1);
        default:
          return std::make_shared<ga_instruction_contraction_opt2_0>
            (t,tc1,tc2, n1, q1);
        }
      }
    }

    switch(n) {
    case  2 : return std::make_shared<ga_instruction_contraction_unrolled< 2>>
                     (t, tc1, tc2);
    case  3 : return std::make_shared<ga_instruction_contraction_unrolled< 3>>
                     (t, tc1, tc2);
    case  4 : return std::make_shared<ga_instruction_contraction_unrolled< 4>>
                     (t, tc1, tc2);
    case  5 : return std::make_shared<ga_instruction_contraction_unrolled< 5>>
                     (t, tc1, tc2);
    case  6 : return std::make_shared<ga_instruction_contraction_unrolled< 6>>
                     (t, tc1, tc2);
    case  7 : return std::make_shared<ga_instruction_contraction_unrolled< 7>>
                     (t, tc1, tc2);
    case  8 : return std::make_shared<ga_instruction_contraction_unrolled< 8>>
                     (t, tc1, tc2);
    case  9 : return std::make_shared<ga_instruction_contraction_unrolled< 9>>
                     (t, tc1, tc2);
    case 10 : return std::make_shared<ga_instruction_contraction_unrolled<10>>
                     (t, tc1, tc2);
    case 11 : return std::make_shared<ga_instruction_contraction_unrolled<11>>
                     (t, tc1, tc2);
    case 12 : return std::make_shared<ga_instruction_contraction_unrolled<12>>
                     (t, tc1, tc2);
    case 13 : return std::make_shared<ga_instruction_contraction_unrolled<13>>
                     (t, tc1, tc2);
    case 14 : return std::make_shared<ga_instruction_contraction_unrolled<14>>
                     (t, tc1, tc2);
    case 15 : return std::make_shared<ga_instruction_contraction_unrolled<15>>
                     (t, tc1, tc2);
    case 16 : return std::make_shared<ga_instruction_contraction_unrolled<16>>
                     (t, tc1, tc2);
    default : return std::make_shared<ga_instruction_contraction>
                     (t, tc1, tc2, n);
    }
  }

  pga_instruction  ga_uniform_instruction_contraction_switch
  (assembly_tensor &t_, assembly_tensor &tc1_, assembly_tensor &tc2_,
   size_type n, bool &to_clear) {
    base_tensor &t = t_.tensor(), &tc1 = tc1_.tensor(), &tc2 = tc2_.tensor();

    if (tc1_.sparsity() == 1 && tc2_.sparsity() == 1 &&
        tc1_.qdim() == n && tc2_.qdim() == n) {
      to_clear = true;
      t_.set_sparsity(10, tc1_.qdim());
      return std::make_shared<ga_instruction_contraction_opt1_1>(t,tc1,tc2,n);
    }
    if (tc2_.sparsity() == 1) {
      switch(n) {
      case 2:
        return std::make_shared<ga_instruction_contraction_opt0_1_unrolled<2>>
          (t, tc1, tc2);
      case 3:
        return std::make_shared<ga_instruction_contraction_opt0_1_unrolled<3>>
          (t, tc1, tc2);
      case 4:
        return std::make_shared<ga_instruction_contraction_opt0_1_unrolled<4>>
          (t, tc1, tc2);
      case 5:
        return std::make_shared<ga_instruction_contraction_opt0_1_unrolled<5>>
          (t, tc1, tc2);
      default:
        return std::make_shared<ga_instruction_contraction_opt0_1>(t,tc1,tc2, n);
      }
    }
    if (tc2_.sparsity() == 2) {
      size_type q2 = tc2.sizes()[1];
      size_type n2 = (tc2.sizes().size() > 2) ? tc2.sizes()[2] : 1;
      if (n2*q2 == n) {
        switch (n2) {
        case 1:
          switch (q2) {
          case 2:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<1,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<1,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<1,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_contraction_opt0_2_unrolled<1>>
              (t, tc1, tc2, q2);
          }
        case 2:
          switch (q2) {
          case 2:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<2,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<2,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<2,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_contraction_opt0_2_unrolled<2>>
              (t, tc1, tc2, q2);
          }
        case 3:
          switch (q2) {
          case 2:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<3,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<3,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_contraction_opt0_2_dunrolled<3,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_contraction_opt0_2_unrolled<3>>
              (t, tc1, tc2, q2);
          }
        case 4:
          return std::make_shared<ga_instruction_contraction_opt0_2_unrolled<4>>
            (t, tc1, tc2, q2);
        case 5:
          return std::make_shared<ga_instruction_contraction_opt0_2_unrolled<5>>
            (t, tc1, tc2, q2);
        default:
          return std::make_shared<ga_instruction_contraction_opt0_2>
            (t,tc1,tc2,n2,q2);
        }
      }
    }
    if (tc1_.sparsity() == 2) {
      size_type q1 = tc1.sizes()[1];
      size_type n1 = (tc1.sizes().size() > 2) ? tc1.sizes()[2] : 1;
      if (n1*q1 == n) {
        switch (n1) {
        case 1:
          switch (q1) {
          case 2:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<1,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<1,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<1,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_contraction_opt2_0_unrolled<1>>
              (t, tc1, tc2, q1);
          }
        case 2:
          switch (q1) {
          case 2:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<2,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<2,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<2,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_contraction_opt2_0_unrolled<2>>
              (t, tc1, tc2, q1);
          }
        case 3:
          switch (q1) {
          case 2:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<3,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<3,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_contraction_opt2_0_dunrolled<3,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_contraction_opt2_0_unrolled<3>>
              (t, tc1, tc2, q1);
          }
          return std::make_shared<ga_instruction_contraction_opt2_0_unrolled<3>>
            (t, tc1, tc2, q1);
        case 4:
          return std::make_shared<ga_instruction_contraction_opt2_0_unrolled<4>>
            (t, tc1, tc2, q1);
        case 5:
          return std::make_shared<ga_instruction_contraction_opt2_0_unrolled<5>>
            (t, tc1, tc2, q1);
        default:
          return std::make_shared<ga_instruction_contraction_opt2_0>
            (t,tc1,tc2, n1, q1);
        }
      }
    }

    // Only specialized for certain values
    size_type s2 = tc2.size()/n;
    switch(s2) {
    case 1 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,1>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,1>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,1>>(t, tc1, tc2);
      default: return ga_instruction_contraction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 2 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,2>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,2>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,2>>(t, tc1, tc2);
      default: return ga_instruction_contraction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 3 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,3>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,3>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,3>>(t, tc1, tc2);
      default: return ga_instruction_contraction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 4 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,4>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,4>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,4>>(t, tc1, tc2);
      default: return ga_instruction_contraction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 5 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,5>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,5>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,5>>(t, tc1, tc2);
      default: return ga_instruction_contraction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 6 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,6>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,6>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,6>>(t, tc1, tc2);
      default: return ga_instruction_contraction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 7 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,7>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,7>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,7>>(t, tc1, tc2);
      default: return ga_instruction_contraction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 8 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,8>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,8>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,8>>(t, tc1, tc2);
      default: return ga_instruction_contraction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 9 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,9>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,9>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,9>>(t, tc1, tc2);
      default: return ga_instruction_contraction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 10:
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,10>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,10>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,10>>(t, tc1, tc2);
      default: return ga_instruction_contraction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    default: return ga_instruction_contraction_switch(t_,tc1_,tc2_,n,to_clear);
    }
  }


  // Performs Amij Bnj -> Cmni. To be optimized.
  struct ga_instruction_spec_contraction : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: specific contraction operation of "
                    "size " << nn);
      size_type s1 = tc1.sizes()[0], s11 = tc1.size() / (s1*nn), s111 = s1*s11;
      size_type s2 = tc2.sizes()[0];
      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < s11; ++i)
        for (size_type n = 0; n < s2; ++n)
          for (size_type m = 0; m < s1; ++m, ++it) {
            *it = scalar_type(0);
            for (size_type j = 0; j < nn; ++j)
              *it += tc1[m+i*s1+j*s111] * tc2[n+j*s2];
          }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_spec_contraction(base_tensor &t_, base_tensor &tc1_,
                                  base_tensor &tc2_, size_type n_)
      : t(t_), tc1(tc1_), tc2(tc2_), nn(n_) {}
  };

  // Performs Amik Bnjk -> Cmnij. To be optimized.
  struct ga_instruction_spec2_contraction : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: second specific contraction operation of "
                    "size " << nn);
      size_type s1 = tc1.sizes()[0], s11 = tc1.size() / (s1*nn), s111 = s1*s11;
      size_type s2 = tc2.sizes()[0], s22 = tc2.size() / (s2*nn), s222 = s2*s22;
      base_tensor::iterator it = t.begin();
      for (size_type j = 0; j < s22; ++j)
        for (size_type i = 0; i < s11; ++i)
          for (size_type m = 0; m < s1; ++m)
            for (size_type n = 0; n < s2; ++n, ++it) {
              *it = scalar_type(0);
              for (size_type k = 0; k < nn; ++k)
                *it += tc1[m+i*s1+k*s111] * tc2[n+j*s2+k*s222];
            }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_spec2_contraction(base_tensor &t_, base_tensor &tc1_,
                                   base_tensor &tc2_, size_type n_)
      : t(t_), tc1(tc1_), tc2(tc2_), nn(n_) {}
  };

  // Performs Aij Bkl -> Cijkl
  struct ga_instruction_simple_tmult : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: simple tensor product");
      size_type s1 = tc1.size();
      GA_DEBUG_ASSERT(t.size() == s1 * tc2.size(), "Wrong sizes");
      base_tensor::iterator it2=tc2.begin(), it1=tc1.begin(), it1end=it1 + s1;
      for (base_tensor::iterator it = t.begin(); it != t.end(); ++it) {
        *it = *(it2) * (*it1);
        ++it1; if (it1 == it1end) { it1 = tc1.begin(), ++it2; }
      }
      return 0;
    }
    ga_instruction_simple_tmult(base_tensor &t_, base_tensor &tc1_,
                                base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  template<int S1> inline void tmult_elem_unrolled__
  (base_tensor::iterator &it, base_tensor::iterator &it1,
   base_tensor::iterator &it2) {
    *it++ = (*it1++)*(*it2);
    tmult_elem_unrolled__<S1-1>(it, it1, it2);
  }
  template<> inline void tmult_elem_unrolled__<0>
  (base_tensor::iterator &/*it*/, base_tensor::iterator &/*it1*/,
   base_tensor::iterator &/*it2*/) { }

  // Performs Aij Bkl -> Cijkl, partially unrolled version
  template<int S1> struct ga_instruction_simple_tmult_unrolled
  : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      size_type s2 = tc2.size();
      GA_DEBUG_ASSERT(tc1.size() == S1,
                      "Wrong sizes " << tc1.size() << " != " << S1);
      GA_DEBUG_INFO("Instruction: simple tensor product, unrolled with "
                    << S1 << " operations");
      GA_DEBUG_ASSERT(t.size() == S1 * s2,
                      "Wrong sizes " << t.size() << " != " << S1 << "*" << s2);
      base_tensor::iterator it = t.begin(), it2 = tc2.begin();
      for (size_type ii = 0; ii < s2; ++ii, ++it2) {
        base_tensor::iterator it1 = tc1.begin();
        tmult_elem_unrolled__<S1>(it, it1, it2);
      }
      GA_DEBUG_ASSERT(it == t.end(), "Internal error");
      return 0;
    }
    ga_instruction_simple_tmult_unrolled(base_tensor &t_, base_tensor &tc1_,
                                         base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  pga_instruction  ga_uniform_instruction_simple_tmult
  (base_tensor &t, base_tensor &tc1, base_tensor &tc2) {
    switch(tc1.size()) {
    case  2 : return std::make_shared<ga_instruction_simple_tmult_unrolled< 2>>
                     (t, tc1, tc2);
    case  3 : return std::make_shared<ga_instruction_simple_tmult_unrolled< 3>>
                     (t, tc1, tc2);
    case  4 : return std::make_shared<ga_instruction_simple_tmult_unrolled< 4>>
                     (t, tc1, tc2);
    case  5 : return std::make_shared<ga_instruction_simple_tmult_unrolled< 5>>
                     (t, tc1, tc2);
    case  6 : return std::make_shared<ga_instruction_simple_tmult_unrolled< 6>>
                     (t, tc1, tc2);
    case  7 : return std::make_shared<ga_instruction_simple_tmult_unrolled< 7>>
                     (t, tc1, tc2);
    case  8 : return std::make_shared<ga_instruction_simple_tmult_unrolled< 8>>
                     (t, tc1, tc2);
    case  9 : return std::make_shared<ga_instruction_simple_tmult_unrolled< 9>>
                     (t, tc1, tc2);
    case 10 : return std::make_shared<ga_instruction_simple_tmult_unrolled<10>>
                     (t, tc1, tc2);
    case 11 : return std::make_shared<ga_instruction_simple_tmult_unrolled<11>>
                     (t, tc1, tc2);
    case 12 : return std::make_shared<ga_instruction_simple_tmult_unrolled<12>>
                     (t, tc1, tc2);
    case 13 : return std::make_shared<ga_instruction_simple_tmult_unrolled<13>>
                     (t, tc1, tc2);
    case 14 : return std::make_shared<ga_instruction_simple_tmult_unrolled<14>>
                     (t, tc1, tc2);
    case 15 : return std::make_shared<ga_instruction_simple_tmult_unrolled<15>>
                     (t, tc1, tc2);
    case 16 : return std::make_shared<ga_instruction_simple_tmult_unrolled<16>>
                     (t, tc1, tc2);
    default : return std::make_shared<ga_instruction_simple_tmult>
                     (t, tc1, tc2);
    }
  }


  // Performs Ami Bnj -> Cmnij. To be optimized.
  struct ga_instruction_spec_tmult : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type s1_2, s2_2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: specific tensor product");
      GA_DEBUG_ASSERT(t.size() == tc1.size() * tc2.size(), "Wrong sizes");
      size_type s1_1 = tc1.size() / s1_2;
      size_type s2_1 = tc2.size() / s2_2;

      base_tensor::iterator it = t.begin();
      for (size_type j = 0; j < s2_2; ++j)
        for (size_type i = 0; i < s1_2; ++i)
          for (size_type n = 0; n < s2_1; ++n)
            for (size_type m = 0; m < s1_1; ++m, ++it)
              *it = tc1[m+i*s1_1] * tc2[n+j*s2_1];
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_spec_tmult(base_tensor &t_, base_tensor &tc1_,
                              base_tensor &tc2_, size_type s1_2_,
                              size_type s2_2_)
      : t(t_), tc1(tc1_), tc2(tc2_), s1_2(s1_2_), s2_2(s2_2_) {}
  };

  // Performs Ai Bmj -> Cmij. To be optimized.
  struct ga_instruction_spec2_tmult : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: second specific tensor product");
      GA_DEBUG_ASSERT(t.size() == tc1.size() * tc2.size(), "Wrong sizes");
      size_type s1 = tc1.size();
      size_type s2_1 = tc2.sizes()[0], s2_2 = tc2.size() / s2_1;

      base_tensor::iterator it = t.begin();
      for (size_type j = 0; j < s2_2; ++j)
        for (size_type i = 0; i < s1; ++i)
          for (size_type m = 0; m < s2_1; ++m, ++it)
            *it = tc1[i] * tc2[m+j*s2_1];
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_spec2_tmult(base_tensor &t_, base_tensor &tc1_,
                              base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };



  struct ga_instruction_simple_c_matrix : public ga_instruction {
    base_tensor &t;
    std::vector<scalar_type *> components;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: gathering components for explicit "
                    "matrix");
      GA_DEBUG_ASSERT(t.size() == components.size(), "Wrong sizes");
      for (size_type i = 0; i < components.size(); ++i)
        t[i] = *(components[i]);
      return 0;
    }
    ga_instruction_simple_c_matrix(base_tensor &t_,
                                   std::vector<scalar_type *> &components_)
      : t(t_), components(components_) {}
  };

  struct ga_instruction_c_matrix_with_tests : public ga_instruction {
    base_tensor &t;
    const std::vector<const base_tensor *> components;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: gathering components for explicit "
                    "matrix with tests functions");
      size_type s = t.size() / components.size();
      GA_DEBUG_ASSERT(s, "Wrong sizes");
      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < components.size(); ++i) {
        const base_tensor &t1 = *(components[i]);
        if (t1.size() > 1) {
          GA_DEBUG_ASSERT(t1.size() == s, "Wrong sizes, " << t1.size()
                          << " != " << s);
          for (size_type j = 0; j < s; ++j) *it++ = t1[j];
        } else {
          for (size_type j = 0; j < s; ++j) *it++ = t1[0];
        }
      }
      return 0;
    }
    ga_instruction_c_matrix_with_tests
    (base_tensor &t_, const std::vector<const base_tensor *> &components_)
      : t(t_), components(components_) {}
  };

  struct ga_instruction_eval_func_1arg_1res : public ga_instruction {
    scalar_type &t;
    const scalar_type &c;
    pscalar_func_onearg f1;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: evaluation of a one argument "
                    "predefined function on a scalar");
      t = (*f1)(c);
      return 0;
    }
    ga_instruction_eval_func_1arg_1res(scalar_type &t_, const scalar_type &c_,
                                       pscalar_func_onearg f1_)
      : t(t_), c(c_), f1(f1_) {}
  };

  struct ga_instruction_eval_func_1arg_1res_expr : public ga_instruction {
    scalar_type &t;
    const scalar_type &c;
    const ga_predef_function &F;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: evaluation of a one argument "
                    "predefined function on a scalar");
      t = F(c);
      return 0;
    }
    ga_instruction_eval_func_1arg_1res_expr(scalar_type &t_,
                                            const scalar_type &c_,
                                            const ga_predef_function &F_)
      : t(t_), c(c_), F(F_) {}
  };

  struct ga_instruction_eval_func_1arg : public ga_instruction {
    base_tensor &t, &tc1;
    pscalar_func_onearg f1;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: evaluation of a one argument "
                    "predefined function on tensor");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      for (size_type i = 0; i < t.size(); ++i)  t[i] = (*f1)(tc1[i]);
      return 0;
    }
    ga_instruction_eval_func_1arg(base_tensor &t_, base_tensor &c_,
                                  pscalar_func_onearg f1_)
      : t(t_), tc1(c_), f1(f1_) {}
  };

  struct ga_instruction_eval_func_1arg_expr : public ga_instruction {
    base_tensor &t, &tc1;
    const ga_predef_function &F;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: evaluation of a one argument "
                    "predefined function on tensor");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      for (size_type i = 0; i < t.size(); ++i)  t[i] = F(tc1[i]);
      return 0;
    }
    ga_instruction_eval_func_1arg_expr(base_tensor &t_, base_tensor &c_,
                                       const ga_predef_function &F_)
      : t(t_), tc1(c_), F(F_) {}
  };

  struct ga_instruction_eval_func_2arg_1res : public ga_instruction {
    scalar_type &t;
    const scalar_type &c, &d;
    pscalar_func_twoargs f2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                    "predefined function on two scalar");
      t = (*f2)(c, d);
      return 0;
    }
    ga_instruction_eval_func_2arg_1res(scalar_type &t_, const scalar_type &c_,
                                       const scalar_type &d_,
                                       pscalar_func_twoargs f2_)
      : t(t_), c(c_), d(d_), f2(f2_) {}
  };

  struct ga_instruction_eval_func_2arg_1res_expr : public ga_instruction {
    scalar_type &t;
    const scalar_type &c, &d;
    const ga_predef_function &F;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                           "predefined function on two scalar");
      t = F(c, d);
      return 0;
    }
    ga_instruction_eval_func_2arg_1res_expr(scalar_type &t_,
                                            const scalar_type &c_,
                                            const scalar_type &d_,
                                            const ga_predef_function &F_)
      : t(t_), c(c_), d(d_), F(F_) {}
  };

  struct ga_instruction_eval_func_2arg_first_scalar : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    pscalar_func_twoargs f2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                    "predefined function on one scalar and one tensor");
      GA_DEBUG_ASSERT(t.size() == tc2.size(), "Wrong sizes");
      for (size_type i = 0; i < t.size(); ++i)  t[i] = (*f2)(tc1[0], tc2[i]);
      return 0;
    }
    ga_instruction_eval_func_2arg_first_scalar
    (base_tensor &t_, base_tensor &c_, base_tensor &d_,
     pscalar_func_twoargs f2_)
      : t(t_), tc1(c_), tc2(d_), f2(f2_) {}
  };

  struct ga_instruction_eval_func_2arg_first_scalar_expr
    : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    const ga_predef_function &F;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                    "predefined function on one scalar and one tensor");
      GA_DEBUG_ASSERT(t.size() == tc2.size(), "Wrong sizes");
      for (size_type i = 0; i < t.size(); ++i)  t[i] = F(tc1[0], tc2[i]);
      return 0;
    }
    ga_instruction_eval_func_2arg_first_scalar_expr
    (base_tensor &t_, base_tensor &c_, base_tensor &d_,
     const ga_predef_function &F_)
      : t(t_), tc1(c_), tc2(d_), F(F_) {}
  };

  struct ga_instruction_eval_func_2arg_second_scalar : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    pscalar_func_twoargs f2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                    "predefined function on one tensor and one scalar");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      for (size_type i = 0; i < t.size(); ++i)  t[i] = (*f2)(tc1[i], tc2[0]);
      return 0;
    }
    ga_instruction_eval_func_2arg_second_scalar(base_tensor &t_,
                                                base_tensor &c_,
                                                base_tensor &d_,
                                                pscalar_func_twoargs f2_)
      : t(t_), tc1(c_), tc2(d_), f2(f2_) {}
  };

  struct ga_instruction_eval_func_2arg_second_scalar_expr
    : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    const ga_predef_function &F;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                    "predefined function on one tensor and one scalar");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      for (size_type i = 0; i < t.size(); ++i)  t[i] = F(tc1[i], tc2[0]);
      return 0;
    }
    ga_instruction_eval_func_2arg_second_scalar_expr
    (base_tensor &t_, base_tensor &c_, base_tensor &d_,
     const ga_predef_function &F_)
      : t(t_), tc1(c_), tc2(d_), F(F_) {}
  };

  struct ga_instruction_eval_func_2arg : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    pscalar_func_twoargs f2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                    "predefined function on two tensors");
      GA_DEBUG_ASSERT(t.size() == tc1.size() && t.size() == tc2.size(),
                      "Wrong sizes");

      for (size_type i = 0; i < t.size(); ++i)  t[i] = (*f2)(tc1[i], tc2[i]);
      return 0;
    }
    ga_instruction_eval_func_2arg(base_tensor &t_, base_tensor &c_,
                                  base_tensor &d_, pscalar_func_twoargs f2_)
      : t(t_), tc1(c_), tc2(d_), f2(f2_) {}
  };

  struct ga_instruction_eval_func_2arg_expr : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    const ga_predef_function &F;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                    "predefined function on two tensors");
      GA_DEBUG_ASSERT(t.size() == tc1.size() && t.size() == tc2.size(),
                      "Wrong sizes");

      for (size_type i = 0; i < t.size(); ++i)  t[i] = F(tc1[i], tc2[i]);
      return 0;
    }
    ga_instruction_eval_func_2arg_expr(base_tensor &t_, base_tensor &c_,
                                       base_tensor &d_,
                                       const ga_predef_function &F_)
      : t(t_), tc1(c_), tc2(d_), F(F_) {}
  };

  struct ga_instruction_eval_OP : public ga_instruction {
    base_tensor &t;
    const ga_nonlinear_operator &OP;
    ga_nonlinear_operator::arg_list args;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: operator evaluation");
      OP.value(args, t);
      return 0;
    }
    ga_instruction_eval_OP(base_tensor &t_, const ga_nonlinear_operator &OP_,
                           ga_nonlinear_operator::arg_list &args_)
      : t(t_), OP(OP_), args(args_) {}
  };

  struct ga_instruction_eval_derivative_OP : public ga_instruction {
    base_tensor &t;
    const ga_nonlinear_operator &OP;
    ga_nonlinear_operator::arg_list args;
    size_type der1;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: operator derivative evaluation");
      OP.derivative(args, der1, t);
      return 0;
    }
    ga_instruction_eval_derivative_OP(base_tensor &t_,
                                      const ga_nonlinear_operator &OP_,
                                      ga_nonlinear_operator::arg_list &args_,
                                      size_type der1_)
      : t(t_), OP(OP_), args(args_), der1(der1_) {}
  };

  struct ga_instruction_eval_second_derivative_OP : public ga_instruction {
    base_tensor &t;
    const ga_nonlinear_operator &OP;
    ga_nonlinear_operator::arg_list args;
    size_type der1, der2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: operator second derivative evaluation");
      OP.second_derivative(args, der1, der2, t);
      return 0;
    }
    ga_instruction_eval_second_derivative_OP
    (base_tensor &t_, const ga_nonlinear_operator &OP_,
     ga_nonlinear_operator::arg_list &args_, size_type der1_, size_type der2_)
      : t(t_), OP(OP_), args(args_), der1(der1_), der2(der2_) {}
  };

  struct ga_instruction_tensor_slice : public ga_instruction {
    base_tensor &t, &tc1;
    bgeot::multi_index mi, indices;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: tensor slice");
      size_type order = t.sizes().size();
      for (bgeot::multi_index mi3(order); !mi3.finished(t.sizes());
           mi3.incrementation(t.sizes())) {
        for (size_type j = 0; j < order; ++j)
          mi[indices[j]] = mi3[j];
        t(mi3) = tc1(mi);
      }
      return 0;
    }
    ga_instruction_tensor_slice(base_tensor &t_, base_tensor &tc1_,
                                bgeot::multi_index &mi_,
                                bgeot::multi_index &indices_)
      : t(t_), tc1(tc1_), mi(mi_), indices(indices_)  {}
  };

  struct ga_instruction_transformation_call : public ga_instruction {
    const ga_workspace &workspace;
    ga_instruction_set::interpolate_info &inin;
    pinterpolate_transformation trans;
    fem_interpolation_context &ctx;
    const base_small_vector &Normal;
    const mesh &m;
    bool compute_der;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: call interpolate transformation");
      base_node P_ref;
      size_type cv;
      short_type face_num;
      inin.pt_type = trans->transform(workspace, m, ctx, Normal, &(inin.m), cv,
                                      face_num, P_ref, inin.Normal,
                                      inin.derivatives, compute_der);
      if (inin.pt_type) {
        if (cv != size_type(-1)) {
          inin.m->points_of_convex(cv, inin.G);
          inin.ctx.change((inin.m)->trans_of_convex(cv),
                          0, P_ref, inin.G, cv, face_num);
          inin.has_ctx = true;
          if (face_num != short_type(-1)) {
            inin.Normal = bgeot::compute_normal(inin.ctx, face_num);
            gmm::scale(inin.Normal, 1.0/gmm::vect_norm2(inin.Normal));
          } else
            inin.Normal.resize(0);
          inin.pt_y = inin.ctx.xreal();
        } else {
          inin.ctx.invalid_convex_num();
          inin.Normal.resize(0);
          inin.pt_y = P_ref;
          inin.has_ctx = false;
        }
      } else {
        inin.ctx.invalid_convex_num();
        inin.Normal.resize(0);
        inin.pt_y.resize(0);
        inin.has_ctx = false;
      }
      GA_DEBUG_INFO("Instruction: end of call interpolate transformation");
      return 0;
    }
    ga_instruction_transformation_call
    (const ga_workspace &w, ga_instruction_set::interpolate_info &i,
     pinterpolate_transformation t, fem_interpolation_context &ctxx,
     const base_small_vector &No, const mesh &mm, bool compute_der_)
      : workspace(w), inin(i), trans(t), ctx(ctxx), Normal(No), m(mm),
        compute_der(compute_der_) {}
  };

  struct ga_instruction_neighbor_transformation_call : public ga_instruction {
    const ga_workspace &workspace;
    ga_instruction_set::interpolate_info &inin;
    pinterpolate_transformation trans;
    fem_interpolation_context &ctx;
    base_small_vector dummy_normal;
    const mesh &m;
    size_type &ipt;
    papprox_integration &pai;
    bgeot::geotrans_precomp_pool &gp_pool;
    std::map<gauss_pt_corresp, bgeot::pstored_point_tab> &neighbor_corresp;

    virtual int exec() {
      bool cancel_optimization = false;
      GA_DEBUG_INFO("Instruction: call interpolate neighbor transformation");
      if (ipt == 0) {
        if (!(ctx.have_pgp()) || !pai || pai->is_built_on_the_fly()
            || cancel_optimization) {
          inin.ctx.invalid_convex_num();
        } else {
          // Test if the situation has already been encountered
          size_type cv = ctx.convex_num();
          short_type f = ctx.face_num();
          auto adj_face = m.adjacent_face(cv, f);
          if (adj_face.cv == size_type(-1)) {
            GMM_WARNING2("Adjacent face not found, "
                         "probably an non-interior face");
            inin.ctx.invalid_convex_num();
          } else {
            gauss_pt_corresp gpc;
            gpc.pgt1 = m.trans_of_convex(cv);
            gpc.pgt2 = m.trans_of_convex(adj_face.cv);
            gpc.pai = pai;
            auto inds_pt1 = m.ind_points_of_face_of_convex(cv, f);
            auto inds_pt2 = m.ind_points_of_face_of_convex(adj_face.cv,
                                                           adj_face.f);
            auto str1 = gpc.pgt1->structure();
            auto str2 = gpc.pgt2->structure();
            size_type nbptf1 = str1->nb_points_of_face(f);
            size_type nbptf2 = str2->nb_points_of_face(adj_face.f);
            gpc.nodes.resize(nbptf1*2);
            for (size_type i = 0; i < nbptf1; ++i)  {
              gpc.nodes[2*i] = str1->ind_points_of_face(f)[i];
              bool found = false;
              for (size_type j = 0; j < nbptf2; ++j) {
                if (inds_pt2[j] == inds_pt1[i]) {
                  gpc.nodes[2*i+1] = str2->ind_points_of_face(adj_face.f)[j];
                  found = true;
                  break;
                }
              }
              GMM_ASSERT1(found, "Internal error");
            }
            bgeot::pstored_point_tab pspt = 0;
            auto itm = neighbor_corresp.find(gpc);
            if (itm != neighbor_corresp.end()) {
              pspt = itm->second;
            } else {
              size_type nbpt = pai->nb_points_on_face(f);
              bgeot::geotrans_inv_convex gic;
              gic.init(m.points_of_convex(adj_face.cv), gpc.pgt2);
              size_type first_ind = pai->ind_first_point_on_face(f);
              const bgeot::stored_point_tab
                &spt = *(pai->pintegration_points());
              base_matrix G;
              m.points_of_convex(cv, G);
              fem_interpolation_context ctx_x(gpc.pgt1, 0, spt[0], G, cv, f);
              std::vector<base_node> P_ref(nbpt);

              for (size_type i = 0; i < nbpt; ++i) {
                ctx_x.set_xref(spt[first_ind+i]);
                bool converged = true;
                gic.invert(ctx_x.xreal(), P_ref[i], converged);
                bool is_in = (gpc.pgt2->convex_ref()->is_in(P_ref[i]) < 1E-4);
                GMM_ASSERT1(is_in && converged,"Geometric transformation "
                            "inversion has failed in neighbor transformation");
              }
              pspt = store_point_tab(P_ref);
              neighbor_corresp[gpc] = pspt;
            }
            m.points_of_convex(adj_face.cv, inin.G);
            bgeot::pgeotrans_precomp pgp = gp_pool(gpc.pgt2, pspt);
            inin.ctx.change(pgp, 0, 0, inin.G, adj_face.cv, adj_face.f);
          }
        }
      }

      if (inin.ctx.have_pgp() && inin.ctx.is_convex_num_valid()) {
        inin.ctx.set_ii(ipt);
        inin.pt_type = 1;
        inin.has_ctx = true;
        inin.pt_y = inin.ctx.xreal();
        inin.Normal = bgeot::compute_normal(inin.ctx, inin.ctx.face_num());
        gmm::scale(inin.Normal, 1.0/gmm::vect_norm2(inin.Normal));
        inin.m = &m;
      } else {
        base_node P_ref;
        size_type cv;
        short_type face_num;
        gmm::clear(inin.Normal);
        inin.pt_type = trans->transform(workspace, m, ctx, dummy_normal,
                                        &(inin.m), cv, face_num, P_ref,
                                        dummy_normal, inin.derivatives,
                                        false);
        if (inin.pt_type) {
          if (cv != size_type(-1)) {
            inin.m->points_of_convex(cv, inin.G);
            inin.ctx.change((inin.m)->trans_of_convex(cv),
                            0, P_ref, inin.G, cv, face_num);
            inin.has_ctx = true;
            if (face_num != short_type(-1)) {
              inin.Normal = bgeot::compute_normal(inin.ctx, face_num);
              gmm::scale(inin.Normal, 1.0/gmm::vect_norm2(inin.Normal));
            } else
              inin.Normal.resize(0);
            inin.pt_y = inin.ctx.xreal();
          } else {
            inin.ctx.invalid_convex_num();
            inin.pt_y = P_ref;
            inin.has_ctx = false;
          }
        } else {
          inin.ctx.invalid_convex_num();
          inin.Normal.resize(0);
          inin.pt_y.resize(0);
          inin.has_ctx = false;
        }
      }
      GA_DEBUG_INFO("Instruction: end of call neighbor interpolate "
                    "transformation");
      return 0;
    }
    ga_instruction_neighbor_transformation_call
    (const ga_workspace &w, ga_instruction_set::interpolate_info &i,
     pinterpolate_transformation t, fem_interpolation_context &ctxx,
     const mesh &mm, size_type &ipt_, papprox_integration &pai_,
     bgeot::geotrans_precomp_pool &gp_pool_,
     std::map<gauss_pt_corresp, bgeot::pstored_point_tab> &neighbor_corresp_)
      : workspace(w), inin(i), trans(t), ctx(ctxx), m(mm),
        ipt(ipt_), pai(pai_), gp_pool(gp_pool_),
        neighbor_corresp(neighbor_corresp_) {}
  };


  struct ga_instruction_scalar_assembly : public ga_instruction {
    const base_tensor &t;
    scalar_type &E, &coeff;
     virtual int exec() {
      GA_DEBUG_INFO("Instruction: scalar term assembly");
      E += t[0] * coeff;
      return 0;
     }
    ga_instruction_scalar_assembly(base_tensor &t_, scalar_type &E_,
                                   scalar_type &coeff_)
      : t(t_), E(E_), coeff(coeff_) {}
  };

  struct ga_instruction_vector_assembly_mf : public ga_instruction
  {
    const base_tensor &t;
    base_vector &VI, &Vi;
    const fem_interpolation_context &ctx;
    const gmm::sub_interval *const&I, *const I__;
    const mesh_fem *const&mf, *const mf__;
    const bool &reduced_mf;
    const scalar_type &coeff;
    const size_type &nbpt, &ipt;
    base_vector elem;
    const bool interpolate;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: vector term assembly for fem variable");
      bool empty_weight = (coeff == scalar_type(0));
      if (ipt == 0 || interpolate) {
        if (empty_weight) elem.resize(0);
        elem.resize(t.size());
        if (!empty_weight)
          copy_scaled_4(t, coeff, elem);
      } else if (!empty_weight)
        // gmm::add(gmm::scaled(t.as_vector(), coeff), elem);
        add_scaled_4(t, coeff, elem);

      if (ipt == nbpt-1 || interpolate) { // finalize
        GA_DEBUG_ASSERT(mf, "Internal error");
        if (!ctx.is_convex_num_valid()) return 0;
        size_type cv_1 = ctx.convex_num();
        size_type qmult = mf->get_qdim();
        if (qmult > 1) qmult /= mf->fem_of_element(cv_1)->target_dim();
        base_vector &V = reduced_mf ? Vi : VI;
        GA_DEBUG_ASSERT(V.size() >= I->first() + mf->nb_basic_dof(),
                        "Bad assembly vector size " << V.size() << ">=" <<
                        I->first() << "+"<< mf->nb_basic_dof());
        auto itr = elem.cbegin();
        auto itw = V.begin() + I->first();
        for (const auto &dof : mf->ind_scalar_basic_dof_of_element(cv_1))
          for (size_type q = 0; q < qmult; ++q)
              *(itw+dof+q) += *itr++;
        GMM_ASSERT1(itr == elem.end(), "Internal error");
      }
      return 0;
    }

    ga_instruction_vector_assembly_mf
    (const base_tensor &t_, base_vector &VI_, base_vector &Vi_,
     const fem_interpolation_context &ctx_,
     const gmm::sub_interval *&I_, const mesh_fem *&mf_,
     const bool &reduced_mf_,
     const scalar_type &coeff_, const size_type &nbpt_, const size_type &ipt_,
     bool interpolate_)
    : t(t_), VI(VI_), Vi(Vi_), ctx(ctx_),
      I(I_), I__(nullptr), mf(mf_), mf__(nullptr), reduced_mf(reduced_mf_),
      coeff(coeff_), nbpt(nbpt_), ipt(ipt_), interpolate(interpolate_) {}

    ga_instruction_vector_assembly_mf
    (const base_tensor &t_, base_vector &V_,
     const fem_interpolation_context &ctx_,
     const gmm::sub_interval &I_, const mesh_fem &mf_,
     const scalar_type &coeff_, const size_type &nbpt_, const size_type &ipt_,
     bool interpolate_)
    : t(t_), VI(V_), Vi(V_), ctx(ctx_),
      I(I__), I__(&I_), mf(mf__), mf__(&mf_), reduced_mf(false_),
      coeff(coeff_), nbpt(nbpt_), ipt(ipt_), interpolate(interpolate_) {}
  protected:
    const bool false_=false;
  };

  struct ga_instruction_vector_assembly_imd : public ga_instruction {
    const base_tensor &t;
    base_vector &V;
    const fem_interpolation_context &ctx;
    const gmm::sub_interval &I;
    const im_data &imd;
    scalar_type &coeff;
    const size_type &ipt;
    const bool initialize;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: vector term assembly for im_data variable");
      size_type cv = ctx.convex_num();
      size_type i = t.size() * imd.filtered_index_of_point(cv, ctx.ii());
      GMM_ASSERT1(i+t.size() <= I.size(),
                  "Internal error "<<i<<"+"<<t.size()<<" <= "<<I.size());
      auto itw = V.begin() + I.first() + i;
      if (initialize)
        for (const auto &val : t.as_vector())
          *itw++ = coeff*val;
      else
        for (const auto &val : t.as_vector())
          *itw++ += coeff*val;
      return 0;
    }
    ga_instruction_vector_assembly_imd
    (const base_tensor &t_, base_vector &V_,
     const fem_interpolation_context &ctx_, const gmm::sub_interval &I_,
     const im_data &imd_, scalar_type &coeff_, const size_type &ipt_,
     bool initialize_=false)
    : t(t_), V(V_), ctx(ctx_), I(I_), imd(imd_), coeff(coeff_), ipt(ipt_),
      initialize(initialize_)
    {}
  };

  struct ga_instruction_vector_assembly : public ga_instruction {
    const base_tensor &t;
    base_vector &V;
    const gmm::sub_interval &I;
    scalar_type &coeff;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: vector term assembly for "
                    "fixed size variable");
      gmm::add(gmm::scaled(t.as_vector(), coeff), gmm::sub_vector(V, I));
      return 0;
    }
    ga_instruction_vector_assembly(const base_tensor &t_, base_vector &V_,
                                   const gmm::sub_interval &I_,
                                   scalar_type &coeff_)
      : t(t_), V(V_), I(I_), coeff(coeff_) {}
  };

  struct ga_instruction_assignment : public ga_instruction {
    const base_tensor &t;
    base_vector &V;
    const fem_interpolation_context &ctx;
    const im_data *imd;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Assignement to im_data");
      imd->set_tensor(V, ctx.convex_num(), ctx.ii(), t);
      return 0;
    }
    ga_instruction_assignment(const base_tensor &t_, base_vector &V_,
                              const fem_interpolation_context &ctx_,
                              const im_data *imd_)
      : t(t_), V(V_), ctx(ctx_), imd(imd_) {}
  };

  struct ga_instruction_extract_residual_on_imd_dofs : public ga_instruction {
    base_tensor &t;
    const base_vector &V;
    const fem_interpolation_context &ctx;
    const gmm::sub_interval &I;
    const im_data &imd;
    const size_type &ipt;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: extract residual for im_data variable");
      size_type ifirst = I.first();
      size_type cv = ctx.convex_num();
      size_type i = t.size() * imd.filtered_index_of_point(cv, ctx.ii());
      GMM_ASSERT1(i+t.size() <= I.size(),
                  "Internal error "<<i<<"+"<<t.size()<<" <= "<<I.size());
      for (auto &&val : t.as_vector())
        val = V[ifirst+(i++)];
      return 0;
    }
    ga_instruction_extract_residual_on_imd_dofs
    (base_tensor &t_, const base_vector &V_,
     const fem_interpolation_context &ctx_, const gmm::sub_interval &I_,
     const im_data &imd_, const size_type &ipt_)
    : t(t_), V(V_), ctx(ctx_), I(I_), imd(imd_), ipt(ipt_)
    {}
  };


  template <class MAT>
  inline void add_elem_matrix
  (MAT &K, const std::vector<size_type> &dofs1,
   const std::vector<size_type> &dofs2, std::vector<size_type> &/*dofs1_sort*/,
   const base_vector &elem, scalar_type threshold, size_type /* N */) {

    base_vector::const_iterator it = elem.cbegin();
    for (const size_type &dof2 : dofs2)
      for (const size_type &dof1 : dofs1) {
        if (gmm::abs(*it) > threshold)
          K(dof1, dof2) += *it;
        ++it;
      }
  }

  // static const std::vector<size_type> *the_indto_sort;
  // int compare_my_indices(const void *a, const void *b) {
  //   size_type aa = *((const size_type *)(a));
  //   size_type bb = *((const size_type *)(b));
  //   return  int((*the_indto_sort)[aa]) - int((*the_indto_sort)[bb]);
  // }

  inline void add_elem_matrix
  (gmm::col_matrix<gmm::rsvector<scalar_type>> &K,
   const std::vector<size_type> &dofs1, const std::vector<size_type> &dofs2,
   std::vector<size_type> &dofs1_sort,
   const base_vector &elem, scalar_type threshold, size_type N) {

    size_type s1 = dofs1.size();

    dofs1_sort.resize(s1);
    for (size_type i = 0; i < s1; ++i) { // insertion sort
      size_type j = i, k = j-1;
      while (j > 0 && dofs1[i] < dofs1[dofs1_sort[k]])
        { dofs1_sort[j] = dofs1_sort[k]; j--; k--; }
      dofs1_sort[j] = i;
    }

    // dofs1_sort.resize(s1); // test with qsort: not faster in the tested cases
    // for (size_type i = 0; i < s1; ++i) dofs1_sort[i] = i;
    // the_indto_sort = &dofs1;
    // qsort(&(dofs1_sort[0]), s1, sizeof(size_type), compare_my_indices);

    gmm::elt_rsvector_<scalar_type> ev;

    size_type maxest = (N+1) * s1;
    base_vector::const_iterator it = elem.cbegin();
    bool first(true);
    for (const size_type &dof2 : dofs2) { // Iteration on columns
      if (first) first = false;
      else it += s1;
      std::vector<gmm::elt_rsvector_<scalar_type>> &col = K[dof2];
      size_type nb = col.size();

      if (nb == 0) {
        col.reserve(maxest);
        for (size_type k : dofs1_sort) {
          ev.e = *(it+k);
          if (gmm::abs(ev.e) > threshold) {
            ev.c=dofs1[k];
            col.push_back(ev);
          }
        }
      } else { // column merge
        size_type ind = 0;
        for (size_type k : dofs1_sort) {
          ev.e = *(it+k);
          if (gmm::abs(ev.e) > threshold) {
            ev.c = dofs1[k];

            size_type count = nb - ind, step, l;
            while (count > 0) {
              step = count / 2;
              l = ind + step;
              if (col[l].c < ev.c) {
                ind = ++l;
                count -= step + 1;
              }
              else
                count = step;
            }

            auto itc = col.begin() + ind;
            if (ind != nb && itc->c == ev.c)
              itc->e += ev.e;
            else {
              if (nb - ind > 1300)
                GMM_WARNING2("Inefficient addition of element in rsvector with "
                             << col.size() - ind << " non-zero entries");
              col.push_back(ev);
              if (ind != nb) {
                itc = col.begin() + ind;
                auto ite = col.end();
                --ite;
                auto itee = ite;
                for (; ite != itc; --ite) { --itee; *ite = *itee; }
                *itc = ev;
              }
              ++nb;
            }
            ++ind;
          }
        }
      }
    }
  }


  inline void add_elem_matrix_contiguous_rows
  (gmm::col_matrix<gmm::rsvector<scalar_type>> &K,
   const size_type &i1, const size_type &s1,
   const std::vector<size_type> &dofs2,
   const base_vector &elem, scalar_type threshold) {

    gmm::elt_rsvector_<scalar_type> ev;

    base_vector::const_iterator it = elem.cbegin();
    bool first(true);
    for (const size_type &dof2 : dofs2) { // Iteration on columns
      if (first) first = false;
      else it += s1;
      std::vector<gmm::elt_rsvector_<scalar_type>> &col = K[dof2];
      size_type nb = col.size();

      if (nb == 0) {
        col.reserve(s1);
        for (size_type i = 0; i < s1; ++i) {
          ev.e = *(it+i);
          if (gmm::abs(ev.e) > threshold) {
            ev.c = i1 + i;
            col.push_back(ev);
          }
        }
      } else { // column merge (can be optimized for a contiguous range)
        size_type ind = 0;
        for (size_type i = 0; i < s1; ++i) {
          ev.e = *(it+i);
          if (gmm::abs(ev.e) > threshold) {
            ev.c = i1 + i;

            size_type count = nb - ind, step, l;
            while (count > 0) {
              step = count / 2;
              l = ind + step;
              if (col[l].c < ev.c) {
                ind = ++l;
                count -= step + 1;
              }
              else
                count = step;
            }

            auto itc = col.begin() + ind;
            if (ind != nb && itc->c == ev.c)
              itc->e += ev.e;
            else {
              if (nb - ind > 1300)
                GMM_WARNING2("Inefficient addition of element in rsvector with "
                             << col.size() - ind << " non-zero entries");
              col.push_back(ev);
              if (ind != nb) {
                itc = col.begin() + ind;
                auto ite = col.end();
                --ite;
                auto itee = ite;
                for (; ite != itc; --ite) { --itee; *ite = *itee; }
                *itc = ev;
              }
              ++nb;
            }
            ++ind;
          }
        }
      }
    }
  }

  inline void populate_dofs_vector
  (std::vector<size_type> &dofs,
   const size_type &size, const size_type &ifirst, const size_type &qmult,
   const getfem::mesh::ind_set &mfdofs)
  {
    dofs.assign(size, ifirst);
    auto itd = dofs.begin();
    if (qmult == 1)
      for (const auto &dof : mfdofs) *itd++ += dof;
    else
      for (const auto &dof : mfdofs)
        for (size_type q = 0; q < qmult; ++q) *itd++ += dof + q;
  }

  inline void populate_dofs_vector // special case for qmult == 1
  (std::vector<size_type> &dofs, const size_type &size, const size_type &ifirst,
   const getfem::mesh::ind_set &mfdofs)
  {
    dofs.assign(size, ifirst);
    auto itd = dofs.begin();
    for (const auto &dof : mfdofs) *itd++ += dof;
  }


  inline void populate_contiguous_dofs_vector
  (std::vector<size_type> &dofs, const size_type &size, const size_type &ifirst)
  {
    dofs.assign(size, ifirst);
    for (size_type i=0; i < size; ++i) dofs[i] += i;
  }

  struct ga_instruction_matrix_assembly_base : public ga_instruction {
    const base_tensor &t;
    const fem_interpolation_context &ctx1, &ctx2;
    const scalar_type &alpha1, &alpha2, &coeff;
    const size_type &nbpt, &ipt;
    base_vector elem;
    bool interpolate;
    std::vector<size_type> dofs1, dofs2, dofs1_sort;
    void add_tensor_to_element_matrix(bool initialize, bool empty_weight) {
      if (initialize) {
        if (empty_weight) elem.resize(0);
        elem.resize(t.size());
        if (!empty_weight)
          copy_scaled_4(t, coeff*alpha1*alpha2, elem);
      } else if (!empty_weight)
        // gmm::add(gmm::scaled(t.as_vector(), coeff*alpha1*alpha2), elem);
        // Faster than a daxpy blas call on my config
        add_scaled_4(t, coeff*alpha1*alpha2, elem);
    }
    ga_instruction_matrix_assembly_base
    (const base_tensor &t_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const scalar_type &a1, const scalar_type &a2, const scalar_type &coeff_,
     const size_type &nbpt_, const size_type &ipt_, bool interpolate_)
      : t(t_), ctx1(ctx1_), ctx2(ctx2_), alpha1(a1), alpha2(a2),
        coeff(coeff_), nbpt(nbpt_), ipt(ipt_), interpolate(interpolate_),
        dofs1(0), dofs2(0), dofs1_sort(0)
    {}
  protected:
    const bool false_=false;
    const size_type zero_=0;
  };


  struct ga_instruction_matrix_assembly_mf_mf
    : public ga_instruction_matrix_assembly_base
  {
    model_real_sparse_matrix &Krr, &Kru, &Kur, &Kuu;
    const gmm::sub_interval *const&I1, *const&I2, *const I1__, *const I2__;
    const mesh_fem *const&mf1, *const&mf2, *const mf1__, *const mf2__;
    const bool &reduced_mf1, &reduced_mf2; // refs to mf1/2->is_reduced()
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: matrix term assembly mf-mf");
      if (!ctx1.is_convex_num_valid() || !ctx2.is_convex_num_valid()) return 0;

      bool initialize = (ipt == 0 || interpolate);
      bool empty_weight = (coeff == scalar_type(0));
      add_tensor_to_element_matrix(initialize, empty_weight); // t --> elem

      if (ipt == nbpt-1 || interpolate) { // finalize
        model_real_sparse_matrix &K = reduced_mf1 ? (reduced_mf2 ? Kuu : Kur)
                                                  : (reduced_mf2 ? Kru : Krr);
        GA_DEBUG_ASSERT(I1->size() && I2->size(), "Internal error");

        scalar_type ninf = gmm::vect_norminf(elem);
        if (ninf == scalar_type(0)) return 0;

        size_type s1 = t.sizes()[0], s2 = t.sizes()[1];
        size_type cv1 = ctx1.convex_num(), cv2 = ctx2.convex_num();
        size_type ifirst1 = I1->first(), ifirst2 = I2->first();

        size_type N = ctx1.N();
        size_type qmult1 = mf1->get_qdim();
        if (qmult1 > 1) qmult1 /= mf1->fem_of_element(cv1)->target_dim();
        populate_dofs_vector(dofs1, s1, ifirst1, qmult1,        // --> dofs1
                             mf1->ind_scalar_basic_dof_of_element(cv1));
        if (mf1 == mf2 && cv1 == cv2) {
          if (ifirst1 == ifirst2) {
            add_elem_matrix(K, dofs1, dofs1, dofs1_sort, elem, ninf*1E-14, N);
          } else {
            populate_dofs_vector(dofs2, dofs1.size(), ifirst2 - ifirst1, dofs1);
            add_elem_matrix(K, dofs1, dofs2, dofs1_sort, elem, ninf*1E-14, N);
          }
        } else {
          N = std::max(N, ctx2.N());
          size_type qmult2 = mf2->get_qdim();
          if (qmult2 > 1) qmult2 /= mf2->fem_of_element(cv2)->target_dim();
          populate_dofs_vector(dofs2, s2, ifirst2, qmult2,        // --> dofs2
                               mf2->ind_scalar_basic_dof_of_element(cv2));
          add_elem_matrix(K, dofs1, dofs2, dofs1_sort, elem, ninf*1E-14, N);
        }
      }
      return 0;
    }

    ga_instruction_matrix_assembly_mf_mf
    (const base_tensor &t_,
     model_real_sparse_matrix &Krr_, model_real_sparse_matrix &Kru_,
     model_real_sparse_matrix &Kur_, model_real_sparse_matrix &Kuu_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const ga_instruction_set::variable_group_info &vgi1,
     const ga_instruction_set::variable_group_info &vgi2,
     const scalar_type &coeff_, const size_type &nbpt_, const size_type &ipt_,
     bool interpolate_)
      : ga_instruction_matrix_assembly_base
        (t_, ctx1_, ctx2_, vgi1.alpha, vgi2.alpha, coeff_, nbpt_, ipt_,
         interpolate_),
      Krr(Krr_), Kru(Kru_), Kur(Kur_), Kuu(Kuu_),
      I1(vgi1.I), I2(vgi2.I), I1__(nullptr), I2__(nullptr),
      mf1(vgi1.mf), mf2(vgi2.mf), mf1__(nullptr), mf2__(nullptr),
      reduced_mf1(vgi1.reduced_mf), reduced_mf2(vgi2.reduced_mf) {}

    ga_instruction_matrix_assembly_mf_mf
    (const base_tensor &t_,
     model_real_sparse_matrix &Kxr_, model_real_sparse_matrix &Kxu_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const gmm::sub_interval &I1_, const mesh_fem &mf1_, const scalar_type &a1,
     const ga_instruction_set::variable_group_info &vgi2,
     const scalar_type &coeff_, const size_type &nbpt_, const size_type &ipt_,
     bool interpolate_)
      : ga_instruction_matrix_assembly_base
        (t_, ctx1_, ctx2_, a1, vgi2.alpha, coeff_, nbpt_, ipt_, interpolate_),
      Krr(Kxr_), Kru(Kxu_), Kur(Kxr_), Kuu(Kxu_),
      I1(I1__), I2(vgi2.I), I1__(&I1_), I2__(nullptr),
      mf1(mf1__), mf2(vgi2.mf), mf1__(&mf1_), mf2__(nullptr),
      reduced_mf1(false_), reduced_mf2(vgi2.reduced_mf) {}

    ga_instruction_matrix_assembly_mf_mf
    (const base_tensor &t_,
     model_real_sparse_matrix &Krx_, model_real_sparse_matrix &Kux_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const ga_instruction_set::variable_group_info &vgi1,
     const gmm::sub_interval &I2_, const mesh_fem &mf2_, const scalar_type &a2,
     const scalar_type &coeff_, const size_type &nbpt_, const size_type &ipt_,
     bool interpolate_)
      : ga_instruction_matrix_assembly_base
        (t_, ctx1_, ctx2_, vgi1.alpha, a2, coeff_, nbpt_, ipt_, interpolate_),
      Krr(Krx_), Kru(Krx_), Kur(Kux_), Kuu(Kux_),
      I1(vgi1.I), I2(I2__), I1__(nullptr), I2__(&I2_),
      mf1(vgi1.mf), mf2(mf2__), mf1__(nullptr), mf2__(&mf2_),
      reduced_mf1(vgi1.reduced_mf), reduced_mf2(false_) {}

    ga_instruction_matrix_assembly_mf_mf
    (const base_tensor &t_, model_real_sparse_matrix &K_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const gmm::sub_interval &I1_, const mesh_fem &mf1_, const scalar_type &a1,
     const gmm::sub_interval &I2_, const mesh_fem &mf2_, const scalar_type &a2,
     const scalar_type &coeff_, const size_type &nbpt_, const size_type &ipt_,
     bool interpolate_)
      : ga_instruction_matrix_assembly_base
        (t_, ctx1_, ctx2_, a1, a2, coeff_, nbpt_, ipt_, interpolate_),
      Krr(K_), Kru(K_), Kur(K_), Kuu(K_),
      I1(I1__), I2(I2__), I1__(&I1_), I2__(&I2_),
      mf1(mf1__), mf2(mf2__), mf1__(&mf1_), mf2__(&mf2_),
      reduced_mf1(false_), reduced_mf2(false_) {}
  };


  struct ga_instruction_matrix_assembly_imd_mf
    : public ga_instruction_matrix_assembly_base
  {
    model_real_sparse_matrix &Kxr, &Kxu;
    const gmm::sub_interval *I1, *I2__, * const &I2;
    const im_data *imd1;
    const mesh_fem * const mf2__, * const &mf2;
    const bool &reduced_mf2; // ref to mf2->is_reduced()
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: matrix term assembly "
                    "(imdata or fixed size)-mf");
      if (!ctx1.is_convex_num_valid() || !ctx2.is_convex_num_valid()) return 0;

      bool empty_weight = (coeff == scalar_type(0));
      add_tensor_to_element_matrix(true, empty_weight); // t --> elem

      scalar_type ninf = gmm::vect_norminf(elem);
      if (ninf == scalar_type(0)) return 0;

      model_real_sparse_matrix &K = reduced_mf2 ? Kxu : Kxr;
      GA_DEBUG_ASSERT(I1->size() && I2->size(), "Internal error");
      size_type s1 = t.sizes()[0], s2 = t.sizes()[1];
      size_type cv1 = ctx1.convex_num(), cv2 = ctx2.convex_num();
      size_type ifirst1 = I1->first(), ifirst2 = I2->first();
      if (imd1) ifirst1 += s1 * imd1->filtered_index_of_point(cv1, ctx1.ii());

      populate_contiguous_dofs_vector(dofs1, s1, ifirst1); // --> dofs1
      size_type qmult2 = mf2->get_qdim();
      if (qmult2 > 1) qmult2 /= mf2->fem_of_element(cv2)->target_dim();
      populate_dofs_vector(dofs2, s2, ifirst2, qmult2,     // --> dofs2
                           mf2->ind_scalar_basic_dof_of_element(cv2));
      add_elem_matrix(K, dofs1, dofs2, dofs1_sort, elem, ninf*1E-14, ctx2.N());
      return 0;
    }

    ga_instruction_matrix_assembly_imd_mf
    (const base_tensor &t_,
     model_real_sparse_matrix &Kxr_, model_real_sparse_matrix &Kxu_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const gmm::sub_interval &I1_, const im_data *imd1_, const scalar_type &a1,
     const ga_instruction_set::variable_group_info &vgi2,
     const scalar_type &coeff_, const size_type &ipt_)
      : ga_instruction_matrix_assembly_base
        (t_, ctx1_, ctx2_, a1, vgi2.alpha, coeff_, zero_, ipt_, false),
        Kxr(Kxr_), Kxu(Kxu_), I1(&I1_), I2__(nullptr), I2(vgi2.I),
        imd1(imd1_), mf2__(nullptr), mf2(vgi2.mf), reduced_mf2(vgi2.reduced_mf)
    {}

    ga_instruction_matrix_assembly_imd_mf
    (const base_tensor &t_, model_real_sparse_matrix &K_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const gmm::sub_interval &I1_, const im_data *imd1_, const scalar_type &a1,
     const gmm::sub_interval &I2_, const mesh_fem &mf2_, const scalar_type &a2,
     const scalar_type &coeff_, const size_type &ipt_)
      : ga_instruction_matrix_assembly_base
        (t_, ctx1_, ctx2_, a1, a2, coeff_, zero_, ipt_, false),
        Kxr(K_), Kxu(K_), I1(&I1_), I2__(&I2_), I2(I2__),
        imd1(imd1_), mf2__(&mf2_), mf2(mf2__), reduced_mf2(false_) {}
  };

  struct ga_instruction_matrix_assembly_mf_imd
    : public ga_instruction_matrix_assembly_base
  {
    model_real_sparse_matrix &Krx, &Kux;
    const gmm::sub_interval * const &I1, *const I1__, *I2;
    const mesh_fem * const &mf1, *const mf1__;
    const bool &reduced_mf1; // ref to mf1->is_reduced()
    const im_data *imd2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: matrix term assembly "
                    "mf-(imdata or fixed size)");
      if (!ctx1.is_convex_num_valid() || !ctx2.is_convex_num_valid()) return 0;

      bool empty_weight = (coeff == scalar_type(0));
      add_tensor_to_element_matrix(true, empty_weight); // t --> elem

      scalar_type ninf = gmm::vect_norminf(elem);
      if (ninf == scalar_type(0)) return 0;

      model_real_sparse_matrix &K = reduced_mf1 ? Kux : Krx;
      GA_DEBUG_ASSERT(I1->size() && I2->size(), "Internal error");
      size_type s1 = t.sizes()[0], s2 = t.sizes()[1];
      size_type cv1 = ctx1.convex_num(), cv2 = ctx2.convex_num();
      size_type ifirst1 = I1->first(), ifirst2 = I2->first();
      if (imd2) ifirst2 += s2 * imd2->filtered_index_of_point(cv2, ctx2.ii());

      size_type qmult1 = mf1->get_qdim();
      if (qmult1 > 1) qmult1 /= mf1->fem_of_element(cv1)->target_dim();
      populate_dofs_vector(dofs1, s1, ifirst1, qmult1,     // --> dofs1
                           mf1->ind_scalar_basic_dof_of_element(cv1));
      populate_contiguous_dofs_vector(dofs2, s2, ifirst2); // --> dofs2
      add_elem_matrix(K, dofs1, dofs2, dofs1_sort, elem, ninf*1E-14, ctx1.N());
      return 0;
    }

    ga_instruction_matrix_assembly_mf_imd
    (const base_tensor &t_,
     model_real_sparse_matrix &Krx_, model_real_sparse_matrix &Kux_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const ga_instruction_set::variable_group_info &vgi1,
     const gmm::sub_interval &I2_, const im_data *imd2_, const scalar_type &a2,
     const scalar_type &coeff_, const size_type &ipt_)
      : ga_instruction_matrix_assembly_base
        (t_, ctx1_, ctx2_, vgi1.alpha, a2, coeff_, zero_, ipt_, false),
        Krx(Krx_), Kux(Kux_), I1(vgi1.I), I1__(nullptr), I2(&I2_),
        mf1(vgi1.mf), mf1__(nullptr), reduced_mf1(vgi1.reduced_mf), imd2(imd2_)
    {}

    ga_instruction_matrix_assembly_mf_imd
    (const base_tensor &t_, model_real_sparse_matrix &K_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const gmm::sub_interval &I1_, const mesh_fem &mf1_, const scalar_type &a1,
     const gmm::sub_interval &I2_, const im_data *imd2_, const scalar_type &a2,
     const scalar_type &coeff_, const size_type &ipt_)
      : ga_instruction_matrix_assembly_base
        (t_, ctx1_, ctx2_, a1, a2, coeff_, zero_, ipt_, false),
        Krx(K_), Kux(K_), I1(I1__), I1__(&I1_), I2(&I2_),
        mf1(mf1__), mf1__(&mf1_), reduced_mf1(false_), imd2(imd2_) {}
  };



  struct ga_instruction_matrix_assembly_imd_imd
    : public ga_instruction_matrix_assembly_base
  {
    model_real_sparse_matrix &K;
    const gmm::sub_interval &I1, &I2;
    const im_data *imd1, *imd2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: matrix term assembly "
                    "(imdata or fixed size)-(imdata or fixed size)");
      GA_DEBUG_ASSERT(I1.size() && I2.size(), "Internal error");

      bool empty_weight = (coeff == scalar_type(0));
      add_tensor_to_element_matrix(true, empty_weight); // t --> elem

      scalar_type ninf = gmm::vect_norminf(elem);
      if (ninf == scalar_type(0)) return 0;

      size_type s1 = t.sizes()[0], s2 = t.sizes()[1];
      size_type ifirst1 = I1.first(), ifirst2 = I2.first();
      if (imd1)
        ifirst1 += s1 * imd1->filtered_index_of_point(ctx1.convex_num(), ctx1.ii());
      if (imd2)
        ifirst2 += s2 * imd2->filtered_index_of_point(ctx2.convex_num(), ctx2.ii());

      populate_contiguous_dofs_vector(dofs2, s2, ifirst2);
      add_elem_matrix_contiguous_rows(K, ifirst1, s1, dofs2, elem, ninf*1E-14);
      return 0;
    }
    ga_instruction_matrix_assembly_imd_imd
    (const base_tensor &t_, model_real_sparse_matrix &K_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const gmm::sub_interval &I1_, const im_data *imd1_, const scalar_type &a1,
     const gmm::sub_interval &I2_, const im_data *imd2_, const scalar_type &a2,
     const scalar_type &coeff_, const size_type &ipt_)
      : ga_instruction_matrix_assembly_base
        (t_, ctx1_, ctx2_, a1, a2, coeff_, zero_, ipt_, false),
        K(K_), I1(I1_), I2(I2_), imd1(imd1_), imd2(imd2_) {}
  };


  struct ga_instruction_matrix_assembly_standard_scalar
    : public ga_instruction_matrix_assembly_base
  {
    model_real_sparse_matrix &K;
    const gmm::sub_interval &I1, &I2;
    const mesh_fem *pmf1, *pmf2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: matrix term assembly for standard "
                    "scalar fems");
      if (ipt == 0) {
        elem.resize(t.size());
        // gmm::copy(gmm::scaled(t.as_vector(), coeff*alpha1*alpha2), elem);
        copy_scaled_4(t, coeff*alpha1*alpha2, elem);
      } else
        // gmm::add(gmm::scaled(t.as_vector(), coeff*alpha1*alpha2), elem);
        // Faster than a daxpy blas call on my config
        add_scaled_4(t, coeff*alpha1*alpha2, elem);

      if (ipt == nbpt-1) { // finalize
        GA_DEBUG_ASSERT(I1.size() && I2.size(), "Internal error");

        scalar_type ninf = gmm::vect_norminf(elem);
        if (ninf == scalar_type(0)) return 0;

        size_type cv1 = ctx1.convex_num(), cv2 = ctx2.convex_num(), N=ctx1.N();
        if (cv1 == size_type(-1)) return 0;
        auto &ct1 = pmf1->ind_scalar_basic_dof_of_element(cv1);
        GA_DEBUG_ASSERT(ct1.size() == t.sizes()[0], "Internal error");
        populate_dofs_vector(dofs1, ct1.size(), I1.first(), ct1);

        if (pmf2 == pmf1 && cv1 == cv2) {
          if (I1.first() == I2.first()) {
            add_elem_matrix(K, dofs1, dofs1, dofs1_sort, elem, ninf*1E-14, N);
          } else {
            populate_dofs_vector(dofs2, dofs1.size(), I2.first() - I1.first(),
                                 dofs1);
            add_elem_matrix(K, dofs1, dofs2, dofs1_sort, elem, ninf*1E-14, N);
          }
        } else {
          if (cv2 == size_type(-1)) return 0;
          auto &ct2 = pmf2->ind_scalar_basic_dof_of_element(cv2);
          GA_DEBUG_ASSERT(ct2.size() == t.sizes()[1], "Internal error");
          populate_dofs_vector(dofs2, ct2.size(), I2.first(), ct2);
          add_elem_matrix(K, dofs1, dofs2, dofs1_sort, elem, ninf*1E-14, N);
        }
      }
      return 0;
    }
    ga_instruction_matrix_assembly_standard_scalar
    (const base_tensor &t_, model_real_sparse_matrix &K_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const gmm::sub_interval &I1_, const gmm::sub_interval &I2_,
     const mesh_fem *mfn1_, const mesh_fem *mfn2_,
     const scalar_type &a1, const scalar_type &a2, const scalar_type &coeff_,
     const size_type &nbpt_, const size_type &ipt_)
      : ga_instruction_matrix_assembly_base
        (t_, ctx1_, ctx2_, a1, a2, coeff_, nbpt_, ipt_, false),
        K(K_), I1(I1_), I2(I2_), pmf1(mfn1_), pmf2(mfn2_) {}
  };

  struct ga_instruction_matrix_assembly_standard_vector
    : public ga_instruction_matrix_assembly_base
  {
    model_real_sparse_matrix &K;
    const gmm::sub_interval &I1, &I2;
    const mesh_fem *pmf1, *pmf2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: matrix term assembly for standard "
                        "vector fems");
      if (ipt == 0) {
        elem.resize(t.size());
        copy_scaled_8(t, coeff*alpha1*alpha2, elem);
        // gmm::copy(gmm::scaled(t.as_vector(), coeff*alpha1*alpha2), elem);
      } else
        // gmm::add(gmm::scaled(t.as_vector(), coeff*alpha1*alpha2), elem);
        // (Far) faster than a daxpy blas call on my config.
        add_scaled_8(t, coeff*alpha1*alpha2, elem);

      if (ipt == nbpt-1) { // finalize
        GA_DEBUG_ASSERT(I1.size() && I2.size(), "Internal error");

        scalar_type ninf = gmm::vect_norminf(elem);
        if (ninf == scalar_type(0)) return 0;
        size_type s1 = t.sizes()[0], s2 = t.sizes()[1], N = ctx1.N();

        size_type cv1 = ctx1.convex_num(), cv2 = ctx2.convex_num();
        if (cv1 == size_type(-1)) return 0;
        size_type qmult1 = pmf1->get_qdim();
        if (qmult1 > 1) qmult1 /= pmf1->fem_of_element(cv1)->target_dim();
        populate_dofs_vector(dofs1, s1, I1.first(), qmult1,         // --> dofs1
                             pmf1->ind_scalar_basic_dof_of_element(cv1));

        if (pmf2 == pmf1 && cv1 == cv2 && I1.first() == I2.first()) {
          add_elem_matrix(K, dofs1, dofs1, dofs1_sort, elem, ninf*1E-14, N);
        } else {
          if (pmf2 == pmf1 && cv1 == cv2) {
            populate_dofs_vector(dofs2, dofs1.size(), I2.first() - I1.first(),
                                 dofs1);
          } else {
            if (cv2 == size_type(-1)) return 0;
            size_type qmult2 = pmf2->get_qdim();
            if (qmult2 > 1) qmult2 /= pmf2->fem_of_element(cv2)->target_dim();
            populate_dofs_vector(dofs2, s2, I2.first(), qmult2,      // --> dofs2
                                 pmf2->ind_scalar_basic_dof_of_element(cv2));
          }
          add_elem_matrix(K, dofs1, dofs2, dofs1_sort, elem, ninf*1E-14, N);
        }
      }
      return 0;
    }
    ga_instruction_matrix_assembly_standard_vector
    (const base_tensor &t_, model_real_sparse_matrix &K_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const gmm::sub_interval &I1_, const gmm::sub_interval &I2_,
     const mesh_fem *mfn1_, const mesh_fem *mfn2_,
     const scalar_type &a1, const scalar_type &a2, const scalar_type &coeff_,
     const size_type &nbpt_, const size_type &ipt_)
      : ga_instruction_matrix_assembly_base
        (t_, ctx1_, ctx2_, a1, a2, coeff_, nbpt_, ipt_, false),
        K(K_), I1(I1_), I2(I2_), pmf1(mfn1_), pmf2(mfn2_) {}
  };

  template<int QQ>
  struct ga_instruction_matrix_assembly_standard_vector_opt10
    : public ga_instruction_matrix_assembly_base
  {
    model_real_sparse_matrix &K;
    const gmm::sub_interval &I1, &I2;
    const mesh_fem *pmf1, *pmf2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: matrix term assembly for standard "
                    "vector fems optimized for format 10 qdim " << QQ);
      size_type s1_q = QQ*t.sizes()[0];
      size_type ss1 = t.sizes()[0]/QQ, ss2 = t.sizes()[1]/QQ;
      scalar_type e = coeff*alpha1*alpha2;
      if (ipt == 0) {
        elem.resize(ss1*ss2);
        auto itel = elem.begin();
        for (size_type j = 0; j < ss2; ++j) {
          auto it = t.begin() + j*s1_q;
          for (size_type i = 0; i < ss1; ++i, it += QQ)
            *itel++ = (*it) * e;
        }
      } else {
        auto itel = elem.begin();
        for (size_type j = 0; j < ss2; ++j) {
          auto it = t.begin() + j*s1_q;
          for (size_type i = 0; i < ss1; ++i, it += QQ)
            *itel++ += (*it) * e;
        }
      }
      if (ipt == nbpt-1) { // finalize
        GA_DEBUG_ASSERT(I1.size() && I2.size(), "Internal error");

        scalar_type ninf = gmm::vect_norminf(elem) * 1E-14;
        if (ninf == scalar_type(0)) return 0;
        size_type N = ctx1.N();
        size_type cv1 = ctx1.convex_num(), cv2 = ctx2.convex_num();
        size_type i1 = I1.first(), i2 = I2.first();
        if (cv1 == size_type(-1)) return 0;
        populate_dofs_vector(dofs1, ss1, i1,
                             pmf1->ind_scalar_basic_dof_of_element(cv1));
        bool same_dofs(pmf2 == pmf1 && cv1 == cv2 && i1 == i2);

        if (!same_dofs) {
          if (cv2 == size_type(-1)) return 0;
          populate_dofs_vector(dofs2, ss2, i2,
                               pmf2->ind_scalar_basic_dof_of_element(cv2));
        }
        std::vector<size_type> &dofs2_ = same_dofs ? dofs1 : dofs2;
        add_elem_matrix(K, dofs1, dofs2_, dofs1_sort, elem, ninf, N);
        for (size_type i = 0; i < ss1; ++i) (dofs1[i])++;
        if (!same_dofs) for (size_type i = 0; i < ss2; ++i) (dofs2[i])++;
        add_elem_matrix(K, dofs1, dofs2_, dofs1_sort, elem, ninf, N);
        if (QQ >= 3) {
          for (size_type i = 0; i < ss1; ++i) (dofs1[i])++;
          if (!same_dofs) for (size_type i = 0; i < ss2; ++i) (dofs2[i])++;
          add_elem_matrix(K, dofs1, dofs2_, dofs1_sort, elem, ninf, N);
        }
      }
      return 0;
    }

    ga_instruction_matrix_assembly_standard_vector_opt10
    (const base_tensor &t_, model_real_sparse_matrix &Kn_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const gmm::sub_interval &In1_, const gmm::sub_interval &In2_,
     const mesh_fem *mfn1_, const mesh_fem *mfn2_,
     const scalar_type &a1, const scalar_type &a2, const scalar_type &coeff_,
     const size_type &nbpt_, const size_type &ipt_)
      : ga_instruction_matrix_assembly_base
        (t_, ctx1_, ctx2_, a1, a2, coeff_, nbpt_, ipt_, false),
        K(Kn_), I1(In1_), I2(In2_), pmf1(mfn1_), pmf2(mfn2_)
    {
      static_assert(QQ >= 2 && QQ <=3,
                    "Template implemented only for QQ=2 and QQ=3");
    }
  };


  struct ga_instruction_condensation_sub : public ga_instruction {
    // one such instruction is used for every cluster of intercoupled
    // condensed variables
    gmm::dense_matrix<base_tensor *> KQJprime;
    std::vector<base_tensor *> RQprime;
    gmm::dense_matrix<base_tensor const *> KQQloc, KQJloc;
    base_tensor invKqqqq, Kqqjj;
    base_vector Rqq;
    std::vector<std::array<size_type,3>> partQ, partJ;
    const scalar_type &coeff; // &alpha1, &alpha2 ?
    virtual int exec() { // implementation can be optimized
      GA_DEBUG_INFO("Instruction: variable cluster subdiagonal condensation");
      // copy from KQQ to invKqqqq
      for (const auto &qqq1 : partQ) {
        size_type q1 = qqq1[0], qq1start = qqq1[1], qq1end = qqq1[2];
        for (const auto &qqq2 : partQ) {
          size_type q2 = qqq2[0], qq2start = qqq2[1], qq2end = qqq2[2];
          if (KQQloc(q1,q2)) {
            auto itr = KQQloc(q1,q2)->cbegin();
            GMM_ASSERT1(KQQloc(q1,q2)->size()
                        == (qq1end-qq1start)*(qq2end-qq2start),
                        "Internal error");
            for (size_type qq2=qq2start; qq2 < qq2end; ++qq2)
              for (size_type qq1=qq1start; qq1 < qq1end; ++qq1)
                invKqqqq(qq1,qq2) = *itr++;
          }
        }
      }
      // calculate inverse matrix invKqqqq
      bgeot::lu_inverse(&(invKqqqq[0]), invKqqqq.size(0));

      // Resize Kqqjj as primary variable sizes may change dynamically
      size_type prev_j(0);
      for (auto &&jjj : partJ) {
        size_type j=jjj[0];
        size_type new_j(0);
        for (const auto &qqq : partQ) {
          size_type q=qqq[0];
          if (KQJloc(q,j)) {
            if (new_j) {
              GMM_ASSERT1(new_j == KQJloc(q,j)->size(1), "Internal error");
            } else
              new_j = KQJloc(q,j)->size(1);
          }
        }
        // Resize KQJprime submatrices to match KQJloc sizes
        for (const auto &qqq : partQ) {
          size_type q=qqq[0];
          KQJprime(q,j)->adjust_sizes(qqq[2]-qqq[1], new_j);
        }
        jjj[1] = prev_j;
        prev_j += new_j;
        jjj[2] = prev_j;
      }

      Kqqjj.adjust_sizes(partQ.back()[2], partJ.back()[2]);
      gmm::clear(Kqqjj.as_vector());
      gmm::clear(Rqq);

      // multiply invKqqqq with all submatrices in KQJloc and RQprime and store
      // the results in Kqqjj and Rqq
      for (const auto &jjj : partJ) {
        size_type j = jjj[0], jjstart = jjj[1], jjend = jjj[2];
        for (const auto &qqq2 : partQ) {
          size_type q2 = qqq2[0], qq2start = qqq2[1], qq2end = qqq2[2];
          if (KQJloc(q2,j)) {
            auto itr = KQJloc(q2,j)->begin(); // auto &mat = KQJloc(q2,j);
            for (size_type jj=jjstart; jj < jjend; ++jj) {
              for (size_type qq2=qq2start; qq2 < qq2end; ++qq2, ++itr) {
                for (size_type qq1=0; qq1 < partQ.back()[2]; ++qq1) {
                  Kqqjj(qq1,jj) += invKqqqq(qq1,qq2)*(*itr);
                  // Kqqjj(qq1,jj) += invKqq(qq1,qq2)*mat(qq2-qqstart,jj-jjstart);
                } // for qq1
              } // for qq2
            } // for jj
            GMM_ASSERT1(itr == KQJloc(q2,j)->cend(), "Internal error");
          }
        } // in partQ
      } // in partJ
      for (const auto &qqq2 : partQ) {
        size_type q2 = qqq2[0], qq2start = qqq2[1], qq2end = qqq2[2];
        if (RQprime[q2]) {
          auto itr = RQprime[q2]->cbegin();
          for (size_type qq2=qq2start; qq2 < qq2end; ++qq2, ++itr) {
            for (size_type qq1=0; qq1 < invKqqqq.size(0); ++qq1)
              Rqq[qq1] += invKqqqq(qq1,qq2)*(*itr);
          } // for qq2
          GMM_ASSERT1(itr == RQprime[q2]->cend(), "Internal error");
        }
      } // in partQ

      // distribute the results from Kqqjj/Rqq to KQJprime/RQprime
      // submatrices/subvectors
      for (const auto &qqq1 : partQ) {
        size_type q1 = qqq1[0], qq1start = qqq1[1], qq1end = qqq1[2];
        { // writing into RQprime
          auto itw = RQprime[q1]->begin();
          for (size_type qq1=qq1start; qq1 < qq1end; ++qq1)
            *itw++ = Rqq[qq1]/coeff;
        }
        for (const auto &jjj2 : partJ) {
          size_type j2 = jjj2[0], jj2start = jjj2[1], jj2end = jjj2[2];
          auto itw = KQJprime(q1,j2)->begin();
          for (size_type jj2=jj2start; jj2 < jj2end; ++jj2)
            for (size_type qq1=qq1start; qq1 < qq1end; ++qq1)
              *itw++ = Kqqjj(qq1,jj2);
        }
      }
      return 0;
    }

    ga_instruction_condensation_sub(gmm::dense_matrix<base_tensor *> &KQJpr,
                                    std::vector<base_tensor *> &RQpr, // input/output
                                    const gmm::dense_matrix<base_tensor *> &KQQ,
                                    const gmm::dense_matrix<base_tensor *> &KQJ,
                                    const std::set<size_type> &Qset,
                                    const scalar_type &coeff_)
      : KQJprime(KQJpr), RQprime(RQpr), coeff(coeff_)
    {
      // * to const *
      KQQloc.resize(KQQ.nrows(), KQQ.ncols());
      KQJloc.resize(KQJ.nrows(), KQJ.ncols());
      for (size_type i=0; i < KQQ.as_vector().size(); ++i) KQQloc[i] = KQQ[i];
      for (size_type i=0; i < KQJ.as_vector().size(); ++i) KQJloc[i] = KQJ[i];

      for (size_type j=0; j < KQJ.ncols(); ++j)
        for (const size_type &q : Qset)
          if (KQJ(q,j)) {
            partJ.push_back(std::array<size_type,3>{j,0,0});
            break;
          }

      partQ.resize(0);
      for (const size_type &q : Qset)
        partQ.push_back(std::array<size_type,3>{q,0,0});
      size_type prev_q(0);
      for (auto &qqq1 : partQ) {
        size_type q1 = qqq1[0];
        size_type new_q(0);
        for (const size_type &q2 : Qset)
          if (new_q) {
            GMM_ASSERT1(new_q == KQQ(q1,q2)->size(0) &&
                        new_q == KQQ(q2,q1)->size(1), "Internal error");
          } else
            new_q = KQQ(q1,q2)->size(0);
        qqq1[1] = prev_q;
        prev_q += new_q;
        qqq1[2] = prev_q;
      }
      invKqqqq.adjust_sizes(partQ.back()[2], partQ.back()[2]);
      Rqq.resize(partQ.back()[2]);
      // Kqqjj will be resized dynamically due to possible changes in j interval
    }
  };


  struct ga_instruction_condensation_super_K : public ga_instruction {
    base_tensor &Kij;
    std::vector<base_tensor *> KiQ, KQj; // indexed wrt q in Q
    size_type Qsize;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: contribution of condensation to kept part");

      size_type m = KiQ[0]->size(0);
      size_type n = KQj[0]->size(1);
      Kij.adjust_sizes(m,n);
      gmm::clear(Kij.as_vector());
      for (size_type k=0; k < Qsize; ++k) {
        const base_tensor &K1 = *KiQ[k], &K2 = *KQj[k];
        size_type qqsize = K1.size(1);
        GMM_ASSERT1(K1.size(0) == m && K2.size(1) == n && K2.size(0) == qqsize,
                    "Internal error");

        base_tensor::iterator it = Kij.begin();
        for (size_type jj = 0; jj < n; ++jj)
          for (size_type ii = 0; ii < m; ++ii, ++it)
            for (size_type qq = 0; qq < qqsize; ++qq)
              *it -= K1[ii+qq*m] * K2[qq+jj*qqsize];
        GA_DEBUG_ASSERT(it == Kij.end(), "Wrong sizes");
      }
      return 0;
    }
    ga_instruction_condensation_super_K(base_tensor &Kij_,
                                        const std::vector<base_tensor *> KiQ_,
                                        const std::vector<base_tensor *> KQj_)
      : Kij(Kij_), KiQ(KiQ_), KQj(KQj_)
    {
      Qsize = KiQ.size();
      GMM_ASSERT1(KiQ.size() == KQj.size(), "Internal error");
    }
  };

  struct ga_instruction_condensation_super_R : public ga_instruction {
    base_tensor &Ri;
    std::vector<base_tensor *> KiQ, RQpr; // indexed wrt q in Q
    size_type Qsize;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: contribution of condensation to primary rhs");

      size_type m = KiQ[0]->size(0);
      Ri.adjust_sizes(m);
      gmm::clear(Ri.as_vector());
      for (size_type k=0; k < Qsize; ++k) {
        const base_tensor &K1 = *KiQ[k], &R2 = *RQpr[k];
        size_type qqsize = K1.size(1);
        GMM_ASSERT1(K1.size(0) == m && R2.size(0) == qqsize, "Internal error");
        base_tensor::iterator it = Ri.begin();
        for (size_type ii = 0; ii < m; ++ii, ++it)
          for (size_type qq = 0; qq < qqsize; ++qq)
            *it -= K1[ii+qq*m] * R2[qq];
        GA_DEBUG_ASSERT(it == Ri.end(), "Wrong sizes");
      }
      return 0;
    }
    ga_instruction_condensation_super_R(base_tensor &Ri_,
                                        const std::vector<base_tensor *> KiQ_,
                                        const std::vector<base_tensor *> RQpr_)
      : Ri(Ri_), KiQ(KiQ_), RQpr(RQpr_)
    {
      Qsize = KiQ.size();
      GMM_ASSERT1(KiQ.size() == RQpr.size(), "Internal error");
    }
  };

  //=========================================================================
  // Compilation of assembly trees into a list of basic instructions
  //=========================================================================

  static void extend_variable_in_gis(const ga_workspace &workspace,
                                     const std::string &varname,
                                     ga_instruction_set &gis) {
    if (workspace.variable_group_exists(varname)) {
      for (const std::string &v : workspace.variable_group(varname))
        extend_variable_in_gis(workspace, v, gis);
    } else if (gis.extended_vars.count(varname) == 0) {
      const mesh_fem *mf = workspace.associated_mf(varname);
      if (mf->is_reduced()) {
        auto n = (mf->get_qdim() == 1) ? workspace.qdim(varname) : 1;
        base_vector &U = gis.really_extended_vars[varname];
        gmm::resize(U, mf->nb_basic_dof() * n);
        mf->extend_vector(workspace.value(varname), U);
        gis.extended_vars[varname] = &(gis.really_extended_vars[varname]);
      } else {
        gis.extended_vars[varname] = &(workspace.value(varname));
      }
    }
  }

  static void ga_clear_node_list
  (pga_tree_node pnode, std::map<scalar_type,
   std::list<pga_tree_node> > &node_list) {
    std::list<pga_tree_node> &loc_node_list = node_list[pnode->hash_value];
    for (std::list<pga_tree_node>::iterator it = loc_node_list.begin();
         it != loc_node_list.end(); ) {
      if (*it == pnode) it = loc_node_list.erase(it); else ++it;
    }
    for (size_type i = 0; i < pnode->children.size(); ++i)
      ga_clear_node_list(pnode->children[i], node_list);
  }

  // workspace argument  is not const because of declaration of temporary
  // unreduced variables
  static void ga_compile_node(const pga_tree_node pnode,
                              ga_workspace &workspace,
                              ga_instruction_set &gis,
                              ga_instruction_set::region_mim_instructions &rmi,
                              const mesh &m, bool function_case,
                              ga_if_hierarchy &if_hierarchy) {

    if (pnode->node_type == GA_NODE_PREDEF_FUNC ||
        pnode->node_type == GA_NODE_OPERATOR ||
        pnode->node_type == GA_NODE_SPEC_FUNC ||
        pnode->node_type == GA_NODE_CONSTANT ||
        pnode->node_type == GA_NODE_ALLINDICES ||
        pnode->node_type == GA_NODE_RESHAPE ||
        pnode->node_type == GA_NODE_SWAP_IND ||
        pnode->node_type == GA_NODE_IND_MOVE_LAST ||
        pnode->node_type == GA_NODE_CONTRACT) return;

    // cout << "compiling "; ga_print_node(pnode, cout); cout << endl;

    pga_instruction pgai;
    ga_if_hierarchy *pif_hierarchy = &if_hierarchy;
    ga_if_hierarchy new_if_hierarchy;

    const mesh_fem *mf1 = 0, *mf2 = 0;
    const mesh_fem **mfg1 = 0, **mfg2 = 0;
    fem_interpolation_context *pctx1 = 0, *pctx2 = 0;
    bool tensor_to_clear = false;
    bool tensor_to_adapt = false;

    if (pnode->test_function_type) {
      if (pnode->name_test1.size())
        mf1 = workspace.associated_mf(pnode->name_test1);
      if (mf1) {
        pctx1 = &(gis.ctx);
        const std::string &intn1 = pnode->interpolate_name_test1;
        if (intn1.size()) {
          if (workspace.secondary_domain_exists(intn1)) {
            pctx1 = &(rmi.secondary_domain_infos.ctx);
          } else {
            tensor_to_adapt = true;
            pctx1 = &(rmi.interpolate_infos[intn1].ctx);
            if (workspace.variable_group_exists(pnode->name_test1)) {
              ga_instruction_set::variable_group_info &vgi =
                rmi.interpolate_infos[intn1].groups_info[pnode->name_test1];
              mfg1 = &(vgi.mf);
              mf1 = 0;
            }
          }
        }
      }
      if (pnode->name_test2.size())
        mf2 = workspace.associated_mf(pnode->name_test2);
      if (mf2) {
        pctx2 = &(gis.ctx);
        const std::string &intn2 = pnode->interpolate_name_test2;
        if (intn2.size()) {
          if (workspace.secondary_domain_exists(intn2)) {
            pctx2 = &(rmi.secondary_domain_infos.ctx);
          } else {
            tensor_to_adapt = true;
            pctx2 = &(rmi.interpolate_infos[intn2].ctx);
            if (workspace.variable_group_exists(pnode->name_test2)) {
              ga_instruction_set::variable_group_info &vgi =
                rmi.interpolate_infos[intn2].groups_info[pnode->name_test2];
              mfg2 = &(vgi.mf);
              mf2 = 0;
            }
          }
        }
      }
    }

    // Produce a resize instruction which is stored if no equivalent node is
    // detected and if the mesh is not uniform.
    pnode->t.set_to_original(); pnode->t.set_sparsity(0, 0);
    bool is_uniform = false;
    if (pnode->test_function_type == 1) {
      if (mf1 || mfg1)
        pgai = std::make_shared<ga_instruction_first_ind_tensor>
          (pnode->tensor(), *pctx1, pnode->qdim1, mf1, mfg1);
      if (mf1 && mf1->is_uniform())
        { is_uniform = true; pctx1->invalid_convex_num(); }
    } else if (pnode->test_function_type == 2) {
      if (mf2 || mfg2)
        pgai = std::make_shared<ga_instruction_first_ind_tensor>
          (pnode->tensor(), *pctx2, pnode->qdim2, mf2, mfg2);
      if (mf2 && mf2->is_uniform())
        { is_uniform = true; pctx2->invalid_convex_num(); }
    } else if (pnode->test_function_type == 3) {
      if ((mf1 || mfg1) && (mf2 || mfg2)) {
        pgai = std::make_shared<ga_instruction_two_first_ind_tensor>
          (pnode->tensor(), *pctx1, *pctx2, pnode->qdim1, mf1, mfg1,
           pnode->qdim2, mf2, mfg2);
        if (mf1 && mf1->is_uniform() && mf2 && mf2->is_uniform()) {
          is_uniform = true;
          pctx1->invalid_convex_num();
          pctx2->invalid_convex_num();
        }
      } else if (mf1 || mfg1) {
        pgai = std::make_shared<ga_instruction_first_ind_tensor>
          (pnode->tensor(), *pctx1, pnode->qdim1, mf1, mfg1);
        if (mf1 && mf1->is_uniform())
          { is_uniform = true; pctx1->invalid_convex_num(); }
      } else if (mf2 || mfg2) {
        pgai = std::make_shared<ga_instruction_second_ind_tensor>
          (pnode->tensor(), *pctx2, pnode->qdim2, mf2, mfg2);
        if (mf2 && mf2->is_uniform())
          { is_uniform = true; pctx2->invalid_convex_num(); }
      }
    }

    // Optimization: detects if an equivalent node has already been compiled
    pnode->t.set_to_original();
    if (rmi.node_list.count(pnode->hash_value) != 0) {
      for (pga_tree_node &pnode1 : rmi.node_list[pnode->hash_value]) {
        // cout << "found potential equivalent nodes ";
        // ga_print_node(pnode, cout);
        // cout << " and "; ga_print_node(pnode1, cout); cout << endl;
        if (sub_tree_are_equal(pnode, pnode1, workspace, 1)) {
          pnode->t.set_to_copy(pnode1->t);
          return;
        }
        if (sub_tree_are_equal(pnode, pnode1, workspace, 2)) {
          // cout << "confirmed with transpose" << endl;
          if (pnode->nb_test_functions() == 2) {
            if (pgai) { // resize instruction if needed
              if (is_uniform)
                { pgai->exec(); }
              else { rmi.instructions.push_back(std::move(pgai)); }
            }
            pgai = std::make_shared<ga_instruction_transpose_test>
              (pnode->tensor(), pnode1->tensor());
            rmi.instructions.push_back(std::move(pgai));
          } else {
            pnode->t.set_to_copy(pnode1->t);
          }
          return;
        }
        // cout << "sub_tree_are_equal = " << int(sub_tree_are_equal(pnode, pnode1, workspace, 1)) << endl;
        std::stringstream ss;
        ss << "Detected wrong equivalent nodes:" << endl;
        ga_print_node(pnode, ss);
        ss << endl << " and " << endl;
        ga_print_node(pnode1, ss);
        ss << endl << "No problem, but hash values could be adapted." << endl;
        GMM_TRACE2(ss.str());
      }
    }

    if (pgai) { // resize instruction if needed and no equivalent node detected
      if (is_uniform) { pgai->exec(); }
      else {
        if (tensor_to_adapt)
          rmi.instructions.push_back(std::move(pgai));
        else
          rmi.elt_instructions.push_back(std::move(pgai));
      }
    }

    size_type interpolate_filter_inst = rmi.instructions.size();
    if (pnode->node_type == GA_NODE_INTERPOLATE_FILTER) {
      pgai = pga_instruction();
      rmi.instructions.push_back(std::move(pgai));
      if_hierarchy.increment();
      new_if_hierarchy.child_of(if_hierarchy);
      pif_hierarchy = &new_if_hierarchy;
    }

    for (size_type i = 0; i < pnode->children.size(); ++i)
      ga_compile_node(pnode->children[i], workspace, gis, rmi, m,
                      function_case, *pif_hierarchy);

    if (pnode->node_type == GA_NODE_INTERPOLATE_FILTER) {
      const std::string &intn = pnode->interpolate_name;
      ga_instruction_set::interpolate_info &inin = rmi.interpolate_infos[intn];
      pgai = std::make_shared<ga_instruction_interpolate_filter>
        (pnode->tensor(), inin, pnode->nbc1,
         int(rmi.instructions.size() - interpolate_filter_inst));
      rmi.instructions[interpolate_filter_inst].swap(pgai);
      pgai = std::make_shared<ga_instruction_copy_tensor>
        (pnode->tensor(), pnode->children[0]->tensor());
      rmi.instructions.push_back(std::move(pgai));
      ga_clear_node_list(pnode->children[0], rmi.node_list);
    }

    static scalar_type minus = -scalar_type(1);
    size_type nbch = pnode->children.size();
    pga_tree_node child0 = (nbch > 0) ? pnode->children[0] : 0;
    pga_tree_node child1 = (nbch > 1) ? pnode->children[1] : 0;
    bgeot::multi_index mi;
    const bgeot::multi_index &size0 = child0 ? child0->t.sizes() : mi;
    // const bgeot::multi_index &size1 = child1 ? child1->t.sizes() : mi;
    size_type dim0 = child0 ? child0->tensor_order() : 0;
    size_type dim1 = child1 ? child1->tensor_order() : 0;

    switch (pnode->node_type) {

    case GA_NODE_PREDEF_FUNC: case GA_NODE_OPERATOR: case GA_NODE_SPEC_FUNC:
    case GA_NODE_CONSTANT: case GA_NODE_ALLINDICES: case GA_NODE_ZERO:
    case GA_NODE_RESHAPE:  case GA_NODE_CROSS_PRODUCT:
    case GA_NODE_SWAP_IND: case GA_NODE_IND_MOVE_LAST:
    case GA_NODE_CONTRACT: case GA_NODE_INTERPOLATE_FILTER:
      break;

    case GA_NODE_X:
      GMM_ASSERT1(!function_case,
                  "No use of X is allowed in scalar functions");
      if (pnode->nbc1) {
        GA_DEBUG_ASSERT(pnode->tensor().size() == 1, "dimensions mismatch");
        GMM_ASSERT1(pnode->nbc1 <= m.dim(),
                    "Bad index for X in expression");
        pgai = std::make_shared<ga_instruction_X_component>
            (pnode->tensor()[0], gis.ctx, pnode->nbc1-1);
      } else {
        if (pnode->tensor().size() != m.dim())
          pnode->init_vector_tensor(m.dim());
        pgai = std::make_shared<ga_instruction_X>(pnode->tensor(), gis.ctx);
      }
      rmi.instructions.push_back(std::move(pgai));
      break;

    case GA_NODE_ELT_SIZE:
      GMM_ASSERT1(!function_case,
                  "No use of element_size is allowed in functions");
      if (pnode->tensor().size() != 1) pnode->init_scalar_tensor(0);
      pgai = std::make_shared<ga_instruction_element_size>
        (pnode->tensor(), gis.elt_size);
      gis.need_elt_size = true;
      rmi.instructions.push_back(std::move(pgai));
      break;

    case GA_NODE_ELT_K:
      GMM_ASSERT1(!function_case,
                  "No use of element_K is allowed in functions");
      pgai = std::make_shared<ga_instruction_element_K>(pnode->tensor(),
                                                        gis.ctx);
      rmi.instructions.push_back(std::move(pgai));
      break;

    case GA_NODE_ELT_B:
      GMM_ASSERT1(!function_case,
                  "No use of element_B is allowed in functions");
      pgai = std::make_shared<ga_instruction_element_B>(pnode->tensor(),
                                                        gis.ctx);
      rmi.instructions.push_back(std::move(pgai));
      break;

    case GA_NODE_NORMAL:
      {
        GMM_ASSERT1(!function_case,
                    "No use of Normal is allowed in functions");
        if (pnode->tensor().size() != m.dim())
          pnode->init_vector_tensor(m.dim());
        const mesh_im_level_set *mimls
          = dynamic_cast<const mesh_im_level_set *>(rmi.im);
        if (mimls && mimls->location()==mesh_im_level_set::INTEGRATE_BOUNDARY) {
          // Appel avec ctx (pt de Gauss)
          pgai = std::make_shared<ga_instruction_level_set_normal_vector>
            (pnode->tensor(), mimls, gis.ctx);
          rmi.instructions.push_back(std::move(pgai));
        } else {
          pgai = std::make_shared<ga_instruction_copy_Normal>
            (pnode->tensor(), gis.Normal);
          rmi.instructions.push_back(std::move(pgai));
        }
      }
      break;

    case GA_NODE_INTERPOLATE_X:
    case GA_NODE_INTERPOLATE_NORMAL:
      GMM_ASSERT1(!function_case,
                  "No use of Interpolate is allowed in functions");
      if (pnode->tensor().size() != m.dim())
        pnode->init_vector_tensor(m.dim());
      if (pnode->node_type == GA_NODE_INTERPOLATE_X)
        pgai = std::make_shared<ga_instruction_copy_interpolated_small_vect>
               (pnode->tensor(),
                rmi.interpolate_infos[pnode->interpolate_name].pt_y,
                rmi.interpolate_infos[pnode->interpolate_name]);
      else if (pnode->node_type == GA_NODE_INTERPOLATE_NORMAL)
        pgai = std::make_shared<ga_instruction_copy_Normal>
               (pnode->tensor(),
                rmi.interpolate_infos[pnode->interpolate_name].Normal);
      rmi.instructions.push_back(std::move(pgai));
      break;

    case GA_NODE_INTERPOLATE_ELT_K:
    case GA_NODE_INTERPOLATE_ELT_B:
      GMM_ASSERT1(!function_case,
                  "No use of Interpolate is allowed in functions");
      if (pnode->node_type == GA_NODE_INTERPOLATE_ELT_K)
        pgai = std::make_shared<ga_instruction_element_K>
               (pnode->tensor(),
                rmi.interpolate_infos[pnode->interpolate_name].ctx);
      else if (pnode->node_type == GA_NODE_INTERPOLATE_ELT_B)
        pgai = std::make_shared<ga_instruction_element_B>
               (pnode->tensor(),
                rmi.interpolate_infos[pnode->interpolate_name].ctx);
      rmi.instructions.push_back(std::move(pgai));
      break;

    case GA_NODE_SECONDARY_DOMAIN_X:
    case GA_NODE_SECONDARY_DOMAIN_NORMAL:
      {
        GMM_ASSERT1(!function_case,
                    "No use of Secondary_domain is allowed in functions");
        auto psd = workspace.secondary_domain(pnode->interpolate_name);
        size_type sddim = psd->mim().linked_mesh().dim();
        if (pnode->tensor().size() != sddim)
          pnode->init_vector_tensor(sddim);
        if (pnode->node_type == GA_NODE_SECONDARY_DOMAIN_X)
          pgai = std::make_shared<ga_instruction_X>
            (pnode->tensor(), rmi.secondary_domain_infos.ctx);
        else if (pnode->node_type == GA_NODE_SECONDARY_DOMAIN_NORMAL)
          pgai = std::make_shared<ga_instruction_copy_Normal>
            (pnode->tensor(), rmi.secondary_domain_infos.Normal);
        rmi.instructions.push_back(std::move(pgai));
      }
      break;

    case GA_NODE_VAL: case GA_NODE_GRAD:
    case GA_NODE_HESS: case GA_NODE_DIVERG:
    case GA_NODE_ELEMENTARY_VAL: case GA_NODE_ELEMENTARY_GRAD:
    case GA_NODE_ELEMENTARY_HESS: case GA_NODE_ELEMENTARY_DIVERG:
    case GA_NODE_XFEM_PLUS_VAL: case GA_NODE_XFEM_PLUS_GRAD:
    case GA_NODE_XFEM_PLUS_HESS: case GA_NODE_XFEM_PLUS_DIVERG:
    case GA_NODE_XFEM_MINUS_VAL: case GA_NODE_XFEM_MINUS_GRAD:
    case GA_NODE_XFEM_MINUS_HESS: case GA_NODE_XFEM_MINUS_DIVERG:
      {
        bool is_elementary = (pnode->node_type == GA_NODE_ELEMENTARY_VAL ||
                              pnode->node_type == GA_NODE_ELEMENTARY_GRAD ||
                              pnode->node_type == GA_NODE_ELEMENTARY_HESS ||
                              pnode->node_type == GA_NODE_ELEMENTARY_DIVERG);
        if (function_case) {
          GMM_ASSERT1(!is_elementary,
                      "No elementary transformation is allowed in functions");
          GMM_ASSERT1(pnode->node_type != GA_NODE_XFEM_PLUS_VAL &&
                      pnode->node_type != GA_NODE_XFEM_PLUS_GRAD &&
                      pnode->node_type != GA_NODE_XFEM_PLUS_HESS &&
                      pnode->node_type != GA_NODE_XFEM_PLUS_DIVERG,
                      "Xfem_plus not allowed in functions");
          GMM_ASSERT1(pnode->node_type != GA_NODE_XFEM_MINUS_VAL &&
                      pnode->node_type != GA_NODE_XFEM_MINUS_GRAD &&
                      pnode->node_type != GA_NODE_XFEM_MINUS_HESS &&
                      pnode->node_type != GA_NODE_XFEM_MINUS_DIVERG,
                      "Xfem_plus not allowed in functions");
          const mesh_fem *mf = workspace.associated_mf(pnode->name);
          const im_data *imd = workspace.associated_im_data(pnode->name);
          GMM_ASSERT1(!mf, "No fem expression is allowed in "
                      "function expression");
          GMM_ASSERT1(!imd, "No integration method data is allowed in "
                      "function expression");
          if (gmm::vect_size(workspace.value(pnode->name)) == 1)
            pgai = std::make_shared<ga_instruction_copy_scalar>
              (pnode->tensor()[0], (workspace.value(pnode->name))[0]);
          else
            pgai = std::make_shared<ga_instruction_copy_vect>
              (pnode->tensor().as_vector(), workspace.value(pnode->name));
          rmi.instructions.push_back(std::move(pgai));
        } else {
          const mesh_fem *mf = workspace.associated_mf(pnode->name), *mfo=mf;
          const im_data *imd = workspace.associated_im_data(pnode->name);
          
          if (is_elementary) {
            mf = workspace.associated_mf(pnode->elementary_target);
            GMM_ASSERT1(mf && mfo,
                        "Wrong context for elementary transformation");
            GMM_ASSERT1(&(mfo->linked_mesh()) == &(m),
                        "The finite element of variable " << pnode->name
                        << " has to be defined on the same mesh as the "
                        << "integration method or interpolation used");
          }
          
          if (imd) {
            GMM_ASSERT1(pnode->node_type == GA_NODE_VAL,
                        "Only values can be extracted on im_data (no " <<
                        "gradient, Hessian, xfem or elementary tranformation" <<
                        " allowed)");
            pgai = std::make_shared<ga_instruction_extract_local_im_data>
              (pnode->tensor(), *imd, workspace.value(pnode->name),
               gis.pai, gis.ctx, workspace.qdim(pnode->name));
            rmi.instructions.push_back(std::move(pgai));
          } else {
            GMM_ASSERT1(mf, "Internal error");
            
            GMM_ASSERT1(&(mf->linked_mesh()) == &(m),
                        "The finite element of variable " <<
                        (is_elementary ? pnode->elementary_target : pnode->name)
                        << " has to be defined on the same mesh as the "
                        << "integration method or interpolation used");

            // An instruction for extracting local dofs of the variable.
            if (rmi.local_dofs.count(pnode->name) == 0) {
              rmi.local_dofs[pnode->name] = base_vector(1);
              extend_variable_in_gis(workspace, pnode->name, gis);
              // cout << "local dof of " << pnode->name << endl;
              size_type qmult2 = mfo->get_qdim();
              if (qmult2 > 1 && !(mfo->is_uniformly_vectorized()))
                qmult2 = size_type(-1);
              pgai = std::make_shared<ga_instruction_slice_local_dofs>
                (*mfo, *(gis.extended_vars[pnode->name]), gis.ctx,
                 rmi.local_dofs[pnode->name],
                 workspace.qdim(pnode->name) / mfo->get_qdim(), qmult2);
              rmi.elt_instructions.push_back(std::move(pgai));
            }
            
            // An instruction for pfp update
            if (mf->is_uniform()) {
              if (rmi.pfps.count(mf) == 0) {
                rmi.pfps[mf] = 0;
                pgai = std::make_shared<ga_instruction_update_pfp>
                  (*mf, rmi.pfps[mf], gis.ctx, gis.fp_pool);
                rmi.begin_instructions.push_back(std::move(pgai));
              }
            } else if (rmi.pfps.count(mf) == 0 ||
                       !if_hierarchy.is_compatible(rmi.pfp_hierarchy[mf])) {
              rmi.pfp_hierarchy[mf].push_back(if_hierarchy);
              rmi.pfps[mf] = 0;
              pgai = std::make_shared<ga_instruction_update_pfp>
                (*mf, rmi.pfps[mf], gis.ctx, gis.fp_pool);
              rmi.instructions.push_back(std::move(pgai));
            }

            // An instruction for the base value
            pgai = pga_instruction();
            switch (pnode->node_type) {
            case GA_NODE_VAL: case GA_NODE_ELEMENTARY_VAL:
              if (rmi.base.count(mf) == 0 ||
                  !if_hierarchy.is_compatible(rmi.base_hierarchy[mf])) {
                rmi.base_hierarchy[mf].push_back(if_hierarchy);
                pgai = std::make_shared<ga_instruction_val_base>
                  (rmi.base[mf], gis.ctx, *mf, rmi.pfps[mf]);
              }
              break;
            case GA_NODE_XFEM_PLUS_VAL:
              if (rmi.xfem_plus_base.count(mf) == 0 ||
                  !if_hierarchy.is_compatible(rmi.xfem_plus_base_hierarchy[mf]))
                {
                  rmi.xfem_plus_base_hierarchy[mf].push_back(if_hierarchy);
                  pgai = std::make_shared<ga_instruction_xfem_plus_val_base>
                    (rmi.xfem_plus_base[mf], gis.ctx, *mf, rmi.pfps[mf]);
                }
              break;
            case GA_NODE_XFEM_MINUS_VAL:
              if (rmi.xfem_minus_base.count(mf) == 0 ||
                  !if_hierarchy.is_compatible(rmi.xfem_minus_base_hierarchy[mf]))
                {
                  rmi.xfem_minus_base_hierarchy[mf].push_back(if_hierarchy);
                  pgai = std::make_shared<ga_instruction_xfem_minus_val_base>
                    (rmi.xfem_minus_base[mf], gis.ctx, *mf, rmi.pfps[mf]);
                }
              break;
            case GA_NODE_GRAD: case GA_NODE_DIVERG:
            case GA_NODE_ELEMENTARY_GRAD: case GA_NODE_ELEMENTARY_DIVERG:
              if (rmi.grad.count(mf) == 0 ||
                  !if_hierarchy.is_compatible(rmi.grad_hierarchy[mf])) {
                rmi.grad_hierarchy[mf].push_back(if_hierarchy);
                pgai = std::make_shared<ga_instruction_grad_base>
                  (rmi.grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
              }
              break;
            case GA_NODE_XFEM_PLUS_GRAD: case GA_NODE_XFEM_PLUS_DIVERG:
              if (rmi.xfem_plus_grad.count(mf) == 0 ||
                  !if_hierarchy.is_compatible(rmi.xfem_plus_grad_hierarchy[mf]))
                {
                  rmi.xfem_plus_grad_hierarchy[mf].push_back(if_hierarchy);
                  pgai = std::make_shared<ga_instruction_xfem_plus_grad_base>
                    (rmi.xfem_plus_grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
                }
              break;
            case GA_NODE_XFEM_MINUS_GRAD: case GA_NODE_XFEM_MINUS_DIVERG:
              if (rmi.xfem_minus_grad.count(mf) == 0 ||
                  !if_hierarchy.is_compatible(rmi.xfem_minus_grad_hierarchy[mf]))
                {
                  rmi.xfem_minus_grad_hierarchy[mf].push_back(if_hierarchy);
                  pgai = std::make_shared<ga_instruction_xfem_minus_grad_base>
                    (rmi.xfem_minus_grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
                }
              break;
            case GA_NODE_HESS: case GA_NODE_ELEMENTARY_HESS:
              if (rmi.hess.count(mf) == 0 ||
                  !if_hierarchy.is_compatible(rmi.hess_hierarchy[mf])) {
                rmi.hess_hierarchy[mf].push_back(if_hierarchy);
                pgai = std::make_shared<ga_instruction_hess_base>
                  (rmi.hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
              }
              break;
            case GA_NODE_XFEM_PLUS_HESS:
              if (rmi.xfem_plus_hess.count(mf) == 0 ||
                  !if_hierarchy.is_compatible(rmi.xfem_plus_hess_hierarchy[mf]))
                {
                  rmi.xfem_plus_hess_hierarchy[mf].push_back(if_hierarchy);
                  pgai = std::make_shared<ga_instruction_xfem_plus_hess_base>
                    (rmi.xfem_plus_hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
                }
              break;
            case GA_NODE_XFEM_MINUS_HESS:
              if (rmi.xfem_minus_hess.count(mf) == 0 ||
                  !if_hierarchy.is_compatible(rmi.xfem_minus_hess_hierarchy[mf]))
                {
                  rmi.xfem_minus_hess_hierarchy[mf].push_back(if_hierarchy);
                  pgai = std::make_shared<ga_instruction_xfem_minus_hess_base>
                    (rmi.xfem_minus_hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
                }
              break;

            default : GMM_ASSERT1(false, "Internal error");
            }
            if (pgai) rmi.instructions.push_back(std::move(pgai));
            
            // The eval instruction
            switch (pnode->node_type) {
            case GA_NODE_VAL: // --> t(target_dim*Qmult)
              pgai = std::make_shared<ga_instruction_val>
                (pnode->tensor(), rmi.base[mf], rmi.local_dofs[pnode->name],
                 workspace.qdim(pnode->name));
              break;
            case GA_NODE_GRAD: // --> t(target_dim*Qmult,N)
              pgai = std::make_shared<ga_instruction_grad>
                (pnode->tensor(), rmi.grad[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
              break;
            case GA_NODE_HESS: // --> t(target_dim*Qmult,N,N)
              pgai = std::make_shared<ga_instruction_hess>
                (pnode->tensor(), rmi.hess[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
              break;
            case GA_NODE_DIVERG: // --> t(1)
              pgai = std::make_shared<ga_instruction_diverg>
                (pnode->tensor(), rmi.grad[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
              break;
            case GA_NODE_XFEM_PLUS_VAL: // --> t(target_dim*Qmult)
              pgai = std::make_shared<ga_instruction_val>
                (pnode->tensor(), rmi.xfem_plus_base[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
              break;
            case GA_NODE_XFEM_PLUS_GRAD: // --> t(target_dim*Qmult,N)
              pgai = std::make_shared<ga_instruction_grad>
                (pnode->tensor(), rmi.xfem_plus_grad[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
              break;
            case GA_NODE_XFEM_PLUS_HESS: // --> t(target_dim*Qmult,N,N)
              pgai = std::make_shared<ga_instruction_hess>
                (pnode->tensor(), rmi.xfem_plus_hess[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
              break;
            case GA_NODE_XFEM_PLUS_DIVERG: // --> t(1)
              pgai = std::make_shared<ga_instruction_diverg>
                (pnode->tensor(), rmi.xfem_plus_grad[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
              break;
            case GA_NODE_XFEM_MINUS_VAL: // --> t(target_dim*Qmult)
              pgai = std::make_shared<ga_instruction_val>
                (pnode->tensor(), rmi.xfem_minus_base[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
              break;
            case GA_NODE_XFEM_MINUS_GRAD: // --> t(target_dim*Qmult,N)
              pgai = std::make_shared<ga_instruction_grad>
                (pnode->tensor(), rmi.xfem_minus_grad[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
              break;
            case GA_NODE_XFEM_MINUS_HESS: // --> t(target_dim*Qmult,N,N)
              pgai = std::make_shared<ga_instruction_hess>
                (pnode->tensor(), rmi.xfem_minus_hess[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
              break;
            case GA_NODE_XFEM_MINUS_DIVERG: // --> t(1)
              pgai = std::make_shared<ga_instruction_diverg>
                (pnode->tensor(), rmi.xfem_minus_grad[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
              break;
            case GA_NODE_ELEMENTARY_VAL:
              { // --> t(target_dim*Qmult)
                ga_instruction_set::elementary_trans_info &eti
                  = rmi.elementary_trans_infos
                  [std::make_tuple(pnode->elementary_name, mfo, mf)];
                pgai =
                  std::make_shared<ga_instruction_elementary_trans_val>
                  (pnode->tensor(), rmi.base[mf],
                   rmi.local_dofs[pnode->name],
                   workspace.qdim(pnode->elementary_target),
                   workspace.elementary_transformation(pnode->elementary_name),
                   *mfo, *mf, gis.ctx, eti.M, eti.icv);
              }
              break;
            case GA_NODE_ELEMENTARY_GRAD:
              { // --> t(target_dim*Qmult,N)
                ga_instruction_set::elementary_trans_info &eti
                  = rmi.elementary_trans_infos
                  [std::make_tuple(pnode->elementary_name, mfo, mf)];
                pgai =
                  std::make_shared<ga_instruction_elementary_trans_grad>
                  (pnode->tensor(), rmi.grad[mf],
                   rmi.local_dofs[pnode->name],
                   workspace.qdim(pnode->elementary_target),
                   workspace.elementary_transformation(pnode->elementary_name),
                   *mfo, *mf, gis.ctx, eti.M, eti.icv);
              }
              break;
            case GA_NODE_ELEMENTARY_HESS:
              { // --> t(target_dim*Qmult,N,N)
                ga_instruction_set::elementary_trans_info &eti
                  = rmi.elementary_trans_infos
                  [std::make_tuple(pnode->elementary_name, mfo, mf)];
                pgai =
                  std::make_shared<ga_instruction_elementary_trans_hess>
                  (pnode->tensor(), rmi.hess[mf],
                   rmi.local_dofs[pnode->name],
                   workspace.qdim(pnode->elementary_target),
                   workspace.elementary_transformation(pnode->elementary_name),
                   *mfo, *mf, gis.ctx, eti.M, eti.icv);
              }
              break;
            case GA_NODE_ELEMENTARY_DIVERG:
              { // --> t(1)
                ga_instruction_set::elementary_trans_info &eti
                  = rmi.elementary_trans_infos
                  [std::make_tuple(pnode->elementary_name, mfo, mf)];
                pgai =
                  std::make_shared<ga_instruction_elementary_trans_diverg>
                  (pnode->tensor(), rmi.grad[mf],
                   rmi.local_dofs[pnode->name],
                   workspace.qdim(pnode->elementary_target),
                   workspace.elementary_transformation(pnode->elementary_name),
                   *mfo, *mf, gis.ctx, eti.M, eti.icv);
              }
              break;
            default: break;
            }
            rmi.instructions.push_back(std::move(pgai));
          }
        }
      }
      break;

    case GA_NODE_SECONDARY_DOMAIN_VAL: case GA_NODE_SECONDARY_DOMAIN_GRAD:
    case GA_NODE_SECONDARY_DOMAIN_HESS: case GA_NODE_SECONDARY_DOMAIN_DIVERG:
      {
        GMM_ASSERT1(!function_case, "internal error");
        const mesh_fem *mf = workspace.associated_mf(pnode->name);
        const im_data *imd = workspace.associated_im_data(pnode->name);
        const std::string &intn = pnode->interpolate_name;
        auto &sdi = rmi.secondary_domain_infos;

        fem_interpolation_context *pctx = &(sdi.ctx);
        papprox_integration pai = sdi.pai;
        psecondary_domain psd = workspace.secondary_domain(intn);

        if (imd) {
          pgai = std::make_shared<ga_instruction_extract_local_im_data>
            (pnode->tensor(), *imd, workspace.value(pnode->name),
             pai, *pctx, workspace.qdim(pnode->name));
          rmi.instructions.push_back(std::move(pgai));
        } else {
          GMM_ASSERT1(mf, "Internal error");
          GMM_ASSERT1(&(mf->linked_mesh()) == &(psd->mim().linked_mesh()),
                      "The finite element of variable " << pnode->name <<
                      " has to be defined on the same mesh as the "
                      "integration method or interpolation used on the "
                      "secondary domain");

          // An instruction for extracting local dofs of the variable.
          if (sdi.local_dofs.count(pnode->name) == 0) {
            sdi.local_dofs[pnode->name] = base_vector(1);
            extend_variable_in_gis(workspace, pnode->name, gis);
            size_type qmult2 = mf->get_qdim();
            if (qmult2 > 1 && !(mf->is_uniformly_vectorized()))
              qmult2 = size_type(-1);
            pgai = std::make_shared<ga_instruction_slice_local_dofs>
              (*mf, *(gis.extended_vars[pnode->name]), *pctx,
               sdi.local_dofs[pnode->name],
               workspace.qdim(pnode->name) / mf->get_qdim(), qmult2);
            rmi.elt_instructions.push_back(std::move(pgai));
          }

          // An instruction for pfp update
          if (mf->is_uniform()) {
            if (sdi.pfps.count(mf) == 0) {
              sdi.pfps[mf] = 0;
              pgai = std::make_shared<ga_instruction_update_pfp>
                (*mf, sdi.pfps[mf], *pctx, gis.fp_pool);
              rmi.begin_instructions.push_back(std::move(pgai));
            }
          } else if (sdi.pfps.count(mf) == 0 ||
                     !if_hierarchy.is_compatible(rmi.pfp_hierarchy[mf])) {
            rmi.pfp_hierarchy[mf].push_back(if_hierarchy);
            sdi.pfps[mf] = 0;
            pgai = std::make_shared<ga_instruction_update_pfp>
              (*mf, sdi.pfps[mf], *pctx, gis.fp_pool);
            rmi.instructions.push_back(std::move(pgai));
          }

          // An instruction for the base value
          pgai = pga_instruction();
          switch (pnode->node_type) {
          case GA_NODE_SECONDARY_DOMAIN_VAL:
            if (sdi.base.count(mf) == 0 ||
               !(if_hierarchy.is_compatible(rmi.base_hierarchy[mf]))) {
              rmi.base_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_val_base>
                (sdi.base[mf], *pctx, *mf, sdi.pfps[mf]);
            }
            break;
          case GA_NODE_SECONDARY_DOMAIN_GRAD:
          case GA_NODE_SECONDARY_DOMAIN_DIVERG:
            if (sdi.grad.count(mf) == 0 ||
                !(if_hierarchy.is_compatible(rmi.grad_hierarchy[mf]))) {
              rmi.grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_grad_base>
                (sdi.grad[mf], *pctx, *mf, sdi.pfps[mf]);
            }
            break;
          case GA_NODE_SECONDARY_DOMAIN_HESS:
            if (sdi.hess.count(mf) == 0 ||
                !(if_hierarchy.is_compatible(rmi.hess_hierarchy[mf]))) {
              rmi.hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_hess_base>
                (sdi.hess[mf], *pctx, *mf, sdi.pfps[mf]);
            }
            break;
          default : GMM_ASSERT1(false, "Internal error");
          }
          if (pgai) rmi.instructions.push_back(std::move(pgai));

          // The eval instruction
          switch (pnode->node_type) {
          case GA_NODE_SECONDARY_DOMAIN_VAL: // --> t(target_dim*Qmult)
            pgai = std::make_shared<ga_instruction_val>
              (pnode->tensor(), sdi.base[mf], sdi.local_dofs[pnode->name],
               workspace.qdim(pnode->name));
            break;
          case GA_NODE_SECONDARY_DOMAIN_GRAD: // --> t(target_dim*Qmult,N)
            pgai = std::make_shared<ga_instruction_grad>
              (pnode->tensor(), sdi.grad[mf],
               sdi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          case GA_NODE_SECONDARY_DOMAIN_HESS: // --> t(target_dim*Qmult,N,N)
            pgai = std::make_shared<ga_instruction_hess>
              (pnode->tensor(), sdi.hess[mf],
               sdi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          case GA_NODE_SECONDARY_DOMAIN_DIVERG: // --> t(1)
            pgai = std::make_shared<ga_instruction_diverg>
              (pnode->tensor(), sdi.grad[mf],
               sdi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          default: break;
          }
          rmi.instructions.push_back(std::move(pgai));
        }
      }
      break;

    case GA_NODE_INTERPOLATE_VAL: case GA_NODE_INTERPOLATE_GRAD:
    case GA_NODE_INTERPOLATE_HESS: case GA_NODE_INTERPOLATE_DIVERG:
      {
        extend_variable_in_gis(workspace, pnode->name, gis);

        const mesh_fem *mfn = workspace.associated_mf(pnode->name), **mfg = 0;
        const std::string &intn = pnode->interpolate_name;
        const base_vector *Un = gis.extended_vars[pnode->name], **Ug = 0;
        fem_interpolation_context *pctx = &(rmi.interpolate_infos[intn].ctx);
        const mesh **m2 = &(rmi.interpolate_infos[intn].m);
        if (workspace.variable_group_exists(pnode->name)) {
          ga_instruction_set::variable_group_info &vgi =
            rmi.interpolate_infos[intn].groups_info[pnode->name];
          mfg = &(vgi.mf); mfn = 0; Ug = &(vgi.U); Un = 0;
        }

        if (pnode->node_type == GA_NODE_INTERPOLATE_VAL) {
          // --> t(target_dim*Qmult)
          pgai = std::make_shared<ga_instruction_interpolate_val>
            (pnode->tensor(), m2, mfn, mfg, Un, Ug, *pctx,
             workspace.qdim(pnode->name),
             gis.ipt, gis.fp_pool, rmi.interpolate_infos[intn]);
        } else if (pnode->node_type == GA_NODE_INTERPOLATE_GRAD) {
          // --> t(target_dim*Qmult,N)
          pgai = std::make_shared<ga_instruction_interpolate_grad>
            (pnode->tensor(), m2, mfn, mfg, Un, Ug, *pctx,
             workspace.qdim(pnode->name),
             gis.ipt, gis.fp_pool, rmi.interpolate_infos[intn]);
        } else if (pnode->node_type == GA_NODE_INTERPOLATE_HESS) {
          // --> t(target_dim*Qmult,N,N)
          pgai = std::make_shared<ga_instruction_interpolate_hess>
            (pnode->tensor(), m2, mfn, mfg, Un, Ug, *pctx,
             workspace.qdim(pnode->name),
             gis.ipt, gis.fp_pool, rmi.interpolate_infos[intn]);
        } else { // --> t(1)
          pgai = std::make_shared<ga_instruction_interpolate_diverg>
            (pnode->tensor(), m2, mfn, mfg, Un, Ug, *pctx,
             workspace.qdim(pnode->name),
             gis.ipt, gis.fp_pool, rmi.interpolate_infos[intn]);
        }
        rmi.instructions.push_back(std::move(pgai));
      }
      break;

    case GA_NODE_INTERPOLATE_DERIVATIVE:
      GMM_ASSERT1(!function_case,
                  "No use of Interpolate is allowed in functions");
      pgai = std::make_shared<ga_instruction_copy_tensor_possibly_void>
        (pnode->tensor(),
         rmi.interpolate_infos[pnode->interpolate_name_der]
         .derivatives[var_trans_pair(pnode->name, pnode->interpolate_name)]);
      rmi.instructions.push_back(std::move(pgai));
      break;

    case GA_NODE_VAL_TEST: case GA_NODE_GRAD_TEST:
    case GA_NODE_HESS_TEST: case GA_NODE_DIVERG_TEST:
    case GA_NODE_ELEMENTARY_VAL_TEST: case GA_NODE_ELEMENTARY_GRAD_TEST:
    case GA_NODE_ELEMENTARY_HESS_TEST: case GA_NODE_ELEMENTARY_DIVERG_TEST:
    case GA_NODE_XFEM_PLUS_VAL_TEST: case GA_NODE_XFEM_PLUS_GRAD_TEST:
    case GA_NODE_XFEM_PLUS_HESS_TEST: case GA_NODE_XFEM_PLUS_DIVERG_TEST:
    case GA_NODE_XFEM_MINUS_VAL_TEST: case GA_NODE_XFEM_MINUS_GRAD_TEST:
    case GA_NODE_XFEM_MINUS_HESS_TEST: case GA_NODE_XFEM_MINUS_DIVERG_TEST:
      // GMM_ASSERT1(!function_case,
      //            "Test functions not allowed in functions");
      {
        bool is_elementary = (pnode->node_type==GA_NODE_ELEMENTARY_VAL_TEST ||
                              pnode->node_type==GA_NODE_ELEMENTARY_GRAD_TEST ||
                              pnode->node_type==GA_NODE_ELEMENTARY_HESS_TEST ||
                              pnode->node_type==GA_NODE_ELEMENTARY_DIVERG_TEST);
        const mesh_fem *mf = workspace.associated_mf(pnode->name), *mfo=mf;
        if (is_elementary) {
          mf = workspace.associated_mf(pnode->elementary_target);
          GMM_ASSERT1(mf && mfo,
                      "Wrong context for elementary transformation");
          GMM_ASSERT1(&(mfo->linked_mesh()) == &(m),
                      "The finite element of variable " << pnode->name
                      << " has to be defined on the same mesh as the "
                      << "integration method or interpolation used");
        }
        
        if (mf) {
          GMM_ASSERT1(&(mf->linked_mesh()) == &(m),
                      "The finite element of variable " <<
                      (is_elementary ? pnode->elementary_target : pnode->name)
                      << " and the applied integration method have to be"
                      << " defined on the same mesh");

          // An instruction for pfp update
          if (is_uniform) {
            if (rmi.pfps.count(mf) == 0) {
              rmi.pfps[mf] = 0;
              pgai = std::make_shared<ga_instruction_update_pfp>
                (*mf, rmi.pfps[mf], gis.ctx, gis.fp_pool);
              rmi.begin_instructions.push_back(std::move(pgai));
            }
          } else if (rmi.pfps.count(mf) == 0 ||
                     !if_hierarchy.is_compatible(rmi.pfp_hierarchy[mf])) {
            rmi.pfp_hierarchy[mf].push_back(if_hierarchy);
            rmi.pfps[mf] = 0;
            pgai = std::make_shared<ga_instruction_update_pfp>
              (*mf, rmi.pfps[mf], gis.ctx, gis.fp_pool);
            rmi.instructions.push_back(std::move(pgai));
          }

          // An instruction for the base value
          pgai = pga_instruction();
          switch (pnode->node_type) {
          case GA_NODE_VAL_TEST: case GA_NODE_ELEMENTARY_VAL_TEST:
             if (rmi.base.count(mf) == 0 ||
                 !if_hierarchy.is_compatible(rmi.base_hierarchy[mf])) {
              rmi.base_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_val_base>
                (rmi.base[mf], gis.ctx, *mf, rmi.pfps[mf]);
             }
             break;
          case GA_NODE_XFEM_PLUS_VAL_TEST:
            if (rmi.xfem_plus_base.count(mf) == 0 ||
                !if_hierarchy.is_compatible(rmi.xfem_plus_base_hierarchy[mf]))
            {
              rmi.xfem_plus_base_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_plus_val_base>
                (rmi.xfem_plus_base[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_VAL_TEST:
            if (rmi.xfem_minus_base.count(mf) == 0 ||
                !if_hierarchy.is_compatible(rmi.xfem_minus_base_hierarchy[mf]))
            {
              rmi.xfem_minus_base_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_minus_val_base>
                (rmi.xfem_minus_base[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_GRAD_TEST: case GA_NODE_DIVERG_TEST:
          case GA_NODE_ELEMENTARY_GRAD_TEST:
          case GA_NODE_ELEMENTARY_DIVERG_TEST:
            if (rmi.grad.count(mf) == 0 ||
                !if_hierarchy.is_compatible(rmi.grad_hierarchy[mf])) {
              rmi.grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_grad_base>
                (rmi.grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_PLUS_GRAD_TEST: case GA_NODE_XFEM_PLUS_DIVERG_TEST:
            if (rmi.xfem_plus_grad.count(mf) == 0 ||
                !if_hierarchy.is_compatible(rmi.xfem_plus_grad_hierarchy[mf]))
            {
              rmi.xfem_plus_grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_plus_grad_base>
                (rmi.xfem_plus_grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_GRAD_TEST:
          case GA_NODE_XFEM_MINUS_DIVERG_TEST:
            if (rmi.xfem_minus_grad.count(mf) == 0 ||
                !if_hierarchy.is_compatible(rmi.xfem_minus_grad_hierarchy[mf]))
            {
              rmi.xfem_minus_grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_minus_grad_base>
                (rmi.xfem_minus_grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_HESS_TEST: case GA_NODE_ELEMENTARY_HESS_TEST:
            if (rmi.hess.count(mf) == 0 ||
                !if_hierarchy.is_compatible(rmi.hess_hierarchy[mf])) {
              rmi.hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_hess_base>
                (rmi.hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_PLUS_HESS_TEST:
            if (rmi.xfem_plus_hess.count(mf) == 0 ||
                !if_hierarchy.is_compatible(rmi.xfem_plus_hess_hierarchy[mf]))
            {
              rmi.xfem_plus_hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_plus_hess_base>
                (rmi.xfem_plus_hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_HESS_TEST:
            if (rmi.xfem_minus_hess.count(mf) == 0 ||
                !if_hierarchy.is_compatible(rmi.xfem_minus_hess_hierarchy[mf]))
            {
              rmi.xfem_minus_hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_minus_hess_base>
                (rmi.xfem_minus_hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;

          default : GMM_ASSERT1(false, "Internal error");
          }
          if (pgai) rmi.instructions.push_back(std::move(pgai));

          // The copy of the real_base_value
          switch(pnode->node_type) {
          case GA_NODE_VAL_TEST:
            // --> t(Qmult*ndof,Qmult*target_dim)
            if (mf->get_qdim() > 1 && mf->is_uniformly_vectorized()) {
              pnode->t.set_sparsity(1, mf->get_qdim());
              tensor_to_clear = true;
              pgai = std::make_shared<ga_instruction_copy_vect_val_base>
                (pnode->tensor(), rmi.base[mf], mf->get_qdim());
            } else {
              pgai = std::make_shared<ga_instruction_copy_val_base>
                (pnode->tensor(), rmi.base[mf], mf->get_qdim());
            }
            break;
          case GA_NODE_GRAD_TEST:
            // --> t(Qmult*ndof,Qmult*target_dim,N)
            if (mf->get_qdim() > 1 && mf->is_uniformly_vectorized()) {
              pnode->t.set_sparsity(2, mf->get_qdim());
              tensor_to_clear = true;
              pgai = std::make_shared<ga_instruction_copy_vect_grad_base>
                (pnode->tensor(), rmi.grad[mf], mf->get_qdim());
            } else {
              pgai = std::make_shared<ga_instruction_copy_grad_base>
                (pnode->tensor(), rmi.grad[mf], mf->get_qdim());
            }
            break;
          case GA_NODE_HESS_TEST:
            // --> t(Qmult*ndof,Qmult*target_dim,N,N)
            pgai = std::make_shared<ga_instruction_copy_hess_base>
              (pnode->tensor(), rmi.hess[mf], mf->get_qdim());
            if (mf->get_qdim() > 1 && mf->is_uniformly_vectorized())
              pnode->t.set_sparsity(3, mf->get_qdim());
            break;
          case GA_NODE_DIVERG_TEST:
            // --> t(Qmult*ndof)
            pgai = std::make_shared<ga_instruction_copy_diverg_base>
              (pnode->tensor(), rmi.grad[mf], mf->get_qdim());
            break;
          case GA_NODE_XFEM_PLUS_VAL_TEST:
            // -->t(Qmult*ndof,Qmult*target_dim)
            pgai = std::make_shared<ga_instruction_copy_val_base>
              (pnode->tensor(), rmi.xfem_plus_base[mf], mf->get_qdim());
            if (mf->get_qdim() > 1 && mf->is_uniformly_vectorized())
              pnode->t.set_sparsity(1, mf->get_qdim());
            break;
          case GA_NODE_XFEM_PLUS_GRAD_TEST:
            // --> t(Qmult*ndof,Qmult*target_dim,N)
            pgai = std::make_shared<ga_instruction_copy_grad_base>
              (pnode->tensor(), rmi.xfem_plus_grad[mf], mf->get_qdim());
            if (mf->get_qdim() > 1 && mf->is_uniformly_vectorized())
              pnode->t.set_sparsity(2, mf->get_qdim());
            break;
          case GA_NODE_XFEM_PLUS_HESS_TEST:
            // --> t(Qmult*ndof,Qmult*target_dim,N,N)
            pgai = std::make_shared<ga_instruction_copy_hess_base>
              (pnode->tensor(), rmi.xfem_plus_hess[mf], mf->get_qdim());
            if (mf->get_qdim() > 1 && mf->is_uniformly_vectorized())
              pnode->t.set_sparsity(3, mf->get_qdim());
            break;
          case GA_NODE_XFEM_PLUS_DIVERG_TEST:
            // --> t(Qmult*ndof)
            pgai = std::make_shared<ga_instruction_copy_diverg_base>
              (pnode->tensor(), rmi.xfem_plus_grad[mf], mf->get_qdim());
            break;
          case GA_NODE_XFEM_MINUS_VAL_TEST:
            // -->t(Qmult*ndof,Qmult*target_dim)
            pgai = std::make_shared<ga_instruction_copy_val_base>
              (pnode->tensor(), rmi.xfem_minus_base[mf], mf->get_qdim());
            if (mf->get_qdim() > 1 && mf->is_uniformly_vectorized())
              pnode->t.set_sparsity(1, mf->get_qdim());
            break;
          case GA_NODE_XFEM_MINUS_GRAD_TEST:
            // --> t(Qmult*ndof,Qmult*target_dim,N)
            pgai = std::make_shared<ga_instruction_copy_grad_base>
              (pnode->tensor(), rmi.xfem_minus_grad[mf], mf->get_qdim());
            if (mf->get_qdim() > 1 && mf->is_uniformly_vectorized())
              pnode->t.set_sparsity(2, mf->get_qdim());
            break;
          case GA_NODE_XFEM_MINUS_HESS_TEST:
            // --> t(Qmult*ndof,Qmult*target_dim,N,N)
            pgai = std::make_shared<ga_instruction_copy_hess_base>
              (pnode->tensor(), rmi.xfem_minus_hess[mf], mf->get_qdim());
            if (mf->get_qdim() > 1 && mf->is_uniformly_vectorized())
              pnode->t.set_sparsity(3, mf->get_qdim());
            break;
          case GA_NODE_XFEM_MINUS_DIVERG_TEST:
            // --> t(Qmult*ndof)
            pgai = std::make_shared<ga_instruction_copy_diverg_base>
              (pnode->tensor(), rmi.xfem_minus_grad[mf], mf->get_qdim());
            break;
          case GA_NODE_ELEMENTARY_VAL_TEST:
            { // --> t(Qmult*ndof,Qmult*target_dim)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos
                [std::make_tuple(pnode->elementary_name, mfo, mf)];
              pgai =
             std::make_shared<ga_instruction_elementary_trans_val_base>
                (pnode->tensor(), rmi.base[mf], mf->get_qdim(),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mfo, *mf, gis.ctx, eti.M, eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_GRAD_TEST:
            { // --> t(Qmult*ndof,Qmult*target_dim,N)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos
                [std::make_tuple(pnode->elementary_name, mfo, mf)];
              pgai =
            std::make_shared<ga_instruction_elementary_trans_grad_base>
                (pnode->tensor(), rmi.grad[mf], mf->get_qdim(),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mfo, *mf, gis.ctx, eti.M, eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_HESS_TEST:
            { // --> t(Qmult*ndof,Qmult*target_dim,N,N)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos
                [std::make_tuple(pnode->elementary_name, mfo, mf)];
              pgai =
            std::make_shared<ga_instruction_elementary_trans_hess_base>
                (pnode->tensor(), rmi.hess[mf], mf->get_qdim(),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mfo, *mf, gis.ctx, eti.M, eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_DIVERG_TEST:
            { // --> t(Qmult*ndof)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos
                [std::make_tuple(pnode->elementary_name, mfo, mf)];
              pgai =
          std::make_shared<ga_instruction_elementary_trans_diverg_base>
                (pnode->tensor(), rmi.grad[mf], mf->get_qdim(),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mfo, *mf, gis.ctx, eti.M, eti.icv);
            }
            break;
          default: break;
          }
          if (pgai) rmi.instructions.push_back(std::move(pgai));
        }
        workspace.add_temporary_interval_for_unreduced_variable(pnode->name);
      }
      break;

    case GA_NODE_SECONDARY_DOMAIN_VAL_TEST:
    case GA_NODE_SECONDARY_DOMAIN_GRAD_TEST:
    case GA_NODE_SECONDARY_DOMAIN_HESS_TEST:
    case GA_NODE_SECONDARY_DOMAIN_DIVERG_TEST:
      {
        GMM_ASSERT1(!function_case, "internal error");
        const mesh_fem *mf = workspace.associated_mf(pnode->name);
        const std::string &intn = pnode->interpolate_name;
        auto &sdi = rmi.secondary_domain_infos;

        fem_interpolation_context *pctx = &(sdi.ctx);
        papprox_integration pai = sdi.pai;
        psecondary_domain psd = workspace.secondary_domain(intn);
        if (mf) {
          GMM_ASSERT1(&(mf->linked_mesh()) == &(psd->mim().linked_mesh()),
                      "The finite element of variable " << pnode->name <<
                      " and the applied integration method have to be"
                      " defined on the same mesh for secondary domain");

          // An instruction for pfp update
          if (is_uniform) {
            if (sdi.pfps.count(mf) == 0) {
              sdi.pfps[mf] = 0;
              pgai = std::make_shared<ga_instruction_update_pfp>
                (*mf, sdi.pfps[mf], *pctx, gis.fp_pool);
              rmi.begin_instructions.push_back(std::move(pgai));
            }
          } else if (sdi.pfps.count(mf) == 0 ||
                     !if_hierarchy.is_compatible(rmi.pfp_hierarchy[mf])) {
            rmi.pfp_hierarchy[mf].push_back(if_hierarchy);
            sdi.pfps[mf] = 0;
            pgai = std::make_shared<ga_instruction_update_pfp>
              (*mf, sdi.pfps[mf], *pctx, gis.fp_pool);
            rmi.instructions.push_back(std::move(pgai));
          }

          // An instruction for the base value
          pgai = pga_instruction();
          switch (pnode->node_type) {
          case GA_NODE_SECONDARY_DOMAIN_VAL_TEST:
             if (sdi.base.count(mf) == 0 ||
                !(if_hierarchy.is_compatible(rmi.base_hierarchy[mf]))) {
              rmi.base_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_val_base>
                (sdi.base[mf], *pctx, *mf, sdi.pfps[mf]);
             }
             break;
          case GA_NODE_SECONDARY_DOMAIN_GRAD_TEST:
          case GA_NODE_SECONDARY_DOMAIN_DIVERG_TEST:
            if (sdi.grad.count(mf) == 0 ||
                !(if_hierarchy.is_compatible(rmi.grad_hierarchy[mf]))) {
              rmi.grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_grad_base>
                (sdi.grad[mf], *pctx, *mf, sdi.pfps[mf]);
            }
            break;
          case GA_NODE_SECONDARY_DOMAIN_HESS_TEST:
            if (sdi.hess.count(mf) == 0 ||
                !(if_hierarchy.is_compatible(rmi.hess_hierarchy[mf]))) {
              rmi.hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_hess_base>
                (sdi.hess[mf], *pctx, *mf, sdi.pfps[mf]);
            }
            break;
          default : GMM_ASSERT1(false, "Internal error");
          }
          if (pgai) rmi.instructions.push_back(std::move(pgai));

          // The copy of the real_base_value
          switch(pnode->node_type) {
          case GA_NODE_SECONDARY_DOMAIN_VAL_TEST:
            // --> t(Qmult*ndof,Qmult*target_dim)
            if (mf->get_qdim() > 1 && mf->is_uniformly_vectorized()) {
              pnode->t.set_sparsity(1, mf->get_qdim());
              tensor_to_clear = true;
              pgai = std::make_shared<ga_instruction_copy_vect_val_base>
                (pnode->tensor(), sdi.base[mf], mf->get_qdim());
            } else {
              pgai = std::make_shared<ga_instruction_copy_val_base>
                (pnode->tensor(), sdi.base[mf], mf->get_qdim());
            }
            break;
          case GA_NODE_SECONDARY_DOMAIN_GRAD_TEST:
            // --> t(Qmult*ndof,Qmult*target_dim,N)
            if (mf->get_qdim() > 1 && mf->is_uniformly_vectorized()) {
              pnode->t.set_sparsity(2, mf->get_qdim());
              tensor_to_clear = true;
              pgai = std::make_shared<ga_instruction_copy_vect_grad_base>
                (pnode->tensor(), sdi.grad[mf], mf->get_qdim());
            } else {
              pgai = std::make_shared<ga_instruction_copy_grad_base>
                (pnode->tensor(), sdi.grad[mf], mf->get_qdim());
            }
            break;
          case GA_NODE_SECONDARY_DOMAIN_HESS_TEST:
            // --> t(Qmult*ndof,Qmult*target_dim,N,N)
            pgai = std::make_shared<ga_instruction_copy_hess_base>
              (pnode->tensor(), sdi.hess[mf], mf->get_qdim());
            if (mf->get_qdim() > 1 && mf->is_uniformly_vectorized())
              pnode->t.set_sparsity(3, mf->get_qdim());
            break;
          case GA_NODE_SECONDARY_DOMAIN_DIVERG_TEST:
            // --> t(Qmult*ndof)
            pgai = std::make_shared<ga_instruction_copy_diverg_base>
              (pnode->tensor(), sdi.grad[mf], mf->get_qdim());
            break;
          default: break;
          }
          if (pgai) rmi.instructions.push_back(std::move(pgai));
        }
        workspace.add_temporary_interval_for_unreduced_variable(pnode->name);
      }
      break;

    case GA_NODE_INTERPOLATE_VAL_TEST: case GA_NODE_INTERPOLATE_GRAD_TEST:
    case GA_NODE_INTERPOLATE_HESS_TEST: case GA_NODE_INTERPOLATE_DIVERG_TEST:
      {
        const mesh_fem *mfn = workspace.associated_mf(pnode->name), **mfg = 0;
        const std::string &intn = pnode->interpolate_name;
        const mesh **m2 = &(rmi.interpolate_infos[intn].m);
        if (workspace.variable_group_exists(pnode->name)) {
          ga_instruction_set::variable_group_info &vgi =
            rmi.interpolate_infos[intn].groups_info[pnode->name];
          mfg = &(vgi.mf); mfn = 0;
        }

        if (pnode->node_type == GA_NODE_INTERPOLATE_VAL_TEST) {
          // --> t(Qmult*ndof,Qmult*target_dim)
          pgai = std::make_shared<ga_instruction_interpolate_val_base>
            (pnode->tensor(), m2, mfn, mfg, gis.ipt,
             workspace.qdim(pnode->name), rmi.interpolate_infos[intn],
             gis.fp_pool);
        } else if (pnode->node_type == GA_NODE_INTERPOLATE_GRAD_TEST) {
           // --> t(Qmult*ndof,Qmult*target_dim,N)
          pgai = std::make_shared<ga_instruction_interpolate_grad_base>
            (pnode->tensor(), m2, mfn, mfg, gis.ipt,
             workspace.qdim(pnode->name),
             rmi.interpolate_infos[intn], gis.fp_pool);
        } else if (pnode->node_type == GA_NODE_INTERPOLATE_HESS_TEST) {
           // --> t(Qmult*ndof,Qmult*target_dim,N,N)
          pgai = std::make_shared<ga_instruction_interpolate_hess_base>
            (pnode->tensor(), m2, mfn, mfg, gis.ipt,
             workspace.qdim(pnode->name),
             rmi.interpolate_infos[intn], gis.fp_pool);
        } else { // if (pnode->node_type == GA_NODE_INTERPOLATE_DIVERG_TEST) {
           // --> t(Qmult*ndof)
          pgai = std::make_shared<ga_instruction_interpolate_diverg_base>
            (pnode->tensor(), m2, mfn, mfg, gis.ipt,
             workspace.qdim(pnode->name),
             rmi.interpolate_infos[intn], gis.fp_pool);
        }
        rmi.instructions.push_back(std::move(pgai));
        workspace.add_temporary_interval_for_unreduced_variable(pnode->name);
      }
      break;

     case GA_NODE_OP:
       switch(pnode->op_type) {

       case GA_PLUS:
         if (pnode->tensor().size() == 1) {
           GA_DEBUG_ASSERT(child0->tensor().size() == 1,
                           "Internal error: child0 not scalar");
           GA_DEBUG_ASSERT(child1->tensor().size() == 1,
                           "Internal error: child1 not scalar");
           pgai = std::make_shared<ga_instruction_scalar_add>
             (pnode->tensor()[0], child0->tensor()[0], child1->tensor()[0]);
         } else {
           pgai = std::make_shared<ga_instruction_add>
             (pnode->tensor(), child0->tensor(), child1->tensor());
         }
         if (child0->t.sparsity() == child1->t.sparsity()
             && child0->t.qdim() == child1->t.qdim())
           pnode->t.set_sparsity(child0->t.sparsity(), child0->t.qdim());
         rmi.instructions.push_back(std::move(pgai));
         break;

       case GA_MINUS:
         if (pnode->tensor().size() == 1) {
           GA_DEBUG_ASSERT(child0->tensor().size() == 1,
                           "Internal error: child0 not scalar");
           GA_DEBUG_ASSERT(child1->tensor().size() == 1,
                           "Internal error: child1 not scalar");
           pgai = std::make_shared<ga_instruction_scalar_sub>
             (pnode->tensor()[0], child0->tensor()[0], child1->tensor()[0]);
         } else {
           pgai = std::make_shared<ga_instruction_sub>
             (pnode->tensor(), child0->tensor(), child1->tensor());
         }
         if (child0->t.sparsity() == child1->t.sparsity()
             && child0->t.qdim() == child1->t.qdim())
           pnode->t.set_sparsity(child0->t.sparsity(), child0->t.qdim());
         rmi.instructions.push_back(std::move(pgai));
         break;

       case GA_UNARY_MINUS:
         if (pnode->tensor().size() == 1) {
           GA_DEBUG_ASSERT(child0->tensor().size() == 1, "Internal error");
           pgai = std::make_shared<ga_instruction_scalar_scalar_mult>
             (pnode->tensor()[0], child0->tensor()[0], minus);
         } else {
           pgai = std::make_shared<ga_instruction_scalar_mult>
             (pnode->tensor(), child0->tensor(), minus);
         }
         pnode->t.set_sparsity(child0->t.sparsity(), child0->t.qdim());
         rmi.instructions.push_back(std::move(pgai));
         break;


       case GA_DOT: case GA_COLON: case GA_MULT:
         {
           size_type tps0 = child0->tensor_proper_size();
           size_type tps1 = child1->tensor_proper_size();
           size_type s1 = (tps0 * tps1) / pnode->tensor_proper_size();
           size_type s2 = size_type(round(sqrt(scalar_type(s1))));

           pgai = pga_instruction();
           if ((pnode->op_type == GA_DOT && dim1 <= 1) ||
               (pnode->op_type == GA_COLON && dim1 <= 2) ||
               (pnode->op_type == GA_MULT && dim0 == 4) ||
               (pnode->op_type == GA_MULT && dim1 <= 1) ||
               child0->tensor().size() == 1 || tps1 == 1) {

             if (child0->tensor().size() == 1 && child1->tensor().size() == 1) {
               pgai = std::make_shared<ga_instruction_scalar_scalar_mult>
                 (pnode->tensor()[0], child0->tensor()[0], child1->tensor()[0]);
             }
             else if (child0->tensor().size() == 1) {
               pnode->t.set_sparsity(child1->t.sparsity(), child1->t.qdim());
               pgai = std::make_shared<ga_instruction_scalar_mult>
                 (pnode->tensor(), child1->tensor(), child0->tensor()[0]);
             }
             else if (child1->tensor().size() == 1) {
               pnode->t.set_sparsity(child0->t.sparsity(), child0->t.qdim());
               pgai = std::make_shared<ga_instruction_scalar_mult>
                 (pnode->tensor(), child0->tensor(), child1->tensor()[0]);
             }
             else if (pnode->test_function_type < 3) {
               if (tps0 == 1) {
                 if (is_uniform) // Unrolled instruction
                   pgai = ga_uniform_instruction_simple_tmult
                     (pnode->tensor(), child0->tensor(), child1->tensor());
                 else
                   pgai = std::make_shared<ga_instruction_simple_tmult>
                     (pnode->tensor(), child0->tensor(), child1->tensor());
               } else {
                 if (tps1 == 1) {
                   if (is_uniform) // Unrolled instruction
                     pgai = ga_uniform_instruction_simple_tmult
                       (pnode->tensor(), child1->tensor(), child0->tensor());
                   else
                     pgai = std::make_shared<ga_instruction_simple_tmult>
                       (pnode->tensor(), child1->tensor(), child0->tensor());
                 } else if (is_uniform) // Unrolled instruction
                   pgai = ga_uniform_instruction_contraction_switch
                     (pnode->t, child0->t, child1->t, s2, tensor_to_clear);
                 else // Unrolled instruction
                   pgai = ga_instruction_contraction_switch
                     (pnode->t, child0->t, child1->t, s2, tensor_to_clear);
               }
             } else {
               if (child1->test_function_type == 1 ||
                   child1->test_function_type == 3) {
                 if (child1->test_function_type == 3 ||
                     child1->tensor_proper_size() <= s2) {
                   if (tps0 == 1) {
                     if (is_uniform) { // Unrolled instruction
                       pgai = ga_uniform_instruction_simple_tmult
                         (pnode->tensor(), child1->tensor(), child0->tensor());
                     } else
                       pgai = std::make_shared<ga_instruction_simple_tmult>
                         (pnode->tensor(), child1->tensor(), child0->tensor());
                   } else if (is_uniform) // Unrolled instruction
                     pgai = ga_uniform_instruction_contraction_switch
                       (pnode->t, child0->t, child1->t, s2, tensor_to_clear);
                   else // Unrolled instruction
                     pgai = ga_instruction_contraction_switch
                       (pnode->t, child0->t, child1->t, s2, tensor_to_clear);
                 } else
                   pgai = std::make_shared<ga_instruction_spec_contraction>
                     (pnode->tensor(), child1->tensor(), child0->tensor(), s2);
               } else if (child1->test_function_type == 0 ||
                          (child0->tensor_proper_size() == s2 &&
                           child1->tensor_proper_size() == s2)) {
                 if (tps0 == 1) {
                   if (is_uniform) { // Unrolled instruction
                     pgai = ga_uniform_instruction_simple_tmult
                       (pnode->tensor(), child0->tensor(), child1->tensor());
                   } else
                     pgai = std::make_shared<ga_instruction_simple_tmult>
                       (pnode->tensor(), child0->tensor(), child1->tensor());
                 } else {
                   if (is_uniform) // Unrolled instruction
                     pgai = ga_uniform_instruction_contraction_switch
                       (pnode->t, child1->t, child0->t, s2, tensor_to_clear);
                   else // Unrolled instruction
                     pgai = ga_instruction_contraction_switch
                       (pnode->t, child1->t, child0->t, s2, tensor_to_clear);
                 }
               } else {
                 if (child0->tensor_proper_size() == s2)
                   pgai = ga_uniform_instruction_contraction_switch
                     (pnode->t, child1->t, child0->t, s2, tensor_to_clear);
                 else if (child1->tensor_proper_size() == s2)
                   pgai = std::make_shared<ga_instruction_spec_contraction>
                     (pnode->tensor(), child0->tensor(), child1->tensor(), s2);
                 else
                   pgai = std::make_shared<ga_instruction_spec2_contraction>
                     (pnode->tensor(), child0->tensor(), child1->tensor(), s2);
               }
             }
           } else { // GA_MULT or GA_DOT for dim1 > 1 or GA_COLON for dim1 > 2
                    // and child1->tensor_proper_size() > 1
             if (pnode->test_function_type < 3) {
               if (tps0 == 1) {
                 if (is_uniform) // Unrolled instruction
                   pgai = ga_uniform_instruction_simple_tmult
                     (pnode->tensor(), child0->tensor(), child1->tensor());
                 else
                   pgai = std::make_shared<ga_instruction_simple_tmult>
                     (pnode->tensor(), child0->tensor(), child1->tensor());
               } else {
                 if (child1->test_function_type == 0)
                   pgai = std::make_shared<ga_instruction_matrix_mult>
                     (pnode->tensor(), child0->tensor(), child1->tensor(), s2);
                 else
                   pgai = std::make_shared<ga_instruction_matrix_mult_spec>
                     (pnode->tensor(), child0->tensor(), child1->tensor(),
                      s2, tps0/s2, tps1/s2);
               }
             } else {
               if (child0->tensor_proper_size() == 1) {
                 if (child0->test_function_type == 0 ||
                     child0->test_function_type == 1) {
                   if (is_uniform) // Unrolled instruction
                     pgai = ga_uniform_instruction_simple_tmult
                       (pnode->tensor(), child0->tensor(), child1->tensor());
                   else
                     pgai = std::make_shared<ga_instruction_simple_tmult>
                       (pnode->tensor(), child0->tensor(), child1->tensor());
                 } else
                   pgai = std::make_shared<ga_instruction_spec_tmult>
                     (pnode->tensor(), child1->tensor(), child0->tensor(),
                      tps1, tps0);
               } else {
                 if (child1->test_function_type == 0)
                   pgai = std::make_shared<ga_instruction_matrix_mult>
                     (pnode->tensor(), child0->tensor(), child1->tensor(), s2);
                 else if (child1->test_function_type == 2)
                   pgai = std::make_shared<ga_instruction_matrix_mult_spec>
                     (pnode->tensor(), child0->tensor(), child1->tensor(),
                      s2, tps0/s2, tps1/s2);
                 else
                   pgai = std::make_shared<ga_instruction_matrix_mult_spec2>
                     (pnode->tensor(), child0->tensor(), child1->tensor(),
                      s2, tps0/s2, tps1/s2);
               }
             }
           }
           rmi.instructions.push_back(std::move(pgai));
         }
         break;

       case GA_DIV:
         if (child0->tensor().size() == 1 && child1->tensor().size() == 1) {
           pgai = std::make_shared<ga_instruction_scalar_scalar_div>
             (pnode->tensor()[0], child0->tensor()[0], child1->tensor()[0]);
         } else if (child1->tensor().size() == 1) {
           pnode->t.set_sparsity(child0->t.sparsity(), child0->t.qdim());
           pgai = std::make_shared<ga_instruction_scalar_div>
             (pnode->tensor(), child0->tensor(), child1->tensor()[0]);
         } else GMM_ASSERT1(false, "Internal error");
         rmi.instructions.push_back(std::move(pgai));
         break;

       case GA_PRINT:
         pnode->t.set_to_copy(child0->t);
         pgai = std::make_shared<ga_instruction_print_tensor>
           (pnode->tensor(), child0, gis.ctx, gis.nbpt, gis.ipt);
         rmi.instructions.push_back(std::move(pgai));
         break;

       case GA_QUOTE:
         if (pnode->tensor_proper_size() > 1) {
           size_type n1 = child0->tensor_proper_size(0);
           size_type n2 = (child0->tensor_order() > 1) ?
             child0->tensor_proper_size(1) : 1;
           size_type nn = 1;
           for (size_type i = 2; i < child0->tensor_order(); ++i)
             nn *= child0->tensor_proper_size(i);
           if (child0->nb_test_functions() == 0)
             pgai = std::make_shared<ga_instruction_transpose_no_test>
               (pnode->tensor(), child0->tensor(), n1, n2, nn);
           else
             pgai = std::make_shared<ga_instruction_transpose>
               (pnode->tensor(), child0->tensor(), n1, n2, nn);
           rmi.instructions.push_back(std::move(pgai));
         } else {
           pnode->t.set_to_copy(child0->t);
         }
         break;

       case GA_SYM:
         if (pnode->tensor_proper_size() != 1) {
           pgai = std::make_shared<ga_instruction_sym>
             (pnode->tensor(), child0->tensor());
           rmi.instructions.push_back(std::move(pgai));
         } else {
           pnode->t.set_to_copy(child0->t);
         }
         break;

       case GA_SKEW:
         {
           pgai = std::make_shared<ga_instruction_skew>
             (pnode->tensor(), child0->tensor());
           rmi.instructions.push_back(std::move(pgai));
         }
         break;

       case GA_TRACE:
         {
           size_type N = (child0->tensor_proper_size() == 1) ? 1:size0.back();
           if (N == 1) {
             pnode->t.set_to_copy(child0->t);
           } else {
             pgai = std::make_shared<ga_instruction_trace>
               (pnode->tensor(), child0->tensor(), N);
             rmi.instructions.push_back(std::move(pgai));
           }
         }
         break;

       case GA_DEVIATOR:
         {
           size_type N = (child0->tensor_proper_size() == 1) ? 1:size0.back();
           pgai = std::make_shared<ga_instruction_deviator>
             (pnode->tensor(), child0->tensor(), N);
           rmi.instructions.push_back(std::move(pgai));
         }
         break;

       case GA_DOTMULT:

         if (child0->tensor().size() == 1 && child1->tensor().size() == 1) {
           pgai = std::make_shared<ga_instruction_scalar_scalar_mult>
             (pnode->tensor()[0], child0->tensor()[0], child1->tensor()[0]);
         } else if (child0->tensor().size() == 1) {
           pnode->t.set_sparsity(child1->t.sparsity(), child1->t.qdim());
           pgai = std::make_shared<ga_instruction_scalar_mult>
             (pnode->tensor(), child1->tensor(), child0->tensor()[0]);
         }
         else if (child1->tensor().size() == 1) {
           pnode->t.set_sparsity(child0->t.sparsity(), child0->t.qdim());
           pgai = std::make_shared<ga_instruction_scalar_mult>
             (pnode->tensor(), child0->tensor(), child1->tensor()[0]);
         }
         else if (child1->test_function_type == 0)
           pgai = std::make_shared<ga_instruction_dotmult>
             (pnode->tensor(), child0->tensor(), child1->tensor());
         else if (child0->test_function_type == 0)
           pgai = std::make_shared<ga_instruction_dotmult>
             (pnode->tensor(), child1->tensor(), child0->tensor());
         else if (child0->test_function_type == 1)
           pgai = std::make_shared<ga_instruction_dotmult_spec>
             (pnode->tensor(), child0->tensor(), child1->tensor());
         else
           pgai = std::make_shared<ga_instruction_dotmult_spec>
             (pnode->tensor(), child1->tensor(), child0->tensor());

         rmi.instructions.push_back(std::move(pgai));
         break;


       case GA_DOTDIV:
         if (child0->tensor().size() == 1 && child1->tensor().size() == 1) {
           pgai = std::make_shared<ga_instruction_scalar_scalar_div>
             (pnode->tensor()[0], child0->tensor()[0], child1->tensor()[0]);
         } else if (child1->tensor().size() == 1) {
           pnode->t.set_sparsity(child0->t.sparsity(), child0->t.qdim());
           pgai = std::make_shared<ga_instruction_scalar_div>
             (pnode->tensor(), child0->tensor(), child1->tensor()[0]);
         } else if (child1->test_function_type == 0) {
           pgai = std::make_shared<ga_instruction_dotdiv>
             (pnode->tensor(), child0->tensor(), child1->tensor());
         } else GMM_ASSERT1(false, "Internal error");
         rmi.instructions.push_back(std::move(pgai));
         break;


       case GA_TMULT:
         if (child0->tensor().size() == 1 && child1->tensor().size() == 1) {
           pgai = std::make_shared<ga_instruction_scalar_scalar_mult>
             (pnode->tensor()[0], child0->tensor()[0], child1->tensor()[0]);
         } else if (child0->tensor().size() == 1) {
           pnode->t.set_sparsity(child1->t.sparsity(), child1->t.qdim());
           pgai = std::make_shared<ga_instruction_scalar_mult>
             (pnode->tensor(), child1->tensor(), child0->tensor()[0]);
         }
         else if (child1->tensor().size() == 1) {
           pnode->t.set_sparsity(child0->t.sparsity(), child0->t.qdim());
           pgai = std::make_shared<ga_instruction_scalar_mult>
             (pnode->tensor(), child0->tensor(), child1->tensor()[0]);
         }
         else if (child1->test_function_type == 0) {
           if (is_uniform) // Unrolled instruction
             pgai = ga_uniform_instruction_simple_tmult
               (pnode->tensor(), child0->tensor(), child1->tensor());
           else
             pgai = std::make_shared<ga_instruction_simple_tmult>
               (pnode->tensor(), child0->tensor(), child1->tensor());
         } else if (child1->tensor_proper_size() == 1)
           pgai = std::make_shared<ga_instruction_spec2_tmult>
             (pnode->tensor(), child0->tensor(), child1->tensor());
         else
           pgai = std::make_shared<ga_instruction_spec_tmult>
             (pnode->tensor(), child0->tensor(), child1->tensor(),
              child0->tensor_proper_size(),
              child1->tensor_proper_size());

         rmi.instructions.push_back(std::move(pgai));
         break;

       default:GMM_ASSERT1(false, "Unexpected operation. Internal error.");
       }
       break;

    case GA_NODE_C_MATRIX:
      {
        if (pnode->test_function_type) {
          std::vector<const base_tensor *> components(pnode->children.size());
          for (size_type i = 0; i < pnode->children.size(); ++i)
            components[i]  = &(pnode->children[i]->tensor());
          pgai = std::make_shared<ga_instruction_c_matrix_with_tests>
            (pnode->tensor(), components);
        } else {
          std::vector<scalar_type *> components(pnode->children.size());
          for (size_type i = 0; i < pnode->children.size(); ++i)
            components[i]  = &(pnode->children[i]->tensor()[0]);
          pgai = std::make_shared<ga_instruction_simple_c_matrix>
            (pnode->tensor(), components);
        }
        rmi.instructions.push_back(std::move(pgai));
      }
      break;

    case GA_NODE_PARAMS:
      if (child0->node_type == GA_NODE_RESHAPE) {
        pgai = std::make_shared<ga_instruction_copy_tensor>(pnode->tensor(),
                                                            child1->tensor());
        rmi.instructions.push_back(std::move(pgai));
      } else if (child0->node_type == GA_NODE_CROSS_PRODUCT) {
        pga_tree_node child2 = pnode->children[2];
        if (child1->test_function_type==2 && child2->test_function_type==1)
          pgai = std::make_shared<ga_instruction_cross_product_tf>
            (pnode->tensor(), child2->tensor(), child1->tensor(), true);
        else if (child1->test_function_type || child2->test_function_type)
          pgai = std::make_shared<ga_instruction_cross_product_tf>
            (pnode->tensor(), child1->tensor(), child2->tensor(), false);
        else
          pgai = std::make_shared<ga_instruction_cross_product>
            (pnode->tensor(), child1->tensor(), child2->tensor());
        rmi.instructions.push_back(std::move(pgai));
      } else if (child0->node_type == GA_NODE_IND_MOVE_LAST) {
        size_type ind;
        ind = size_type(round(pnode->children[2]->tensor()[0])-1);
        size_type ii2 = 1;
        for (size_type i = 0; i < child1->tensor_order(); ++i)
          if (i>ind) ii2 *= child1->tensor_proper_size(i);
        size_type nn = child1->tensor_proper_size(ind);
        pgai = std::make_shared<ga_instruction_index_move_last>
          (pnode->tensor(), child1->tensor(), nn, ii2);
        rmi.instructions.push_back(std::move(pgai));
      } else if (child0->node_type == GA_NODE_SWAP_IND) {
        size_type ind[4];
        for (size_type i = 2; i < 4; ++i)
          ind[i] = size_type(round(pnode->children[i]->tensor()[0])-1);
        if (ind[2] > ind[3]) std::swap(ind[2], ind[3]);
        size_type ii2 = 1, ii3 = 1;
        for (size_type i = 0; i < child1->tensor_order(); ++i) {
          if (i>ind[2] && i<ind[3]) ii2 *= child1->tensor_proper_size(i);
          if (i>ind[3]) ii3 *= child1->tensor_proper_size(i);
        }
        size_type nn1 = child1->tensor_proper_size(ind[2]);
        size_type nn2 = child1->tensor_proper_size(ind[3]);

        pgai = std::make_shared<ga_instruction_swap_indices>
          (pnode->tensor(), child1->tensor(), nn1, nn2, ii2, ii3);
        rmi.instructions.push_back(std::move(pgai));
      } else if (child0->node_type == GA_NODE_CONTRACT) {
        std::vector<size_type> ind(2), indsize(2);
        pga_tree_node child2(0);
        if (pnode->children.size() == 4)
          { ind[0] = 2; ind[1] = 3; }
        else if (pnode->children.size() == 5)
          { ind[0] = 2; ind[1] = 4; child2 = pnode->children[3]; }
        else if (pnode->children.size() == 7) {
          ind.resize(4); indsize.resize(4);
          ind[0] = 2; ind[1] = 3; ind[2] = 5; ind[3] = 6;
          child2 = pnode->children[4];
        }
        size_type kk = 0, ll = 1;
        for (size_type i = 1; i < pnode->children.size(); ++i) {
          if (i == ind[kk]) {
            ind[kk] = size_type(round(pnode->children[i]->tensor()[0])-1);
            indsize[kk] =  pnode->children[ll]->tensor_proper_size(ind[kk]);
            ++kk;
          } else ll = i;
        }

        if (pnode->children.size() == 4) {
          size_type i1 = ind[0], i2 = ind[1];
          if (i1 > i2) std::swap(i1, i2);
          size_type ii2 = 1, ii3 = 1;
            for (size_type i = 0; i < child1->tensor_order(); ++i) {
              if (i > i1 && i < i2) ii2 *= child1->tensor_proper_size(i);
              if (i > i2) ii3 *= child1->tensor_proper_size(i);
            }
            pgai = std::make_shared<ga_instruction_contract_1_1>
              (pnode->tensor(), child1->tensor(), indsize[0], ii2, ii3);
        }
        else if (pnode->children.size() == 5) {
          // Particular cases should be detected (ii2=ii3=1 in particular).
          size_type i1 = ind[0], i2 = ind[1];
          size_type ii1 = 1, ii2 = 1, ii3 = 1, ii4 = 1;
          for (size_type i = 0; i < child1->tensor_order(); ++i) {
            if (i < i1) ii1 *= child1->tensor_proper_size(i);
            if (i > i1) ii2 *= child1->tensor_proper_size(i);
          }
          for (size_type i = 0; i < child2->tensor_order(); ++i) {
            if (i < i2) ii3 *= child2->tensor_proper_size(i);
            if (i > i2) ii4 *= child2->tensor_proper_size(i);
          }
          if (child1->test_function_type==1 && child2->test_function_type==2)
            pgai = std::make_shared<ga_instruction_contract_2_1_rev>
              (pnode->tensor(), child1->tensor(), child2->tensor(),
               indsize[0], ii1, ii2, ii3, ii4);
          else
            pgai = std::make_shared<ga_instruction_contract_2_1>
              (pnode->tensor(), child1->tensor(), child2->tensor(),
               indsize[0], ii1, ii2, ii3, ii4);
        }
        else if (pnode->children.size() == 7) {
          // Particular cases should be detected (ii2=ii3=1 in particular).
          size_type i1 = ind[0], i2 = ind[1], i3 = ind[2], i4 = ind[3];
          size_type nn1 = indsize[0], nn2 = indsize[1];
          size_type ii1 = 1, ii2 = 1, ii3 = 1, ii4 = 1, ii5 = 1, ii6 = 1;
          if (i1 > i2)
            { std::swap(i1, i2); std::swap(i3, i4); std::swap(nn1, nn2); }
          for (size_type i = 0; i < child1->tensor_order(); ++i) {
            if (i < i1) ii1 *= child1->tensor_proper_size(i);
            if (i > i1 && i < i2) ii2 *= child1->tensor_proper_size(i);
            if (i > i2) ii3 *= child1->tensor_proper_size(i);
          }
          for (size_type i = 0; i < child2->tensor_order(); ++i) {
            if (i < i3 && i < i4) ii4 *= child2->tensor_proper_size(i);
            if ((i > i3 && i < i4) || (i > i4 && i < i3))
              ii5 *= child2->tensor_proper_size(i);
            if (i > i3 && i > i4) ii6 *= child2->tensor_proper_size(i);
          }
          if (child1->test_function_type==1 && child2->test_function_type==2)
            pgai = std::make_shared<ga_instruction_contract_2_2_rev>
              (pnode->tensor(), child1->tensor(), child2->tensor(),
               nn1, nn2, ii1, ii2, ii3, ii4, ii5, ii6, i4 < i3);
          else
            pgai = std::make_shared<ga_instruction_contract_2_2>
              (pnode->tensor(), child1->tensor(), child2->tensor(),
               nn1, nn2, ii1, ii2, ii3, ii4, ii5, ii6, i4 < i3);
        }
        rmi.instructions.push_back(std::move(pgai));
      } else if (child0->node_type == GA_NODE_PREDEF_FUNC) {

        std::string name = child0->name;
        const ga_predef_function_tab &PREDEF_FUNCTIONS
          = dal::singleton<ga_predef_function_tab>::instance(0);
        ga_predef_function_tab::const_iterator it = PREDEF_FUNCTIONS.find(name);
        const ga_predef_function &F = it->second;
        size_type nbargs = F.nbargs();
        pga_tree_node child2 = (nbargs == 2) ? pnode->children[2] : child1;

        if (nbargs == 1) {
          if (child1->tensor().size() == 1) {
            if (F.ftype() == 0)
              pgai = std::make_shared<ga_instruction_eval_func_1arg_1res>
                (pnode->tensor()[0], child1->tensor()[0], F.f1());
            else
              pgai = std::make_shared<ga_instruction_eval_func_1arg_1res_expr>
                (pnode->tensor()[0], child1->tensor()[0], F);
          } else {
            if (F.ftype() == 0)
              pgai = std::make_shared<ga_instruction_eval_func_1arg>
                (pnode->tensor(), child1->tensor(), F.f1());
            else
              pgai = std::make_shared<ga_instruction_eval_func_1arg_expr>
                (pnode->tensor(), child1->tensor(), F);
          }
        } else {
          if (child1->tensor().size() == 1 && child2->tensor().size() == 1) {
            if (F.ftype() == 0)
              pgai = std::make_shared<ga_instruction_eval_func_2arg_1res>
                (pnode->tensor()[0], child1->tensor()[0], child2->tensor()[0],
                 F.f2());
            else
              pgai = std::make_shared<ga_instruction_eval_func_2arg_1res_expr>
                (pnode->tensor()[0], child1->tensor()[0], child2->tensor()[0],
                 F);
          } else if (child1->tensor().size() == 1) {
            if (F.ftype() == 0)
              pgai =
                std::make_shared<ga_instruction_eval_func_2arg_first_scalar>
                (pnode->tensor(), child1->tensor(), child2->tensor(), F.f2());
            else
              pgai =
              std::make_shared<ga_instruction_eval_func_2arg_first_scalar_expr>
                (pnode->tensor(), child1->tensor(), child2->tensor(), F);
          } else if (child2->tensor().size() == 1) {
            if (F.ftype() == 0)
              pgai =
                std::make_shared<ga_instruction_eval_func_2arg_second_scalar>
                (pnode->tensor(), child1->tensor(), child2->tensor(), F.f2());
            else
              pgai =
              std::make_shared<ga_instruction_eval_func_2arg_second_scalar_expr>
                (pnode->tensor(), child1->tensor(), child2->tensor(), F);
          } else {
            if (F.ftype() == 0)
              pgai = std::make_shared<ga_instruction_eval_func_2arg>
                (pnode->tensor(), child1->tensor(), child2->tensor(), F.f2());
            else
              pgai = std::make_shared<ga_instruction_eval_func_2arg_expr>
                (pnode->tensor(), child1->tensor(), child2->tensor(), F);
          }
        }
        rmi.instructions.push_back(std::move(pgai));

      } else if (child0->node_type == GA_NODE_SPEC_FUNC) {

        GMM_ASSERT1(false, "Internal error");

      } else if (child0->node_type == GA_NODE_OPERATOR) {

        ga_predef_operator_tab &PREDEF_OPERATORS
          = dal::singleton<ga_predef_operator_tab>::instance(0);
        ga_predef_operator_tab::T::iterator it
          = PREDEF_OPERATORS.tab.find(child0->name);
        const ga_nonlinear_operator &OP = *(it->second);
        ga_nonlinear_operator::arg_list args;
        for (size_type i = 1; i < pnode->children.size(); ++i)
          args.push_back(&(pnode->children[i]->tensor()));

        if (child0->der1 && child0->der2 == 0) {
          pgai = std::make_shared<ga_instruction_eval_derivative_OP>
             (pnode->tensor(), OP, args, child0->der1);
        } else if (child0->der1 && child0->der2) {
          pgai = std::make_shared<ga_instruction_eval_second_derivative_OP>
             (pnode->tensor(), OP, args, child0->der1, child0->der2);
        } else {
          pgai = std::make_shared<ga_instruction_eval_OP>(pnode->tensor(),
                                                          OP, args);
        }
        rmi.instructions.push_back(std::move(pgai));

      } else { // Access to a component of the tensor
        bgeot::multi_index mi1(size0.size()), indices;
        size_type nb_test = pnode->nb_test_functions();
        if (pnode->tensor().size() == 1) {
          for (size_type i = 0; i < child0->tensor_order(); ++i)
            mi1[i+nb_test] = size_type(round(pnode->children[i+1]->tensor()[0])-1);
          pgai = std::make_shared<ga_instruction_copy_scalar>
            (pnode->tensor()[0], child0->tensor()(mi1));
        } else {
          for (size_type i = 0; i < nb_test; ++i) indices.push_back(i);
          for (size_type i = 0; i < child0->tensor_order(); ++i) {
            if (pnode->children[i+1]->node_type != GA_NODE_ALLINDICES)
              mi1[i+nb_test]
                = size_type(round(pnode->children[i+1]->tensor()[0])- 1);
            else
              indices.push_back(i+nb_test);
          }
          pgai = std::make_shared<ga_instruction_tensor_slice>
            (pnode->tensor(), child0->tensor(), mi1, indices);
        }
        rmi.instructions.push_back(std::move(pgai));
      }

      break;

    default:GMM_ASSERT1(false, "Unexpected node type " << pnode->node_type
                        << " in compilation. Internal error.");
    }
    if (tensor_to_clear) {
      gmm::clear(pnode->tensor().as_vector());
      if (!is_uniform) {
        pgai = std::make_shared<ga_instruction_clear_tensor>(pnode->tensor());
        rmi.elt_instructions.push_back(std::move(pgai));
      }
    }
    rmi.node_list[pnode->hash_value].push_back(pnode);
  } // ga_compile_node

  void ga_compile_function(ga_workspace &workspace,
                           ga_instruction_set &gis, bool scalar) {
    for (size_type i = 0; i < workspace.nb_trees(); ++i) {
      const ga_workspace::tree_description &td = workspace.tree_info(i);

      gis.trees.push_back(*(td.ptree));
      pga_tree_node root = gis.trees.back().root;
      if (root) {
        GMM_ASSERT1(!scalar || (root->tensor().size() == 1),
                    "The result of the given expression is not a scalar");
        ga_instruction_set::region_mim rm(td.mim, td.rg, 0);
        gis.all_instructions[rm].m = td.m;
        ga_if_hierarchy if_hierarchy;
        ga_compile_node(root, workspace, gis, gis.all_instructions[rm],
                        *(td.m), true, if_hierarchy);

        gis.coeff = scalar_type(1);
        pga_instruction pgai;
        workspace.assembled_tensor() = root->tensor();
        pgai = std::make_shared<ga_instruction_add_to_coeff>
          (workspace.assembled_tensor(), root->tensor(), gis.coeff);
        gis.all_instructions[rm].instructions.push_back(std::move(pgai));
      }
    }
  }

  static bool ga_node_used_interpolates
  (const pga_tree_node pnode, const ga_workspace &workspace,
   std::map<std::string, std::set<std::string> > &interpolates,
   std::set<std::string> &interpolates_der) {
    bool found = false;
    bool intrpl(pnode->node_type == GA_NODE_INTERPOLATE_VAL ||
                pnode->node_type == GA_NODE_INTERPOLATE_GRAD ||
                pnode->node_type == GA_NODE_INTERPOLATE_HESS ||
                pnode->node_type == GA_NODE_INTERPOLATE_DIVERG);
    bool intrpl_test(pnode->node_type == GA_NODE_INTERPOLATE_VAL_TEST ||
                     pnode->node_type == GA_NODE_INTERPOLATE_GRAD_TEST ||
                     pnode->node_type == GA_NODE_INTERPOLATE_HESS_TEST ||
                     pnode->node_type == GA_NODE_INTERPOLATE_DIVERG_TEST);

    if (intrpl || intrpl_test ||
        pnode->node_type == GA_NODE_INTERPOLATE_FILTER ||
        pnode->node_type == GA_NODE_INTERPOLATE_X ||
        pnode->node_type == GA_NODE_INTERPOLATE_NORMAL) {
      interpolates[pnode->interpolate_name].size();
      if (intrpl || intrpl_test) {
        if (workspace.variable_group_exists(pnode->name))
          interpolates[pnode->interpolate_name].insert(pnode->name);
      }
      found = true;
    }
    if (pnode->node_type == GA_NODE_INTERPOLATE_DERIVATIVE) {
      interpolates_der.insert(pnode->interpolate_name_der);
      interpolates[pnode->interpolate_name_der].size();
      if (workspace.variable_group_exists(pnode->name))
        interpolates[pnode->interpolate_name_der].insert(pnode->name);
    }
    for (size_type i = 0; i < pnode->children.size(); ++i)
      found = ga_node_used_interpolates(pnode->children[i], workspace,
                                        interpolates, interpolates_der)
        || found;
    return found;
  }


  static void ga_compile_interpolate_trans
  (const pga_tree_node pnode, const ga_workspace &workspace,
   ga_instruction_set &gis, ga_instruction_set::region_mim_instructions &rmi,
   const mesh &m) {

    std::set<std::string> interpolates_der;
    std::map<std::string, std::set<std::string> > transformations;
    ga_node_used_interpolates(pnode, workspace, transformations,
                              interpolates_der);

    for (const auto &transformation : transformations) {
      const std::string &transname = transformation.first;
      bool compute_der = (interpolates_der.count(transname) != 0);
      if (rmi.transformations.count(transname) == 0 ||
          (compute_der && rmi.transformations_der.count(transname) == 0)) {
        rmi.transformations[transname].size();
        gis.transformations.insert(transname);
        if (compute_der) rmi.transformations_der.insert(transname);
        pga_instruction pgai;
        if (transname.compare("neighbor_element") == 0 ||
            transname.compare("neighbour_elt") == 0) {
          pgai = std::make_shared<ga_instruction_neighbor_transformation_call>
            (workspace, rmi.interpolate_infos[transname],
             workspace.interpolate_transformation(transname), gis.ctx,
             m, gis.ipt, gis.pai, gis.gp_pool, gis.neighbor_corresp);
        } else {
          pgai = std::make_shared<ga_instruction_transformation_call>
            (workspace, rmi.interpolate_infos[transname],
             workspace.interpolate_transformation(transname), gis.ctx,
             gis.Normal, m, compute_der);
        }
        if (pgai) rmi.instructions.push_back(std::move(pgai));
      }

      for (const std::string &nodename : transformation.second) {
        if (rmi.transformations[transname].count(nodename) == 0) {
          auto&& inin = rmi.interpolate_infos[transname];
          pga_instruction pgai =
            std::make_shared<ga_instruction_update_group_info>
            (workspace, gis, inin, nodename, inin.groups_info[nodename]);
          rmi.instructions.push_back(std::move(pgai));
          rmi.transformations[transname].insert(nodename);
        }
      }
    }
  }

  void ga_compile_interpolation(ga_workspace &workspace,
                                ga_instruction_set &gis) {
    gis.transformations.clear();
    gis.all_instructions.clear();
    for (size_type i = 0; i < workspace.nb_trees(); ++i) {
      const ga_workspace::tree_description &td = workspace.tree_info(i);
      if (td.operation != ga_workspace::ASSEMBLY) {
        gis.trees.push_back(*(td.ptree));

        // Semantic analysis mainly to evaluate fixed size variables and data
        const mesh *m = td.m;
        GMM_ASSERT1(m, "Internal error");
        ga_semantic_analysis(gis.trees.back(), workspace, *m,
                             ref_elt_dim_of_mesh(*m, *(td.rg)), true, false);
        pga_tree_node root = gis.trees.back().root;
        if (root) {
          // Compile tree
          ga_instruction_set::region_mim rm(td.mim, td.rg, 0);
          auto &rmi = gis.all_instructions[rm];
          rmi.m = td.m;
          rmi.im = td.mim;
          // rmi.interpolate_infos.clear();
          ga_compile_interpolate_trans(root, workspace, gis, rmi, *(td.m));
          ga_compile_node(root, workspace, gis,rmi, *(td.m), false,
                          rmi.current_hierarchy);

          // After compile tree
          workspace.assembled_tensor() = root->tensor();
          pga_instruction pgai = std::make_shared<ga_instruction_add_to>
            (workspace.assembled_tensor(), root->tensor());
          rmi.instructions.push_back(std::move(pgai));
        }
      }
    }
  }


  struct var_set : std::map<std::string,size_type> {
    // This class indexes variable names in the order of their addition
    size_type operator[](const std::string &name) {
      if (name.empty()) return size_type(-1);
      size_type id = size();
      auto it = find(name);
      if (it == end()) {
        emplace(name, id);
        return id;
      }
      return it->second;
    }
    std::string operator[](const size_type &id) const {
      for (const auto &key_value : *this) // brute force reverse search
        if (key_value.second == id)
          return key_value.first;
      return std::string("");
    }
  };


  struct condensation_description {
    var_set Ivars, Jvars, Qvars; // sets of variables involved in condensation
    // Clusters of intercoupled condensed variables and subdiagonally coupled
    // primary variables for each cluster
    std::vector<std::set<size_type>> Qclusters, Jclusters;
    // Each element of Qclusters contains a group of intercoupled condensed
    // variables. Due to the couplings within each group, all variables of the
    // same group need to be condensed out simultaneously. Per definition two
    // clusters cannot share a common variable.
    // indexing of groups
    std::vector<size_type> cluster_of_Qvar;
    // Matrices of pointers to submatrices for all coupling terms
    gmm::dense_matrix<base_tensor *> KQQ,        // diagonal
                                     KQJ, KQJpr, // subdiagonal
                                     KIQ,        // superdiagonal
                                     KIJ;        // outcome
    std::vector<base_tensor *> RI,   // res. vector of coupled primary variables
                               RQpr; // partial solution for condensed variables (initially stores residuals)
  };

  void ga_compile(ga_workspace &workspace,
                  ga_instruction_set &gis, size_type order, bool condensation) {
    gis.transformations.clear();
    gis.all_instructions.clear();
    gis.unreduced_terms.clear();
    workspace.clear_temporary_variable_intervals();

    std::map<const ga_instruction_set::region_mim, condensation_description>
      condensations;

    if (condensation && order == 2) {
      for (size_type i = 0; i < workspace.nb_trees(); ++i) {
        ga_workspace::tree_description &td = workspace.tree_info(i);
        if (td.order != 2 && td.order != size_type(-1))
          continue;
        ga_tree tree(*(td.ptree)); // temporary tree (not used later)
        ga_semantic_analysis(tree, workspace, td.mim->linked_mesh(),
                            ref_elt_dim_of_mesh(td.mim->linked_mesh(),*(td.rg)),
                            true, false);
        pga_tree_node root = tree.root;
        if (root) {
          const bool
            v1_is_intern = workspace.is_internal_variable(root->name_test1),
            v2_is_intern = workspace.is_internal_variable(root->name_test2);
          if (v1_is_intern || v2_is_intern) {
            GMM_ASSERT1(tree.secondary_domain.empty(),
                        "Condensed variable cannot be used in secondary domain");

            for (const auto &key_val : condensations) {
              const ga_instruction_set::region_mim rm0 = key_val.first;
              const condensation_description &CC0 = key_val.second;
              if (rm0.mim() == td.mim && rm0.region() != td.rg
                  && (CC0.Qvars.count(root->name_test1) ||
                      CC0.Qvars.count(root->name_test2))) {
                mesh_region intrsct = getfem::mesh_region::intersection
                                      (*(rm0.region()), *(td.rg));
                GMM_ASSERT1(intrsct.is_empty(),
                            "Cannot condense coupled variables between "
                            "intersecting regions");
              }
            }
            const ga_instruction_set::region_mim rm(td.mim, td.rg, nullptr);

            condensation_description &CC = condensations[rm];
            size_type
              q1 = v1_is_intern ? CC.Qvars[root->name_test1] : size_type(-1),
              q2 = v2_is_intern ? CC.Qvars[root->name_test2] : size_type(-1);
            GMM_ASSERT1(q1 != size_type(-1) || q2 != size_type(-1), "Error");
            std::vector<size_type> selected_clusters;
            for (size_type j=0; j < CC.Qclusters.size(); ++j)
              if (CC.Qclusters[j].count(q1) || CC.Qclusters[j].count(q2))
                selected_clusters.push_back(j);

            if (selected_clusters.empty()) { // create new cluster
              CC.Qclusters.push_back(std::set<size_type>());
              if (q1 != size_type(-1)) CC.Qclusters.back().insert(q1);
              if (q2 != size_type(-1)) CC.Qclusters.back().insert(q2);
            } else { // add into existing cluster / merge clusters together
              auto &target = CC.Qclusters[selected_clusters[0]];
              if (q1 != size_type(-1)) target.insert(q1);
              if (q2 != size_type(-1)) target.insert(q2);
              for (size_type j=selected_clusters.size()-1; j > 1; --j) {
                auto &source = CC.Qclusters[selected_clusters[j]];
                target.insert(source.begin(), source.end());
                CC.Qclusters.erase(CC.Qclusters.begin() + selected_clusters[j]);
              }
            }
          } // is_internal_variable
        } // if (root)
      } // for (size_type i = 0; i < workspace.nb_trees(); ++i)

      for (auto &key_value : condensations) {
        condensation_description &CC = key_value.second;
        //for (const auto &cluster : CC.Qclusters) {
        //  cout << "Clusters of coupled variables:" << endl;
        //  for (const auto &varid : cluster) cout << "/" << CC.Qvars[varid];
        //  cout << "/" << endl;
        //}
        size_type Qsize = CC.Qvars.size();

        // Jclusters will hold all J variables each cluster is coupled to
        CC.Jclusters.resize(CC.Qclusters.size());

        CC.cluster_of_Qvar.resize(Qsize);
        for (size_type i=0; i < CC.Qclusters.size(); ++i)
          for (const size_type &var : CC.Qclusters[i])
            CC.cluster_of_Qvar[var] = i;

        // Qvars: all condensed variables
        // Qclusters: definition of clusters of intercoupled variables of Qvars
        // cluster_of_Qvar: dictionary for which cluster each variable belongs to
        CC.KQQ.resize(Qsize, Qsize);
        CC.RQpr.resize(Qsize);
        for (size_type q=0; q < Qsize; ++q) {
          bgeot::multi_index mi(1);
          mi[0] = workspace.associated_im_data(CC.Qvars[q]) ->nb_tensor_elem();
          gis.condensation_tensors.push_back // memory allocation
                                   (std::make_shared<base_tensor>(mi));
          CC.RQpr[q] = gis.condensation_tensors.back().get();
        }
      }
    } // if (condensation && order == 2)

    std::array<ga_workspace::operation_type,3>
      phases{ga_workspace::PRE_ASSIGNMENT,
             ga_workspace::ASSEMBLY,
             ga_workspace::POST_ASSIGNMENT};
    for (const auto &phase : phases) {

      for (size_type i = 0; i < workspace.nb_trees(); ++i) {
        ga_workspace::tree_description &td = workspace.tree_info(i);
        if (td.operation != phase)
          continue; // skip this tree in this phase

        if (td.order == order || td.order == size_type(-1)) {
          std::list<ga_tree> &trees = (phase == ga_workspace::ASSEMBLY)
                                    ? gis.trees
                                    : gis.interpolation_trees;
          trees.push_back(*(td.ptree));
          // Semantic analysis mainly to evaluate fixed size variables and data
          ga_semantic_analysis(trees.back(), workspace, td.mim->linked_mesh(),
                            ref_elt_dim_of_mesh(td.mim->linked_mesh(),*(td.rg)),
                            true, false);
          pga_tree_node root = trees.back().root;
          if (root) {
            // Compile tree
            // cout << "Will compile "; ga_print_node(root, cout); cout << endl;

            psecondary_domain psd(0);
            if (trees.back().secondary_domain.size())
              psd = workspace.secondary_domain(trees.back().secondary_domain);
            ga_instruction_set::region_mim rm(td.mim, td.rg, psd);
            auto &rmi = gis.all_instructions[rm];
            rmi.m = td.m;
            rmi.im = td.mim;
            // rmi.interpolate_infos.clear();
            ga_compile_interpolate_trans(root, workspace, gis, rmi, *(td.m));
            ga_compile_node(root, workspace, gis, rmi, *(td.m), false,
                            rmi.current_hierarchy);
            // cout << "compilation finished "; ga_print_node(root, cout);
            // cout << endl;

            if (phase != ga_workspace::ASSEMBLY) { // Assignment/interpolation
              if (!td.varname_interpolation.empty()) {
                auto *imd
                  = workspace.associated_im_data(td.varname_interpolation);
                auto &V = const_cast<model_real_plain_vector &>
                  (workspace.value(td.varname_interpolation));
                GMM_ASSERT1(imd, "Internal error");
                auto pgai = std::make_shared<ga_instruction_assignment>
                  (root->tensor(), V, gis.ctx, imd);
                rmi.instructions.push_back(std::move(pgai));
              }
            } else { // Addition of an assembly instruction
              pga_instruction pgai;
              switch(order) {
              case 0: {
                workspace.assembled_tensor() = root->tensor();
                pgai = std::make_shared<ga_instruction_add_to_coeff>
                  (workspace.assembled_tensor(), root->tensor(), gis.coeff);
                break;
              }
              case 1: {
                GMM_ASSERT1(root->tensor_proper_size() == 1,
                            "Invalid vector or tensor quantity. An order 1 "
                            "weak form has to be a scalar quantity");
                const mesh_fem * const
                  mf = workspace.associated_mf(root->name_test1);
                const im_data * const
                  imd = workspace.associated_im_data(root->name_test1);
                workspace.add_temporary_interval_for_unreduced_variable
                  (root->name_test1);

                base_vector &Vu = workspace.unreduced_vector(),
                            &Vr = workspace.assembled_vector();
                if (mf) {
                  const std::string &intn1 = root->interpolate_name_test1;
                  bool secondary = !intn1.empty() &&
                                   workspace.secondary_domain_exists(intn1);
                  fem_interpolation_context
                    &ctx = intn1.empty() ? gis.ctx
                         : (secondary ? rmi.secondary_domain_infos.ctx
                                      : rmi.interpolate_infos[intn1].ctx);
                  bool interpolate =
                    !(intn1.empty() || intn1 == "neighbor_element"
                                    || intn1 == "neighbour_elt" || secondary);

                  if (intn1.size() && !secondary &&
                      workspace.variable_group_exists(root->name_test1)) {
                    ga_instruction_set::variable_group_info
                      &vgi = rmi.interpolate_infos[intn1]
                                .groups_info[root->name_test1];
                    pgai = std::make_shared<ga_instruction_vector_assembly_mf>
                           (root->tensor(), Vr, Vu, ctx,
                            vgi.I, vgi.mf, vgi.reduced_mf,
                            gis.coeff, gis.nbpt, gis.ipt, interpolate);
                    for (const std::string &name
                         : workspace.variable_group(root->name_test1))
                      gis.unreduced_terms.emplace(name, "");
                  } else {
                    base_vector &V = mf->is_reduced() ? Vu : Vr;
                    const gmm::sub_interval
                      &I = mf->is_reduced()
                         ? workspace.temporary_interval_of_variable
                                     (root->name_test1)
                         : workspace.interval_of_variable(root->name_test1);
                    pgai = std::make_shared<ga_instruction_vector_assembly_mf>
                           (root->tensor(), V, ctx, I, *mf,
                            gis.coeff, gis.nbpt, gis.ipt, interpolate);
                    if (mf->is_reduced())
                      gis.unreduced_terms.emplace(root->name_test1, "");
                  }
                } else if (imd) {
                  GMM_ASSERT1(root->interpolate_name_test1.size() == 0,
                              "Interpolate transformation on integration "
                              "point variable");
                  if (!workspace.is_internal_variable(root->name_test1) ||
                      condensation)
                    pgai = std::make_shared<ga_instruction_vector_assembly_imd>
                           (root->tensor(), Vr, gis.ctx,
                            workspace.interval_of_variable(root->name_test1),
                            *imd, gis.coeff, gis.ipt);
                    // Variable root->name_test1 can be internal or not
                } else {
                  pgai = std::make_shared<ga_instruction_vector_assembly>
                         (root->tensor(), Vr,
                          workspace.interval_of_variable(root->name_test1),
                          gis.coeff);
                }
                break;
              }
              case 2: {
                GMM_ASSERT1(root->tensor_proper_size() == 1,
                            "Invalid vector or tensor quantity. An order 2 "
                            "weak form has to be a scalar quantity");
                const mesh_fem *mf1=workspace.associated_mf(root->name_test1),
                               *mf2=workspace.associated_mf(root->name_test2);
                const im_data
                  *imd1 = workspace.associated_im_data(root->name_test1),
                  *imd2 = workspace.associated_im_data(root->name_test2);
                const std::string &intn1 = root->interpolate_name_test1,
                                  &intn2 = root->interpolate_name_test2;
                bool secondary1 = intn1.size() &&
                                  workspace.secondary_domain_exists(intn1);
                bool secondary2 = intn2.size() &&
                                  workspace.secondary_domain_exists(intn2);
                fem_interpolation_context
                  &ctx1 = intn1.empty() ? gis.ctx
                        : (secondary1 ? rmi.secondary_domain_infos.ctx
                                      : rmi.interpolate_infos[intn1].ctx),
                  &ctx2 = intn2.empty() ? gis.ctx
                        : (secondary2 ? rmi.secondary_domain_infos.ctx
                                      : rmi.interpolate_infos[intn2].ctx);
                bool interpolate = !(intn1.empty() || intn1 == "neighbor_element"
                                     || intn1 == "neighbour_elt"
                                     || secondary1) ||
                                   !(intn2.empty() || intn2 == "neighbor_element"
                                     || intn2 == "neighbour_elt"
                                     || secondary2);

                workspace.add_temporary_interval_for_unreduced_variable
                  (root->name_test1);
                workspace.add_temporary_interval_for_unreduced_variable
                  (root->name_test2);

                bool has_var_group1 = (!intn1.empty() && !secondary1 &&
                                       workspace.variable_group_exists
                                                 (root->name_test1));
                bool has_var_group2 = (!intn2.empty() && !secondary2 &&
                                       workspace.variable_group_exists
                                                 (root->name_test2));
                bool simple = !interpolate &&
                              !has_var_group1 && !has_var_group2 &&
                              mf1 && !(mf1->is_reduced()) &&
                              mf2 && !(mf2->is_reduced());

                // ga instructions write into one of the following matrices
                auto &Krr = workspace.assembled_matrix();
                auto &Kru = workspace.col_unreduced_matrix();
                auto &Kur = workspace.row_unreduced_matrix();
                auto &Kuu = workspace.row_col_unreduced_matrix();

                if (simple) { // --> Krr
                  const gmm::sub_interval
                    &I1 = workspace.interval_of_variable(root->name_test1),
                    &I2 = workspace.interval_of_variable(root->name_test2);
                  const scalar_type
                    &alpha1 = workspace.factor_of_variable(root->name_test1),
                    &alpha2 = workspace.factor_of_variable(root->name_test2);
                  if (mf1->get_qdim() == 1 && mf2->get_qdim() == 1)
                    pgai = std::make_shared
                      <ga_instruction_matrix_assembly_standard_scalar>
                      (root->tensor(), Krr, ctx1, ctx2, I1, I2, mf1, mf2,
                       alpha1, alpha2, gis.coeff, gis.nbpt, gis.ipt);
                  else if (root->sparsity() == 10 && root->t.qdim() == 2)
                    pgai = std::make_shared
                      <ga_instruction_matrix_assembly_standard_vector_opt10<2>>
                      (root->tensor(), Krr, ctx1, ctx2, I1, I2, mf1, mf2,
                       alpha1, alpha2, gis.coeff, gis.nbpt, gis.ipt);
                  else if (root->sparsity() == 10 && root->t.qdim() == 3)
                    pgai = std::make_shared
                      <ga_instruction_matrix_assembly_standard_vector_opt10<3>>
                      (root->tensor(), Krr, ctx1, ctx2, I1, I2, mf1, mf2,
                       alpha1, alpha2, gis.coeff, gis.nbpt, gis.ipt);
                  else
                    pgai = std::make_shared
                      <ga_instruction_matrix_assembly_standard_vector>
                      (root->tensor(), Krr, ctx1, ctx2, I1, I2, mf1, mf2,
                       alpha1, alpha2, gis.coeff, gis.nbpt, gis.ipt);
                } else if (condensation &&
                           workspace.is_internal_variable(root->name_test1) &&
                           workspace.is_internal_variable(root->name_test2)) {
                  // diagonal condensation matrix KQQ
                  // Only memory allocation, gathering of relevant pointers
                  // and data summation instructions
                  GMM_ASSERT1(imd1 && imd2, "Internal error");
                  GMM_ASSERT1(!interpolate, "Internal error");
                  size_type s1 = imd1->nb_tensor_elem();
                  size_type s2 = imd2->nb_tensor_elem();

                  condensation_description &CC = condensations[rm];
                  GMM_ASSERT1(CC.Qvars.count(root->name_test1) > 0 &&
                              CC.Qvars.count(root->name_test2) > 0,
                              "Internal error");
                  size_type q1 = CC.Qvars[root->name_test1],
                            q2 = CC.Qvars[root->name_test2];
                  if (!CC.KQQ(q1,q2)) {
                    // allocate a new matrix
                    gis.condensation_tensors.push_back
                      (std::make_shared<base_tensor>(s1,s2));
                    CC.KQQ(q1,q2) = gis.condensation_tensors.back().get();
                    pgai = std::make_shared<ga_instruction_copy_vect>
                      (CC.KQQ(q1,q2)->as_vector(), root->tensor().as_vector());
                  } else {
                    // addition instruction to the previously allocated matrix
                    pgai = std::make_shared<ga_instruction_add_to>
                      (*CC.KQQ(q1,q2), root->tensor());
                  }
                  rmi.instructions.push_back(std::move(pgai));
                } else if (condensation &&
                           workspace.is_internal_variable(root->name_test1)) {
                  // subdiagonal condensation matrix KQJ
                  // Only memory allocation, gathering of relevant pointers
                  // and data summation instructions
                  GMM_ASSERT1(imd1, "Internal error");
                  GMM_ASSERT1(!interpolate, "Internal error");
                  size_type s1 = imd1->nb_tensor_elem();

                  condensation_description &CC = condensations[rm];
                  GMM_ASSERT1(CC.Qvars.count(root->name_test1),
                              "Internal error");
                  size_type q1 = CC.Qvars[root->name_test1],
                            j2 = CC.Jvars[root->name_test2];
                  CC.Jclusters[CC.cluster_of_Qvar[q1]].insert(j2);
                  if (q1 >= CC.KQJ.nrows() || j2 >= CC.KQJ.ncols())
                    CC.KQJ.resize(std::max(CC.KQJ.nrows(), q1+1),
                                  std::max(CC.KQJ.ncols(), j2+1));
                  if (!CC.KQJ(q1,j2)) {
                    // allocate a new matrix. Here we do not know the size as
                    // it may change dynamically, but for now, just use the
                    // size of root->tensor()
                    gis.condensation_tensors.push_back
                      (std::make_shared<base_tensor>(root->tensor()));
                    GMM_ASSERT1(root->tensor().size(0) == s1, "Internal error");
                    CC.KQJ(q1,j2) = gis.condensation_tensors.back().get();
                    pgai = std::make_shared<ga_instruction_copy_vect>
                      (CC.KQJ(q1,j2)->as_vector(), root->tensor().as_vector());
                  } else {
                    // an extra matrix for this entry has already been
                    // allocated, so just add the current tensor to it
                    pgai = std::make_shared<ga_instruction_add_to>
                      (*CC.KQJ(q1,j2), root->tensor());
                  }
                  rmi.instructions.push_back(std::move(pgai));
                } else if (condensation &&
                           workspace.is_internal_variable(root->name_test2)) {
                  // superdiagonal condensation matrix KIQ
                  // Only memory allocation, gathering of relevant pointers
                  // and data summation instructions
                  GMM_ASSERT1(imd2, "Internal error");
                  GMM_ASSERT1(!interpolate, "Internal error");
                  size_type s2 = imd2->nb_tensor_elem();

                  condensation_description &CC = condensations[rm];
                  GMM_ASSERT1(CC.Qvars.count(root->name_test2),
                              "Internal error");
                  size_type i1 = CC.Ivars[root->name_test1],
                            q2 = CC.Qvars[root->name_test2];
                  if (i1 >= CC.KIQ.nrows() || q2 >= CC.KIQ.ncols())
                    CC.KIQ.resize(std::max(CC.KIQ.nrows(), i1+1),
                                  std::max(CC.KIQ.ncols(), q2+1));
                  if (!CC.KIQ(i1,q2)) {
                    // allocate a new matrix. Here we do not know the size as
                    // it may change dynamically, but for now, just use the
                    // size of root->tensor()
                    gis.condensation_tensors.push_back
                      (std::make_shared<base_tensor>(root->tensor()));
                    GMM_ASSERT1(root->tensor().size(1) == s2,
                                "Internal error");
                    CC.KIQ(i1,q2) = gis.condensation_tensors.back().get();
                    pgai = std::make_shared<ga_instruction_copy_vect>
                      (CC.KIQ(i1,q2)->as_vector(), root->tensor().as_vector());
                  } else {
                    // an extra matrix for this entry has already been
                    // allocated, so just add the current tensor to it
                    pgai = std::make_shared<ga_instruction_add_to>
                      (*CC.KIQ(i1,q2), root->tensor());
                  }
                  rmi.instructions.push_back(std::move(pgai));
                } else if (!workspace.is_internal_variable(root->name_test1) &&
                           !workspace.is_internal_variable(root->name_test2)) {

                  if ((mf1 && mf1->is_reduced()) || (mf2 && mf2->is_reduced())
                      || has_var_group1 || has_var_group2)
                    gis.unreduced_terms.emplace(root->name_test1,
                                                root->name_test2);

                  auto &Kxu = (mf1 && mf1->is_reduced()) ? Kuu : Kru;
                  auto &Kxr = (mf1 && mf1->is_reduced()) ? Kur : Krr;
                  auto &Kux = (mf2 && mf2->is_reduced()) ? Kuu : Kur;
                  auto &Krx = (mf2 && mf2->is_reduced()) ? Kru : Krr;
                  auto &Kxx = (mf2 && mf2->is_reduced()) ? Kxu : Kxr;

                  const scalar_type
                    &alpha1 = workspace.factor_of_variable(root->name_test1),
                    &alpha2 = workspace.factor_of_variable(root->name_test2);

                  if (has_var_group1) {
                    ga_instruction_set::variable_group_info
                      &vgi1 = rmi.interpolate_infos[intn1]
                                 .groups_info[root->name_test1];
                    if (has_var_group2) {
                      ga_instruction_set::variable_group_info
                        &vgi2 = rmi.interpolate_infos[intn2]
                                   .groups_info[root->name_test2];
                      pgai = std::make_shared
                             <ga_instruction_matrix_assembly_mf_mf>
                             (root->tensor(), Krr, Kru, Kur, Kuu, ctx1, ctx2,
                              vgi1, vgi2,
                              gis.coeff, gis.nbpt, gis.ipt, interpolate);
                    } else {
                      const gmm::sub_interval &I2 = mf2 && mf2->is_reduced()
                        ? workspace.temporary_interval_of_variable
                                    (root->name_test2)
                        : workspace.interval_of_variable(root->name_test2);
                      if (mf2)
                        pgai = std::make_shared
                               <ga_instruction_matrix_assembly_mf_mf>
                               (root->tensor(), Krx, Kux,  ctx1, ctx2,
                                vgi1, I2, *mf2, alpha2,
                                gis.coeff, gis.nbpt, gis.ipt, interpolate);
                      else // for global variable imd2 == 0
                        pgai = std::make_shared
                               <ga_instruction_matrix_assembly_mf_imd>
                               (root->tensor(), Krr, Kur, ctx1, ctx2,
                                vgi1, I2, imd2, alpha2, gis.coeff, gis.ipt);
                    }
                  } else { // !has_var_group1
                    const gmm::sub_interval &I1 = mf1 && mf1->is_reduced()
                      ? workspace.temporary_interval_of_variable
                                  (root->name_test1)
                      : workspace.interval_of_variable(root->name_test1);
                    if (has_var_group2) {
                      ga_instruction_set::variable_group_info
                        &vgi2 = rmi.interpolate_infos[intn2]
                                   .groups_info[root->name_test2];
                      if (mf1)
                        pgai = std::make_shared
                               <ga_instruction_matrix_assembly_mf_mf>
                               (root->tensor(), Kxr, Kxu, ctx1, ctx2,
                                I1, *mf1, alpha1, vgi2,
                                gis.coeff, gis.nbpt, gis.ipt, interpolate);
                      else // for global variable imd1 == 0
                        pgai = std::make_shared
                               <ga_instruction_matrix_assembly_imd_mf>
                               (root->tensor(), Krr, Kru, ctx1, ctx2,
                                I1, imd1, alpha1, vgi2, gis.coeff, gis.ipt);
                    } else { // !has_var_group2
                      const gmm::sub_interval &I2 = mf2 && mf2->is_reduced()
                        ? workspace.temporary_interval_of_variable
                                    (root->name_test2)
                        : workspace.interval_of_variable(root->name_test2);
                      if (mf1 && mf2)
                        pgai = std::make_shared
                               <ga_instruction_matrix_assembly_mf_mf>
                               (root->tensor(), Kxx, ctx1, ctx2,
                                I1, *mf1, alpha1, I2, *mf2, alpha2,
                                gis.coeff, gis.nbpt, gis.ipt, interpolate);
                      else if (mf1) // for global variable imd2 == 0
                        pgai = std::make_shared
                               <ga_instruction_matrix_assembly_mf_imd>
                               (root->tensor(), Kxr, ctx1, ctx2,
                                I1, *mf1, alpha1, I2, imd2, alpha2,
                                gis.coeff, gis.ipt);
                      else if (mf2)
                        pgai = std::make_shared
                               <ga_instruction_matrix_assembly_imd_mf>
                               (root->tensor(), Krx, ctx1, ctx2,
                                I1, imd1, alpha1, I2, *mf2, alpha2,
                                gis.coeff, gis.ipt);
                      else
                        pgai = std::make_shared
                               <ga_instruction_matrix_assembly_imd_imd>
                               (root->tensor(), Krr, ctx1, ctx2,
                                I1, imd1, alpha1, I2, imd2, alpha2,
                                gis.coeff, gis.ipt);
                    }
                  }
                } // if (!simple)
                break;
              } // case 2
              } // switch(order)
              if (pgai)
                rmi.instructions.push_back(std::move(pgai));
            }
          } // if (root)
        } // if (td.order == order || td.order == size_type(-1))
      } // for (const ga_workspace::tree_description &td : trees_of_current_phase)

      if (condensation && order == 2 && phase == ga_workspace::ASSEMBLY) {

        auto &Krr = workspace.assembled_matrix();
        auto &Kru = workspace.col_unreduced_matrix();
        auto &Kur = workspace.row_unreduced_matrix();
        auto &Kuu = workspace.row_col_unreduced_matrix();

        for (auto &&key_val : condensations) {
          const ga_instruction_set::region_mim rm = key_val.first;
          condensation_description &CC = key_val.second;
          auto &rmi = gis.all_instructions[rm];

          CC.KQJpr.resize(CC.KQJ.nrows(), CC.KQJ.ncols());
          for (size_type k=0; k < CC.KQJpr.size(); ++k) {
            gis.condensation_tensors.push_back // memory allocation
                                     (std::make_shared<base_tensor>(2,2));
            CC.KQJpr[k] = gis.condensation_tensors.back().get();
          }

          pga_instruction pgai;

          // Add one diagonal/subdiagonal condensation instruction per cluster
          for (size_type k=0; k < CC.Qclusters.size(); ++k) {
            // extract condensed variables residuals from
            // workspace.assembled_vector() into RQpr
            for (size_type q1 : CC.Qclusters[k]) {
              std::string name_test1 = CC.Qvars[q1];
              const im_data *imd1 = workspace.associated_im_data(name_test1);
              const gmm::sub_interval
                &I1 = workspace.interval_of_variable(name_test1);
              pgai =
                std::make_shared<ga_instruction_extract_residual_on_imd_dofs>
                (*(CC.RQpr[q1]), workspace.cached_vector(), // cached_V --> CC.RQpr[q1]
                 gis.ctx, I1, *imd1, gis.ipt);
              rmi.instructions.push_back(std::move(pgai));
            }

            // the exec() of this instruction calculates KQJpr including any
            // necessary size update to match the sizes of KQJ, upon size change
            // of primary variables J
            pgai = std::make_shared<ga_instruction_condensation_sub>
              (CC.KQJpr, CC.RQpr, CC.KQQ, CC.KQJ, CC.Qclusters[k], gis.coeff); // factor_of_variable()?
            rmi.instructions.push_back(std::move(pgai));

            // assemble/store KQJpr/RQpr matrices/vectors into the
            // corresponding global matrix/vector
            for (size_type q1 : CC.Qclusters[k]) {
              std::string name_test1 = CC.Qvars[q1];
              const im_data *imd1 = workspace.associated_im_data(name_test1);
//              const scalar_type
//                &alpha1 = workspace.factor_of_variable(name_test1); // TODO
              const gmm::sub_interval
                &I1 = workspace.interval_of_variable(name_test1);
              GMM_ASSERT1(imd1, "Internal error");
              for (size_type j2 : CC.Jclusters[k]) {
                std::string name_test2 = CC.Jvars[j2];
                const mesh_fem *mf2 = workspace.associated_mf(name_test2); // TODO: name_test2 variable group
                const im_data *imd2 = workspace.associated_im_data(name_test2);
//                const std::string &intn2 = root->interpolate_name_test2;
//                GMM_ASSERT1(intn2.empty(), "Coupling of internal variables "
//                                           "with interpolated variables not "
//                                           "implemented yet");
//                const scalar_type
//                  &alpha2 = workspace.factor_of_variable(name_test2); // TODO
                const gmm::sub_interval
                  &I2 = mf2 && mf2->is_reduced()
                      ? workspace.temporary_interval_of_variable(name_test2)
                      : workspace.interval_of_variable(name_test2);
                const base_tensor &Kq1j2pr = *(CC.KQJpr(q1,j2)); // <- input
                model_real_sparse_matrix
                  &KQJpr = mf2 && mf2->is_reduced()
                         ? workspace.col_unreduced_matrix()
                         : workspace.internal_coupling_matrix(); // <- output
                if (mf2) {
                  pgai =
                    std::make_shared<ga_instruction_matrix_assembly_imd_mf>
                    (Kq1j2pr, KQJpr, gis.ctx, gis.ctx,
                     I1, imd1, gis.ONE, I2, *mf2, gis.ONE, gis.ONE, gis.ipt); // without gis.coeff
                    // TODO: name_test2 variable group
                    if (mf2->is_reduced())
                      gis.unreduced_terms.emplace(name_test1, name_test2);
                } else // for global variable imd2 == 0
                  pgai =
                    std::make_shared<ga_instruction_matrix_assembly_imd_imd>
                    (Kq1j2pr, KQJpr, gis.ctx, gis.ctx,
                     I1, imd1, gis.ONE, I2, imd2, gis.ONE, gis.ONE, gis.ipt); // without gis.coeff
                rmi.instructions.push_back(std::move(pgai));
              } // for j2
              const bool initialize = true;
              pgai = std::make_shared<ga_instruction_vector_assembly_imd>
                (*(CC.RQpr[q1]), workspace.assembled_vector(), // <- overwriting internal variables residual with internal solution
                 gis.ctx, I1, *imd1, gis.ONE, gis.ipt, initialize); // without gis.coeff
              rmi.instructions.push_back(std::move(pgai));
            } // for q1
          }

          // Add superdiagonal condensation instructions
          for (size_type i1=0; i1 < CC.Ivars.size(); ++i1) {

            std::string name_test1 = CC.Ivars[i1];
            const mesh_fem *mf1 = workspace.associated_mf(name_test1); // TODO: name_test1 variable group
            const im_data *imd1 = workspace.associated_im_data(name_test1);
            const scalar_type
              &alpha1 = workspace.factor_of_variable(name_test1);
            const gmm::sub_interval
              &I1 = mf1 && mf1->is_reduced()
                  ? workspace.temporary_interval_of_variable(name_test1)
                  : workspace.interval_of_variable(name_test1);

            // Q_of_J[j2] will hold all condensed variables q that couple
            // variable i1 to each variable j2
            std::vector<std::set<size_type>> Q_of_J(CC.Jvars.size());
            for (size_type q=0; q < CC.Qvars.size(); ++q)
              if (CC.KIQ(i1,q)) {
                size_type cid = CC.cluster_of_Qvar[q];
                for (size_type j : CC.Jclusters[cid])
                  Q_of_J[j].insert(q);
              }

            for (size_type j2=0; j2 < CC.Jvars.size(); ++j2) {
              if (Q_of_J[j2].size()) { // a coupling between i1 and j2 exists
                std::vector<base_tensor *> Ki1Q, KQj2;
                for (size_type q : Q_of_J[j2]) {
                  Ki1Q.push_back(CC.KIQ(i1,q));
                  KQj2.push_back(CC.KQJpr(q,j2));
                }
                // allocate a tensor for storing the coupling between i1 and j2
                gis.condensation_tensors.push_back
                                         (std::make_shared<base_tensor>());
                base_tensor &Kij = *gis.condensation_tensors.back();
                pgai = std::make_shared<ga_instruction_condensation_super_K>
                       (Kij, Ki1Q, KQj2);
                rmi.instructions.push_back(std::move(pgai));
                // add assembly instruction
                std::string name_test2 = CC.Jvars[j2];
                const mesh_fem *mf2 = workspace.associated_mf(name_test2); // TODO: name_test2 variable group
                const im_data *imd2 = workspace.associated_im_data(name_test2);
                // Here assuming interpolate_name_test1.empty() &&
                //               interpolate_name_test2.empty() &&
                //               !(secondary1 || secondary2) && !interpolate;
                const scalar_type
                  &alpha2 = workspace.factor_of_variable(name_test2);
                const gmm::sub_interval
                  &I2 = mf2 && mf2->is_reduced()
                      ? workspace.temporary_interval_of_variable(name_test2)
                      : workspace.interval_of_variable(name_test2);

                auto &Kxu = (mf1 && mf1->is_reduced()) ? Kuu : Kru;
                auto &Kxr = (mf1 && mf1->is_reduced()) ? Kur : Krr;
                auto &Krx = (mf2 && mf2->is_reduced()) ? Kru : Krr;
                auto &Kxx = (mf2 && mf2->is_reduced()) ? Kxu : Kxr;

                if ((mf1 && mf1->is_reduced()) || (mf2 && mf2->is_reduced()))
                  gis.unreduced_terms.emplace(name_test1, name_test2);

                if (mf1 && mf2) // TODO: name_test1 or name_test2 variable group
                  pgai = std::make_shared
                         <ga_instruction_matrix_assembly_mf_mf>
                         (Kij, Kxx, gis.ctx, gis.ctx,
                          I1, *mf1, alpha1, I2, *mf2, alpha2,
                          gis.coeff, gis.nbpt, gis.ipt, false);
                else if (mf1) // for global variable imd2 == 0
                  pgai = std::make_shared
                         <ga_instruction_matrix_assembly_mf_imd>
                         (Kij, Kxr, gis.ctx, gis.ctx,
                          I1, *mf1, alpha1, I2, imd2, alpha2,
                          gis.coeff, gis.ipt);
                else if (mf2)
                  pgai = std::make_shared
                         <ga_instruction_matrix_assembly_imd_mf>
                         (Kij, Krx, gis.ctx, gis.ctx,
                          I1, imd1, alpha1, I2, *mf2, alpha2,
                          gis.coeff, gis.ipt);
                else
                  pgai = std::make_shared
                         <ga_instruction_matrix_assembly_imd_imd>
                         (Kij, Krr, gis.ctx, gis.ctx,
                          I1, imd1, alpha1, I2, imd2, alpha2,
                          gis.coeff, gis.ipt);
                rmi.instructions.push_back(std::move(pgai));
              } // if (Q_of_J[j2].size())
            } // for j2

            // RHS condensation instructions
            std::vector<base_tensor *> Ki1Q, RQpr;
            for (size_type q=0; q < CC.Qvars.size(); ++q)
              if (CC.KIQ(i1,q)) {
                Ki1Q.push_back(CC.KIQ(i1,q));
                RQpr.push_back(CC.RQpr[q]);
              }
            gis.condensation_tensors.push_back
                                     (std::make_shared<base_tensor>());
            base_tensor &Ri = *gis.condensation_tensors.back();
            pgai = std::make_shared<ga_instruction_condensation_super_R>
              (Ri, Ki1Q, RQpr);
            rmi.instructions.push_back(std::move(pgai));

            base_vector &R = mf1->is_reduced() ? workspace.unreduced_vector()
                                               : workspace.assembled_vector();
            if (mf1)
              pgai = std::make_shared<ga_instruction_vector_assembly_mf>
                (Ri, R, gis.ctx, I1, *mf1, gis.coeff, gis.nbpt, gis.ipt, false);
            else if (imd1)
              pgai = std::make_shared<ga_instruction_vector_assembly_imd>
                (Ri, R, gis.ctx, I1, *imd1, gis.coeff, gis.ipt);
            else
              pgai = std::make_shared<ga_instruction_vector_assembly>
                (Ri, R, I1, gis.coeff);
            rmi.instructions.push_back(std::move(pgai));
          } // for i1
        } // for (const auto &key_val : condensations)
      } // if (phase == ga_workspace::ASSEMBLY)
    } // for (const auto &phase : phases)

  } // ga_compile(...)



  //=========================================================================
  // Execution of a compiled set of assembly terms
  //=========================================================================


  void ga_function_exec(ga_instruction_set &gis) {

    for (auto &&instr : gis.all_instructions) {
      const auto &gil = instr.second.instructions;
      for (size_type j = 0; j < gil.size(); ++j) j += gil[j]->exec();
    }
  }

  void ga_interpolation_exec(ga_instruction_set &gis,
                             ga_workspace &workspace,
                             ga_interpolation_context &gic) {
    base_matrix G;
    base_small_vector un, up;

    for (const std::string &t : gis.transformations)
      workspace.interpolate_transformation(t)->init(workspace);

    for (auto &&instr : gis.all_instructions) {

      const getfem::mesh_im &mim = *(instr.first.mim());
      const mesh_region &region = *(instr.first.region());
      const getfem::mesh &m = *(instr.second.m);
      GMM_ASSERT1(&m == &(gic.linked_mesh()),
                  "Incompatibility of meshes in interpolation");
      const auto &gilb = instr.second.begin_instructions;
      const auto &gile = instr.second.elt_instructions;
      const auto &gil = instr.second.instructions;

      // iteration on elements (or faces of elements)
      std::vector<size_type> ind;
      auto pai_old = papprox_integration{};
      for (getfem::mr_visitor v(region, m, true); !v.finished(); ++v) {
        if (gic.use_mim()) {
          if (!mim.convex_index().is_in(v.cv())) continue;
          gis.pai = mim.int_method_of_element(v.cv())->approx_method();
        } else
          gis.pai = 0;

        ind.resize(0);
        bgeot::pstored_point_tab pspt
          = gic.ppoints_for_element(v.cv(), v.f(), ind);

        if (pspt.get() && ind.size() && pspt->size()) {
          m.points_of_convex(v.cv(), G);
          bgeot::pgeometric_trans pgt = m.trans_of_convex(v.cv());
          up.resize(G.nrows());
          un.resize(pgt->dim());

          if (gis.ctx.have_pgp() && gis.ctx.pgt() == pgt && pai_old == gis.pai) {
            gis.ctx.change(gis.ctx.pgp(), 0, 0, G, v.cv(), v.f());
          } else {
            if (!(gic.use_pgp(v.cv()))) {
              gis.ctx.change(pgt, 0, (*pspt)[0], G, v.cv(), v.f());
            } else {
              gis.ctx.change(gis.gp_pool(pgt, pspt), 0, 0, G, v.cv(), v.f());
            }
          }
          pai_old = gis.pai;

          if (gis.need_elt_size)
            gis.elt_size = m.convex_radius_estimate(v.cv()) * scalar_type(2);

          // iterations on interpolation points
          gis.nbpt = pspt->size();
          for (size_type ii = 0; ii < ind.size(); ++ii) {
            gis.ipt = ii;
            if (gis.ctx.have_pgp()) gis.ctx.set_ii(ind[ii]);
            else gis.ctx.set_xref((*pspt)[gis.ipt]);

            if (ii == 0 || !(pgt->is_linear())) {
              // Computation of unit normal vector in case of a boundary
              if (v.f() != short_type(-1)) {
                const base_matrix& B = gis.ctx.B();
                gmm::copy(pgt->normals()[v.f()], un);
                gmm::mult(B, un, up);
                scalar_type nup = gmm::vect_norm2(up);
                gmm::scale(up,1.0/nup);
                gmm::clean(up, 1e-13);
                gis.Normal = up;
              } else gis.Normal.resize(0);
            }
            gmm::clear(workspace.assembled_tensor().as_vector());
            if (ii == 0) {
              for (size_type j = 0; j < gilb.size(); ++j) j += gilb[j]->exec();
              for (size_type j = 0; j < gile.size(); ++j) j += gile[j]->exec();
            }
            for (size_type j = 0; j < gil.size(); ++j) j += gil[j]->exec();
            gic.store_result(v.cv(), ind[ii], workspace.assembled_tensor());
          }
        }
      }
    }
    for (const std::string &t : gis.transformations)
      workspace.interpolate_transformation(t)->finalize();

    gic.finalize();
  }

  void ga_exec(ga_instruction_set &gis, ga_workspace &workspace) {
    base_matrix G1, G2;
    base_small_vector un;
    scalar_type J1(0), J2(0);

    for (const std::string &t : gis.transformations)
      workspace.interpolate_transformation(t)->init(workspace);

    for (auto &instr : gis.all_instructions) {
      const getfem::mesh_im &mim = *(instr.first.mim());
      psecondary_domain psd = instr.first.psd();
      const getfem::mesh &m = *(instr.second.m);
      GMM_ASSERT1(&m == &(mim.linked_mesh()), "Incompatibility of meshes");
      const auto &gilb = instr.second.begin_instructions;
      const auto &gile = instr.second.elt_instructions;
      const auto &gil = instr.second.instructions;

      // if (gilb.size()) cout << "Begin instructions\n";
      // for (size_type j = 0; j < gilb.size(); ++j)
      //   cout << typeid(*(gilb[j])).name() << endl;
      // if (gile.size()) cout << "\nElement instructions\n";
      // for (size_type j = 0; j < gile.size(); ++j)
      //   cout << typeid(*(gile[j])).name() << endl;
      // cout << "\nGauss pt instructions\n";
      // for (size_type j = 0; j < gil.size(); ++j)
      //   cout << typeid(*(gil[j])).name() << endl;

      if (!psd) { // standard integration on a single domain

        const mesh_region &region = *(instr.first.region());

        // iteration on elements (or faces of elements)
        size_type old_cv = size_type(-1);
        bgeot::pgeometric_trans pgt = 0, pgt_old = 0;
        pintegration_method pim = 0;
        papprox_integration pai = 0;
        bgeot::pstored_point_tab pspt = 0, old_pspt = 0;
        bgeot::pgeotrans_precomp pgp = 0;
        bool first_gp = true;
        for (getfem::mr_visitor v(region, m, true); !v.finished(); ++v) {
          if (mim.convex_index().is_in(v.cv())) {
            // cout << "proceed with elt " << v.cv() << " face " << v.f()<<endl;
            if (v.cv() != old_cv) {
              pgt = m.trans_of_convex(v.cv());
              pim = mim.int_method_of_element(v.cv());
              m.points_of_convex(v.cv(), G1);

              if (pim->type() == IM_NONE) continue;
              GMM_ASSERT1(pim->type() == IM_APPROX, "Sorry, exact methods "
                          "cannot be used in high level generic assembly");
              pai = pim->approx_method();
              pspt = pai->pintegration_points();
              if (pspt->size()) {
                if (pgp && gis.pai == pai && pgt_old == pgt) {
                  gis.ctx.change(pgp, 0, 0, G1, v.cv(), v.f());
                } else {
                  if (pai->is_built_on_the_fly()) {
                    gis.ctx.change(pgt, 0, (*pspt)[0], G1, v.cv(), v.f());
                    pgp = 0;
                  } else {
                    pgp = gis.gp_pool(pgt, pspt);
                    gis.ctx.change(pgp, 0, 0, G1, v.cv(), v.f());
                  }
                  pgt_old = pgt; gis.pai = pai;
                }
                if (gis.need_elt_size)
                  gis.elt_size = convex_radius_estimate(pgt, G1)*scalar_type(2);
              }
              old_cv = v.cv();
            } else {
              if (pim->type() == IM_NONE) continue;
              gis.ctx.set_face_num(v.f());
            }
            if (pspt != old_pspt) { first_gp = true; old_pspt = pspt; }
            if (pspt->size()) {
              // Iterations on Gauss points
              size_type first_ind = 0;
              if (v.f() != short_type(-1)) {
                gis.nbpt = pai->nb_points_on_face(v.f());
                first_ind = pai->ind_first_point_on_face(v.f());
              } else {
                gis.nbpt = pai->nb_points_on_convex();
              }
              for (gis.ipt = 0; gis.ipt < gis.nbpt; ++(gis.ipt)) {
                if (pgp) gis.ctx.set_ii(first_ind+gis.ipt);
                else gis.ctx.set_xref((*pspt)[first_ind+gis.ipt]);
                if (gis.ipt == 0 || !(pgt->is_linear())) {
                  J1 = gis.ctx.J();
                  // Computation of unit normal vector in case of a boundary
                  if (v.f() != short_type(-1)) {
                    gis.Normal.resize(G1.nrows());
                    un.resize(pgt->dim());
                    gmm::copy(pgt->normals()[v.f()], un);
                    gmm::mult(gis.ctx.B(), un, gis.Normal);
                    scalar_type nup = gmm::vect_norm2(gis.Normal);
                    J1 *= nup;
                    gmm::scale(gis.Normal, 1.0/nup);
                    gmm::clean(gis.Normal, 1e-13);
                  } else gis.Normal.resize(0);
                }
                auto ipt_coeff = pai->coeff(first_ind+gis.ipt);
                gis.coeff = J1 * ipt_coeff;
                bool enable_ipt = (gmm::abs(ipt_coeff) > 0.0 ||
                                   workspace.include_empty_int_points());
                if (!enable_ipt) gis.coeff = scalar_type(0);
                if (first_gp) {
                  for (size_type j=0; j < gilb.size(); ++j) j+=gilb[j]->exec();
                  first_gp = false;
                }
                if (gis.ipt == 0) {
                  for (size_type j=0; j < gile.size(); ++j) j+=gile[j]->exec();
                }
                if (enable_ipt || gis.ipt == 0 || gis.ipt == gis.nbpt-1) {
                  for (size_type j=0; j < gil.size(); ++j) j+=gil[j]->exec();
                }
                GA_DEBUG_INFO("");
              }
            }
          }
        }
        GA_DEBUG_INFO("-----------------------------");

      } else { // Integration on the product of two domains (secondary domain)

        auto &sdi = instr.second.secondary_domain_infos;
        const mesh_region &region1 = *(instr.first.region());

        // iteration on elements (or faces of elements)
        size_type old_cv1=size_type(-1), old_cv2=size_type(-1);
        size_type nbpt1 = 0, nbpt2 = 0;
        bgeot::pgeometric_trans pgt1 = 0, pgt1_old = 0, pgt2 = 0, pgt2_old = 0;
        pintegration_method pim1 = 0, pim2 = 0;
        papprox_integration pai1 = 0, pai2 = 0;
        bgeot::pstored_point_tab pspt1=0, old_pspt1=0, pspt2=0, old_pspt2=0;
        bgeot::pgeotrans_precomp pgp1 = 0, pgp2 = 0;
        bool first_gp = true;
        for (getfem::mr_visitor v1(region1, m, true); !v1.finished(); ++v1) {
          if (mim.convex_index().is_in(v1.cv())) {
            // cout << "proceed with elt " << v1.cv()<<" face " << v1.f()<<endl;
            if (v1.cv() != old_cv1) {
              pgt1 = m.trans_of_convex(v1.cv());
              pim1 = mim.int_method_of_element(v1.cv());
              m.points_of_convex(v1.cv(), G1);

              if (pim1->type() == IM_NONE) continue;
              GMM_ASSERT1(pim1->type() == IM_APPROX, "Sorry, exact methods "
                          "cannot be used in high level generic assembly");
              pai1 = pim1->approx_method();
              pspt1 = pai1->pintegration_points();
              if (pspt1->size()) {
                if (pgp1 && gis.pai == pai1 && pgt1_old == pgt1) {
                  gis.ctx.change(pgp1, 0, 0, G1, v1.cv(), v1.f());
                } else {
                  if (pai1->is_built_on_the_fly()) {
                    gis.ctx.change(pgt1, 0, (*pspt1)[0], G1, v1.cv(), v1.f());
                    pgp1 = 0;
                  } else {
                    pgp1 = gis.gp_pool(pgt1, pspt1);
                    gis.ctx.change(pgp1, 0, 0, G1, v1.cv(), v1.f());
                  }
                  pgt1_old = pgt1; gis.pai = pai1;
                }
                if (gis.need_elt_size)
                  gis.elt_size = convex_radius_estimate(pgt1,G1)*scalar_type(2);
              }
              old_cv1 = v1.cv();
            } else {
              if (pim1->type() == IM_NONE) continue;
              gis.ctx.set_face_num(v1.f());
            }
            if (pspt1 != old_pspt1) { first_gp = true; old_pspt1 = pspt1; }
            if (pspt1->size()) {
              // iterations on Gauss points
              size_type first_ind1 = 0;
              if (v1.f() != short_type(-1)) {
                nbpt1 = pai1->nb_points_on_face(v1.f());
                first_ind1 = pai1->ind_first_point_on_face(v1.f());
              } else {
                nbpt1 = pai1->nb_points_on_convex();
              }

              const mesh &m2 = psd->mim().linked_mesh();
              const mesh_region &region2 = psd->give_region(m, v1.cv(), v1.f());
              for (getfem::mr_visitor v2(region2, m2, true);
                   !v2.finished(); ++v2) {
                if (v2.cv() != old_cv2) {
                  pgt2 = m2.trans_of_convex(v2.cv());
                  pim2 = psd->mim().int_method_of_element(v2.cv());
                  m2.points_of_convex(v2.cv(), G2);

                  if (pim2->type() == IM_NONE) continue;
                  GMM_ASSERT1(pim2->type() == IM_APPROX, "Sorry, exact methods "
                              "cannot be used in high level generic assembly");
                  pai2 = pim2->approx_method();
                  pspt2 = pai2->pintegration_points();
                  if (pspt2->size()) {
                    if (pgp2 && sdi.pai == pai2 && pgt2_old == pgt2) {
                      sdi.ctx.change(pgp2, 0, 0, G2, v2.cv(), v2.f());
                    } else {
                      if (pai2->is_built_on_the_fly()) {
                        sdi.ctx.change(pgt2, 0, (*pspt2)[0], G2,v2.cv(),v2.f());
                        pgp2 = 0;
                      } else {
                        pgp2 = gis.gp_pool(pgt2, pspt2);
                        sdi.ctx.change(pgp2, 0, 0, G2, v2.cv(), v2.f());
                      }
                      pgt2_old = pgt2; sdi.pai = pai2;
                    }
                  }
                  old_cv2 = v2.cv();
                } else {
                  if (pim2->type() == IM_NONE) continue;
                  sdi.ctx.set_face_num(v2.f());
                }
                if (pspt2 != old_pspt2) { first_gp = true; old_pspt2 = pspt2; }
                if (pspt2->size()) {
                  // iterations on Gauss points
                  size_type first_ind2 = 0;
                  if (v2.f() != short_type(-1)) {
                    nbpt2 = pai2->nb_points_on_face(v2.f());
                    first_ind2 = pai2->ind_first_point_on_face(v2.f());
                  } else {
                    nbpt2 = gis.nbpt = pai2->nb_points_on_convex();
                  }
                  gis.nbpt = nbpt1 * nbpt2;
                  gis.ipt = 0;
                  for (size_type ipt1=0; ipt1 < nbpt1; ++ipt1) {
                    for (size_type ipt2=0; ipt2 < nbpt2; ++ipt2, ++(gis.ipt)) {

                      if (pgp1) gis.ctx.set_ii(first_ind1+ipt1);
                      else gis.ctx.set_xref((*pspt1)[first_ind1+ipt1]);
                      if (pgp2) sdi.ctx.set_ii(first_ind2+ipt2);
                      else sdi.ctx.set_xref((*pspt2)[first_ind2+ipt2]);

                      if (gis.ipt == 0 || !(pgt1->is_linear())) {
                        J1 = gis.ctx.J();
                        if (v1.f() != short_type(-1)) {
                          gis.Normal.resize(G1.nrows());
                          un.resize(pgt1->dim());
                          gmm::copy(pgt1->normals()[v1.f()], un);
                          gmm::mult(gis.ctx.B(), un, gis.Normal);
                          scalar_type nup = gmm::vect_norm2(gis.Normal);
                          J1 *= nup;
                          gmm::scale(gis.Normal, 1.0/nup);
                          gmm::clean(gis.Normal, 1e-13);
                        } else gis.Normal.resize(0);
                      }

                      if (gis.ipt == 0 || !(pgt2->is_linear())) {
                        J2 = sdi.ctx.J();
                        if (v2.f() != short_type(-1)) {
                          sdi.Normal.resize(G2.nrows());
                          un.resize(pgt2->dim());
                          gmm::copy(pgt2->normals()[v2.f()], un);
                          gmm::mult(sdi.ctx.B(), un, sdi.Normal);
                          scalar_type nup = gmm::vect_norm2(sdi.Normal);
                          J2 *= nup;
                          gmm::scale(sdi.Normal, 1.0/nup);
                          gmm::clean(sdi.Normal, 1e-13);
                        } else sdi.Normal.resize(0);
                      }

                      auto ipt_coeff = pai1->coeff(first_ind1+ipt1)
                                     * pai2->coeff(first_ind2+ipt2);
                      gis.coeff = J1 * J2 * ipt_coeff;
                      bool enable_ipt = (gmm::abs(ipt_coeff) > 0.0 ||
                                         workspace.include_empty_int_points());
                      if (!enable_ipt) gis.coeff = scalar_type(0);

                      if (first_gp) {
                        for (size_type j=0; j < gilb.size(); ++j)
                          j+=gilb[j]->exec();
                        first_gp = false;
                      }
                      if (gis.ipt == 0) {
                        for (size_type j=0; j < gile.size(); ++j)
                          j+=gile[j]->exec();
                      }
                      if (enable_ipt || gis.ipt == 0 || gis.ipt == gis.nbpt-1) {
                        for (size_type j=0; j < gil.size(); ++j)
                          j+=gil[j]->exec();
                      }
                      GA_DEBUG_INFO("");
                    }
                  }
                }
              }
            }
          }
        }
        GA_DEBUG_INFO("-----------------------------");
      }

    }

    for (const std::string &t : gis.transformations)
      workspace.interpolate_transformation(t)->finalize();
  }


} /* end of namespace */
