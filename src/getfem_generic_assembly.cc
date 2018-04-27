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

namespace getfem {

  extern bool predef_operators_nonlinear_elasticity_initialized;
  extern bool predef_operators_plasticity_initialized;
  extern bool predef_operators_contact_initialized;

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
        cv_old(-1) {}
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
      GA_DEBUG_INFO("Instruction: Normal");
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

  struct ga_instruction_elementary_transformation {
    const base_vector &coeff_in;
    base_vector coeff_out;
    pelementary_transformation elemtrans;
    const mesh_fem &mf;
    const fem_interpolation_context &ctx;
    base_matrix &M;
    const mesh_fem **mf_M;
    size_type &icv;

    void do_transformation() {
      size_type nn = gmm::vect_size(coeff_in);
      if (M.size() == 0 || icv != ctx.convex_num() || &mf != *mf_M) {
        M.base_resize(nn, nn);
        *mf_M = &mf; icv = ctx.convex_num();
        elemtrans->give_transformation(mf, icv, M);
      }
      coeff_out.resize(nn);
      gmm::mult(M, coeff_in, coeff_out); // remember: coeff == coeff_out
    }

    ga_instruction_elementary_transformation
    (const base_vector &co, pelementary_transformation e,
     const mesh_fem &mf_, const fem_interpolation_context &ctx_,
     base_matrix &M_, const mesh_fem **mf_M_, size_type &icv_)
      : coeff_in(co), elemtrans(e), mf(mf_), ctx(ctx_),
        M(M_), mf_M(mf_M_), icv(icv_) {}
    ~ga_instruction_elementary_transformation() {};
  };

  struct ga_instruction_elementary_transformation_val
    : public ga_instruction_val, ga_instruction_elementary_transformation {
    // Z(ndof,target_dim), coeff_in(Qmult,ndof) --> t(target_dim*Qmult)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: variable value with elementary "
                    "transformation");
      do_transformation();
      return ga_instruction_val::exec();
    }

    ga_instruction_elementary_transformation_val
    (base_tensor &tt, const base_tensor &Z_, const base_vector &co, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_val(tt, Z_, coeff_out, q),
        ga_instruction_elementary_transformation(co, e, mf_, ctx_, M_,
                                                 mf_M_, icv_) {}
  };

  struct ga_instruction_elementary_transformation_grad
    : public ga_instruction_grad, ga_instruction_elementary_transformation {
    // Z(ndof,target_dim,N), coeff_in(Qmult,ndof) --> t(target_dim*Qmult,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: gradient with elementary transformation");
      do_transformation();
      return ga_instruction_grad::exec();
    }

    ga_instruction_elementary_transformation_grad
    (base_tensor &tt, const base_tensor &Z_, const base_vector &co, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_grad(tt, Z_, coeff_out, q),
        ga_instruction_elementary_transformation(co, e, mf_, ctx_, M_,
                                                 mf_M_, icv_) {}
  };

  struct ga_instruction_elementary_transformation_hess
    : public ga_instruction_hess, ga_instruction_elementary_transformation {
    // Z(ndof,target_dim,N,N), coeff_in(Qmult,ndof) --> t(target_dim*Qmult,N,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Hessian with elementary transformation");
      do_transformation();
      return ga_instruction_hess::exec();
    }

    ga_instruction_elementary_transformation_hess
    (base_tensor &tt, const base_tensor &Z_, const base_vector &co, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_hess(tt, Z_, coeff_out, q),
        ga_instruction_elementary_transformation(co, e, mf_, ctx_, M_,
                                                 mf_M_, icv_) {}
  };

  struct ga_instruction_elementary_transformation_diverg
    : public ga_instruction_diverg, ga_instruction_elementary_transformation {
    // Z(ndof,target_dim,N), coeff_in(Qmult,ndof) --> t(1)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: divergence with elementary transformation");
      do_transformation();
      return ga_instruction_diverg::exec();
    }

    ga_instruction_elementary_transformation_diverg
    (base_tensor &tt, const base_tensor &Z_, const base_vector &co, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_diverg(tt, Z_, coeff_out, q),
        ga_instruction_elementary_transformation(co, e, mf_, ctx_, M_,
                                                 mf_M_, icv_) {}
  };

  struct ga_instruction_update_group_info : public ga_instruction {
    const ga_workspace &workspace;
    ga_instruction_set &gis;
    ga_instruction_set::interpolate_info &inin;
    const std::string gname;
    ga_instruction_set::variable_group_info &vgi;

    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Update group info for "+gname);
      if (vgi.varname &&
          &(workspace.associated_mf(*(vgi.varname))->linked_mesh())==inin.m)
        return 0;
      const std::string &varname
        = inin.m ? workspace.variable_in_group(gname, *(inin.m))
                 : workspace.first_variable_of_group(gname);
      vgi.mf = workspace.associated_mf(varname);
      vgi.Ir = gis.var_intervals[varname];
      vgi.In = workspace.interval_of_variable(varname);
      vgi.alpha = workspace.factor_of_variable(varname);
      vgi.U = gis.extended_vars[varname];
      vgi.varname = &varname;
      return 0;
    }

    ga_instruction_update_group_info
    (const ga_workspace &workspace_, ga_instruction_set &gis_,
     ga_instruction_set::interpolate_info &inin_, const std::string &gname_,
     ga_instruction_set::variable_group_info &vgi_) :
      workspace(workspace_), gis(gis_), inin(inin_), gname(gname_),
      vgi(vgi_) {}
  };

  struct ga_instruction_interpolate_filter : public ga_instruction {
    base_tensor &t;
    const ga_instruction_set::interpolate_info &inin;
    size_type pt_type;
    int nb;

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


  struct ga_instruction_elementary_transformation_base {
    base_tensor t_in;
    base_tensor &t_out;
    pelementary_transformation elemtrans;
    const mesh_fem &mf;
    const fem_interpolation_context &ctx;
    base_matrix &M;
    const mesh_fem **mf_M;
    size_type &icv;

    void do_transformation(size_type n) {
      if (M.size() == 0 || icv != ctx.convex_num() || &mf != *mf_M) {
        M.base_resize(n, n);
        *mf_M = &mf; icv = ctx.convex_num();
        elemtrans->give_transformation(mf, icv, M);
      }
      t_out.mat_reduction(t_in, M, 0);
    }

    ga_instruction_elementary_transformation_base
    (base_tensor &t_, pelementary_transformation e, const mesh_fem &mf_,
     const fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : t_out(t_), elemtrans(e), mf(mf_), ctx(ctx_),
        M(M_), mf_M(mf_M_), icv(icv_) {}
  };

  struct ga_instruction_elementary_transformation_val_base
    : public ga_instruction_copy_val_base,
             ga_instruction_elementary_transformation_base {
    // Z(ndof,target_dim) --> t_in --> t_out(Qmult*ndof,Qmult*target_dim)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: value of test functions with elementary "
                    "transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      t_in.adjust_sizes(t_out.sizes());
      ga_instruction_copy_val_base::exec();
      do_transformation(ndof*Qmult);
      return 0;
    }

    ga_instruction_elementary_transformation_val_base
    (base_tensor &t_, const base_tensor &Z_, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_copy_val_base(t_in, Z_, q),
        ga_instruction_elementary_transformation_base(t_, e, mf_, ctx_, M_,
                                                      mf_M_, icv_) {}
  };

  struct ga_instruction_elementary_transformation_grad_base
    : public ga_instruction_copy_grad_base,
             ga_instruction_elementary_transformation_base {
    // Z(ndof,target_dim,N) --> t_in --> t_out(Qmult*ndof,Qmult*target_dim,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: gradient of test functions with elementary "
                    "transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      t_in.adjust_sizes(t_out.sizes());
      ga_instruction_copy_grad_base::exec();
      do_transformation(ndof*Qmult);
      return 0;
    }

    ga_instruction_elementary_transformation_grad_base
    (base_tensor &t_, const base_tensor &Z_, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_copy_grad_base(t_in, Z_, q),
        ga_instruction_elementary_transformation_base(t_, e, mf_, ctx_, M_,
                                                      mf_M_, icv_) {}
  };

  struct ga_instruction_elementary_transformation_hess_base
    : public ga_instruction_copy_hess_base,
             ga_instruction_elementary_transformation_base {
    // Z(ndof,target_dim,N*N) --> t_out(Qmult*ndof,Qmult*target_dim,N,N)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Hessian of test functions with elementary "
                    "transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      t_in.adjust_sizes(t_out.sizes());
      ga_instruction_copy_hess_base::exec();
      do_transformation(ndof*Qmult);
      return 0;
    }

    ga_instruction_elementary_transformation_hess_base
    (base_tensor &t_, const base_tensor &Z_, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_copy_hess_base(t_in, Z_, q),
        ga_instruction_elementary_transformation_base(t_, e, mf_, ctx_, M_,
                                                      mf_M_, icv_) {}
  };

  struct ga_instruction_elementary_transformation_diverg_base
    : public ga_instruction_copy_diverg_base,
             ga_instruction_elementary_transformation_base {
    // Z(ndof,target_dim,N) --> t_out(Qmult*ndof)
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: divergence of test functions with elementary "
                    "transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      t_in.adjust_sizes(t_out.sizes());
      ga_instruction_copy_diverg_base::exec();
      do_transformation(ndof*Qmult);
      return 0;
    }

    ga_instruction_elementary_transformation_diverg_base
    (base_tensor &t_, const base_tensor &Z_, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_copy_diverg_base(t_in, Z_, q),
        ga_instruction_elementary_transformation_base(t_, e, mf_, ctx_, M_,
                                                      mf_M_, icv_) {}
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

  struct ga_instruction_transpose : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: transpose");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      size_type order = t.sizes().size();
      size_type s1 = t.sizes()[order-2], s2 = t.sizes()[order-1];
      size_type s = t.size() / (s1*s2);
      for (size_type i = 0; i < s1;  ++i)
        for (size_type j = 0; j < s2;  ++j) {
          base_tensor::iterator it = t.begin() + s*(i + s1*j);
          base_tensor::const_iterator it1 = tc1.begin() + s*(j + s2*i);
          for (size_type k = 0; k < s; ++k) *it++ = *it1++;
        }
      return 0;
    }
    ga_instruction_transpose(base_tensor &t_, const base_tensor &tc1_)
      : t(t_), tc1(tc1_) {}
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
      GA_DEBUG_INFO("Instruction: specific componentwise "
                           "multiplication");
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

  // Performs Amij Bjk -> Cmik. To be optimized
  struct ga_instruction_matrix_mult : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: matrix multiplication");
      size_type order = tc2.sizes().size();
      size_type s2_1 = tc2.sizes()[order-2];
      size_type s2_2 = tc2.sizes()[order-1];
      size_type s1 = tc1.size() / s2_1;
      size_type s2 = tc2.size() / (s2_1*s2_2);

      base_tensor::iterator it = t.begin();
      for (size_type k = 0; k < s2_2; ++k)
        for (size_type i = 0; i < s1; ++i)
          for (size_type m = 0; m < s2; ++m, ++it) {
            *it = scalar_type(0);
            for (size_type j = 0; j < s2_1; ++j)
              *it += tc1[i+j*s1] * tc2[m+j*s2+k*s2_1*s2];
          }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_matrix_mult(base_tensor &t_, base_tensor &tc1_,
                               base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  // Performs Amij Bnjk -> Cmnik. To be optimized
  struct ga_instruction_matrix_mult_spec : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: specific matrix multiplication");
      size_type s1_1 = tc1.sizes()[0];
      size_type s1_2 = tc1.sizes()[1];
      // size_type s1_3 = tc1.sizes()[2];
      size_type s2_1 = tc2.sizes()[0];
      size_type s2_2 = tc2.sizes()[1];
      size_type s2_3 = tc2.sizes()[2];

      base_tensor::iterator it = t.begin();
      for (size_type k = 0; k < s2_3; ++k)
        for (size_type i = 0; i < s1_2; ++i)
          for (size_type n = 0; n < s2_1; ++n)
            for (size_type m = 0; m < s1_1; ++m, ++it) {
              *it = scalar_type(0);
              for (size_type j = 0; j < s2_2; ++j)
                *it += tc1[m+i*s1_1+j*s1_1*s1_2] * tc2[n+j*s2_1+k*s2_1*s2_2];
            }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_matrix_mult_spec(base_tensor &t_, base_tensor &tc1_,
                               base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };


  // Performs Ani Bmi -> Cmn
  struct ga_instruction_reduction : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: reduction operation of size " << nn);
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
    ga_instruction_reduction(base_tensor &t_, base_tensor &tc1_,
                             base_tensor &tc2_, size_type n_)
      : t(t_), tc1(tc1_), tc2(tc2_), nn(n_) {}
  };

  // Performs Ani Bmi -> Cmn
  struct ga_instruction_reduction_opt0_2 : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type n, q;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: reduction operation of size " << n*q <<
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
      // ga_instruction_reduction toto(t, tc1, tc2, n*q);
      // toto.exec();
      // GMM_ASSERT1(gmm::vect_dist2(t.as_vector(), u.as_vector()) < 1E-9, "Erroneous");
      return 0;
    }
    ga_instruction_reduction_opt0_2(base_tensor &t_, base_tensor &tc1_,
                                    base_tensor &tc2_, size_type n_,
                                    size_type q_)
      : t(t_), tc1(tc1_), tc2(tc2_), n(n_), q(q_) {}
  };

  // Performs Ani Bmi -> Cmn
  template <int N>
  struct ga_instruction_reduction_opt0_2_unrolled : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type q;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: unrolled reduction operation of size " << N*q
                    << " optimized for vectorized second tensor of type 2");
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
    ga_instruction_reduction_opt0_2_unrolled(base_tensor &t_, base_tensor &tc1_,
                                             base_tensor &tc2_, size_type q_)
      : t(t_), tc1(tc1_), tc2(tc2_), q(q_) {}
  };

  // Performs Ani Bmi -> Cmn
  template <int N, int Q>
  struct ga_instruction_reduction_opt0_2_dunrolled : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: unrolled reduction operation of size " << N*Q
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
    ga_instruction_reduction_opt0_2_dunrolled
    (base_tensor &t_, base_tensor &tc1_, base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  // Performs Ani Bmi -> Cmn
  struct ga_instruction_reduction_opt2_0 : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type n, q;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: reduction operation of size " << n*q <<
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
    ga_instruction_reduction_opt2_0(base_tensor &t_, base_tensor &tc1_,
                                    base_tensor &tc2_, size_type n_,
                                    size_type q_)
      : t(t_), tc1(tc1_), tc2(tc2_), n(n_), q(q_) { }
  };

  // Performs Ani Bmi -> Cmn
  template <int N>
  struct ga_instruction_reduction_opt2_0_unrolled : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type q;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: unrolled reduction operation of size " << N*q
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
    ga_instruction_reduction_opt2_0_unrolled(base_tensor &t_, base_tensor &tc1_,
                                             base_tensor &tc2_, size_type q_)
      : t(t_), tc1(tc1_), tc2(tc2_), q(q_) {}
  };

  // Performs Ani Bmi -> Cmn
  template <int N, int Q>
  struct ga_instruction_reduction_opt2_0_dunrolled : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: unrolled reduction operation of size " << N*Q
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
    ga_instruction_reduction_opt2_0_dunrolled
    (base_tensor &t_, base_tensor &tc1_, base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  // Performs Ani Bmi -> Cmn
  struct ga_instruction_reduction_opt0_1 : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: reduction operation of size " << nn <<
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
    ga_instruction_reduction_opt0_1(base_tensor &t_, base_tensor &tc1_,
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
  struct ga_instruction_reduction_opt0_1_unrolled : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: unrolled reduction operation of size " << N
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
    ga_instruction_reduction_opt0_1_unrolled(base_tensor &t_, base_tensor &tc1_,
                                             base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  // Performs Ani Bmi -> Cmn
  struct ga_instruction_reduction_opt1_1 : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: reduction operation of size " << nn <<
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
    ga_instruction_reduction_opt1_1(base_tensor &t_, base_tensor &tc1_,
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
  template<int N> struct ga_instruction_reduction_unrolled
    : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: unrolled reduction operation of size " << N);
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
    ga_instruction_reduction_unrolled(base_tensor &t_, base_tensor &tc1_,
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
      GA_DEBUG_INFO("Instruction: doubly unrolled reduction operation of size "
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


  pga_instruction ga_instruction_reduction_switch
  (assembly_tensor &t_, assembly_tensor &tc1_, assembly_tensor &tc2_,
   size_type n, bool &to_clear) {
    base_tensor &t = t_.tensor(), &tc1 = tc1_.tensor(), &tc2 = tc2_.tensor();

    if (tc1_.sparsity() == 1 && tc2_.sparsity() == 1 &&
        tc1_.qdim() == n && tc2_.qdim() == n) {
      to_clear = true;
      t_.set_sparsity(10, tc1_.qdim());
      return std::make_shared<ga_instruction_reduction_opt1_1>(t, tc1, tc2, n);
    }

    if (tc2_.sparsity() == 1) {
      switch(n) {
      case 2:
        return std::make_shared<ga_instruction_reduction_opt0_1_unrolled<2>>
          (t, tc1, tc2);
      case 3:
        return std::make_shared<ga_instruction_reduction_opt0_1_unrolled<3>>
          (t, tc1, tc2);
      case 4:
        return std::make_shared<ga_instruction_reduction_opt0_1_unrolled<4>>
          (t, tc1, tc2);
      case 5:
        return std::make_shared<ga_instruction_reduction_opt0_1_unrolled<5>>
          (t, tc1, tc2);
      default:
        return std::make_shared<ga_instruction_reduction_opt0_1>(t,tc1,tc2, n);
      }
    }
        if (tc2_.sparsity() == 2) {
      size_type q2 = tc2.sizes()[1];
      size_type n2 = (tc2.sizes().size() > 2) ? tc2.sizes()[1] : 1;
      if (n2*q2 == n) {
        switch (n2) {
        case 1:
          switch (q2) {
          case 2:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<1,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<1,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<1,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_reduction_opt0_2_unrolled<1>>
              (t, tc1, tc2, q2);
          }
        case 2:
          switch (q2) {
          case 2:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<2,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<2,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<2,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_reduction_opt0_2_unrolled<2>>
              (t, tc1, tc2, q2);
          }
        case 3:
          switch (q2) {
          case 2:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<3,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<3,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<3,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_reduction_opt0_2_unrolled<3>>
              (t, tc1, tc2, q2);
          }
        case 4:
          return std::make_shared<ga_instruction_reduction_opt0_2_unrolled<4>>
            (t, tc1, tc2, q2);
        case 5:
          return std::make_shared<ga_instruction_reduction_opt0_2_unrolled<5>>
            (t, tc1, tc2, q2);
        default:
          return std::make_shared<ga_instruction_reduction_opt0_2>
            (t,tc1,tc2,n2,q2);
        }
      }
    }
    if (tc1_.sparsity() == 2) {
      size_type q1 = tc1.sizes()[1];
      size_type n1 = (tc1.sizes().size() > 2) ? tc1.sizes()[1] : 1;
      if (n1*q1 == n) {
        switch (n1) {
        case 1:
          switch (q1) {
          case 2:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<1,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<1,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<1,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_reduction_opt2_0_unrolled<1>>
              (t, tc1, tc2, q1);
          }
        case 2:
          switch (q1) {
          case 2:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<2,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<2,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<2,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_reduction_opt2_0_unrolled<2>>
              (t, tc1, tc2, q1);
          }
        case 3:
          switch (q1) {
          case 2:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<3,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<3,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<3,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_reduction_opt2_0_unrolled<3>>
              (t, tc1, tc2, q1);
          }
          return std::make_shared<ga_instruction_reduction_opt2_0_unrolled<3>>
            (t, tc1, tc2, q1);
        case 4:
          return std::make_shared<ga_instruction_reduction_opt2_0_unrolled<4>>
            (t, tc1, tc2, q1);
        case 5:
          return std::make_shared<ga_instruction_reduction_opt2_0_unrolled<5>>
            (t, tc1, tc2, q1);
        default:
          return std::make_shared<ga_instruction_reduction_opt2_0>
            (t,tc1,tc2, n1, q1);
        }
      }
    }

    switch(n) {
    case  2 : return std::make_shared<ga_instruction_reduction_unrolled< 2>>
                     (t, tc1, tc2);
    case  3 : return std::make_shared<ga_instruction_reduction_unrolled< 3>>
                     (t, tc1, tc2);
    case  4 : return std::make_shared<ga_instruction_reduction_unrolled< 4>>
                     (t, tc1, tc2);
    case  5 : return std::make_shared<ga_instruction_reduction_unrolled< 5>>
                     (t, tc1, tc2);
    case  6 : return std::make_shared<ga_instruction_reduction_unrolled< 6>>
                     (t, tc1, tc2);
    case  7 : return std::make_shared<ga_instruction_reduction_unrolled< 7>>
                     (t, tc1, tc2);
    case  8 : return std::make_shared<ga_instruction_reduction_unrolled< 8>>
                     (t, tc1, tc2);
    case  9 : return std::make_shared<ga_instruction_reduction_unrolled< 9>>
                     (t, tc1, tc2);
    case 10 : return std::make_shared<ga_instruction_reduction_unrolled<10>>
                     (t, tc1, tc2);
    case 11 : return std::make_shared<ga_instruction_reduction_unrolled<11>>
                     (t, tc1, tc2);
    case 12 : return std::make_shared<ga_instruction_reduction_unrolled<12>>
                     (t, tc1, tc2);
    case 13 : return std::make_shared<ga_instruction_reduction_unrolled<13>>
                     (t, tc1, tc2);
    case 14 : return std::make_shared<ga_instruction_reduction_unrolled<14>>
                     (t, tc1, tc2);
    case 15 : return std::make_shared<ga_instruction_reduction_unrolled<15>>
                     (t, tc1, tc2);
    case 16 : return std::make_shared<ga_instruction_reduction_unrolled<16>>
                     (t, tc1, tc2);
    default : return std::make_shared<ga_instruction_reduction>
                     (t, tc1, tc2, n);
    }
  }

  pga_instruction  ga_uniform_instruction_reduction_switch
  (assembly_tensor &t_, assembly_tensor &tc1_, assembly_tensor &tc2_,
   size_type n, bool &to_clear) {
    base_tensor &t = t_.tensor(), &tc1 = tc1_.tensor(), &tc2 = tc2_.tensor();

    if (tc1_.sparsity() == 1 && tc2_.sparsity() == 1 &&
        tc1_.qdim() == n && tc2_.qdim() == n) {
      to_clear = true;
      t_.set_sparsity(10, tc1_.qdim());
      return std::make_shared<ga_instruction_reduction_opt1_1>(t,tc1,tc2,n);
    }
    if (tc2_.sparsity() == 1) {
      switch(n) {
      case 2:
        return std::make_shared<ga_instruction_reduction_opt0_1_unrolled<2>>
          (t, tc1, tc2);
      case 3:
        return std::make_shared<ga_instruction_reduction_opt0_1_unrolled<3>>
          (t, tc1, tc2);
      case 4:
        return std::make_shared<ga_instruction_reduction_opt0_1_unrolled<4>>
          (t, tc1, tc2);
      case 5:
        return std::make_shared<ga_instruction_reduction_opt0_1_unrolled<5>>
          (t, tc1, tc2);
      default:
        return std::make_shared<ga_instruction_reduction_opt0_1>(t,tc1,tc2, n);
      }
    }
    if (tc2_.sparsity() == 2) {
      size_type q2 = tc2.sizes()[1];
      size_type n2 = (tc2.sizes().size() > 2) ? tc2.sizes()[1] : 1;
      if (n2*q2 == n) {
        switch (n2) {
        case 1:
          switch (q2) {
          case 2:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<1,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<1,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<1,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_reduction_opt0_2_unrolled<1>>
              (t, tc1, tc2, q2);
          }
        case 2:
          switch (q2) {
          case 2:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<2,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<2,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<2,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_reduction_opt0_2_unrolled<2>>
              (t, tc1, tc2, q2);
          }
        case 3:
          switch (q2) {
          case 2:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<3,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<3,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_reduction_opt0_2_dunrolled<3,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_reduction_opt0_2_unrolled<3>>
              (t, tc1, tc2, q2);
          }
        case 4:
          return std::make_shared<ga_instruction_reduction_opt0_2_unrolled<4>>
            (t, tc1, tc2, q2);
        case 5:
          return std::make_shared<ga_instruction_reduction_opt0_2_unrolled<5>>
            (t, tc1, tc2, q2);
        default:
          return std::make_shared<ga_instruction_reduction_opt0_2>
            (t,tc1,tc2,n2,q2);
        }
      }
    }
    if (tc1_.sparsity() == 2) {
      size_type q1 = tc1.sizes()[1];
      size_type n1 = (tc1.sizes().size() > 2) ? tc1.sizes()[1] : 1;
      if (n1*q1 == n) {
        switch (n1) {
        case 1:
          switch (q1) {
          case 2:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<1,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<1,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<1,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_reduction_opt2_0_unrolled<1>>
              (t, tc1, tc2, q1);
          }
        case 2:
          switch (q1) {
          case 2:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<2,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<2,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<2,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_reduction_opt2_0_unrolled<2>>
              (t, tc1, tc2, q1);
          }
        case 3:
          switch (q1) {
          case 2:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<3,2>>
              (t, tc1, tc2);
          case 3:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<3,3>>
              (t, tc1, tc2);
          case 4:
            return
              std::make_shared<ga_instruction_reduction_opt2_0_dunrolled<3,4>>
              (t, tc1, tc2);
          default :
            return std::make_shared<ga_instruction_reduction_opt2_0_unrolled<3>>
              (t, tc1, tc2, q1);
          }
          return std::make_shared<ga_instruction_reduction_opt2_0_unrolled<3>>
            (t, tc1, tc2, q1);
        case 4:
          return std::make_shared<ga_instruction_reduction_opt2_0_unrolled<4>>
            (t, tc1, tc2, q1);
        case 5:
          return std::make_shared<ga_instruction_reduction_opt2_0_unrolled<5>>
            (t, tc1, tc2, q1);
        default:
          return std::make_shared<ga_instruction_reduction_opt2_0>
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
      default: return ga_instruction_reduction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 2 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,2>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,2>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,2>>(t, tc1, tc2);
      default: return ga_instruction_reduction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 3 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,3>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,3>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,3>>(t, tc1, tc2);
      default: return ga_instruction_reduction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 4 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,4>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,4>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,4>>(t, tc1, tc2);
      default: return ga_instruction_reduction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 5 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,5>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,5>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,5>>(t, tc1, tc2);
      default: return ga_instruction_reduction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 6 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,6>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,6>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,6>>(t, tc1, tc2);
      default: return ga_instruction_reduction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 7 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,7>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,7>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,7>>(t, tc1, tc2);
      default: return ga_instruction_reduction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 8 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,8>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,8>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,8>>(t, tc1, tc2);
      default: return ga_instruction_reduction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 9 :
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,9>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,9>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,9>>(t, tc1, tc2);
      default: return ga_instruction_reduction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    case 10:
      switch(n) {
      case 2: return std::make_shared<ga_ins_red_d_unrolled<2,10>>(t, tc1, tc2);
      case 3: return std::make_shared<ga_ins_red_d_unrolled<3,10>>(t, tc1, tc2);
      case 4: return std::make_shared<ga_ins_red_d_unrolled<4,10>>(t, tc1, tc2);
      default: return ga_instruction_reduction_switch(t_,tc1_,tc2_,n,to_clear);
      }
    default: return ga_instruction_reduction_switch(t_,tc1_,tc2_,n,to_clear);
    }
  }


  // Performs Amij Bnj -> Cmni. To be optimized.
  struct ga_instruction_spec_reduction : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: specific reduction operation of "
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
    ga_instruction_spec_reduction(base_tensor &t_, base_tensor &tc1_,
                                  base_tensor &tc2_, size_type n_)
      : t(t_), tc1(tc1_), tc2(tc2_), nn(n_) {}
  };

  // Performs Amik Bnjk -> Cmnij. To be optimized.
  struct ga_instruction_spec2_reduction : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: second specific reduction operation of "
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
    ga_instruction_spec2_reduction(base_tensor &t_, base_tensor &tc1_,
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
      GA_DEBUG_INFO("Instruction: simple tensor product, unrolled with "
                    << tc1.size() << " operations");
      GA_DEBUG_ASSERT(t.size() == tc1.size() * s2, "Wrong sizes");
      GA_DEBUG_ASSERT(tc1.size() == S1, "Wrong sizes");

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
    (base_tensor &t_, const std::vector<const base_tensor *>  &components_)
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
      gmm::clear(inin.Normal);
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

  struct ga_instruction_neighbour_transformation_call : public ga_instruction {
    const ga_workspace &workspace;
    ga_instruction_set::interpolate_info &inin;
    pinterpolate_transformation trans;
    fem_interpolation_context &ctx;
    base_small_vector &Normal;
    const mesh &m;
    size_type &ipt;
    papprox_integration &pai;
    bgeot::geotrans_precomp_pool &gp_pool;
    std::map<gauss_pt_corresp, bgeot::pstored_point_tab> &neighbour_corresp;

    virtual int exec() {
      bool cancel_optimization = false;
      GA_DEBUG_INFO("Instruction: call interpolate transformation");
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
            auto itm = neighbour_corresp.find(gpc);
            if (itm != neighbour_corresp.end()) {
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
                bool is_in = gic.invert(ctx_x.xreal(), P_ref[i],converged,1E-4);
                GMM_ASSERT1(is_in && converged,"Geometric transformation "
                            "inversion has failed in neighbour transformation");
              }
              pspt = store_point_tab(P_ref);
              neighbour_corresp[gpc] = pspt;
            }
            m.points_of_convex(adj_face.cv, inin.G);
            bgeot::pgeotrans_precomp pgp = gp_pool(gpc.pgt2, pspt);
            inin.ctx.change(pgp, 0, 0, inin.G, adj_face.cv, adj_face.f);
          }
        }
      }

      if (inin.ctx.have_pgp()) {
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
        inin.pt_type = trans->transform(workspace, m, ctx, Normal, &(inin.m),
                                        cv, face_num, P_ref, inin.Normal,
                                        inin.derivatives, false);
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
      GA_DEBUG_INFO("Instruction: end of call interpolate transformation");
      return 0;
    }
    ga_instruction_neighbour_transformation_call
    (const ga_workspace &w, ga_instruction_set::interpolate_info &i,
     pinterpolate_transformation t, fem_interpolation_context &ctxx,
     base_small_vector &No, const mesh &mm, size_type &ipt_,
     papprox_integration &pai_, bgeot::geotrans_precomp_pool &gp_pool_,
     std::map<gauss_pt_corresp, bgeot::pstored_point_tab> &neighbour_corresp_)
      : workspace(w), inin(i), trans(t), ctx(ctxx), Normal(No), m(mm),
        ipt(ipt_), pai(pai_), gp_pool(gp_pool_),
        neighbour_corresp(neighbour_corresp_) {}
  };


  struct ga_instruction_scalar_assembly : public ga_instruction {
    base_tensor &t;
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

  struct ga_instruction_fem_vector_assembly : public ga_instruction {
    base_tensor &t;
    base_vector &Vr, &Vn;
    const fem_interpolation_context &ctx;
    const gmm::sub_interval &Ir, &In;
    const mesh_fem *mfn, **mfg;
    scalar_type &coeff;
    const size_type &nbpt, &ipt;
    base_vector elem;
    bool interpolate;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: vector term assembly for fem variable");
      if (ipt == 0 || interpolate) {
        elem.resize(t.size());
        auto itt = t.begin(); auto it = elem.begin(), ite = elem.end();
        size_type nd = ((t.size()) >> 2);
        for (size_type i = 0; i < nd; ++i) {
          *it++ = (*itt++) * coeff; *it++ = (*itt++) * coeff;
          *it++ = (*itt++) * coeff; *it++ = (*itt++) * coeff;
        }
        for (; it != ite;) *it++ = (*itt++) * coeff;
        // gmm::copy(gmm::scaled(t.as_vector(), coeff), elem);
      } else {
        auto itt = t.begin(); auto it = elem.begin(), ite = elem.end();
        size_type nd = ((t.size()) >> 2);
        for (size_type i = 0; i < nd; ++i) {
          *it++ += (*itt++) * coeff; *it++ += (*itt++) * coeff;
          *it++ += (*itt++) * coeff; *it++ += (*itt++) * coeff;
        }
        for (; it != ite;) *it++ += (*itt++) * coeff;
        // gmm::add(gmm::scaled(t.as_vector(), coeff), elem);
      }
      if (ipt == nbpt-1 || interpolate) {
        const mesh_fem &mf = *(mfg ? *mfg : mfn);
        GMM_ASSERT1(mfg ? *mfg : mfn, "Internal error");
        const gmm::sub_interval &I = mf.is_reduced() ? Ir : In;
        base_vector &V = mf.is_reduced() ? Vr : Vn;
        if (!(ctx.is_convex_num_valid())) return 0;
        size_type cv_1 = ctx.convex_num();
        // size_type cv_1 = ctx.is_convex_num_valid()
        //   ? ctx.convex_num() : mf.convex_index().first_true();
        GA_DEBUG_ASSERT(V.size() >= I.first() + mf.nb_basic_dof(),
                        "Bad assembly vector size");
        auto &ct = mf.ind_scalar_basic_dof_of_element(cv_1);
        size_type qmult = mf.get_qdim();
        if (qmult > 1) qmult /= mf.fem_of_element(cv_1)->target_dim();
        size_type ifirst = I.first();
        auto ite = elem.begin();
        for (auto itc = ct.begin(); itc != ct.end(); ++itc)
          for (size_type q = 0; q < qmult; ++q)
            V[ifirst+(*itc)+q] += *ite++;
        GMM_ASSERT1(ite == elem.end(), "Internal error");
      }
      return 0;
    }
    ga_instruction_fem_vector_assembly
    (base_tensor &t_, base_vector &Vr_, base_vector &Vn_,
     const fem_interpolation_context &ctx_,
     const gmm::sub_interval &Ir_, const gmm::sub_interval &In_,
     const mesh_fem *mfn_, const mesh_fem **mfg_,
     scalar_type &coeff_,
     const size_type &nbpt_, const size_type &ipt_, bool interpolate_)
    : t(t_), Vr(Vr_), Vn(Vn_), ctx(ctx_), Ir(Ir_), In(In_), mfn(mfn_),
      mfg(mfg_), coeff(coeff_), nbpt(nbpt_), ipt(ipt_),
      interpolate(interpolate_) {}
  };

  struct ga_instruction_vector_assembly : public ga_instruction {
    base_tensor &t;
    base_vector &V;
    const gmm::sub_interval &I;
    scalar_type &coeff;
     virtual int exec() {
      GA_DEBUG_INFO("Instruction: vector term assembly for "
                    "fixed size variable");
      gmm::add(gmm::scaled(t.as_vector(), coeff), gmm::sub_vector(V, I));
      return 0;
     }
    ga_instruction_vector_assembly(base_tensor &t_, base_vector &V_,
                                   const gmm::sub_interval &I_,
                                   scalar_type &coeff_)
      : t(t_), V(V_), I(I_), coeff(coeff_) {}
  };

  struct ga_instruction_assignment : public ga_instruction {
    base_tensor &t;
    base_vector &V;
    const fem_interpolation_context &ctx;
    const im_data *imd;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: Assignement to im_data");
      imd->set_tensor(V, ctx.convex_num(), ctx.ii(), t);
      return 0;
     }
    ga_instruction_assignment(base_tensor &t_, base_vector &V_,
                              const fem_interpolation_context &ctx_,
                              const im_data *imd_)
      : t(t_), V(V_), ctx(ctx_), imd(imd_) {}
  };

  template <class MAT>
  inline void add_elem_matrix_
  (MAT &K, const std::vector<size_type> &dofs1,
   const std::vector<size_type> &dofs2, std::vector<size_type> &/*dofs1_sort*/,
   base_vector &elem, scalar_type threshold, size_type /* N */) {
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

  inline void add_elem_matrix_
  (gmm::col_matrix<gmm::rsvector<scalar_type>> &K,
   const std::vector<size_type> &dofs1, const std::vector<size_type> &dofs2,
   std::vector<size_type> &dofs1_sort,
   base_vector &elem, scalar_type threshold, size_type N) {
    size_type maxest = (N+1) * std::max(dofs1.size(), dofs2.size());
    size_type s1 = dofs1.size(), s2 = dofs2.size();
    gmm::elt_rsvector_<scalar_type> ev;

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

    base_vector::const_iterator it = elem.cbegin();
    for (size_type j = 0; j < s2; ++j) { // Iteration on columns
      if (j) it += s1;
      std::vector<gmm::elt_rsvector_<scalar_type>> &col = K[dofs2[j]];
      size_type nb = col.size();

      if (nb == 0) {
        col.reserve(maxest);
        for (size_type i = 0; i < s1; ++i) {
          size_type k = dofs1_sort[i]; ev.e = *(it+k);
          if (gmm::abs(ev.e) > threshold) { ev.c=dofs1[k]; col.push_back(ev); }
        }
      } else { // column merge
        size_type ind = 0;
        for (size_type i = 0; i < s1; ++i) {
          size_type k = dofs1_sort[i]; ev.e = *(it+k);
          if (gmm::abs(ev.e) > threshold) {
            ev.c = dofs1[k];

            size_type count = nb - ind, step, l;
            while (count > 0) {
              step = count / 2; l = ind + step;
              if (col[l].c < ev.c) { ind = ++l; count -= step + 1; }
              else count = step;
            }

            auto itc = col.begin() + ind;
            if (ind != nb && itc->c == ev.c) itc->e += ev.e;
            else {
              if (nb - ind > 1100)
                GMM_WARNING2("Inefficient addition of element in rsvector with "
                             << col.size() - ind << " non-zero entries");
              col.push_back(ev);
              if (ind != nb) {
                itc = col.begin() + ind;
                auto ite = col.end(); --ite; auto itee = ite;
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


  template <class MAT = model_real_sparse_matrix>
  struct ga_instruction_matrix_assembly : public ga_instruction {
    const base_tensor &t;
    MAT &Kr, &Kn;
    const fem_interpolation_context &ctx1, &ctx2;
    const gmm::sub_interval &Ir1, &Ir2;
    const gmm::sub_interval &In1, &In2;
    const mesh_fem *mfn1, *mfn2;
    const mesh_fem **mfg1, **mfg2;
    const scalar_type &coeff, &alpha1, &alpha2;
    const size_type &nbpt, &ipt;
    base_vector elem;
    bool interpolate;
    std::vector<size_type> dofs1, dofs2, dofs1_sort;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: matrix term assembly");
      if (ipt == 0 || interpolate) {
        elem.resize(t.size());
        auto itt = t.begin(); auto it = elem.begin(), ite = elem.end();
        scalar_type e = coeff*alpha1*alpha2;
        size_type nd = ((t.size()) >> 2);
        for (size_type i = 0; i < nd; ++i) {
          *it++ = (*itt++) * e; *it++ = (*itt++) * e;
          *it++ = (*itt++) * e; *it++ = (*itt++) * e;
        }
        for (; it != ite;) *it++ = (*itt++) * e;
        // gmm::copy(gmm::scaled(t.as_vector(), coeff*alpha1*alpha2), elem);
      } else {
        // Faster than a daxpy blas call on my config
        auto itt = t.begin(); auto it = elem.begin(), ite = elem.end();
        scalar_type e = coeff*alpha1*alpha2;
        size_type nd = ((t.size()) >> 2);
        for (size_type i = 0; i < nd; ++i) {
          *it++ += (*itt++) * e; *it++ += (*itt++) * e;
          *it++ += (*itt++) * e; *it++ += (*itt++) * e;
        }
        for (; it != ite;) *it++ += (*itt++) * e;
        // gmm::add(gmm::scaled(t.as_vector(), coeff*alpha1*alpha2), elem);
      }
      if (ipt == nbpt-1 || interpolate) {
        const mesh_fem *pmf1 = mfg1 ? *mfg1 : mfn1;
        const mesh_fem *pmf2 = mfg2 ? *mfg2 : mfn2;
        bool reduced = (pmf1 && pmf1->is_reduced())
          || (pmf2 && pmf2->is_reduced());
        MAT &K = reduced ? Kr : Kn;
        const gmm::sub_interval &I1 = reduced ? Ir1 : In1;
        const gmm::sub_interval &I2 = reduced ? Ir2 : In2;
        GA_DEBUG_ASSERT(I1.size() && I2.size(), "Internal error");

        scalar_type ninf = gmm::vect_norminf(elem);
        if (ninf == scalar_type(0)) return 0;

        size_type s1 = t.sizes()[0], s2 = t.sizes()[1];
        size_type cv1 = pmf1 ? ctx1.convex_num() : s1;
        size_type cv2 = pmf2 ? ctx2.convex_num() : s2;
        size_type N = 1;

        dofs1.assign(s1, I1.first());
        if (pmf1) {
          if (!(ctx1.is_convex_num_valid())) return 0;
          N = ctx1.N();
          auto &ct1 = pmf1->ind_scalar_basic_dof_of_element(cv1);
          size_type qmult1 = pmf1->get_qdim();
          if (qmult1 > 1) qmult1 /= pmf1->fem_of_element(cv1)->target_dim();
          auto itd = dofs1.begin();
          if (qmult1 == 1) {
            for (auto itt = ct1.begin(); itt != ct1.end(); ++itt)
              *itd++ += *itt;
          } else {
            for (auto itt = ct1.begin(); itt != ct1.end(); ++itt)
              for (size_type q = 0; q < qmult1; ++q)
                  *itd++ += *itt + q;
          }
        } else
          for (size_type i=0; i < s1; ++i) dofs1[i] += i;

        if (pmf1 == pmf2 && cv1 == cv2) {
          if (I1.first() == I2.first()) {
            add_elem_matrix_(K, dofs1, dofs1, dofs1_sort, elem, ninf*1E-14, N);
          } else {
            dofs2.resize(dofs1.size());
            for (size_type i = 0; i < dofs1.size(); ++i)
              dofs2[i] =  dofs1[i] + I2.first() - I1.first();
            add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf*1E-14, N);
          }
        } else {
          dofs2.assign(s2, I2.first());
          if (pmf2) {
            if (!(ctx2.is_convex_num_valid())) return 0;
            N = std::max(N, ctx2.N());
            auto &ct2 = pmf2->ind_scalar_basic_dof_of_element(cv2);
            size_type qmult2 = pmf2->get_qdim();
            if (qmult2 > 1) qmult2 /= pmf2->fem_of_element(cv2)->target_dim();
            auto itd = dofs2.begin();
            if (qmult2 == 1) {
              for (auto itt = ct2.begin(); itt != ct2.end(); ++itt)
                *itd++ += *itt;
            } else {
              for (auto itt = ct2.begin(); itt != ct2.end(); ++itt)
                for (size_type q = 0; q < qmult2; ++q)
                  *itd++ += *itt + q;
            }
          } else
            for (size_type i=0; i < s2; ++i) dofs2[i] += i;

          add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf*1E-14, N);
        }
      }
      return 0;
    }
    ga_instruction_matrix_assembly
    (const base_tensor &t_, MAT &Kr_, MAT &Kn_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const gmm::sub_interval &Ir1_, const gmm::sub_interval &In1_,
     const gmm::sub_interval &Ir2_, const gmm::sub_interval &In2_,
     const mesh_fem *mfn1_, const mesh_fem **mfg1_,
     const mesh_fem *mfn2_, const mesh_fem **mfg2_,
     const scalar_type &coeff_,
     const scalar_type &alpha2_, const scalar_type &alpha1_,
     const size_type &nbpt_, const size_type &ipt_, bool interpolate_)
      : t(t_), Kr(Kr_), Kn(Kn_), ctx1(ctx1_), ctx2(ctx2_),
        Ir1(Ir1_), Ir2(Ir2_), In1(In1_), In2(In2_),
        mfn1(mfn1_), mfn2(mfn2_), mfg1(mfg1_), mfg2(mfg2_),
        coeff(coeff_), alpha1(alpha1_), alpha2(alpha2_),
        nbpt(nbpt_), ipt(ipt_), interpolate(interpolate_),
        dofs1(0), dofs2(0) {}
  };

  template <class MAT = model_real_sparse_matrix>
  struct ga_instruction_matrix_assembly_standard_scalar: public ga_instruction {
    const base_tensor &t;
    MAT &K;
    const fem_interpolation_context &ctx1, &ctx2;
    const gmm::sub_interval &I1, &I2;
    const mesh_fem *pmf1, *pmf2;
    const scalar_type &coeff, &alpha1, &alpha2;
    const size_type &nbpt, &ipt;
    base_vector elem;
    std::vector<size_type> dofs1, dofs2, dofs1_sort;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: matrix term assembly for standard "
                    "scalar fems");
      if (ipt == 0) {
        elem.resize(t.size());
        auto itt = t.begin(); auto it = elem.begin(), ite = elem.end();
        scalar_type e = coeff*alpha1*alpha2;
        size_type nd = ((t.size()) >> 2);
        for (size_type i = 0; i < nd; ++i) {
          *it++ = (*itt++) * e; *it++ = (*itt++) * e;
          *it++ = (*itt++) * e; *it++ = (*itt++) * e;
        }
        for (; it != ite;) *it++ = (*itt++) * e;
        // gmm::copy(gmm::scaled(t.as_vector(), coeff*alpha1*alpha2), elem);
      } else {
        // Faster than a daxpy blas call on my config
        auto itt = t.begin(); auto it = elem.begin(), ite = elem.end();
        scalar_type e = coeff*alpha1*alpha2;
        size_type nd = ((t.size()) >> 2);
        for (size_type i = 0; i < nd; ++i) {
          *it++ += (*itt++) * e; *it++ += (*itt++) * e;
          *it++ += (*itt++) * e; *it++ += (*itt++) * e;
        }
        for (; it != ite;) *it++ += (*itt++) * e;
        // gmm::add(gmm::scaled(t.as_vector(), coeff*alpha1*alpha2), elem);
      }
      if (ipt == nbpt-1) {
        GA_DEBUG_ASSERT(I1.size() && I2.size(), "Internal error");

        scalar_type ninf = gmm::vect_norminf(elem);
        if (ninf == scalar_type(0)) return 0;

        size_type cv1 = ctx1.convex_num(), cv2 = ctx2.convex_num(), N=ctx1.N();
        if (cv1 == size_type(-1)) return 0;
        auto &ct1 = pmf1->ind_scalar_basic_dof_of_element(cv1);
        GA_DEBUG_ASSERT(ct1.size() == t.sizes()[0], "Internal error");
        dofs1.resize(ct1.size());
        for (size_type i = 0; i < ct1.size(); ++i)
          dofs1[i] = ct1[i] + I1.first();

        if (pmf2 == pmf1 && cv1 == cv2) {
          if (I1.first() == I2.first()) {
            add_elem_matrix_(K, dofs1, dofs1, dofs1_sort, elem, ninf*1E-14, N);
          } else {
            dofs2.resize(dofs1.size());
            for (size_type i = 0; i < dofs1.size(); ++i)
              dofs2[i] =  dofs1[i] + I2.first() - I1.first();
            add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf*1E-14, N);
          }
        } else {
          if (cv2 == size_type(-1)) return 0;
          auto &ct2 = pmf2->ind_scalar_basic_dof_of_element(cv2);
          GA_DEBUG_ASSERT(ct2.size() == t.sizes()[1], "Internal error");
          dofs2.resize(ct2.size());
          for (size_type i = 0; i < ct2.size(); ++i)
            dofs2[i] = ct2[i] + I2.first();
          add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf*1E-14, N);
        }
      }
      return 0;
    }
    ga_instruction_matrix_assembly_standard_scalar
    (const base_tensor &t_, MAT &Kn_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const gmm::sub_interval &In1_, const gmm::sub_interval &In2_,
     const mesh_fem *mfn1_, const mesh_fem *mfn2_,
     const scalar_type &coeff_, const scalar_type &alpha2_,
     const scalar_type &alpha1_,
     const size_type &nbpt_, const size_type &ipt_)
      : t(t_), K(Kn_), ctx1(ctx1_), ctx2(ctx2_),
        I1(In1_), I2(In2_),  pmf1(mfn1_), pmf2(mfn2_),
        coeff(coeff_), alpha1(alpha1_), alpha2(alpha2_),
        nbpt(nbpt_), ipt(ipt_) {}
  };

  template <class MAT = model_real_sparse_matrix>
  struct ga_instruction_matrix_assembly_standard_vector: public ga_instruction {
    const base_tensor &t;
    MAT &K;
    const fem_interpolation_context &ctx1, &ctx2;
    const gmm::sub_interval &I1, &I2;
    const mesh_fem *pmf1, *pmf2;
    const scalar_type &coeff, &alpha1, &alpha2;
    const size_type &nbpt, &ipt;
    mutable base_vector elem;
    mutable std::vector<size_type> dofs1, dofs2, dofs1_sort;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: matrix term assembly for standard "
                        "vector fems");
      if (ipt == 0) {
        elem.resize(t.size());
        auto itt = t.begin(); auto it = elem.begin(), ite = elem.end();
        scalar_type e = coeff*alpha1*alpha2;
        size_type nd = ((t.size()) >> 3);
        for (size_type i = 0; i < nd; ++i) {
          *it++ = (*itt++) * e; *it++ = (*itt++) * e;
          *it++ = (*itt++) * e; *it++ = (*itt++) * e;
          *it++ = (*itt++) * e; *it++ = (*itt++) * e;
          *it++ = (*itt++) * e; *it++ = (*itt++) * e;
        }
        for (; it != ite;) *it++ = (*itt++) * e;
        // gmm::copy(gmm::scaled(t.as_vector(), coeff*alpha1*alpha2), elem);
      } else {
        // (Far) faster than a daxpy blas call on my config.
        auto itt = t.begin(); auto it = elem.begin(), ite = elem.end();
        scalar_type e = coeff*alpha1*alpha2;
        size_type nd = ((t.size()) >> 3);
        for (size_type i = 0; i < nd; ++i) {
          *it++ += (*itt++) * e; *it++ += (*itt++) * e;
          *it++ += (*itt++) * e; *it++ += (*itt++) * e;
          *it++ += (*itt++) * e; *it++ += (*itt++) * e;
          *it++ += (*itt++) * e; *it++ += (*itt++) * e;
        }
        for (; it != ite;) *it++ += (*itt++) * e;
        // gmm::add(gmm::scaled(t.as_vector(), coeff*alpha1*alpha2), elem);
      }
      if (ipt == nbpt-1) {
        GA_DEBUG_ASSERT(I1.size() && I2.size(), "Internal error");

        scalar_type ninf = gmm::vect_norminf(elem);
        if (ninf == scalar_type(0)) return 0;
        size_type s1 = t.sizes()[0], s2 = t.sizes()[1], N = ctx1.N();

        size_type cv1 = ctx1.convex_num(), cv2 = ctx2.convex_num();
        if (cv1 == size_type(-1)) return 0;
        auto &ct1 = pmf1->ind_scalar_basic_dof_of_element(cv1);
        size_type qmult1 = pmf1->get_qdim();
        if (qmult1 > 1) qmult1 /= pmf1->fem_of_element(cv1)->target_dim();
        dofs1.assign(s1, I1.first());
        auto itd = dofs1.begin();
        for (auto itt = ct1.begin(); itt != ct1.end(); ++itt)
          for (size_type q = 0; q < qmult1; ++q)
            *itd++ += *itt + q;

        if (pmf2 == pmf1 && cv1 == cv2) {
          if (I1.first() == I2.first()) {
            add_elem_matrix_(K, dofs1, dofs1, dofs1_sort, elem, ninf*1E-14, N);
          } else {
            dofs2.resize(dofs1.size());
            for (size_type i = 0; i < dofs1.size(); ++i)
              dofs2[i] =  dofs1[i] + I2.first() - I1.first();
            add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf*1E-14, N);
          }
        } else {
          if (cv2 == size_type(-1)) return 0;
          auto &ct2 = pmf2->ind_scalar_basic_dof_of_element(cv2);
          size_type qmult2 = pmf2->get_qdim();
          if (qmult2 > 1) qmult2 /= pmf2->fem_of_element(cv2)->target_dim();
          dofs2.assign(s2, I2.first());
          itd = dofs2.begin();
          for (auto itt = ct2.begin(); itt != ct2.end(); ++itt)
            for (size_type q = 0; q < qmult2; ++q)
              *itd++ += *itt + q;

          add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf*1E-14, N);
        }
      }
      return 0;
    }
    ga_instruction_matrix_assembly_standard_vector
    (const base_tensor &t_, MAT &Kn_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const gmm::sub_interval &In1_, const gmm::sub_interval &In2_,
     const mesh_fem *mfn1_, const mesh_fem *mfn2_,
     const scalar_type &coeff_, const scalar_type &alpha2_,
     const scalar_type &alpha1_, const size_type &nbpt_,
     const size_type &ipt_)
      : t(t_), K(Kn_), ctx1(ctx1_), ctx2(ctx2_),
        I1(In1_), I2(In2_),  pmf1(mfn1_), pmf2(mfn2_),
        coeff(coeff_), alpha1(alpha1_), alpha2(alpha2_),
        nbpt(nbpt_), ipt(ipt_), dofs1(0), dofs2(0) {}
  };

  struct ga_instruction_matrix_assembly_standard_vector_opt10_2
    : public ga_instruction {
    const base_tensor &t;
    model_real_sparse_matrix &K;
    const fem_interpolation_context &ctx1, &ctx2;
    const gmm::sub_interval &I1, &I2;
    const mesh_fem *pmf1, *pmf2;
    const scalar_type &coeff, &alpha1, &alpha2;
    const size_type &nbpt, &ipt;
    mutable base_vector elem;
    mutable std::vector<size_type> dofs1, dofs2, dofs1_sort;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: matrix term assembly for standard "
                    "vector fems optimized for format 10 qdim 2");
      size_type s1 = t.sizes()[0], s2 = t.sizes()[1], s1_q = 2*s1;
      size_type ss1 = s1/2, ss2 = s2/2;
      scalar_type e = coeff*alpha1*alpha2;
      if (ipt == 0) {
        elem.resize(ss1*ss2);
        auto itel = elem.begin();
        for (size_type j = 0; j < ss2; ++j) {
          auto it = t.begin() + j*s1_q;
          for (size_type i = 0; i < ss1; ++i, it += 2)
            *itel++ = (*it) * e;
        }
      } else {
        auto itel = elem.begin();
        for (size_type j = 0; j < ss2; ++j) {
          auto it = t.begin() + j*s1_q;
          for (size_type i = 0; i < ss1; ++i, it += 2)
            *itel++ += (*it) * e;
        }
      }
      if (ipt == nbpt-1) {
        GA_DEBUG_ASSERT(I1.size() && I2.size(), "Internal error");

        scalar_type ninf = gmm::vect_norminf(elem) * 1E-14;
        if (ninf == scalar_type(0)) return 0;
        size_type N = ctx1.N();
        size_type cv1 = ctx1.convex_num(), cv2 = ctx2.convex_num();
        size_type i1 = I1.first(), i2 = I2.first();
        if (cv1 == size_type(-1)) return 0;
        auto &ct1 = pmf1->ind_scalar_basic_dof_of_element(cv1);
        dofs1.resize(ss1);
        for (size_type i = 0; i < ss1; ++i) dofs1[i] = i1 + ct1[i];

        if (pmf2 == pmf1 && cv1 == cv2) {
          if (i1 == i2) {
            add_elem_matrix_(K, dofs1, dofs1, dofs1_sort, elem, ninf, N);
            for (size_type i = 0; i < ss1; ++i) (dofs1[i])++;
            add_elem_matrix_(K, dofs1, dofs1, dofs1_sort, elem, ninf, N);
          } else {
            dofs2.resize(ss2);
            for (size_type i = 0; i < ss2; ++i) dofs2[i] = i2 + ct1[i];
            add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf, N);
            for (size_type i = 0; i < ss1; ++i) (dofs1[i])++;
            for (size_type i = 0; i < ss2; ++i) (dofs2[i])++;
            add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf, N);
          }
        } else {
          if (cv2 == size_type(-1)) return 0;
          auto &ct2 = pmf2->ind_scalar_basic_dof_of_element(cv2);
          dofs2.resize(ss2);
          for (size_type i = 0; i < ss2; ++i) dofs2[i] = i2 + ct2[i];
          add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf, N);
          for (size_type i = 0; i < ss1; ++i) (dofs1[i])++;
          for (size_type i = 0; i < ss2; ++i) (dofs2[i])++;
          add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf, N);
        }
      }
      return 0;
    }
    ga_instruction_matrix_assembly_standard_vector_opt10_2
    (const base_tensor &t_, model_real_sparse_matrix &Kn_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const gmm::sub_interval &In1_, const gmm::sub_interval &In2_,
     const mesh_fem *mfn1_, const mesh_fem *mfn2_,
     const scalar_type &coeff_, const scalar_type &alpha2_,
     const scalar_type &alpha1_, const size_type &nbpt_,
     const size_type &ipt_)
      : t(t_), K(Kn_), ctx1(ctx1_), ctx2(ctx2_),
        I1(In1_), I2(In2_),  pmf1(mfn1_), pmf2(mfn2_),
        coeff(coeff_), alpha1(alpha1_), alpha2(alpha2_),
        nbpt(nbpt_), ipt(ipt_), dofs1(0), dofs2(0) {}
  };

  struct ga_instruction_matrix_assembly_standard_vector_opt10_3
    : public ga_instruction {
    const base_tensor &t;
    model_real_sparse_matrix &K;
    const fem_interpolation_context &ctx1, &ctx2;
    const gmm::sub_interval &I1, &I2;
    const mesh_fem *pmf1, *pmf2;
    const scalar_type &coeff, &alpha1, &alpha2;
    const size_type &nbpt, &ipt;
    mutable base_vector elem;
    mutable std::vector<size_type> dofs1, dofs2, dofs1_sort;
    virtual int exec() {
      GA_DEBUG_INFO("Instruction: matrix term assembly for standard "
                    "vector fems optimized for format 10 qdim 3");
      size_type s1 = t.sizes()[0], s2 = t.sizes()[1], s1_q = 3*s1;
      size_type ss1 = s1/3, ss2 = s2/3;
      scalar_type e = coeff*alpha1*alpha2;
      if (ipt == 0) {
        elem.resize(ss1*ss2);
        auto itel = elem.begin();
        for (size_type j = 0; j < ss2; ++j) {
          auto it = t.begin() + j*s1_q;
          for (size_type i = 0; i < ss1; ++i, it += 3)
            *itel++ = (*it) * e;
        }
      } else {
        auto itel = elem.begin();
        for (size_type j = 0; j < ss2; ++j) {
          auto it = t.begin() + j*s1_q;
          for (size_type i = 0; i < ss1; ++i, it += 3)
            *itel++ += (*it) * e;
        }
      }
      if (ipt == nbpt-1) {
        GA_DEBUG_ASSERT(I1.size() && I2.size(), "Internal error");

        scalar_type ninf = gmm::vect_norminf(elem)*1E-14;
        if (ninf == scalar_type(0)) return 0;
        size_type N = ctx1.N();
        size_type cv1 = ctx1.convex_num(), cv2 = ctx2.convex_num();
        size_type i1 = I1.first(), i2 = I2.first();
        if (cv1 == size_type(-1)) return 0;
        auto &ct1 = pmf1->ind_scalar_basic_dof_of_element(cv1);
        dofs1.resize(ss1);
        for (size_type i = 0; i < ss1; ++i) dofs1[i] = i1 + ct1[i];

        if (pmf2 == pmf1 && cv1 == cv2) {
          if (i1 == i2) {
            add_elem_matrix_(K, dofs1, dofs1, dofs1_sort, elem, ninf, N);
            for (size_type i = 0; i < ss1; ++i) (dofs1[i])++;
            add_elem_matrix_(K, dofs1, dofs1, dofs1_sort, elem, ninf, N);
            for (size_type i = 0; i < ss1; ++i) (dofs1[i])++;
            add_elem_matrix_(K, dofs1, dofs1, dofs1_sort, elem, ninf, N);
          } else {
            dofs2.resize(ss2);
            for (size_type i = 0; i < ss2; ++i) dofs2[i] = i2 + ct1[i];
            add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf, N);
            for (size_type i = 0; i < ss1; ++i) (dofs1[i])++;
            for (size_type i = 0; i < ss2; ++i) (dofs2[i])++;
            add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf, N);
            for (size_type i = 0; i < ss1; ++i) (dofs1[i])++;
            for (size_type i = 0; i < ss2; ++i) (dofs2[i])++;
            add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf, N);
          }
        } else {
          if (cv2 == size_type(-1)) return 0;
          auto &ct2 = pmf2->ind_scalar_basic_dof_of_element(cv2);
          dofs2.resize(ss2);
          for (size_type i = 0; i < ss2; ++i) dofs2[i] = i2 + ct2[i];
          add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf, N);
          for (size_type i = 0; i < ss1; ++i) (dofs1[i])++;
          for (size_type i = 0; i < ss2; ++i) (dofs2[i])++;
          add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf, N);
          for (size_type i = 0; i < ss1; ++i) (dofs1[i])++;
          for (size_type i = 0; i < ss2; ++i) (dofs2[i])++;
          add_elem_matrix_(K, dofs1, dofs2, dofs1_sort, elem, ninf, N);
        }
      }
      return 0;
    }
    ga_instruction_matrix_assembly_standard_vector_opt10_3
    (const base_tensor &t_, model_real_sparse_matrix &Kn_,
     const fem_interpolation_context &ctx1_,
     const fem_interpolation_context &ctx2_,
     const gmm::sub_interval &In1_, const gmm::sub_interval &In2_,
     const mesh_fem *mfn1_, const mesh_fem *mfn2_,
     const scalar_type &coeff_, const scalar_type &alpha2_,
     const scalar_type &alpha1_, const size_type &nbpt_,
     const size_type &ipt_)
      : t(t_), K(Kn_), ctx1(ctx1_), ctx2(ctx2_),
        I1(In1_), I2(In2_),  pmf1(mfn1_), pmf2(mfn2_),
        coeff(coeff_), alpha1(alpha1_), alpha2(alpha2_),
        nbpt(nbpt_), ipt(ipt_), dofs1(0), dofs2(0) {}
  };




  //=========================================================================
  // Structure dealing with user defined environment : constant, variables,
  // functions, operators.
  //=========================================================================

  void ga_workspace::init() {
    // allocate own storage for K an V to be used unless/until external
    // storage is provided with set_assembled_matrix/vector
    K = std::make_shared<model_real_sparse_matrix>(2,2);
    V = std::make_shared<base_vector>(2);
    // add default transformations
    add_interpolate_transformation
      ("neighbour_elt", interpolate_transformation_neighbour_instance());
  }

  // variables and variable groups
  void ga_workspace::add_fem_variable
  (const std::string &name, const mesh_fem &mf,
   const gmm::sub_interval &I, const model_real_plain_vector &VV) {
    variables[name] = var_description(true, true, &mf, I, &VV, 0, 1);
  }

  void ga_workspace::add_fixed_size_variable
  (const std::string &name,
   const gmm::sub_interval &I, const model_real_plain_vector &VV) {
    variables[name] = var_description(true, false, 0, I, &VV, 0,
                                      dim_type(gmm::vect_size(VV)));
  }

  void ga_workspace::add_fem_constant
  (const std::string &name, const mesh_fem &mf,
   const model_real_plain_vector &VV) {
    GMM_ASSERT1(mf.nb_dof(), "The provided mesh_fem of variable" << name
                             << "has zero degrees of freedom.");
    size_type Q = gmm::vect_size(VV)/mf.nb_dof();
    if (Q == 0) Q = size_type(1);
    variables[name] = var_description(false, true, &mf,
                                      gmm::sub_interval(), &VV, 0, Q);
  }

  void ga_workspace::add_fixed_size_constant
  (const std::string &name, const model_real_plain_vector &VV) {
    variables[name] = var_description(false, false, 0,
                                      gmm::sub_interval(), &VV, 0,
                                      gmm::vect_size(VV));
  }

  void ga_workspace::add_im_data(const std::string &name, const im_data &imd,
                                 const model_real_plain_vector &VV) {
    variables[name] = var_description
      (false, false, 0, gmm::sub_interval(), &VV, &imd,
       gmm::vect_size(VV)/(imd.nb_filtered_index() * imd.nb_tensor_elem()));
  }


  bool ga_workspace::variable_exists(const std::string &name) const {
    return (md && md->variable_exists(name)) ||
      (parent_workspace && parent_workspace->variable_exists(name)) ||
      (variables.find(name) != variables.end());
  }

  bool ga_workspace::variable_group_exists(const std::string &name) const {
    return (variable_groups.find(name) != variable_groups.end()) ||
      (md && md->variable_group_exists(name)) ||
      (parent_workspace && parent_workspace->variable_group_exists(name));
  }

  const std::vector<std::string>&
  ga_workspace::variable_group(const std::string &group_name) const {
    std::map<std::string, std::vector<std::string> >::const_iterator
      it = variable_groups.find(group_name);
    if (it != variable_groups.end())
      return (variable_groups.find(group_name))->second;
    if (md && md->variable_group_exists(group_name))
      return md->variable_group(group_name);
    if (parent_workspace &&
        parent_workspace->variable_group_exists(group_name))
      return parent_workspace->variable_group(group_name);
    GMM_ASSERT1(false, "Undefined variable group " << group_name);
  }

  const std::string&
  ga_workspace::first_variable_of_group(const std::string &name) const {
    const std::vector<std::string> &t = variable_group(name);
    GMM_ASSERT1(t.size(), "Variable group " << name << " is empty");
    return t[0];
  }

  bool ga_workspace::is_constant(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) return !(it->second.is_variable);
    if (variable_group_exists(name))
      return is_constant(first_variable_of_group(name));
    if (md && md->variable_exists(name)) {
      if (enable_all_md_variables) return md->is_true_data(name);
      return md->is_data(name);
    }
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->is_constant(name);
    GMM_ASSERT1(false, "Undefined variable " << name);
  }

  bool ga_workspace::is_disabled_variable(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) return false;
    if (variable_group_exists(name))
      return is_disabled_variable(first_variable_of_group(name));
    if (md && md->variable_exists(name)) {
      if (enable_all_md_variables) return false;
      return md->is_disabled_variable(name);
    }
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->is_disabled_variable(name);
    GMM_ASSERT1(false, "Undefined variable " << name);
  }

  const scalar_type &
  ga_workspace::factor_of_variable(const std::string &name) const {
    static const scalar_type one(1);
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) return one;
    if (variable_group_exists(name))
      return one;
    if (md && md->variable_exists(name)) return md->factor_of_variable(name);
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->factor_of_variable(name);
    GMM_ASSERT1(false, "Undefined variable " << name);
  }

  const gmm::sub_interval &
  ga_workspace::interval_of_disabled_variable(const std::string &name) const {
    std::map<std::string, gmm::sub_interval>::const_iterator
      it1 = int_disabled_variables.find(name);
    if (it1 != int_disabled_variables.end()) return it1->second;
    if (md->is_affine_dependent_variable(name))
      return interval_of_disabled_variable(md->org_variable(name));

    size_type first = md->nb_dof();
    for (const std::pair<std::string, gmm::sub_interval> &p
         : int_disabled_variables)
      first = std::max(first, p.second.last());

    int_disabled_variables[name]
      = gmm::sub_interval(first, gmm::vect_size(value(name)));
    return int_disabled_variables[name];
  }

  const gmm::sub_interval &
  ga_workspace::interval_of_variable(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) return it->second.I;
    if (md && md->variable_exists(name)) {
      if (enable_all_md_variables && md->is_disabled_variable(name))
        return interval_of_disabled_variable(name);
      return md->interval_of_variable(name);
    }
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->interval_of_variable(name);
    GMM_ASSERT1(false, "Undefined variable " << name);
  }

  const mesh_fem *
  ga_workspace::associated_mf(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end())
      return it->second.is_fem_dofs ? it->second.mf : 0;
    if (md && md->variable_exists(name))
      return md->pmesh_fem_of_variable(name);
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->associated_mf(name);
    if (variable_group_exists(name))
      return associated_mf(first_variable_of_group(name));
    GMM_ASSERT1(false, "Undefined variable or group " << name);
  }

  const im_data *
  ga_workspace::associated_im_data(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) return it->second.imd;
    if (md && md->variable_exists(name))
      return md->pim_data_of_variable(name);
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->associated_im_data(name);
    if (variable_group_exists(name)) return 0;
    GMM_ASSERT1(false, "Undefined variable " << name);
  }

  size_type ga_workspace::qdim(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) {
      const mesh_fem *mf =  it->second.is_fem_dofs ? it->second.mf : 0;
      const im_data *imd = it->second.imd;
      size_type n = it->second.qdim();
      if (mf) {
        return n * mf->get_qdim();
      } else if (imd) {
        return n * imd->tensor_size().total_size();
      }
      return n;
    }
    if (md && md->variable_exists(name))
      return md->qdim_of_variable(name);
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->qdim(name);
    if (variable_group_exists(name))
      return qdim(first_variable_of_group(name));
    GMM_ASSERT1(false, "Undefined variable or group " << name);
  }

  bgeot::multi_index
  ga_workspace::qdims(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) {
      const mesh_fem *mf =  it->second.is_fem_dofs ? it->second.mf : 0;
      const im_data *imd = it->second.imd;
      size_type n = it->second.qdim();
      if (mf) {
        bgeot::multi_index mi = mf->get_qdims();
        if (n > 1 || it->second.qdims.size() > 1) {
          size_type i = 0;
          if (mi.back() == 1) { mi.back() *= it->second.qdims[0]; ++i; }
          for (; i < it->second.qdims.size(); ++i)
            mi.push_back(it->second.qdims[i]);
        }
        return mi;
      } else if (imd) {
        bgeot::multi_index mi = imd->tensor_size();
        size_type q = n / imd->nb_filtered_index();
        GMM_ASSERT1(q % imd->nb_tensor_elem() == 0,
                    "Invalid mesh im data vector");
        if (n > 1 || it->second.qdims.size() > 1) {
          size_type i = 0;
          if (mi.back() == 1) { mi.back() *= it->second.qdims[0]; ++i; }
          for (; i < it->second.qdims.size(); ++i)
            mi.push_back(it->second.qdims[i]);
        }
        return mi;
      }
      return it->second.qdims;
    }
    if (md && md->variable_exists(name))
      return md->qdims_of_variable(name);
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->qdims(name);
    if (variable_group_exists(name))
      return qdims(first_variable_of_group(name));
    GMM_ASSERT1(false, "Undefined variable or group " << name);
  }

  const model_real_plain_vector &
  ga_workspace::value(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end())
      return *(it->second.V);
    if (md && md->variable_exists(name))
      return md->real_variable(name);
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->value(name);
    if (variable_group_exists(name))
      return value(first_variable_of_group(name));
    GMM_ASSERT1(false, "Undefined variable or group " << name);
  }

  scalar_type ga_workspace::get_time_step() const {
    if (md) return md->get_time_step();
    if (parent_workspace) return parent_workspace->get_time_step();
    GMM_ASSERT1(false, "No time step defined here");
  }


  // Macros
  bool ga_workspace::macro_exists(const std::string &name) const {
    if (macros.find(name) != macros.end()) return true;
    if (md && md->macro_exists(name)) return true;
    if (parent_workspace &&
        parent_workspace->macro_exists(name)) return true;
    return false;
  }

  const std::string&
  ga_workspace::get_macro(const std::string &name) const {
    std::map<std::string, std::string>::const_iterator it=macros.find(name);
    if (it != macros.end()) return it->second;
    if (md && md->macro_exists(name)) return md->get_macro(name);
    if (parent_workspace &&
        parent_workspace->macro_exists(name))
      return parent_workspace->get_macro(name);
    GMM_ASSERT1(false, "Undefined macro");
  }

  ga_tree &
  ga_workspace::macro_tree(const std::string &name, size_type meshdim,
                           size_type ref_elt_dim, bool ignore_X) const {
    GMM_ASSERT1(macro_exists(name), "Undefined macro");
    auto it = macro_trees.find(name);
    bool to_be_analyzed = false;
    m_tree *mt = 0;

    if (it == macro_trees.end()) {
      mt = &(macro_trees[name]);
      to_be_analyzed = true;
    } else {
      mt = &(it->second);
      GMM_ASSERT1(mt->ptree, "Recursive definition of macro " << name);
      if (mt->meshdim != meshdim || mt->ignore_X != ignore_X) {
        to_be_analyzed = true;
        delete mt->ptree; mt->ptree = 0;
      }
    }
    if (to_be_analyzed) {
      ga_tree tree;
      ga_read_string(get_macro(name), tree);
      ga_semantic_analysis(get_macro(name), tree, *this, meshdim, ref_elt_dim,
                           false, ignore_X, 3);
      GMM_ASSERT1(tree.root, "Invalid macro");
      mt->ptree = new ga_tree(tree);
      mt->meshdim = meshdim;
      mt->ignore_X = ignore_X;
    }
    return *(mt->ptree);
  }

  void ga_workspace::add_interpolate_transformation
  (const std::string &name, pinterpolate_transformation ptrans) {
    if (transformations.find(name) != transformations.end())
      GMM_ASSERT1(name.compare("neighbour_elt"), "neighbour_elt is a "
                  "reserved interpolate transformation name");
    transformations[name] = ptrans;
  }

  bool ga_workspace::interpolate_transformation_exists
  (const std::string &name) const {
    return (md && md->interpolate_transformation_exists(name)) ||
      (parent_workspace &&
       parent_workspace->interpolate_transformation_exists(name)) ||
      (transformations.find(name) != transformations.end());
  }

  pinterpolate_transformation
  ga_workspace::interpolate_transformation(const std::string &name) const {
    std::map<std::string, pinterpolate_transformation>::const_iterator
      it = transformations.find(name);
    if (it != transformations.end()) return it->second;
    if (md && md->interpolate_transformation_exists(name))
      return md->interpolate_transformation(name);
    if (parent_workspace &&
       parent_workspace->interpolate_transformation_exists(name))
      return parent_workspace->interpolate_transformation(name);
    GMM_ASSERT1(false, "Inexistent transformation " << name);
  }

  bool ga_workspace::elementary_transformation_exists
  (const std::string &name) const {
    return (md && md->elementary_transformation_exists(name)) ||
      (parent_workspace &&
       parent_workspace->elementary_transformation_exists(name)) ||
      (elem_transformations.find(name) != elem_transformations.end());
  }

  pelementary_transformation
  ga_workspace::elementary_transformation(const std::string &name) const {
    std::map<std::string, pelementary_transformation>::const_iterator
      it = elem_transformations.find(name);
    if (it != elem_transformations.end()) return it->second;
    if (md && md->elementary_transformation_exists(name))
      return md->elementary_transformation(name);
    if (parent_workspace &&
       parent_workspace->elementary_transformation_exists(name))
      return parent_workspace->elementary_transformation(name);
    GMM_ASSERT1(false, "Inexistent elementary transformation " << name);
  }

  const mesh_region &
  ga_workspace::register_region(const mesh &m, const mesh_region &region) {
    if (&m == &dummy_mesh())
      return dummy_mesh_region();

    std::list<mesh_region> &lmr = registred_mesh_regions[&m];
    for (const mesh_region &rg : lmr)
      if (rg.compare(m, region, m)) return rg;
    lmr.push_back(region);
    return lmr.back();
  }


  void ga_workspace::add_tree(ga_tree &tree, const mesh &m,
                              const mesh_im &mim, const mesh_region &rg,
                              const std::string &expr,
                              size_type add_derivative_order,
                              bool function_expr, size_type for_interpolation,
                              const std::string varname_interpolation) {
    if (tree.root) {

      // Eliminate the term if it corresponds to disabled variables
      if ((tree.root->test_function_type >= 1 &&
           is_disabled_variable(tree.root->name_test1)) ||
          (tree.root->test_function_type >= 2 &&
           is_disabled_variable(tree.root->name_test2))) {
        // cout<<"disabling term ";  ga_print_node(tree.root, cout); cout<<endl;
        return;
      }
      // cout << "add tree with tests functions of " <<  tree.root->name_test1
      //      << " and " << tree.root->name_test2 << endl;
      //      ga_print_node(tree.root, cout); cout << endl;
      bool remain = true;
      size_type order = 0, ind_tree = 0;

      if (for_interpolation)
        order = size_type(-1) - add_derivative_order;
      else {
        switch(tree.root->test_function_type) {
        case 0: order = 0; break;
        case 1: order = 1; break;
        case 3: order = 2; break;
        default: GMM_ASSERT1(false, "Inconsistent term "
                             << tree.root->test_function_type);
        }
      }

      bool found = false;
      for (size_type i = 0; i < trees.size(); ++i) {
        if (trees[i].mim == &mim && trees[i].m == &m &&
            trees[i].order == order &&
            trees[i].name_test1.compare(tree.root->name_test1) == 0 &&
            trees[i].interpolate_name_test1.compare
            (tree.root->interpolate_name_test1) == 0 &&
            trees[i].name_test2.compare(tree.root->name_test2) == 0 &&
            trees[i].interpolate_name_test2.compare
            (tree.root->interpolate_name_test2) == 0 &&
            trees[i].rg == &rg && trees[i].interpolation == for_interpolation &&
            trees[i].varname_interpolation.compare(varname_interpolation)==0) {
          ga_tree &ftree = *(trees[i].ptree);

          ftree.insert_node(ftree.root, GA_NODE_OP);
          ftree.root->op_type = GA_PLUS;
          ftree.root->children.resize(2, nullptr);
          ftree.copy_node(tree.root, ftree.root, ftree.root->children[1]);
          ga_semantic_analysis("", ftree, *this, m.dim(),
                               ref_elt_dim_of_mesh(m), false, function_expr);
          found = true;
          break;
        }
      }

      if (!found) {
        ind_tree = trees.size(); remain = false;
        trees.push_back(tree_description());
        trees.back().mim = &mim; trees.back().m = &m;
        trees.back().rg = &rg;
        trees.back().ptree = new ga_tree;
        trees.back().ptree->swap(tree);
        pga_tree_node root = trees.back().ptree->root;
        trees.back().name_test1 = root->name_test1;
        trees.back().name_test2 = root->name_test2;
        trees.back().interpolate_name_test1 = root->interpolate_name_test1;
        trees.back().interpolate_name_test2 = root->interpolate_name_test2;
        trees.back().order = order;
        trees.back().interpolation = for_interpolation;
        trees.back().varname_interpolation = varname_interpolation;
       }

      if (for_interpolation == 0 && order < add_derivative_order) {
        std::set<var_trans_pair> expr_variables;
        ga_extract_variables((remain ? tree : *(trees[ind_tree].ptree)).root,
                             *this, m, expr_variables, true);
        for (const var_trans_pair &var : expr_variables) {
          if (!(is_constant(var.varname))) {
            ga_tree dtree = (remain ? tree : *(trees[ind_tree].ptree));
            // cout << "Derivation with respect to " << var.first << " : "
            //     << var.second << " of " << ga_tree_to_string(dtree) << endl;
            GA_TIC;
            ga_derivative(dtree, *this, m, var.varname, var.transname, 1+order);
            // cout << "Result : " << ga_tree_to_string(dtree) << endl;
            GA_TOCTIC("Derivative time");
            ga_semantic_analysis(expr, dtree, *this, m.dim(),
                                 ref_elt_dim_of_mesh(m), false, function_expr);
            GA_TOCTIC("Analysis after Derivative time");
            // cout << "after analysis "  << ga_tree_to_string(dtree) << endl;
            add_tree(dtree, m, mim, rg, expr, add_derivative_order,
                     function_expr, for_interpolation, varname_interpolation);
          }
        }
      }
    }
  }

  ga_workspace::m_tree::~m_tree() { if (ptree) delete ptree; }
  ga_workspace::m_tree::m_tree(const m_tree& o)
    : ptree(o.ptree), meshdim(o.meshdim), ignore_X(o.ignore_X)
  { if (o.ptree) ptree = new ga_tree(*(o.ptree)); }
  ga_workspace::m_tree &ga_workspace::m_tree::operator =(const m_tree& o) {
    if (ptree) delete ptree;
    ptree = o.ptree; meshdim = o.meshdim; ignore_X = o.ignore_X;
    if (o.ptree) ptree = new ga_tree(*(o.ptree));
    return *this;
  }

  size_type ga_workspace::add_expression(const std::string &expr,
                                         const mesh_im &mim,
                                         const mesh_region &rg_,
                                         size_type add_derivative_order) {
    const mesh_region &rg = register_region(mim.linked_mesh(), rg_);
    // cout << "adding expression " << expr << endl;
    GA_TIC;
    size_type max_order = 0;
    std::vector<ga_tree> ltrees(1);
    ga_read_string(expr, ltrees[0]);
    // cout << "read : " << ga_tree_to_string(ltrees[0])  << endl;
    ga_semantic_analysis(expr, ltrees[0], *this, mim.linked_mesh().dim(),
                         ref_elt_dim_of_mesh(mim.linked_mesh()),
                         false, false, 1);
    // cout << "analysed : " << ga_tree_to_string(ltrees[0]) << endl;
    GA_TOC("First analysis time");
    if (ltrees[0].root) {
      if (test1.size() > 1 || test2.size() > 1) {
        size_type ntest2 = std::max(size_type(1), test2.size());
        size_type nb_ltrees = test1.size()*ntest2;
        ltrees.resize(nb_ltrees);
        for (size_type i = 1; i < nb_ltrees; ++i) ltrees[i] = ltrees[0];
        std::set<var_trans_pair>::iterator it1 = test1.begin();
        for (size_type i = 0; i < test1.size(); ++i, ++it1) {
          std::set<var_trans_pair>::iterator it2 = test2.begin();
          for (size_type j = 0; j < ntest2; ++j) {
            selected_test1 = *it1;
            if (test2.size()) selected_test2 = *it2++;
            // cout << "analysis with " << selected_test1.first << endl;
            ga_semantic_analysis(expr, ltrees[i*ntest2+j], *this,
                                 mim.linked_mesh().dim(),
                                 ref_elt_dim_of_mesh(mim.linked_mesh()),
                                 false, false, 2);
            // cout <<"split: "<< ga_tree_to_string(ltrees[i*ntest2+j]) << endl;
          }
        }
      }

      for (size_type i = 0; i < ltrees.size(); ++i) {
        if (ltrees[i].root) {
          // cout << "adding tree " << ga_tree_to_string(ltrees[i]) << endl;
          max_order = std::max(ltrees[i].root->nb_test_functions(), max_order);
          add_tree(ltrees[i], mim.linked_mesh(), mim, rg, expr,
                   add_derivative_order, true, 0, "");
        }
      }
    }
    GA_TOC("Time for add expression");
    return max_order;
  }

  void ga_workspace::add_function_expression(const std::string &expr) {
    ga_tree tree;
    ga_read_string(expr, tree);
    ga_semantic_analysis(expr, tree, *this, 1, 1, false, true);
    if (tree.root) {
      // GMM_ASSERT1(tree.root->nb_test_functions() == 0,
      //            "Invalid function expression");
      add_tree(tree, dummy_mesh(), dummy_mesh_im(), dummy_mesh_region(),
               expr, 0, true, 0, "");
    }
  }

  void ga_workspace::add_interpolation_expression(const std::string &expr,
                                                  const mesh &m,
                                                  const mesh_region &rg_) {
    const mesh_region &rg = register_region(m, rg_);
    ga_tree tree;
    ga_read_string(expr, tree);
    ga_semantic_analysis(expr, tree, *this, m.dim(), ref_elt_dim_of_mesh(m),
                         false, false);
    if (tree.root) {
      // GMM_ASSERT1(tree.root->nb_test_functions() == 0,
      //            "Invalid expression containing test functions");
      add_tree(tree, m, dummy_mesh_im(), rg, expr, 0, false, 1, "");
    }
  }

  void ga_workspace::add_interpolation_expression(const std::string &expr,
                                                  const mesh_im &mim,
                                                  const mesh_region &rg_) {
    const mesh &m = mim.linked_mesh();
    const mesh_region &rg = register_region(m, rg_);
    ga_tree tree;
    ga_read_string(expr, tree);
    ga_semantic_analysis(expr, tree, *this, m.dim(), ref_elt_dim_of_mesh(m),
                         false, false);
    if (tree.root) {
      GMM_ASSERT1(tree.root->nb_test_functions() == 0,
                  "Invalid expression containing test functions");
      add_tree(tree, m, mim, rg, expr, 0, false, 1, "");
    }
  }

  void ga_workspace::add_assignment_expression
  (const std::string &varname, const std::string &expr, const mesh_region &rg_,
   size_type order, bool before) {
    const im_data *imd = associated_im_data(varname);
    GMM_ASSERT1(imd != 0, "Only applicable to im_data");
    const mesh_im &mim = imd->linked_mesh_im();
    const mesh &m = mim.linked_mesh();
    const mesh_region &rg = register_region(m, rg_);
    ga_tree tree;
    ga_read_string(expr, tree);
    ga_semantic_analysis(expr, tree, *this, m.dim(), ref_elt_dim_of_mesh(m),
                         false, false);
    if (tree.root) {
      GMM_ASSERT1(tree.root->nb_test_functions() == 0,
                  "Invalid expression containing test functions");
      add_tree(tree, m, mim, rg, expr, order+1, false, (before ? 1 : 2),
               varname);
    }
  }

  size_type ga_workspace::nb_trees() const { return trees.size(); }

  ga_workspace::tree_description &ga_workspace::tree_info(size_type i)
  { return trees[i]; }

  bool ga_workspace::used_variables(std::vector<std::string> &vl,
                                    std::vector<std::string> &vl_test1,
                                    std::vector<std::string> &vl_test2,
                                    std::vector<std::string> &dl,
                                    size_type order) {
    bool islin = true;
    std::set<var_trans_pair> vll, dll;
    for (size_type i = 0; i < vl.size(); ++i)
      vll.insert(var_trans_pair(vl[i], ""));
    for (size_type i = 0; i < dl.size(); ++i)
      dll.insert(var_trans_pair(dl[i], ""));

    for (size_type i = 0; i < trees.size(); ++i) {
      ga_workspace::tree_description &td =  trees[i];
      std::set<var_trans_pair> dllaux;
      bool fv = ga_extract_variables(td.ptree->root, *this, *(td.m),
                                     dllaux, false);

      if (td.order == order) {
        for (std::set<var_trans_pair>::iterator it = dllaux.begin();
             it!=dllaux.end(); ++it)
          dll.insert(*it);
      }
      switch (td.order) {
      case 0:  break;
      case 1:
        if (td.order == order) {
          if (variable_group_exists(td.name_test1)) {
            for (const std::string &t : variable_group(td.name_test1))
              vll.insert(var_trans_pair(t, td.interpolate_name_test1));
          } else {
            vll.insert(var_trans_pair(td.name_test1,
                                      td.interpolate_name_test1));
          }
          bool found = false;
          for (const std::string &t : vl_test1)
            if (td.name_test1.compare(t) == 0)
              found = true;
          if (!found)
            vl_test1.push_back(td.name_test1);
        }
        break;
      case 2:
        if (td.order == order) {
          if (variable_group_exists(td.name_test1)) {
            for (const std::string &t : variable_group(td.name_test1))
              vll.insert(var_trans_pair(t, td.interpolate_name_test1));
          } else {
            vll.insert(var_trans_pair(td.name_test1,
                                      td.interpolate_name_test1));
          }
          if (variable_group_exists(td.name_test2)) {
            for (const std::string &t : variable_group(td.name_test2))
              vll.insert(var_trans_pair(t, td.interpolate_name_test2));
          } else {
            vll.insert(var_trans_pair(td.name_test2,
                                      td.interpolate_name_test2));
          }
          bool found = false;
          for (size_type j = 0; j < vl_test1.size(); ++j)
            if ((td.name_test1.compare(vl_test1[j]) == 0) &&
                (td.name_test2.compare(vl_test2[j]) == 0))
              found = true;
          if (!found) {
            vl_test1.push_back(td.name_test1);
            vl_test2.push_back(td.name_test2);
          }
        }
        if (fv) islin = false;
        break;
      }
    }
    vl.clear();
    for (const auto &var : vll)
      if (vl.size() == 0 || var.varname.compare(vl.back()))
        vl.push_back(var.varname);
    dl.clear();
    for (const auto &var : dll)
      if (dl.size() == 0 || var.varname.compare(dl.back()))
        dl.push_back(var.varname);

    return islin;
  }

  void ga_workspace::define_variable_group(const std::string &group_name,
                                           const std::vector<std::string> &nl) {
    GMM_ASSERT1(!(variable_exists(group_name)), "The name of a group of "
                "variables cannot be the same as a variable name");

    std::set<const mesh *> ms;
    bool is_data_ = false;
    for (size_type i = 0; i < nl.size(); ++i) {
      if (i == 0)
        is_data_ = is_constant(nl[i]);
      else {
        GMM_ASSERT1(is_data_ == is_constant(nl[i]),
                    "It is not possible to mix variables and data in a group");
      }
      GMM_ASSERT1(variable_exists(nl[i]),
                  "All variables in a group have to exist in the model");
      const mesh_fem *mf = associated_mf(nl[i]);
      GMM_ASSERT1(mf, "Variables in a group should be fem variables");
      GMM_ASSERT1(ms.find(&(mf->linked_mesh())) == ms.end(),
                  "Two variables in a group cannot share the same mesh");
      ms.insert(&(mf->linked_mesh()));
    }
    variable_groups[group_name] = nl;
  }


  const std::string &ga_workspace::variable_in_group
  (const std::string &group_name, const mesh &m) const {
    if (variable_group_exists(group_name)) {
      for (const std::string &t : variable_group(group_name))
        if (&(associated_mf(t)->linked_mesh()) == &m)
          return t;
      GMM_ASSERT1(false, "No variable in this group for the given mesh");
    } else
      return group_name;
  }


  void ga_workspace::assembly(size_type order) {
    size_type ndof;
    const ga_workspace *w = this;
    while (w->parent_workspace) w = w->parent_workspace;
    if (w->md) ndof = w->md->nb_dof(); // To eventually call actualize_sizes()

    GA_TIC;
    ga_instruction_set gis;
    ga_compile(*this, gis, order);
    ndof = gis.nb_dof;
    size_type max_dof =  gis.max_dof;
    GA_TOCTIC("Compile time");

    if (order == 2) {
      if (K.use_count()) {
        gmm::clear(*K);
        gmm::resize(*K, max_dof, max_dof);
      }
      gmm::clear(unreduced_K);
      gmm::resize(unreduced_K, ndof, ndof);
    }
    if (order == 1) {
      if (V.use_count()) {
        gmm::clear(*V);
        gmm::resize(*V, max_dof);
      }
      gmm::clear(unreduced_V);
      gmm::resize(unreduced_V, ndof);
    }
    E = 0;
    GA_TOCTIC("Init time");

    ga_exec(gis, *this);
    GA_TOCTIC("Exec time");

    if (order == 1) {
      MPI_SUM_VECTOR(assembled_vector());
      MPI_SUM_VECTOR(unreduced_V);
    } else if (order == 0) {
      assembled_potential() = MPI_SUM_SCALAR(assembled_potential());
    }

    // Deal with reduced fems.
    if (order > 0) {
      std::set<std::string> vars_vec_done;
      std::set<std::pair<std::string, std::string> > vars_mat_done;
      for (ga_tree &tree : gis.trees) {
        if (tree.root) {
          if (order == 1) {
            const std::string &name = tree.root->name_test1;
            const std::vector<std::string> vnames_(1,name);
            const std::vector<std::string> &vnames
              = variable_group_exists(name) ? variable_group(name)
                                            : vnames_;
            for (const std::string &vname : vnames) {
              const mesh_fem *mf = associated_mf(vname);
              if (mf && mf->is_reduced() &&
                  vars_vec_done.find(vname) == vars_vec_done.end()) {
                gmm::mult_add(gmm::transposed(mf->extension_matrix()),
                              gmm::sub_vector(unreduced_V,
                                              gis.var_intervals[vname]),
                              gmm::sub_vector(*V,
                                              interval_of_variable(vname)));
                vars_vec_done.insert(vname);
              }
            }
          } else {
            std::string &name1 = tree.root->name_test1;
            std::string &name2 = tree.root->name_test2;
            const std::vector<std::string> vnames1_(1,name1),
                                           vnames2_(2,name2);
            const std::vector<std::string> &vnames1
              = variable_group_exists(name1) ? variable_group(name1)
                                             : vnames1_;
            const std::vector<std::string> &vnames2
              = variable_group_exists(name2) ? variable_group(name2)
                                             : vnames2_;
            for (const std::string &vname1 : vnames1) {
              for (const std::string &vname2 : vnames2) {
                const mesh_fem *mf1 = associated_mf(vname1);
                const mesh_fem *mf2 = associated_mf(vname2);
                if (((mf1 && mf1->is_reduced())
                     || (mf2 && mf2->is_reduced()))) {
                  std::pair<std::string, std::string> p(vname1, vname2);
                  if (vars_mat_done.find(p) == vars_mat_done.end()) {
                    gmm::sub_interval uI1 = gis.var_intervals[vname1];
                    gmm::sub_interval uI2 = gis.var_intervals[vname2];
                    gmm::sub_interval I1 = interval_of_variable(vname1);
                    gmm::sub_interval I2 = interval_of_variable(vname2);
                    if ((mf1 && mf1->is_reduced()) &&
                        (mf2 && mf2->is_reduced())) {
                      model_real_sparse_matrix aux(I1.size(), uI2.size());
                      model_real_row_sparse_matrix M(I1.size(), I2.size());
                      gmm::mult(gmm::transposed(mf1->extension_matrix()),
                                gmm::sub_matrix(unreduced_K, uI1, uI2), aux);
                      gmm::mult(aux, mf2->extension_matrix(), M);
                      gmm::add(M, gmm::sub_matrix(*K, I1, I2));
                    } else if (mf1 && mf1->is_reduced()) {
                      model_real_sparse_matrix M(I1.size(), I2.size());
                      gmm::mult(gmm::transposed(mf1->extension_matrix()),
                                gmm::sub_matrix(unreduced_K, uI1, uI2), M);
                      gmm::add(M, gmm::sub_matrix(*K, I1, I2));
                    } else {
                      model_real_row_sparse_matrix M(I1.size(), I2.size());
                      gmm::mult(gmm::sub_matrix(unreduced_K, uI1, uI2),
                                mf2->extension_matrix(), M);
                      gmm::add(M, gmm::sub_matrix(*K, I1, I2));
                    }
                    vars_mat_done.insert(p);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  void ga_workspace::clear_expressions() {
    trees.clear();
    macro_trees.clear();
  }

  void ga_workspace::print(std::ostream &str) {
    for (size_type i = 0; i < trees.size(); ++i)
      if (trees[i].ptree->root) {
        cout << "Expression tree " << i << " of order " <<
                trees[i].ptree->root->nb_test_functions() << " :" << endl;
        ga_print_node(trees[i].ptree->root, str);
        cout << endl;
      }
  }

  void ga_workspace::tree_description::copy(const tree_description& td) {
    order = td.order;
    interpolation = td.interpolation;
    varname_interpolation = td.varname_interpolation;
    name_test1 = td.name_test1;
    name_test2 = td.name_test2;
    interpolate_name_test1 = td.interpolate_name_test1;
    interpolate_name_test2 = td.interpolate_name_test2;
    mim = td.mim;
    m = td.m;
    rg = td.rg;
    ptree = 0;
    if (td.ptree) ptree = new ga_tree(*(td.ptree));
  }

  ga_workspace::tree_description &ga_workspace::tree_description::operator =
  (const ga_workspace::tree_description& td)
  { if (ptree) delete ptree; ptree = 0; copy(td); return *this; }
  ga_workspace::tree_description::~tree_description()
  { if (ptree) delete ptree; ptree = 0; }



  //=========================================================================
  // Extract the constant term of degree 1 expressions
  //=========================================================================

  std::string ga_workspace::extract_constant_term(const mesh &m) {
    std::string constant_term;
    for (size_type i = 0; i < trees.size(); ++i) {
      ga_workspace::tree_description &td =  trees[i];

      if (td.order == 1) {
        ga_tree local_tree = *(td.ptree);
        if (local_tree.root)
          ga_node_extract_constant_term(local_tree, local_tree.root, *this, m);
        if (local_tree.root)
          ga_semantic_analysis("", local_tree, *this, m.dim(),
                               ref_elt_dim_of_mesh(m), false, false);
        if (local_tree.root && local_tree.root->node_type != GA_NODE_ZERO) {
          constant_term += "-("+ga_tree_to_string(local_tree)+")";
        }
      }
    }
    return constant_term;
  }

  //=========================================================================
  // Extract the order zero term
  //=========================================================================

  std::string ga_workspace::extract_order0_term() {
    std::string term;
    for (size_type i = 0; i < trees.size(); ++i) {
      ga_workspace::tree_description &td =  trees[i];
      if (td.order == 0) {
        ga_tree &local_tree = *(td.ptree);
        if (term.size())
          term += "+("+ga_tree_to_string(local_tree)+")";
        else
          term = "("+ga_tree_to_string(local_tree)+")";
      }
    }
    return term;
  }


  //=========================================================================
  // Extract the order one term corresponding to a certain test function
  //=========================================================================

  std::string ga_workspace::extract_order1_term(const std::string &varname) {
    std::string term;
    for (size_type i = 0; i < trees.size(); ++i) {
      ga_workspace::tree_description &td =  trees[i];
      if (td.order == 1 && td.name_test1.compare(varname) == 0) {
        ga_tree &local_tree = *(td.ptree);
        if (term.size())
          term += "+("+ga_tree_to_string(local_tree)+")";
        else
          term = "("+ga_tree_to_string(local_tree)+")";
      }
    }
    return term;
  }

  //=========================================================================
  // Extract Neumann terms
  //=========================================================================

  std::string ga_workspace::extract_Neumann_term(const std::string &varname) {
    std::string result;
    for (size_type i = 0; i < trees.size(); ++i) {
      ga_workspace::tree_description &td =  trees[i];
      if (td.order == 1 && td.name_test1.compare(varname) == 0) {
        ga_tree &local_tree = *(td.ptree);
        if (local_tree.root)
          ga_extract_Neumann_term(local_tree, varname, *this,
                                      local_tree.root, result);
      }
    }
    return result;
  }

  //=========================================================================
  // Compilation of assembly trees into a list of basic instructions
  //=========================================================================

  static void add_interval_to_gis(const ga_workspace &workspace,
                                  const std::string &varname,
                                  ga_instruction_set &gis) {
    if (workspace.variable_group_exists(varname)) {
      for (const std::string &v : workspace.variable_group(varname))
        add_interval_to_gis(workspace, v, gis);
    } else {
      if (gis.var_intervals.find(varname) == gis.var_intervals.end()) {
	const mesh_fem *mf = workspace.associated_mf(varname);
	size_type nd = mf ? mf->nb_basic_dof() :
	  gmm::vect_size(workspace.value(varname));
	gis.var_intervals[varname]=gmm::sub_interval(gis.nb_dof, nd);
	gis.nb_dof += nd;
      }
      gis.max_dof = std::max(gis.max_dof,
			     workspace.interval_of_variable(varname).last());
    }
  }

  static void extend_variable_in_gis(const ga_workspace &workspace,
                                     const std::string &varname,
                                     ga_instruction_set &gis) {
    if (workspace.variable_group_exists(varname)) {
      for (const std::string &v : workspace.variable_group(varname))
        extend_variable_in_gis(workspace, v, gis);
    } else if (gis.extended_vars.find(varname)==gis.extended_vars.end()) {
      const mesh_fem *mf = workspace.associated_mf(varname);
      if (mf->is_reduced()) {
        auto n = (mf->get_qdim() == 1) ? workspace.qdim(varname) : 1;
        base_vector U(mf->nb_basic_dof() * n);
        mf->extend_vector(workspace.value(varname), U);
        gis.really_extended_vars[varname] = U;
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

  static void ga_compile_node(const pga_tree_node pnode,
                              const ga_workspace &workspace,
                              ga_instruction_set &gis,
                              ga_instruction_set::region_mim_instructions &rmi,
                              const mesh &m, bool function_case,
                              ga_if_hierarchy &if_hierarchy) {

    if (pnode->node_type == GA_NODE_PREDEF_FUNC ||
        pnode->node_type == GA_NODE_OPERATOR ||
        pnode->node_type == GA_NODE_SPEC_FUNC ||
        pnode->node_type == GA_NODE_CONSTANT ||
        pnode->node_type == GA_NODE_ALLINDICES ||
        pnode->node_type == GA_NODE_RESHAPE) return;

    // cout << "compiling "; ga_print_node(pnode, cout); cout << endl;

    pga_instruction pgai;
    ga_if_hierarchy *pif_hierarchy = &if_hierarchy;
    ga_if_hierarchy new_if_hierarchy;

    const mesh_fem *mf1 = 0, *mf2 = 0;
    const mesh_fem **mfg1 = 0, **mfg2 = 0;
    fem_interpolation_context *pctx1 = 0, *pctx2 = 0;
    bool tensor_to_clear = false;

    if (pnode->test_function_type) {
      if (pnode->name_test1.size())
        mf1 = workspace.associated_mf(pnode->name_test1);
      if (mf1) {
        pctx1 = &(gis.ctx);
        const std::string &intn1 = pnode->interpolate_name_test1;
        if (intn1.size()) {
          pctx1 = &(rmi.interpolate_infos[intn1].ctx);
          if (workspace.variable_group_exists(pnode->name_test1)) {
            ga_instruction_set::variable_group_info &vgi =
              rmi.interpolate_infos[intn1].groups_info[pnode->name_test1];
            mfg1 = &(vgi.mf);
            mf1 = 0;
          }
        }
      }
      if (pnode->name_test2.size())
        mf2 = workspace.associated_mf(pnode->name_test2);
      if (mf2) {
        pctx2 = &(gis.ctx);
        const std::string &intn2 = pnode->interpolate_name_test2;
        if (intn2.size()) {
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

    // Produce a resize instruction which is stored if no equivalent node is
    // detected and is the mesh is not uniform.
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
    if (rmi.node_list.find(pnode->hash_value) != rmi.node_list.end()) {
      std::list<pga_tree_node> &node_list = rmi.node_list[pnode->hash_value];
      for (std::list<pga_tree_node>::iterator it = node_list.begin();
           it != node_list.end(); ++it) {
        // cout << "found potential equivalent nodes ";
        // ga_print_node(pnode, cout);
        // cout << " and "; ga_print_node(*it, cout); cout << endl;
        if (sub_tree_are_equal(pnode, *it, workspace, 1)) {
          pnode->t.set_to_copy((*it)->t);
          return;
        }
        if (sub_tree_are_equal(pnode, *it, workspace, 2)) {
          // cout << "confirmed with transpose" << endl;
          if (pnode->nb_test_functions() == 2) {
            if (pgai) { // resize instruction if needed
              if (is_uniform)
                { pgai->exec(); }
              else { rmi.instructions.push_back(std::move(pgai)); }
            }
            pgai = std::make_shared<ga_instruction_transpose_test>
              (pnode->tensor(), (*it)->tensor());
            rmi.instructions.push_back(std::move(pgai));
          } else {
            pnode->t.set_to_copy((*it)->t);
          }
          return;
        }
        std::stringstream ss;
        ss << "Detected wrong equivalent nodes: ";
        ga_print_node(pnode, ss);
        ss << " and "; ga_print_node(*it, ss);
        ss << " (no problem, but hash code would be adapted) " << endl;
        GMM_TRACE2(ss.str());
      }
    }

    if (pgai) { // resize instruction if needed and no equivalent node detected
      if (is_uniform) { pgai->exec(); }
      else {
        if (mfg1 || mfg2)
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
    case GA_NODE_RESHAPE: case GA_NODE_INTERPOLATE_FILTER:
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
      GMM_ASSERT1(!function_case,
                  "No use of Normal is allowed in functions");
      if (pnode->tensor().size() != m.dim())
        pnode->init_vector_tensor(m.dim());
      pgai = std::make_shared<ga_instruction_copy_Normal>
             (pnode->tensor(), gis.Normal);
      rmi.instructions.push_back(std::move(pgai));
      break;

    case GA_NODE_INTERPOLATE_X:
    case GA_NODE_INTERPOLATE_NORMAL:
      GMM_ASSERT1(!function_case,
                  "No use of Interpolate is allowed in functions");
      if (pnode->tensor().size() != m.dim())
        pnode->init_vector_tensor(m.dim());
      if (pnode->node_type == GA_NODE_INTERPOLATE_X)
        pgai = std::make_shared<ga_instruction_copy_small_vect>
               (pnode->tensor(),
                rmi.interpolate_infos[pnode->interpolate_name].pt_y);
      else if (pnode->node_type == GA_NODE_INTERPOLATE_NORMAL)
        pgai = std::make_shared<ga_instruction_copy_Normal>
               (pnode->tensor(),
                rmi.interpolate_infos[pnode->interpolate_name].Normal);
      rmi.instructions.push_back(std::move(pgai));
      break;

    case GA_NODE_VAL: case GA_NODE_GRAD:
    case GA_NODE_HESS: case GA_NODE_DIVERG:
    case GA_NODE_ELEMENTARY_VAL: case GA_NODE_ELEMENTARY_GRAD:
    case GA_NODE_ELEMENTARY_HESS: case GA_NODE_ELEMENTARY_DIVERG:
    case GA_NODE_XFEM_PLUS_VAL: case GA_NODE_XFEM_PLUS_GRAD:
    case GA_NODE_XFEM_PLUS_HESS: case GA_NODE_XFEM_PLUS_DIVERG:
    case GA_NODE_XFEM_MINUS_VAL: case GA_NODE_XFEM_MINUS_GRAD:
    case GA_NODE_XFEM_MINUS_HESS: case GA_NODE_XFEM_MINUS_DIVERG:
      if (function_case) {
        GMM_ASSERT1(pnode->node_type != GA_NODE_ELEMENTARY_VAL &&
                    pnode->node_type != GA_NODE_ELEMENTARY_GRAD &&
                    pnode->node_type != GA_NODE_ELEMENTARY_HESS &&
                    pnode->node_type != GA_NODE_ELEMENTARY_DIVERG,
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
        GMM_ASSERT1(!mf,"No fem expression is allowed in function expression");
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
        const mesh_fem *mf = workspace.associated_mf(pnode->name);
        const im_data *imd = workspace.associated_im_data(pnode->name);

        if (imd) {
          pgai = std::make_shared<ga_instruction_extract_local_im_data>
            (pnode->tensor(), *imd, workspace.value(pnode->name),
             gis.pai, gis.ctx, workspace.qdim(pnode->name));
          rmi.instructions.push_back(std::move(pgai));
        } else {
          GMM_ASSERT1(mf, "Internal error");

          GMM_ASSERT1(&(mf->linked_mesh()) == &(m),
                      "The finite element of variable " << pnode->name <<
                      " has to be defined on the same mesh than the "
                      "integration method or interpolation used");

          // An instruction for extracting local dofs of the variable.
          if (rmi.local_dofs.count(pnode->name) == 0) {
            rmi.local_dofs[pnode->name] = base_vector(1);
            extend_variable_in_gis(workspace, pnode->name, gis);
            // cout << "local dof of " << pnode->name << endl;
            size_type qmult2 = mf->get_qdim();
            if (qmult2 > 1 && !(mf->is_uniformly_vectorized()))
              qmult2 = size_type(-1);
            pgai = std::make_shared<ga_instruction_slice_local_dofs>
              (*mf, *(gis.extended_vars[pnode->name]), gis.ctx,
               rmi.local_dofs[pnode->name],
               workspace.qdim(pnode->name) / mf->get_qdim(), qmult2);
            rmi.elt_instructions.push_back(std::move(pgai));
          }

          // An instruction for pfp update
          if (rmi.pfps.count(mf) == 0) {
            rmi.pfps[mf] = 0;
            pgai = std::make_shared<ga_instruction_update_pfp>
              (*mf, rmi.pfps[mf], gis.ctx, gis.fp_pool);
            if (mf->is_uniform())
              rmi.begin_instructions.push_back(std::move(pgai));
            else
              rmi.instructions.push_back(std::move(pgai));
          }

          // An instruction for the base value
          pgai = pga_instruction();
          switch (pnode->node_type) {
          case GA_NODE_VAL: case GA_NODE_ELEMENTARY_VAL:
            if (rmi.base.count(mf) == 0 ||
               !(if_hierarchy.is_compatible(rmi.base_hierarchy[mf]))) {
              rmi.base_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_val_base>
                (rmi.base[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_PLUS_VAL:
            if (rmi.xfem_plus_base.count(mf) == 0 ||
             !(if_hierarchy.is_compatible(rmi.xfem_plus_base_hierarchy[mf]))) {
              rmi.xfem_plus_base_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_plus_val_base>
                (rmi.xfem_plus_base[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_VAL:
            if (rmi.xfem_minus_base.count(mf) == 0 ||
            !(if_hierarchy.is_compatible(rmi.xfem_minus_base_hierarchy[mf]))) {
              rmi.xfem_minus_base_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_minus_val_base>
                (rmi.xfem_minus_base[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_GRAD: case GA_NODE_DIVERG:
          case GA_NODE_ELEMENTARY_GRAD: case GA_NODE_ELEMENTARY_DIVERG:
            if (rmi.grad.count(mf) == 0 ||
                !(if_hierarchy.is_compatible(rmi.grad_hierarchy[mf]))) {
              rmi.grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_grad_base>
                (rmi.grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_PLUS_GRAD: case GA_NODE_XFEM_PLUS_DIVERG:
            if (rmi.xfem_plus_grad.count(mf) == 0 ||
             !(if_hierarchy.is_compatible(rmi.xfem_plus_grad_hierarchy[mf]))) {
              rmi.xfem_plus_grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_plus_grad_base>
                (rmi.xfem_plus_grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_GRAD: case GA_NODE_XFEM_MINUS_DIVERG:
            if (rmi.xfem_minus_grad.count(mf) == 0 ||
            !(if_hierarchy.is_compatible(rmi.xfem_minus_grad_hierarchy[mf]))) {
              rmi.xfem_minus_grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_minus_grad_base>
                (rmi.xfem_minus_grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_HESS: case GA_NODE_ELEMENTARY_HESS:
            if (rmi.hess.count(mf) == 0 ||
                !(if_hierarchy.is_compatible(rmi.hess_hierarchy[mf]))) {
              rmi.hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_hess_base>
                (rmi.hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_PLUS_HESS:
            if (rmi.xfem_plus_hess.count(mf) == 0 ||
             !(if_hierarchy.is_compatible(rmi.xfem_plus_hess_hierarchy[mf]))) {
              rmi.xfem_plus_hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_plus_hess_base>
                (rmi.xfem_plus_hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_HESS:
            if (rmi.xfem_minus_hess.count(mf) == 0 ||
            !(if_hierarchy.is_compatible(rmi.xfem_minus_hess_hierarchy[mf]))) {
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
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai =
                std::make_shared<ga_instruction_elementary_transformation_val>
                (pnode->tensor(), rmi.base[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_GRAD:
            { // --> t(target_dim*Qmult,N)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai =
                std::make_shared<ga_instruction_elementary_transformation_grad>
                (pnode->tensor(), rmi.grad[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_HESS:
            { // --> t(target_dim*Qmult,N,N)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai =
                std::make_shared<ga_instruction_elementary_transformation_hess>
                (pnode->tensor(), rmi.hess[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_DIVERG:
            { // --> t(1)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai =
               std::make_shared<ga_instruction_elementary_transformation_diverg>
                (pnode->tensor(), rmi.grad[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
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
        const mesh_fem *mf = workspace.associated_mf(pnode->name);
        if (mf) {
          GMM_ASSERT1(&(mf->linked_mesh()) == &(m),
                      "The finite element of variable " << pnode->name <<
                      " and the applied integration method have to be"
                      " defined on the same mesh");

          // An instruction for pfp update
          if (rmi.pfps.count(mf) == 0) {
            rmi.pfps[mf] = 0;
            pgai = std::make_shared<ga_instruction_update_pfp>
              (*mf, rmi.pfps[mf], gis.ctx, gis.fp_pool);
            if (is_uniform)
              rmi.begin_instructions.push_back(std::move(pgai));
            else
              rmi.instructions.push_back(std::move(pgai));
          }

          // An instruction for the base value
          pgai = pga_instruction();
          switch (pnode->node_type) {
          case GA_NODE_VAL_TEST: case GA_NODE_ELEMENTARY_VAL_TEST:
             if (rmi.base.find(mf) == rmi.base.end() ||
                !(if_hierarchy.is_compatible(rmi.base_hierarchy[mf]))) {
              rmi.base_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_val_base>
                (rmi.base[mf], gis.ctx, *mf, rmi.pfps[mf]);
             }
             break;
          case GA_NODE_XFEM_PLUS_VAL_TEST:
            if (rmi.xfem_plus_base.find(mf) == rmi.xfem_plus_base.end() ||
             !(if_hierarchy.is_compatible(rmi.xfem_plus_base_hierarchy[mf]))) {
              rmi.xfem_plus_base_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_plus_val_base>
                (rmi.xfem_plus_base[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_VAL_TEST:
            if (rmi.xfem_minus_base.find(mf) == rmi.xfem_minus_base.end() ||
            !(if_hierarchy.is_compatible(rmi.xfem_minus_base_hierarchy[mf]))) {
              rmi.xfem_minus_base_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_minus_val_base>
                (rmi.xfem_minus_base[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_GRAD_TEST: case GA_NODE_DIVERG_TEST:
          case GA_NODE_ELEMENTARY_GRAD_TEST:
          case GA_NODE_ELEMENTARY_DIVERG_TEST:
            if (rmi.grad.find(mf) == rmi.grad.end() ||
                !(if_hierarchy.is_compatible(rmi.grad_hierarchy[mf]))) {
              rmi.grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_grad_base>
                (rmi.grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_PLUS_GRAD_TEST: case GA_NODE_XFEM_PLUS_DIVERG_TEST:
            if (rmi.xfem_plus_grad.find(mf) == rmi.xfem_plus_grad.end() ||
             !(if_hierarchy.is_compatible(rmi.xfem_plus_grad_hierarchy[mf]))) {
              rmi.xfem_plus_grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_plus_grad_base>
                (rmi.xfem_plus_grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_GRAD_TEST:
          case GA_NODE_XFEM_MINUS_DIVERG_TEST:
            if (rmi.xfem_minus_grad.find(mf) == rmi.xfem_minus_grad.end() ||
            !(if_hierarchy.is_compatible(rmi.xfem_minus_grad_hierarchy[mf]))) {
              rmi.xfem_minus_grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_minus_grad_base>
                (rmi.xfem_minus_grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_HESS_TEST: case GA_NODE_ELEMENTARY_HESS_TEST:
            if (rmi.hess.count(mf) == 0 ||
                !(if_hierarchy.is_compatible(rmi.hess_hierarchy[mf]))) {
              rmi.hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_hess_base>
                (rmi.hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_PLUS_HESS_TEST:
            if (rmi.xfem_plus_hess.count(mf) == 0 ||
                !(if_hierarchy.is_compatible(rmi.xfem_plus_hess_hierarchy[mf]))
                ) {
              rmi.xfem_plus_hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = std::make_shared<ga_instruction_xfem_plus_hess_base>
                (rmi.xfem_plus_hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_HESS_TEST:
            if (rmi.xfem_minus_hess.find(mf) == rmi.xfem_minus_hess.end() ||
            !(if_hierarchy.is_compatible(rmi.xfem_minus_hess_hierarchy[mf]))) {
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
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai =
             std::make_shared<ga_instruction_elementary_transformation_val_base>
                (pnode->tensor(), rmi.base[mf], mf->get_qdim(),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_GRAD_TEST:
            { // --> t(Qmult*ndof,Qmult*target_dim,N)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai =
            std::make_shared<ga_instruction_elementary_transformation_grad_base>
                (pnode->tensor(), rmi.grad[mf], mf->get_qdim(),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_HESS_TEST:
            { // --> t(Qmult*ndof,Qmult*target_dim,N,N)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai =
            std::make_shared<ga_instruction_elementary_transformation_hess_base>
                (pnode->tensor(), rmi.hess[mf], mf->get_qdim(),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_DIVERG_TEST:
            { // --> t(Qmult*ndof)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai =
          std::make_shared<ga_instruction_elementary_transformation_diverg_base>
                (pnode->tensor(), rmi.grad[mf], mf->get_qdim(),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          default: break;
          }
          if (pgai) rmi.instructions.push_back(std::move(pgai));
        }
        add_interval_to_gis(workspace, pnode->name, gis);
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
        add_interval_to_gis(workspace, pnode->name, gis);
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
           size_type tps1 = child0->tensor_proper_size();
           size_type tps2 = child1->tensor_proper_size();
           size_type s1 = (tps1 * tps2) / pnode->tensor_proper_size();
           size_type s2 = size_type(round(sqrt(scalar_type(s1))));

           pgai = pga_instruction();
           if (pnode->op_type == GA_DOT || pnode->op_type == GA_COLON ||
               (pnode->op_type == GA_MULT && dim0 == 4) ||
               (pnode->op_type == GA_MULT && dim1 <= 1) ||
               child0->tensor().size() == 1 || child1->tensor().size() == 1) {

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
               if (child0->tensor_proper_size() == 1) {
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
                   pgai = ga_uniform_instruction_reduction_switch
                     (pnode->t, child0->t, child1->t, s2, tensor_to_clear);
                 else // Unrolled instruction
                   pgai = ga_instruction_reduction_switch
                     (pnode->t, child0->t, child1->t, s2, tensor_to_clear);
               }
             } else {
               if (child1->test_function_type == 1 ||
                   child1->test_function_type == 3) {
                 if (child1->test_function_type == 3 ||
                     child1->tensor_proper_size() <= s2) {
                   if (tps1 == 1) {
                     if (is_uniform) { // Unrolled instruction
                       pgai = ga_uniform_instruction_simple_tmult
                         (pnode->tensor(), child1->tensor(), child0->tensor());
                     } else
                       pgai = std::make_shared<ga_instruction_simple_tmult>
                         (pnode->tensor(), child1->tensor(), child0->tensor());
                   } else if (is_uniform) // Unrolled instruction
                     pgai = ga_uniform_instruction_reduction_switch
                       (pnode->t, child0->t, child1->t, s2, tensor_to_clear);
                   else // Unrolled instruction
                     pgai = ga_instruction_reduction_switch
                       (pnode->t, child0->t, child1->t, s2, tensor_to_clear);
                 } else
                   pgai = std::make_shared<ga_instruction_spec_reduction>
                     (pnode->tensor(), child1->tensor(), child0->tensor(), s2);
               } else if (child1->test_function_type == 0 ||
                          (child0->tensor_proper_size() == s2 &&
                           child1->tensor_proper_size() == s2)) {
                 if (tps1 == 1) {
                   if (is_uniform) { // Unrolled instruction
                     pgai = ga_uniform_instruction_simple_tmult
                       (pnode->tensor(), child0->tensor(), child1->tensor());
                   } else
                     pgai = std::make_shared<ga_instruction_simple_tmult>
                       (pnode->tensor(), child0->tensor(), child1->tensor());
                 } else {
                   if (is_uniform) // Unrolled instruction
                     pgai = ga_uniform_instruction_reduction_switch
                       (pnode->t, child1->t, child0->t, s2, tensor_to_clear);
                   else // Unrolled instruction
                     pgai = ga_instruction_reduction_switch
                       (pnode->t, child1->t, child0->t, s2, tensor_to_clear);
                 }
               } else {
                 if (child0->tensor_proper_size() == s2)
                   pgai = ga_uniform_instruction_reduction_switch
                     (pnode->t, child1->t, child0->t, s2, tensor_to_clear);
                 else if (child1->tensor_proper_size() == s2)
                   pgai = std::make_shared<ga_instruction_spec_reduction>
                     (pnode->tensor(), child0->tensor(), child1->tensor(), s2);
                 else
                   pgai = std::make_shared<ga_instruction_spec2_reduction>
                     (pnode->tensor(), child0->tensor(), child1->tensor(), s2);
               }
             }


           } else { // GA_MULT

             if (pnode->test_function_type < 3) {
               if (child1->tensor_proper_size() == 1) {
                 if (is_uniform) // Unrolled instruction
                   pgai = ga_uniform_instruction_simple_tmult
                     (pnode->tensor(), child1->tensor(), child0->tensor());
                 else
                   pgai = std::make_shared<ga_instruction_simple_tmult>
                     (pnode->tensor(), child1->tensor(), child0->tensor());
               } else if (child0->tensor_proper_size() == 1) {
                 if (is_uniform) // Unrolled instruction
                   pgai = ga_uniform_instruction_simple_tmult
                     (pnode->tensor(), child0->tensor(), child1->tensor());
                 else
                   pgai = std::make_shared<ga_instruction_simple_tmult>
                     (pnode->tensor(), child0->tensor(), child1->tensor());
               } else {
                 if (dim0 == 2)
                   pgai = std::make_shared<ga_instruction_matrix_mult>
                     (pnode->tensor(), child0->tensor(), child1->tensor());
               }
             } else {
               if (child1->tensor_proper_size() == 1) {
                 if (child1->test_function_type == 0 ||
                     child1->test_function_type == 1) {
                   if (is_uniform) // Unrolled instruction
                     pgai = ga_uniform_instruction_simple_tmult
                       (pnode->tensor(), child1->tensor(), child0->tensor());
                   else
                     pgai = std::make_shared<ga_instruction_simple_tmult>
                       (pnode->tensor(), child1->tensor(), child0->tensor());
                 } else
                   pgai = std::make_shared<ga_instruction_spec_tmult>
                     (pnode->tensor(), child0->tensor(), child1->tensor(),
                      child0->tensor_proper_size(),
                      child1->tensor_proper_size());
               } else if (child0->tensor_proper_size() == 1) {
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
                      child1->tensor_proper_size(),
                      child0->tensor_proper_size());
               } else if (dim0 == 2) {
                 if (child1->test_function_type != 2)
                   pgai = std::make_shared<ga_instruction_matrix_mult>
                     (pnode->tensor(), child0->tensor(), child1->tensor());
                 else
                   pgai = std::make_shared<ga_instruction_matrix_mult_spec>
                     (pnode->tensor(), child0->tensor(), child1->tensor());
               }
             }
           }
           GMM_ASSERT1(pgai.get(), "Internal error");
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
         if (pnode->tensor_proper_size() != 1) {
           pgai = std::make_shared<ga_instruction_transpose>
             (pnode->tensor(), child0->tensor());
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
        size_type nbc1 = pnode->nbc1, nbc2 = pnode->nbc2, nbc3 = pnode->nbc3;
        size_type nbl = pnode->children.size() / (nbc1*nbc2*nbc3);
        if (pnode->test_function_type) {
          std::vector<const base_tensor *> components(pnode->children.size());

          if (nbc1 == 1 && nbc2 == 1 && nbc3 == 1) {
            for (size_type i = 0; i < pnode->children.size(); ++i)
              components[i]  = &(pnode->children[i]->tensor());
          } else if (nbc2 == 1 && nbc3 == 1) {
            for (size_type i = 0; i < nbl; ++i)
              for (size_type j = 0; j < nbc1; ++j)
                components[i+j*nbl] = &(pnode->children[i*nbc1+j]->tensor());
          } else {
            size_type n = 0;
            for (size_type i = 0; i < nbl; ++i)
              for (size_type j = 0; j < nbc3; ++j)
                for (size_type k = 0; k < nbc2; ++k)
                  for (size_type l = 0; l < nbc1; ++l)
                    components[i+j*nbl+k*nbl*nbc3+l*nbc2*nbc3*nbl]
                      = &(pnode->children[n++]->tensor());
          }
          pgai = std::make_shared<ga_instruction_c_matrix_with_tests>
            (pnode->tensor(), components);
        } else {
          std::vector<scalar_type *> components(pnode->children.size());
          if (nbc1 == 1 && nbc2 == 1 && nbc3 == 1) {
            for (size_type i = 0; i < pnode->children.size(); ++i)
              components[i]  = &(pnode->children[i]->tensor()[0]);
          } else if (nbc2 == 1 && nbc3 == 1) {
            for (size_type i = 0; i < nbl; ++i)
              for (size_type j = 0; j < nbc1; ++j)
                components[i+j*nbl] = &(pnode->children[i*nbc1+j]->tensor()[0]);
          } else {
            size_type n = 0;
            for (size_type i = 0; i < nbl; ++i)
              for (size_type j = 0; j < nbc3; ++j)
                for (size_type k = 0; k < nbc2; ++k)
                  for (size_type l = 0; l < nbc1; ++l)
                    components[i+j*nbl+k*nbl*nbc3+l*nbc2*nbc3*nbl]
                      = &(pnode->children[n++]->tensor()[0]);
          }
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
        if (pnode->tensor().size() == 1) {
          for (size_type i = 0; i < child0->tensor_order(); ++i)
            mi1[i] = size_type(round(pnode->children[i+1]->tensor()[0])-1);
          pgai = std::make_shared<ga_instruction_copy_scalar>
            (pnode->tensor()[0], child0->tensor()(mi1));
        } else {
          size_type nb_test = pnode->nb_test_functions();
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
  }

  void ga_compile_function(ga_workspace &workspace,
                                  ga_instruction_set &gis, bool scalar) {
    for (size_type i = 0; i < workspace.nb_trees(); ++i) {
      const ga_workspace::tree_description &td = workspace.tree_info(i);

      gis.trees.push_back(*(td.ptree));
      pga_tree_node root = gis.trees.back().root;
      if (root) {
        GMM_ASSERT1(!scalar || (root->tensor().size() == 1),
                    "The result of the given expression is not a scalar");
        ga_instruction_set::region_mim rm(td.mim, td.rg);
        gis.whole_instructions[rm].m = td.m;
        ga_if_hierarchy if_hierarchy;
        ga_compile_node(root, workspace, gis,
                        gis.whole_instructions[rm],*(td.m),true,if_hierarchy);

        gis.coeff = scalar_type(1);
        pga_instruction pgai;
        if (scalar) {
          pgai = std::make_shared<ga_instruction_scalar_assembly>
            (root->tensor(), workspace.assembled_potential(), gis.coeff);

        } else {
          workspace.assembled_tensor() = root->tensor();
          pgai = std::make_shared<ga_instruction_add_to>
            (workspace.assembled_tensor(), root->tensor());
        }
        gis.whole_instructions[rm].instructions.push_back(std::move(pgai));
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
        if (transname.compare("neighbour_elt") == 0) {
          pgai = std::make_shared<ga_instruction_neighbour_transformation_call>
            (workspace, rmi.interpolate_infos[transname],
             workspace.interpolate_transformation(transname), gis.ctx,
             gis.Normal, m, gis.ipt, gis.pai, gis.gp_pool,
             gis.neighbour_corresp);
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

  static void ga_compile_interpolation(ga_workspace &workspace,
                                       ga_instruction_set &gis) {
    gis.transformations.clear();
    gis.whole_instructions.clear();
    for (size_type i = 0; i < workspace.nb_trees(); ++i) {
      const ga_workspace::tree_description &td = workspace.tree_info(i);
      if (td.interpolation > 0) {
        gis.trees.push_back(*(td.ptree));

        // Semantic analysis mainly to evaluate fixed size variables and data
        const mesh *m = td.m;
        GMM_ASSERT1(m, "Internal error");
        ga_semantic_analysis("", gis.trees.back(), workspace, m->dim(),
                             ref_elt_dim_of_mesh(*m), true, false);
        pga_tree_node root = gis.trees.back().root;
        if (root) {
          // Compile tree
          ga_instruction_set::region_mim rm(td.mim, td.rg);
          ga_instruction_set::region_mim_instructions &rmi
            = gis.whole_instructions[rm];
          rmi.m = td.m;
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

  void ga_compile(ga_workspace &workspace,
                         ga_instruction_set &gis, size_type order) {
    gis.transformations.clear();
    gis.whole_instructions.clear();
    for (size_type version : std::array<size_type, 3>{1, 0, 2}) {
      for (size_type i = 0; i < workspace.nb_trees(); ++i) {
        ga_workspace::tree_description &td = workspace.tree_info(i);

        if ((version == td.interpolation) &&
            ((version == 0 && td.order == order) || // Assembly
             ((version > 0 && (td.order == size_type(-1) || // Assignment
                                td.order == size_type(-2) - order))))) {
          ga_tree *added_tree = 0;
          if (td.interpolation) {
            gis.interpolation_trees.push_back(*(td.ptree));
            added_tree = &(gis.interpolation_trees.back());
          } else {
            gis.trees.push_back(*(td.ptree));
            added_tree = &(gis.trees.back());
          }

          // Semantic analysis mainly to evaluate fixed size variables and data
          ga_semantic_analysis("", *added_tree, workspace,
                               td.mim->linked_mesh().dim(),
                               ref_elt_dim_of_mesh(td.mim->linked_mesh()),
                               true, false);
          pga_tree_node root = added_tree->root;
          if (root) {
            // Compile tree
            // cout << "Will compile "; ga_print_node(root, cout); cout << endl;

            ga_instruction_set::region_mim rm(td.mim, td.rg);
            ga_instruction_set::region_mim_instructions &rmi
              = gis.whole_instructions[rm];
            rmi.m = td.m;
            // rmi.interpolate_infos.clear();
            ga_compile_interpolate_trans(root, workspace, gis, rmi, *(td.m));
            ga_compile_node(root, workspace, gis, rmi, *(td.m), false,
                            rmi.current_hierarchy);
            // cout << "compilation finished "; ga_print_node(root, cout);
            // cout << endl;

            if (version > 0) { // Assignment OR interpolation
              if(td.varname_interpolation.size() != 0) {// assignment
                auto *imd
                  = workspace.associated_im_data(td.varname_interpolation);
                auto &V = const_cast<model_real_plain_vector &>
            (workspace.value(td.varname_interpolation));
                GMM_ASSERT1(imd, "Internal error");
                auto pgai = std::make_shared<ga_instruction_assignment>
            (root->tensor(), V, gis.ctx, imd);
                rmi.instructions.push_back(std::move(pgai));
        }
            } else { // assembly
              // Addition of an assembly instruction
              pga_instruction pgai;
              switch(order) {
              case 0:
                pgai = std::make_shared<ga_instruction_scalar_assembly>
                  (root->tensor(), workspace.assembled_potential(), gis.coeff);
                break;
              case 1:
                {
                  const mesh_fem *mf=workspace.associated_mf(root->name_test1);
                  const mesh_fem **mfg = 0;
                  add_interval_to_gis(workspace, root->name_test1, gis);

                  if (mf) {
                    const std::string &intn1 = root->interpolate_name_test1;
                    const gmm::sub_interval *Ir = 0, *In = 0;
                    if (intn1.size() &&
                        workspace.variable_group_exists(root->name_test1)) {
                      ga_instruction_set::variable_group_info &vgi =
                        rmi.interpolate_infos[intn1]
                        .groups_info[root->name_test1];
                      Ir = &(vgi.Ir);
                      In = &(vgi.In);
                      mfg = &(vgi.mf);
                      mf = 0;
                    } else {
                      Ir = &(gis.var_intervals[root->name_test1]);
                      In = &(workspace.interval_of_variable(root->name_test1));
                    }
                    fem_interpolation_context &ctx
                      = intn1.size() ? rmi.interpolate_infos[intn1].ctx
                      : gis.ctx;
                    bool interpolate
                      = (!intn1.empty() && intn1.compare("neighbour_elt")!=0);
                    pgai = std::make_shared<ga_instruction_fem_vector_assembly>
                      (root->tensor(), workspace.unreduced_vector(),
                       workspace.assembled_vector(), ctx, *Ir, *In, mf, mfg,
                       gis.coeff, gis.nbpt, gis.ipt, interpolate);
                  } else {
                    pgai = std::make_shared<ga_instruction_vector_assembly>
                      (root->tensor(), workspace.assembled_vector(),
                       workspace.interval_of_variable(root->name_test1),
                       gis.coeff);
                  }
                }
                break;
              case 2:
                {
                  const mesh_fem *mf1=workspace.associated_mf(root->name_test1);
                  const mesh_fem *mf2=workspace.associated_mf(root->name_test2);
                  const mesh_fem **mfg1 = 0, **mfg2 = 0;
                  const std::string &intn1 = root->interpolate_name_test1;
                  const std::string &intn2 = root->interpolate_name_test2;
                  fem_interpolation_context &ctx1
                    = intn1.empty() ? gis.ctx
                    : rmi.interpolate_infos[intn1].ctx;
                  fem_interpolation_context &ctx2
                    = intn2.empty() ? gis.ctx
                    : rmi.interpolate_infos[intn2].ctx;
                  bool interpolate
                    = (!intn1.empty() && intn1.compare("neighbour_elt")!=0)
                    || (!intn2.empty() && intn2.compare("neighbour_elt")!=0);

                  add_interval_to_gis(workspace, root->name_test1, gis);
                  add_interval_to_gis(workspace, root->name_test2, gis);

                  const gmm::sub_interval *Ir1 = 0, *In1 = 0, *Ir2 = 0, *In2=0;
                  const scalar_type *alpha1 = 0, *alpha2 = 0;

                  if (!intn1.empty() &&
                      workspace.variable_group_exists(root->name_test1)) {
                    ga_instruction_set::variable_group_info &vgi =
                      rmi.interpolate_infos[intn1]
                      .groups_info[root->name_test1];
                    Ir1 = &(vgi.Ir);
                    In1 = &(vgi.In);
                    mfg1 = &(vgi.mf);
                    mf1 = 0;
                    alpha1 = &(vgi.alpha);
                  } else {
                    alpha1 = &(workspace.factor_of_variable(root->name_test1));
                    Ir1 = &(gis.var_intervals[root->name_test1]);
                    In1 = &(workspace.interval_of_variable(root->name_test1));
                  }

                  if (!intn2.empty() &&
                      workspace.variable_group_exists(root->name_test2)) {
                    ga_instruction_set::variable_group_info &vgi =
                      rmi.interpolate_infos[intn2]
                      .groups_info[root->name_test2];
                    Ir2 = &(vgi.Ir);
                    In2 = &(vgi.In);
                    mfg2 = &(vgi.mf);
                    mf2 = 0;
                    alpha2 = &(vgi.alpha);
                  } else {
                    alpha2 = &(workspace.factor_of_variable(root->name_test2));
                    Ir2 = &(gis.var_intervals[root->name_test2]);
                    In2 = &(workspace.interval_of_variable(root->name_test2));
                  }

                  if (!interpolate && mfg1 == 0 && mfg2 == 0 && mf1 && mf2
                      && mf1->get_qdim() == 1 && mf2->get_qdim() == 1
                      && !(mf1->is_reduced()) && !(mf2->is_reduced())) {
                    pgai = std::make_shared
                      <ga_instruction_matrix_assembly_standard_scalar<>>
                      (root->tensor(), workspace.assembled_matrix(), ctx1, ctx2,
                       *In1, *In2, mf1, mf2,
                       gis.coeff, *alpha1, *alpha2, gis.nbpt, gis.ipt);
                  } else if (!interpolate && mfg1 == 0 && mfg2==0 && mf1 && mf2
                             && !(mf1->is_reduced()) && !(mf2->is_reduced())) {
                    if (root->sparsity() == 10 && root->t.qdim()==2)
                      pgai = std::make_shared
                        <ga_instruction_matrix_assembly_standard_vector_opt10_2>
                        (root->tensor(), workspace.assembled_matrix(),ctx1,ctx2,
                         *In1, *In2, mf1, mf2,
                         gis.coeff, *alpha1, *alpha2, gis.nbpt, gis.ipt);
                    else if (root->sparsity() == 10 && root->t.qdim()==3)
                      pgai = std::make_shared
                        <ga_instruction_matrix_assembly_standard_vector_opt10_3>
                        (root->tensor(), workspace.assembled_matrix(),ctx1,ctx2,
                         *In1, *In2, mf1, mf2,
                         gis.coeff, *alpha1, *alpha2, gis.nbpt, gis.ipt);
                    else
                      pgai = std::make_shared
                        <ga_instruction_matrix_assembly_standard_vector<>>
                        (root->tensor(), workspace.assembled_matrix(),ctx1,ctx2,
                         *In1, *In2, mf1, mf2,
                         gis.coeff, *alpha1, *alpha2, gis.nbpt, gis.ipt);

                  } else {
                    pgai = std::make_shared<ga_instruction_matrix_assembly<>>
                      (root->tensor(), workspace.unreduced_matrix(),
                       workspace.assembled_matrix(), ctx1, ctx2,
                       *Ir1, *In1, *Ir2, *In2, mf1, mfg1, mf2, mfg2,
                       gis.coeff, *alpha1, *alpha2, gis.nbpt, gis.ipt,
                       interpolate);
                  }
                  break;
                }
              }
              if (pgai)
                gis.whole_instructions[rm].instructions.push_back
                  (std::move(pgai));
            }
          }
        }
      }
    }
  }


  //=========================================================================
  // Execution of a compiled set of assembly terms
  //=========================================================================


  void ga_function_exec(ga_instruction_set &gis) {

    ga_instruction_set::instructions_set::iterator it
      = gis.whole_instructions.begin();
    for (; it != gis.whole_instructions.end(); ++it) {
      ga_instruction_list &gil = it->second.instructions;
      for (size_type j = 0; j < gil.size(); ++j) j += gil[j]->exec();
    }
  }

  static void ga_interpolation_exec(ga_instruction_set &gis,
                                    ga_workspace &workspace,
                                    ga_interpolation_context &gic) {
    base_matrix G;
    base_small_vector un, up;

    for (const std::string &t : gis.transformations)
      workspace.interpolate_transformation(t)->init(workspace);

    ga_instruction_set::instructions_set::iterator it
      = gis.whole_instructions.begin();
    for (; it != gis.whole_instructions.end(); ++it) {

      const getfem::mesh_im &mim = *(it->first.mim());
      const mesh_region &region = *(it->first.region());
      const getfem::mesh &m = *(it->second.m);
      GMM_ASSERT1(&m == &(gic.linked_mesh()),
                  "Incompatibility of meshes in interpolation");
      ga_instruction_list &gilb = it->second.begin_instructions;
      ga_instruction_list &gile = it->second.elt_instructions;
      ga_instruction_list &gil = it->second.instructions;

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

  static void ga_interpolation_single_point_exec
  (ga_instruction_set &gis, ga_workspace &workspace,
   const fem_interpolation_context &ctx_x, const base_small_vector &Normal,
   const mesh &interp_mesh) {
    gis.ctx = ctx_x;
    gis.Normal = Normal;
    gmm::clear(workspace.assembled_tensor().as_vector());
    gis.nbpt = 1;
    gis.ipt = 0;
    gis.pai = 0;

    ga_instruction_set::instructions_set::iterator
      it = gis.whole_instructions.begin();
    for (; it != gis.whole_instructions.end(); ++it) {
      const getfem::mesh &m = *(it->second.m);
      GMM_ASSERT1(&m == &interp_mesh,
                  "Incompatibility of meshes in interpolation");
      ga_instruction_list &gilb = it->second.begin_instructions;
      for (size_type j = 0; j < gilb.size(); ++j) j += gilb[j]->exec();
      ga_instruction_list &gile = it->second.elt_instructions;
      for (size_type j = 0; j < gile.size(); ++j) j+=gile[j]->exec();
      ga_instruction_list &gil = it->second.instructions;
      for (size_type j = 0; j < gil.size(); ++j) j += gil[j]->exec();
    }
  }

  void ga_exec(ga_instruction_set &gis, ga_workspace &workspace) {
    base_matrix G;
    base_small_vector un;
    scalar_type J(0);

    for (const std::string &t : gis.transformations)
      workspace.interpolate_transformation(t)->init(workspace);

    for (const auto &instr : gis.whole_instructions) {
      const getfem::mesh_im &mim = *(instr.first.mim());
      const getfem::mesh &m = *(instr.second.m);
      GMM_ASSERT1(&m == &(mim.linked_mesh()), "Incompatibility of meshes");
      const ga_instruction_list &gilb = instr.second.begin_instructions;
      const ga_instruction_list &gile = instr.second.elt_instructions;
      const ga_instruction_list &gil = instr.second.instructions;

      // if (gilb.size()) cout << "Begin instructions\n";
      // for (size_type j = 0; j < gilb.size(); ++j)
      //   cout << typeid(*(gilb[j])).name() << endl;
      // if (gile.size()) cout << "\nElement instructions\n";
      // for (size_type j = 0; j < gile.size(); ++j)
      //   cout << typeid(*(gile[j])).name() << endl;
      // cout << "\nGauss pt instructions\n";
      // for (size_type j = 0; j < gil.size(); ++j)
      //   cout << typeid(*(gil[j])).name() << endl;

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
          // cout << "proceed with elt " << v.cv() << " face " << v.f() << endl;
          if (v.cv() != old_cv) {
            pgt = m.trans_of_convex(v.cv());
            pim = mim.int_method_of_element(v.cv());
            m.points_of_convex(v.cv(), G);

            if (pim->type() == IM_NONE) continue;
            GMM_ASSERT1(pim->type() == IM_APPROX, "Sorry, exact methods cannot "
                        "be used in high level generic assembly");
            pai = pim->approx_method();
            pspt = pai->pintegration_points();
            if (pspt->size()) {
              if (pgp && gis.pai == pai && pgt_old == pgt) {
                gis.ctx.change(pgp, 0, 0, G, v.cv(), v.f());
              } else {
                if (pai->is_built_on_the_fly()) {
                  gis.ctx.change(pgt, 0, (*pspt)[0], G, v.cv(), v.f());
                  pgp = 0;
                } else {
                  pgp = gis.gp_pool(pgt, pspt);
                  gis.ctx.change(pgp, 0, 0, G, v.cv(), v.f());
                }
                pgt_old = pgt; gis.pai = pai;
              }
              if (gis.need_elt_size)
                gis.elt_size = convex_radius_estimate(pgt, G)*scalar_type(2);
            }
            old_cv = v.cv();
          } else {
            if (pim->type() == IM_NONE) continue;
            gis.ctx.set_face_num(v.f());
          }
          if (pspt != old_pspt) { first_gp = true; old_pspt = pspt; }
          if (pspt->size()) {
            // iterations on Gauss points
            gis.nbpt = pai->nb_points_on_convex();
            size_type first_ind = 0;
            if (v.f() != short_type(-1)) {
              gis.nbpt = pai->nb_points_on_face(v.f());
              first_ind = pai->ind_first_point_on_face(v.f());
            }
            for (gis.ipt = 0; gis.ipt < gis.nbpt; ++(gis.ipt)) {
	      // cout << "Gauss pt " << gis.ipt << endl;
              if (pgp) gis.ctx.set_ii(first_ind+gis.ipt);
              else gis.ctx.set_xref((*pspt)[first_ind+gis.ipt]);
              if (gis.ipt == 0 || !(pgt->is_linear())) {
                J = gis.ctx.J();
                // Computation of unit normal vector in case of a boundary
                if (v.f() != short_type(-1)) {
                  gis.Normal.resize(G.nrows());
                  un.resize(pgt->dim());
                  gmm::copy(pgt->normals()[v.f()], un);
                  gmm::mult(gis.ctx.B(), un, gis.Normal);
                  scalar_type nup = gmm::vect_norm2(gis.Normal);
                  J *= nup;
                  gmm::scale(gis.Normal, 1.0/nup);
                  gmm::clean(gis.Normal, 1e-13);
                } else gis.Normal.resize(0);
              }
              gis.coeff = J * pai->coeff(first_ind+gis.ipt);
              if (first_gp) {
                for (size_type j = 0; j < gilb.size(); ++j) j+=gilb[j]->exec();
                first_gp = false;
              }
              if (gis.ipt == 0) {
                for (size_type j = 0; j < gile.size(); ++j) j+=gile[j]->exec();
              }
              for (size_type j = 0; j < gil.size(); ++j) j+=gil[j]->exec();
              GA_DEBUG_INFO("");
            }
          }
        }
      }
      GA_DEBUG_INFO("-----------------------------");
    }
    for (const std::string &t : gis.transformations)
      workspace.interpolate_transformation(t)->finalize();
  }

  //=========================================================================
  // User defined functions
  //=========================================================================

  ga_function::ga_function(const ga_workspace &workspace_,
                           const std::string &e)
    : local_workspace(true, workspace_), expr(e), gis(0) {}

  ga_function::ga_function(const model &md, const std::string &e)
    : local_workspace(md), expr(e), gis(0) {}

  ga_function::ga_function(const std::string &e)
    : local_workspace(), expr(e), gis(0) {}

  ga_function::ga_function(const ga_function &gaf)
    : local_workspace(gaf.local_workspace), expr(gaf.expr), gis(0)
  { if (gaf.gis) compile(); }

  void ga_function::compile() const {
    if (gis) delete gis;
    gis = new ga_instruction_set;
    local_workspace.clear_expressions();
    local_workspace.add_function_expression(expr);
    ga_compile_function(local_workspace, *gis, false);
  }

  ga_function &ga_function::operator =(const ga_function &gaf) {
    if (gis) delete gis;
    gis = 0;
    local_workspace = gaf.local_workspace;
    expr = gaf.expr;
    if (gaf.gis) compile();
    return *this;
  }

  ga_function::~ga_function() { if (gis) delete gis; gis = 0; }

  const base_tensor &ga_function::eval() const {
    GMM_ASSERT1(gis, "Uncompiled function");
    gmm::clear(local_workspace.assembled_tensor().as_vector());
    ga_function_exec(*gis);
    return local_workspace.assembled_tensor();
  }

  void ga_function::derivative(const std::string &var) {
    GMM_ASSERT1(gis, "Uncompiled function");
    if (local_workspace.nb_trees()) {
      ga_tree tree = *(local_workspace.tree_info(0).ptree);
      ga_derivative(tree, local_workspace, *((const mesh *)(0)), var, "", 1);
      if (tree.root) {
        ga_semantic_analysis(expr, tree, local_workspace, 1, 1, false, true);
        // To be improved to suppress test functions in the expression ...
        // ga_replace_test_by_cte do not work in all operations like
        // vector components x(1)
        // ga_replace_test_by_cte(tree.root, false);
        // ga_semantic_analysis(expr, tree, local_workspace, 1, 1,
        //                      false, true);
      }
      expr = ga_tree_to_string(tree);
    }
    if (gis) delete gis;
    gis = 0;
    compile();
  }

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

  class  interpolate_transformation_expression
    : public virtual_interpolate_transformation, public context_dependencies {

    struct workspace_gis_pair : public std::pair<ga_workspace, ga_instruction_set> {
      inline ga_workspace &workspace() { return this->first; }
      inline ga_instruction_set &gis() { return this->second; }
    };

    const mesh &source_mesh;
    const mesh &target_mesh;
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
            ga_semantic_analysis(expr, tree, local_workspace, 1, 1,
                                 false, true);
          ga_compile_interpolation(pwi.workspace(), pwi.gis());
        }
      }

      // Element_boxes update (if necessary)
      if (recompute_elt_boxes) {

        element_boxes.clear();
        base_node bmin(N), bmax(N);
        for (dal::bv_visitor cv(target_mesh.convex_index());
             !cv.finished(); ++cv) {

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
                  "Wrong dimension of the tranformation expression");
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
      if (compute_derivatives) { // To be tested both with the computation
                                 // of derivative. Could be optimized ?
        for (auto &&d : derivatives) {
          workspace_gis_pair &pwi = compiled_derivatives[d.first];

          gmm::clear(pwi.workspace().assembled_tensor().as_vector());
          ga_function_exec(pwi.gis());
          d.second = pwi.workspace().assembled_tensor();
        }
      }
      return ret_type;
    }

    interpolate_transformation_expression(const mesh &sm, const mesh &tm,
                                          const std::string &expr_)
      : source_mesh(sm), target_mesh(tm), expr(expr_),
        recompute_elt_boxes(true), extract_variable_done(false),
        extract_data_done(false)
    { this->add_dependency(tm); }

  };


  void add_interpolate_transformation_from_expression
  (model &md, const std::string &name, const mesh &sm, const mesh &tm,
   const std::string &expr) {
    pinterpolate_transformation
      p = std::make_shared<interpolate_transformation_expression>(sm, tm, expr);
    md.add_interpolate_transformation(name, p);
  }

  void add_interpolate_transformation_from_expression
  (ga_workspace &workspace, const std::string &name, const mesh &sm,
   const mesh &tm, const std::string &expr) {
    pinterpolate_transformation
      p = std::make_shared<interpolate_transformation_expression>(sm, tm, expr);
    workspace.add_interpolate_transformation(name, p);
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
        bool is_in = gic.invert(ctx_x.xreal(), P_ref, converged, 1E-4);
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


  pinterpolate_transformation interpolate_transformation_neighbour_instance()
  {
    pinterpolate_transformation
      p = std::make_shared<interpolate_transformation_neighbour>();
    return p;
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

} /* end of namespace */
