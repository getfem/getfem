/*===========================================================================

 Copyright (C) 2000-2020 Yves Renard

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


#include <deque>
#include "getfem/dal_singleton.h"
#include "getfem/getfem_fem.h"
#include "getfem/getfem_mat_elem.h"
#include "getfem/getfem_omp.h"

namespace getfem {
  /* ********************************************************************* */
  /*       Elementary matrices computation.                                */
  /* ********************************************************************* */

  struct emelem_comp_key_  : virtual public dal::static_stored_object_key {
    pmat_elem_type pmt;
    pintegration_method ppi;
    bgeot::pgeometric_trans pgt;
    /* prefer_comp_on_real_element: compute elementary matrices on the real
       element if possible (i.e. if no exact integration is used); this allow
       using inline reduction during the integration */
    bool prefer_comp_on_real_element;
    bool compare(const static_stored_object_key &oo) const override{
      auto &o = dynamic_cast<const emelem_comp_key_ &>(oo);
      if (pmt < o.pmt) return true;
      if (o.pmt < pmt) return false;
      if (ppi < o.ppi) return true;
      if (o.ppi < ppi) return false;
      if (pgt < o.pgt) return true;
      if (o.pgt < pgt) return false;
      return prefer_comp_on_real_element < o.prefer_comp_on_real_element;
    }
    bool equal(const static_stored_object_key &oo) const override{
      auto &o = dynamic_cast<const emelem_comp_key_ &>(oo);

      if (pmt == o.pmt && ppi == o.ppi && pgt == o.pgt) return true;

      auto pmat_key = dal::key_of_stored_object(pmt);
      auto poo_mat_key = dal::key_of_stored_object(o.pmt);
      if (*pmat_key != *poo_mat_key) return false;

      auto pint_key = dal::key_of_stored_object(ppi);
      auto poo_int_key = dal::key_of_stored_object(o.ppi);
      if (*pint_key != *poo_int_key) return false;

      auto pgt_key = dal::key_of_stored_object(pgt);
      auto poo_gt_key = dal::key_of_stored_object(o.pgt);
      if (*pgt_key != *poo_gt_key) return false;

      return true;
    }
    emelem_comp_key_(pmat_elem_type pm, pintegration_method pi,
                       bgeot::pgeometric_trans pg, bool on_relt)
    { pmt = pm; ppi = pi; pgt = pg; prefer_comp_on_real_element = on_relt; }
    emelem_comp_key_(void) { }
  };

  struct emelem_comp_structure_ : public mat_elem_computation {
    bgeot::pgeotrans_precomp pgp;
    ppoly_integration ppi;
    papprox_integration pai;
    bool is_ppi;
    mutable std::vector<base_tensor> mref;
    mutable std::vector<pfem_precomp> pfp;
    mutable std::vector<base_tensor> elmt_stored;
    short_type nbf, dim;
    std::deque<short_type> grad_reduction, hess_reduction, trans_reduction;
    std::deque<short_type> K_reduction;
    std::deque<pfem> trans_reduction_pfi;
    mutable base_vector un, up;
    mutable bool faces_computed;
    mutable bool volume_computed;
    bool is_linear;
    bool computed_on_real_element;
    size_type memsize() const {
      size_type sz = sizeof(emelem_comp_structure_) +
        mref.capacity()*sizeof(base_tensor) +
        grad_reduction.size()*sizeof(short_type) +
        K_reduction.size()*sizeof(short_type) +
        hess_reduction.size()*sizeof(short_type) +
        trans_reduction.size()*sizeof(short_type) +
        trans_reduction_pfi.size()*sizeof(pfem);

      for (size_type i=0; i < mref.size(); ++i) sz += mref[i].memsize();
      return sz;
    }

    emelem_comp_structure_(pmat_elem_type pm, pintegration_method pi,
                           bgeot::pgeometric_trans pg,
                           bool prefer_comp_on_real_element) {

      pgt = pg;
      pgp = bgeot::geotrans_precomp(pg, pi->pintegration_points(), pi);
      pme = pm;
      switch (pi->type()) {
      case IM_EXACT:
        ppi = pi->exact_method(); pai = 0;  is_ppi = true; break;
      case IM_APPROX:
        ppi = 0; pai = pi->approx_method(); is_ppi = false; break;
      case IM_NONE:
        GMM_ASSERT1(false, "Attempt to use IM_NONE integration method "
                    "in assembly!\n");
      }

      faces_computed = volume_computed = false;
      is_linear = pgt->is_linear();
      computed_on_real_element = !is_linear || (prefer_comp_on_real_element && !is_ppi);
      // computed_on_real_element = true;
      nbf = pgt->structure()->nb_faces();
      dim = pgt->structure()->dim();
      mat_elem_type::const_iterator it = pme->begin(), ite = pme->end();
      //      size_type d = pgt->dim();

      for (short_type k = 0; it != ite; ++it, ++k) {
        if ((*it).pfi) {
          if ((*it).pfi->is_on_real_element()) computed_on_real_element = true;
          GMM_ASSERT1(!is_ppi || (((*it).pfi->is_polynomial()) && is_linear
                                  && !computed_on_real_element),
                      "Exact integration not allowed in this context");

          if ((*it).t != GETFEM_NONLINEAR_ && !((*it).pfi->is_equivalent())) {
            // TODO : le numero d'indice à reduire peut changer ...
            trans_reduction.push_back(k);
            trans_reduction_pfi.push_back((*it).pfi);
          }
        }
        switch ((*it).t) {
          case GETFEM_BASE_    :
            if ((*it).pfi->target_dim() > 1) {
              ++k;
              switch((*it).pfi->vectorial_type()) {
              case virtual_fem::VECTORIAL_PRIMAL_TYPE:
                K_reduction.push_back(k); break;
              case virtual_fem::VECTORIAL_DUAL_TYPE:
                grad_reduction.push_back(k); break; // reduction with B
              default: break;
              }
            }
            break;
          case GETFEM_UNIT_NORMAL_ :
            computed_on_real_element = true; break;
          case GETFEM_GRAD_GEOTRANS_ :
          case GETFEM_GRAD_GEOTRANS_INV_ :
            ++k; computed_on_real_element = true; break;
          case GETFEM_GRAD_    : {
            ++k;
            switch((*it).pfi->vectorial_type()) {
            case virtual_fem::VECTORIAL_PRIMAL_TYPE:
              K_reduction.push_back(k); break;
            case virtual_fem::VECTORIAL_DUAL_TYPE:
              grad_reduction.push_back(k); break; // reduction with B
            default: break;
            }
            if ((*it).pfi->target_dim() > 1) ++k;
            if (!((*it).pfi->is_on_real_element()))
              grad_reduction.push_back(k);
          } break;
          case GETFEM_HESSIAN_ : {
            ++k;
            switch((*it).pfi->vectorial_type()) {
            case virtual_fem::VECTORIAL_PRIMAL_TYPE:
              K_reduction.push_back(k); break;
            case virtual_fem::VECTORIAL_DUAL_TYPE:
              grad_reduction.push_back(k); break;
            default: break;
            }

            if ((*it).pfi->target_dim() > 1) ++k;
            if (!((*it).pfi->is_on_real_element()))
              hess_reduction.push_back(k);
          } break;
          case GETFEM_NONLINEAR_ : {
            if ((*it).nl_part == 0) {
              k = short_type(k+(*it).nlt->sizes(size_type(-1)).size()-1);
              GMM_ASSERT1(!is_ppi, "For nonlinear terms you have "
                          "to use approximated integration");
              computed_on_real_element = true;
            }
          } break;
        }
      }

      if (!is_ppi) {
        pfp.resize(pme->size());
        it = pme->begin(), ite = pme->end();
        for (size_type k = 0; it != ite; ++it, ++k)
          if ((*it).pfi)
            pfp[k] = fem_precomp((*it).pfi, pai->pintegration_points(), pi);
          else pfp[k] = 0;
        elmt_stored.resize(pme->size());
      }
      if (!computed_on_real_element) mref.resize(nbf + 1);
    }

    void add_elem(base_tensor &t, fem_interpolation_context& ctx,
                  scalar_type J, bool first, bool trans,
                  mat_elem_integration_callback *icb,
                  bgeot::multi_index sizes) const {
      mat_elem_type::const_iterator it = pme->begin(), ite = pme->end();
      bgeot::multi_index aux_ind;

      for (size_type k = 0; it != ite; ++it, ++k) {
        if ((*it).t == GETFEM_NONLINEAR_)
          (*it).nlt->term_num() = size_type(-1);
      }
      it = pme->begin();

      // incrementing "mit" should match increments of "j" in mat_elem_type::sizes
      bgeot::multi_index::iterator mit = sizes.begin();
      for (size_type k = 0; it != ite; ++it, ++k, ++mit) {
        if (pfp[k]) ctx.set_pfp(pfp[k]);

        switch ((*it).t) {
          case GETFEM_BASE_    :
            if ((*it).pfi && (*it).pfi->target_dim() > 1) ++mit;
            if (trans)
              (*it).pfi->real_base_value(ctx, elmt_stored[k], icb != 0);
            else
              elmt_stored[k] = pfp[k]->val(ctx.ii());
            break;
          case GETFEM_GRAD_    :
            ++mit;
            if ((*it).pfi && (*it).pfi->target_dim() > 1) ++mit;
            if (trans) {
              (*it).pfi->real_grad_base_value(ctx, elmt_stored[k], icb != 0);
              *mit = short_type(ctx.N());
            }
            else
              elmt_stored[k] = pfp[k]->grad(ctx.ii());
            break;
          case GETFEM_HESSIAN_ :
            ++mit;
            if ((*it).pfi && (*it).pfi->target_dim() > 1) ++mit;
            if (trans) {
              (*it).pfi->real_hess_base_value(ctx, elmt_stored[k], icb != 0);
              *mit = short_type(gmm::sqr(ctx.N()));
            }
            else {
              base_tensor tt = pfp[k]->hess(ctx.ii());
              aux_ind.resize(3);
              aux_ind[2] = gmm::sqr(tt.sizes()[2]); aux_ind[1] = tt.sizes()[1];
              aux_ind[0] = tt.sizes()[0];
              tt.adjust_sizes(aux_ind);
              elmt_stored[k] = tt;
            }
            break;
          case GETFEM_UNIT_NORMAL_ :
            *mit = short_type(ctx.N());
            {
              aux_ind.resize(1); aux_ind[0] = short_type(ctx.N());
              elmt_stored[k].adjust_sizes(aux_ind);
            }
            std::copy(up.begin(), up.end(), elmt_stored[k].begin());
            break;
          case GETFEM_GRAD_GEOTRANS_ :
          case GETFEM_GRAD_GEOTRANS_INV_ : {
            size_type P = gmm::mat_ncols(ctx.K()), N=ctx.N();
            base_matrix Bt;
            if (it->t == GETFEM_GRAD_GEOTRANS_INV_) {
              Bt.resize(P,N); gmm::copy(gmm::transposed(ctx.B()),Bt);
            }
            const base_matrix &A = (it->t==GETFEM_GRAD_GEOTRANS_) ? ctx.K():Bt;
            aux_ind.resize(2);
            *mit++ = aux_ind[0] = short_type(gmm::mat_nrows(A));
            *mit = aux_ind[1] = short_type(gmm::mat_ncols(A));
            elmt_stored[k].adjust_sizes(aux_ind);
            std::copy(A.begin(), A.end(), elmt_stored[k].begin());
          } break;
          case GETFEM_NONLINEAR_ :
            if ((*it).nl_part != 0) { /* for auxiliary fem of nonlinear_term,*/
              /* the "prepare" method is called           */
              if ((*it).nlt->term_num() == size_type(-1)) {
                (*it).nlt->prepare(ctx, (*it).nl_part);
                /* the dummy assistant multiplies everybody by 1
                   -> not efficient ! */
              }
              aux_ind.resize(1); aux_ind[0] = 1;
              elmt_stored[k].adjust_sizes(aux_ind); elmt_stored[k][0] = 1.;
            } else {
              if ((*it).nlt->term_num() == size_type(-1)) {
                const bgeot::multi_index &nltsizes
                  = (*it).nlt->sizes(ctx.convex_num());
                elmt_stored[k].adjust_sizes(nltsizes);
                (*it).nlt->compute(ctx, elmt_stored[k]);
                (*it).nlt->term_num() = k;
                for (dim_type ii = 0; ii < nltsizes.size(); ++ii)
                  *mit++ = nltsizes[ii];
                --mit;
              } else {
                elmt_stored[k] = elmt_stored[(*it).nlt->term_num()];
                const bgeot::multi_index &nltsizes = elmt_stored[k].sizes();
                for (dim_type ii = 0; ii < nltsizes.size(); ++ii)
                  *mit++ = nltsizes[ii];
                --mit;
              }
            }
            break;
        }
      }

      GMM_ASSERT1(mit == sizes.end(), "internal error");

      //expand_product_old(t,J*pai->coeff(ctx.ii()), first);
      scalar_type c = J*pai->coeff(ctx.ii());
      if (!icb) {
        if (first) { t.adjust_sizes(sizes); }
        expand_product_daxpy(t, c, first);
      } else {
        icb->eltm.resize(0);
        for (unsigned k=0; k != pme->size(); ++k) {
          if (icb && !((*pme)[k].t == GETFEM_NONLINEAR_
                       && (*pme)[k].nl_part != 0))
            icb->eltm.push_back(&elmt_stored[k]);
        }
        icb->exec(t, first, c);
      }
    }


    void expand_product_old(base_tensor &t, scalar_type J, bool first) const {
      scalar_type V;
      size_type k;
      if (first) std::fill(t.begin(), t.end(), 0.0);
      base_tensor::iterator pt = t.begin();
      std::vector<base_tensor::const_iterator> pts(pme->size());
      std::vector<scalar_type> Vtab(pme->size());
      for (k = 0; k < pme->size(); ++k)
        pts[k] = elmt_stored[k].begin();

      size_type k0 = 0;
      unsigned n0 = unsigned(elmt_stored[0].size());
      /*while (elmt_stored[k0].size() == 1 && k0+1 < pme->size()) {
        J *= elmt_stored[k0][0];
        ++k0; n0 = elmt_stored[k0].size();
        }*/
      base_tensor::const_iterator pts0 = pts[k0];


      k = pme->size()-1; Vtab[k] = J;
      /* very heavy expansion .. takes much time */
      do {
        for (V = Vtab[k]; k!=k0; --k)
          Vtab[k-1] = V = *pts[k] * V;
        for (k=0; k < n0; ++k)
          *pt++ += V * pts0[k];
        for (k=k0+1; k != pme->size() && ++pts[k] == elmt_stored[k].end(); ++k)
          pts[k] = elmt_stored[k].begin();
      } while (k != pme->size());
      GMM_ASSERT1(pt == t.end(), "Internal error");
    }

    /* do the tensorial product using the blas function daxpy (much more
       efficient than a loop).

       efficiency is maximized when the first tensor has a large dimension
     */
    void expand_product_daxpy(base_tensor &t, scalar_type J, bool first)const {
      size_type k;
      base_tensor::iterator pt = t.begin();
      THREAD_SAFE_STATIC std::vector<base_tensor::const_iterator> pts;
      THREAD_SAFE_STATIC std::vector<base_tensor::const_iterator> es_beg;
      THREAD_SAFE_STATIC std::vector<base_tensor::const_iterator> es_end;
      THREAD_SAFE_STATIC std::vector<scalar_type> Vtab;

      pts.resize(0); pts.resize(pme->size()); // resize(0) necessary, do not remove
      es_beg.resize(0); es_beg.resize(pme->size());
      es_end.resize(0); es_end.resize(pme->size());
      Vtab.resize(pme->size());
      size_type nm = 0;
      if (first) memset(&(*t.begin()), 0, t.size()*sizeof(*t.begin())); //std::fill(t.begin(), t.end(), 0.0);
      for (k = 0, nm = 0; k < pme->size(); ++k) {
        if (elmt_stored[k].size() != 1) {
          es_beg[nm] = elmt_stored[k].begin();
          es_end[nm] = elmt_stored[k].end();
          pts[nm] = elmt_stored[k].begin();
          ++nm;
        } else J *= elmt_stored[k][0];
      }
      if (nm == 0) {
        t[0] += J;
      } else {
        BLAS_INT n0 = BLAS_INT(es_end[0] - es_beg[0]);
        base_tensor::const_iterator pts0 = pts[0];

        /* very heavy reduction .. takes much time */
        k = nm-1; Vtab[k] = J;
        BLAS_INT one = BLAS_INT(1);
        scalar_type V;
        do {
          for (V = Vtab[k]; k; --k)
            Vtab[k-1] = V = *pts[k] * V;
          GMM_ASSERT1(pt+n0 <= t.end(), "Internal error");
	  gmm::daxpy_(&n0, &V, const_cast<double*>(&(pts0[0])), &one,
		      (double*)&(*pt), &one);
          pt+=n0;
          for (k=1; k != nm && ++pts[k] == es_end[k]; ++k)
            pts[k] = es_beg[k];
        } while (k != nm);
        GMM_ASSERT1(pt == t.end(), "Internal error");
      }
    }


    void pre_tensors_for_linear_trans(bool volumic) const {

      if ((volumic && volume_computed) || (!volumic && faces_computed)) return;
      // scalar_type exectime = ftool::uclock_sec();

      bgeot::multi_index sizes = pme->sizes(0), mi(sizes.size());
      bgeot::multi_index::iterator mit = sizes.begin(), mite = sizes.end();
      size_type f = 1;
      for ( ; mit != mite; ++mit, ++f) f *= *mit;
      if (f > 1000000)
        GMM_WARNING2("Warning, very large elementary computations.\n"
                    << "Be sure you need to compute this elementary matrix.\n"
                    << "(sizes = " << sizes << " )\n");

      base_tensor aux(sizes);
      std::fill(aux.begin(), aux.end(), 0.0);
      if (volumic) {
        volume_computed = true;
        mref[0] = aux;
      }
      else {
        faces_computed = true;
        std::fill(mref.begin()+1, mref.end(), aux);
      }

      if (is_ppi) // pour accelerer, il faudrait précalculer les dérivées
      {
        base_poly P(dim, 0), Q(dim, 0), R(dim, 0);
        size_type old_ind = size_type(-1), ind;
        for ( ; !mi.finished(sizes); mi.incrementation(sizes)) {

          mat_elem_type::const_iterator it = pme->begin(), ite = pme->end();
          mit = mi.begin();

          ind = *mit; ++mit;

          if ((*it).pfi) {
            if ((*it).pfi->target_dim() > 1)
              { ind += (*it).pfi->nb_base(0) * (*mit); ++mit; }

            Q = ((ppolyfem)((*it).pfi).get())->base()[ind];
          }

          switch ((*it).t) {
          case GETFEM_GRAD_    : Q.derivative(short_type(*mit)); ++mit; break;
          case GETFEM_HESSIAN_ :
            Q.derivative(short_type(*mit % dim));
            Q.derivative(short_type(*mit / dim));
            ++mit; break;
          case GETFEM_BASE_ : break;
          case GETFEM_GRAD_GEOTRANS_:
          case GETFEM_GRAD_GEOTRANS_INV_:
          case GETFEM_UNIT_NORMAL_ :
          case GETFEM_NONLINEAR_ :
            GMM_ASSERT1(false,
                        "Normals, gradients of geotrans and non linear "
                        "terms are not compatible with exact integration, "
                        "use an approximate method instead");
          }
          ++it;

          if (it != ite && *mit != old_ind) {
            old_ind = *mit;
            P.one();
            for (; it != ite; ++it) {
              ind = *mit; ++mit;

              if ((*it).pfi->target_dim() > 1)
                { ind += (*it).pfi->nb_base(0) * (*mit); ++mit; }
              R = ((ppolyfem)((*it).pfi).get())->base()[ind];

              switch ((*it).t) {
              case GETFEM_GRAD_    :
                R.derivative(short_type(*mit)); ++mit;
                break;
              case GETFEM_HESSIAN_ :
                R.derivative(short_type(*mit % dim));
                R.derivative(short_type(*mit / dim));
                ++mit; break;
              case GETFEM_BASE_ : break;
              case GETFEM_UNIT_NORMAL_ :
              case GETFEM_GRAD_GEOTRANS_:
              case GETFEM_GRAD_GEOTRANS_INV_ :
              case GETFEM_NONLINEAR_ :
                GMM_ASSERT1(false, "No nonlinear term allowed here");
              }
              P *= R;
            }
          }
          R = P * Q;
          if (volumic) mref[0](mi) = bgeot::to_scalar(ppi->int_poly(R));
          for (f = 0; f < nbf && !volumic; ++f)
            mref[f+1](mi) = bgeot::to_scalar(ppi->int_poly_on_face(R, short_type(f)));
        }
      }
      else {
        bool first = true;
        fem_interpolation_context ctx;
        size_type ind_l = 0, nb_ptc = pai->nb_points_on_convex(),
          nb_pt_l = nb_ptc, nb_pt_tot =(volumic ? nb_ptc : pai->nb_points());
        for (size_type ip = (volumic ? 0:nb_ptc); ip < nb_pt_tot; ++ip) {
          while (ip == nb_pt_l && ind_l < nbf)
            { nb_pt_l += pai->nb_points_on_face(short_type(ind_l)); ind_l++; }
          ctx.set_ii(ip);
          add_elem(mref[ind_l], ctx, 1.0, first, false, NULL, sizes);
          first = false;
        }
      }
      // cout << "precompute Mat elem computation time : "
      //   << ftool::uclock_sec() - exectime << endl;
    }


    void compute(base_tensor &t, const base_matrix &G, short_type ir,
                 size_type elt, mat_elem_integration_callback *icb = 0) const {
      dim_type P = dim_type(dim), N = dim_type(G.nrows());
      short_type NP = short_type(pgt->nb_points());
      fem_interpolation_context ctx(pgp, 0, 0, G, elt,
				    short_type(ir-1));

      GMM_ASSERT1(G.ncols() == NP, "dimensions mismatch");
      if (ir > 0) {
        up.resize(N); un.resize(P);
        //un = pgt->normals()[ir-1];
        gmm::copy(pgt->normals()[ir-1],un);
      }
      base_tensor taux;
      bool flag = false;

      if (!computed_on_real_element) {
        pre_tensors_for_linear_trans(ir == 0);
        const base_matrix& B = ctx.B(); // compute B and J
        scalar_type J=ctx.J();
        if (ir > 0) {
          gmm::mult(B, un, up);
          scalar_type nup = gmm::vect_norm2(up);
          J *= nup; //up /= nup;
          gmm::scale(up,1.0/nup);
        }

        t = mref[ir]; gmm::scale(t.as_vector(), J);

        if (grad_reduction.size() > 0) {
          std::deque<short_type>::const_iterator it = grad_reduction.begin(),
            ite = grad_reduction.end();
          for ( ; it != ite; ++it) {
            (flag ? t:taux).mat_transp_reduction(flag ? taux:t, B, *it);
            flag = !flag;
          }
        }

        if (K_reduction.size() > 0) {
          std::deque<short_type>::const_iterator it = K_reduction.begin(),
            ite = K_reduction.end();
          for ( ; it != ite; ++it) {
            (flag ? t:taux).mat_transp_reduction(flag ? taux:t, ctx.K(), *it);
            // (flag ? t:taux).mat_transp_reduction(flag ? taux:t, B, *it);
            flag = !flag;
          }
        }

        if (hess_reduction.size() > 0) {
          std::deque<short_type>::const_iterator it = hess_reduction.begin(),
            ite = hess_reduction.end();
          for (short_type l = 1; it != ite; ++it, l = short_type(l*2)) {
            (flag ? t:taux).mat_transp_reduction(flag ? taux:t, ctx.B3(), *it);
            flag = !flag;
          }
        }

      } else { // non linear transformation and methods defined on real elements
        bgeot::multi_index sizes = pme->sizes(elt);

        bool first = true;
        for (size_type ip=(ir == 0) ? 0 : pai->repart()[ir-1];
             ip < pai->repart()[ir]; ++ip, first = false) {
          ctx.set_ii(ip);
          const base_matrix& B = ctx.B(); // J computed as side-effect
          scalar_type J = ctx.J();
          if (ir > 0) {
            gmm::mult(B, un, up);
            scalar_type nup = gmm::vect_norm2(up);
            J *= nup; /*up /= nup;*/gmm::scale(up,1.0/nup);
          }
          add_elem(t, ctx, J, first, true, icb, sizes);
        }

        // GMM_ASSERT1(!first, "No integration point on this element.");
        if (first) {
          GMM_WARNING3("No integration point on this element. "
                       "Caution, returning a null tensor");
          t.adjust_sizes(sizes); gmm::clear(t.as_vector());
        }
      }

      /* Applying linear transformation for non tau-equivalent elements.   */

      if (trans_reduction.size() > 0 && !icb) {
        std::deque<short_type>::const_iterator it = trans_reduction.begin(),
          ite = trans_reduction.end();
        std::deque<pfem>::const_iterator iti = trans_reduction_pfi.begin();
        for ( ; it != ite; ++it, ++iti) {
          ctx.set_pf(*iti); // cout << "M = " << ctx.M() << endl;
          (flag ? t:taux).mat_transp_reduction(flag ? taux:t, ctx.M(), *it);
          flag = !flag;
        }
      }
      if (flag) t = taux;
    }

    void compute(base_tensor &t, const base_matrix &G, size_type elt,
                 mat_elem_integration_callback *icb) const
    { compute(t, G, 0, elt, icb); }

    void compute_on_face(base_tensor &t, const base_matrix &G,
                         short_type f, size_type elt,
                         mat_elem_integration_callback *icb) const
    { compute(t, G, short_type(f+1), elt, icb); }
  };

  pmat_elem_computation mat_elem(pmat_elem_type pm, pintegration_method pi,
                                 bgeot::pgeometric_trans pg,
                                 bool prefer_comp_on_real_element) {
    dal::pstatic_stored_object_key
      pk = std::make_shared<emelem_comp_key_>(pm, pi, pg,
					      prefer_comp_on_real_element);
    dal::pstatic_stored_object o = dal::search_stored_object(pk);
    if (o) return std::dynamic_pointer_cast<const mat_elem_computation>(o);
    pmat_elem_computation
      p = std::make_shared<emelem_comp_structure_>
      (pm, pi, pg, prefer_comp_on_real_element);
    dal::add_stored_object(pk, p, pm, pi, pg);
    return p;
  }


}  /* end of namespace getfem.                                            */

