/*===========================================================================
 
 Copyright (C) 2006-2015 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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

#include "getfem/getfem_partial_mesh_fem.h"

namespace getfem {


  void partial_mesh_fem::clear(void)
  { mesh_fem::clear(); is_adapted = false; }

  partial_mesh_fem::partial_mesh_fem(const mesh_fem &mef)
    : mesh_fem(mef.linked_mesh()), mf(mef)
  { is_adapted = false; }

  static getfem::mesh void_mesh__;

  partial_mesh_fem::partial_mesh_fem(const mesh_fem *mef)
    : mesh_fem(*(mef ? &(mef->linked_mesh()) : &(void_mesh__))), mf(*mef)
  { is_adapted = false; }

  DAL_SIMPLE_KEY(special_partialmf_key, pfem);
  void partial_mesh_fem::adapt(const dal::bit_vector &kept_dofs,
                               const dal::bit_vector &rejected_elt) {
    mf.context_check();

    if (!(mi.is_equal(mf.get_qdims()))) {
      mi = mf.get_qdims();
      Qdim = mf.get_qdim();
      dof_enumeration_made = false; touch(); v_num = act_counter();
    }

    fe_convex = mf.convex_index();
    fe_convex.setminus(rejected_elt);

    gmm::row_matrix<gmm::rsvector<scalar_type> >
      RR(kept_dofs.card(), mf.nb_dof());
    size_type j = 0;
    for (dal::bv_visitor i(kept_dofs); !i.finished(); ++i, ++j)
      RR(j, i) = scalar_type(1);

    R_ = REDUCTION_MATRIX(kept_dofs.card(), mf.nb_basic_dof());
    E_ = EXTENSION_MATRIX(mf.nb_basic_dof(), kept_dofs.card());

    if (mf.is_reduced()) {
      gmm::row_matrix<gmm::rsvector<scalar_type> >
        A(kept_dofs.card(), mf.nb_basic_dof());
      gmm::mult(RR, mf.reduction_matrix(), A);
      gmm::copy(A, R_);
      gmm::row_matrix<gmm::rsvector<scalar_type> >
        B(mf.nb_basic_dof(), kept_dofs.card());
      gmm::mult(mf.extension_matrix(), gmm::transposed(RR), B);
      gmm::copy(B, E_);
    }
    else {
      gmm::copy(RR, R_); gmm::copy(gmm::transposed(RR), E_);
    }
    use_reduction = true;

    is_adapted = true; touch(); v_num = act_counter();
  }

  // invalid function for a mesh change.
  // dal::bit_vector partial_mesh_fem::retrieve_kept_dofs() const
  // {
  //   base_vector full(nb_basic_dof());
  //   for (size_type i = 0; i < full.size(); ++i) full[i] = i;
  //   base_vector reduced(nb_dof());
  // 
  //   if (R_.ncols() > 0) gmm::mult(R_, full, reduced);
  //   else reduced = full;
  // 
  //   dal::bit_vector kept_dofs;
  //   for (size_type i=0; i < reduced.size(); ++i) kept_dofs.add(reduced[i]); 
  // 
  //   return kept_dofs;
  // }

  void partial_mesh_fem::write_to_file(std::ostream &ost) const
  { context_check(); mf.context_check();
    gmm::stream_standard_locale sl(ost);
    ost << '\n' << "BEGIN MESH_FEM" << '\n' << '\n';
    mf.write_basic_to_file(ost);
    write_reduction_matrices_to_file(ost);
    ost << "END MESH_FEM" << '\n';
  }

  void partial_mesh_fem::write_to_file(const std::string &name,
                                       bool with_mesh) const {
    std::ofstream o(name.c_str());
    GMM_ASSERT1(o, "impossible to open file '" << name << "'");
    o << "% GETFEM MESH_FEM FILE " << '\n';
    o << "% GETFEM VERSION " << GETFEM_VERSION << '\n' << '\n' << '\n';
    if (with_mesh) mf.linked_mesh().write_to_file(o);
    write_to_file(o);
  }

  dal::bit_vector select_dofs_from_im(const mesh_fem &mf, const mesh_im &mim,
                                      unsigned P) {
    const mesh &m = mf.linked_mesh();
    unsigned N = m.dim();
    if (P == unsigned(-1)) P = N;
    base_matrix G;
    bgeot::pgeometric_trans pgt_old = 0;
    bgeot::pgeotrans_precomp pgp2 = 0;
    getfem::pfem pf_old = 0;
    getfem::pfem_precomp pfp = 0;
    pintegration_method pim1 = 0;

    std::vector<scalar_type> areas(mf.nb_basic_dof());
    std::vector<scalar_type> area_supports(mf.nb_basic_dof());
    dal::bit_vector kept_dofs;

    for (dal::bv_visitor cv(mim.convex_index()); !cv.finished(); ++cv) {
      bgeot::vectors_to_base_matrix(G, m.points_of_convex(cv));
      bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
      pintegration_method pim = mim.int_method_of_element(cv);
      if (pim == im_none()) continue;
      getfem::pfem pf = mf.fem_of_element(cv);
      GMM_ASSERT1(pim->type() == IM_APPROX,
                  "Works only with approximate integration");
      papprox_integration pai2= pim->approx_method();
      static papprox_integration pai2_old = 0;
      if (pgt_old != pgt || pai2 != pai2_old) {
        pim1 = getfem::classical_approx_im(pgt, 2);
              pgp2 = bgeot::geotrans_precomp(pgt,&(pai2->integration_points()),pim);
      }
      if (pai2 != pai2_old || pf != pf_old) {
        pf_old = pf;
        pfp = getfem::fem_precomp(pf, &(pai2->integration_points()), pim);
      }
      pai2_old = pai2;
      pgt_old = pgt;

      bgeot::geotrans_interpolation_context c2(pgp2, 0, G);
      scalar_type area1 = convex_area_estimate(pgt, G, pim1);

      size_type tdim = mf.get_qdim() / pf->target_dim();

      for (size_type i = 0; i < pai2->nb_points_on_convex(); ++i) {
        for (unsigned d = 0; d < pf->nb_dof(cv); ++d) {
          for (size_type j = 0; j < tdim; ++j) {
            if (i == 0)
              areas[mf.ind_basic_dof_of_element(cv)[d*tdim+j]] += area1;
            c2.set_ii(i);
            area_supports[mf.ind_basic_dof_of_element(cv)[d*tdim+j]]
              += pai2->coeff(i) * c2.J() * gmm::sqr(pfp->val(i)[d]);
          }
          //            * ((gmm::abs(pfp->val(i)[d]) < 1e-10) ? 0.0 : 1.0);
        }
      }
    }


    std::vector<scalar_type> areas2(mf.nb_dof());
    std::vector<scalar_type> area_supports2(mf.nb_dof());

    if (mf.is_reduced()) {
      gmm::mult(gmm::transposed(mf.extension_matrix()), areas, areas2);
      gmm::mult(gmm::transposed(mf.extension_matrix()), area_supports,
                area_supports2);
    }
    else {
      gmm::copy(areas, areas2);
      gmm::copy(area_supports, area_supports2);
    }

    for (size_type i = 0; i < mf.nb_dof(); ++i) {
      if (area_supports2[i] > pow(1e-14 * areas2[i], scalar_type(P) / N))
        kept_dofs.add(i);
    }

    return kept_dofs;
  }



}  /* end of namespace getfem.                                            */

