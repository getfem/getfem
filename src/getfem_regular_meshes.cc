/*===========================================================================

 Copyright (C) 1999-2012 Yves Renard

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


#include "getfem/getfem_regular_meshes.h"

namespace getfem
{
  void parallelepiped_regular_simplex_mesh_
  (mesh &me, dim_type N, const base_node &org,
   const base_small_vector *ivect, const size_type *iref) {
    bgeot::convex<base_node>
      pararef = *(bgeot::parallelepiped_of_reference(N));

    if (N >= 5) GMM_WARNING1("CAUTION : Simplexification in dimension >= 5 "
                             "has not been tested and the resulting mesh "
                             "should be not conformal");

    const bgeot::mesh_structure &sl
      = *(bgeot::parallelepiped_of_reference(N)->simplexified_convex());

    base_node a = org;
    size_type i, nbpt = pararef.nb_points();

    for (i = 0; i < nbpt; ++i) {
      gmm::clear(a);
      for (dim_type n = 0; n < N; ++n)
        gmm::add(gmm::scaled(ivect[n],pararef.points()[i][n]),a);
      pararef.points()[i] = a;
    }

    // bgeot::simplexify(cvt, sl, pararef.points(), N, me.eps());

    size_type nbs = sl.nb_convex();
    std::vector<size_type> tab1(N+1), tab(N), tab3(nbpt);
    size_type total = 0;
    std::fill(tab.begin(), tab.end(), 0);
    while (tab[N-1] != iref[N-1]) {
      for (a = org, i = 0; i < N; i++)
        gmm::add(gmm::scaled(ivect[i],scalar_type(tab[i])),a);
        //a.addmul(scalar_type(tab[i]), ivect[i]);

      for (i = 0; i < nbpt; i++)
        tab3[i] = me.add_point(a + pararef.points()[i]);

      for (i = 0; i < nbs; i++) {
        const mesh::ind_cv_ct &tab2 = sl.ind_points_of_convex(i);
        for (dim_type l = 0; l <= N; l++)
          // tab1[l] = tab3[tab2[l]];
          tab1[l] = tab3[(tab2[l]
                          + (((total & 1) && N != 3) ? (nbpt/2) : 0)) % nbpt];
        me.add_simplex(N, tab1.begin());
      }

      for (dim_type l = 0; l < N; l++) {
        tab[l]++; total++;
        if (l < N-1 && tab[l] >= iref[l]) { total -= tab[l]; tab[l] = 0; }
        else break;
      }
    }
  }


  void parallelepiped_regular_prism_mesh_
  (mesh &me, dim_type N, const base_node &org,
   const base_small_vector *ivect, const size_type *iref) {
    mesh aux;
    parallelepiped_regular_simplex_mesh_(aux, dim_type(N-1), org, ivect, iref);
    std::vector<base_node> ptab(2 * N);

    for (dal::bv_visitor cv(aux.convex_index()); !cv.finished(); ++cv) {
      std::copy(aux.points_of_convex(cv).begin(),
                aux.points_of_convex(cv).end(), ptab.begin());

      for (size_type k = 0; k < iref[N-1]; ++k) {

        for (dim_type j = 0; j < N; ++j) ptab[j+N] = ptab[j] + ivect[N-1];
        me.add_prism_by_points(N, ptab.begin());

        std::copy(ptab.begin()+N, ptab.end(), ptab.begin());
      }
    }
  }



  void parallelepiped_regular_mesh_
  (mesh &me, dim_type N, const base_node &org,
   const base_small_vector *ivect, const size_type *iref, bool linear_gt) {
    bgeot::convex<base_node>
      pararef = *(bgeot::parallelepiped_of_reference(N));
    base_node a = org;
    size_type i, nbpt = pararef.nb_points();

    for (i = 0; i < nbpt; ++i) {
      gmm::clear(a);
      for (dim_type n = 0; n < N; ++n)
        gmm::add(gmm::scaled(ivect[n],pararef.points()[i][n]),a);
      //a.addmul(pararef.points()[i][n], ivect[n]);
      pararef.points()[i] = a;
    }

    std::vector<size_type> tab1(N+1), tab(N), tab3(nbpt);
    size_type total = 0;
    std::fill(tab.begin(), tab.end(), 0);
    while (tab[N-1] != iref[N-1]) {
      for (a = org, i = 0; i < N; i++)
        gmm::add(gmm::scaled(ivect[i], scalar_type(tab[i])),a);
      //a.addmul(scalar_type(tab[i]), ivect[i]);

      for (i = 0; i < nbpt; i++)
        tab3[i] = me.add_point(a + pararef.points()[i]);
      me.add_convex(linear_gt ?
                    bgeot::parallelepiped_linear_geotrans(N) :
                    bgeot::parallelepiped_geotrans(N, 1), tab3.begin());

      for (dim_type l = 0; l < N; l++) {
        tab[l]++; total++;
        if (l < N-1 && tab[l] >= iref[l]) { total -= tab[l]; tab[l] = 0; }
        else break;
      }
    }
  }

  /* deformation inside a unit square  -- ugly */
  static base_node shake_func(const base_node& x) {
    base_node z(x.size());
    scalar_type c1 = 1., c2 = 1.;
    for (size_type i=0; i < x.size(); ++i) {
      c1*=(x[i]*(1.-x[i]));
      c2*=(.5 - gmm::abs(x[i]-.5));
    }
    z[0] = x[0] + c1;  // 1.
    for (size_type i=1; i < x.size(); ++i) {
      z[i] = x[i] + c2/3.; // /10
    }
    return z;
  }

  static base_node radial_deformation(const base_node& x) {
    GMM_ASSERT1(x.size() == 2, "For two-dimensional meshes only. \n");
    base_node z(x.size());
    z[0] = x[0] - 0.5 ;
    z[1] = x[1] - 0.5 ;
    scalar_type r = sqrt( z[0] * z[0] + z[1] * z[1] ) ;
    scalar_type theta = atan2(z[1], z[0]);
    if ( r < 0.5 - 1.e-6)
//        theta += 1000. * gmm::sqrt(r) * (0.5 - r) * (0.5 - r) * (0.5 - r) * (0.5 - r) * gmm::sqrt(gmm::abs(0.1 - r)) * gmm::sqrt(gmm::abs(0.15 - r))  ;
    theta += 10000. * gmm::sqrt(r) * (0.5 - r) * (0.5 - r) * (0.5 - r) * (0.5 - r) * (0.1 - r) * (0.15 - r)  ;
    z[0] = r * cos(theta) + 0.5;
    z[1] = r * sin(theta) + 0.5;
    return z;
  }

  static void noise_unit_mesh(mesh& m, std::vector<size_type> nsubdiv,
                              bgeot::pgeometric_trans pgt) {
    size_type N = nsubdiv.size();
    for (dal::bv_visitor ip(m.points().index()); !ip.finished(); ++ip) {
      bool is_border = false;
      base_node& P = m.points()[ip];
      for (size_type i=0; i < N; ++i) {
        if (gmm::abs(P[i]) < 1e-10 || gmm::abs(P[i]-1.) < 1e-10)
          is_border = true;
      }
      if (!is_border) {
        P = shake_func(P);
        if (N == 2) P = radial_deformation(P) ;
        for (size_type i=0; i < N; ++i)
          P[i] += 0.*(double(1)/double(nsubdiv[i]* pgt->complexity()))
            * gmm::random(double());
      }
    }
  }


  void regular_unit_mesh(mesh& m, std::vector<size_type> nsubdiv,
                         bgeot::pgeometric_trans pgt, bool noised) {
    mesh msh;
    dim_type N = dim_type(nsubdiv.size());
    base_node org(N); gmm::clear(org);
    std::vector<base_small_vector> vtab(N);
    for (dim_type i = 0; i < N; i++) {
      vtab[i] = base_small_vector(N); gmm::clear(vtab[i]);
      (vtab[i])[i] = 1. / scalar_type(nsubdiv[i]) * 1.;
    }
    if (pgt->basic_structure() == bgeot::simplex_structure(N)) {
      getfem::parallelepiped_regular_simplex_mesh
        (msh, N, org, vtab.begin(), nsubdiv.begin());
    } else if (pgt->basic_structure() == bgeot::parallelepiped_structure(N)) {
      getfem::parallelepiped_regular_mesh
        (msh, N, org, vtab.begin(), nsubdiv.begin());
    } else if (pgt->basic_structure() == bgeot::prism_structure(N)) {
      getfem::parallelepiped_regular_prism_mesh
        (msh, N, org, vtab.begin(), nsubdiv.begin());
    } else {
      GMM_ASSERT1(false, "cannot build a regular mesh for "
                  << bgeot::name_of_geometric_trans(pgt));
    }

    m.clear();
    /* build a mesh with a geotrans of degree K */
    for (dal::bv_visitor cv(msh.convex_index()); !cv.finished(); ++cv) {
      if (pgt == msh.trans_of_convex(cv)) {
        m.add_convex_by_points(msh.trans_of_convex(cv),
                               msh.points_of_convex(cv).begin());
      } else {
        std::vector<base_node> pts(pgt->nb_points());
        for (size_type i=0; i < pgt->nb_points(); ++i) {
          pts[i] = msh.trans_of_convex(cv)->transform
            (pgt->convex_ref()->points()[i], msh.points_of_convex(cv));
        }
        m.add_convex_by_points(pgt, pts.begin());
      }
    }

    /* apply a continuous deformation + some noise */
    if (noised) noise_unit_mesh(m, nsubdiv, pgt);

    m.optimize_structure();
  }



  void regular_mesh(mesh& m, const std::string &st) {
    std::stringstream s(st);
    bgeot::md_param PARAM;
    PARAM.read_param_file(s);

    std::string GT = PARAM.string_value("GT");
    GMM_ASSERT1(!GT.empty(), "regular mesh : you have at least to "
                "specify the geometric transformation");
    bgeot::pgeometric_trans pgt =
      bgeot::geometric_trans_descriptor(GT);

    size_type N = pgt->dim();
    base_small_vector org(N); gmm::clear(org);

    const std::vector<bgeot::md_param::param_value> &o
      = PARAM.array_value("ORG");
    if (o.size() > 0) {
      GMM_ASSERT1(o.size() == N, "ORG parameter should be an array of size "
                  << N);
      for (size_type i = 0; i < N; ++i) {
        GMM_ASSERT1(o[i].type_of_param() == bgeot::md_param::REAL_VALUE,
                    "ORG should be a real array.");
        org[i] = o[i].real();
      }
    }

    bool noised = (PARAM.int_value("NOISED") != 0);

    std::vector<size_type> nsubdiv(N);
    gmm::fill(nsubdiv, 2);
    const std::vector<bgeot::md_param::param_value> &ns
      = PARAM.array_value("NSUBDIV");
    if (ns.size() > 0) {
      GMM_ASSERT1(ns.size() == N,
                  "NSUBDIV parameter should be an array of size " << N);
      for (size_type i = 0; i < N; ++i) {
        GMM_ASSERT1(ns[i].type_of_param() == bgeot::md_param::REAL_VALUE,
                    "NSUBDIV should be an integer array.");
        nsubdiv[i] = size_type(ns[i].real()+0.5);
      }
    }

    base_small_vector sizes(N);
    gmm::fill(sizes, 1.0);

    const std::vector<bgeot::md_param::param_value> &si
      = PARAM.array_value("SIZES");
    if (si.size() > 0) {
      GMM_ASSERT1(si.size() == N,
                  "SIZES parameter should be an array of size " << N);
      for (size_type i = 0; i < N; ++i) {
        GMM_ASSERT1(si[i].type_of_param() == bgeot::md_param::REAL_VALUE,
                    "SIZES should be a real array.");
        sizes[i] = si[i].real();
      }
    }

    regular_unit_mesh(m, nsubdiv, pgt, noised);

    base_matrix M(N,N);
    for (size_type i=0; i < N; ++i) M(i,i) = sizes[i];
    m.transformation(M);
    m.translation(org);

  }


  void regular_ball_mesh(mesh& m, const std::string &st) {
    std::stringstream s(st);
    bgeot::md_param PARAM;
    PARAM.read_param_file(s);

    std::string GT = PARAM.string_value("GT");
    GMM_ASSERT1(!GT.empty(), "regular ball mesh : you have at least to "
                "specify the geometric transformation");
    bgeot::pgeometric_trans pgt =
      bgeot::geometric_trans_descriptor(GT);

    size_type N = pgt->dim();
    base_small_vector org(N);

    const std::vector<bgeot::md_param::param_value> &o
      = PARAM.array_value("ORG");
    if (o.size() > 0) {
      GMM_ASSERT1(o.size() == N, "ORG parameter should be an array of size "
                  << N);
      for (size_type i = 0; i < N; ++i) {
        GMM_ASSERT1(o[i].type_of_param() == bgeot::md_param::REAL_VALUE,
                    "ORG should be a real array.");
        org[i] = o[i].real();
      }
    }

    // "NOISED" applies only to the interior of all merged sub-meshes
    bool noised = (PARAM.int_value("NOISED") != 0);

    size_type nsubdiv0(2), nsubdiv1(2);
    const std::vector<bgeot::md_param::param_value> &ns
      = PARAM.array_value("NSUBDIV");
    if (ns.size() > 0) {
      GMM_ASSERT1(ns.size() == 2,
                  "NSUBDIV parameter should be an array of size " << 2);
      for (size_type i = 0; i < 2; ++i)
        GMM_ASSERT1(ns[i].type_of_param() == bgeot::md_param::REAL_VALUE,
                    "NSUBDIV should be an integer array.");
      nsubdiv0 = size_type(ns[0].real()+0.5);
      nsubdiv1 = size_type(ns[1].real()+0.5);
    }

    scalar_type radius(1), core_ratio(M_SQRT1_2);
    const std::vector<bgeot::md_param::param_value> &si
      = PARAM.array_value("SIZES");
    if (si.size() > 0) {
      GMM_ASSERT1(si.size() == 1,
                  "SIZES parameter should be an array of size " << 1);
      GMM_ASSERT1(si[0].type_of_param() == bgeot::md_param::REAL_VALUE,
                  "SIZES should be a real array.");
      radius = si[0].real();
    }

    
    std::vector<size_type> nsubdiv(N);
    gmm::fill(nsubdiv, nsubdiv0);
    regular_unit_mesh(m, nsubdiv, pgt, noised);
    std::vector<mesh> mm(N);
    for (size_type i = 0; i < N; ++i) {
      gmm::fill(nsubdiv, nsubdiv0);
      nsubdiv[i] = nsubdiv1;
      regular_unit_mesh(mm[i], nsubdiv, pgt, noised);
    }

    base_matrix M(N,N);
    for (size_type i=0; i < N; ++i) M(i,i) = core_ratio;
    m.transformation(M);

    scalar_type peel_ratio = scalar_type(1)-core_ratio; 
    for (size_type i = 0; i < N; ++i) {
      base_matrix MM(M);
      MM(i,i) = peel_ratio;
      mm[i].transformation(MM);

      std::vector<base_node> pts(mm[i].points().card(), base_node(N));
      size_type j(0);
      for (dal::bv_visitor pt(mm[i].points().index()); !pt.finished(); ++pt, ++j) {
        pts[j] = mm[i].points()[pt];
        for (unsigned k=0; k < N; ++k)
          if (k != i) pts[j][k] += (pts[j][k]/core_ratio) * pts[j][i];
      }
      j = size_type(0);
      for (dal::bv_visitor pt(mm[i].points().index()); !pt.finished(); ++pt, ++j)
          mm[i].points()[pt] = pts[j];

      base_small_vector trsl(N);
      trsl[i] = core_ratio;
      mm[i].translation(trsl);
      for (dal::bv_visitor cv(mm[i].convex_index()); !cv.finished(); ++cv)
	     m.add_convex_by_points(mm[i].trans_of_convex(cv), mm[i].points_of_convex(cv).begin());
    }

    std::vector<base_node> pts(m.points().card(), base_node(N));
    size_type j(0);
    for (dal::bv_visitor pt(m.points().index()); !pt.finished(); ++pt, ++j) {
      pts[j] = m.points()[pt];
      scalar_type maxcoord(0);
      size_type kmax(0);
      for (size_type k=0; k < N; ++k)
        if (gmm::abs(pts[j][k]) > maxcoord) {
          maxcoord = gmm::abs(pts[j][k]);
          kmax = k;
        }
      if (maxcoord > 1e-10) {
        scalar_type l(0), l0(0);
        for (size_type k=0; k < N; ++k) {
          if (k != kmax) {
            scalar_type theta = M_PI_4 * pts[j][k] / maxcoord;
            scalar_type c0 = std::min(scalar_type(1), maxcoord);
            pts[j][k] = c0*tan(theta)* maxcoord + (scalar_type(1)-c0)*pts[j][k];
          }
          l += pts[j][k] * pts[j][k];
        }
        l = sqrt(l);
        l0 = l/maxcoord;
        scalar_type scale(radius);
        scalar_type c(core_ratio);
        c *= std::max(scalar_type(0.3),
                      (scalar_type(1) - sqrt(l0*l0 - scalar_type(1))));
        if (l > c) {
          scale -= (l - c) / (l0 - c) * (radius - radius/(l/maxcoord));
        }
        for (size_type k=0; k < N; ++k)
          pts[j][k] *= scale;
      }
    }
    j = size_type(0);
    for (dal::bv_visitor pt(m.points().index()); !pt.finished(); ++pt, ++j)
      m.points()[pt] = pts[j];

    size_type symmetries(PARAM.int_value("SYMMETRIES"));
    symmetries = std::min(symmetries,N);

    for (size_type sym=0; sym < N-symmetries; ++sym) {
      size_type sym0 = (sym+1) % N;
      if (sym0 != 0) {
        gmm::clear(M);
        M(sym,sym0) = scalar_type(-1);
        M(sym0,sym) = scalar_type(1);
        for (size_type i=0; i < N; ++i)
          if (i != sym && i != sym0) M(i,i) = scalar_type(1);
      } else {
        base_matrix M1(M), M2(M);
        gmm::mult(M1,M2,M);
      }
      mesh m0;
      m0.copy_from(m);
      m0.transformation(M);
      for (dal::bv_visitor cv(m0.convex_index()); !cv.finished(); ++cv)
	     m.add_convex_by_points(m0.trans_of_convex(cv), m0.points_of_convex(cv).begin());
    }

    m.translation(org);

  }

}  /* end of namespace getfem.                                             */
