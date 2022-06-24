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

#include "getfem/bgeot_torus.h"
#include "getfem/dal_singleton.h"
#include "gmm/gmm_condition_number.h"
#include "getfem/getfem_mesh.h"
#include "getfem/getfem_integration.h"

#if GETFEM_HAVE_METIS_OLD_API
extern "C" void METIS_PartGraphKway(int *, int *, int *, int *, int *, int *,
                                    int *, int *, int *, int *, int *);
#elif GETFEM_HAVE_METIS
#  include <metis.h>
#endif

namespace getfem {

  gmm::uint64_type act_counter(void) {
    static gmm::uint64_type c = gmm::uint64_type(1);
    return ++c;
  }

  void mesh::sup_convex_from_regions(size_type c) {
    for (dal::bv_visitor i(valid_cvf_sets); !i.finished(); ++i)
      cvf_sets[i].sup_all(c);
    touch();
  }

  void mesh::swap_convex_in_regions(size_type c1, size_type c2) {
    for (dal::bv_visitor i(valid_cvf_sets); !i.finished(); ++i)
      cvf_sets[i].swap_convex(c1, c2);
    touch();
  }

  void mesh::handle_region_refinement(size_type ic,
                                      const std::vector<size_type> &icv,
                                      bool refine) {

    bgeot::pgeometric_trans pgt = trans_of_convex(ic);
    ind_set s;

    for (dal::bv_visitor ir(valid_cvf_sets); !ir.finished(); ++ir) {
      mesh_region &r = cvf_sets[ir];


      if (refine && r[ic].any()) {
        if (r[ic][0])
          for (size_type jc = 0; jc < icv.size(); ++jc)
            r.add(icv[jc]);

        bgeot::geotrans_inv_convex giv;
        for (short_type f = 0; f < pgt->structure()->nb_faces(); ++f) {
          if (r[ic][f+1]) {

            for (size_type jc = 0; jc < icv.size(); ++jc) {
              bgeot::pgeometric_trans pgtsub = trans_of_convex(icv[jc]);
              for (short_type fsub = 0; fsub < pgtsub->structure()->nb_faces();
                   ++fsub) {
                neighbors_of_convex(icv[jc], fsub, s);
                ind_set::const_iterator it = s.begin(), ite = s.end();
                bool found = false;
                for (; it != ite; ++it)
                  if (std::find(icv.begin(), icv.end(), *it) != icv.end())
                    { found = true; break; }
                if (found) continue;

                base_node pt, barycentre
                  =gmm::mean_value(pgtsub->convex_ref()->points_of_face(fsub));
                pt = pgtsub->transform(barycentre, points_of_convex(icv[jc]));

                giv.init(points_of_convex(ic), pgt);
                giv.invert(pt, barycentre);


                if (gmm::abs(pgt->convex_ref()->is_in_face(f,barycentre)) < 0.001)
                  r.add(icv[jc], fsub);
              }
            }
          }
        }
      }

      for (size_type jc = 0; jc < icv.size(); ++jc)
        if (!refine && r[icv[jc]].any()) {
          if (r[icv[jc]][0])
            r.add(ic);
          bgeot::geotrans_inv_convex giv;
          bgeot::pgeometric_trans pgtsub = trans_of_convex(icv[jc]);
          for (short_type fsub = 0; fsub < pgtsub->structure()->nb_faces();
               ++fsub)
            if (r[icv[jc]][fsub+1]) {
              base_node pt, barycentre
                = gmm::mean_value(pgtsub->convex_ref()->points_of_face(fsub));
              pt = pgtsub->transform(barycentre, points_of_convex(icv[jc]));

              giv.init(points_of_convex(ic), pgt);
              giv.invert(pt, barycentre);

              for (short_type f = 0; f < pgt->structure()->nb_faces(); ++f)
                if (gmm::abs(pgt->convex_ref()->is_in_face(f,barycentre)) < 0.001)
                  { r.add(ic, f); break; }
          }
        }
    }
  }

  void mesh::init(void) {
#if GETFEM_PARA_LEVEL > 1
    modified = true;
#endif
    cuthill_mckee_uptodate = false;
  }

  mesh::mesh(const std::string name) : name_(name)  { init(); }

  mesh::mesh(const bgeot::basic_mesh &m, const std::string name)
    : bgeot::basic_mesh(m), name_(name)  { init(); }

  mesh::mesh(const mesh &m) : context_dependencies(),
                              std::enable_shared_from_this<getfem::mesh>()
  { copy_from(m); }

  mesh &mesh::operator=(const mesh &m) {
    clear_dependencies();
    clear();
    copy_from(m);
    return *this;
  }

#if GETFEM_PARA_LEVEL > 1

  void mesh::compute_mpi_region() const {
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size < 2) {
      mpi_region = mesh_region::all_convexes();
      mpi_region.from_mesh(*this);
    } else {
      int ne = int(nb_convex());
      std::vector<int> xadj(ne+1), adjncy, numelt(ne), npart(ne);
      std::vector<int> indelt(nb_allocated_convex());

      double t_ref = MPI_Wtime();

      int j = 0, k = 0;
      ind_set s;
      for (dal::bv_visitor ic(convex_index()); !ic.finished(); ++ic, ++j) {
        numelt[j] = ic;
        indelt[ic] = j;
      }
      j = 0;
      for (dal::bv_visitor ic(convex_index()); !ic.finished(); ++ic, ++j) {
        xadj[j] = k;
        neighbors_of_convex(ic, s);
        for (ind_set::iterator it = s.begin();
             it != s.end(); ++it) { adjncy.push_back(indelt[*it]); ++k; }
      }
      xadj[j] = k;

#ifdef GETFEM_HAVE_METIS_OLD_API
      int wgtflag = 0, numflag = 0, edgecut;
      int options[5] = {0,0,0,0,0};
      METIS_PartGraphKway(&ne, &(xadj[0]), &(adjncy[0]), 0, 0, &wgtflag,
                          &numflag, &size, options, &edgecut, &(npart[0]));
#else
      int ncon = 1, edgecut;
      int options[METIS_NOPTIONS] = { 0 };
      METIS_SetDefaultOptions(options);
      METIS_PartGraphKway(&ne, &ncon, &(xadj[0]), &(adjncy[0]), 0, 0, 0,
                          &size, 0, 0, options, &edgecut, &(npart[0]));
#endif

      for (size_type i = 0; i < size_type(ne); ++i)
        if (npart[i] == rank) mpi_region.add(numelt[i]);

      if (MPI_IS_MASTER())
        cout << "Partition time "<< MPI_Wtime()-t_ref << endl;
    }
    modified = false;
    valid_sub_regions.clear();
  }

  void mesh::compute_mpi_sub_region(size_type n) const {
    if (valid_cvf_sets.is_in(n)) {
      mpi_sub_region[n] = mesh_region::intersection(cvf_sets[n], mpi_region);
    }
    else
      mpi_sub_region[n] = mesh_region();
    valid_sub_regions.add(n);
  }

  void mesh::intersect_with_mpi_region(mesh_region &rg) const {
    if (rg.id() == mesh_region::all_convexes().id()) {
      rg = get_mpi_region();
    } else if (int(rg.id()) >= 0) {
      rg = get_mpi_sub_region(rg.id());
    } else
      rg = mesh_region::intersection(rg, get_mpi_region());
  }
#endif

  void mesh::optimize_structure(bool with_renumbering) {
    pts.resort();
    size_type i, j = nb_convex(), nbc = j;
    for (i = 0; i < j; i++)
      if (!convex_tab.index_valid(i))
        swap_convex(i, convex_tab.ind_last());
    if (pts.size())
      for (i = 0, j = pts.size()-1;
           i < j && j != ST_NIL; ++i, --j) {
        while (i < j && j != ST_NIL && pts.index()[i]) ++i;
        while (i < j && j != ST_NIL && !(pts.index()[j])) --j;
        if (i < j && j != ST_NIL ) swap_points(i, j);
      }
    if (with_renumbering) { // Could be optimized no using only swap_convex
      std::vector<size_type> cmk, iord(nbc), iordinv(nbc);
      for (i = 0; i < nbc; ++i) iord[i] = iordinv[i] = i;

      bgeot::cuthill_mckee_on_convexes(*this, cmk);
      for (i = 0; i < nbc; ++i) {
        j = iordinv[cmk[i]];
        if (i != j) {
          swap_convex(i, j);
          std::swap(iord[i], iord[j]);
          std::swap(iordinv[iord[i]], iordinv[iord[j]]);
        }
      }
    }
  }

  void mesh::translation(const base_small_vector &V)
  { pts.translation(V); touch(); }

  void mesh::transformation(const base_matrix &M) {
    pts.transformation(M);
    Bank_info = std::unique_ptr<Bank_info_struct>();
    touch();
  }

  void mesh::bounding_box(base_node& Pmin, base_node& Pmax) const {
    bool is_first = true;
    Pmin.clear(); Pmax.clear();
    for (dal::bv_visitor i(pts.index()); !i.finished(); ++i) {
      if (is_first) { Pmin = Pmax = pts[i]; is_first = false; }
      else for (dim_type j=0; j < dim(); ++j) {
        Pmin[j] = std::min(Pmin[j], pts[i][j]);
        Pmax[j] = std::max(Pmax[j], pts[i][j]);
      }
    }
  }

  void mesh::clear(void) {
    mesh_structure::clear();
    pts.clear();
    gtab.clear(); trans_exists.clear();
    cvf_sets.clear(); valid_cvf_sets.clear();
    cvs_v_num.clear();
    Bank_info = nullptr;
    touch();
  }

  size_type mesh::add_segment(size_type a, size_type b) {
    size_type ipt[2]; ipt[0] = a; ipt[1] = b;
    return add_convex(bgeot::simplex_geotrans(1, 1), &(ipt[0]));
  }

  size_type mesh::add_triangle(size_type a,
                               size_type b, size_type c) {
    size_type ipt[3]; ipt[0] = a; ipt[1] = b; ipt[2] = c;
    return add_simplex(2, &(ipt[0]));
  }

  size_type mesh::add_triangle_by_points
    (const base_node &p1, const base_node &p2, const base_node &p3)
  { return add_triangle(add_point(p1), add_point(p2), add_point(p3)); }

  size_type mesh::add_tetrahedron(size_type a, size_type b,
                                  size_type c, size_type d) {
    size_type ipt[4]; ipt[0] = a; ipt[1] = b; ipt[2] = c; ipt[3] = d;
    return add_simplex(3, &(ipt[0]));
  }

  size_type mesh::add_pyramid(size_type a, size_type b,
                              size_type c, size_type d, size_type e) {
    size_type ipt[5] = {a, b, c, d, e};
    return add_convex(bgeot::pyramid_QK_geotrans(1), &(ipt[0]));
  }

  size_type mesh::add_tetrahedron_by_points
  (const base_node &p1, const base_node &p2,
   const base_node &p3, const base_node &p4) {
    return add_tetrahedron(add_point(p1), add_point(p2),
                           add_point(p3), add_point(p4));
  }

  void mesh::merge_convexes_from_mesh(const mesh &msource, size_type rg,
                                      scalar_type tol) {

    size_type nbpt = points_index().last()+1;
    GMM_ASSERT1(nbpt == nb_points(),
                "Please call the optimize_structure() function before "
                "merging elements from another mesh");
    GMM_ASSERT1(rg == size_type(-1) || msource.region(rg).is_only_convexes(),
                "The provided mesh region should only contain convexes");

    const dal::bit_vector &convexes = (rg == size_type(-1))
                                    ? msource.convex_index()
                                    : msource.region(rg).index();
    std::vector<size_type> old2new(msource.points_index().last()+1, size_type(-1));
    for (dal::bv_visitor cv(convexes); !cv.finished(); ++cv) {

      bgeot::pgeometric_trans pgt = msource.trans_of_convex(cv);
      short_type nb = short_type(pgt->nb_points());
      const ind_cv_ct &rct = msource.ind_points_of_convex(cv);
      GMM_ASSERT1(nb == rct.size(), "Internal error");
      std::vector<size_type> ind(nb);
      for (short_type i = 0; i < nb; ++i) {
        size_type old_pid = rct[i];
        size_type new_pid = old2new[old_pid];
        if (new_pid == size_type(-1)) {
          size_type next_pid = points_index().last()+1;
          base_node pt = msource.points()[old_pid];
          new_pid = pts.add_node(pt, tol);
          if (new_pid < next_pid && new_pid >= nbpt) {
            // do not allow internal merging of nodes in the source mesh
            new_pid = pts.add_node(pt, -1.);
            GMM_ASSERT1(new_pid == next_pid, "Internal error");
          }
          old2new[old_pid] = new_pid;
        }
        ind[i] = new_pid;
      }
      add_convex(pgt, ind.begin());
    }
  }

  void mesh::sup_convex(size_type ic, bool sup_points) {
    static std::vector<size_type> ipt;
    if (sup_points) {
      const ind_cv_ct &ct = ind_points_of_convex(ic);
      ipt.assign(ct.begin(), ct.end());
    }
    bgeot::mesh_structure::sup_convex(ic);
    if (sup_points)
      for (const size_type &ip : ipt) sup_point(ip);
    trans_exists.sup(ic);
    sup_convex_from_regions(ic);
    if (Bank_info.get()) Bank_sup_convex_from_green(ic);
    touch();
  }

  void mesh::swap_convex(size_type i, size_type j) {
    if (i != j) {
      bgeot::mesh_structure::swap_convex(i,j);
      trans_exists.swap(i, j);
      gtab.swap(i,j);
      swap_convex_in_regions(i, j);
      if (Bank_info.get()) Bank_swap_convex(i,j);
      cvs_v_num[i] = cvs_v_num[j] = act_counter(); touch();
    }
  }

  base_small_vector mesh::normal_of_face_of_convex(size_type ic, short_type f,
                                                   const base_node &pt) const {
    bgeot::pgeometric_trans pgt = trans_of_convex(ic);
    base_matrix G(dim(),pgt->nb_points());
    vectors_to_base_matrix(G, points_of_convex(ic));
    bgeot::geotrans_interpolation_context c(trans_of_convex(ic), pt, G);
    return bgeot::compute_normal(c, f);
  }


  base_small_vector mesh::normal_of_face_of_convex(size_type ic, short_type f,
                                                   size_type n) const {
    bgeot::pgeometric_trans pgt = trans_of_convex(ic);
    bgeot::pgeotrans_precomp pgp
      = bgeot::geotrans_precomp(pgt, pgt->pgeometric_nodes(), 0);
    base_matrix G;
    vectors_to_base_matrix(G, points_of_convex(ic));
    bgeot::geotrans_interpolation_context
      c(pgp,pgt->structure()->ind_points_of_face(f)[n], G);
    return bgeot::compute_normal(c, f);
  }

  base_small_vector mesh::mean_normal_of_face_of_convex(size_type ic,
                                                        short_type f) const {
    bgeot::pgeometric_trans pgt = trans_of_convex(ic);
    base_matrix G; vectors_to_base_matrix(G,points_of_convex(ic));
    base_small_vector ptmean(dim());
    size_type nbpt = pgt->structure()->nb_points_of_face(f);
    for (size_type i = 0; i < nbpt; ++i)
      gmm::add(pgt->geometric_nodes()[pgt->structure()->ind_points_of_face(f)[i]], ptmean);
    ptmean /= scalar_type(nbpt);
    bgeot::geotrans_interpolation_context c(pgt,ptmean, G);
    base_small_vector n = bgeot::compute_normal(c, f);
    n /= gmm::vect_norm2(n);
    return n;
  }

  base_matrix mesh::local_basis_of_face_of_convex(size_type ic, short_type f,
                                                  const base_node &pt) const {
    bgeot::pgeometric_trans pgt = trans_of_convex(ic);
    base_matrix G(dim(),pgt->nb_points());
    vectors_to_base_matrix(G,points_of_convex(ic));
    bgeot::geotrans_interpolation_context c(trans_of_convex(ic), pt, G);
    return bgeot::compute_local_basis(c, f);
  }

  base_matrix mesh::local_basis_of_face_of_convex(size_type ic, short_type f,
                                                         size_type n) const {
    bgeot::pgeometric_trans pgt = trans_of_convex(ic);
    bgeot::pgeotrans_precomp pgp
      = bgeot::geotrans_precomp(pgt, pgt->pgeometric_nodes(), 0);
    base_matrix G(dim(),pgt->nb_points());
    vectors_to_base_matrix(G,points_of_convex(ic));
    bgeot::geotrans_interpolation_context
      c(pgp,pgt->structure()->ind_points_of_face(f)[n], G);
    return bgeot::compute_local_basis(c, f);
  }

  scalar_type  mesh::convex_area_estimate(size_type ic, size_type deg) const {
    base_matrix G;
    bgeot::vectors_to_base_matrix(G, points_of_convex(ic));
    return getfem::convex_area_estimate
      (trans_of_convex(ic), G, classical_approx_im(trans_of_convex(ic),
                                                   dim_type(deg)));
  }

  scalar_type  mesh::convex_quality_estimate(size_type ic) const {
    base_matrix G;
    bgeot::vectors_to_base_matrix(G, points_of_convex(ic));
    auto pgt = trans_of_convex(ic);
    if (auto pgt_torus = dynamic_cast<const bgeot::torus_geom_trans*>(pgt.get())) {
      pgt = pgt_torus->get_original_transformation();
      G.resize(pgt->dim(), G.ncols());
    }
    return getfem::convex_quality_estimate(pgt, G);
  }

  scalar_type  mesh::convex_radius_estimate(size_type ic) const {
    base_matrix G;
    bgeot::vectors_to_base_matrix(G, points_of_convex(ic));
    return getfem::convex_radius_estimate(trans_of_convex(ic), G);
  }

  scalar_type mesh::minimal_convex_radius_estimate() const {
    if (convex_index().empty()) return 1;
    scalar_type r = convex_radius_estimate(convex_index().first_true());
    for (dal::bv_visitor cv(convex_index()); !cv.finished(); ++cv) {
      r = std::min(r, convex_radius_estimate(cv));
    }
    return r;
  }

  scalar_type mesh::maximal_convex_radius_estimate() const {
    if (convex_index().empty()) return 1;
    scalar_type r = convex_radius_estimate(convex_index().first_true());
    for (dal::bv_visitor cv(convex_index()); !cv.finished(); ++cv) {
      r = std::max(r, convex_radius_estimate(cv));
    }
    return r;
  }

  void mesh::set_name(const std::string& name){name_=name;}

  void mesh::copy_from(const mesh& m) {
    clear();
    set_name(m.name_);
    bgeot::basic_mesh::operator=(m);
    for (const auto &kv : m.cvf_sets) {
      if (kv.second.get_parent_mesh() != 0)
        cvf_sets[kv.first].set_parent_mesh(this);
      cvf_sets[kv.first] = kv.second;
    }
    valid_cvf_sets = m.valid_cvf_sets;
    cvs_v_num.clear();
    gmm::uint64_type d = act_counter();
    for (dal::bv_visitor i(convex_index()); !i.finished(); ++i)
      cvs_v_num[i] = d;
    Bank_info = std::unique_ptr<Bank_info_struct>();
    if (m.Bank_info.get())
      Bank_info = std::make_unique<Bank_info_struct>(*(m.Bank_info));
  }

  struct mesh_convex_structure_loc {
    bgeot::pgeometric_trans cstruct;
    std::vector<size_type> pts;
  };

  void mesh::read_from_file(std::istream &ist) {
    gmm::stream_standard_locale sl(ist);
    dal::bit_vector npt;
    dal::dynamic_array<double> tmpv;
    std::string tmp;
    bool te = false, please_get = true;

    ist.precision(16);
    clear();
    ist.seekg(0);ist.clear();
    bgeot::read_until(ist, "BEGIN POINTS LIST");

    while (!te) {
      if (please_get) bgeot::get_token(ist, tmp); else please_get = true;

      if (!bgeot::casecmp(tmp, "END"))
      { te = true; }
      else if (!bgeot::casecmp(tmp, "POINT")) {
        bgeot::get_token(ist, tmp);
        if (!bgeot::casecmp(tmp, "COUNT")) {
          bgeot::get_token(ist, tmp); // Ignored. Used in some applications
        } else {
          size_type ip = atoi(tmp.c_str());
          dim_type d = 0;
          GMM_ASSERT1(!npt.is_in(ip),
                      "Two points with the same index. loading aborted.");
          npt.add(ip);
          bgeot::get_token(ist, tmp);
          while (isdigit(tmp[0]) || tmp[0] == '-' || tmp[0] == '+'
                 || tmp[0] == '.')
            { tmpv[d++] = atof(tmp.c_str()); bgeot::get_token(ist, tmp); }
          please_get = false;
          base_node v(d);
          for (size_type i = 0; i < d; i++) v[i] = tmpv[i];
          size_type ipl = add_point(v);
          if (ip != ipl) {
            GMM_ASSERT1(!npt.is_in(ipl), "Two points [#" << ip << " and #"
                        << ipl << "] with the same coords "<< v
                        << ". loading aborted.");
            swap_points(ip, ipl);
          }
        }
      } else if (tmp.size()) {
        GMM_ASSERT1(false, "Syntax error in file, at token '" << tmp
                    << "', pos=" << std::streamoff(ist.tellg()));
      } else if (ist.eof()) {
        GMM_ASSERT1(false, "Unexpected end of stream while reading mesh");
      }
    }

    bool tend = false;
    dal::dynamic_array<mesh_convex_structure_loc> cv;
    dal::bit_vector ncv;

    ist.seekg(0);
    if (!bgeot::read_until(ist, "BEGIN MESH STRUCTURE DESCRIPTION"))
      GMM_ASSERT1(false, "This seems not to be a mesh file");

    while (!tend) {
      tend = !bgeot::get_token(ist, tmp);
      if (!bgeot::casecmp(tmp, "END"))
      { tend = true; }
      else if (!bgeot::casecmp(tmp, "CONVEX")) {
        size_type ic;
        bgeot::get_token(ist, tmp);
        if (!bgeot::casecmp(tmp, "COUNT")) {
          bgeot::get_token(ist, tmp); // Ignored. Used in some applications
        } else {
          ic = gmm::abs(atoi(tmp.c_str()));
          GMM_ASSERT1(!ncv.is_in(ic),
                      "Negative or repeated index, loading aborted.");
          ncv.add(ic);

          int rgt = bgeot::get_token(ist, tmp);
          if (rgt != 3) { // for backward compatibility with version 1.7
            char c; ist.get(c);
            while (!isspace(c)) { tmp.push_back(c); ist.get(c); }
          }

          bgeot::pgeometric_trans pgt = bgeot::geometric_trans_descriptor(tmp);
          size_type nb = pgt->nb_points();

          cv[ic].cstruct = pgt;
          cv[ic].pts.resize(nb);
          for (size_type i = 0; i < nb; i++) {
            bgeot::get_token(ist, tmp);
            cv[ic].pts[i] = gmm::abs(atoi(tmp.c_str()));
          }
        }
      }
      else if (tmp.size()) {
        GMM_ASSERT1(false, "Syntax error reading a mesh file "
                    " at pos " << std::streamoff(ist.tellg())
                    << "(expecting 'CONVEX' or 'END', found '"<< tmp << "')");
      } else if (ist.eof()) {
        GMM_ASSERT1(false, "Unexpected end of stream "
                    << "(missing BEGIN MESH/END MESH ?)");
      }
    }
    ist >> bgeot::skip("MESH STRUCTURE DESCRIPTION");

    for (dal::bv_visitor ic(ncv); !ic.finished(); ++ic) {
      size_type i = add_convex(cv[ic].cstruct, cv[ic].pts.begin());
      if (i != ic) swap_convex(i, ic);
    }

    tend = false;
    while (!tend) {
      tend = !bgeot::get_token(ist, tmp);
      // bool error = false;
      if (bgeot::casecmp(tmp, "BEGIN")==0) {
        bgeot::get_token(ist, tmp);
        if (bgeot::casecmp(tmp, "BOUNDARY")==0 ||
            bgeot::casecmp(tmp, "REGION")==0) {
          bgeot::get_token(ist, tmp);
          size_type bnum = atoi(tmp.c_str());
          bgeot::get_token(ist, tmp);
          while (true) {
            if (bgeot::casecmp(tmp, "END")!=0) {
              size_type ic = atoi(tmp.c_str());
              bgeot::get_token(ist, tmp);
              if (tmp[0] == '/') {
                bgeot::get_token(ist, tmp);
                if (!bgeot::casecmp(tmp, "END")) break;
                int f = atoi(tmp.c_str());
                region(bnum).add(ic, short_type(f));
                bgeot::get_token(ist, tmp);
              }
              else {
                region(bnum).add(ic);
                if (!bgeot::casecmp(tmp, "END")) break;
              }
            } else break;
          }
          bgeot::get_token(ist, tmp);
          bgeot::get_token(ist, tmp);
        } else tend = true;
        /*else GMM_ASSERT1(false, "Syntax error in file at token '"
          << tmp << "' [pos=" << std::streamoff(ist.tellg())
          << "]");*/
      } else tend=true;
    }
  }

  void mesh::read_from_file(const std::string &name) {
    std::ifstream o(name.c_str());
    GMM_ASSERT1(o, "Mesh file '" << name << "' does not exist");
    read_from_file(o);
    o.close();
  }

  template<class ITER>
    void write_tab_to_file_(std::ostream &ost, const ITER& b_, const ITER& e)
  { for (ITER b(b_) ; b != e; ++b) ost << "  " << *b; }

  template<class ITER>
    static void write_convex_to_file_(const mesh &ms,
                                      std::ostream &ost,
                                      ITER b, ITER e) {
    for ( ; b != e; ++b) {
      size_type i = b.index();
      ost << "  CONVEX " << i << "    \'"
          << bgeot::name_of_geometric_trans(ms.trans_of_convex(i)).c_str()
          << "\'    ";
      write_tab_to_file_(ost, ms.ind_points_of_convex(i).begin(),
                         ms.ind_points_of_convex(i).end()  );
      ost << '\n';
    }
  }

  template<class ITER> static void write_point_to_file_(std::ostream &ost,
                                                  ITER b, ITER e)
  { for ( ; b != e; ++b) ost << "  " << *b; ost << '\n'; }

  void mesh::write_to_file(std::ostream &ost) const {
    ost.precision(16);
    gmm::stream_standard_locale sl(ost);
    ost << '\n' << "BEGIN POINTS LIST" << '\n' << '\n';
    ost << "  POINT COUNT " << points().index().last_true()+1 << '\n';
    for (size_type i = 0; i < points_tab.size(); ++i)
      if (is_point_valid(i) ) {
        ost << "  POINT  " << i;
        write_point_to_file_(ost, pts[i].begin(), pts[i].end());
      }
    ost << '\n' << "END POINTS LIST" << '\n' << '\n' << '\n';

    ost << '\n' << "BEGIN MESH STRUCTURE DESCRIPTION" << '\n' << '\n';
    ost << "  CONVEX COUNT " << convex_index().last_true()+1 << '\n';
    write_convex_to_file_(*this, ost, convex_tab.tas_begin(),
                                      convex_tab.tas_end());
    ost << '\n' << "END MESH STRUCTURE DESCRIPTION" << '\n';

    for (dal::bv_visitor bnum(valid_cvf_sets); !bnum.finished(); ++bnum) {
      ost << "BEGIN REGION " << bnum << "\n" << region(bnum) << "\n"
          << "END REGION " << bnum << "\n";
    }
  }

  void mesh::write_to_file(const std::string &name) const {
    std::ofstream o(name.c_str());
    GMM_ASSERT1(o, "impossible to write to file '" << name << "'");
    o << "% GETFEM MESH FILE " << '\n';
    o << "% GETFEM VERSION " << GETFEM_VERSION << '\n' << '\n' << '\n';
    write_to_file(o);
    o.close();
  }

  size_type mesh::memsize(void) const {
    return bgeot::mesh_structure::memsize() - sizeof(bgeot::mesh_structure)
      + pts.memsize() + (pts.index().last_true()+1)*dim()*sizeof(scalar_type)
      + sizeof(mesh) + trans_exists.memsize() + gtab.memsize()
      + valid_cvf_sets.card()*sizeof(mesh_region) + valid_cvf_sets.memsize();
  }

  struct equilateral_to_GT_PK_grad_aux : public std::vector<base_matrix> {};
  static const base_matrix &equilateral_to_GT_PK_grad(dim_type N) {
    std::vector<base_matrix>
      &pbm = dal::singleton<equilateral_to_GT_PK_grad_aux >::instance();
    if (N > pbm.size()) pbm.resize(N);
    if (pbm[N-1].empty()) {
      bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(N,1);
      base_matrix Gr(N,N);
      base_matrix G(N,N+1);
      vectors_to_base_matrix
        (G, bgeot::equilateral_simplex_of_reference(N)->points());
      gmm::mult(G, bgeot::geotrans_precomp
                (pgt, pgt->pgeometric_nodes(), 0)->grad(0), Gr);
      gmm::lu_inverse(Gr);
      pbm[N-1].swap(Gr);
    }
    return pbm[N-1];
  }


  scalar_type convex_area_estimate(bgeot::pgeometric_trans pgt,
                                   const base_matrix& G,
                                   pintegration_method pi) {
    double area(0);
    static bgeot::pgeometric_trans pgt_old = 0;
    static bgeot::pgeotrans_precomp pgp = 0;
    static pintegration_method pim_old = 0;
    papprox_integration pai = get_approx_im_or_fail(pi);
    if (pgt_old != pgt || pim_old != pi) {
      pgt_old = pgt;
      pim_old = pi;
      pgp = bgeot::geotrans_precomp
        (pgt, pai->pintegration_points(), pi);
    }
    bgeot::geotrans_interpolation_context gic(pgp, 0, G);
    for (size_type i = 0; i < pai->nb_points_on_convex(); ++i) {
      gic.set_ii(i);
      area += pai->coeff(i) * gic.J();
    }
    return area;
  }

  /* TODO : use the geotrans from an "equilateral" reference element to
     the real element
     check if the sign of the determinants does change
     (=> very very bad quality of the convex)
  */
  scalar_type convex_quality_estimate(bgeot::pgeometric_trans pgt,
                                      const base_matrix& G) {
    static bgeot::pgeometric_trans pgt_old = 0;
    static bgeot::pgeotrans_precomp pgp = 0;
    if (pgt_old != pgt) {
      pgt_old=pgt;
      pgp=bgeot::geotrans_precomp(pgt, pgt->pgeometric_nodes(), 0);
    }

    size_type n = (pgt->is_linear()) ? 1 : pgt->nb_points();
    scalar_type q = 1;
    dim_type N = dim_type(G.nrows()), P = pgt->structure()->dim();
    base_matrix K(N,P);
    for (size_type ip=0; ip < n; ++ip) {
      gmm::mult(G, pgp->grad(ip), K);
      /* TODO : this is an ugly fix for simplexes only.. there should be
         a transformation of any pgt to the equivalent equilateral pgt
         (for prisms etc) */
      if (bgeot::basic_structure(pgt->structure())
          == bgeot::simplex_structure(P))
        gmm::mult(base_matrix(K),equilateral_to_GT_PK_grad(P),K);
      q = std::max(q, gmm::condition_number(K));
    }
    return 1./q;
  }

  scalar_type convex_radius_estimate(bgeot::pgeometric_trans pgt,
                                     const base_matrix& G) {
    static bgeot::pgeometric_trans pgt_old = 0;
    static bgeot::pgeotrans_precomp pgp = 0;
    if (pgt_old != pgt) {
      pgt_old=pgt;
      pgp=bgeot::geotrans_precomp(pgt, pgt->pgeometric_nodes(), 0);
    }
    size_type N = G.nrows();
    size_type n = (pgt->is_linear()) ? 1 : pgt->nb_points();
    scalar_type q = 0;
    base_matrix K(pgp->grad(0).ncols(),G.nrows());
    for (size_type ip=0; ip < n; ++ip) {
      gmm::mult(gmm::transposed(pgp->grad(ip)), gmm::transposed(G), K);
      scalar_type emax, emin; gmm::condition_number(K,emax,emin);
      q = std::max(q, emax);
    }
    return q * sqrt(scalar_type(N)) / scalar_type(2);
  }

  /* extract faces of convexes which are not shared
      + convexes whose dimension is smaller that m.dim()
  */
  void
  outer_faces_of_mesh(const getfem::mesh &m, const dal::bit_vector& cvlst,
                      convex_face_ct& flist) {
    for (dal::bv_visitor ic(cvlst); !ic.finished(); ++ic) {
      if (m.structure_of_convex(ic)->dim() == m.dim()) {
        for (short_type f = 0; f < m.structure_of_convex(ic)->nb_faces(); f++) {
          if (m.neighbor_of_convex(ic,f) == size_type(-1)) {
            flist.push_back(convex_face(ic,f));
          }
        }
      } else {
        flist.push_back(convex_face(ic, short_type(-1)));
      }
    }
  }

  void  outer_faces_of_mesh(const mesh &m,
                            const mesh_region &cvlst,
                            mesh_region &flist) {
    cvlst.error_if_not_convexes();
    for (mr_visitor i(cvlst); !i.finished(); ++i) {
      if (m.structure_of_convex(i.cv())->dim() == m.dim()) {
        for (short_type f = 0; f < m.structure_of_convex(i.cv())->nb_faces();
             f++) {
          size_type cv2 = m.neighbor_of_convex(i.cv(), f);
          if (cv2 == size_type(-1) || !cvlst.is_in(cv2)) {
            flist.add(i.cv(),f);
          }
        }
      }
      else flist.add(i.cv());
    }
  }

  /* Select all the faces of the given mesh region (counted twice if they
     are shared by two neighbor elements)
  */
  mesh_region all_faces_of_mesh(const mesh &m, const mesh_region &mr) {
    mesh_region mrr;
    mr.from_mesh(m);
    mr.error_if_not_convexes();

    for (mr_visitor i(mr); !i.finished(); ++i) {
      size_type cv1 = i.cv();
      short_type nbf = m.structure_of_convex(i.cv())->nb_faces();
      for (short_type f = 0; f < nbf; ++f)
        mrr.add(cv1, f);
    }
    return mrr;
  }
  
  /* Select all the faces sharing at least two element of the given mesh
      region. Each face is represented only once and is arbitrarily chosen
      between the two neighbor elements. Try to minimize the number of
      elements.
  */
  mesh_region inner_faces_of_mesh(const mesh &m, const mesh_region &mr) {
    mesh_region mrr;
    mr.from_mesh(m);
    mr.error_if_not_convexes();
    dal::bit_vector visited;
    bgeot::mesh_structure::ind_set neighbors;

    for (mr_visitor i(mr); !i.finished(); ++i) {
      size_type cv1 = i.cv();
      short_type nbf = m.structure_of_convex(i.cv())->nb_faces();
      bool neighbor_visited = false;
      for (short_type f = 0; f < nbf; ++f) {
        neighbors.resize(0); m.neighbors_of_convex(cv1, f, neighbors);
        for (size_type j = 0; j < neighbors.size(); ++j)
          if (visited.is_in(neighbors[j]))
            { neighbor_visited = true; break; }
      }
      if (!neighbor_visited) {
        for (short_type f = 0; f < nbf; ++f) {
          size_type cv2 = m.neighbor_of_convex(cv1, f);
          if (cv2 != size_type(-1) && mr.is_in(cv2)) mrr.add(cv1,f);
        }
        visited.add(cv1);
      }
    }

    for (mr_visitor i(mr); !i.finished(); ++i) {
      size_type cv1 = i.cv();
      if (!(visited.is_in(cv1))) {
        short_type nbf = m.structure_of_convex(i.cv())->nb_faces();
        for (short_type f = 0; f < nbf; ++f) {
          neighbors.resize(0); m.neighbors_of_convex(cv1, f, neighbors);
          bool ok = false;
          for (size_type j = 0; j < neighbors.size(); ++j)  {
            if (visited.is_in(neighbors[j])) { ok = false; break; }
            if (mr.is_in(neighbors[j])) { ok = true; }
          }
          if (ok) { mrr.add(cv1,f); }
        }
        visited.add(cv1);
      }
    }
    return mrr;
  }


  mesh_region select_faces_of_normal(const mesh &m, const mesh_region &mr,
                                     const base_small_vector &V,
                                     scalar_type angle) {
    mesh_region mrr;
    scalar_type threshold = gmm::vect_norm2(V)*cos(angle);
    for (getfem::mr_visitor i(mr); !i.finished(); ++i)
      if (i.f() != short_type(-1)) {
        base_node un = m.mean_normal_of_face_of_convex(i.cv(), i.f());
        if (gmm::vect_sp(V, un) - threshold >= -1E-8)
          mrr.add(i.cv(), i.f());
      }
    return mrr;
  }

  mesh_region select_faces_in_box(const mesh &m, const mesh_region &mr,
                                  const base_node &pt1, const base_node &pt2) {
    mesh_region mrr;
    size_type N = m.dim();
    GMM_ASSERT1(pt1.size() == N && pt2.size() == N, "Wrong dimensions");
    for (getfem::mr_visitor i(mr, m); !i.finished(); ++i)
      if (i.f() != short_type(-1)) {
        bgeot::mesh_structure::ind_pt_face_ct pt
          = m.ind_points_of_face_of_convex(i.cv(), i.f());

        bool is_in = true;
        for (bgeot::mesh_structure::ind_pt_face_ct::iterator it = pt.begin();
             it != pt.end(); ++it) {
          for (size_type j = 0; j < N; ++j)
            if (m.points()[*it][j] < pt1[j] || m.points()[*it][j] > pt2[j])
              { is_in = false; break; }
          if (!is_in) break;
        }
        if (is_in) mrr.add(i.cv(), i.f());
      }
    return mrr;
  }

  mesh_region select_faces_in_ball(const mesh &m, const mesh_region &mr,
                                   const base_node &center, scalar_type radius) {
    mesh_region mrr;
    size_type N = m.dim();
    GMM_ASSERT1(center.size() == N, "Wrong dimensions");
    for (getfem::mr_visitor i(mr, m); !i.finished(); ++i)
      if (i.f() != short_type(-1)) {
        bgeot::mesh_structure::ind_pt_face_ct pt
          = m.ind_points_of_face_of_convex(i.cv(), i.f());

        bool is_in = true;
        for (bgeot::mesh_structure::ind_pt_face_ct::iterator it = pt.begin();
             it != pt.end(); ++it) {
          scalar_type checked_radius = scalar_type(0.0);
          for (size_type j = 0; j < N; ++j)
            checked_radius += pow(m.points()[*it][j] - center[j], 2);
          checked_radius = std::sqrt(checked_radius);
          if (checked_radius > radius) { is_in = false; break; }
        }
        if (is_in) mrr.add(i.cv(), i.f());
      }
    return mrr;
  }

  mesh_region select_convexes_in_box(const mesh &m, const mesh_region &mr,
                                     const base_node &pt1, const base_node &pt2) {
    mesh_region mrr;
    size_type N = m.dim();
    GMM_ASSERT1(pt1.size() == N && pt2.size() == N, "Wrong dimensions");
    for (getfem::mr_visitor i(mr, m); !i.finished(); ++i)
      if (i.f() == short_type(-1)) {
        bgeot::mesh_structure::ind_cv_ct pt = m.ind_points_of_convex(i.cv());

        bool is_in = true;
        for (bgeot::mesh_structure::ind_cv_ct::iterator it = pt.begin();
             it != pt.end(); ++it) {
          for (size_type j = 0; j < N; ++j)
            if (m.points()[*it][j] < pt1[j] || m.points()[*it][j] > pt2[j])
              { is_in = false; break; }
          if (!is_in) break;
        }
        if (is_in) mrr.add(i.cv());
      }
    return mrr;
  }

  void extrude(const mesh& in, mesh& out, size_type nb_layers, short_type degree) {
    dim_type dim = in.dim();
    base_node pt(dim+1);
    out.clear();
    size_type nbpt = in.points().index().last()+1;
    GMM_ASSERT1(nbpt == in.points().index().card(),
                "please call the optimize_structure() method before using "
                "the mesh as a base for prismatic mesh");
    size_type nb_layers_total = nb_layers * degree;
    for (const base_node &inpt : in.points()) {
      std::copy(inpt.begin(), inpt.end(), pt.begin());
      pt[dim] = 0.0;
      for (size_type j = 0; j <= nb_layers_total;
           ++j, pt[dim] += scalar_type(1) / scalar_type(nb_layers_total))
        out.add_point(pt);
    }

    std::vector<size_type> tab;
    for (dal::bv_visitor cv(in.convex_index()); !cv.finished(); ++cv) {
      size_type nbp = in.nb_points_of_convex(cv);
      tab.resize((degree+1)*nbp);
      for (size_type j = 0; j < nb_layers; ++j) {
        for (short_type d = 0; d < degree+1; ++d)
          for (size_type k = 0; k < nbp; ++k)
            tab[k+nbp*d] = (nb_layers_total+1)*in.ind_points_of_convex(cv)[k] + j*degree + d;
        bgeot::pgeometric_trans pgt =
          bgeot::product_geotrans(in.trans_of_convex(cv),
                                  bgeot::simplex_geotrans(1,degree));
        out.add_convex(pgt, tab.begin());
      }
    }
  }
}


//
// Bank refinement
//


namespace bgeot {
  size_type refinement_simplexe_tab(size_type n, size_type **tab);
}

namespace getfem {

  bool mesh::edge::operator <(const edge &e) const {
    if (i0 < e.i0) return true;
    if (i0 > e.i0) return false;
    if (i1 < e.i1) return true;
    if (i1 > e.i1) return false;
    if (i2 < e.i2) return true;
    return false;
  }

  void mesh::Bank_sup_convex_from_green(size_type i) {
    if (Bank_info.get() && Bank_info->is_green_simplex.is_in(i)) {
      size_type igs = Bank_info->num_green_simplex[i];
      green_simplex &gs = Bank_info->green_simplices[igs];
      for (size_type j = 0; j < gs.sub_simplices.size(); ++j) {
        Bank_info->num_green_simplex.erase(gs.sub_simplices[j]);
        Bank_info->is_green_simplex.sup(gs.sub_simplices[j]);
      }
      Bank_info->green_simplices.sup(igs);
    }
  }

  void mesh::Bank_swap_convex(size_type i, size_type j) {
    if (Bank_info.get()) {
      Bank_info->is_green_simplex.swap(i, j);
      std::map<size_type, size_type>::iterator
        iti = Bank_info->num_green_simplex.find(i);
      std::map<size_type, size_type>::iterator
        ite = Bank_info->num_green_simplex.end();
      std::map<size_type, size_type>::iterator
        itj = Bank_info->num_green_simplex.find(j);
      size_type numi(0), numj(0);
      if (iti != ite)
        { numi = iti->second; Bank_info->num_green_simplex.erase(i); }
      if (itj != ite)
        { numj = itj->second; Bank_info->num_green_simplex.erase(j); }
      if (iti != ite) {
        Bank_info->num_green_simplex[j] = numi;
        green_simplex &gs = Bank_info->green_simplices[numi];
        for (size_type k = 0; k < gs.sub_simplices.size(); ++k)
          if (gs.sub_simplices[k] == i) gs.sub_simplices[k] = j;
          else if (gs.sub_simplices[k] == j) gs.sub_simplices[k] = i;
      }
      if (itj != ite) {
        Bank_info->num_green_simplex[i] = numj;
        if (iti == ite || numi != numj) {
          green_simplex &gs = Bank_info->green_simplices[numj];
          for (size_type k = 0; k < gs.sub_simplices.size(); ++k)
            if (gs.sub_simplices[k] == i) gs.sub_simplices[k] = j;
            else if (gs.sub_simplices[k] == j) gs.sub_simplices[k] = i;
        }
      }
    }
  }

  void mesh::Bank_build_first_mesh(mesh &m, size_type n) {
    bgeot::pconvex_ref pcr = bgeot::simplex_of_reference(dim_type(n), 2);
    m.clear();
    for (size_type ip = 0; ip < pcr->nb_points(); ++ip)
      m.add_point(pcr->points()[ip]);
    size_type *tab;
    size_type nbc = bgeot::refinement_simplexe_tab(n, &tab);
    for (size_type ic = 0; ic < nbc; ++ic, tab+=(n+1))
      m.add_simplex(dim_type(n), tab);
  }

  struct mesh_cache_for_Bank_basic_refine_convex : public mesh {};

  void mesh::Bank_basic_refine_convex(size_type i) {
    bgeot::pgeometric_trans pgt = trans_of_convex(i);
    size_type n = pgt->basic_structure()->dim();

    static bgeot::pgeometric_trans pgt1 = 0;

    mesh &mesh2 = dal::singleton<mesh_cache_for_Bank_basic_refine_convex>::instance();

    static bgeot::pstored_point_tab pspt = 0;
    static bgeot::pgeotrans_precomp pgp = 0;
    static std::vector<size_type> ipt, ipt2, icl;

    if (pgt != pgt1) {
      pgt1 = pgt;
      mesh mesh1;
      Bank_build_first_mesh(mesh1, n);

      mesh2.clear();
      ipt.resize(pgt->nb_points());
      for (size_type ic = 0; ic < mesh1.nb_convex(); ++ic) {
        assert(mesh1.convex_index().is_in(ic));
        bgeot::pgeometric_trans pgt2 = mesh1.trans_of_convex(ic);
        for (size_type ip = 0; ip < pgt->nb_points(); ++ip)
          ipt[ip] =
            mesh2.add_point(pgt2->transform(pgt->geometric_nodes()[ip],
                                            mesh1.points_of_convex(ic)));
        mesh2.add_convex(pgt, &ipt[0]);
      }

      // if (pspt) dal::del_stored_object(pspt);
      pspt = bgeot::store_point_tab(mesh2.points());
      pgp = bgeot::geotrans_precomp(pgt, pspt, 0);
    }

    base_node pt(dim());
    ipt.resize(pspt->size());
    for (size_type ip = 0; ip < pspt->size(); ++ip) {
      pgp->transform(points_of_convex(i), ip, pt);
      ipt[ip] = add_point(pt);
    }

    ipt2.resize(pgt->nb_points()); icl.resize(0);
    for (size_type ic = 0; ic < mesh2.nb_convex(); ++ic) {
      for (size_type j = 0; j < pgt->nb_points(); ++j)
        ipt2[j] = ipt[mesh2.ind_points_of_convex(ic)[j]];
      icl.push_back(add_convex(pgt, ipt2.begin()));
    }
    handle_region_refinement(i, icl, true);
    sup_convex(i, true);
  }

  void mesh::Bank_convex_with_edge(size_type i1, size_type i2,
                                   std::vector<size_type> &ipt) {
    ipt.resize(0);
    for (size_type k = 0; k < points_tab[i1].size(); ++k) {
      size_type cv = points_tab[i1][k], found = 0;
      const std::vector<size_type> &loc_ind = trans_of_convex(cv)->vertices();
      for (size_type i = 0; i < loc_ind.size(); ++i) {
        size_type ipp = convex_tab[cv].pts[loc_ind[i]];
        if (ipp == i1) ++found;
        if (ipp == i2) ++found;
      }
      GMM_ASSERT1(found <= 2, "Invalid convex with repeated nodes ");
      if (found == 2) ipt.push_back(cv);
    }
  }

  bool mesh::Bank_is_convex_having_points(size_type cv,
                                          const std::vector<size_type> &ipt) {
    size_type found = 0;
    const std::vector<size_type> &loc_ind = trans_of_convex(cv)->vertices();
    for (size_type i = 0; i < loc_ind.size(); ++i)
      if (std::find(ipt.begin(), ipt.end(),
                    convex_tab[cv].pts[loc_ind[i]]) != ipt.end()) ++found;
    return (found >= ipt.size());
  }

  void mesh::Bank_refine_normal_convex(size_type i) {
    bgeot::pgeometric_trans pgt = trans_of_convex(i);
    GMM_ASSERT1(pgt->basic_structure() == bgeot::simplex_structure(pgt->dim()),
                "Sorry, refinement is only working with simplices.");

    const std::vector<size_type> &loc_ind = pgt->vertices();
    for (size_type ip1 = 0; ip1 < loc_ind.size(); ++ip1)
      for (size_type ip2 = ip1+1; ip2 < loc_ind.size(); ++ip2)
        Bank_info->edges.insert(edge(ind_points_of_convex(i)[loc_ind[ip1]],
                                     ind_points_of_convex(i)[loc_ind[ip2]]));
    Bank_basic_refine_convex(i);
  }

  size_type mesh::Bank_test_and_refine_convex(size_type i,
                                         dal::bit_vector &b, bool ref) {
    if (Bank_info->is_green_simplex[i]) {
      size_type igs = Bank_info->num_green_simplex[i];
      green_simplex &gs = Bank_info->green_simplices[igs];

      size_type icc = add_convex_by_points(gs.pgt, gs.cv.points().begin());
      handle_region_refinement(icc, gs.sub_simplices, false);
      for (size_type ic = 0; ic < gs.sub_simplices.size(); ++ic) {
        sup_convex(gs.sub_simplices[ic], true);
        b.sup(gs.sub_simplices[ic]);
      }
      if (!ref)
        for (size_type ip = 0; ip < gs.ipt_loc.size(); ++ip)
          for (size_type jp = ip+1; jp < gs.ipt_loc.size(); ++jp)
            Bank_info->edges.insert
              (edge(ind_points_of_convex(icc)[gs.ipt_loc[ip]],
                    ind_points_of_convex(icc)[gs.ipt_loc[jp]]));
      Bank_sup_convex_from_green(i);
      if (ref) { Bank_refine_normal_convex(icc); return size_type(-1); }
      else return icc;
    }
    else if (ref) Bank_refine_normal_convex(i);

    return size_type(-1);
  }

  struct mesh_cache_for_Bank_build_green_simplexes : public mesh {};

  void mesh::Bank_build_green_simplexes(size_type ic,
                                        std::vector<size_type> &ipt) {
    GMM_ASSERT1(convex_index().is_in(ic), "Internal error");
    size_type igs = Bank_info->green_simplices.add(green_simplex());
    green_simplex &gs(Bank_info->green_simplices[igs]);
    std::vector<base_node> pt_tab(nb_points_of_convex(ic));
    ref_mesh_pt_ct ptab = points_of_convex(ic);
    pt_tab.assign(ptab.begin(), ptab.end());
    gs.cv = bgeot::convex<base_node>(structure_of_convex(ic), pt_tab);

    bgeot::pgeometric_trans pgt = gs.pgt = trans_of_convex(ic);

    size_type d = ipt.size() - 1, n = structure_of_convex(ic)->dim();
    static size_type d0 = 0;
    static bgeot::pstored_point_tab pspt1 = 0;
    mesh &mesh1
      = dal::singleton<mesh_cache_for_Bank_build_green_simplexes>::instance();
    if (d0 != d) {
      d0 = d;
      Bank_build_first_mesh(mesh1, d);
      pspt1 = bgeot::store_point_tab(mesh1.points());
    }

    const std::vector<size_type> &loc_ind = pgt->vertices();
    const bgeot::mesh_structure::ind_cv_ct &ct = ind_points_of_convex(ic);

    gs.ipt_loc.resize(ipt.size());
    std::vector<size_type> ipt_other;
    size_type nb_found = 0;
    for (size_type ip = 0; ip < loc_ind.size(); ++ip) {
      bool found = false;
      for (size_type i = 0; i < ipt.size(); ++i)
        if (ct[loc_ind[ip]] == ipt[i])
          { gs.ipt_loc[i] = ip; found = true; ++nb_found; break; }
      if (!found) ipt_other.push_back(ip);
    }
    assert(nb_found == ipt.size());

    mesh mesh2;
    for (size_type ip = 0; ip < loc_ind.size(); ++ip)
      mesh2.add_point(pgt->geometric_nodes()[loc_ind[ip]]);
    size_type ic1 = mesh2.add_simplex(dim_type(d), gs.ipt_loc.begin());
    bgeot::pgeometric_trans pgt1 = mesh2.trans_of_convex(ic1);
    for (size_type i = 0; i < ipt.size(); ++i)
      gs.ipt_loc[i] = loc_ind[gs.ipt_loc[i]];

    bgeot::pgeotrans_precomp pgp = bgeot::geotrans_precomp(pgt1, pspt1, 0);

    std::vector<size_type> ipt1(pspt1->size());
    base_node pt(dim());
    for (size_type i = 0; i < pspt1->size(); ++i) {
      pgp->transform(mesh2.points_of_convex(ic1), i, pt);
      ipt1[i] = mesh2.add_point(pt);
    }
    mesh2.sup_convex(ic1);

    std::vector<size_type> ipt2(n+1);
    for (size_type i = 0; i < mesh1.nb_convex(); ++i) {
      for (size_type j = 0; j <= d; ++j)
        ipt2[j] = ipt1[mesh1.ind_points_of_convex(i)[j]];
      for (size_type j = d+1; j <= n; ++j)
        ipt2[j] = ipt_other[j-d-1];
      mesh2.add_simplex(dim_type(n), ipt2.begin());
    }

    mesh mesh3;
    ipt1.resize(pgt->nb_points());
    for (dal::bv_visitor i(mesh2.convex_index()); !i.finished(); ++i) {
      bgeot::pgeometric_trans pgt2 = mesh2.trans_of_convex(i);
      for (size_type ip = 0; ip < pgt->nb_points(); ++ip)
        ipt1[ip] =
          mesh3.add_point(pgt2->transform(pgt->geometric_nodes()[ip],
                                          mesh2.points_of_convex(i)));
      mesh3.add_convex(pgt, ipt1.begin());
    }


    bgeot::pstored_point_tab pspt3 = bgeot::store_point_tab(mesh3.points());
    pgp = bgeot::geotrans_precomp(pgt, pspt3, 0);

    ipt1.resize(pspt3->size());
    for (size_type ip = 0; ip < pspt3->size(); ++ip) {
      pgp->transform(points_of_convex(ic), ip, pt);
      ipt1[ip] = add_point(pt);
    }
    // dal::del_stored_object(pspt3);

    ipt2.resize(pgt->nb_points());
    for (size_type icc = 0; icc < mesh3.nb_convex(); ++icc) {
      for (size_type j = 0; j < pgt->nb_points(); ++j)
        ipt2[j] = ipt1[mesh3.ind_points_of_convex(icc)[j]];
      size_type i = add_convex(pgt, ipt2.begin());
      gs.sub_simplices.push_back(i);
      Bank_info->is_green_simplex.add(i);
      Bank_info->num_green_simplex[i] = igs;
    }

    for (size_type ip1 = 0; ip1 < ipt.size(); ++ip1)
      for (size_type ip2 = ip1+1; ip2 < ipt.size(); ++ip2)
        Bank_info->edges.insert(edge(ipt[ip1], ipt[ip2]));

    handle_region_refinement(ic, gs.sub_simplices, true);
    sup_convex(ic, true);
  }

  void mesh::Bank_refine(dal::bit_vector b) {
    if (!(Bank_info.get())) Bank_info = std::make_unique<Bank_info_struct>();

    b &= convex_index();
    if (b.card() == 0) return;

    Bank_info->edges.clear();
    while (b.card() > 0)
      Bank_test_and_refine_convex(b.take_first(), b);

    std::vector<size_type> ipt;
    edge_set marked_convexes;
    while (Bank_info->edges.size()) {
      marked_convexes.clear();
      b = convex_index();
      edge_set::const_iterator it = Bank_info->edges.begin();
      edge_set::const_iterator ite = Bank_info->edges.end(), it2=it;
      assert(!Bank_info->edges.empty());
      for (; it != ite; it = it2) {
        if (it2 != ite) ++it2;
        Bank_convex_with_edge(it->i1, it->i2, ipt);
        if (ipt.size() == 0) ; // Bank_info->edges.erase(it);
        else for (size_type ic = 0; ic < ipt.size(); ++ic)
          marked_convexes.insert(edge(ipt[ic], it->i1, it->i2));
      }

      it = marked_convexes.begin(); ite = marked_convexes.end();
      if (it == ite) break;

      while (it != ite) {
        size_type ic = it->i0;
        ipt.resize(0);
        while (it != ite && it->i0 == ic) {
          bool found1 = false, found2 = false;
          for (size_type j = 0; j < ipt.size(); ++j) {
            if (ipt[j] == it->i1) found1 = true;
            if (ipt[j] == it->i2) found2 = true;
          }
          if (!found1) ipt.push_back(it->i1);
          if (!found2) ipt.push_back(it->i2);
          ++it;
        }
        if (b.is_in(ic)) {
          if (ipt.size() > structure_of_convex(ic)->dim())
            Bank_test_and_refine_convex(ic, b);
          else if (Bank_info->is_green_simplex[ic]) {
            size_type icc = Bank_test_and_refine_convex(ic, b, false);
            if (!Bank_is_convex_having_points(icc, ipt)) {
              Bank_test_and_refine_convex(icc, b);
            }
          }
          else Bank_build_green_simplexes(ic, ipt);
        }
      }
    }
    Bank_info->edges.clear();
  }

  struct dummy_mesh_ {
    mesh m;
    dummy_mesh_() : m() {}
  };

  const mesh &dummy_mesh()
  { return dal::singleton<dummy_mesh_>::instance().m; }

}  /* end of namespace getfem.                                             */
