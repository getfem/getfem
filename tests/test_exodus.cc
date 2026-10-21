/*===========================================================================

 Copyright (C) 2026 GetFEM contributors

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
 along  with  this program. If not, see https://www.gnu.org/licenses/.

===========================================================================*/

// Round-trip tests for the Exodus II import/export (getfem_exodus.{h,cc}):
// write a mesh, a field and regions, read them back through the generic
// import_mesh() dispatch, and check the geometry, field values and regions
// are recovered -- for a representative set of element types (including
// serendipity ones) and for an incrementally appended transient series.

#include "getfem/getfem_exodus.h"
#include "getfem/getfem_import.h"
#include "getfem/getfem_regular_meshes.h"
#include <algorithm>
#include <set>

using std::cout;
using std::endl;
using getfem::size_type;
using getfem::scalar_type;
using getfem::base_node;
using bgeot::short_type;
using bgeot::dim_type;

// an arbitrary smooth scalar field, sampled at node coordinates
static scalar_type field(const base_node &p, scalar_type t = 1.) {
  scalar_type v = 1. + 2.*p[0] + 0.3*p[0]*p[0];
  if (p.size() > 1) v += 3.*p[1] - 0.2*p[1]*p[1];
  if (p.size() > 2) v += 4.*p[2];
  return t*v;
}

// max |value at node i - field(node coord i)| over the re-read file
static scalar_type
node_value_error(const getfem::exodus_import &imp,
                 const std::vector<scalar_type> &U, scalar_type t) {
  const std::vector<base_node> &P = imp.node_points();
  scalar_type e = 0.;
  for (size_type i=0; i < U.size(); ++i)
    e = std::max(e, gmm::abs(U[i] - field(P[i], t)));
  return e;
}

template<class F>
static void expect_gmm_error(const std::string &what, F f,
                             const std::string &needle) {
  try { f(); }
  catch (const gmm::gmm_error &e) {
    GMM_ASSERT1(std::string(e.what()).find(needle) != std::string::npos,
                what << ": error message should contain '" << needle << "'");
    return;
  }
  GMM_ASSERT1(false, what << ": expected a GetFEM error");
}

// write a mesh_fem (classical of `order`, or the named serendipity `fem` if
// given, qdim components) + a field, read it back, and check the geometry
// and values; qdim>1 exercises the name_x/name_y component split
static void test_roundtrip(const std::string &gt, short_type order,
                           const std::vector<size_type> &nsub,
                           const char *fem = nullptr, size_type qdim = 1) {
  std::string tag = fem ? std::string(fem) : gt + " (order " + char('0'+order) + ")";
  const char *fname = "test_exodus_tmp.exo";

  getfem::mesh m;
  getfem::regular_unit_mesh(m, nsub, bgeot::geometric_trans_descriptor(gt));
  getfem::mesh_fem mf(m); mf.set_qdim(dim_type(qdim));
  if (fem)
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv)
      mf.set_finite_element(cv, getfem::fem_descriptor(fem));
  else
    mf.set_classical_finite_element(order);

  std::vector<scalar_type> U(mf.nb_dof());
  for (size_type i=0; i < mf.nb_dof(); ++i) U[i] = field(mf.point_of_basic_dof(i));

  { getfem::exodus_export ex(fname);
    ex.declare_point_data("u", qdim);
    ex.exporting(mf); ex.write_mesh();
    ex.write_point_data(mf, U, "u"); }              // closed by destructor

  getfem::mesh m2;
  getfem::exodus_import imp(fname);
  imp.read_mesh(m2);
  GMM_ASSERT1(m2.convex_index().card() == m.convex_index().card(),
              tag << ": number of elements");

  scalar_type err = 0.; size_type nnod = 0;
  for (size_type q=0; q < qdim; ++q) {              // read each component back
    std::string nm = qdim == 1 ? "u" : "u_" + std::string(1, char('x'+q));
    std::vector<scalar_type> Ub;
    imp.read_nodal_var(nm, 0, Ub);
    nnod = Ub.size();
    err = std::max(err, node_value_error(imp, Ub, 1.));
  }
  GMM_ASSERT1(nnod == imp.node_points().size(), tag << ": variable length");
  GMM_ASSERT1(err < 1e-10, tag << ": field values (err=" << err << ")");
  cout << "  " << tag << ": " << m2.convex_index().card() << " elements, "
       << nnod << " nodes, value error " << err << endl;
}

// export a mesh with a boundary face region AND a volume region, read it
// back (through the generic import_mesh() dispatch, and directly for
// node_set()), and check both regions and the node set survived
static void test_region_roundtrip() {
  const char *fname = "test_exodus_region.exo";
  getfem::mesh m;
  getfem::regular_unit_mesh(m, {2,2,2}, bgeot::geometric_trans_descriptor("GT_PK(3,1)"));
  getfem::outer_faces_of_mesh(m, m.region(7));          // boundary face region
  size_type n_orig = 0;
  for (getfem::mr_visitor i(m.region(7)); !i.finished(); ++i)
    if (i.is_face()) ++n_orig;
  { size_type c = 0;                                    // + a volume region
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv, ++c)
      if (c % 2 == 0) m.region(5).add(cv); }

  // Identify the volume region's convexes by centroid, not convex id: each
  // volume region becomes its own Exodus element block, numbered block by
  // block, so imported convex ids are legitimately renumbered.
  auto region_centroids = [](const getfem::mesh &mm, size_type rid) {
    std::set<std::vector<double> > s;
    for (getfem::mr_visitor i(mm.region(rid)); !i.finished(); ++i)
      if (!i.is_face()) {
        auto pts = mm.points_of_convex(i.cv());
        std::vector<double> g(mm.dim(), 0.);
        for (size_type k=0; k < pts.size(); ++k)
          for (size_type d=0; d < size_type(mm.dim()); ++d) g[d] += pts[k][d];
        for (double &x : g) x /= double(pts.size());
        s.insert(g);
      }
    return s;
  };
  std::set<std::vector<double> > orig = region_centroids(m, 5);

  { getfem::exodus_export ex(fname); ex.exporting(m); ex.write_mesh(); }

  getfem::mesh m2;                       // via the generic import_mesh() dispatch
  getfem::import_mesh(fname, "exodus", m2);
  GMM_ASSERT1(m2.convex_index().card() == m.convex_index().card(),
              "region roundtrip: element count");

  getfem::exodus_import imp(fname);      // direct, needed for node_set()
  getfem::mesh m3; imp.read_mesh(m3);

  GMM_ASSERT1(m3.regions_index().is_in(7), "region roundtrip: boundary region present");
  size_type n_imp = 0;
  for (getfem::mr_visitor i(m3.region(7)); !i.finished(); ++i)
    if (i.is_face()) ++n_imp;
  GMM_ASSERT1(n_imp == n_orig, "region roundtrip: boundary face count");

  GMM_ASSERT1(m3.regions_index().is_in(5), "region roundtrip: volume region present");
  std::set<std::vector<double> > got = region_centroids(m3, 5);
  GMM_ASSERT1(got == orig, "region roundtrip: volume convex set (by centroid)");
  GMM_ASSERT1(!imp.node_set(5).empty(), "region roundtrip: node set present");

  cout << "  region roundtrip: " << n_imp << "/" << n_orig << " boundary faces, "
       << got.size() << "/" << orig.size() << " volume convexes" << endl;
}

// a transient field written incrementally: create the file at step 0, then
// reopen with append=true to add each following step (the real path used by
// the scripting interface's 'append' option)
static void test_time_series_append() {
  const char *fname = "test_exodus_series.exo";
  std::vector<scalar_type> tvals = {0., 0.5, 1.0};

  getfem::mesh m;
  getfem::regular_unit_mesh(m, {3,3}, bgeot::geometric_trans_descriptor("GT_QK(2,1)"));
  getfem::mesh_fem mf(m);
  mf.set_classical_finite_element(2);              // -> QUAD9

  for (size_type s=0; s < tvals.size(); ++s) {
    std::vector<scalar_type> U(mf.nb_dof());
    for (size_type i=0; i < mf.nb_dof(); ++i)
      U[i] = field(mf.point_of_basic_dof(i), tvals[s]);
    getfem::exodus_export ex(fname, /* append = */ s > 0);
    if (s == 0) ex.declare_point_data("u");
    ex.exporting(mf);
    if (s == 0) ex.write_mesh();
    ex.set_time(tvals[s]);
    ex.write_point_data(mf, U, "u");
  }

  getfem::mesh m2;
  getfem::exodus_import imp(fname);
  imp.read_mesh(m2);
  GMM_ASSERT1(imp.nb_steps() == tvals.size(), "time series: number of steps");
  scalar_type maxerr = 0.;
  for (size_type s=0; s < tvals.size(); ++s) {
    GMM_ASSERT1(gmm::abs(imp.times()[s] - tvals[s]) < 1e-12,
                "time series: time value at step " << s);
    std::vector<scalar_type> Ub;
    imp.read_nodal_var("u", s, Ub);
    maxerr = std::max(maxerr, node_value_error(imp, Ub, tvals[s]));
  }
  GMM_ASSERT1(maxerr < 1e-10, "time series: field values per step");
  cout << "  time series (append): " << tvals.size() << " steps, max value error "
       << maxerr << endl;
}

// appending must be refused when the mesh/block layout changed since the
// file was created (the fingerprint guard)
static void test_append_fingerprint_mismatch() {
  const char *fname = "test_exodus_append_mismatch.exo";
  getfem::mesh m;
  getfem::regular_unit_mesh(m, {3,3}, bgeot::geometric_trans_descriptor("GT_QK(2,1)"));
  getfem::mesh_fem mf(m); mf.set_classical_finite_element(1);
  { getfem::exodus_export ex(fname); ex.exporting(mf); ex.write_mesh(); }

  size_type c = 0;                             // change the block layout
  for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv, ++c)
    m.region(c % 2 ? 6 : 5).add(cv);
  expect_gmm_error("append fingerprint mismatch", [&]() {
    getfem::exodus_export ex(fname, true);
    ex.exporting(mf);
  }, "fingerprint");
  cout << "  append fingerprint mismatch: refused append to changed element "
       << "blocks" << endl;
}

int main() {
  try {
    cout << "Exodus round-trip tests" << endl;
    test_roundtrip("GT_PK(2,1)", 2, {3,3});                           // TRI6
    test_roundtrip("GT_QK(3,1)", 1, {2,2,2});                         // HEX8
    test_roundtrip("GT_QK(2,1)", 1, {3,3}, nullptr, 2);                // vector field
    test_roundtrip("GT_QK(3,1)", 0, {2,2,2}, "FEM_Q2_INCOMPLETE(3)");  // HEX20 (serendipity)
    test_region_roundtrip();
    test_time_series_append();
    test_append_fingerprint_mismatch();
    cout << "all Exodus round-trip tests passed" << endl;
  } GMM_STANDARD_CATCH_ERROR;
  return 0;
}
