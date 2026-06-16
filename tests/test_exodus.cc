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

// Round-trip tests for the Exodus II import/export (getfem_exodus.{h,cc}).
//   * for several element types: write a mesh + a nodal field, read it back,
//     and check the geometry and the field values are recovered;
//   * a transient (time-series) field written to a single Exodus file and
//     read back step by step.

#include "getfem/getfem_exodus.h"
#include "getfem/getfem_regular_meshes.h"
#include <algorithm>
#include <set>

#ifdef GETFEM_HAVE_EXODUS

#include <netcdf.h>      // raw NetCDF reads to check Exodus entities directly
#include <sstream>

using std::cout;
using std::endl;
using getfem::size_type;
using getfem::scalar_type;
using getfem::base_node;
using bgeot::dim_type;
using bgeot::short_type;

static int nb_errors = 0;
static void check(bool ok, const std::string &what) {
  if (!ok) { ++nb_errors; cout << "  *** error has been detected: " << what << endl; }
}

template<class F>
static void expect_gmm_error(const std::string &what, F f,
                             const std::string &needle) {
  bool got = false;
  try { f(); }
  catch (const gmm::gmm_error &e) {
    got = true;
    std::string msg = e.what();
    check(needle.empty() || msg.find(needle) != std::string::npos,
          what + ": error message contains '" + needle + "'");
  }
  catch (...) {
    got = true;
    check(false, what + ": unexpected exception type");
  }
  check(got, what + ": expected GetFEM error");
}

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

static void test_roundtrip(const std::string &gt, short_type order,
                           const std::vector<size_type> &nsub) {
  std::string tag = gt + " (order " + char('0'+order) + ")";
  const char *fname = "test_exodus_tmp.exo";

  getfem::mesh m;
  getfem::regular_unit_mesh(m, nsub, bgeot::geometric_trans_descriptor(gt));
  getfem::mesh_fem mf(m);
  mf.set_classical_finite_element(order);

  std::vector<scalar_type> U(mf.nb_dof());
  for (size_type i=0; i < mf.nb_dof(); ++i) U[i] = field(mf.point_of_basic_dof(i));

  { getfem::exodus_export ex(fname);
    ex.exporting(mf); ex.write_mesh();
    ex.write_point_data(mf, U, "u"); }              // closed by destructor

  getfem::mesh m2;
  getfem::exodus_import imp(fname);
  imp.read_mesh(m2);

  check(m2.convex_index().card() == m.convex_index().card(),
        tag + ": number of elements");

  std::vector<scalar_type> Ub;
  imp.read_nodal_var("u", 0, Ub);
  check(Ub.size() == imp.node_points().size(), tag + ": variable length");
  scalar_type err = node_value_error(imp, Ub, 1.);
  check(err < 1e-10, tag + ": field values (err=" + std::to_string(err) + ")");
  cout << "  " << tag << ": " << m2.convex_index().card() << " elements, "
       << Ub.size() << " nodes, value error " << err << endl;
}

static void test_time_series() {
  const char *fname = "test_exodus_series.exo";
  const size_type nstep = 5;
  std::vector<scalar_type> tvals = {0., 0.25, 0.5, 0.75, 1.0};

  getfem::mesh m;
  getfem::regular_unit_mesh(m, {3,3}, bgeot::geometric_trans_descriptor("GT_QK(2,1)"));
  getfem::mesh_fem mf(m);
  mf.set_classical_finite_element(2);              // -> QUAD9

  // write the transient field to a single file
  { getfem::exodus_export ex(fname);
    ex.exporting(mf); ex.write_mesh();
    for (size_type s=0; s < nstep; ++s) {
      std::vector<scalar_type> U(mf.nb_dof());
      for (size_type i=0; i < mf.nb_dof(); ++i)
        U[i] = field(mf.point_of_basic_dof(i), tvals[s]);
      ex.set_time(tvals[s]);
      ex.write_point_data(mf, U, "u");
      ex.sync();
      getfem::exodus_import live(fname);
      check(live.nb_steps() == s+1, "time series live sync: visible step count");
    }
  }

  getfem::mesh m2;
  getfem::exodus_import imp(fname);
  imp.read_mesh(m2);
  check(imp.nb_steps() == nstep, "time series: number of steps");
  scalar_type terr = 0.;
  for (size_type s=0; s < imp.times().size(); ++s)
    terr = std::max(terr, gmm::abs(imp.times()[s] - tvals[s]));
  check(terr < 1e-12, "time series: time values");

  scalar_type maxerr = 0.;
  for (size_type s=0; s < nstep; ++s) {
    std::vector<scalar_type> Ub;
    imp.read_nodal_var("u", s, Ub);
    maxerr = std::max(maxerr, node_value_error(imp, Ub, tvals[s]));
  }
  check(maxerr < 1e-10, "time series: field values per step");
  cout << "  time series: " << nstep << " steps, max value error "
       << maxerr << endl;
}

// export a mesh whose outer faces form region 7, read it back, and check the
// boundary side set survived as a face region
static void test_region_roundtrip(const std::string &gt,
                                   const std::vector<size_type> &nsub) {
  std::string tag = gt + " regions";
  const char *fname = "test_exodus_region.exo";
  getfem::mesh m;
  getfem::regular_unit_mesh(m, nsub, bgeot::geometric_trans_descriptor(gt));
  getfem::outer_faces_of_mesh(m, m.region(7));
  size_type n_orig = 0;
  for (getfem::mr_visitor i(m.region(7)); !i.finished(); ++i)
    if (i.is_face()) ++n_orig;

  { getfem::exodus_export ex(fname); ex.exporting(m); ex.write_mesh(); }

  getfem::mesh m2;
  getfem::import_mesh_exodus(fname, m2);
  check(m2.regions_index().is_in(7), tag + ": region present");

  size_type n_imp = 0;
  for (getfem::mr_visitor i(m2.region(7)); !i.finished(); ++i)
    if (i.is_face()) ++n_imp;
  check(n_imp == n_orig, tag + ": boundary face count");

  // every face of the imported region must be a genuine boundary face of m2
  getfem::outer_faces_of_mesh(m2, m2.region(8));
  std::set<std::pair<size_type,short_type> > outer;
  for (getfem::mr_visitor i(m2.region(8)); !i.finished(); ++i)
    if (i.is_face()) outer.insert(std::make_pair(i.cv(), i.f()));
  bool all_outer = true;
  for (getfem::mr_visitor i(m2.region(7)); !i.finished(); ++i)
    if (i.is_face() && !outer.count(std::make_pair(i.cv(), i.f())))
      all_outer = false;
  check(all_outer, tag + ": region faces are genuine boundary faces");
  cout << "  " << tag << ": " << n_imp << "/" << n_orig
       << " boundary faces round-tripped" << endl;
}

// export a mesh with a convex (volume) region, read it back, and check the
// region (Exodus element set) and its node set survived
static void test_volume_region_roundtrip(const std::string &gt,
                                          const std::vector<size_type> &nsub) {
  std::string tag = gt + " volume region";
  const char *fname = "test_exodus_vol.exo";
  getfem::mesh m;
  getfem::regular_unit_mesh(m, nsub, bgeot::geometric_trans_descriptor(gt));
  { size_type c = 0;
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv, ++c)
      if (c % 2 == 0) m.region(5).add(cv); }

  // Identify the region's convexes by their centroid, not by convex id: the
  // export gives each volume region its own Exodus element block, and Exodus
  // numbers elements block by block, so the imported convex ids are legitimately
  // renumbered. Coordinates round-trip exactly, so centroid sets must match.
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

  getfem::mesh m2;
  getfem::exodus_import imp(fname);
  imp.read_mesh(m2);
  check(m2.regions_index().is_in(5), tag + ": region present");
  std::set<std::vector<double> > got = region_centroids(m2, 5);
  check(got == orig, tag + ": convex set (by centroid)");
  check(!imp.node_set(5).empty(), tag + ": node set present");
  cout << "  " << tag << ": " << got.size() << "/" << orig.size()
       << " convexes, " << imp.node_set(5).size() << " nodes in node set" << endl;
}

// reopen with the raw NetCDF API and check the "region" element variable, the
// identity element/node number maps, and named element sets
static void test_cell_var_and_names(const std::string &gt,
                                    const std::vector<size_type> &nsub) {
  std::string tag = gt + " cell var/names";
  const char *fname = "test_exodus_cell.exo";
  getfem::mesh m;
  getfem::regular_unit_mesh(m, nsub, bgeot::geometric_trans_descriptor(gt));
  { size_type c = 0;                       // two disjoint volume regions
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv, ++c)
      m.region(c % 2 ? 6 : 5).add(cv); }
  getfem::mesh_fem mf(m);
  mf.set_classical_finite_element(1);
  std::vector<scalar_type> U(mf.nb_dof());
  for (size_type i=0; i < mf.nb_dof(); ++i) U[i] = field(mf.point_of_basic_dof(i));
  { getfem::exodus_export ex(fname);
    ex.enable_region_field();                 // opt in to the "region" cell var
    ex.set_region_name(5, "core"); ex.set_region_name(6, "shell");
    ex.exporting(mf); ex.write_mesh(); ex.write_point_data(mf, U, "u"); }

  int nc = -1;
  check(nc_open(fname, NC_NOWRITE, &nc) == NC_NOERR, tag + ": reopen file");
  if (nc < 0) return;

  size_t lname = 32;
  { int dd; if (nc_inq_dimid(nc, "len_name", &dd) == NC_NOERR)
              nc_inq_dimlen(nc, dd, &lname); }
  auto read_name = [&](int v, size_t row) {
    std::vector<char> buf(lname + 1, '\0');
    size_t st[2] = { row, 0 }, ct[2] = { 1, lname };
    nc_get_vara_text(nc, v, st, ct, buf.data());
    return std::string(buf.data());
  };

  // (1) exactly one element variable, named "region"
  int d; size_t nev = 0;
  bool has_ev = (nc_inq_dimid(nc, "num_elem_var", &d) == NC_NOERR);
  if (has_ev) nc_inq_dimlen(nc, d, &nev);
  check(has_ev && nev == 1, tag + ": one element variable");
  { int v; if (nc_inq_varid(nc, "name_elem_var", &v) == NC_NOERR)
      check(read_name(v, 0) == "region", tag + ": element var named 'region'"); }

  // (2) every element's region value equals its block id
  int d_blk; size_t nblk = 0;
  nc_inq_dimid(nc, "num_el_blk", &d_blk); nc_inq_dimlen(nc, d_blk, &nblk);
  std::vector<int> ebp(nblk);
  { int v; nc_inq_varid(nc, "eb_prop1", &v); nc_get_var_int(nc, v, ebp.data()); }
  bool vals_ok = true;
  for (size_t b=0; b < nblk; ++b) {
    std::ostringstream sv, sd;
    sv << "vals_elem_var1eb" << (b+1); sd << "num_el_in_blk" << (b+1);
    int vv, dd; size_t ne = 0;
    nc_inq_dimid(nc, sd.str().c_str(), &dd); nc_inq_dimlen(nc, dd, &ne);
    nc_inq_varid(nc, sv.str().c_str(), &vv);
    std::vector<double> vals(ne);
    size_t st[2] = { 0, 0 }, ct[2] = { 1, ne };
    nc_get_vara_double(nc, vv, st, ct, vals.data());
    for (size_t e=0; e < ne; ++e) if (int(vals[e]) != ebp[b]) vals_ok = false;
  }
  check(vals_ok, tag + ": region value == block id");

  // (3) identity element/node number maps
  auto identity = [&](const char *name, const char *dimn) {
    int v, dd; size_t n = 0;
    if (nc_inq_varid(nc, name, &v) != NC_NOERR) return false;
    nc_inq_dimid(nc, dimn, &dd); nc_inq_dimlen(nc, dd, &n);
    std::vector<int> a(n); nc_get_var_int(nc, v, a.data());
    for (size_t i=0; i < n; ++i) if (a[i] != int(i)+1) return false;
    return true;
  };
  check(identity("elem_num_map", "num_elem"),  tag + ": elem_num_map identity");
  check(identity("node_num_map", "num_nodes"), tag + ": node_num_map identity");

  // (4) element-set names (region 5 -> "core", region 6 -> "shell")
  { int dd; size_t nes = 0;
    nc_inq_dimid(nc, "num_elem_sets", &dd); nc_inq_dimlen(nc, dd, &nes);
    std::vector<int> esp(nes); int v, vn;
    nc_inq_varid(nc, "els_prop1", &v); nc_get_var_int(nc, v, esp.data());
    bool names_ok = (nc_inq_varid(nc, "els_names", &vn) == NC_NOERR);
    for (size_t i=0; names_ok && i < nes; ++i) {
      std::string want = (esp[i]==5) ? "core" : (esp[i]==6) ? "shell" : "";
      if (read_name(vn, i) != want) names_ok = false;
    }
    check(names_ok, tag + ": element-set names");
  }
  nc_close(nc);

  // (5) a mesh-only export (no field -> no time step) has no element variable
  { const char *f2 = "test_exodus_meshonly.exo";
    { getfem::exodus_export ex(f2); ex.exporting(m); ex.write_mesh(); }
    int n2; nc_open(f2, NC_NOWRITE, &n2);
    int dd;
    check(nc_inq_dimid(n2, "num_elem_var", &dd) == NC_EBADDIM,
          tag + ": mesh-only export has no element variable");
    nc_close(n2);
  }

  // (6) mesh-only export with explicit region_field writes one static step
  { const char *f3 = "test_exodus_meshonly_region.exo";
    { getfem::exodus_export ex(f3);
      ex.enable_region_field(); ex.exporting(m); ex.write_mesh(); ex.close(); }
    int n3; nc_open(f3, NC_NOWRITE, &n3);
    int dd; size_t nev_region = 0, nt = 0;
    check(nc_inq_dimid(n3, "num_elem_var", &dd) == NC_NOERR,
          tag + ": mesh-only region field has elem var dim");
    nc_inq_dimlen(n3, dd, &nev_region);
    check(nev_region == 1, tag + ": mesh-only region field has one elem var");
    nc_inq_dimid(n3, "time_step", &dd); nc_inq_dimlen(n3, dd, &nt);
    check(nt == 1, tag + ": mesh-only region field has one time step");
    nc_close(n3);
  }
  cout << "  " << tag << ": region cell var, id maps, set names OK" << endl;
}

// import a hand-written legacy file using the packed coord(num_dim,num_nodes)
// layout (our exporter only writes split coordx/y/z, so this needs raw NetCDF)
static void test_packed_coord_import() {
  std::string tag = "packed coord import";
  const char *fname = "test_exodus_packed.exo";
  int nc, st = 0;          // deterministic setup; accumulate any NetCDF error
  st |= nc_create(fname, NC_CLOBBER | NC_64BIT_OFFSET, &nc);
  int d_dim, d_nodes, d_blk, d_ne, d_npe, d_time;
  st |= nc_def_dim(nc, "num_dim", 2, &d_dim);
  st |= nc_def_dim(nc, "num_nodes", 4, &d_nodes);
  st |= nc_def_dim(nc, "num_el_blk", 1, &d_blk);
  st |= nc_def_dim(nc, "num_el_in_blk1", 1, &d_ne);
  st |= nc_def_dim(nc, "num_nod_per_el1", 4, &d_npe);
  st |= nc_def_dim(nc, "time_step", NC_UNLIMITED, &d_time);
  int v_coord, v_conn, v_ebp;
  { int dd[2] = { d_dim, d_nodes };
    st |= nc_def_var(nc, "coord", NC_DOUBLE, 2, dd, &v_coord); }
  { int dd[2] = { d_ne, d_npe };
    st |= nc_def_var(nc, "connect1", NC_INT, 2, dd, &v_conn); }
  st |= nc_put_att_text(nc, v_conn, "elem_type", 5, "QUAD4");
  st |= nc_def_var(nc, "eb_prop1", NC_INT, 1, &d_blk, &v_ebp);
  st |= nc_enddef(nc);
  double coord[8] = { 0.,1.,1.,0.,  0.,0.,1.,1. };  // packed: all x, then all y
  st |= nc_put_var_double(nc, v_coord, coord);
  int conn[4] = { 1, 2, 3, 4 };
  st |= nc_put_var_int(nc, v_conn, conn);
  int ebp[1] = { 1 };
  st |= nc_put_var_int(nc, v_ebp, ebp);
  st |= nc_close(nc);
  check(st == NC_NOERR, tag + ": wrote legacy fixture");

  getfem::mesh m;
  getfem::exodus_import imp(fname);
  imp.read_mesh(m);
  const std::vector<base_node> &P = imp.node_points();
  check(P.size() == 4 && m.convex_index().card() == 1,
        tag + ": 1 quad, 4 nodes");
  check(gmm::abs(P[0][0]-0.) < 1e-12 && gmm::abs(P[0][1]-0.) < 1e-12 &&
        gmm::abs(P[2][0]-1.) < 1e-12 && gmm::abs(P[2][1]-1.) < 1e-12,
        tag + ": coordinates recovered from packed layout");
  expect_gmm_error(tag + ": append refused without fingerprint", [&]() {
    getfem::mesh_fem mf(m); mf.set_classical_finite_element(1);
    getfem::exodus_export ex(fname, true);
    ex.exporting(mf);
  }, "fingerprint");
  cout << "  " << tag << ": " << m.convex_index().card() << " element, "
       << P.size() << " nodes" << endl;
}

static void test_import_build_mesh_fem_node_order() {
  std::string tag = "import build_mesh_fem node order";
  const char *fname = "test_exodus_node_order.exo";
  int nc, st = 0;
  st |= nc_create(fname, NC_CLOBBER | NC_64BIT_OFFSET, &nc);
  int d_dim, d_nodes, d_blk, d_ne, d_npe, d_time, d_name, d_nv;
  st |= nc_def_dim(nc, "num_dim", 2, &d_dim);
  st |= nc_def_dim(nc, "num_nodes", 4, &d_nodes);
  st |= nc_def_dim(nc, "num_el_blk", 1, &d_blk);
  st |= nc_def_dim(nc, "num_el_in_blk1", 1, &d_ne);
  st |= nc_def_dim(nc, "num_nod_per_el1", 4, &d_npe);
  st |= nc_def_dim(nc, "time_step", NC_UNLIMITED, &d_time);
  st |= nc_def_dim(nc, "len_name", 33, &d_name);
  st |= nc_def_dim(nc, "num_nod_var", 1, &d_nv);
  int v_x, v_y, v_conn, v_ebp, v_time, v_nvn, v_u;
  st |= nc_def_var(nc, "coordx", NC_DOUBLE, 1, &d_nodes, &v_x);
  st |= nc_def_var(nc, "coordy", NC_DOUBLE, 1, &d_nodes, &v_y);
  { int dd[2] = { d_ne, d_npe };
    st |= nc_def_var(nc, "connect1", NC_INT, 2, dd, &v_conn); }
  st |= nc_put_att_text(nc, v_conn, "elem_type", 5, "QUAD4");
  st |= nc_def_var(nc, "eb_prop1", NC_INT, 1, &d_blk, &v_ebp);
  st |= nc_def_var(nc, "time_whole", NC_DOUBLE, 1, &d_time, &v_time);
  { int dd[2] = { d_nv, d_name };
    st |= nc_def_var(nc, "name_nod_var", NC_CHAR, 2, dd, &v_nvn); }
  { int dd[2] = { d_time, d_nodes };
    st |= nc_def_var(nc, "vals_nod_var1", NC_DOUBLE, 2, dd, &v_u); }
  st |= nc_enddef(nc);
  // Exodus node order intentionally differs from the QUAD4 connectivity order:
  // nodes 1..4 are top-left, bottom-left, top-right, bottom-right.
  double x[4] = { 0., 0., 1., 1. }, y[4] = { 1., 0., 1., 0. };
  int conn[4] = { 2, 4, 3, 1 }, ebp[1] = { 1 };
  double t[1] = { 0. }, vals[4] = { 10., 20., 30., 40. };
  char name[33]; std::fill(name, name+33, '\0'); name[0] = 'u';
  size_t ns[2] = { 0, 0 }, nc2[2] = { 1, 33 };
  size_t us[2] = { 0, 0 }, uc[2] = { 1, 4 };
  st |= nc_put_var_double(nc, v_x, x);
  st |= nc_put_var_double(nc, v_y, y);
  st |= nc_put_var_int(nc, v_conn, conn);
  st |= nc_put_var_int(nc, v_ebp, ebp);
  st |= nc_put_var_double(nc, v_time, t);
  st |= nc_put_vara_text(nc, v_nvn, ns, nc2, name);
  st |= nc_put_vara_double(nc, v_u, us, uc, vals);
  st |= nc_close(nc);
  check(st == NC_NOERR, tag + ": wrote fixture");

  getfem::mesh m;
  getfem::exodus_import imp(fname);
  imp.read_mesh(m);
  getfem::mesh_fem mf;
  imp.build_mesh_fem(m, mf);
  std::vector<scalar_type> U;
  imp.read_nodal_var("u", 0, U);
  check(U.size() == mf.nb_dof() && mf.nb_dof() == 4,
        tag + ": nodal variable length matches reduced mesh_fem dofs");
  std::vector<scalar_type> Ub(mf.nb_basic_dof());
  mf.extend_vector(U, Ub);
  size_type cv = m.convex_index().first_true();
  getfem::mesh_fem::ind_dof_ct dofs = mf.ind_basic_dof_of_element(cv);
  const bgeot::mesh_structure::ind_cv_ct &pts = m.ind_points_of_convex(cv);
  bool order_ok = true;
  for (size_type j=0; j < pts.size(); ++j)
    if (gmm::abs(Ub[dofs[j]] - vals[pts[j]]) > 1e-12) order_ok = false;
  check(order_ok, tag + ": reduced dofs follow Exodus node order");
  cout << "  " << tag << ": reduced dofs match Exodus node order" << endl;
}

static void test_append_block_layout_mismatch() {
  std::string tag = "append block layout mismatch";
  const char *fname = "test_exodus_append_mismatch.exo";
  getfem::mesh m;
  getfem::regular_unit_mesh(m, {3,3}, bgeot::geometric_trans_descriptor("GT_QK(2,1)"));
  getfem::mesh_fem mf(m); mf.set_classical_finite_element(1);
  std::vector<scalar_type> U(mf.nb_dof());
  for (size_type i=0; i < mf.nb_dof(); ++i) U[i] = field(mf.point_of_basic_dof(i));
  { getfem::exodus_export ex(fname);
    ex.enable_region_field(); ex.exporting(mf); ex.write_mesh();
    ex.write_point_data(mf, U, "u"); ex.close(); }

  size_type c = 0;
  for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv, ++c)
    m.region(c % 2 ? 6 : 5).add(cv);
  expect_gmm_error(tag, [&]() {
    getfem::exodus_export ex(fname, true);
    ex.exporting(mf);
  }, "fingerprint");
  cout << "  " << tag << ": refused append to changed element blocks" << endl;
}

static void test_incomplete_step_error_timing() {
  std::string tag = "incomplete transient step";
  const char *fname = "test_exodus_incomplete.exo";
  getfem::mesh m;
  getfem::regular_unit_mesh(m, {2,2}, bgeot::geometric_trans_descriptor("GT_QK(2,1)"));
  getfem::mesh_fem mf(m); mf.set_classical_finite_element(1);
  std::vector<scalar_type> U(mf.nb_dof()), V(mf.nb_dof());
  for (size_type i=0; i < mf.nb_dof(); ++i) {
    U[i] = field(mf.point_of_basic_dof(i));
    V[i] = 2.*U[i];
  }

  { getfem::exodus_export ex(fname);
    ex.exporting(mf); ex.declare_point_data("u"); ex.declare_point_data("v");
    ex.write_mesh();
    ex.set_time(0.); ex.write_point_data(mf, U, "u");
    ex.write_point_data(mf, V, "v");
    ex.set_time(1.); ex.write_point_data(mf, U, "u");
    expect_gmm_error(tag + ": next set_time", [&]() { ex.set_time(2.); },
                     "every field");
  }

  { getfem::exodus_export ex(fname);
    ex.exporting(mf); ex.declare_point_data("u"); ex.declare_point_data("v");
    ex.write_mesh();
    ex.set_time(0.); ex.write_point_data(mf, U, "u");
    expect_gmm_error(tag + ": close", [&]() { ex.close(); }, "every field");
  }
  cout << "  " << tag << ": rejected before syncing/closing incomplete steps" << endl;
}

static void test_shell_sideset_import() {
  std::string tag = "shell side-set import";
  const char *fname = "test_exodus_shell_sideset.exo";
  int nc, st = 0;
  st |= nc_create(fname, NC_CLOBBER | NC_64BIT_OFFSET, &nc);
  int d_dim, d_nodes, d_blk, d_ne, d_npe, d_nss, d_nside;
  st |= nc_def_dim(nc, "num_dim", 2, &d_dim);
  st |= nc_def_dim(nc, "num_nodes", 4, &d_nodes);
  st |= nc_def_dim(nc, "num_el_blk", 1, &d_blk);
  st |= nc_def_dim(nc, "num_el_in_blk1", 1, &d_ne);
  st |= nc_def_dim(nc, "num_nod_per_el1", 4, &d_npe);
  st |= nc_def_dim(nc, "num_side_sets", 1, &d_nss);
  st |= nc_def_dim(nc, "num_side_ss1", 1, &d_nside);
  int v_x, v_y, v_conn, v_ebp, v_ssp, v_e, v_s;
  st |= nc_def_var(nc, "coordx", NC_DOUBLE, 1, &d_nodes, &v_x);
  st |= nc_def_var(nc, "coordy", NC_DOUBLE, 1, &d_nodes, &v_y);
  { int dd[2] = { d_ne, d_npe };
    st |= nc_def_var(nc, "connect1", NC_INT, 2, dd, &v_conn); }
  st |= nc_put_att_text(nc, v_conn, "elem_type", 6, "SHELL4");
  st |= nc_def_var(nc, "eb_prop1", NC_INT, 1, &d_blk, &v_ebp);
  st |= nc_def_var(nc, "ss_prop1", NC_INT, 1, &d_nss, &v_ssp);
  st |= nc_def_var(nc, "elem_ss1", NC_INT, 1, &d_nside, &v_e);
  st |= nc_def_var(nc, "side_ss1", NC_INT, 1, &d_nside, &v_s);
  st |= nc_enddef(nc);
  double x[4] = {0.,1.,1.,0.}, y[4] = {0.,0.,1.,1.};
  int conn[4] = {1,2,3,4}, one[1] = {1};
  st |= nc_put_var_double(nc, v_x, x); st |= nc_put_var_double(nc, v_y, y);
  st |= nc_put_var_int(nc, v_conn, conn);
  st |= nc_put_var_int(nc, v_ebp, one); st |= nc_put_var_int(nc, v_ssp, one);
  st |= nc_put_var_int(nc, v_e, one); st |= nc_put_var_int(nc, v_s, one);
  st |= nc_close(nc);
  check(st == NC_NOERR, tag + ": wrote shell fixture");

  expect_gmm_error(tag, [&]() {
    getfem::mesh m;
    getfem::exodus_import imp(fname);
    imp.read_mesh(m);
  }, "SHELL");
  cout << "  " << tag << ": unsupported shell side sets are rejected" << endl;
}

// serendipity (incomplete) elements QUAD8 / HEX20 / WEDGE15 round-trip
static void test_serendipity_roundtrip(const std::string &gt, const char *fem,
                                       const std::vector<size_type> &nsub) {
  std::string tag = std::string(fem);
  const char *fname = "test_exodus_ser.exo";
  getfem::mesh m;
  getfem::regular_unit_mesh(m, nsub, bgeot::geometric_trans_descriptor(gt));
  getfem::mesh_fem mf(m);
  for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv)
    mf.set_finite_element(cv, getfem::fem_descriptor(fem));
  std::vector<scalar_type> U(mf.nb_dof());
  for (size_type i=0; i < mf.nb_dof(); ++i) U[i] = field(mf.point_of_basic_dof(i));
  { getfem::exodus_export ex(fname); ex.exporting(mf); ex.write_mesh();
    ex.write_point_data(mf, U, "u"); }
  getfem::mesh m2; getfem::exodus_import imp(fname); imp.read_mesh(m2);
  check(m2.convex_index().card() == m.convex_index().card(), tag + ": element count");
  std::vector<scalar_type> Ub; imp.read_nodal_var("u", 0, Ub);
  scalar_type err = node_value_error(imp, Ub, 1.);
  check(err < 1e-10, tag + ": field values (err=" + std::to_string(err) + ")");
  cout << "  " << tag << ": " << m2.convex_index().card() << " elements, "
       << Ub.size() << " nodes, value error " << err << endl;
}

// a vector (qdim>1) field is written as components name_x / name_y
static void test_vector_field_roundtrip() {
  const char *fname = "test_exodus_vec.exo";
  getfem::mesh m;
  getfem::regular_unit_mesh(m, {3,3}, bgeot::geometric_trans_descriptor("GT_QK(2,1)"));
  getfem::mesh_fem mf(m); mf.set_qdim(2); mf.set_classical_finite_element(1);
  std::vector<scalar_type> U(mf.nb_dof());
  for (size_type i=0; i < mf.nb_dof(); ++i) U[i] = field(mf.point_of_basic_dof(i));
  { getfem::exodus_export ex(fname); ex.exporting(mf); ex.write_mesh();
    ex.write_point_data(mf, U, "disp"); }
  getfem::mesh m2; getfem::exodus_import imp(fname); imp.read_mesh(m2);
  std::vector<scalar_type> Ux, Uy;
  imp.read_nodal_var("disp_x", 0, Ux);
  imp.read_nodal_var("disp_y", 0, Uy);
  scalar_type err = std::max(node_value_error(imp, Ux, 1.),
                             node_value_error(imp, Uy, 1.));
  check(err < 1e-10, "vector field: disp_x/_y values (err="
        + std::to_string(err) + ")");
  cout << "  vector field: disp_x/_y round-trip, value error " << err << endl;
}

// default files are compressed NetCDF4 classic-model when the NetCDF build
// supports deflate; 'uncompressed' is always classic 64-bit offset
static void test_compression_format() {
  getfem::mesh m;
  getfem::regular_unit_mesh(m, {3,3}, bgeot::geometric_trans_descriptor("GT_QK(2,1)"));
  getfem::mesh_fem mf(m); mf.set_classical_finite_element(1);
  std::vector<scalar_type> U(mf.nb_dof());
  for (size_type i=0; i < mf.nb_dof(); ++i) U[i] = field(mf.point_of_basic_dof(i));
  { getfem::exodus_export ex("test_exodus_c.exo");
    ex.exporting(mf); ex.write_mesh(); ex.write_point_data(mf, U, "u"); }
  { getfem::exodus_export ex("test_exodus_u.exo");
    ex.enable_compression(0); ex.exporting(mf); ex.write_mesh();
    ex.write_point_data(mf, U, "u"); }
  int nc = -1, fmt = 0;
  if (nc_open("test_exodus_c.exo", NC_NOWRITE, &nc) == NC_NOERR) {
    nc_inq_format(nc, &fmt); nc_close(nc);
#ifdef GETFEM_HAVE_EXODUS_NETCDF4
    check(fmt == NC_FORMAT_NETCDF4_CLASSIC,
          "compression: default file is NetCDF4 classic-model");
#else
    check(fmt == NC_FORMAT_64BIT_OFFSET,
          "compression: default file is 64-bit offset without NetCDF4 deflate");
#endif
  }
  if (nc_open("test_exodus_u.exo", NC_NOWRITE, &nc) == NC_NOERR) {
    nc_inq_format(nc, &fmt); nc_close(nc);
    check(fmt == NC_FORMAT_64BIT_OFFSET,
          "compression: 'uncompressed' file is 64-bit offset");
  }
  getfem::mesh mc;
  { getfem::exodus_import imp("test_exodus_c.exo"); imp.read_mesh(mc); }
  check(mc.convex_index().card() == m.convex_index().card(),
        "compression: default file re-imports");
#ifndef GETFEM_HAVE_EXODUS_NETCDF4
  expect_gmm_error("compression unavailable", [&]() {
    getfem::exodus_export ex("test_exodus_need_netcdf4.exo");
    ex.enable_compression();
  }, "NetCDF4");
#endif
#ifdef GETFEM_HAVE_EXODUS_NETCDF4
  cout << "  compression: default NetCDF4-classic, 'uncompressed' 64-bit offset"
       << endl;
#else
  cout << "  compression: default and 'uncompressed' both use 64-bit offset "
       << "without NetCDF4 deflate" << endl;
#endif
}

// malformed import files must be rejected with a clean error, not crash/UB
static void test_import_rejections() {
  auto write_fixture = [](const char *fname, size_t num_dim, int last_conn) {
    int nc, st = 0;
    st |= nc_create(fname, NC_CLOBBER | NC_64BIT_OFFSET, &nc);
    int d_dim, d_nodes, d_blk, d_ne, d_npe;
    st |= nc_def_dim(nc, "num_dim", num_dim, &d_dim);
    st |= nc_def_dim(nc, "num_nodes", 4, &d_nodes);
    st |= nc_def_dim(nc, "num_el_blk", 1, &d_blk);
    st |= nc_def_dim(nc, "num_el_in_blk1", 1, &d_ne);
    st |= nc_def_dim(nc, "num_nod_per_el1", 4, &d_npe);
    int v_coord, v_conn;
    { int dd[2] = { d_dim, d_nodes };
      st |= nc_def_var(nc, "coord", NC_DOUBLE, 2, dd, &v_coord); }
    { int dd[2] = { d_ne, d_npe };
      st |= nc_def_var(nc, "connect1", NC_INT, 2, dd, &v_conn); }
    st |= nc_put_att_text(nc, v_conn, "elem_type", 5, "QUAD4");
    st |= nc_enddef(nc);
    std::vector<double> coord(num_dim*4, 0.);
    st |= nc_put_var_double(nc, v_coord, coord.data());
    int conn[4] = { 1, 2, 3, last_conn };
    st |= nc_put_var_int(nc, v_conn, conn);
    st |= nc_close(nc);
    return st == NC_NOERR;
  };
  write_fixture("test_exodus_badconn.exo", 2, 99);     // node id out of range
  expect_gmm_error("import rejection (connectivity)", [&]() {
    getfem::mesh m; getfem::exodus_import imp("test_exodus_badconn.exo");
    imp.read_mesh(m);
  }, "out of range");
  write_fixture("test_exodus_baddim.exo", 4, 4);       // unsupported num_dim
  expect_gmm_error("import rejection (num_dim)", [&]() {
    getfem::mesh m; getfem::exodus_import imp("test_exodus_baddim.exo");
    imp.read_mesh(m);
  }, "dimension");

  auto write_set_fixture = [](const char *fname, const std::string &kind) {
    int nc, st = 0;
    st |= nc_create(fname, NC_CLOBBER | NC_64BIT_OFFSET, &nc);
    int d_dim, d_nodes, d_blk, d_ne, d_npe;
    st |= nc_def_dim(nc, "num_dim", 2, &d_dim);
    st |= nc_def_dim(nc, "num_nodes", 4, &d_nodes);
    st |= nc_def_dim(nc, "num_el_blk", 1, &d_blk);
    st |= nc_def_dim(nc, "num_el_in_blk1", 1, &d_ne);
    st |= nc_def_dim(nc, "num_nod_per_el1", 4, &d_npe);
    int v_x, v_y, v_conn;
    st |= nc_def_var(nc, "coordx", NC_DOUBLE, 1, &d_nodes, &v_x);
    st |= nc_def_var(nc, "coordy", NC_DOUBLE, 1, &d_nodes, &v_y);
    { int dd[2] = { d_ne, d_npe };
      st |= nc_def_var(nc, "connect1", NC_INT, 2, dd, &v_conn); }
    st |= nc_put_att_text(nc, v_conn, "elem_type", 5, "QUAD4");
    int v_bad_a = -1, v_bad_b = -1, d_set = -1, d_entries = -1;
    if (kind == "side_elem" || kind == "side_side") {
      st |= nc_def_dim(nc, "num_side_sets", 1, &d_set);
      st |= nc_def_dim(nc, "num_side_ss1", 1, &d_entries);
      st |= nc_def_var(nc, "elem_ss1", NC_INT, 1, &d_entries, &v_bad_a);
      st |= nc_def_var(nc, "side_ss1", NC_INT, 1, &d_entries, &v_bad_b);
    } else if (kind == "element") {
      st |= nc_def_dim(nc, "num_elem_sets", 1, &d_set);
      st |= nc_def_dim(nc, "num_ele_els1", 1, &d_entries);
      st |= nc_def_var(nc, "elem_els1", NC_INT, 1, &d_entries, &v_bad_a);
    } else if (kind == "node") {
      st |= nc_def_dim(nc, "num_node_sets", 1, &d_set);
      st |= nc_def_dim(nc, "num_nod_ns1", 1, &d_entries);
      st |= nc_def_var(nc, "node_ns1", NC_INT, 1, &d_entries, &v_bad_a);
    } else return false;
    st |= nc_enddef(nc);
    double x[4] = { 0., 1., 1., 0. };
    double y[4] = { 0., 0., 1., 1. };
    int conn[4] = { 1, 2, 3, 4 };
    st |= nc_put_var_double(nc, v_x, x);
    st |= nc_put_var_double(nc, v_y, y);
    st |= nc_put_var_int(nc, v_conn, conn);
    if (kind == "side_elem") {
      int elem = 99, side = 1;
      st |= nc_put_var_int(nc, v_bad_a, &elem);
      st |= nc_put_var_int(nc, v_bad_b, &side);
    } else if (kind == "side_side") {
      int elem = 1, side = 99;
      st |= nc_put_var_int(nc, v_bad_a, &elem);
      st |= nc_put_var_int(nc, v_bad_b, &side);
    } else if (kind == "element") {
      int elem = 99;
      st |= nc_put_var_int(nc, v_bad_a, &elem);
    } else {
      int node = 99;
      st |= nc_put_var_int(nc, v_bad_a, &node);
    }
    st |= nc_close(nc);
    return st == NC_NOERR;
  };
  check(write_set_fixture("test_exodus_bad_side_elem.exo", "side_elem"),
        "import rejection fixture: side-set element id");
  expect_gmm_error("import rejection (side-set element id)", [&]() {
    getfem::mesh m; getfem::exodus_import imp("test_exodus_bad_side_elem.exo");
    imp.read_mesh(m);
  }, "side set");
  check(write_set_fixture("test_exodus_bad_side_side.exo", "side_side"),
        "import rejection fixture: side-set side id");
  expect_gmm_error("import rejection (side-set side id)", [&]() {
    getfem::mesh m; getfem::exodus_import imp("test_exodus_bad_side_side.exo");
    imp.read_mesh(m);
  }, "side set");
  check(write_set_fixture("test_exodus_bad_elemset.exo", "element"),
        "import rejection fixture: element-set element id");
  expect_gmm_error("import rejection (element-set element id)", [&]() {
    getfem::mesh m; getfem::exodus_import imp("test_exodus_bad_elemset.exo");
    imp.read_mesh(m);
  }, "element set");
  check(write_set_fixture("test_exodus_bad_nodeset.exo", "node"),
        "import rejection fixture: node-set node id");
  expect_gmm_error("import rejection (node-set node id)", [&]() {
    getfem::mesh m; getfem::exodus_import imp("test_exodus_bad_nodeset.exo");
    imp.read_mesh(m);
  }, "node set");
  cout << "  import rejection: bad dimensions, connectivity and set ids are "
       << "rejected" << endl;
}

int main() {
  cout << "Exodus round-trip tests" << endl;
  test_roundtrip("GT_PK(1,1)",  1, {4});       // BAR2
  test_roundtrip("GT_PK(1,1)",  2, {4});       // BAR3
  test_roundtrip("GT_PK(2,1)",  1, {3,3});     // TRI3
  test_roundtrip("GT_PK(2,1)",  2, {3,3});     // TRI6
  test_roundtrip("GT_QK(2,1)",  1, {3,3});     // QUAD4
  test_roundtrip("GT_QK(2,1)",  2, {3,3});     // QUAD9
  test_roundtrip("GT_PK(3,1)",  1, {2,2,2});   // TETRA4
  test_roundtrip("GT_PK(3,1)",  2, {2,2,2});   // TETRA10
  test_roundtrip("GT_QK(3,1)",  1, {2,2,2});   // HEX8
  test_roundtrip("GT_QK(3,1)",  2, {2,2,2});   // HEX27
  test_roundtrip("GT_PRISM(3,1)", 1, {2,2,2}); // WEDGE6
  test_region_roundtrip("GT_PK(2,1)", {3,3});
  test_region_roundtrip("GT_QK(2,1)", {3,3});
  test_region_roundtrip("GT_PK(3,1)", {2,2,2});
  test_region_roundtrip("GT_QK(3,1)", {2,2,2});
  test_volume_region_roundtrip("GT_PK(2,1)", {3,3});
  test_volume_region_roundtrip("GT_QK(3,1)", {2,2,2});
  test_cell_var_and_names("GT_QK(2,1)", {3,3});
  test_cell_var_and_names("GT_PK(3,1)", {2,2,2});
  test_packed_coord_import();
  test_import_build_mesh_fem_node_order();
  test_append_block_layout_mismatch();
  test_incomplete_step_error_timing();
  test_shell_sideset_import();
  test_serendipity_roundtrip("GT_QK(2,1)", "FEM_Q2_INCOMPLETE(2)", {2,2});       // QUAD8
  test_serendipity_roundtrip("GT_QK(3,1)", "FEM_Q2_INCOMPLETE(3)", {2,2,2});     // HEX20
  test_serendipity_roundtrip("GT_PRISM(3,1)", "FEM_PRISM_INCOMPLETE_P2", {2,2,2}); // WEDGE15
  test_vector_field_roundtrip();
  test_compression_format();
  test_import_rejections();
  test_time_series();
  cout << (nb_errors ? "FAILED" : "all Exodus round-trip tests passed") << endl;
  return nb_errors ? 1 : 0;
}

#else
int main() { cout << "GetFEM built without Exodus support, test skipped\n"; return 0; }
#endif
