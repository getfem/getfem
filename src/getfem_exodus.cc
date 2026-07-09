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

#include "getfem/getfem_exodus.h"

#ifdef GETFEM_HAVE_EXODUS

#include <netcdf.h>
#include <ctime>
#include <cstdint>
#include <iomanip>
#include <limits>

namespace getfem {

  /* helper that turns a NetCDF status code into a GetFEM exception */
#define EXNC(call) do { int e_nc_ = (call);                                  \
    GMM_ASSERT1(e_nc_ == NC_NOERR, "NetCDF/Exodus error: "                    \
                << nc_strerror(e_nc_)); } while (0)

  /* maximum string length used by Exodus (MAX_STR_LENGTH + 1) */
  static const size_t EXO_STRLEN = 33;
  static const size_t EXO_IO_CHUNK = 65536;

  static void fnv_add(uint64_t &h, const void *p, size_t n) {
    const unsigned char *b = static_cast<const unsigned char *>(p);
    for (size_t i=0; i < n; ++i) { h ^= uint64_t(b[i]); h *= 1099511628211ULL; }
  }

  static void fnv_add_int(uint64_t &h, int64_t v) {
    fnv_add(h, &v, sizeof(v));
  }

  static void fnv_add_size(uint64_t &h, size_type v) {
    uint64_t u = uint64_t(v); fnv_add(h, &u, sizeof(u));
  }

  static void fnv_add_scalar(uint64_t &h, scalar_type x) {
    double d = double(x); fnv_add(h, &d, sizeof(d));
  }

  static void fnv_add_string(uint64_t &h, const std::string &s) {
    fnv_add_size(h, s.size());
    if (!s.empty()) fnv_add(h, s.data(), s.size());
  }

  static std::string fnv_hex(uint64_t h) {
    std::ostringstream s;
    s << std::hex << std::setw(16) << std::setfill('0') << h;
    return s.str();
  }

  static std::string exo_get_att_text(int ncid, int varid,
                                      const std::string &name) {
    size_t len = 0;
    if (nc_inq_attlen(ncid, varid, name.c_str(), &len) != NC_NOERR)
      return std::string();
    std::vector<char> buf(len + 1, '\0');
    EXNC(nc_get_att_text(ncid, varid, name.c_str(), buf.data()));
    return std::string(buf.data(), len);
  }

  /* read exactly `n` ints from `varid` into `out`. Using nc_get_vara_int with an
     explicit count bounds the read to the caller's buffer, so a malformed file
     whose variable is larger than its declared set/length cannot overflow it
     (and a too-small variable yields a clean NC_EEDGE error). */
  static void exo_get_ints(int ncid, int varid, size_t n, int *out) {
    if (n == 0) return;
    size_t start = 0, count = n;
    EXNC(nc_get_vara_int(ncid, varid, &start, &count, out));
  }

  static bool exo_is_shell_type(const std::string &t) {
    std::string u;
    for (char ch : t) u.push_back(char(std::toupper(ch)));
    return u.compare(0,5,"SHELL") == 0 || u.find("SHELL") != std::string::npos;
  }

  /* ********************************************************************* */
  /*  Element-type mapping between GetFEM and Exodus.                      */
  /*                                                                       */
  /*  The permutation between GetFEM's local node ordering and Exodus'     */
  /*  PATRAN ordering is derived by matching reference-element             */
  /*  coordinates, so the same machinery serves both writing and reading   */
  /*  and there is no hand-transcribed permutation table to get wrong.     */
  /* ********************************************************************* */

  enum exo_topo { EXO_BAR, EXO_TRI, EXO_QUAD, EXO_TET, EXO_HEX, EXO_WEDGE };

  static base_node exo_nd(scalar_type x) {
    base_node p(1); p[0] = x; return p; }
  static base_node exo_nd(scalar_type x, scalar_type y) {
    base_node p(2); p[0] = x; p[1] = y; return p; }
  static base_node exo_nd(scalar_type x, scalar_type y, scalar_type z) {
    base_node p(3); p[0] = x; p[1] = y; p[2] = z; return p; }

  static base_node exo_mid(const base_node &a, const base_node &b) {
    base_node r(a);
    for (size_type i=0; i < a.size(); ++i) r[i] = 0.5*(a[i]+b[i]);
    return r;
  }
  static base_node exo_center(const std::vector<base_node> &c,
                              std::initializer_list<int> idx) {
    base_node r(c[0].size()); r.fill(0.);
    for (int i : idx) for (size_type k=0; k < r.size(); ++k) r[k] += c[i][k];
    for (size_type k=0; k < r.size(); ++k) r[k] /= scalar_type(idx.size());
    return r;
  }

  /* Exodus corner nodes, in Exodus order, in GetFEM's unit reference frame. */
  static std::vector<base_node> exo_corners(exo_topo topo) {
    switch (topo) {
    case EXO_BAR: return { exo_nd(0.), exo_nd(1.) };
    case EXO_TRI: return { exo_nd(0,0), exo_nd(1,0), exo_nd(0,1) };
    case EXO_QUAD: return { exo_nd(0,0), exo_nd(1,0), exo_nd(1,1), exo_nd(0,1) };
    case EXO_TET: return { exo_nd(0,0,0), exo_nd(1,0,0),
                           exo_nd(0,1,0), exo_nd(0,0,1) };
    case EXO_HEX: return { exo_nd(0,0,0), exo_nd(1,0,0), exo_nd(1,1,0),
                           exo_nd(0,1,0), exo_nd(0,0,1), exo_nd(1,0,1),
                           exo_nd(1,1,1), exo_nd(0,1,1) };
    case EXO_WEDGE: return { exo_nd(0,0,0), exo_nd(1,0,0), exo_nd(0,1,0),
                             exo_nd(0,0,1), exo_nd(1,0,1), exo_nd(0,1,1) };
    }
    GMM_ASSERT1(false, "internal error");
  }

  /* Exodus edges (corner-index pairs), in Exodus edge order. */
  static std::vector<std::pair<int,int> > exo_edges(exo_topo topo) {
    switch (topo) {
    case EXO_BAR: return { {0,1} };
    case EXO_TRI: return { {0,1},{1,2},{2,0} };
    case EXO_QUAD: return { {0,1},{1,2},{2,3},{3,0} };
    case EXO_TET: return { {0,1},{1,2},{2,0},{0,3},{1,3},{2,3} };
    case EXO_HEX: return { {0,1},{1,2},{2,3},{3,0},   // bottom-face edges
                           {0,4},{1,5},{2,6},{3,7},   // vertical edges
                           {4,5},{5,6},{6,7},{7,4} }; // top-face edges
    case EXO_WEDGE: return { {0,1},{1,2},{2,0},   // bottom-triangle edges
                             {0,3},{1,4},{2,5},   // vertical edges
                             {3,4},{4,5},{5,3} }; // top-triangle edges
    }
    GMM_ASSERT1(false, "internal error");
  }

  /* Full ordered Exodus reference coordinates for the given topology / node
     count, in GetFEM's unit reference frame. */
  static std::vector<base_node> exo_reference(exo_topo topo, size_type nbd) {
    std::vector<base_node> c = exo_corners(topo);
    std::vector<base_node> r = c;
    if (r.size() == nbd) return r;                  // linear element
    for (const auto &e : exo_edges(topo))           // quadratic edge nodes
      r.push_back(exo_mid(c[e.first], c[e.second]));
    if (r.size() == nbd) return r;                  // serendipity (8,15,20)
    if (topo == EXO_QUAD && nbd == 9) {             // QUAD9 centre
      r.push_back(exo_center(c, {0,1,2,3})); return r;
    }
    if (topo == EXO_HEX && nbd == 27) {  // HEX27: body centre, then 6 face
      // centres, in the Exodus/PATRAN order verified against libMesh:
      // body, -z(bottom), +z(top), -x(left), +x(right), -y(front), +y(back).
      r.push_back(exo_center(c, {0,1,2,3,4,5,6,7})); // body
      r.push_back(exo_center(c, {0,1,2,3}));         // -z bottom
      r.push_back(exo_center(c, {4,5,6,7}));         // +z top
      r.push_back(exo_center(c, {0,3,4,7}));         // -x left
      r.push_back(exo_center(c, {1,2,5,6}));         // +x right
      r.push_back(exo_center(c, {0,1,4,5}));         // -y front
      r.push_back(exo_center(c, {2,3,6,7}));         // +y back
      return r;
    }
    GMM_ASSERT1(false, "unsupported Exodus node count " << nbd
                << " for this element topology");
  }

  static std::string exo_typestring(exo_topo topo, size_type nbd) {
    switch (topo) {
    case EXO_BAR:   return nbd==2 ? "BAR2" : "BAR3";
    case EXO_TRI:   return nbd==3 ? "TRI3" : "TRI6";
    case EXO_QUAD:  return nbd==4 ? "QUAD4" : (nbd==8 ? "QUAD8" : "QUAD9");
    case EXO_TET:   return nbd==4 ? "TETRA4" : "TETRA10";
    case EXO_HEX:   return nbd==8 ? "HEX8" : (nbd==20 ? "HEX20" : "HEX27");
    case EXO_WEDGE: return nbd==6 ? "WEDGE6" : "WEDGE15";
    }
    GMM_ASSERT1(false, "internal error");
  }

  /* classify a convex's linear topology from its (basic) structure */
  static exo_topo exo_topo_of_structure(bgeot::pconvex_structure cvs) {
    bgeot::pconvex_structure b = bgeot::basic_structure(cvs);
    dim_type d = b->dim();
    short_type nv = b->nb_points();
    if (d == 1) return EXO_BAR;
    if (d == 2) { if (nv == 3) return EXO_TRI; if (nv == 4) return EXO_QUAD; }
    if (d == 3) {
      if (nv == 4) return EXO_TET;
      if (nv == 6) return EXO_WEDGE;
      if (nv == 8) return EXO_HEX;
    }
    GMM_ASSERT1(false, "this element (dim " << int(d) << ", " << nv
                << " vertices) has no Exodus equivalent");
  }

  /* parse an Exodus element-type string into a topology */
  static exo_topo exo_topo_of_name(const std::string &t) {
    std::string u; for (char ch : t) u.push_back(char(std::toupper(ch)));
    if (u.compare(0,3,"BAR")==0 || u.compare(0,4,"BEAM")==0 ||
        u.compare(0,5,"TRUSS")==0) return EXO_BAR;
    if (u.compare(0,5,"TETRA")==0 || u.compare(0,3,"TET")==0) return EXO_TET;
    if (u.compare(0,3,"HEX")==0) return EXO_HEX;
    if (u.compare(0,5,"WEDGE")==0 || u.compare(0,5,"PRISM")==0) return EXO_WEDGE;
    // SHELL* is read as a 2D QUAD surface element. NB: a shell's side numbering
    // (top/bottom faces + edges) is NOT the 2D quad-edge numbering, so side sets
    // on shell elements from external files will not be translated correctly.
    if (u.compare(0,4,"QUAD")==0 || u.compare(0,5,"SHELL")==0) return EXO_QUAD;
    if (u.compare(0,3,"TRI")==0) return EXO_TRI; // after TETRA/TRUSS
    GMM_ASSERT1(false, "unsupported Exodus element type '" << t << "'");
  }

  static bgeot::pgeometric_trans exo_pgt(exo_topo topo, size_type nbd) {
    switch (topo) {
    case EXO_BAR:   return bgeot::simplex_geotrans(1, nbd==2 ? 1 : 2);
    case EXO_TRI:   return bgeot::simplex_geotrans(2, nbd==3 ? 1 : 2);
    case EXO_QUAD:
      if (nbd==4) return bgeot::parallelepiped_geotrans(2,1);
      if (nbd==8) return bgeot::Q2_incomplete_geotrans(2);
      return bgeot::parallelepiped_geotrans(2,2);
    case EXO_TET:   return bgeot::simplex_geotrans(3, nbd==4 ? 1 : 2);
    case EXO_HEX:
      if (nbd==8)  return bgeot::parallelepiped_geotrans(3,1);
      if (nbd==20) return bgeot::Q2_incomplete_geotrans(3);
      return bgeot::parallelepiped_geotrans(3,2);
    case EXO_WEDGE:
      if (nbd==6) return bgeot::prism_geotrans(3,1);
      return bgeot::prism_incomplete_P2_geotrans();
    }
    GMM_ASSERT1(false, "internal error");
  }

  /* perm[exodus_local] = gf_local index, by reference-coordinate matching */
  static std::vector<unsigned>
  exo_build_perm(const std::vector<base_node> &gfref,
                 const std::vector<base_node> &exoref) {
    GMM_ASSERT1(gfref.size() == exoref.size(),
                "Exodus node-count mismatch (" << gfref.size() << " vs "
                << exoref.size() << ")");
    std::vector<unsigned> perm(exoref.size());
    std::vector<bool> used(gfref.size(), false);
    for (size_type e=0; e < exoref.size(); ++e) {
      size_type best = size_type(-1); scalar_type bd = 1e30;
      for (size_type g=0; g < gfref.size(); ++g) if (!used[g]) {
        scalar_type d = gmm::vect_dist2(gfref[g], exoref[e]);
        if (d < bd) { bd = d; best = g; }
      }
      GMM_ASSERT1(best != size_type(-1) && bd < 1e-6,
                  "could not match GetFEM and Exodus reference nodes");
      perm[e] = unsigned(best); used[best] = true;
    }
    return perm;
  }


  static size_type exo_corner_count(exo_topo topo) {
    switch (topo) {
    case EXO_BAR: return 2; case EXO_TRI: return 3; case EXO_QUAD: return 4;
    case EXO_TET: return 4; case EXO_HEX: return 8; case EXO_WEDGE: return 6;
    }
    GMM_ASSERT1(false, "internal error");
  }

  /* Exodus sides (element faces) as corner-index lists, in Exodus side order. */
  static std::vector<std::vector<int> > exo_sides(exo_topo topo) {
    switch (topo) {
    case EXO_BAR:   return { {0}, {1} };
    case EXO_TRI:   return { {0,1},{1,2},{2,0} };
    case EXO_QUAD:  return { {0,1},{1,2},{2,3},{3,0} };
    case EXO_TET:   return { {0,1,3},{1,2,3},{0,3,2},{0,2,1} };
    case EXO_HEX:   return { {0,1,5,4},{1,2,6,5},{2,3,7,6},
                             {0,4,7,3},{0,3,2,1},{4,5,6,7} };
    case EXO_WEDGE: return { {0,1,4,3},{1,2,5,4},{0,3,5,2},{0,2,1},{3,4,5} };
    }
    GMM_ASSERT1(false, "internal error");
  }

  static base_node exo_centroid(const std::vector<base_node> &pts,
                                const std::vector<int> &idx) {
    base_node c(pts[0].size()); c.fill(0.);
    for (int i : idx) for (size_type k=0; k < c.size(); ++k) c[k] += pts[i][k];
    for (size_type k=0; k < c.size(); ++k) c[k] /= scalar_type(idx.size());
    return c;
  }

  /* face_perm[getfem local face] = Exodus side index, matched by comparing the
     face centroids on the reference element (same idea as the node mapping, so
     no hand-maintained face table). */
  static std::vector<unsigned> exo_face_perm(exo_topo topo) {
    bgeot::pgeometric_trans lpgt = exo_pgt(topo, exo_corner_count(topo));
    bgeot::pconvex_structure cvs = lpgt->structure();
    std::vector<base_node> gc(lpgt->geometric_nodes().begin(),
                              lpgt->geometric_nodes().end());
    std::vector<base_node> ec = exo_corners(topo);
    std::vector<std::vector<int> > sides = exo_sides(topo);
    short_type nf = cvs->nb_faces();
    GMM_ASSERT1(size_type(nf) == sides.size(),
                "Exodus side-count mismatch for this topology");
    std::vector<unsigned> fperm(nf);
    for (short_type f=0; f < nf; ++f) {
      const bgeot::convex_ind_ct &ip = cvs->ind_points_of_face(f);
      base_node cg = exo_centroid(gc, std::vector<int>(ip.begin(), ip.end()));
      size_type best = size_type(-1); scalar_type bd = 1e30;
      for (size_type s=0; s < sides.size(); ++s) {
        scalar_type d = gmm::vect_dist2(cg, exo_centroid(ec, sides[s]));
        if (d < bd) { bd = d; best = s; }
      }
      GMM_ASSERT1(best != size_type(-1) && bd < 1e-6,
                  "could not match GetFEM face to Exodus side");
      fperm[f] = unsigned(best);
    }
    return fperm;
  }


  /* ********************************************************************* */
  /*  Exodus export                                                        */
  /* ********************************************************************* */

  exodus_export::exodus_export(const std::string &fname, bool append)
    : fname_(fname), ncid_(-1), file_open_(false) {
    init();
    if (append) open_for_append_();
  }

  exodus_export::~exodus_export() { try { close(); } catch (...) {} }

  void exodus_export::init() {
    strcpy(title_, "Exported by GetFEM");
    nb_nodes_ = 0; dim_ = dim_type(-1); cur_step_ = -1; state_ = EMPTY;
    v_time_ = -1; v_nod_names_ = -1; v_elem_names_ = -1; v_elem_var_tab_ = -1;
    vars_declared_ = false; transient_definitions_ready_ = false;
    append_mode_ = false;
    file_nb_nodes_ = 0; region_field_ = false; append_base_step_ = -1;
    cur_time_ = scalar_type(0); cur_time_valid_ = false; cur_time_written_ = false;
    compression_explicit_ = false;
#ifdef GETFEM_HAVE_EXODUS_NETCDF4
    compression_ = true;  // compressed NetCDF4 by default when available
#else
    compression_ = false; // classic 64-bit offset when NetCDF4 is unavailable
#endif
    compression_level_ = 1;
    file_mesh_fingerprint_.clear();
  }

  /* Reopen an existing Exodus file to append further time steps: read the node
     count, the current number of steps and the already-declared transient
     variables, so set_time/write_point_data continue where the file left off. */
  void exodus_export::open_for_append_() {
    EXNC(nc_open(fname_.c_str(), NC_WRITE, &ncid_));
    file_open_ = true;
    append_mode_ = true;
    int d; size_t n = 0;
    EXNC(nc_inq_dimid(ncid_, "num_nodes", &d));
    EXNC(nc_inq_dimlen(ncid_, d, &n));
    file_nb_nodes_ = n; d_num_nodes_ = d;
    file_mesh_fingerprint_ =
      exo_get_att_text(ncid_, NC_GLOBAL, "getfem_mesh_fingerprint");
    EXNC(nc_inq_dimid(ncid_, "time_step", &d_time_));
    size_t nsteps = 0; nc_inq_dimlen(ncid_, d_time_, &nsteps);
    cur_step_ = int(nsteps) - 1;
    append_base_step_ = cur_step_;        // steps up to here belong to the file
    int dln;
    if (nc_inq_dimid(ncid_, "len_name", &dln) == NC_NOERR) d_len_name_ = dln;
    int d_nv;
    if (nc_inq_dimid(ncid_, "num_nod_var", &d_nv) == NC_NOERR) {
      size_t nv = 0; EXNC(nc_inq_dimlen(ncid_, d_nv, &nv));
      int v_names; EXNC(nc_inq_varid(ncid_, "name_nod_var", &v_names));
      EXNC(nc_inq_varid(ncid_, "time_whole", &v_time_));
      for (size_t k=0; k < nv; ++k) {
        std::vector<char> buf(EXO_STRLEN+1, '\0');
        size_t start[2] = { k, 0 }, count[2] = { 1, EXO_STRLEN };
        EXNC(nc_get_vara_text(ncid_, v_names, start, count, buf.data()));
        std::string nm(buf.data());
        var_names_.push_back(nm);
        std::ostringstream s; s << "vals_nod_var" << (k+1);
        int vid; EXNC(nc_inq_varid(ncid_, s.str().c_str(), &vid));
        var_id_[nm] = vid;
      }
      vars_declared_ = true;
      transient_definitions_ready_ = true;   // already defined in the file
    }
    // rediscover the "region" element variable so appended steps keep writing it
    int d_ev;
    if (nc_inq_dimid(ncid_, "num_elem_var", &d_ev) == NC_NOERR) {
      size_t nev = 0; EXNC(nc_inq_dimlen(ncid_, d_ev, &nev));
      GMM_ASSERT1(nev == 1, "Exodus append: GetFEM only supports appending "
                  "to files whose sole element variable is 'region'");
      int v_enames;
      EXNC(nc_inq_varid(ncid_, "name_elem_var", &v_enames));
      std::vector<char> buf(EXO_STRLEN+1, '\0');
      size_t start[2] = { 0, 0 }, count[2] = { 1, EXO_STRLEN };
      EXNC(nc_get_vara_text(ncid_, v_enames, start, count, buf.data()));
      GMM_ASSERT1(std::string(buf.data()) == "region",
                  "Exodus append: existing element variable is not 'region'");
      int d_blk; EXNC(nc_inq_dimid(ncid_, "num_el_blk", &d_blk));
      size_t nblk = 0; EXNC(nc_inq_dimlen(ncid_, d_blk, &nblk));
      elem_var_id_.resize(nblk);
      for (size_t b=0; b < nblk; ++b) {
        std::ostringstream sv; sv << "vals_elem_var1eb" << (b+1);
        EXNC(nc_inq_varid(ncid_, sv.str().c_str(), &elem_var_id_[b]));
      }
      // the file already holds (region) transient definitions: mark them ready
      // so appended steps stream instead of re-entering define mode (which would
      // collide with the existing variables) -- needed for region-only files
      region_field_ = true;
      vars_declared_ = true;
      transient_definitions_ready_ = true;
    }
  }

  void exodus_export::set_header(const std::string &s) {
    strncpy(title_, s.c_str(), 255); title_[255] = 0;
  }

  void exodus_export::set_region_name(int id, const std::string &name) {
    region_names_[id] = name;
  }

  void exodus_export::enable_region_field(bool on) { region_field_ = on; }

  void exodus_export::enable_compression(int level) {
    GMM_ASSERT1(state_ < MESH_WRITTEN && !append_mode_,
                "Exodus compression must be (de)selected before the file is "
                "created, and cannot be changed on append");
    if (level > 0) {
#ifndef GETFEM_HAVE_EXODUS_NETCDF4
      GMM_ASSERT1(false, "Exodus compression requires NetCDF4/HDF5 deflate "
                  "support; rebuild GetFEM with such a NetCDF library or pass "
                  "level 0 / 'uncompressed'");
#endif
    }
    compression_explicit_ = true;
    compression_ = (level > 0);
    compression_level_ = std::max(0, std::min(level, 9));
  }

  void exodus_export::declare_point_data(const std::string &name,
                                         size_type qdim) {
    GMM_ASSERT1(!append_mode_, "Exodus append uses the variables already "
                "declared in the existing file");
    GMM_ASSERT1(!vars_declared_ && state_ < MESH_WRITTEN,
                "declare_point_data() must be called before write_mesh()");
    append_point_data_names_(name, qdim);
  }

  const mesh_fem &exodus_export::get_exported_mesh_fem() const {
    GMM_ASSERT1(pmf_.get(), "no mesh_fem to export"); return *pmf_;
  }

  void exodus_export::exporting(const mesh &m) {
    dim_ = m.dim();
    GMM_ASSERT1(dim_ <= 3, "Exodus export: dimension " << int(dim_)
                << " not supported");
    pmf_ = std::make_unique<mesh_fem>(const_cast<mesh&>(m), dim_type(1));
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
      pfem pf = classical_fem(pgt, pgt->complexity() > 1 ? 2 : 1);
      pmf_->set_finite_element(cv, pf);
    }
    exporting(*pmf_);
  }

  void exodus_export::exporting(const mesh_fem &mf) {
    dim_ = mf.linked_mesh().dim();
    GMM_ASSERT1(dim_ <= 3, "Exodus export: dimension " << int(dim_)
                << " not supported");
    build_export_mesh_fem(mf);
    compute_blocks();
    if (append_mode_) {
      GMM_ASSERT1(nb_nodes_ == file_nb_nodes_, "Exodus append: the mesh_fem "
                  "yields " << nb_nodes_ << " export nodes but the existing file "
                  "has " << file_nb_nodes_ << " (must be the same mesh_fem)");
      GMM_ASSERT1(!file_mesh_fingerprint_.empty(),
                  "Exodus append: existing file has no GetFEM mesh "
                  "fingerprint; unsafe append is refused");
      std::string fp = mesh_fingerprint_();
      GMM_ASSERT1(fp == file_mesh_fingerprint_,
                  "Exodus append: mesh/block/connectivity fingerprint mismatch; "
                  "append requires the same exported mesh that created the file");
      state_ = MESH_WRITTEN;   // mesh is already in the file; do not rewrite it
    } else
      state_ = EXPORTING;
  }

  /* Build an export mesh_fem of classical (order 1 or 2) elements suitable for
     Exodus. This mirrors vtk_export::exporting(const mesh_fem&). */
  void exodus_export::build_export_mesh_fem(const mesh_fem &mf) {
    if (&mf != pmf_.get())
      pmf_ = std::make_unique<mesh_fem>(mf.linked_mesh());
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
      pfem pf = mf.fem_of_element(cv);
      if (pf == fem_descriptor("FEM_Q2_INCOMPLETE(2)") ||
          pf == fem_descriptor("FEM_Q2_INCOMPLETE(3)") ||
          pf == fem_descriptor("FEM_PRISM_INCOMPLETE_P2") ||
          pf == fem_descriptor("FEM_PRISM_INCOMPLETE_P2_DISCONTINUOUS"))
        pmf_->set_finite_element(cv, pf);
      else {
        bool discontinuous = false;
        for (unsigned i=0; i < pf->nb_dof(cv); ++i)
          if (!dof_linkable(pf->dof_types()[i])) { discontinuous = true; break; }
        pfem classical_pf1 = discontinuous ? classical_discontinuous_fem(pgt, 1)
                                           : classical_fem(pgt, 1);
        short_type degree = 1;
        if ((pf != classical_pf1 && pf->estimated_degree() > 1) ||
            pgt->structure() != pgt->basic_structure())
          degree = 2;
        pmf_->set_finite_element(cv, discontinuous ?
                                 classical_discontinuous_fem(pgt, degree, 0, true) :
                                 classical_fem(pgt, degree, true));
      }
    }
  }

  /* Determine, per convex, the Exodus block it belongs to (one block per
     element type), the GetFEM->Exodus node permutation, and the set of dofs
     actually written as Exodus nodes. */
  void exodus_export::compute_blocks() {
    const mesh &m = pmf_->linked_mesh();
    blocks_.clear();
    pmf_dof_used_.sup(0, pmf_->nb_basic_dof());

    /* A convex's "primary" volume region = the smallest-id region that contains
       it as a whole convex (face entries -> side sets, never blocks). When such
       regions exist each one becomes its own Exodus block, so tools like
       ParaView show/colour the regions natively (by block / ObjectId); convexes
       in no volume region fall in a default block keyed by element type. Element
       sets are still written too, so the GetFEM region round-trip is unchanged. */
    // indexed by convex id (contiguous): -1 = convex in no volume region
    std::vector<int> primary_region(m.convex_index().last_true() + 1, -1);
    for (dal::bv_visitor r(m.regions_index()); !r.finished(); ++r)
      for (getfem::mr_visitor i(m.region(r)); !i.finished(); ++i) {
        if (i.is_face()) continue;
        int &pr = primary_region[i.cv()];
        if (pr < 0 || int(r) < pr) pr = int(r);
      }

    std::map<std::pair<int, std::string>, size_type> block_of_key;
    for (dal::bv_visitor cv(pmf_->convex_index()); !cv.finished(); ++cv) {
      pfem pf = pmf_->fem_of_element(cv);
      size_type nbd = pf->nb_dof(cv);
      exo_topo topo = exo_topo_of_structure(m.structure_of_convex(cv));
      std::string et = exo_typestring(topo, nbd);
      int rid = primary_region[cv];          // -1 when in no volume region

      size_type bi;
      std::pair<int, std::string> key(rid, et);
      auto it = block_of_key.find(key);
      if (it == block_of_key.end()) {
        bi = blocks_.size();
        block_of_key[key] = bi;
        block_info b;
        b.ex_type = et;
        b.nodes_per_elem = nbd;
        b.region_id = rid;
        b.block_id = -1;              // assigned below
        // build the permutation once per block from this convex's reference dofs
        std::vector<base_node> gfref(nbd);
        for (size_type i=0; i < nbd; ++i) gfref[i] = pf->node_of_dof(cv, i);
        b.perm = exo_build_perm(gfref, exo_reference(topo, nbd));
        blocks_.push_back(b);
      } else bi = it->second;

      blocks_[bi].convexes.push_back(cv);
      for (size_type i=0; i < nbd; ++i)
        pmf_dof_used_.add(pmf_->ind_basic_dof_of_element(cv)[i]);
    }

    /* Exodus block ids: a region-based block claims its region id (so ParaView's
       ObjectId equals the GetFEM region id); any extra block -- a second element
       type within the same region, or a non-region block -- gets the next free
       positive id not already claimed by a region. */
    std::set<int> used;
    for (const auto &b : blocks_) if (b.region_id >= 0) used.insert(b.region_id);
    std::set<int> region_taken;
    int nextfree = 1;
    for (auto &b : blocks_) {
      if (b.region_id >= 0 && !region_taken.count(b.region_id)) {
        b.block_id = b.region_id;
        region_taken.insert(b.region_id);
      } else {
        while (used.count(nextfree)) ++nextfree;
        b.block_id = nextfree;
        used.insert(nextfree);
      }
    }

    // compact the used dofs to a 0-based contiguous node numbering
    dofmap_.assign(pmf_->nb_basic_dof(), -1);
    int cnt = 0;
    for (dal::bv_visitor d(pmf_dof_used_); !d.finished(); ++d) dofmap_[d] = cnt++;
    nb_nodes_ = size_type(cnt);
  }

  std::string exodus_export::mesh_fingerprint_() const {
    uint64_t h = 1469598103934665603ULL;
    fnv_add_string(h, "getfem-exodus-mesh-v1");
    fnv_add_int(h, int64_t(dim_));
    fnv_add_size(h, nb_nodes_);
    fnv_add_size(h, blocks_.size());
    for (dal::bv_visitor d(pmf_dof_used_); !d.finished(); ++d) {
      base_node P = pmf_->point_of_basic_dof(d);
      for (size_type k=0; k < size_type(dim_); ++k) fnv_add_scalar(h, P[k]);
    }
    for (size_type b=0; b < blocks_.size(); ++b) {
      const block_info &bk = blocks_[b];
      fnv_add_int(h, bk.block_id);
      fnv_add_int(h, bk.region_id);
      fnv_add_string(h, bk.ex_type);
      fnv_add_size(h, bk.nodes_per_elem);
      fnv_add_size(h, bk.convexes.size());
      for (size_type cv : bk.convexes) {
        fnv_add_size(h, cv);
        auto ind = pmf_->ind_basic_dof_of_element(cv);
        for (size_type e=0; e < bk.nodes_per_elem; ++e)
          fnv_add_int(h, dofmap_[ind[bk.perm[e]]] + 1);
      }
    }
    return std::string("getfem-exodus-mesh-v1:") + fnv_hex(h);
  }

  void exodus_export::deflate_var_(int varid) const {
#ifndef GETFEM_HAVE_EXODUS_NETCDF4
    (void)varid;
#endif
    if (!compression_) return;
#ifdef GETFEM_HAVE_EXODUS_NETCDF4
    EXNC(nc_def_var_deflate(ncid_, varid, 1, 1, compression_level_));
#else
    GMM_ASSERT1(false, "Exodus compression requires NetCDF4/HDF5 deflate "
                "support");
#endif
  }

  /* write a fixed-length, NUL-padded string at row `row` of a (n, len) char var */
  static void exo_put_string(int ncid, int varid, size_t row, size_t len,
                             const std::string &s) {
    std::vector<char> buf(len, '\0');
    for (size_t i=0; i < len-1 && i < s.size(); ++i) buf[i] = s[i];
    size_t start[2] = { row, 0 }, count[2] = { 1, len };
    EXNC(nc_put_vara_text(ncid, varid, start, count, buf.data()));
  }

  void exodus_export::write_mesh() {
    GMM_ASSERT1(state_ >= EXPORTING, "call exporting() first");
    if (state_ >= MESH_WRITTEN) return;

    size_type nb_elem = 0;
    for (const auto &b : blocks_) nb_elem += b.convexes.size();
    GMM_ASSERT1(nb_nodes_ > 0 && nb_elem > 0, "nothing to export");

    /* Assign 1-based global Exodus element ids in block order and collect side
       sets from the GetFEM face regions (a region's face entries -> a side
       set; whole-convex region entries have no Exodus equivalent and are
       skipped). */
    const mesh &m = pmf_->linked_mesh();
    // indexed by convex id (contiguous): 0 = convex not exported as an element
    std::vector<int> cv_to_elem(m.convex_index().last_true() + 1, 0);
    { int gid = 0;
      for (const auto &b : blocks_)
        for (size_type cv : b.convexes) cv_to_elem[cv] = ++gid; }

    /* Each GetFEM region maps (by its id) to: a side set (its face entries),
       an element set (its whole-convex entries), and a node set (all nodes it
       touches). Side/element sets round-trip back to face/convex regions;
       node sets are extra (useful for nodal BCs) and read back via node_set(). */
    std::map<exo_topo, std::vector<unsigned> > fperm_cache;
    std::vector<int> ss_ids, es_ids, ns_ids;
    std::vector<std::vector<int> > ss_elem, ss_side, es_elem, ns_node;
    for (dal::bv_visitor r(m.regions_index()); !r.finished(); ++r) {
      std::vector<int> sel, ssd, eel;
      std::set<int> nds;
      for (getfem::mr_visitor i(m.region(r)); !i.finished(); ++i) {
        int eid = cv_to_elem[i.cv()];
        if (eid == 0) continue;             // convex not exported as an element
        if (i.is_face()) {
          exo_topo topo = exo_topo_of_structure(m.structure_of_convex(i.cv()));
          auto fit = fperm_cache.find(topo);
          if (fit == fperm_cache.end())
            fit = fperm_cache.emplace(topo, exo_face_perm(topo)).first;
          sel.push_back(eid);
          ssd.push_back(int(fit->second[i.f()]) + 1);
          auto fd = pmf_->ind_basic_dof_of_face_of_element(i.cv(), i.f());
          for (size_type j=0; j < fd.size(); ++j) nds.insert(dofmap_[fd[j]] + 1);
        } else {
          eel.push_back(eid);
          auto cd = pmf_->ind_basic_dof_of_element(i.cv());
          for (size_type j=0; j < cd.size(); ++j) nds.insert(dofmap_[cd[j]] + 1);
        }
      }
      if (!sel.empty()) {
        ss_ids.push_back(int(r)); ss_elem.push_back(sel); ss_side.push_back(ssd); }
      if (!eel.empty()) { es_ids.push_back(int(r)); es_elem.push_back(eel); }
      if (!nds.empty()) {
        ns_ids.push_back(int(r));
        ns_node.push_back(std::vector<int>(nds.begin(), nds.end())); }
    }
    size_type n_ss = ss_ids.size(), n_es = es_ids.size(), n_ns = ns_ids.size();
    std::vector<int> v_ss_elem(n_ss, -1), v_ss_side(n_ss, -1);
    std::vector<int> v_es_elem(n_es, -1), v_ns_node(n_ns, -1);
    int v_ssp = -1, v_sss = -1, v_esp = -1, v_ess = -1, v_nsp = -1, v_nss = -1;
    int v_ssn = -1, v_esn = -1, v_nsn = -1;   // side/elem/node-set name char vars

    int cmode = NC_CLOBBER | NC_64BIT_OFFSET;
    int e_create = NC_NOERR;
    if (compression_) {
#ifdef GETFEM_HAVE_EXODUS_NETCDF4
      cmode = NC_CLOBBER | NC_NETCDF4 | NC_CLASSIC_MODEL;
      e_create = nc_create(fname_.c_str(), cmode, &ncid_);
      if (e_create != NC_NOERR && !compression_explicit_) {
        compression_ = false;
        cmode = NC_CLOBBER | NC_64BIT_OFFSET;
        e_create = nc_create(fname_.c_str(), cmode, &ncid_);
      }
#else
      compression_ = false;
      e_create = nc_create(fname_.c_str(), cmode, &ncid_);
#endif
    } else {
      e_create = nc_create(fname_.c_str(), cmode, &ncid_);
    }
    GMM_ASSERT1(e_create == NC_NOERR, "NetCDF/Exodus error: "
                << nc_strerror(e_create));
    file_open_ = true;

    /* --- dimensions --- */
    int d_dim, d_nodes, d_elem, d_blk, d_str, d_name, d_four, d_time, d_qa;
    EXNC(nc_def_dim(ncid_, "num_dim",    size_t(dim_),  &d_dim));
    EXNC(nc_def_dim(ncid_, "num_nodes",  nb_nodes_,     &d_nodes));
    EXNC(nc_def_dim(ncid_, "num_elem",   nb_elem,       &d_elem));
    EXNC(nc_def_dim(ncid_, "num_el_blk", blocks_.size(),&d_blk));
    EXNC(nc_def_dim(ncid_, "len_string", EXO_STRLEN,    &d_str));
    EXNC(nc_def_dim(ncid_, "len_name",   EXO_STRLEN,    &d_name));
    EXNC(nc_def_dim(ncid_, "four",       4,             &d_four));
    EXNC(nc_def_dim(ncid_, "num_qa_rec", 1,             &d_qa));
    EXNC(nc_def_dim(ncid_, "time_step",  NC_UNLIMITED,  &d_time));
    d_num_nodes_ = d_nodes; d_time_ = d_time; d_len_name_ = d_name;

    /* --- global attributes --- */
    float api = 8.03f, ver = 8.03f;
    int   fpws = int(sizeof(scalar_type)), fsz = 1, i64 = 0;
    EXNC(nc_put_att_float(ncid_, NC_GLOBAL, "api_version", NC_FLOAT, 1, &api));
    EXNC(nc_put_att_float(ncid_, NC_GLOBAL, "version",     NC_FLOAT, 1, &ver));
    EXNC(nc_put_att_int(ncid_, NC_GLOBAL, "floating_point_word_size",
                        NC_INT, 1, &fpws));
    EXNC(nc_put_att_int(ncid_, NC_GLOBAL, "file_size",    NC_INT, 1, &fsz));
    EXNC(nc_put_att_int(ncid_, NC_GLOBAL, "int64_status", NC_INT, 1, &i64));
    EXNC(nc_put_att_text(ncid_, NC_GLOBAL, "title",
                         strlen(title_), title_));
    { std::string fp = mesh_fingerprint_();
      EXNC(nc_put_att_text(ncid_, NC_GLOBAL, "getfem_mesh_fingerprint",
                           fp.size(), fp.c_str())); }
    { const char writer[] = "GetFEM direct NetCDF Exodus writer";
      EXNC(nc_put_att_text(ncid_, NC_GLOBAL, "getfem_exodus_writer",
                           strlen(writer), writer)); }

    /* --- coordinate, block and connectivity variables --- */
    int v_x=-1, v_y=-1, v_z=-1, v_cn, v_ebp, v_ebs, v_ebn, v_qa;
    EXNC(nc_def_var(ncid_, "coordx", NC_DOUBLE, 1, &d_nodes, &v_x));
    deflate_var_(v_x);
    if (dim_ >= 2) { EXNC(nc_def_var(ncid_, "coordy", NC_DOUBLE, 1, &d_nodes, &v_y));
      deflate_var_(v_y); }
    if (dim_ >= 3) { EXNC(nc_def_var(ncid_, "coordz", NC_DOUBLE, 1, &d_nodes, &v_z));
      deflate_var_(v_z); }
    { int dd[2] = { d_dim, d_str };
      EXNC(nc_def_var(ncid_, "coor_names", NC_CHAR, 2, dd, &v_cn)); }

    EXNC(nc_def_var(ncid_, "eb_prop1", NC_INT, 1, &d_blk, &v_ebp));
    deflate_var_(v_ebp);
    EXNC(nc_put_att_text(ncid_, v_ebp, "name", 2, "ID"));
    EXNC(nc_def_var(ncid_, "eb_status", NC_INT, 1, &d_blk, &v_ebs));
    deflate_var_(v_ebs);
    { int dd[2] = { d_blk, d_name };
      EXNC(nc_def_var(ncid_, "eb_names", NC_CHAR, 2, dd, &v_ebn)); }
    { int dd[3] = { d_qa, d_four, d_str };
      EXNC(nc_def_var(ncid_, "qa_records", NC_CHAR, 3, dd, &v_qa)); }

    std::vector<int> v_connect(blocks_.size());
    for (size_type b=0; b < blocks_.size(); ++b) {
      std::ostringstream sd1, sd2;
      sd1 << "num_el_in_blk" << (b+1); sd2 << "num_nod_per_el" << (b+1);
      int de, dn;
      EXNC(nc_def_dim(ncid_, sd1.str().c_str(), blocks_[b].convexes.size(), &de));
      EXNC(nc_def_dim(ncid_, sd2.str().c_str(), blocks_[b].nodes_per_elem, &dn));
      std::ostringstream sc; sc << "connect" << (b+1);
      int dd[2] = { de, dn };
      EXNC(nc_def_var(ncid_, sc.str().c_str(), NC_INT, 2, dd, &v_connect[b]));
      deflate_var_(v_connect[b]);
      EXNC(nc_put_att_text(ncid_, v_connect[b], "elem_type",
                           blocks_[b].ex_type.size(), blocks_[b].ex_type.c_str()));
    }

    // identity 1-based id maps (ParaView GlobalElementId / GlobalNodeId)
    int v_emap, v_nmap;
    EXNC(nc_def_var(ncid_, "elem_num_map", NC_INT, 1, &d_elem,  &v_emap));
    deflate_var_(v_emap);
    EXNC(nc_def_var(ncid_, "node_num_map", NC_INT, 1, &d_nodes, &v_nmap));
    deflate_var_(v_nmap);

    if (n_ss > 0) {
      int d_nss;
      EXNC(nc_def_dim(ncid_, "num_side_sets", n_ss, &d_nss));
      EXNC(nc_def_var(ncid_, "ss_prop1", NC_INT, 1, &d_nss, &v_ssp));
      deflate_var_(v_ssp);
      EXNC(nc_put_att_text(ncid_, v_ssp, "name", 2, "ID"));
      EXNC(nc_def_var(ncid_, "ss_status", NC_INT, 1, &d_nss, &v_sss));
      deflate_var_(v_sss);
      { int dd[2] = { d_nss, d_name };
        EXNC(nc_def_var(ncid_, "ss_names", NC_CHAR, 2, dd, &v_ssn)); }
      for (size_type i=0; i < n_ss; ++i) {
        std::ostringstream sd; sd << "num_side_ss" << (i+1);
        int dss;
        EXNC(nc_def_dim(ncid_, sd.str().c_str(), ss_elem[i].size(), &dss));
        std::ostringstream se, sf;
        se << "elem_ss" << (i+1); sf << "side_ss" << (i+1);
        EXNC(nc_def_var(ncid_, se.str().c_str(), NC_INT, 1, &dss, &v_ss_elem[i]));
        deflate_var_(v_ss_elem[i]);
        EXNC(nc_def_var(ncid_, sf.str().c_str(), NC_INT, 1, &dss, &v_ss_side[i]));
        deflate_var_(v_ss_side[i]);
      }
    }

    if (n_es > 0) {
      int d_nes;
      EXNC(nc_def_dim(ncid_, "num_elem_sets", n_es, &d_nes));
      EXNC(nc_def_var(ncid_, "els_prop1", NC_INT, 1, &d_nes, &v_esp));
      deflate_var_(v_esp);
      EXNC(nc_put_att_text(ncid_, v_esp, "name", 2, "ID"));
      EXNC(nc_def_var(ncid_, "els_status", NC_INT, 1, &d_nes, &v_ess));
      deflate_var_(v_ess);
      { int dd[2] = { d_nes, d_name };
        EXNC(nc_def_var(ncid_, "els_names", NC_CHAR, 2, dd, &v_esn)); }
      for (size_type i=0; i < n_es; ++i) {
        std::ostringstream sd; sd << "num_ele_els" << (i+1);
        int des;
        EXNC(nc_def_dim(ncid_, sd.str().c_str(), es_elem[i].size(), &des));
        std::ostringstream se; se << "elem_els" << (i+1);
        EXNC(nc_def_var(ncid_, se.str().c_str(), NC_INT, 1, &des, &v_es_elem[i]));
        deflate_var_(v_es_elem[i]);
      }
    }

    if (n_ns > 0) {
      int d_nns;
      EXNC(nc_def_dim(ncid_, "num_node_sets", n_ns, &d_nns));
      EXNC(nc_def_var(ncid_, "ns_prop1", NC_INT, 1, &d_nns, &v_nsp));
      deflate_var_(v_nsp);
      EXNC(nc_put_att_text(ncid_, v_nsp, "name", 2, "ID"));
      EXNC(nc_def_var(ncid_, "ns_status", NC_INT, 1, &d_nns, &v_nss));
      deflate_var_(v_nss);
      { int dd[2] = { d_nns, d_name };
        EXNC(nc_def_var(ncid_, "ns_names", NC_CHAR, 2, dd, &v_nsn)); }
      for (size_type i=0; i < n_ns; ++i) {
        std::ostringstream sd; sd << "num_nod_ns" << (i+1);
        int dns;
        EXNC(nc_def_dim(ncid_, sd.str().c_str(), ns_node[i].size(), &dns));
        std::ostringstream sn; sn << "node_ns" << (i+1);
        EXNC(nc_def_var(ncid_, sn.str().c_str(), NC_INT, 1, &dns, &v_ns_node[i]));
        deflate_var_(v_ns_node[i]);
      }
    }
    bool predeclared_transient = !var_names_.empty();
    if (predeclared_transient) define_transient_vars_();
    EXNC(nc_enddef(ncid_));
    if (predeclared_transient) {
      write_transient_metadata_();
      vars_declared_ = true;
    }

    /* --- coordinate data (gathered over the compacted nodes) --- */
    std::vector<double> X(nb_nodes_,0.), Y(nb_nodes_,0.), Z(nb_nodes_,0.);
    for (dal::bv_visitor d(pmf_dof_used_); !d.finished(); ++d) {
      base_node P = pmf_->point_of_basic_dof(d);
      int i = dofmap_[d];
      X[i] = P[0];
      if (dim_ >= 2) Y[i] = P[1];
      if (dim_ >= 3) Z[i] = P[2];
    }
    EXNC(nc_put_var_double(ncid_, v_x, X.data()));
    if (dim_ >= 2) EXNC(nc_put_var_double(ncid_, v_y, Y.data()));
    if (dim_ >= 3) EXNC(nc_put_var_double(ncid_, v_z, Z.data()));
    { const char *cn[3] = { "x", "y", "z" };
      for (size_type k=0; k < size_type(dim_); ++k)
        exo_put_string(ncid_, v_cn, k, EXO_STRLEN, cn[k]); }

    /* --- block ids / status / names --- */
    std::vector<int> ebp(blocks_.size()), ebs(blocks_.size(), 1);
    for (size_type b=0; b < blocks_.size(); ++b) ebp[b] = blocks_[b].block_id;
    EXNC(nc_put_var_int(ncid_, v_ebp, ebp.data()));
    EXNC(nc_put_var_int(ncid_, v_ebs, ebs.data()));
    for (size_type b=0; b < blocks_.size(); ++b) {
      // a named region uses its name; an unnamed region block -> "region_<id>";
      // a non-region block keeps the element-type name.
      std::string bname = blocks_[b].ex_type;
      if (blocks_[b].region_id >= 0) {
        auto it = region_names_.find(blocks_[b].region_id);
        if (it != region_names_.end() && !it->second.empty()) bname = it->second;
        else { std::ostringstream nm; nm << "region_" << blocks_[b].region_id;
               bname = nm.str(); }
      }
      exo_put_string(ncid_, v_ebn, b, EXO_STRLEN, bname);
    }

    /* --- one QA record --- */
    char date[16] = "00000000", tim[16] = "000000";
    std::time_t tt = std::time(nullptr);
    if (tt != std::time_t(-1)) {
      std::tm *lt = std::localtime(&tt);
      if (lt) { std::strftime(date, sizeof date, "%Y/%m/%d", lt);
                std::strftime(tim, sizeof tim, "%H:%M:%S", lt); }
    }
    { const char *qa[4] = { "GetFEM", "exodus", date, tim };
      size_t st[3], ct[3] = { 1, 1, EXO_STRLEN };
      for (int j=0; j < 4; ++j) {
        std::vector<char> buf(EXO_STRLEN, '\0');
        for (size_t i=0; i+1 < EXO_STRLEN && i < strlen(qa[j]); ++i)
          buf[i] = qa[j][i];
        st[0] = 0; st[1] = size_t(j); st[2] = 0;
        EXNC(nc_put_vara_text(ncid_, v_qa, st, ct, buf.data()));
      } }

    /* --- connectivity (1-based, Exodus node ordering) --- */
    for (size_type b=0; b < blocks_.size(); ++b) {
      const block_info &bk = blocks_[b];
      size_type chunk_elems =
        std::max<size_type>(1, EXO_IO_CHUNK / std::max<size_type>(1, bk.nodes_per_elem));
      std::vector<int> conn(chunk_elems * bk.nodes_per_elem);
      for (size_type first=0; first < bk.convexes.size(); first += chunk_elems) {
        size_type ne = std::min(chunk_elems, bk.convexes.size() - first);
        size_type p = 0;
        for (size_type j=0; j < ne; ++j) {
          size_type cv = bk.convexes[first + j];
          auto ind = pmf_->ind_basic_dof_of_element(cv);
          for (size_type e=0; e < bk.nodes_per_elem; ++e)
            conn[p++] = dofmap_[ind[bk.perm[e]]] + 1;
        }
        size_t start[2] = { first, 0 }, count[2] = { ne, bk.nodes_per_elem };
        EXNC(nc_put_vara_int(ncid_, v_connect[b], start, count, conn.data()));
      }
    }

    /* --- identity element/node number maps --- */
    { std::vector<int> emap(nb_elem);
      for (size_type i=0; i < nb_elem; ++i) emap[i] = int(i) + 1;
      EXNC(nc_put_var_int(ncid_, v_emap, emap.data())); }
    { std::vector<int> nmap(nb_nodes_);
      for (size_type i=0; i < nb_nodes_; ++i) nmap[i] = int(i) + 1;
      EXNC(nc_put_var_int(ncid_, v_nmap, nmap.data())); }

    /* --- side sets --- */
    if (n_ss > 0) {
      EXNC(nc_put_var_int(ncid_, v_ssp, ss_ids.data()));
      std::vector<int> sstat(n_ss, 1);
      EXNC(nc_put_var_int(ncid_, v_sss, sstat.data()));
      for (size_type i=0; i < n_ss; ++i) {
        EXNC(nc_put_var_int(ncid_, v_ss_elem[i], ss_elem[i].data()));
        EXNC(nc_put_var_int(ncid_, v_ss_side[i], ss_side[i].data()));
      }
      for (size_type i=0; i < n_ss; ++i) {  // optional set names (empty if unset)
        auto it = region_names_.find(ss_ids[i]);
        exo_put_string(ncid_, v_ssn, i, EXO_STRLEN,
                       it == region_names_.end() ? std::string() : it->second);
      }
    }

    /* --- element sets (from whole-convex region entries) --- */
    if (n_es > 0) {
      EXNC(nc_put_var_int(ncid_, v_esp, es_ids.data()));
      std::vector<int> estat(n_es, 1);
      EXNC(nc_put_var_int(ncid_, v_ess, estat.data()));
      for (size_type i=0; i < n_es; ++i)
        EXNC(nc_put_var_int(ncid_, v_es_elem[i], es_elem[i].data()));
      for (size_type i=0; i < n_es; ++i) {  // optional set names (empty if unset)
        auto it = region_names_.find(es_ids[i]);
        exo_put_string(ncid_, v_esn, i, EXO_STRLEN,
                       it == region_names_.end() ? std::string() : it->second);
      }
    }

    /* --- node sets (nodes touched by each region) --- */
    if (n_ns > 0) {
      EXNC(nc_put_var_int(ncid_, v_nsp, ns_ids.data()));
      std::vector<int> nstat(n_ns, 1);
      EXNC(nc_put_var_int(ncid_, v_nss, nstat.data()));
      for (size_type i=0; i < n_ns; ++i)
        EXNC(nc_put_var_int(ncid_, v_ns_node[i], ns_node[i].data()));
      for (size_type i=0; i < n_ns; ++i) {  // optional set names (empty if unset)
        auto it = region_names_.find(ns_ids[i]);
        exo_put_string(ncid_, v_nsn, i, EXO_STRLEN,
                       it == region_names_.end() ? std::string() : it->second);
      }
    }

    state_ = MESH_WRITTEN;
  }

  static std::string exo_comp_name(const std::string &name, size_type q,
                                   size_type Q) {
    if (Q == 1) return name;
    static const char *xyz[3] = { "x", "y", "z" };
    std::ostringstream s; s << name << "_";
    if (Q <= 3) s << xyz[q]; else s << q;
    return s.str();
  }

  void exodus_export::append_point_data_names_(const std::string &name,
                                               size_type Q) {
    GMM_ASSERT1(Q >= 1, "empty data");
    for (size_type q=0; q < Q; ++q) {
      std::string nm = remove_spaces(exo_comp_name(name, q, Q));
      if (nm.size() >= EXO_STRLEN) nm.resize(EXO_STRLEN - 1);
      bool seen = false;
      for (const auto &n : var_names_) if (n == nm) { seen = true; break; }
      GMM_ASSERT1(!seen, "Exodus transient export: duplicate variable name '"
                  << nm << "'");
      var_names_.push_back(nm);
    }
  }

  void exodus_export::set_time(scalar_type t) {
    if (state_ < MESH_WRITTEN) write_mesh();
    finish_current_step_(true);
    step_written_.clear();
    ++cur_step_;                          // continues from the file's last step
    cur_time_ = t;
    cur_time_valid_ = true;
    cur_time_written_ = false;
  }

  void exodus_export::add_nodal_data_(const std::vector<scalar_type> &V,
                                      size_type Q, const std::string &name) {
    if (cur_step_ < 0) set_time(0.);
    size_type step = size_type(cur_step_);
    for (size_type q=0; q < Q; ++q) {
      std::string nm = remove_spaces(exo_comp_name(name, q, Q));
      // Exodus names are fixed-length; truncate here so the in-memory name
      // matches the (truncated) name written to / read back from the file --
      // otherwise a long name would be rejected on append.
      if (nm.size() >= EXO_STRLEN) nm.resize(EXO_STRLEN - 1);
      if (!vars_declared_) {
        // step 0: record the variable (preserving order) and gather its values
        // straight into the step-0 buffer it will be flushed from
        bool seen = false;
        for (const auto &n : var_names_) if (n == nm) { seen = true; break; }
        if (!seen) var_names_.push_back(nm);
        std::vector<scalar_type> &buf = step0_[nm];
        buf.resize(nb_nodes_);
        for (size_type i=0; i < nb_nodes_; ++i) buf[i] = V[i*Q+q];
      } else {
        auto it = var_id_.find(nm);
        GMM_ASSERT1(it != var_id_.end(), "Exodus transient export: variable '"
                    << nm << "' appeared after the first time step; the set of "
                    "written fields must be the same at every step");
        col_.resize(nb_nodes_);              // reused across components/steps
        for (size_type i=0; i < nb_nodes_; ++i) col_[i] = V[i*Q+q];
        size_t start[2] = { step, 0 }, count[2] = { 1, nb_nodes_ };
        EXNC(nc_put_vara_double(ncid_, it->second, start, count, col_.data()));
        bool already = false;            // record that this field got this step
        for (const auto &n : step_written_) if (n == nm) { already = true; break; }
        if (!already) step_written_.push_back(nm);
      }
    }
  }

  void exodus_export::define_transient_vars_() {
    if (transient_definitions_ready_) return;
    GMM_ASSERT1(!var_names_.empty() || region_field_,
                "Exodus transient export: no variables to declare");
    EXNC(nc_def_var(ncid_, "time_whole", NC_DOUBLE, 1, &d_time_, &v_time_));
    deflate_var_(v_time_);
    if (!var_names_.empty()) {
      int d_nv;
      EXNC(nc_def_dim(ncid_, "num_nod_var", var_names_.size(), &d_nv));
      { int dd[2] = { d_nv, d_len_name_ };
        EXNC(nc_def_var(ncid_, "name_nod_var", NC_CHAR, 2, dd, &v_nod_names_)); }
      for (size_type k=0; k < var_names_.size(); ++k) {
        std::ostringstream s; s << "vals_nod_var" << (k+1);
        int dd[2] = { d_time_, d_num_nodes_ }, vid;
        EXNC(nc_def_var(ncid_, s.str().c_str(), NC_DOUBLE, 2, dd, &vid));
        deflate_var_(vid);
        var_id_[var_names_[k]] = vid;
      }
    }

    if (region_field_) {
      int d_nev, d_blk;
      EXNC(nc_def_dim(ncid_, "num_elem_var", 1, &d_nev));
      { int dd[2] = { d_nev, d_len_name_ };
        EXNC(nc_def_var(ncid_, "name_elem_var", NC_CHAR, 2, dd, &v_elem_names_)); }
      EXNC(nc_inq_dimid(ncid_, "num_el_blk", &d_blk));
      { int dd[2] = { d_blk, d_nev };
        EXNC(nc_def_var(ncid_, "elem_var_tab", NC_INT, 2, dd, &v_elem_var_tab_));
        deflate_var_(v_elem_var_tab_); }
      elem_var_id_.resize(blocks_.size());
      for (size_type b=0; b < blocks_.size(); ++b) {
        std::ostringstream sd; sd << "num_el_in_blk" << (b+1);
        int de; EXNC(nc_inq_dimid(ncid_, sd.str().c_str(), &de));
        std::ostringstream sv; sv << "vals_elem_var1eb" << (b+1);
        int dd[2] = { d_time_, de }, vid;
        EXNC(nc_def_var(ncid_, sv.str().c_str(), NC_DOUBLE, 2, dd, &vid));
        deflate_var_(vid);
        elem_var_id_[b] = vid;
      }
    }
    transient_definitions_ready_ = true;
  }

  void exodus_export::write_transient_metadata_() {
    if (v_nod_names_ >= 0)
      for (size_type k=0; k < var_names_.size(); ++k)
        exo_put_string(ncid_, v_nod_names_, k, EXO_STRLEN, var_names_[k]);
    if (v_elem_names_ >= 0) {
      exo_put_string(ncid_, v_elem_names_, 0, EXO_STRLEN, "region");
      std::vector<int> tab(blocks_.size(), 1);   // defined on every block
      EXNC(nc_put_var_int(ncid_, v_elem_var_tab_, tab.data()));
    }
  }

  /* Declare the transient NetCDF variables (once step 0's field set is known)
     and flush the buffered step 0. From then on steps stream directly. */
  void exodus_export::declare_transient_vars_() {
    EXNC(nc_redef(ncid_));
    define_transient_vars_();
    EXNC(nc_enddef(ncid_));
    write_transient_metadata_();
    vars_declared_ = true;

    for (const auto &kv : step0_) {
      size_t start[2] = { size_t(cur_step_), 0 }, count[2] = { 1, nb_nodes_ };
      EXNC(nc_put_vara_double(ncid_, var_id_[kv.first], start, count,
                              kv.second.data()));
      bool already = false;
      for (const auto &n : step_written_) if (n == kv.first) already = true;
      if (!already) step_written_.push_back(kv.first);
    }
    step0_.clear();
    write_step_time_and_region_();
    EXNC(nc_sync(ncid_));
  }

  void exodus_export::check_current_step_complete_() const {
    if (cur_step_ < 0 || !vars_declared_) return;
    if (append_mode_ && cur_step_ <= append_base_step_) return;
    GMM_ASSERT1(step_written_.size() == var_names_.size(),
                "Exodus transient export: time step " << cur_step_ << " wrote "
                << step_written_.size() << " of " << var_names_.size()
                << " field(s); every field must be written at every step");
  }

  void exodus_export::write_step_time_and_region_() {
    if (cur_step_ < 0 || cur_time_written_) return;
    GMM_ASSERT1(vars_declared_ && cur_time_valid_,
                "Exodus transient export: no current time step");
    size_t st = size_t(cur_step_), ct = 1;
    EXNC(nc_put_vara_double(ncid_, v_time_, &st, &ct, &cur_time_));
    if (!elem_var_id_.empty()) write_region_slab_(size_type(cur_step_));
    cur_time_written_ = true;
  }

  void exodus_export::finish_current_step_(bool do_sync) {
    if (cur_step_ < 0 || !cur_time_valid_) return;
    if (!vars_declared_) declare_transient_vars_();
    else {
      check_current_step_complete_();
      write_step_time_and_region_();
      if (do_sync) EXNC(nc_sync(ncid_));
    }
  }

  /* Write the (time-constant) "region" element variable for one step: every
     element of block b gets that block's id. Cheap: O(num_elem) per step. */
  void exodus_export::write_region_slab_(size_type step) {
    for (size_type b=0; b < blocks_.size(); ++b) {
      col_.assign(blocks_[b].convexes.size(), scalar_type(blocks_[b].block_id));
      size_t start[2] = { step, 0 }, count[2] = { 1, blocks_[b].convexes.size() };
      EXNC(nc_put_vara_double(ncid_, elem_var_id_[b], start, count, col_.data()));
    }
  }

  void exodus_export::sync() {
    if (!file_open_) return;
    if (state_ < MESH_WRITTEN) write_mesh();
    finish_current_step_(true);
  }

  void exodus_export::close() {
    if (!file_open_) return;
    try {
      if (region_field_ && cur_step_ < 0) set_time(0.);
      finish_current_step_(false);
    } catch (...) {
      nc_close(ncid_);
      file_open_ = false;
      throw;
    }
    EXNC(nc_close(ncid_));
    file_open_ = false;
  }


  /* ********************************************************************* */
  /*  Exodus import                                                        */
  /* ********************************************************************* */

  exodus_import::exodus_import(const std::string &fname)
    : ncid_(-1), file_open_(false) {
    int comp = 0;
    EXNC(nc_open(fname.c_str(), NC_NOWRITE, &ncid_));
    file_open_ = true; (void)comp;

    /* transient variable names */
    int d_nv;
    if (nc_inq_dimid(ncid_, "num_nod_var", &d_nv) == NC_NOERR) {
      size_t nv; EXNC(nc_inq_dimlen(ncid_, d_nv, &nv));
      int v_names, d_name; size_t lname = EXO_STRLEN;
      EXNC(nc_inq_varid(ncid_, "name_nod_var", &v_names));
      if (nc_inq_dimid(ncid_, "len_name", &d_name) == NC_NOERR)
        EXNC(nc_inq_dimlen(ncid_, d_name, &lname));
      for (size_t k=0; k < nv; ++k) {
        std::vector<char> buf(lname+1, '\0');
        size_t start[2] = { k, 0 }, count[2] = { 1, lname };
        EXNC(nc_get_vara_text(ncid_, v_names, start, count, buf.data()));
        nodal_var_names_.push_back(std::string(buf.data()));
      }
    }
    /* time values */
    int d_t;
    if (nc_inq_dimid(ncid_, "time_step", &d_t) == NC_NOERR) {
      size_t nt = 0; nc_inq_dimlen(ncid_, d_t, &nt);
      int v_t;
      if (nt > 0 && nc_inq_varid(ncid_, "time_whole", &v_t) == NC_NOERR) {
        times_.resize(nt);
        EXNC(nc_get_var_double(ncid_, v_t, times_.data()));
      }
    }
  }

  exodus_import::~exodus_import() { try { close(); } catch (...) {} }

  void exodus_import::close() {
    if (file_open_) { nc_close(ncid_); file_open_ = false; }
  }

  void exodus_import::read_mesh(mesh &m) {
    GMM_ASSERT1(file_open_, "Exodus file is closed");
    int d_dim, d_nodes, d_blk;
    size_t num_dim, num_nodes, num_blk;
    EXNC(nc_inq_dimid(ncid_, "num_dim", &d_dim));
    EXNC(nc_inq_dimlen(ncid_, d_dim, &num_dim));
    GMM_ASSERT1(num_dim >= 1 && num_dim <= 3,
                "Exodus import: unsupported spatial dimension " << num_dim
                << " (1, 2 or 3 expected)");
    EXNC(nc_inq_dimid(ncid_, "num_nodes", &d_nodes));
    EXNC(nc_inq_dimlen(ncid_, d_nodes, &num_nodes));
    EXNC(nc_inq_dimid(ncid_, "num_el_blk", &d_blk));
    EXNC(nc_inq_dimlen(ncid_, d_blk, &num_blk));

    /* coordinates */
    std::vector<double> X(num_nodes,0.), Y(num_nodes,0.), Z(num_nodes,0.);
    int v;
    if (nc_inq_varid(ncid_, "coordx", &v) == NC_NOERR) {       // split layout (modern)
      EXNC(nc_get_var_double(ncid_, v, X.data()));
      if (num_dim >= 2 && nc_inq_varid(ncid_, "coordy", &v) == NC_NOERR)
        EXNC(nc_get_var_double(ncid_, v, Y.data()));
      if (num_dim >= 3 && nc_inq_varid(ncid_, "coordz", &v) == NC_NOERR)
        EXNC(nc_get_var_double(ncid_, v, Z.data()));
    } else if (nc_inq_varid(ncid_, "coord", &v) == NC_NOERR) { // packed (legacy)
      std::vector<double> C(num_dim * num_nodes);              // coord(num_dim,num_nodes)
      EXNC(nc_get_var_double(ncid_, v, C.data()));
      for (size_t i=0; i < num_nodes; ++i) {
        X[i] = C[i];
        if (num_dim >= 2) Y[i] = C[num_nodes + i];
        if (num_dim >= 3) Z[i] = C[2*num_nodes + i];
      }
    } else GMM_ASSERT1(false, "Exodus import: file has no nodal coordinates "
                              "(neither 'coordx' nor packed 'coord')");

    /* add every node so that getfem point id == Exodus node id (0-based) */
    std::vector<size_type> node2pid(num_nodes);
    node_pts_.assign(num_nodes, base_node(dim_type(num_dim)));
    for (size_t i=0; i < num_nodes; ++i) {
      base_node P(num_dim);
      P[0] = X[i];
      if (num_dim >= 2) P[1] = Y[i];
      if (num_dim >= 3) P[2] = Z[i];
      node_pts_[i] = P;
      node2pid[i] = m.add_point(P, scalar_type(-1)); // no merge
    }

    /* block connectivity -> convexes; track the 1-based global Exodus element
       id of each convex (and its topology) so side sets can be resolved. */
    std::vector<size_type> elem_to_cv;
    std::vector<exo_topo> elem_to_topo;
    std::vector<bool> elem_to_shell;
    for (size_t b=0; b < num_blk; ++b) {
      std::ostringstream sc, sd1, sd2;
      sc << "connect" << (b+1); sd1 << "num_el_in_blk" << (b+1);
      sd2 << "num_nod_per_el" << (b+1);
      int v_cn, de, dn; size_t nel, npe;
      EXNC(nc_inq_varid(ncid_, sc.str().c_str(), &v_cn));
      EXNC(nc_inq_dimid(ncid_, sd1.str().c_str(), &de));
      EXNC(nc_inq_dimlen(ncid_, de, &nel));
      EXNC(nc_inq_dimid(ncid_, sd2.str().c_str(), &dn));
      EXNC(nc_inq_dimlen(ncid_, dn, &npe));

      // read elem_type via the length-checked helper (the attribute can be
      // longer than any fixed buffer in a foreign/malformed file)
      std::string etype = exo_get_att_text(ncid_, v_cn, "elem_type");
      GMM_ASSERT1(!etype.empty(),
                  "Exodus import: block " << (b+1) << " has no 'elem_type'");
      bool is_shell = exo_is_shell_type(etype);
      exo_topo topo = exo_topo_of_name(etype);
      bgeot::pgeometric_trans pgt = exo_pgt(topo, npe);
      GMM_ASSERT1(pgt->nb_points() == npe,
                  "Exodus block element type '" << etype << "' (" << npe
                  << " nodes) is not supported for import");

      /* inverse permutation: for each GetFEM geometric node, which Exodus
         local position supplies it */
      std::vector<base_node> gn(pgt->geometric_nodes().begin(),
                                pgt->geometric_nodes().end());
      std::vector<unsigned> perm = exo_build_perm(gn, exo_reference(topo, npe));
      std::vector<unsigned> inv(npe);
      for (size_type e=0; e < npe; ++e) inv[perm[e]] = unsigned(e);

      std::vector<size_type> cvnodes(npe);
      size_t chunk_elems =
        std::max<size_t>(1, EXO_IO_CHUNK / std::max<size_t>(1, npe));
      std::vector<int> conn(chunk_elems*npe);
      for (size_t first=0; first < nel; first += chunk_elems) {
        size_t ne = std::min(chunk_elems, nel - first);
        size_t start[2] = { first, 0 }, count[2] = { ne, npe };
        EXNC(nc_get_vara_int(ncid_, v_cn, start, count, conn.data()));
        for (size_t el=0; el < ne; ++el) {
          for (size_type g=0; g < npe; ++g) {
            int cnode = conn[el*npe + inv[g]];        // 1-based Exodus node id
            GMM_ASSERT1(cnode >= 1 && size_t(cnode) <= num_nodes,
                        "Exodus import: connectivity node id " << cnode
                        << " out of range [1," << num_nodes << "]");
            cvnodes[g] = node2pid[size_t(cnode) - 1];
          }
          size_type ic = m.add_convex(pgt, cvnodes.begin());
          elem_to_cv.push_back(ic);
          elem_to_topo.push_back(topo);
          elem_to_shell.push_back(is_shell);
        }
      }
    }

    /* side sets -> face regions (Exodus side-set id -> GetFEM region id) */
    int d_nss;
    if (nc_inq_dimid(ncid_, "num_side_sets", &d_nss) == NC_NOERR) {
      size_t n_ss; EXNC(nc_inq_dimlen(ncid_, d_nss, &n_ss));
      std::vector<int> ss_ids(n_ss);
      if (nc_inq_varid(ncid_, "ss_prop1", &v) == NC_NOERR)
        exo_get_ints(ncid_, v, n_ss, ss_ids.data());
      else for (size_t i=0; i < n_ss; ++i) ss_ids[i] = int(i+1);
      std::map<exo_topo, std::vector<unsigned> > finv_cache;
      for (size_t i=0; i < n_ss; ++i) {
        std::ostringstream se, sf, sd;
        se << "elem_ss" << (i+1); sf << "side_ss" << (i+1);
        sd << "num_side_ss" << (i+1);
        int v_e, v_s, dss; size_t nside;
        EXNC(nc_inq_dimid(ncid_, sd.str().c_str(), &dss));
        EXNC(nc_inq_dimlen(ncid_, dss, &nside));
        EXNC(nc_inq_varid(ncid_, se.str().c_str(), &v_e));
        EXNC(nc_inq_varid(ncid_, sf.str().c_str(), &v_s));
        std::vector<int> els(nside), sds(nside);
        exo_get_ints(ncid_, v_e, nside, els.data());
        exo_get_ints(ncid_, v_s, nside, sds.data());
        mesh_region reg = m.region(size_type(ss_ids[i]));
        for (size_t k=0; k < nside; ++k) {
          int elem_id = els[k];                 // 1-based global element id
          GMM_ASSERT1(elem_id >= 1 && size_t(elem_id) <= elem_to_cv.size(),
                      "Exodus import: side set " << ss_ids[i]
                      << " references element id " << elem_id
                      << " out of range [1," << elem_to_cv.size() << "]");
          size_t gid = size_t(elem_id);
          GMM_ASSERT1(!elem_to_shell[gid-1],
                      "Exodus import: side sets on SHELL elements are not "
                      "supported because Exodus shell side numbering is not "
                      "the same as GetFEM surface-element face numbering");
          exo_topo topo = elem_to_topo[gid-1];
          auto it = finv_cache.find(topo);
          if (it == finv_cache.end()) {
            std::vector<unsigned> fp = exo_face_perm(topo), finv(fp.size());
            for (size_type f=0; f < fp.size(); ++f) finv[fp[f]] = unsigned(f);
            it = finv_cache.emplace(topo, finv).first;
          }
          int side_id = sds[k];                 // 1-based Exodus side
          GMM_ASSERT1(side_id >= 1 && size_t(side_id) <= it->second.size(),
                      "Exodus import: side set " << ss_ids[i]
                      << " references side " << side_id << " out of range "
                      << "[1," << it->second.size() << "] for element "
                      << elem_id);
          size_t side = size_t(side_id);
          reg.add(elem_to_cv[gid-1], short_type(it->second[side-1]));
        }
      }
    }

    /* element sets -> convex regions (Exodus element-set id -> region id) */
    int d_nes;
    if (nc_inq_dimid(ncid_, "num_elem_sets", &d_nes) == NC_NOERR) {
      size_t n_es; EXNC(nc_inq_dimlen(ncid_, d_nes, &n_es));
      std::vector<int> es_ids(n_es);
      if (nc_inq_varid(ncid_, "els_prop1", &v) == NC_NOERR)
        exo_get_ints(ncid_, v, n_es, es_ids.data());
      else for (size_t i=0; i < n_es; ++i) es_ids[i] = int(i+1);
      for (size_t i=0; i < n_es; ++i) {
        std::ostringstream se, sd;
        se << "elem_els" << (i+1); sd << "num_ele_els" << (i+1);
        int v_e, des; size_t nel;
        EXNC(nc_inq_dimid(ncid_, sd.str().c_str(), &des));
        EXNC(nc_inq_dimlen(ncid_, des, &nel));
        EXNC(nc_inq_varid(ncid_, se.str().c_str(), &v_e));
        std::vector<int> els(nel);
        exo_get_ints(ncid_, v_e, nel, els.data());
        mesh_region reg = m.region(size_type(es_ids[i]));
        for (size_t k=0; k < nel; ++k) {
          int elem_id = els[k];
          GMM_ASSERT1(elem_id >= 1 && size_t(elem_id) <= elem_to_cv.size(),
                      "Exodus import: element set " << es_ids[i]
                      << " references element id " << elem_id
                      << " out of range [1," << elem_to_cv.size() << "]");
          size_t gid = size_t(elem_id);
          reg.add(elem_to_cv[gid-1]);
        }
      }
    }

    /* node sets -> recorded for node_set() (GetFEM has no node region) */
    int d_nns;
    if (nc_inq_dimid(ncid_, "num_node_sets", &d_nns) == NC_NOERR) {
      size_t n_ns; EXNC(nc_inq_dimlen(ncid_, d_nns, &n_ns));
      std::vector<int> ns_ids(n_ns);
      if (nc_inq_varid(ncid_, "ns_prop1", &v) == NC_NOERR)
        exo_get_ints(ncid_, v, n_ns, ns_ids.data());
      else for (size_t i=0; i < n_ns; ++i) ns_ids[i] = int(i+1);
      for (size_t i=0; i < n_ns; ++i) {
        std::ostringstream sn, sd;
        sn << "node_ns" << (i+1); sd << "num_nod_ns" << (i+1);
        int v_n, dns; size_t nn;
        EXNC(nc_inq_dimid(ncid_, sd.str().c_str(), &dns));
        EXNC(nc_inq_dimlen(ncid_, dns, &nn));
        EXNC(nc_inq_varid(ncid_, sn.str().c_str(), &v_n));
        std::vector<int> nodes(nn);
        exo_get_ints(ncid_, v_n, nn, nodes.data());
        std::vector<size_type> pts(nn);
        for (size_t k=0; k < nn; ++k) {
          int node_id = nodes[k];
          GMM_ASSERT1(node_id >= 1 && size_t(node_id) <= node2pid.size(),
                      "Exodus import: node set " << ns_ids[i]
                      << " references node id " << node_id
                      << " out of range [1," << node2pid.size() << "]");
          pts[k] = node2pid[size_t(node_id) - 1]; // 1-based -> 0-based node
        }
        node_set_ids_.push_back(ns_ids[i]);
        node_sets_[ns_ids[i]] = pts;
      }
    }
  }

  const std::vector<size_type> &exodus_import::node_set(int id) const {
    static const std::vector<size_type> empty;
    auto it = node_sets_.find(id);
    return (it == node_sets_.end()) ? empty : it->second;
  }

  void exodus_import::build_mesh_fem(const mesh &m, mesh_fem &mf) const {
    mf = mesh_fem(m);
    mf.set_qdim(1);
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
      pfem pf;
      // keep the matching SERENDIPITY element for incomplete geotrans, so the
      // mesh_fem has exactly the file's nodes (dof i <-> Exodus node i); a
      // complete classical_fem would add interior dofs absent from the file
      if (pgt == bgeot::Q2_incomplete_geotrans(2))
        pf = fem_descriptor("FEM_Q2_INCOMPLETE(2)");
      else if (pgt == bgeot::Q2_incomplete_geotrans(3))
        pf = fem_descriptor("FEM_Q2_INCOMPLETE(3)");
      else if (pgt == bgeot::prism_incomplete_P2_geotrans())
        pf = fem_descriptor("FEM_PRISM_INCOMPLETE_P2");
      else {
        short_type degree =
          (pgt->structure() != pgt->basic_structure()) ? 2 : 1;
        pf = classical_fem(pgt, degree, true);
      }
      mf.set_finite_element(cv, pf);
    }

    /* The imported mesh was built with point id == Exodus node id.  GetFEM's
       natural mesh_fem enumeration is not point-id ordered, so expose a reduced
       mesh_fem whose dof i is exactly Exodus node i. */
    size_type nn = node_pts_.size();
    GMM_ASSERT1(nn > 0 || m.convex_index().card() == 0,
                "Exodus import: read_mesh() must be called before build_mesh_fem()");
    std::vector<size_type> node_to_basic(nn, size_type(-1));
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv) {
      const bgeot::mesh_structure::ind_cv_ct &pts = m.ind_points_of_convex(cv);
      mesh_fem::ind_dof_ct dofs = mf.ind_basic_dof_of_element(cv);
      GMM_ASSERT1(pts.size() == dofs.size(),
                  "Exodus import: imported finite element does not have one "
                  "basic dof per Exodus node");
      for (size_type j=0; j < pts.size(); ++j) {
        size_type n = pts[j], d = dofs[j];
        GMM_ASSERT1(n < nn, "Exodus import: mesh point id " << n
                    << " is outside the Exodus node range [0," << nn << ")");
        if (node_to_basic[n] == size_type(-1)) node_to_basic[n] = d;
        else GMM_ASSERT1(node_to_basic[n] == d,
                         "Exodus import: node " << n
                         << " maps to inconsistent mesh_fem dofs");
      }
    }
    GMM_ASSERT1(mf.nb_basic_dof() == nn,
                "Exodus import: the mesh_fem has " << mf.nb_basic_dof()
                << " basic dofs but the file has " << nn
                << " nodes; isolated Exodus nodes cannot be represented by "
                "this mesh_fem");

    gmm::row_matrix<gmm::rsvector<scalar_type> > R(nn, mf.nb_basic_dof());
    gmm::row_matrix<gmm::rsvector<scalar_type> > E(mf.nb_basic_dof(), nn);
    for (size_type n=0; n < nn; ++n) {
      GMM_ASSERT1(node_to_basic[n] != size_type(-1),
                  "Exodus import: node " << n
                  << " is not attached to any imported element");
      R(n, node_to_basic[n]) = scalar_type(1);
      E(node_to_basic[n], n) = scalar_type(1);
    }
    mf.set_reduction_matrices(R, E);
  }

  void exodus_import::read_nodal_var_(const std::string &name, size_type step,
                                      std::vector<scalar_type> &U) const {
    GMM_ASSERT1(file_open_, "Exodus file is closed");
    size_type k = size_type(-1);
    for (size_type i=0; i < nodal_var_names_.size(); ++i)
      if (nodal_var_names_[i] == name) { k = i; break; }
    GMM_ASSERT1(k != size_type(-1),
                "no nodal variable named '" << name << "' in the Exodus file");
    GMM_ASSERT1(step < times_.size() || (times_.empty() && step == 0),
                "time step " << step << " out of range");
    int d_nodes; size_t num_nodes;
    EXNC(nc_inq_dimid(ncid_, "num_nodes", &d_nodes));
    EXNC(nc_inq_dimlen(ncid_, d_nodes, &num_nodes));
    std::ostringstream s; s << "vals_nod_var" << (k+1);
    int v; EXNC(nc_inq_varid(ncid_, s.str().c_str(), &v));
    U.assign(num_nodes, 0.);
    size_t start[2] = { step, 0 }, count[2] = { 1, num_nodes };
    EXNC(nc_get_vara_double(ncid_, v, start, count, U.data()));
  }

  void import_mesh_exodus(const std::string &fname, mesh &m) {
    exodus_import imp(fname);
    imp.read_mesh(m);
  }

} /* end of namespace getfem */

#endif /* GETFEM_HAVE_EXODUS */
