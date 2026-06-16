/* -*- c++ -*- (enables emacs c++ mode) */
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

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/**@file getfem_exodus.h
   @brief Export/import solutions to/from the Exodus II file format.

   Exodus II is the finite-element database format defined by Sandia's SEACAS
   project. It is a naming convention layered on top of NetCDF, so this
   implementation writes/reads the documented Exodus entities directly through
   the NetCDF C API (no dependency on the SEACAS library itself).

   This whole file is only active when GetFEM is configured with
   --enable-exodus (which defines GETFEM_HAVE_EXODUS and links NetCDF).
*/
#ifndef GETFEM_EXODUS_H__
#define GETFEM_EXODUS_H__

#include "getfem_export.h"

namespace getfem {

#ifdef GETFEM_HAVE_EXODUS

  /** @brief Exodus II export.

      Mirrors the logic of getfem::vtk_export: an internal classical mesh_fem
      (of order 1 or 2, suitable for Exodus element types) is built from the
      exported mesh/mesh_fem, fields are interpolated onto it, and only the
      degrees of freedom that map to Exodus nodes are written.

      A single Exodus file natively stores a transient (time-dependent) series:
      call set_time() to open a new time step, then write_point_data() for each
      field at that step. The file is finalised (and the buffered transient data
      flushed) when close() is called or the object is destroyed.
  */
  class exodus_export {
  protected:
    std::string fname_;
    int ncid_;                 // NetCDF file id (valid while file_open_)
    bool file_open_;
    char title_[256];
    int d_num_nodes_, d_time_, d_len_name_; // NetCDF dim ids reused on close()

    std::unique_ptr<mesh_fem> pmf_; // classical export mesh_fem
    dal::bit_vector pmf_dof_used_;  // basic dofs actually exported as nodes
    std::vector<int> dofmap_;       // basic dof -> 0-based compacted node index
    size_type nb_nodes_;
    dim_type dim_;
    std::map<int, std::string> region_names_; // region id -> Exodus block/set name

    // One Exodus element block per (volume region, element type). When no
    // volume region covers a convex it falls in a default block keyed by type.
    struct block_info {
      std::string ex_type;            // e.g. "TETRA10"
      size_type nodes_per_elem;
      std::vector<unsigned> perm;     // perm[exodus_local] = getfem local dof
      std::vector<size_type> convexes;// convexes (of pmf_) in this block
      int region_id;                  // GetFEM volume-region id, or -1 if none
      int block_id;                   // Exodus block id (eb_prop1); == region_id
                                      // when possible, else a free positive id
    };
    std::vector<block_info> blocks_;

    // Transient streaming: each step's slab is written straight to the open
    // file. The variable set is discovered from step 0 (buffered), declared
    // once, then later steps stream with no in-memory history.
    std::vector<std::string> var_names_;    // scalar nodal variable names (order)
    std::map<std::string, int> var_id_;     // name -> NetCDF varid (once declared)
    std::map<std::string, std::vector<scalar_type> > step0_; // step-0 buffer
    int v_time_;                            // NetCDF varid of time_whole
    int v_nod_names_, v_elem_names_, v_elem_var_tab_;
    bool vars_declared_;
    bool transient_definitions_ready_;
    bool append_mode_;                      // reopened an existing file to append
    size_type file_nb_nodes_;               // num_nodes read from a reopened file
    int append_base_step_;                  // file's last step on append (-1 else)
    int cur_step_;                          // -1 until first set_time/write
    scalar_type cur_time_;
    bool cur_time_valid_, cur_time_written_;
    std::vector<int> elem_var_id_;          // varid of vals_elem_var1eb{b}, per block
    std::vector<scalar_type> col_;          // reused gather buffer for one component
    std::vector<std::string> step_written_; // fields written in the current step
    bool region_field_;                     // write the "region" element variable?
    bool compression_;                      // create a compressed NetCDF4 file?
    bool compression_explicit_;             // user explicitly requested mode?
    int compression_level_;
    std::string file_mesh_fingerprint_;     // fingerprint read from append target

    enum { EMPTY, EXPORTING, MESH_WRITTEN } state_;

  public:
    /** Open `fname` for export. If `append` is true the file is reopened and
        further time steps are appended to it (it must already contain the mesh
        written for the same mesh_fem); otherwise it is created. */
    exodus_export(const std::string &fname, bool append = false);
    ~exodus_export();

    /** should be called before write_mesh / write_point_data */
    void exporting(const mesh &m);
    void exporting(const mesh_fem &mf);

    /** the Exodus database title (truncated to 80 chars by the format) */
    void set_header(const std::string &s);

    /** give GetFEM region `id` a human-readable name, used as the name of the
        matching Exodus element block / side set / element set / node set.
        Call before write_mesh(). */
    void set_region_name(int id, const std::string &name);

    /** also write a "region" element (cell) variable holding each element's
        block id, so a viewer can colour regions by a discrete per-cell value.
        Off by default: it is a time-constant field, so storing it at every
        transient step would bloat the file. Call before write_mesh(). */
    void enable_region_field(bool on = true);

    /** Set the deflate level (1..9) of the NetCDF4/CLASSIC_MODEL compression
        used for newly created files. Compression is ON by default when GetFEM
        was built with usable NetCDF4/HDF5 deflate support; otherwise files
        default to classic 64-bit offset. Pass level 0 to force classic output.
        A positive level requires NetCDF4/HDF5 support. Must be called before
        the file is created (not on append). */
    void enable_compression(int level = 1);

    /** Predeclare a nodal variable before write_mesh(), avoiding a later
        NetCDF redefine. Vector fields use the same component names as
        write_point_data(): name_x, name_y, name_z or name_0... */
    void declare_point_data(const std::string &name, size_type qdim = 1);

    /** write the mesh (coordinates + element blocks + connectivity). */
    void write_mesh();

    /** open a new transient time step at value t. Subsequent write_point_data
        calls attach to this step. */
    void set_time(scalar_type t);

    /** append a scalar/vector field defined on mf to the current time step.
        Vector fields are stored as scalar components name_x, name_y, name_z
        (or name_0 .. for more than 3 components). NO SPACE ALLOWED in name. */
    template<class VECT>
    void write_point_data(const mesh_fem &mf, const VECT &U,
                          const std::string &name);

    /** flush all complete data written so far. A live reader such as ParaView
        can reload the file after this call and see completed time steps. */
    void sync();

    /** finalise and close the file (also done by the destructor). */
    void close();

    const mesh_fem &get_exported_mesh_fem() const;

  private:
    void init();
    void build_export_mesh_fem(const mesh_fem &mf);
    void compute_blocks();
    void open_for_append_();          // reopen an existing file to append steps
    std::string mesh_fingerprint_() const;
    void append_point_data_names_(const std::string &name, size_type qdim);
    void define_transient_vars_();    // define time/vars while already in define mode
    void write_transient_metadata_(); // write var names/tables after enddef
    void write_step_time_and_region_();
    void check_current_step_complete_() const;
    void finish_current_step_(bool do_sync);
    void deflate_var_(int varid) const;
    void declare_transient_vars_();   // define time/vars + flush buffered step 0
    void write_region_slab_(size_type step); // (constant) "region" element var
    // stream already-interpolated, compacted nodal values (size nb_nodes_*Q)
    void add_nodal_data_(const std::vector<scalar_type> &V, size_type Q,
                         const std::string &name);
  };

  /** @brief Exodus II import / field read-back.

      read_mesh() rebuilds a getfem mesh (and turns Exodus side sets into face
      regions). The transient nodal variables can then be read back into a
      classical mesh_fem matching the file (build_mesh_fem + read_nodal_var),
      which is what makes an export -> import round-trip verifiable.
  */
  class exodus_import {
  protected:
    int ncid_;
    bool file_open_;
    std::vector<std::string> nodal_var_names_;
    std::vector<scalar_type> times_;
    std::vector<base_node> node_pts_; // filled by read_mesh(), Exodus node order
    std::vector<int> node_set_ids_;
    std::map<int, std::vector<size_type> > node_sets_; // id -> getfem point ids

  public:
    exodus_import(const std::string &fname);
    ~exodus_import();

    /** rebuild the mesh from the file. Exodus side sets become GetFEM face
        regions and element sets become convex regions (both keyed by the set
        id); node sets are recorded separately (see node_set()), as GetFEM has
        no native node region. */
    void read_mesh(mesh &m);

    /** ids of the node sets in the file (valid after read_mesh()). */
    const std::vector<int> &node_set_ids() const { return node_set_ids_; }
    /** node ids (== getfem point ids) of node set `id`; empty if absent. */
    const std::vector<size_type> &node_set(int id) const;

    const std::vector<std::string> &nodal_var_names() const
    { return nodal_var_names_; }
    const std::vector<scalar_type> &times() const { return times_; }
    size_type nb_steps() const { return times_.size(); }
    /** node coordinates in Exodus node order (valid after read_mesh()); node i
        here corresponds to value i returned by read_nodal_var(). */
    const std::vector<base_node> &node_points() const { return node_pts_; }

    /** build a classical mesh_fem on m whose nodes coincide with the Exodus
        nodes (order chosen so dof i <-> Exodus node i). */
    void build_mesh_fem(const mesh &m, mesh_fem &mf) const;

    /** read nodal variable `name` at time step `step` (0-based) into U, indexed
        by Exodus node number (== dof of the mesh_fem built by build_mesh_fem). */
    template<class VECT>
    void read_nodal_var(const std::string &name, size_type step, VECT &U) const;

    void close();

  private:
    void read_nodal_var_(const std::string &name, size_type step,
                         std::vector<scalar_type> &U) const;
  };

  /** Hook used by getfem::import_mesh(..., "exodus", ...). */
  void import_mesh_exodus(const std::string &fname, mesh &m);

  /* --- template implementations --- */

  template<class VECT>
  void exodus_export::write_point_data(const mesh_fem &mf, const VECT &U,
                                       const std::string &name) {
    GMM_ASSERT1(state_ >= EXPORTING, "call exporting() first");
    size_type Q = (gmm::vect_size(U) / mf.nb_dof()) * mf.get_qdim();
    GMM_ASSERT1(Q >= 1, "empty data");
    // interpolate the field onto the export mesh_fem, then compact to the
    // exported nodes (identical bookkeeping to vtk_export::write_point_data).
    std::vector<scalar_type> V(pmf_->nb_dof() * Q);
    if (&mf != pmf_.get()) interpolation(mf, *pmf_, U, V);
    else gmm::copy(U, V);
    size_type cnt = 0;
    for (dal::bv_visitor d(pmf_dof_used_); !d.finished(); ++d, ++cnt) {
      if (cnt != d)
        for (size_type q=0; q < Q; ++q) V[cnt*Q + q] = V[d*Q + q];
    }
    V.resize(Q*nb_nodes_);
    add_nodal_data_(V, Q, name);
  }

  template<class VECT>
  void exodus_import::read_nodal_var(const std::string &name, size_type step,
                                     VECT &U) const {
    std::vector<scalar_type> tmp;
    read_nodal_var_(name, step, tmp);
    gmm::resize(U, tmp.size());
    gmm::copy(tmp, U);
  }

#endif /* GETFEM_HAVE_EXODUS */

} /* end of namespace getfem */

#endif /* GETFEM_EXODUS_H__ */
