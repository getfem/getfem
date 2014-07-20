/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2001-2012 Yves Renard, Julien Pommier
 
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
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/

/**@file getfem_export.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>,
   @author  Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date October 15, 2001.
   @brief Export solutions to various formats.
*/
#ifndef GETFEM_EXPORT_H__
#define GETFEM_EXPORT_H__

#include "getfem_interpolation.h"
#include "getfem_mesh_slice.h"
#include <list>

namespace getfem {

  /* ********************************************************************* */
  /*                                                                       */
  /*  Save a solution in a file with a Pk interpolation.                   */
  /*                                                                       */
  /* ********************************************************************* */

  inline std::string remove_spaces(const std::string &s) {
    std::string s2(s);
    for (unsigned i=0; i < s.size(); ++i)
      if (s2[i] <= ' ') s2[i] = '_';
    return s2;
  }

  /** @brief VTK export.

      export class to VTK ( http://www.kitware.com/vtk.html ) file format
      (not the XML format, but the old format)

      A vtk_export can store multiple scalar/vector fields.
  */
  class vtk_export {
  protected:
    std::ostream &os;
    char header[256]; // hard limit in vtk
    bool ascii;
    const stored_mesh_slice *psl;
    std::unique_ptr<mesh_fem> pmf;
    dal::bit_vector pmf_dof_used;
    std::vector<unsigned> pmf_cell_type;
    std::ofstream real_os;
    dim_type dim_;
    bool reverse_endian;
    enum { EMPTY, HEADER_WRITTEN, STRUCTURE_WRITTEN, IN_CELL_DATA,
           IN_POINT_DATA } state;
  public:
    typedef enum { VTK_VERTEX = 1, VTK_LINE = 3, VTK_QUADRATIC_EDGE = 21,
                   VTK_TRIANGLE = 5, VTK_QUADRATIC_TRIANGLE = 22,
                   VTK_PIXEL = 8, VTK_QUAD = 9, VTK_QUADRATIC_QUAD = 23,
                   VTK_TETRA = 10, VTK_QUADRATIC_TETRA = 24,
                   VTK_WEDGE = 13, /*VTK_QUADRATIC_WEDGE = 26,*/
                   VTK_VOXEL = 11, VTK_HEXAHEDRON = 12,
                   VTK_QUADRATIC_HEXAHEDRON = 25 } vtk_cell_type;
    vtk_export(const std::string& fname, bool ascii_ = false);
    vtk_export(std::ostream &os_, bool ascii_ = false);

    /** should be called before write_*_data */
    void exporting(const mesh& m);
    void exporting(const mesh_fem& mf);
    void exporting(const stored_mesh_slice& sl);

    /** the header is the second line of text in the exported file,
       you can put whatever you want -- call this before any write_dataset
       or write_mesh */
    void set_header(const std::string& s);
    void write_mesh();

    /** append a new scalar or vector field defined on mf to the .vtk file.  If
        you are exporting a slice, or if mf != get_exported_mesh_fem(), U will
        be interpolated on the slice, or on get_exported_mesh_fem().

        Note that vectors should be written AFTER scalars, and tensors
        after vectors

        NO SPACE ALLOWED in 'name' */
    template<class VECT> void write_point_data(const getfem::mesh_fem &mf,
                                               const VECT& U0,
                                               const std::string& name);

    /** append a new scalar or vector field to .vtk file. The Uslice vector is
        the field interpolated on the exported mesh_slice This function should
        not be used if you are not exporting a slice!  NO SPACE ALLOWED in
        'name' */
    template<class VECT> void write_sliced_point_data(const VECT& Uslice,
                                                      const std::string& name,
                                                      getfem::size_type qdim=1);
    /** export data which is constant over each element. You should not use
        this function if you are exporting a slice.  U should have
        convex_index().card() elements.  */

    template<class VECT> void write_cell_data(const VECT& U,
                                              const std::string& name, 
                                              getfem::size_type qdim = 1);
    /** export a data_set correspounding to measures of quality for each convex
        of the supplied mesh (which should have the same number of convex than
        the one used in the vtk_export)

        If a slice is being exported, the convex quality is written as
        point_data (TO IMPROVE ONEDAY), if a mesh/mesh_fem is being exported,
        it is written as cell_data
    */
    void write_mesh_quality(const mesh &m);
    void write_normals();
    const stored_mesh_slice& get_exported_slice() const;
    const mesh_fem& get_exported_mesh_fem() const;
  private:
    void init();
    void check_header();
    void write_mesh_structure_from_slice();
    void write_mesh_structure_from_mesh_fem();
    void switch_to_cell_data();
    void switch_to_point_data();
    template<class T> void write_val(T v);
    template<class V> void write_vec(V p, getfem::size_type qdim);
    template<class IT> void write_3x3tensor(IT p);
    void write_separ();
    template<class VECT> void write_dataset_(const VECT& U,
                                             const std::string& name,
                                             getfem::size_type qdim,
                                             bool cell_data=false);
  };

  template<class T> void vtk_export::write_val(T v) {
    if (ascii) os << " " << v;
    else {
      char *p = (char*)&v;
      if (reverse_endian)
        for (size_type i=0; i < sizeof(v)/2; ++i)
          std::swap(p[i], p[sizeof(v)-i-1]);
      os.write(p, sizeof(T));
    }
  }

  template<class IT> void vtk_export::write_vec(IT p, getfem::size_type qdim) {
    float v[3];
    for (size_type i=0; i < qdim; ++i) {
      v[i] = float(p[i]);
    }
    for (size_type i=qdim; i < 3; ++i) v[i] = 0.0f;
    write_val(v[0]);write_val(v[1]);write_val(v[2]);
  }

  template<class IT> void vtk_export::write_3x3tensor(IT p) {
    float v[3][3];
    memset(v, 0, sizeof v);
    for (size_type i=0; i < dim_; ++i) {
      for (size_type j=0; j < dim_; ++j)
        v[i][j] = float(p[i + j*dim_]);
    }
    for (size_type i=0; i < 3; ++i) {
      for (size_type j=0; j < 3; ++j) {
        write_val(v[i][j]);
      }
      if (ascii) os << "\n";
    }
  }

  template<class VECT>
  void vtk_export::write_point_data(const getfem::mesh_fem &mf, const VECT& U,
                                    const std::string& name) {
    size_type Q = (gmm::vect_size(U) / mf.nb_dof()) * mf.get_qdim();
    size_type qdim = mf.get_qdim();
    if (psl) {
      std::vector<scalar_type> Uslice(Q*psl->nb_points());
      psl->interpolate(mf, U, Uslice);
      write_dataset_(Uslice, name, qdim);
    } else {
      std::vector<scalar_type> V(pmf->nb_dof() * Q);
      if (&mf != &(*pmf)) {
        interpolation(mf, *pmf, U, V);
      } else gmm::copy(U,V);
      size_type cnt = 0;
      for (dal::bv_visitor d(pmf_dof_used); !d.finished(); ++d, ++cnt) {
        if (cnt != d)
          for (size_type q=0; q < Q; ++q) {
            V[cnt*Q + q] = V[d*Q + q];
          }
      }
      V.resize(Q*pmf_dof_used.card());
      write_dataset_(V, name, qdim);
    }
  }

  template<class VECT>
  void vtk_export::write_cell_data(const VECT& U, const std::string& name, getfem::size_type qdim) {
    write_dataset_(U, name, qdim, true);
  }

  template<class VECT>
  void vtk_export::write_sliced_point_data(const VECT& U,
                                           const std::string& name,
                                           getfem::size_type qdim) {
    write_dataset_(U, name, qdim, false);
  }

  template<class VECT>
  void vtk_export::write_dataset_(const VECT& U, const std::string& name,
                                  getfem::size_type qdim,
                                  bool cell_data) {
    write_mesh();
    size_type nb_val = 0;
    if (cell_data) {
      switch_to_cell_data();
      nb_val = psl ? psl->linked_mesh().convex_index().card()
                   : pmf->linked_mesh().convex_index().card();
    } else {
      switch_to_point_data();
      nb_val = psl ? psl->nb_points() : pmf_dof_used.card();
    }
    //size_type Q = gmm::vect_size(U) / nb_val;
    size_type Q = qdim;
    GMM_ASSERT1(gmm::vect_size(U) == nb_val*Q,
                "inconsistency in the size of the dataset: "
                << gmm::vect_size(U) << " != " << nb_val << "*" << Q);
    write_separ();
    if (Q == 1) {
      os << "SCALARS " << remove_spaces(name) << " float 1\n";
      os << "LOOKUP_TABLE default\n";
      for (size_type i=0; i < nb_val; ++i) {
        write_val(float(U[i]));
      }
    } else if (Q <= 3) {
      os << "VECTORS " << remove_spaces(name) << " float\n";
      for (size_type i=0; i < nb_val; ++i) {
        write_vec(U.begin() + i*Q, Q);
      }
    } else if (Q == gmm::sqr(dim_)) {
      /* tensors : coef are supposed to be stored in FORTRAN order
         in the VTK file, they are written with C (row major) order
       */
      os << "TENSORS " << remove_spaces(name) << " float\n";
      for (size_type i=0; i < nb_val; ++i) {
        write_3x3tensor(U.begin() + i*Q);
      }
    } else GMM_ASSERT1(false, "vtk does not accept vectors of dimension > 3");
    write_separ();
  }


  /** @brief A (quite large) class for exportation of data to IBM OpenDX.

                     http://www.opendx.org/

      This class is more capable than the VTK export, as it is
      possible to export many different meshes/slices, with their
      edges, datasets, and create series of dataset for animations
      etc, in a single '.dx' file.

      Moreover, it is able to reopen a '.dx' file and append new data
      into it.  Hence it is possible, if many time-steps are to be
      saved, to view intermediate results in OpenDX during the
      computation.
   */
  class dx_export {
    std::ostream &os;
    char header[256];
    bool ascii;
    const stored_mesh_slice *psl;
    bool psl_use_merged; /* flag enabled if we merge the points of
                            psl before export */
    std::unique_ptr<mesh_fem> pmf;
    dal::bit_vector pmf_dof_used;
    std::vector<unsigned> pmf_cell_type;
    std::fstream real_os;
    dim_type dim_, connections_dim;
    struct dxSeries {
      std::string name;
      std::list<std::string> members;
    };
    struct dxObject {
      std::string name;
      std::string mesh;
    };
    struct dxMesh {
      unsigned flags;
      typedef enum { NONE=0, WITH_EDGES=1, STRUCTURE_WRITTEN=2 } flags_t;
      std::string name;
      dxMesh() : flags(NONE) {}
    };
    std::list<dxObject> objects;
    std::list<dxSeries> series;
    std::list<dxMesh> meshes;
    bool header_written;
  public:
    dx_export(const std::string& fname, bool ascii_ = false,
              bool append_ = false);
    dx_export(std::ostream &os_, bool ascii_ = false);
    ~dx_export(); /* the file is not complete until the destructor
                     has been executed */
    void exporting(const mesh& m, std::string name = std::string());
    void exporting(const mesh_fem& mf, std::string name = std::string());
    void exporting(const stored_mesh_slice& sl, bool merge_points = true,
                   std::string name = std::string());
    /** append edges information (if you want to draw the mesh and are
        using a refined slice. Should be called just after exporting(..)  */
    void exporting_mesh_edges(bool with_slice_edge = true);

    /** the header is the second line of text in the exported file,
       you can put whatever you want -- call this before any write_dataset
       or write_mesh */
    void set_header(const std::string& s);
    void write_mesh();
    /** add an object (typically the name of a data field) to a
       'series', i.e.  an aggregate of consecutive objects. Using
       'series' is useful for animations in opendx

       If 'field_name' corresponds to a data_set whose mesh edges have
       been exported, a second series called serie_name + '_edges'
       will be filled, which will allow you to view the mesh edges.
    */
    void serie_add_object(const std::string& serie_name,
                          const std::string& object_name);
    void serie_add_object(const std::string& serie_name)
    { serie_add_object(serie_name, current_data_name()); }
    /** return the name of current mesh (use exporting(...) to change
        the current mesh) */
    std::string current_mesh_name() { return current_mesh().name; }
    /** return the name of last written data_set */
    std::string current_data_name() { return current_data().name; }
    template<class VECT> void
    write_point_data(const getfem::mesh_fem &mf,
                     const VECT& U0, std::string name = std::string());
    template<class VECT> void
    write_sliced_point_data(const VECT& Uslice,
                            std::string name = std::string());
    /* TOBEDONE !!!!!!!!!!!
       template<class VECT> void
    write_cell_data(const VECT& U, std::string name = std::string());
    void write_mesh_quality(const mesh &m);*/
    void write_normals();
    const stored_mesh_slice& get_exported_slice() const;
    const mesh_fem& get_exported_mesh_fem() const;

  private:
    void init();
    void reread_metadata();
    void update_metadata(std::ios::pos_type);
    void write_series();
    void serie_add_object_(const std::string &serie_name,
                           const std::string &object_name);
    void write_separ();
    std::string default_name(std::string s, int count,
                             const char *default_prefix) {
      if (s.size() == 0) {
        std::stringstream ss; ss << default_prefix << count; return ss.str();
      } else return s;
    }
    template<class T> void write_val(T v) {
      if (ascii) os << " " << v;
      else os.write((char*)&v, sizeof(T));
    }
    static const char* endianness() {
      static int i=0x12345678;
      char *p = (char*)&i;
      if (*p == 0x12) return "msb";
      else if (*p == 0x78) return "lsb";
      else return "this is very strange..";
    }
    bool new_mesh(std::string &name);
    std::list<dxMesh>::iterator get_mesh(const std::string& name,
                                         bool raise_error = true);
    std::list<dxObject>::iterator get_object(const std::string& name,
                                             bool raise_error = true);
    dxMesh &current_mesh() {
      if (meshes.size()) return meshes.back();
      else GMM_ASSERT1(false, "no mesh!");
    }
    dxObject &current_data() {
      if (objects.size()) return objects.back();
      else GMM_ASSERT1(false, "no data!");
    }

    std::string name_of_pts_array(const std::string &meshname)
    { return meshname + std::string("_pts"); }
    std::string name_of_conn_array(const std::string &meshname)
    { return meshname + std::string("_conn"); }
    std::string name_of_edges_array(const std::string &meshname)
    { return meshname + std::string("_edges"); }
    void check_header();
    const char *dxname_of_convex_structure(bgeot::pconvex_structure cvs);
    void write_convex_attributes(bgeot::pconvex_structure cvs);
    void write_mesh_structure_from_slice();
    void write_mesh_structure_from_mesh_fem();
    void write_mesh_edges_from_slice(bool with_slice_edge);
    void write_mesh_edges_from_mesh();
    template <class VECT>
    void smooth_field(const VECT& U, base_vector &sU);
    template<class VECT>
    void write_dataset_(const VECT& U, std::string name, bool cell_data=false);
  };

  template <class VECT>
  void dx_export::smooth_field(const VECT& U, base_vector &sU) {
    size_type Q = gmm::vect_size(U) / psl->nb_points();
    sU.clear(); sU.resize(Q*psl->nb_merged_nodes());
    for (size_type i=0; i < psl->nb_merged_nodes(); ++i) {
      for (size_type j=0; j < psl->merged_point_cnt(i); ++j)
        for (size_type q=0; q < Q; ++q)
          sU[i*Q+q] += U[psl->merged_point_nodes(i)[j].pos*Q+q];
      for (size_type q=0; q < Q; ++q)
        sU[i*Q+q] /= double(psl->merged_point_cnt(i));
    }
  }

  template<class VECT>
  void dx_export::write_point_data(const getfem::mesh_fem &mf, const VECT& U,
                                   std::string name) {
    size_type Q = (gmm::vect_size(U) / mf.nb_dof())*mf.get_qdim();
    if (psl) {
      std::vector<scalar_type> Uslice(Q*psl->nb_points());
      psl->interpolate(mf, U, Uslice);
      write_sliced_point_data(Uslice,name);
    } else {
      std::vector<scalar_type> V(pmf->nb_dof() * Q);
      if (&mf != &(*pmf)) {
        interpolation(mf, *pmf, U, V);
      } else gmm::copy(U,V);
      size_type cnt = 0;
      for (dal::bv_visitor d(pmf_dof_used); !d.finished(); ++d, ++cnt) {
        if (cnt != d)
          for (size_type q=0; q < Q; ++q) {
            V[cnt*Q + q] = V[d*Q + q];
          }
      }
      V.resize(Q*pmf_dof_used.card());
      write_dataset_(V, name);
    }
  }

  template<class VECT> void
  dx_export::write_sliced_point_data(const VECT& Uslice, std::string name) {
    if (!psl_use_merged)
      write_dataset_(Uslice, name, false);
    else {
      base_vector Umerged; smooth_field(Uslice,Umerged);
      write_dataset_(Umerged, name, false);
    }
  }

  template<class VECT> void
  dx_export::write_dataset_(const VECT& U, std::string name, bool cell_data) {
    write_mesh();
    objects.push_back(dxObject());
    name = default_name(name, int(objects.size()), "gf_field");
    objects.back().name = name;
    objects.back().mesh = current_mesh_name();
    size_type nb_val = 0;
    if (cell_data) {
      nb_val = psl ? psl->linked_mesh().convex_index().card()
                   : pmf->linked_mesh().convex_index().card();
    } else {
      nb_val = psl ? (psl_use_merged ? psl->nb_merged_nodes() : psl->nb_points())
                   : pmf_dof_used.card();
    }
    size_type Q = gmm::vect_size(U) / nb_val;
    GMM_ASSERT1(gmm::vect_size(U) == nb_val*Q,
                "inconsistency in the size of the dataset: "
                << gmm::vect_size(U) << " != " << nb_val << "*" << Q);

    os << "\nobject \"" << name << "_data\" class array type float rank ";
    if (Q == 1) { os << "0"; } /* scalar data */
    else if (Q == 4) { os << "2 shape 2 2"; } /* or 2x2 tensor data */
    else if (Q == 9) { os << "2 shape 3 3"; } /* or 2x2 tensor data */
    else { os << "1 shape " << Q; } /* fallback: vector data */
    os << " items " << nb_val;
    if (!ascii) os << " " << endianness() << " binary";
    os << " data follows" << endl;
    for (size_type i=0; i < nb_val*Q; ++i) {
      write_val(float(U[i]));
      if (((i+1) % (Q > 1 ? Q : 10)) == 0) write_separ();
    }
    write_separ();

    if (!cell_data)
      os << "\n  attribute \"dep\" string \"positions\"\n";
    else os << "\n  attribute \"dep\" string \"connections\"\n";
    os << "\n";

    if (current_mesh().flags & dxMesh::WITH_EDGES) {
      os << "\nobject \"" << name << "_edges\" class field\n"
         << "  component \"positions\" value \""
         << name_of_pts_array(current_mesh_name()) << "\"\n"
         << "  component \"connections\" value \""
         << name_of_conn_array(name_of_edges_array(current_mesh_name()))
         << "\"\n"
         << "  component \"data\" value \"" << name << "_data\"\n";
    }

    /* write footer */
    os << "\nobject \"" << name << "\" class field\n"
       << "  component \"positions\" value \""
       << name_of_pts_array(current_mesh_name()) << "\"\n"
       << "  component \"connections\" value \""
       << name_of_conn_array(current_mesh_name()) << "\"\n"
       << "  component \"data\" value \"" << name << "_data\"\n";
  }

  /***************************************************************
      @brief POS export.

      export class to Gmsh post-processing file format.

      ( http://geuz.org/gmsh )

      A pos_export can store multiple scalar/vector/tensor fields.
  ****************************************************************/

  class pos_export {
  protected:
    std::ostream& os;
    char header[256];

    std::vector<std::vector<float> > pos_pts;
    std::vector<unsigned> pos_cell_type;
    std::vector<std::vector<unsigned> > pos_cell_dof;

    std::unique_ptr<mesh_fem> pmf;
    const stored_mesh_slice *psl;

    size_type view;
    dim_type dim;
    enum { EMPTY, HEADER_WRITTEN, STRUCTURE_WRITTEN, IN_CELL_DATA} state;
    std::ofstream real_os;

  public:
    typedef enum {
                   POS_PT = 0, //point
                   POS_LN = 1, //line
                   POS_TR = 2, //triangles
                   POS_QU = 3, //quadrangles
                   POS_SI = 4, //tetrahedra
                   POS_HE = 5, //hexahedra
                   POS_PR = 6  //prisms
    } pos_cell_types;

    pos_export(const std::string& fname);
    pos_export(std::ostream& osname);

    void set_header(const std::string& s);

    void exporting(const mesh& m);
    void exporting(const mesh_fem& mf);
    void exporting(const stored_mesh_slice& sl);

    void write(const mesh& m, const std::string& name="");
    void write(const mesh_fem& mf, const std::string& name="");
    void write(const stored_mesh_slice& sl, const std::string& name="");

    template <class VECT>
    void write(const mesh_fem& mf,const VECT& U, const std::string& name);
    template <class VECT>
    void write(const stored_mesh_slice& sl,const VECT& U, const std::string& name);

  private:
    void init();
    void check_header();

    template <class VECT>
    void write(const VECT& V, const size_type qdim_v);

    template <class VECT>
    void write_cell(const int& t, const std::vector<unsigned>& dof,
                                  const VECT& val);
  };

  template <class VECT>
  void pos_export::write(const mesh_fem& mf,const VECT& U,
                         const std::string& name){
    check_header();
    exporting(mf);

    os << "View \"" << name.c_str() <<"\" {\n";

    size_type nb_points = mf.nb_dof()/mf.get_qdim();
    size_type qdim_u = gmm::vect_size(U)/nb_points;
    if (psl){
      std::vector<scalar_type> Uslice(psl->nb_points()*qdim_u);
      psl->interpolate(mf, U, Uslice);
      qdim_u = gmm::vect_size(Uslice)/psl->nb_points();
      write(Uslice, qdim_u);
    }else {
      std::vector<scalar_type> V(pmf->nb_dof()*qdim_u);
      if (&mf != &(*pmf)) {
        interpolation(mf, *pmf, U, V);
      } else gmm::copy(U,V);
      /*for (dal::bv_visitor d(pmf_dof_used); !d.finished(); ++d, ++cnt) {
        if (cnt != d)
          for (size_type q=0; q < Q; ++q) {
            V[cnt*Q + q] = V[d*Q + q];
          }
      }
      V.resize(Q*pmf_dof_used.card());*/
      nb_points = pmf->nb_dof()/pmf->get_qdim();
      qdim_u = gmm::vect_size(V)/nb_points;
      write(V, qdim_u);
    }

    os << "};\n";
    os << "View[" << view << "].ShowScale = 1;\n";
    os << "View[" << view << "].ShowElement = 0;\n";
    os << "View[" << view << "].DrawScalars = 1;\n";
    os << "View[" << view << "].DrawVectors = 1;\n";
    os << "View[" << view++ << "].DrawTensors = 1;\n";
  }

  template <class VECT>
  void pos_export::write(const stored_mesh_slice& sl,const VECT& V,
                         const std::string& name){
    check_header();
    exporting(sl);

    os << "View \"" << name.c_str() <<"\" {\n";

    size_type qdim_v = gmm::vect_size(V)/psl->nb_points();
    write(V, qdim_v);

    os << "};\n";
    os << "View[" << view << "].ShowScale = 1;\n";
    os << "View[" << view << "].ShowElement = 0;\n";
    os << "View[" << view << "].DrawScalars = 1;\n";
    os << "View[" << view << "].DrawVectors = 1;\n";
    os << "View[" << view++ << "].DrawTensors = 1;\n";
  }

  template <class VECT>
  void pos_export::write(const VECT& V, const size_type qdim_v){
    int t;
    std::vector<unsigned> cell_dof;
    std::vector<scalar_type> cell_dof_val;
    for (size_type cell = 0; cell < pos_cell_type.size(); ++cell) {
      t = pos_cell_type[cell];
      cell_dof = pos_cell_dof[cell];
      cell_dof_val.resize(cell_dof.size()*qdim_v, scalar_type(0));
      for (size_type i=0; i< cell_dof.size(); ++i)
        for (size_type j=0; j< qdim_v; ++j)
          cell_dof_val[i*qdim_v+j] = scalar_type(V[cell_dof[i]*qdim_v+j]);
      write_cell(t,cell_dof,cell_dof_val);
    }
  }

  template <class VECT>
  void pos_export::write_cell(const int& t, const std::vector<unsigned>& dof,
                                            const VECT& val){
    size_type qdim_cell = val.size()/dof.size();
    size_type dim3D = size_type(-1);
    if (1==qdim_cell){
      dim3D = size_type(1);
      os << "S";
    } else if (2==qdim_cell || 3==qdim_cell){
      dim3D = size_type(3);
      os << "V";
    } else if (4<=qdim_cell && qdim_cell<=9){
      dim3D = size_type(9);
      os << "T";
    }
    switch (t){
      case POS_PT: os << "P("; break; // point
      case POS_LN: os << "L("; break; // line
      case POS_TR: os << "T("; break; // triangle
      case POS_QU: os << "Q("; break; // quadrangle
      case POS_SI: os << "S("; break; // tetrahedra (simplex)
      case POS_HE: os << "H("; break; // hexahedra
      case POS_PR: os << "I("; break; // prism
    }
    for (size_type i=0; i<dof.size(); ++i){
      for(size_type j=0; j<dim; ++j){
        if(0!=i || 0!=j) os << ",";
        os << pos_pts[dof[i]][j];
      }
      for (size_type j=dim; j<3; ++j){
        os << ",0.00";
      }
    }

    os << "){";
    for (size_type i=0; i<dof.size(); ++i){
      for(size_type j=0; j<qdim_cell; ++j){
        if(0!=i || 0!=j) os << ",";
        os << val[i*qdim_cell+j];
      }
      for (size_type j=qdim_cell; j< dim3D; ++j){
        os << ",0.00";
      }
    }
    os << "};\n";
  }
}  /* end of namespace getfem. */

#endif /* GETFEM_EXPORT_H__  */
