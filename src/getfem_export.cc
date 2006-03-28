// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2001-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

#include <iomanip>
#include <dal_singleton.h>
#include <bgeot_comma_init.h>
#include <getfem_export.h>

namespace getfem
{
  /* -------------------------------------------------------------
   * VTK export 
   * ------------------------------------------------------------- */

  struct gf2vtk_dof_mapping : public std::vector<std::vector<unsigned> > {};

  static const std::vector<unsigned> &getfem_to_vtk_dof_mapping(int t) {
    gf2vtk_dof_mapping &dm = dal::singleton<gf2vtk_dof_mapping>::instance();
    if (dm.size() == 0) {
      dm.resize(30);
      bgeot::sc(dm[vtk_export::VTK_VERTEX]) = 0;
      bgeot::sc(dm[vtk_export::VTK_LINE]) = 0, 1;
      bgeot::sc(dm[vtk_export::VTK_QUADRATIC_EDGE]) = 0, 2, 1;
      bgeot::sc(dm[vtk_export::VTK_TRIANGLE]) = 0, 1, 2;
      bgeot::sc(dm[vtk_export::VTK_QUADRATIC_TRIANGLE]) = 0, 2, 5, 1, 4, 3;
      bgeot::sc(dm[vtk_export::VTK_QUAD]) = 0, 1, 3, 2;
      bgeot::sc(dm[vtk_export::VTK_PIXEL]) = 0, 1, 2, 3;
      bgeot::sc(dm[vtk_export::VTK_QUADRATIC_QUAD]) = 0, 2, 8, 6, 1, 5, 7, 3;
      bgeot::sc(dm[vtk_export::VTK_TETRA]) = 0, 1, 2, 3;
      bgeot::sc(dm[vtk_export::VTK_QUADRATIC_TETRA]) = 0, 2, 5, 9, 1, 4, 3, 6, 7, 8;
      bgeot::sc(dm[vtk_export::VTK_WEDGE]) = 0, 1, 2, 3, 4, 5;
      //bgeot::sc(dm[vtk_export::VTK_QUADRATIC_WEDGE]) = 0, 6, 1, 7, 2, 8, 3, 4, 5;
      bgeot::sc(dm[vtk_export::VTK_VOXEL]) = 0, 1, 2, 3, 4, 5, 6, 7;
      bgeot::sc(dm[vtk_export::VTK_HEXAHEDRON]) = 0, 1, 3, 2, 4, 5, 7, 6;
      bgeot::sc(dm[vtk_export::VTK_QUADRATIC_HEXAHEDRON]) = 0, 2, 8, 6, 18, 20, 26, 24, 1, 5, 7, 3, 19, 23, 25, 21, 9, 11, 17, 15;
    }
    return dm[t];
  }
  
  vtk_export::vtk_export(std::ostream &os_, bool ascii_)
    : os(os_), ascii(ascii_) { init(); } 
  
  vtk_export::vtk_export(const std::string& fname, bool ascii_)
    : os(real_os), ascii(ascii_), real_os(fname.c_str()) {
    if (!real_os) DAL_THROW(failure_error, "impossible to write to vtk file '"
			    << fname << "'");
    init();
  }

  void vtk_export::init() {
    static int test_endian = 0x01234567;
    strcpy(header, "Exported by getfem++");
    psl = 0; dim_ = dim_type(-1);
    if (*((char*)&test_endian) == 0x67)
      reverse_endian = true;
    else reverse_endian = false;
    state = EMPTY;
  }
  
  void vtk_export::switch_to_cell_data() {
    if (state != IN_CELL_DATA) {
      state = IN_CELL_DATA;
      write_separ();
      if (psl) {
        os << "CELL_DATA " << psl->nb_simplexes(0) + psl->nb_simplexes(1) + psl->nb_simplexes(2) + psl->nb_simplexes(3) << "\n";
      } else {
        os << "CELL_DATA " << pmf->convex_index().card() << "\n";
      }
      write_separ();
    }
  }

  void vtk_export::switch_to_point_data() {
    if (state != IN_POINT_DATA) {
      state = IN_POINT_DATA;
      write_separ();
      if (psl) {
        write_separ(); os << "POINT_DATA " << psl->nb_points() << "\n";
      } else {
        os << "POINT_DATA " << pmf_dof_used.card() << "\n";
      }
      write_separ();
    }
  }
  

  /* try to check if a quad or hexahedric cell is "rectangular" and oriented
     along the axes */
  template<typename C> static bool check_voxel(const C& c) {
    scalar_type h[3];
    unsigned N = c[0].size();
    if (c.size() != (1U << N)) return false;
    const base_node P0 = c[0];
    h[0] = c[1][0] - P0[0];
    h[1] = c[2][0] - P0[0];
    if (c.size() != 4) h[2] = c[4][0] - P0[0];
    for (unsigned i=1; i < c.size(); ++i) {
      const base_node d = c[i] - P0;
      for (unsigned j=0; j < N; ++j) 
        if (gmm::abs(d[j]) > 1e-7*h[j] && gmm::abs(d[j] - h[j]) > 1e-7*h[j]) 
          return false;
    }
    return true;
  }


  void vtk_export::exporting(const stored_mesh_slice& sl) {
    psl = &sl; dim_ = sl.dim();
    if (psl->dim() > 3) DAL_THROW(dal::failure_error, "attempt to export a " << int(dim_) << "D slice (not supported)");
  }
  
  void vtk_export::exporting(const mesh& m) {
    dim_ = m.dim();
    if (dim_ > 3) DAL_THROW(dal::failure_error, "attempt to export a " << int(dim_) << "D slice (not supported)");
    pmf.reset(new mesh_fem(const_cast<mesh&>(m),1));
    pmf->set_classical_finite_element(1);
    exporting(*pmf);
  }

  void vtk_export::exporting(const mesh_fem& mf) {
    dim_ = mf.linked_mesh().dim();
    if (dim_ > 3) DAL_THROW(dal::failure_error, "attempt to export a " << int(dim_) << "D slice (not supported)");
    if (&mf != pmf.get())
      pmf.reset(new mesh_fem(mf.linked_mesh(),1));
    /* initialize pmf with finite elements suitable for VTK (which only knows
       isoparametric FEMs of order 1 and 2 (with inner points missing)) */
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
      pfem pf = mf.fem_of_element(cv);
      
      bool discontinuous = false;
      for (unsigned i=0; i < pf->nb_dof(cv); ++i) {
        /* could be a better test for discontinuity .. */
        if (!dof_linkable(pf->dof_types()[i])) { discontinuous = true; break; }
      }
      pfem classical_pf1 = discontinuous ? classical_discontinuous_fem(pgt, 1) : classical_fem(pgt, 1);
      int degree = ((pf != classical_pf1 && pf->estimated_degree() > 1) || 
                     pgt->structure() != pgt->basic_structure()) ? 2 : 1;
      pmf->set_finite_element(cv, discontinuous ? 
                              classical_discontinuous_fem(pgt, degree) : 
                              classical_fem(pgt, degree));
    }
    /* find out which dof will be exported to VTK */

    const mesh &m = pmf->linked_mesh();
    pmf_cell_type.resize(pmf->convex_index().last_true() + 1, unsigned(-1));
    pmf_dof_used.sup(0, pmf->nb_dof());
    for (dal::bv_visitor cv(pmf->convex_index()); !cv.finished(); ++cv) {
      int t = -1;
      size_type nbd = pmf->fem_of_element(cv)->nb_dof(cv);
      switch (pmf->fem_of_element(cv)->dim()) {
        case 0: t = VTK_VERTEX; break;
        case 1: 
          if (nbd == 2) t = VTK_LINE;
          else if (nbd == 3) t = VTK_QUADRATIC_EDGE; break;
        case 2: 
          if (nbd == 3) t = VTK_TRIANGLE;
          else if (nbd == 4) t = check_voxel(m.points_of_convex(cv)) ? VTK_PIXEL : VTK_QUAD;
          else if (nbd == 6) t = VTK_QUADRATIC_TRIANGLE;
          else if (nbd == 9) t = VTK_QUADRATIC_QUAD; break;
        case 3:
          if (nbd == 4) t = VTK_TETRA;
          else if (nbd == 10) t = VTK_QUADRATIC_TETRA;
          else if (nbd == 8) t = check_voxel(m.points_of_convex(cv)) ? VTK_VOXEL : VTK_HEXAHEDRON;
          else if (nbd == 27) t = VTK_QUADRATIC_HEXAHEDRON;
          else if (nbd == 6) t = VTK_WEDGE; break;        
      }
      if (t == -1)
        DAL_THROW(dal::failure_error, "semi internal error.. could not map " << 
                  name_of_fem(pmf->fem_of_element(cv)) << " to a VTK cell type");
      pmf_cell_type[cv] = t;
      const std::vector<unsigned> &dmap = getfem_to_vtk_dof_mapping(t);
      //cout << "nbd = " << nbd << ", t = " << t << ", dmap = " << dmap << "\n";
      if (dmap.size() > pmf->nb_dof_of_element(cv)) DAL_INTERNAL_ERROR("inconsistency in vtk_dof_mapping");
      for (unsigned i=0; i < dmap.size(); ++i)
        pmf_dof_used.add(pmf->ind_dof_of_element(cv)[dmap[i]]);
    }
    //cout << "mf.nb_dof = " << mf.nb_dof() << ", pmf->nb_dof=" << pmf->nb_dof() << ", dof_used = " << pmf_dof_used.card() << "\n";
  }


  const stored_mesh_slice& vtk_export::get_exported_slice() const { 
    if (!psl) DAL_THROW(dal::failure_error,"no slice!")
    else return *psl; 
  }

  const mesh_fem& vtk_export::get_exported_mesh_fem() const { 
    if (!pmf.get()) DAL_THROW(dal::failure_error,"no mesh_fem!")
    else return *pmf; 
  }

  void vtk_export::set_header(const std::string& s)
  { strncpy(header, s.c_str(), 256); header[255] = 0; }

  void vtk_export::check_header() {
    if (state >= HEADER_WRITTEN) return;
    os << "# vtk DataFile Version 2.0\n";
    os << header << "\n";
    if (ascii) os << "ASCII\n"; else os << "BINARY\n";
    state = HEADER_WRITTEN;
  }

  void vtk_export::write_separ()
  { if (ascii) os << "\n"; }

  void vtk_export::write_mesh() {
    if (psl) write_mesh_structure_from_slice();
    else write_mesh_structure_from_mesh_fem();
  }

  /* export the slice data as an unstructured mesh composed of simplexes */
  void vtk_export::write_mesh_structure_from_slice() {
    /* element type code for (linear) simplexes of dimensions 0,1,2,3 in VTK */
    static int vtk_simplex_code[4] = { VTK_VERTEX, VTK_LINE, VTK_TRIANGLE, VTK_TETRA }; 
    if (state >= STRUCTURE_WRITTEN) return;
    check_header();
    /* possible improvement: detect structured grids */
    os << "DATASET UNSTRUCTURED_GRID\n";
    os << "POINTS " << psl->nb_points() << " float\n";
    /* 
       points are not merge, vtk is mostly fine with that (except for
       transparency and normals at vertices of 3D elements: all simplex faces
       are considered to be "boundary" faces by vtk)
    */
    for (size_type ic=0; ic < psl->nb_convex(); ++ic) {
      for (size_type i=0; i < psl->nodes(ic).size(); ++i)
	write_vec(psl->nodes(ic)[i].pt.begin());
      write_separ();
    }
    /* count total number of simplexes, and total number of entries
     * in the CELLS section */
    size_type cells_cnt = 0, splx_cnt = 0;
    for (size_type ic=0; ic < psl->nb_convex(); ++ic) {
      for (size_type i=0; i < psl->simplexes(ic).size(); ++i)
	cells_cnt += psl->simplexes(ic)[i].dim() + 2;
      splx_cnt += psl->simplexes(ic).size();
    }
    size_type nodes_cnt = 0;
    write_separ(); os << "CELLS " << splx_cnt << " " << cells_cnt << "\n";
    for (size_type ic=0; ic < psl->nb_convex(); ++ic) {
      const getfem::mesh_slicer::cs_simplexes_ct& s = psl->simplexes(ic);
      for (size_type i=0; i < s.size(); ++i) {
	write_val(int(s[i].dim()+1));
	for (size_type j=0; j < s[i].dim()+1; ++j)
	  write_val(int(s[i].inodes[j] + nodes_cnt));
	write_separ();
      }
      nodes_cnt += psl->nodes(ic).size();
    }
    assert(nodes_cnt == psl->nb_points()); // sanity check
    write_separ(); os << "CELL_TYPES " << splx_cnt << "\n";
    for (size_type ic=0; ic < psl->nb_convex(); ++ic) {
      const getfem::mesh_slicer::cs_simplexes_ct& s = psl->simplexes(ic);
      for (size_type i=0; i < s.size(); ++i) {
	write_val(int(vtk_simplex_code[s[i].dim()]));
      }
      write_separ();
      splx_cnt -= s.size();
    }
    assert(splx_cnt == 0); // sanity check
    state = STRUCTURE_WRITTEN;
  }


  void vtk_export::write_mesh_structure_from_mesh_fem() {
    if (state >= STRUCTURE_WRITTEN) return;
    check_header();
    /* possible improvement: detect structured grids */
    os << "DATASET UNSTRUCTURED_GRID\n";
    os << "POINTS " << pmf_dof_used.card() << " float\n";
    std::vector<int> dofmap(pmf->nb_dof());
    size_type cnt = 0;
    for (dal::bv_visitor d(pmf_dof_used); !d.finished(); ++d) {
      dofmap[d] = cnt++;
      base_node P = pmf->point_of_dof(d);
      write_vec(P.const_begin());
      write_separ();
    }

    size_type nb_cell_values = 0;
    for (dal::bv_visitor cv(pmf->convex_index()); !cv.finished(); ++cv)
      nb_cell_values += (1 + getfem_to_vtk_dof_mapping(pmf_cell_type[cv]).size());

    write_separ(); os << "CELLS " << pmf->convex_index().card() << " " << nb_cell_values << "\n";

    for (dal::bv_visitor cv(pmf->convex_index()); !cv.finished(); ++cv) {
      const std::vector<unsigned> &dmap = getfem_to_vtk_dof_mapping(pmf_cell_type[cv]);
      write_val(int(dmap.size()));
      for (size_type i=0; i < dmap.size(); ++i)
        write_val(int(dofmap[pmf->ind_dof_of_element(cv)[dmap[i]]]));
      write_separ();
    }

    write_separ(); os << "CELL_TYPES " << pmf->convex_index().card() << "\n";
    for (dal::bv_visitor cv(pmf->convex_index()); !cv.finished(); ++cv) {
      write_val(int(pmf_cell_type[cv]));
      write_separ();
    }

    state = STRUCTURE_WRITTEN;
  }

  void vtk_export::write_mesh_quality(const mesh &m) {
    if (psl) {
      mesh_fem mf(const_cast<mesh&>(m),1);
      mf.set_classical_finite_element(0);
      std::vector<scalar_type> q(mf.nb_dof());
      for (size_type d=0; d < mf.nb_dof(); ++d) {
        q[d] = m.convex_quality_estimate(mf.first_convex_of_dof(d));
      }
      write_point_data(mf, q, "convex_quality");
    } else {
      std::vector<scalar_type> q(pmf->convex_index().card());
      for (dal::bv_visitor cv(pmf->convex_index()); !cv.finished(); ++cv) {
        q[cv] = m.convex_quality_estimate(cv);
      }
      write_cell_data(q, "convex_quality");
    }
  }


  /* -------------------------------------------------------------
   * OPENDX export 
   * ------------------------------------------------------------- */

  dx_export::dx_export(std::ostream &os_, bool ascii_)
    : os(os_), ascii(ascii_) { init(); } 
  
  dx_export::dx_export(const std::string& fname, bool ascii_, bool append_)
    : os(real_os), ascii(ascii_) {
    real_os.open(fname.c_str(), 
                 std::ios_base::openmode(std::ios_base::in | std::ios_base::out | 
                                         (append_ ? std::ios_base::ate : std::ios_base::trunc)));
    if (!real_os.good()) DAL_THROW(failure_error, "impossible to write to dx file '"
			    << fname << "'");
    init();
    if (append_) { reread_metadata(); header_written = true; }
  }

  dx_export::~dx_export() { 
    std::ios::pos_type p = os.tellp();
    write_series();
    os << "\n# --end of getfem export\nend\n"; 
    update_metadata(p);
  }

  void dx_export::init() {    
    strcpy(header, "Exported by getfem++");
    psl = 0; dim_ = dim_type(-1); connections_dim = dim_type(-1);  psl_use_merged = false;
    header_written = false;
  }

  void dx_export::write_separ()
  { if (ascii) os << "\n"; }

  template<typename T> static typename std::list<T>::iterator
  get_from_name(std::list<T> &c, 
		const std::string& name, bool raise_error) {
    for (typename std::list<T>::iterator it = c.begin(); 
	 it != c.end(); ++it) {
      if (it->name == name) return it;
    }
    if (raise_error) {
      DAL_THROW(dal::failure_error, 
		"object not found in dx file: " << name);
    } else return c.end();
  }

  std::list<dx_export::dxMesh>::iterator
  dx_export::get_mesh(const std::string& name, bool raise_error) {
    return get_from_name(meshes,name,raise_error);
  }
  std::list<dx_export::dxObject>::iterator
  dx_export::get_object(const std::string& name, bool raise_error) {
    return get_from_name(objects,name,raise_error);
  }


  bool dx_export::new_mesh(std::string &name) {
    name = default_name(name, meshes.size(), "mesh");
    std::list<dxMesh>::iterator it = get_mesh(name, false);
    if (it != meshes.end()) {
      if (&(*it) != &current_mesh())
	std::swap(current_mesh(),*it);
      return false;
    } else {
      meshes.push_back(dxMesh()); meshes.back().name = name;
      return true;
    }
  }

  void dx_export::exporting(const stored_mesh_slice& sl, bool merge_points, std::string name) {
    if (!new_mesh(name)) return;
    psl_use_merged = merge_points;
    if (merge_points) sl.merge_nodes();
    psl = &sl; dim_ = sl.dim();
    if (psl->dim() > 3) DAL_THROW(dal::failure_error, "4D slices and more are not supported");
    for (dim_type d = 0; d <= psl->dim(); ++d) {
      if (psl->nb_simplexes(d)) {
        if (connections_dim == dim_type(-1)) connections_dim = d;
        else DAL_THROW(dal::failure_error, "Cannot export a slice containing simplexes of different dimensions");
      }
    }
    if (connections_dim == dim_type(-1)) 
      DAL_THROW(dal::failure_error, "empty slice!");
  }

 
  void dx_export::exporting(const mesh_fem& mf, std::string name) {
    name = default_name(name, meshes.size(), "mesh");
    if (!new_mesh(name)) return;
    const mesh &m = mf.linked_mesh();
    if (mf.linked_mesh().convex_index().card() == 0) 
      DAL_THROW(dal::failure_error, "won't export an empty mesh");
    
    dim_ = m.dim();
    if (dim_ > 3) DAL_THROW(dal::failure_error, "4D meshes and more are not supported");
    if (&mf != pmf.get())
      pmf.reset(new mesh_fem(const_cast<mesh&>(m),1));
    bgeot::pgeometric_trans pgt = m.trans_of_convex(m.convex_index().first_true());
    if (dxname_of_convex_structure(pgt->structure()->basic_structure()) == 0) 
      DAL_THROW(dal::failure_error, "DX Cannot handle " << 
		bgeot::name_of_geometric_trans(pgt) << ", use slices");
    /* initialize pmf with finite elements suitable for OpenDX */
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt2 = mf.linked_mesh().trans_of_convex(cv);
      if (pgt->structure()->basic_structure() != pgt2->structure()->basic_structure()) {
	DAL_THROW(dal::failure_error, 
		  "Cannot export this mesh to opendx, it contains "
		  "different convex types. Slice it first.");
      }
      pfem pf = mf.fem_of_element(cv);      
      bool discontinuous = false;
      for (unsigned i=0; i < pf->nb_dof(cv); ++i) {
        /* could be a better test for discontinuity .. */
        if (!dof_linkable(pf->dof_types()[i])) { discontinuous = true; break; }
      }
      pfem classical_pf1 = discontinuous ? classical_discontinuous_fem(pgt, 1) : classical_fem(pgt, 1);
      pmf->set_finite_element(cv, classical_pf1);
    }
    pmf_dof_used.add(0, pmf->nb_dof());
    connections_dim = pmf->nb_dof_of_element(m.convex_index().first_true());
  }

  void dx_export::exporting(const mesh& m, std::string name) {
    dim_ = m.dim();
    if (dim_ > 3) DAL_THROW(dal::failure_error, "4D meshes and more are not supported");
    pmf.reset(new mesh_fem(const_cast<mesh&>(m),1));
    pmf->set_classical_finite_element(1);
    exporting(*pmf, name);
  }

  void dx_export::write_series() {
    for (std::list<dxSeries>::const_iterator it = series.begin(); 
	 it != series.end(); ++it) {
      if (it->members.size() == 0) continue;
      size_type count = 0;
      os << "\nobject  \"" << it->name << "\" class series\n";
      for (std::list<std::string>::const_iterator ito = it->members.begin();
	   ito != it->members.end(); ++ito, ++count) {
	os << "  member  " << count << " \"" << (*ito) << "\"\n";
      }
    }
  }

  void dx_export::serie_add_object_(const std::string &serie_name, 
				   const std::string &object_name) {
    std::list<dxSeries>::iterator it = series.begin();
    while (it != series.end() && it->name != serie_name) ++it;
    if (it == series.end()) { 
      series.push_back(dxSeries()); it = series.end(); --it; 
      it->name = serie_name;
    }
    it->members.push_back(object_name);
  }

  void dx_export::serie_add_object(const std::string &serie_name, 
				   const std::string &object_name) {
    /* create a series for edge data if possible (the cost is null 
     and it may be useful) */
    std::list<dxObject>::iterator ito = get_object(object_name, false);
    if (ito != objects.end()) {
      std::list<dxMesh>::iterator itm = get_mesh(ito->mesh);
      if (itm != meshes.end() && (itm->flags & dxMesh::WITH_EDGES)) {
	serie_add_object_(serie_name + "_edges", 
			 object_name + "_edges");
      }
    }
    /* fill the real serie */
    serie_add_object_(serie_name, object_name);
  }

  void dx_export::set_header(const std::string& s)
  { strncpy(header, s.c_str(), 256); header[255] = 0; }

  void dx_export::check_header() {
    if (header_written) return; header_written = true;
    os << "# data file for IBM OpenDX, generated by GetFem++ v " 
       << GETFEM_VERSION << "\n";
    os << "# " << header << "\n";
  }

  void dx_export::update_metadata(std::ios::pos_type pos_series) {
    os.seekp(0,std::ios::end);
    os << "# This file contains the following objects\n";
    std::ios::pos_type pos_end = os.tellp();
    for (std::list<dxSeries>::const_iterator it = series.begin(); 
	 it != series.end(); ++it) {
      os << "#S \"" << it->name << "\" which contains:\n";
      for (std::list<std::string>::const_iterator its = it->members.begin(); 
           its != it->members.end(); ++its)
        os << "#+   \"" << *its << "\"\n";
    }
    for (std::list<dxObject>::const_iterator it = objects.begin(); 
	 it != objects.end(); ++it) {
      os << "#O \"" << it->name << "\" \"" << it->mesh << "\"\n";
    }
    for (std::list<dxMesh>::const_iterator it = meshes.begin(); 
	 it != meshes.end(); ++it) {
      os << "#M \"" << it->name << "\" " << it->flags << "\n";
    }
    os << "#E \"THE_END\" " << std::setw(20) << pos_series << std::setw(20) << pos_end << "\n";
  }

  void dx_export::reread_metadata() {
    char line[512];
    real_os.seekg(0, std::ios::end);
    int count=0; char c;
    unsigned long lu_end, lu_series;
    do { 
      real_os.seekg(-1, std::ios::cur); 
      c = real_os.peek();
    } while (++count < 512 && c != '#');
    real_os.getline(line, sizeof line);
    if (sscanf(line, "#E \"THE_END\" %lu %lu", &lu_series, &lu_end) != 2)
      DAL_THROW(dal::failure_error, "this file was not generated by getfem, "
                "cannot append data to it!\n");
    real_os.seekg(lu_end, std::ios::beg);
    do {
      char name[512]; unsigned n;
      int pos;
      real_os.getline(line, sizeof line);
      if (sscanf(line, "#%c \"%512[^\"]\"%n", &c, name, &pos) < 1) 
        DAL_THROW(dal::failure_error, "corrupted file! your .dx file is broken\n");
      if (c == 'S') {
        series.push_back(dxSeries()); series.back().name = name;
      } else if (c == '+') {
        series.back().members.push_back(name);
      } else if (c == 'O') {
        objects.push_back(dxObject()); objects.back().name = name; 
        sscanf(line+pos, " \"%512[^\"]\"", name); objects.back().mesh = name;
      } else if (c == 'M') {
        meshes.push_back(dxMesh()); meshes.back().name = name;
        sscanf(line+pos, "%u", &n); meshes.back().flags = n;
      } else if (c == 'E') {
        break;
      } else DAL_THROW(dal::failure_error, "corrupted file! your .dx file is broken\n");
    } while (1);
    real_os.seekp(lu_series, std::ios::beg);
  }
  
  void dx_export::write_convex_attributes(bgeot::pconvex_structure cvs) {
    const char *s_elem_type = dxname_of_convex_structure(cvs);
    if (!s_elem_type) 
      DAL_WARNING1("OpenDX won't handle this kind of convexes");
    os << "\n  attribute \"element type\" string \"" << s_elem_type << "\"\n"
       << "  attribute \"ref\" string \"positions\"\n\n";
  }

  const char *dx_export::dxname_of_convex_structure(bgeot::pconvex_structure cvs) {
    const char *s_elem_type = 0;
    switch (cvs->dim()) {
      /* TODO: do something for point data */
      case 1: s_elem_type = "lines"; break;
      case 2: 
	if (cvs->nb_points() == 3) 
	  s_elem_type = "triangles"; 
	else if (cvs->nb_points() == 4) 
	  s_elem_type = "quads"; 
	break;
      case 3: 
	if (cvs->nb_points() == 4)
	  s_elem_type = "tetrahedra"; 
	else if (cvs->nb_points() == 8)
	  s_elem_type = "cubes";
	break;
    }
    return s_elem_type;
  }

  void dx_export::write_mesh() {
    check_header();
    if (current_mesh().flags & dxMesh::STRUCTURE_WRITTEN) return;
    if (psl) write_mesh_structure_from_slice();
    else write_mesh_structure_from_mesh_fem();

    os << "\nobject \"" << current_mesh_name() << "\" class field\n"
       << "  component \"positions\" value \"" 
       << name_of_pts_array(current_mesh_name()) << "\"\n"
       << "  component \"connections\" value \"" 
       << name_of_conn_array(current_mesh_name()) << "\"\n";
    current_mesh().flags |= dxMesh::STRUCTURE_WRITTEN;
  }

  /* export the slice data as an unstructured mesh composed of simplexes */
  void dx_export::write_mesh_structure_from_slice() {
    os << "\nobject \"" << name_of_pts_array(current_mesh_name()) 
       << "\" class array type float rank 1 shape " 
      << int(psl->dim()) 
       << " items " << (psl_use_merged ? psl->nb_merged_nodes() : psl->nb_points());
    if (!ascii) os << " " << endianness() << " binary";
    os << " data follows\n";
    if (psl_use_merged) {
      for (size_type i=0; i < psl->nb_merged_nodes(); ++i) {
	for (size_type k=0; k < psl->dim(); ++k) 
	  write_val(float(psl->merged_point(i)[k]));
	write_separ();
      }
    } else {
      for (size_type ic=0; ic < psl->nb_convex(); ++ic) {
	for (size_type i=0; i < psl->nodes(ic).size(); ++i)
	  for (size_type k=0; k < psl->dim(); ++k)
	    write_val(float(psl->nodes(ic)[i].pt[k]));
	write_separ();
      }    
    }

    os << "\nobject \"" << name_of_conn_array(current_mesh_name()) 
       << "\" class array type int rank 1 shape " 
       << int(connections_dim+1)
       << " items " << psl->nb_simplexes(connections_dim);
    if (!ascii) os << " " << endianness() << " binary";
    os << " data follows\n";

    size_type nodes_cnt = 0; /* <- a virer , global_index le remplace */
    for (size_type ic=0; ic < psl->nb_convex(); ++ic) {
      const getfem::mesh_slicer::cs_simplexes_ct& s = psl->simplexes(ic);
      for (size_type i=0; i < s.size(); ++i) {
        if (s[i].dim() == connections_dim) {
          for (size_type j=0; j < s[i].dim()+1; ++j) {
	    size_type k;
	    if (psl_use_merged)
	      k = psl->merged_index(ic, s[i].inodes[j]);
	    else k = psl->global_index(ic, s[i].inodes[j]);
            write_val(int(k));
	  }
          write_separ();
        }
      }
      nodes_cnt += psl->nodes(ic).size();
    }

    write_convex_attributes(bgeot::simplex_structure(connections_dim));
    assert(nodes_cnt == psl->nb_points()); // sanity check
  }



  void dx_export::write_mesh_structure_from_mesh_fem() {
    os << "\nobject \"" << name_of_pts_array(current_mesh_name()) 
       << "\" class array type float rank 1 shape " 
       << int(pmf->linked_mesh().dim())
       << " items " << pmf->nb_dof();
    if (!ascii) os << " " << endianness() << " binary";
    os << " data follows\n";

    /* possible improvement: detect structured grids */
    for (size_type d = 0; d < pmf->nb_dof(); ++d) {
      const base_node P = pmf->point_of_dof(d);
      for (size_type k=0; k < dim_; ++k)
	write_val(float(P[k]));
      write_separ();
    }

    os << "\nobject \"" << name_of_conn_array(current_mesh_name()) 
       << "\" class array type int rank 1 shape " 
       << int(connections_dim)
       << " items " << pmf->convex_index().card();
    if (!ascii) os << " " << endianness() << " binary";
    os << " data follows\n";

    for (dal::bv_visitor cv(pmf->convex_index()); !cv.finished(); ++cv) {
      for (size_type i=0; i < connections_dim; ++i)
        write_val(int(pmf->ind_dof_of_element(cv)[i]));
      write_separ();
    }
    write_convex_attributes(pmf->linked_mesh().structure_of_convex(pmf->convex_index().first_true())->basic_structure());
  }

  void dx_export::exporting_mesh_edges(bool with_slice_edges) {
    write_mesh();
    if (current_mesh().flags & dxMesh::WITH_EDGES) return;
    if (psl) write_mesh_edges_from_slice(with_slice_edges);
    else write_mesh_edges_from_mesh();
    current_mesh().flags |= dxMesh::WITH_EDGES;
    os << "\nobject \"" << name_of_edges_array(current_mesh_name()) << "\" class field\n"
       << "  component \"positions\" value \"" 
       << name_of_pts_array(current_mesh_name()) << "\"\n"
       << "  component \"connections\" value \"" 
       << name_of_conn_array(name_of_edges_array(current_mesh_name())) << "\"\n";
  }

  void dx_export::write_mesh_edges_from_slice(bool with_slice_edges) {
    std::vector<size_type> edges;
    dal::bit_vector slice_edges;
    psl->get_edges(edges, slice_edges, psl_use_merged);
    if (with_slice_edges) slice_edges.clear();
    os << "\nobject \"" << name_of_conn_array(name_of_edges_array(current_mesh_name())) 
       << "\" class array type int rank 1 shape 2"
       << " items " << edges.size()/2 - slice_edges.card();
    if (!ascii) os << " " << endianness() << " binary";
    os << " data follows\n";
    for (size_type i=0; i < edges.size()/2; ++i) {
      if (!slice_edges.is_in(i)) {
	write_val(int(edges[2*i]));
	write_val(int(edges[2*i+1]));
      }
      if ((i+1)%10 == 0) write_separ();
    }
    write_separ();
    write_convex_attributes(bgeot::simplex_structure(1));
  }

  void dx_export::write_mesh_edges_from_mesh() {
    bgeot::mesh_structure ms(pmf->linked_mesh()); ms.to_edges();
    os << "\nobject \"" << name_of_conn_array(name_of_edges_array(current_mesh_name())) 
       << "\" class array type int rank 1 shape 2"
       << " items " << ms.convex_index().card();
    if (!ascii) os << " " << endianness() << " binary";
    os << " data follows\n";
    for (dal::bv_visitor cv(ms.convex_index()); !cv.finished(); ++cv) {
      write_val(int(ms.ind_points_of_convex(cv)[0]));
      write_val(int(ms.ind_points_of_convex(cv)[1]));
      if ((cv+1)%20 == 0) write_separ();
    }
    write_separ();
    write_convex_attributes(bgeot::simplex_structure(1));
  }
  
}  /* end of namespace getfem.                                             */
