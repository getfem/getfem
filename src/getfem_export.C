/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_export.C : export and import data.                    */
/*     									   */
/* Date : October 15, 2001.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001-2002  Yves Renard.                                   */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */


#include <getfem_export.h>

namespace getfem
{
  void classical_mesh_fem(mesh_fem& mf, short_type K)
  {
    for (dal::bv_visitor cv(mf.linked_mesh().convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
      mf.set_finite_element(cv, classical_fem(pgt,K), exact_classical_im(pgt));
    }
  }

  vtk_export::vtk_export(std::ostream &os_, bool ascii_) : os(os_), ascii(ascii_) {
    init();
  } 
  
  vtk_export::vtk_export(const std::string& fname, bool ascii_) : os(real_os), ascii(ascii_), 
								  real_os(fname.c_str()) {
    if (!real_os)
      DAL_THROW(failure_error, "impossible to write to vtk file '" << fname << "'");
    init();
  }

  void vtk_export::init() {
    static int test_endian = 0x01234567;
    strcpy(header, "Exported by getfem++");
    header_done = false; mesh_structure_done = false; psl = 0;
    if (*((char*)&test_endian) == 0x67)
      reverse_endian = true;
    else reverse_endian = false;
    //cout << "reverse_endian = " << reverse_endian << "\n";
  }

  void vtk_export::check_mesh(const getfem_mesh &m) {
    if (&psl->linked_mesh() != &m) {
      DAL_THROW(dal::failure_error, "the meshes are different.");
    }
  }

  void vtk_export::set_mesh(const getfem_mesh &m, unsigned nrefine) {
    if (psl == 0) {
      owned_psl.reset(new stored_mesh_slice());
      psl = owned_psl.get();
      owned_psl.get()->build(m, slicer_none(), nrefine);
    } else check_mesh(m);
  }
  
  void vtk_export::set_slice(const stored_mesh_slice& sl) {
    if (psl) { DAL_THROW(dal::failure_error, "a slice was already set."); }
    else psl = &sl;
  }

  const stored_mesh_slice& vtk_export::get_slice() const { 
    if (!psl) DAL_THROW(dal::failure_error,"no slice!")
    else return *psl; 
  }

  void vtk_export::set_header(const std::string& s) {
    strncpy(header, s.c_str(), 256); header[255] = 0;
  }

  void vtk_export::check_header() {
    if (header_done) return;
    os << "# vtk DataFile Version 2.0\n";
    os << header << "\n";
    if (ascii) os << "ASCII\n"; else os << "BINARY\n";
    header_done = true;
  }

  void vtk_export::write_separ() {
    if (ascii) os << "\n";
  }

  /* export the slice data as an unstructured mesh composed of simplexes */
  void vtk_export::write_mesh_structure() {
    static int vtk_simplex_code[4] = { 1, 3, 5, 10 }; /* element type code for simplexes of dimensions 0,1,2,3 in VTK */
    if (mesh_structure_done) return;
    check_header();
    if (psl == 0) 
      DAL_THROW(dal::failure_error, "cannot export the mesh data: you did not give a mesh!\n");
    /* possible improvement: detect structured grids */
    os << "DATASET UNSTRUCTURED_GRID\n";
    os << "POINTS " << psl->nb_points() << " float\n";
    /* no need to merge the points, vtk is fine with that */
    for (size_type ic=0; ic < psl->nb_convex(); ++ic) {
      for (size_type i=0; i < psl->nodes(ic).size(); ++i)
	write_vec(psl->nodes(ic)[i].pt.begin());
      write_separ();
    }
    /* count total number of simplexes, and total number of entries in the CELLS section */
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
    write_separ(); os << "POINT_DATA " << psl->nb_points() << "\n";
    mesh_structure_done = true;
  }

  void vtk_export::write_mesh_quality(const getfem_mesh &m, unsigned nrefine) {
    set_mesh(m,nrefine);
    mesh_fem mf(const_cast<getfem_mesh&>(m),1); mf.set_classical_finite_element(0);
    std::vector<scalar_type> q(mf.nb_dof());
    for (size_type d=0; d < mf.nb_dof(); ++d) {
      q[d] = m.convex_quality_estimate(mf.first_convex_of_dof(d));
    }
    std::vector<scalar_type> Uslice(psl->nb_points());
    psl->interpolate(mf, q, Uslice);
    write_dataset(Uslice, "convex_quality");
  }
}  /* end of namespace getfem.                                             */


