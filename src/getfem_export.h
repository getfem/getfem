/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_export.h : export and import data.                    */
/*     									   */
/* Date : October 15, 2001.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*          Julien Pommier, pommier@gmm.insa-tlse.fr                       */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001-2004  Yves Renard, Julien Pommier.                   */
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


#ifndef GETFEM_EXPORT_H__
#define GETFEM_EXPORT_H__

#include <getfem_interpolation.h>

namespace getfem {

  /* ********************************************************************* */
  /*                                                                       */
  /*  Save a solution in a file with a Pk interpolation.                   */
  /*                                                                       */
  /* ********************************************************************* */

  /** Semi-obsoleted function: In general, it is preferable to avoid this
      function, and directly save the mesh_fem (mf.write_to_file) in one file,
      and the fields in another file. */
  template<typename VECT>
    void save_solution(const std::string &filename, const mesh_fem &mf,
		       const VECT &U, short_type K)
  { // a corriger
    dim_type N = mf.linked_mesh().dim();
    dim_type P = mf.get_qdim();

    std::ofstream o(filename.c_str());
    if (!o) DAL_THROW(internal_error, "impossible to open file");
    base_node pt1(N), pt3(P);
    base_matrix G;

    o << "% GETFEM++ DATA FILE\nBEGIN DATA ELEMENT\nN = " << int(N)
      << "\nP = " << int(P) << "\nK = " << K << endl << endl;
    o.precision(14);
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
      pfem pfe = classical_fem(pgt, K);
      pfem pf1 = mf.fem_of_element(cv);
      // size_type nbd1 = mf.nb_dof_of_element(cv);
      size_type nbd2 = pfe->nb_dof();

      o << "DIM = " << int(pgt->dim()) << endl;

      if (pf1->need_G()) 
	bgeot::vectors_to_base_matrix(G,mf.linked_mesh().points_of_convex(cv));

      if (pf1->target_dim() != 1)
	DAL_THROW(to_be_done_error, "to be done ... ");
      fem_interpolation_context ctx(pgt,pf1,base_node(),G,cv);
      for (size_type i = 0; i < nbd2; ++i)
      {
	/* point corresponding to the i th node.                           */
	pt1 = pgt->transform(pfe->node_of_dof(i),
			     mf.linked_mesh().points_of_convex(cv));
	for (size_type j = 0; j < N; ++j) o << pt1[j] << " ";
	o << "  ";

	/* interpolation of the solution.                                  */
	/* faux dans le cas des éléments vectoriel.                        */
	ctx.set_xref(pfe->node_of_dof(i));
	pf1->interpolation(ctx, 
			   gmm::sub_vector(U,
	   gmm::sub_index(mf.ind_dof_of_element(cv))), pt3, mf.get_qdim());
	for (size_type j = 0; j < P; ++j) o << pt3[j] << " ";
	o << endl;
      }
      o << endl;
    }
    o << "END DATA ELEMENT" << endl;
    o.close();
  }

  /** Semi-obsoleted function */
  template<typename VECT>
    void load_solution(const std::string &fi, getfem_mesh &mesh,
		       mesh_fem &mef, VECT &U, short_type &K) {
    std::ifstream ist(fi.c_str());
    size_type token, DIM, Np, N;
    char tmp[100], c;
    bool error = false, began = false;
    dal::dynamic_array<base_node> ptab;
    dal::dynamic_array<scalar_type> vtab;
    size_type nbvtab = 0;
    dim_type P;
    
    if (!ist) DAL_THROW(file_not_found_error, "File " << fi << " not found");

    N = 1; P = 1; K = 1; DIM = 1; // Default parameters
    
    if (!(ftool::read_until(ist, "BEGIN DATA ELEMENT"))) error = true;

    while (!error) {
      if (ist.eof()) { error = true; break; }
      ftool::get_token(ist, tmp, 99); token = 0;
      if (!strcmp(tmp, "N")) token = 1;
      if (!strcmp(tmp, "P")) token = 2;
      if (!strcmp(tmp, "K")) token = 3;
      if (!strcmp(tmp, "DIM")) token = 4;
      if (!strcmp(tmp, "END")) break;
      if (token) {
	if (token >= 1 && token <= 2 && began)
	  cerr << "WARNING : Ignoring principal dimensions redefinition"<<endl;
	else {
	  ftool::get_token(ist, tmp, 99);
	  if (strcmp(tmp, "=")) error = true;
	  ftool::get_token(ist, tmp, 99);
	  int k = atoi(tmp); if (k < 0 || k > 254) error = true;
	  // cout << "token " << token << " val = " << k << endl;
	  switch (token) {
	  case 1 : DIM = N = k; if (k == 0) error = true; break;
	  case 2 : P = k; mef.set_qdim(P); if (k == 0) error = true; break;
	  case 3 : K = k; break;
	  case 4 : DIM = k; break;
	  }
	}
      }
      else {
	began = true;
	// cout << "tmp =" << tmp << endl;
	ist.putback(' ');
	for(int k = strlen(tmp)-1; k >= 0; --k) ist.putback(tmp[k]);
	size_type nbpt = 0;
	// cout << "N = " << N << " P = " << int(P) << " K = "
	//      << K << " DIM = " << DIM << endl;
	while(nbpt < 1000) {
	  if (ptab[nbpt].size() != N) ptab[nbpt].resize(N);
	  for (dim_type i = 0; i < N; ++i) ist >> (ptab[nbpt])[i];
	  // cout << "Ptab[" << nbpt << "] = " << ptab[nbpt] << endl;
	  for (short_type i = 0; i < P; ++i) ist >> vtab[nbvtab * P + i];
	  nbpt++; nbvtab++;
	  while(ist.get(c))
	    { if (!(isspace(c)) || c == 10) { ist.putback(c); break; } }
	  ist.get(c); if (c != 10 && !(ist.eof())) { error = true; break; }
	  while(ist.get(c))
	    { if (!(isspace(c)) || c == 10) { ist.putback(c); break; } }
	  ist.get(c); if (c == 10 || ist.eof()) break; else ist.putback(c);
	}
	if (nbpt >= 1000) { error = true; break; }
	if (nbpt == 0) break;

	Np = 1; for (dim_type i = 0; i < DIM; ++i) Np *= K + 1;
	
	if (nbpt == bgeot::alpha(DIM, K)) {
	  std::vector<getfem::base_node> pnode(DIM+1);
	  for (dim_type i = 0; i <= DIM; ++i) 
	    pnode[i] = ptab[bgeot::alpha(i, K) - 1];
	  size_type i = mesh.add_simplex_by_points(DIM, pnode.begin());
	  mef.set_finite_element(i, getfem::PK_fem(DIM, K),
				 getfem::exact_simplex_im(DIM));
	}
	else if (nbpt == Np) {
	  std::vector<getfem::base_node> pnode(1 << DIM);
	  size_type j = 0;
	  for (size_type i = 0; i < (size_type(1) << DIM); ++i) {
	    pnode[i] = ptab[j];
	    size_type k = i + 1, l = K-1;
	    while (!(k & 1)) { l *= K+1; k >>= 1; }
	    j += l + 1;
	  }
	  j = mesh.add_parallelepiped_by_points(DIM, pnode.begin());
	  mef.set_finite_element(j, getfem::QK_fem(DIM, K),
				 getfem::exact_parallelepiped_im(DIM));
	}
	else if (nbpt == bgeot::alpha(DIM - 1, K) * bgeot::alpha(1, K)) {
	  std::vector<getfem::base_node> pnode(2*DIM);
	  for (dim_type i = 0; i < DIM; ++i)
	    pnode[i] = ptab[bgeot::alpha(i, K) - 1];
	  for (dim_type i = 0; i < DIM; ++i)
	    pnode[i+DIM] = ptab[bgeot::alpha(i, K)-1 + bgeot::alpha(DIM-1, K)];
	  int i = mesh.add_prism_by_points(DIM, pnode.begin());
	  mef.set_finite_element(i, getfem::PK_prism_fem(DIM, K),
				 getfem::exact_prism_im(DIM));
	}
	else 
	  DAL_THROW(failure_error, "Unknown element in file " << fi);
      }
    }
    
    ist.close();
    if (error) DAL_THROW(failure_error, "Format error in file " << fi);
    
    // Data repartition
    // à corriger
    
    U.resize(mef.nb_dof());
    std::fill(U.begin(), U.end(), 0.0);
    std::vector<int> cp(mef.nb_dof() / P);
    std::fill(cp.begin(), cp.end(), 0);
    size_type l = 0;
    for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i) {
      size_type nbd = mef.nb_dof_of_element(i);
      for (size_type j = 0; j < nbd / P; ++j, ++l) {
	size_type dof = mef.ind_dof_of_element(i)[j*P];
	for (short_type k = 0; k < P; ++k) U[dof +k] += vtab[l*P+k];
	(cp[dof])++;
      }
    }
    for (size_type i = 0; i < mef.nb_dof() / P; ++i) {
      if (cp[i] == 0) DAL_THROW(internal_error, "Internal error");
      for (short_type k = 0; k < P; ++k)
	U[i*P +k] /= getfem::scalar_type(cp[i]);
    }
    
  }


  /** 
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
    std::auto_ptr<mesh_fem> pmf;
    dal::bit_vector pmf_dof_used;
    std::vector<unsigned> pmf_cell_type;
    std::ofstream real_os;
    dim_type dim_;
    bool reverse_endian;
    enum { EMPTY, HEADER_WRITTEN, STRUCTURE_WRITTEN, IN_CELL_DATA, IN_POINT_DATA } state;
  public:
    typedef enum { VTK_VERTEX = 1, VTK_LINE = 3, VTK_QUADRATIC_EDGE = 21, 
                   VTK_TRIANGLE = 5, VTK_QUADRATIC_TRIANGLE = 22, 
                   VTK_PIXEL = 8, VTK_QUAD = 9, VTK_QUADRATIC_QUAD = 23,
                   VTK_TETRA = 10, VTK_QUADRATIC_TETRA = 24, 
                   VTK_WEDGE = 13, /*VTK_QUADRATIC_WEDGE = 26,*/
                   VTK_VOXEL = 11, VTK_HEXAHEDRON = 12, VTK_QUADRATIC_HEXAHEDRON = 25 } vtk_cell_type;
    vtk_export(const std::string& fname, bool ascii_ = false);
    vtk_export(std::ostream &os_, bool ascii_ = false);

    
    void exporting(const getfem_mesh& m);
    void exporting(const mesh_fem& mf);
    void exporting(const stored_mesh_slice& sl);

    /** the header is the second line of text in the exported file,
       you can put whatever you want -- call this before any write_dataset or write_mesh */
    void set_header(const std::string& s);
    void write_mesh();
    /** append a new scalar or vector field defined on mf to the .vtk file.  If
        you are exporting a slice, or if mf != get_exported_mesh_fem(), U will
        be interpolated on the slice, or on get_exported_mesh_fem().

        NO SPACE ALLOWED in 'name' */
    template<class VECT> void write_point_data(const getfem::mesh_fem &mf, const VECT& U0, const std::string& name);
    /** append a new scalar or vector field to .vtk file. The Uslice vector is
        the field interpolated on the exported mesh_slice This function should
        not be used if you are not exporting a slice!  NO SPACE ALLOWED in
        'name' */
    template<class VECT> void write_sliced_point_data(const VECT& Uslice, const std::string& name);
    /** export data which is constant over each element. You should not use
        this function if you are exporting a slice.  U should have
        convex_index().card() elements.  */
    
    template<class VECT> void write_cell_data(const VECT& U, const std::string& name); 
    /** export a data_set correspounding to measures of quality for each convex
        of the supplied mesh (which should have the same number of convex than
        the one used in the vtk_export)

        If a slice is being exported, the convex quality is written as
        point_data (TO IMPROVE ONEDAY), if a mesh/mesh_fem is being exported,
        it is written as cell_data
    */
    void write_mesh_quality(const getfem_mesh &m);
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
    template<class V> void write_vec(V p);
    template<class IT> void write_3x3tensor(IT p);
    void write_separ();
    template<class VECT> void write_dataset_(const VECT& U, const std::string& name, bool cell_data=false); 
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

  template<class IT> void vtk_export::write_vec(IT p) {
    float v[3];
    for (size_type i=0; i < dim_; ++i) {
      v[i] = p[i];
    }
    for (size_type i=dim_; i < 3; ++i) v[i] = 0.0f;
    write_val(v[0]);write_val(v[1]);write_val(v[2]);
  }

  template<class IT> void vtk_export::write_3x3tensor(IT p) {
    float v[3][3];
    memset(v, 0, sizeof v);
    for (size_type i=0; i < dim_; ++i) {
      for (size_type j=0; j < dim_; ++j)
        v[i][j] = p[i + j*dim_];
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
    size_type Q = (gmm::vect_size(U) / mf.nb_dof())*mf.get_qdim();
    if (psl) {
      std::vector<scalar_type> Uslice(Q*psl->nb_points());
      psl->interpolate(mf, U, Uslice);
      write_dataset_(Uslice,name);
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

  template<class VECT>
  void vtk_export::write_cell_data(const VECT& U, const std::string& name) {
    write_dataset_(U, name, true);
  }

  template<class VECT>
  void vtk_export::write_sliced_point_data(const VECT& U, const std::string& name) {
    write_dataset_(U, name, false);
  }

  template<class VECT>
  void vtk_export::write_dataset_(const VECT& U, const std::string& name, bool cell_data) {
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
    size_type Q = gmm::vect_size(U) / nb_val;
    if (gmm::vect_size(U) != nb_val*Q) 
      DAL_THROW(dal::failure_error, "inconsistency in the size of the dataset: " 
                << gmm::vect_size(U) << " != " << nb_val << "*" << Q);
    write_separ();
    if (Q == 1) {
      os << "SCALARS " << name << " float 1\n";
      os << "LOOKUP_TABLE default\n";
      for (size_type i=0; i < nb_val; ++i) {
	write_val(float(U[i]));
      }
    } else if (Q <= 3) {
      os << "VECTORS " << name << " float\n";
      for (size_type i=0; i < nb_val; ++i) {
	write_vec(U.begin() + i*Q);
      }
    } else if (Q == dal::sqr(dim_)) {
      /* tensors : coef are supposed to be stored in FORTRAN order 
         in the VTK file, they are written with C (row major) order
       */
      os << "TENSORS " << name << " float\n";
      for (size_type i=0; i < nb_val; ++i) {
        write_3x3tensor(U.begin() + i*Q);
      }
    } else DAL_THROW(dal::dimension_error,
		     "vtk does not accept vectors of dimension > 3");
    write_separ();
  }

  class dx_export {
    std::ostream &os;
    char header[256]; 
    bool ascii;
    const stored_mesh_slice *psl;
    std::auto_ptr<mesh_fem> pmf;
    dal::bit_vector pmf_dof_used;
    std::vector<unsigned> pmf_cell_type;
    std::fstream real_os;
    dim_type dim_, connections_dim;
    std::string series_name;
    struct log_names {
      unsigned count;
      std::list<std::string> names;
      log_names() : count(0) {}
      void add(const std::string &s) { ++count; names.push_back(s); }
    };
    enum { MESHES, SCALAR_FIELDS, VECTOR_FIELDS, TENSOR_FIELDS, SERIES_OBJECTS, NB_LOG };
    std::vector<log_names> log;
    enum { EMPTY, HEADER_WRITTEN, STRUCTURE_WRITTEN } state;
  public:
    dx_export(const std::string& fname, bool ascii_ = false, bool append_ = false);
    dx_export(std::ostream &os_, bool ascii_ = false);
    ~dx_export();
    void exporting(const getfem_mesh& m, const char *name = 0);
    void exporting(const mesh_fem& mf, const char *name = 0);
    void exporting(const stored_mesh_slice& sl, const char *name = 0);

    /** the header is the second line of text in the exported file,
       you can put whatever you want -- call this before any write_dataset or write_mesh */
    void set_header(const std::string& s);
    void write_mesh();
    void begin_series(const char *name = 0);
    void end_series();
    template<class VECT> void write_point_data(const getfem::mesh_fem &mf, const VECT& U0, const std::string& name);
    template<class VECT> void write_sliced_point_data(const VECT& Uslice, const std::string& name);
    template<class VECT> void write_cell_data(const VECT& U, const std::string& name); 
    void write_mesh_quality(const getfem_mesh &m);
    void write_normals();
    const stored_mesh_slice& get_exported_slice() const;
    const mesh_fem& get_exported_mesh_fem() const;

  private:
    void init();
    void write_separ();
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
    std::string name_of_pts_array() { return log[MESHES].names.back() + std::string("_pts"); }
    std::string name_of_conn_array() { return log[MESHES].names.back() + std::string("_conn"); }
    void check_header();
    const char *dxname_of_convex_structure(bgeot::pconvex_structure cvs);
    void write_convex_attributes(bgeot::pconvex_structure cvs);
    void write_mesh_structure_from_slice();
    void write_mesh_structure_from_mesh_fem();
    template<class VECT> void write_dataset_(const VECT& U, const std::string& name, bool cell_data=false);
  };


  template<class VECT>
  void dx_export::write_point_data(const getfem::mesh_fem &mf, const VECT& U,
                                   const std::string& name) {
    size_type Q = (gmm::vect_size(U) / mf.nb_dof())*mf.get_qdim();
    if (psl) {
      std::vector<scalar_type> Uslice(Q*psl->nb_points());
      psl->interpolate(mf, U, Uslice);
      write_dataset_(Uslice,name);
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

  template<class VECT>
  void dx_export::write_sliced_point_data(const VECT& U, const std::string& name) {
    write_dataset_(U, name, false);
  }

  template<class VECT> void dx_export::write_dataset_(const VECT& U, const std::string& name, bool cell_data) {
    write_mesh();
    size_type nb_val = 0;
    log[SERIES_OBJECTS].add(name);
    if (cell_data) {
      nb_val = psl ? psl->linked_mesh().convex_index().card() 
        : pmf->linked_mesh().convex_index().card(); 
    } else {
      nb_val = psl ? psl->nb_points() : pmf_dof_used.card();
    }
    size_type Q = gmm::vect_size(U) / nb_val;
    if (gmm::vect_size(U) != nb_val*Q) 
      DAL_THROW(dal::failure_error, "inconsistency in the size of the dataset: " 
                << gmm::vect_size(U) << " != " << nb_val << "*" << Q);

    os << "object \"" << name << "_data\" class array type float rank ";
    if (Q == 1) { os << "0"; log[SCALAR_FIELDS].add(name); } /* scalar data */
    else if (Q == 4) { os << "2 shape 2 2"; log[TENSOR_FIELDS].add(name); } /* or 2x2 tensor data */
    else if (Q == 9) { os << "2 shape 3 3"; log[TENSOR_FIELDS].add(name); } /* or 2x2 tensor data */
    else { os << "1 shape " << Q; log[VECTOR_FIELDS].add(name); } /* fallback: vector data */
    os << " items " << nb_val;
    if (!ascii) os << " " << endianness() << " binary";
    os << " data follows" << endl;
    for (size_type i=0; i < nb_val*Q; ++i) { 
      write_val(float(U[i])); 
      if (((i+1) % (Q > 1 ? Q : 10)) == 0) write_separ();
    }
    write_separ();

    if (!cell_data)
      os << "attribute \"dep\" string \"positions\"\n";
    else os << "attribute \"dep\" string \"connections\"\n";
    os << "\n";
    /* write footer */
    os << "object \"" << name << "\" class field\n"
       << "component \"positions\" value \"" << name_of_pts_array() << "\"\n"
       << "component \"connections\" value \"" << name_of_conn_array() << "\"\n"
       << "component \"data\" value \"" << name << "_data\"\n";
  }

}  /* end of namespace getfem. */


#endif /* GETFEM_EXPORT_H__  */
