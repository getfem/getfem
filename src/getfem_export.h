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


  /* this function has nothing to do here .. */
  void classical_mesh_fem(mesh_fem& mf, short_type K) IS_DEPRECATED;

  /** 
      export class to VTK ( http://www.kitware.com/vtk.html ) file format 
      (not the XML format, but the old format)

      A vtk_export can store multiple scalar/vector fields.
  */
  class vtk_export {
    std::ostream &os;
    char header[256]; // hard limit in vtk
    bool ascii;
    bool header_done, mesh_structure_done;
    const stored_mesh_slice *psl;
    std::auto_ptr<stored_mesh_slice> owned_psl;
    std::ofstream real_os;
    bool reverse_endian;
  public:
    vtk_export(const std::string& fname, bool ascii_ = false);
    vtk_export(std::ostream &os_, bool ascii_ = false);
    /** choose the mesh that will be exported -- if you are juste calling 
       write_dataset(mf), you are not required to call this function */
    void set_mesh(const getfem_mesh &m, unsigned nrefine=1);
    /** if you do want to export a mesh slice instead of a whole mesh,
       use set_slice instead of set_mesh */
    void set_slice(const stored_mesh_slice& sl);
    /** the header is the second line of text in the exported file,
       you can put whatever you want */
    void set_header(const std::string& s);
    /** append a new scalar or vector field defined on mf to the .vtk file. 
        NO SPACE ALLOWED in 'name' */
    template<class VECT> void write_dataset(const getfem::mesh_fem &mf, const VECT& U0, const std::string& name);
    /** append a new scalar or vector field to .vtk file. The Uslice
     vector is the field interpolated on the exported mesh_slice 
     NO SPACE ALLOWED in 'name' */
    template<class VECT> void write_dataset(const VECT& U, const std::string& name);
    /** export a data_set correspounding to measures of quality for each convex
        of the mesh */
    void write_mesh_quality(const getfem_mesh &m, unsigned nrefine=1);
    const stored_mesh_slice& get_slice() const;
  private:
    void init();
    void check_header();
    void write_mesh_structure();
    void check_mesh(const getfem_mesh &m);
    template<class T> void write_val(T v);
    template<class V> void write_vec(V p);
    void write_separ();
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
    for (size_type i=0; i < psl->dim(); ++i) {
      v[i] = p[i];
    }
    for (size_type i=psl->dim(); i < 3; ++i) v[i] = 0.0f;
    write_val(v[0]);write_val(v[1]);write_val(v[2]);
  }

  template<class VECT>
  void vtk_export::write_dataset(const getfem::mesh_fem &mf, const VECT& U,
				 const std::string& name) {
    set_mesh(mf.linked_mesh());
    size_type Q = (U.size() / mf.nb_dof())*mf.get_qdim();
    std::vector<scalar_type> Uslice(Q*psl->nb_points());
    psl->interpolate(mf, U, Uslice);
    write_dataset(Uslice,name);
  }
  
  template<class VECT>
  void vtk_export::write_dataset(const VECT& Uslice, const std::string& name) {
    write_mesh_structure();
    size_type Q = Uslice.size() / psl->nb_points();
    write_separ();
    if (Q == 1) {
      os << "SCALARS " << name << " float\n";
      os << "LOOKUP_TABLE default\n";
      for (size_type i=0; i < psl->nb_points(); ++i) {
	write_val(float(Uslice[i]));
      }
    } else if (Q <= 3) {
      os << "VECTORS " << name << " float\n";
      for (size_type i=0; i < psl->nb_points(); ++i) {
	write_vec(Uslice.begin() + i*Q);
      }
    } else DAL_THROW(dal::dimension_error,
		     "vtk does not accept vectors of dimension > 3");
    write_separ();
  }
}  /* end of namespace getfem.                                             */


#endif /* GETFEM_EXPORT_H__  */
