/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_export.h : export and import data.                    */
/*     									   */
/* Date : 2004/01/12.                                                      */
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

#include <getfem_mesh_fem.h>
//#include <getfem_mat_elem.h>
#include <bgeot_geotrans_inv.h>
#include <dal_tree_sorted.h>
#include <getfem_mesh_slice.h>

namespace getfem
{

  /* ********************************************************************* */
  /*                                                                       */
  /*  Save a solution in a file with a Pk interpolation.                   */
  /*                                                                       */
  /* ********************************************************************* */

  template<class VECT>
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
      size_type nbd1 = mf.nb_dof_of_element(cv);
      size_type nbd2 = pfe->nb_dof();

      o << "DIM = " << int(pgt->dim()) << endl;

      if (pf1->need_G()) 
	bgeot::vectors_to_base_matrix(G, mf.linked_mesh().points_of_convex(cv));

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
	   gmm::sub_vector(U, gmm::sub_index(mf.ind_dof_of_element(cv))), pt3,cv);
	for (size_type j = 0; j < P; ++j) o << pt3[j] << " ";
	o << endl;
      }
      o << endl;
    }
    o << "END DATA ELEMENT" << endl;
    o.close();
  }

  template<class VECT>
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
	  for (size_type i = 0; i <= (size_type(1) << DIM); ++i) {
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
     interpolation of a solution on same mesh.
     - &mf_target.linked_mesh() == &mf_source.linked_mesh()
     - mf_target must be of lagrange type.
     - mf_target's qdim should be equal to mf_source qdim, or equal to 1
     - U.size() >= mf_source.get_qdim()
     - V.size() >= (mf_target.nb_dof() / mf_target.get_qdim()) * mf_source.get_qdim()
  */

  template<class VECT>
    void interpolation_solution_same_mesh(const mesh_fem &mf_source,
					  const mesh_fem &mf_target,
					  const VECT &U, VECT &V)
  {
    base_matrix G;
    size_type qdim = mf_source.get_qdim();
    base_vector coeff, val(qdim);
    if ( &(mf_source.linked_mesh()) != &(mf_target.linked_mesh()))
      DAL_THROW(failure_error, "Meshes should be the same in this function.");

    if (qdim != mf_target.get_qdim() && mf_target.get_qdim() != 1)
      DAL_THROW(failure_error, "Attempt to interpolate a field of dimension " << 
		qdim << " on a mesh_fem whose Qdim is " << 
		int(mf_target.get_qdim()));
    size_type qmult = mf_source.get_qdim()/mf_target.get_qdim();
    fem_precomp_pool fppool;
    /* we should sort convex by their fem */
    for (dal::bv_visitor cv(mf_source.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt = mf_source.linked_mesh().trans_of_convex(cv);
      pfem pf_s = mf_source.fem_of_element(cv);
      pfem pf_t = mf_target.fem_of_element(cv);
      size_type nbd_s = pf_s->nb_dof();
      size_type nbd_t = pf_t->nb_dof();
      coeff.resize(nbd_s*qdim);
      ref_mesh_dof_ind_ct::iterator itdof;
      itdof = mf_source.ind_dof_of_element(cv).begin();
      for (size_type k = 0; k < mf_source.nb_dof_of_element(cv); ++k, ++itdof)
	coeff[k] = U[*itdof];
      if (pf_s->need_G()) 
	bgeot::vectors_to_base_matrix(G, mf_source.linked_mesh().points_of_convex(cv));

      if (pf_s->target_dim() != 1 || pf_t->target_dim() != 1)
	DAL_THROW(to_be_done_error, "vector FEM interpolation still to be done ... ");
      pfem_precomp pfp = fppool(pf_s, pf_t->node_tab());
      fem_interpolation_context ctx(pgt,pfp,size_type(-1),G,cv);
      itdof = mf_target.ind_dof_of_element(cv).begin();
      for (size_type i = 0; i < nbd_t; ++i, itdof+=mf_target.get_qdim()) {
	size_type dof_t = *itdof*qmult;
	/* faux dans le cas des éléments vectoriel.                        */
	ctx.set_ii(i);
	pf_s->interpolation(ctx, coeff, val, qdim);
	for (size_type k=0; k < qdim; ++k) 
	  V[dof_t + k] = val[k];
      }
    }
  }


  /**
     interpolation of a solution on another mesh.
     - mf_target must be of lagrange type.
     - the solution should be continuous..
   */
  template<class VECT>
    void interpolation_solution(const mesh_fem &mf_source, const mesh_fem &mf_target,
				const VECT &U, VECT &V) {
    bgeot::geotrans_inv gti;
    dal::dynamic_array<base_node> ptab;
    dal::dynamic_array<size_type> itab;    
    base_matrix G;
    
    if (&mf_source.linked_mesh() == &mf_target.linked_mesh()) {
      interpolation_solution_same_mesh(mf_source, mf_target, U, V);
      return;
    }
    size_type qdim_s = mf_source.get_qdim();
    size_type qdim_t = mf_target.get_qdim();
    if (qdim_s != qdim_t && qdim_t != 1)
      DAL_THROW(failure_error, "Attempt to interpolate a field of dimension " << 
		qdim_s << " on a mesh_fem whose Qdim is " << qdim_t);
    /* initialisation of the geotrans_inv */
    {
      dal::bit_vector dof_done; dof_done.sup(0,mf_target.nb_dof());
      /* store all dof nodes into the geotrans_inv */
      for (dal::bv_visitor cv(mf_target.convex_index()); !cv.finished(); ++cv) {
	pfem pf_t = mf_target.fem_of_element(cv);
	if (pf_t->target_dim() != 1) DAL_THROW(failure_error, "still some work to do on vector FEMs!");  
	for (size_type j=0; j < pf_t->nb_dof(); ++j) {
	  size_type dof_t = mf_target.ind_dof_of_element(cv)[j*qdim_t];
	  if (!dof_done[dof_t]) {
	    // TODO: add a function in getfem_fem for inquiry about the lagrangitude of a dof ..
	    //if (dof_is_lagrange(pf_t->dof_types()[j])) {
	    gti.add_point_with_id(mf_target.point_of_dof(dof_t), dof_t);
	    dof_done.add(dof_t);
	    //} else DAL_WARNING(1, "ignoring non lagrange dof\n");
	  }
	}
      }
    }

    /* interpolation */
    dal::bit_vector ddl_touched; ddl_touched.add(0, mf_target.nb_dof());

    for (dal::bv_visitor cv(mf_source.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt=mf_source.linked_mesh().trans_of_convex(cv);
      size_type nb = gti.points_in_convex(mf_source.linked_mesh().convex(cv),
					  pgt, ptab, itab);
      if (!nb) continue;
      pfem pf_s = mf_source.fem_of_element(cv);
      if (pf_s->need_G()) 
	bgeot::vectors_to_base_matrix(G, mf_source.linked_mesh().points_of_convex(cv));

      fem_interpolation_context ctx(pgt,pf_s,base_node(),G,cv);
      base_vector coeff(mf_source.nb_dof_of_element(cv));
      gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf_source.ind_dof_of_element(cv))), coeff);
      for (size_type i = 0; i < nb; ++i) {
	size_type dof_t = itab[i];
	if (ddl_touched[dof_t]) {
	  ddl_touched.sup(dof_t, qdim_t);
	  ctx.set_xref(ptab[i]);
	  size_type pos = dof_t*(qdim_s/qdim_t);
	  typename gmm::sub_vector_type<VECT*, gmm::sub_interval>::vector_type dest = 
	    gmm::sub_vector(V,gmm::sub_interval(pos,qdim_s));
	  pf_s->interpolation(ctx, coeff, dest, qdim_s);
	  pos += qdim_s;
	}
      }
    }
    if (ddl_touched.card() != 0) {
      cerr << "WARNING : in interpolation_solution (different meshes),"
	   << ddl_touched.card() << " dof of the target mesh_fem have not been touched\n";
      for (dal::bv_visitor d(ddl_touched); !d.finished(); ++d) {
        cerr << "ddl_touched[" << d << "]=" << mf_target.point_of_dof(d) << "\n";
      }
    }
  }

  /** 
      interpolation of a solution on a set of sparse points filled in
      the provided geotrans_inv object.  The gradient is also
      interpolated of PVGRAD is non-null.
   */
  template<class VECT>
    void interpolation_solution(const mesh_fem &mf_source, bgeot::geotrans_inv &gti,
				const VECT &U, VECT &V, VECT* PVGRAD = 0) {
    size_type mdim = mf_source.linked_mesh().dim();
    size_type qdim = mf_source.get_qdim();
    dal::dynamic_array<base_node> ptab;
    dal::dynamic_array<size_type> itab;
    base_vector coeff;
    base_matrix G;

    dal::bit_vector ddl_touched; ddl_touched.add(0, gti.nb_points());

    for (dal::bv_visitor cv(mf_source.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt=mf_source.linked_mesh().trans_of_convex(cv);
      size_type nb = gti.points_in_convex(mf_source.linked_mesh().convex(cv),
					  pgt, ptab, itab);
      if (nb == 0) continue;
      pfem pf_s = mf_source.fem_of_element(cv);

      if (pf_s->need_G())
	bgeot::vectors_to_base_matrix(G, mf_source.linked_mesh().points_of_convex(cv));
      // cerr << "is_equiv:" << pf_s->is_equivalent() << ",inerp: G=" << G << ",nrow=" << G.nrows() << ", ncols=" << G.ncols() << endl;
      size_type nbd_s = pf_s->nb_dof();
      coeff.resize(nbd_s);

      fem_interpolation_context ctx(pgt,pf_s,base_node(),G,cv);
      gmm::resize(coeff, mf_source.nb_dof_of_element(cv));
      gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf_source.ind_dof_of_element(cv))), coeff);

      for (size_type i = 0; i < nb; ++i) {
	size_type dof_t = itab[i];
	if (ddl_touched[dof_t])
	{ // inverser les deux boucles pour gagner du temps ?
	  // Il faut verifier que le ddl est bien de Lagrange ...
	  ctx.set_xref(ptab[i]);
	  // cerr << "cv=" << cv << ", ptab[" << i << "]=" << ptab[i] << ", coeff=" << coeff << endl;
	  typename gmm::sub_vector_type<VECT*, gmm::sub_interval>::vector_type dest = 
	    gmm::sub_vector(V,gmm::sub_interval(dof_t*qdim,qdim));
	  pf_s->interpolation(ctx, coeff, dest, qdim);
	  if (PVGRAD) {
	    base_matrix grad(mdim, qdim);
	    pf_s->interpolation_grad(ctx, coeff, grad, qdim);
	    std::copy(grad.begin(), grad.end(), V.begin() + dof_t*qdim*mdim);
	  }
	  ddl_touched.sup(dof_t, qdim);
	}
      }
    }
    if (ddl_touched.card() != 0)
      cerr << "WARNING : in interpolation_solution (set of sparse points),"
	   << ddl_touched.card() << " points have not been touched\n";
  }



  /* ugly remplacement for "pre-qdim era" interpolation ( contrib/compare_solutions uses that.. ) */
  template<class VECT>
    void interpolation_solution(mesh_fem &mf, const mesh_fem &mf_target,
				const VECT &U, VECT &V, dim_type P) {
    if (mf.get_qdim() == 1) { mf.set_qdim(P); interpolation_solution(mf,mf_target,U,V); mf.set_qdim(1); }
    else interpolation_solution(mf,mf_target,U,V);
  }

  /* this function has nothing to do here .. */
  void classical_mesh_fem(mesh_fem& mf, short_type K);

  /** 
      export class to VTK ( http://www.kitware.com/vtk.html ) file format 
      (not the XML format, but the old format)
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
    /* choose the mesh that will be exported -- if you are juste calling 
       write_dataset(mf), you are not required to call this function */
    void set_mesh(const getfem_mesh &m, unsigned nrefine=1);
    /* if you do want to export a mesh slice instead of a whole mesh,
       use set_slice instead of set_mesh */
    void set_slice(const stored_mesh_slice& sl);
    /* the header is the second line of text in the exported file,
       you can put whatever you want */
    void set_header(const std::string& s);
    /* write a field described on mf */
    template<class VECT> void write_dataset(const getfem::mesh_fem &mf, const VECT& U0, const std::string& name);
    /* write a field already interpolated on the mesh_slice */
    template<class VECT> void write_dataset(const VECT& U, const std::string& name);

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
  void vtk_export::write_dataset(const getfem::mesh_fem &mf, const VECT& U, const std::string& name) {
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
    } else DAL_THROW(dal::dimension_error, "vtk does not accept vectors of dimension > 3");
    write_separ();
  }
}  /* end of namespace getfem.                                             */


#endif /* GETFEM_EXPORT_H__  */
