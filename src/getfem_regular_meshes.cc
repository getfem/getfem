// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Elements Methods (getfem)
// File    : getfem_regular_meshes.cc : basic refinement of meshes.
//           
// Date    : December 20, 1999.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1999-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================


#include <getfem_regular_meshes.h>

namespace getfem
{
  void parallelepiped_regular_simplex_mesh_
  (mesh &me, dim_type N, const base_node &org,
   const base_small_vector *ivect, const size_type *iref) {
    bgeot::convex<base_node>
      pararef = *(bgeot::parallelepiped_of_reference(N));

    if (N >= 5) DAL_WARNING1("CAUTION : Simplexification in dimension >= 5 "
			    "has not been tested and the resulting mesh "
			    "should be not conformal");
    
    const bgeot::mesh_structure &sl
      = *(bgeot::parallelepiped_of_reference(N)->simplexified_convex());
    
    base_node a = org;
    size_type i, nbpt = pararef.nb_points();
    
    for (i = 0; i < nbpt; ++i) {
      gmm::clear(a);
      for (dim_type n = 0; n < N; ++n)
	gmm::add(gmm::scaled(ivect[n],pararef.points()[i][n]),a);
      pararef.points()[i] = a;
    }

    // bgeot::simplexify(cvt, sl, pararef.points(), N, me.eps());

    size_type nbs = sl.nb_convex();
    std::vector<size_type> tab1(N+1), tab(N), tab3(nbpt);
    size_type total = 0;
    std::fill(tab.begin(), tab.end(), 0);
    while (tab[N-1] != iref[N-1]) {
      for (a = org, i = 0; i < N; i++) 
	gmm::add(gmm::scaled(ivect[i],scalar_type(tab[i])),a);
        //a.addmul(scalar_type(tab[i]), ivect[i]);

      for (i = 0; i < nbpt; i++)
	tab3[i] = me.add_point(a + pararef.points()[i]);
	
      for (i = 0; i < nbs; i++) {
	const mesh::ind_cv_ct &tab2 = sl.ind_points_of_convex(i);
	for (dim_type l = 0; l <= N; l++)
	  // tab1[l] = tab3[tab2[l]];
	  tab1[l] = tab3[(tab2[l]
			  + (((total & 1) && N != 3) ? (nbpt/2) : 0)) % nbpt];
	me.add_simplex(N, tab1.begin());
      }

      for (dim_type l = 0; l < N; l++) {
	tab[l]++; total++;
	if (l < N-1 && tab[l] >= iref[l]) { total -= tab[l]; tab[l] = 0; }
	else break;
      }
    }
  }


  void parallelepiped_regular_prism_mesh_
  (mesh &me, dim_type N, const base_node &org,
   const base_small_vector *ivect, const size_type *iref) {
    mesh aux;
    parallelepiped_regular_simplex_mesh_(aux, N-1, org, ivect, iref);
    std::vector<base_node> ptab(2 * N);
    
    for (dal::bv_visitor cv(aux.convex_index()); !cv.finished(); ++cv) {
      std::copy(aux.points_of_convex(cv).begin(),
		aux.points_of_convex(cv).end(), ptab.begin());

      for (size_type k = 0; k < iref[N-1]; ++k) {
	
	for (dim_type j = 0; j < N; ++j) ptab[j+N] = ptab[j] + ivect[N-1];
	me.add_prism_by_points(N, ptab.begin());
	
	std::copy(ptab.begin()+N, ptab.end(), ptab.begin());
      }
    }
  }



  void parallelepiped_regular_mesh_
  (mesh &me, dim_type N, const base_node &org,
   const base_small_vector *ivect, const size_type *iref, bool linear_gt) {
    bgeot::convex<base_node>
      pararef = *(bgeot::parallelepiped_of_reference(N));
    base_node a = org;
    size_type i, nbpt = pararef.nb_points();
    
    for (i = 0; i < nbpt; ++i) {
      gmm::clear(a);
      for (dim_type n = 0; n < N; ++n)
	gmm::add(gmm::scaled(ivect[n],pararef.points()[i][n]),a);
      //a.addmul(pararef.points()[i][n], ivect[n]);
      pararef.points()[i] = a;
    }

    std::vector<size_type> tab1(N+1), tab(N), tab3(nbpt);
    size_type total = 0;
    std::fill(tab.begin(), tab.end(), 0);
    while (tab[N-1] != iref[N-1]) {
      for (a = org, i = 0; i < N; i++) 
	gmm::add(gmm::scaled(ivect[i], scalar_type(tab[i])),a);
      //a.addmul(scalar_type(tab[i]), ivect[i]);
      
      for (i = 0; i < nbpt; i++)
	tab3[i] = me.add_point(a + pararef.points()[i]);
      me.add_convex(linear_gt ? 
                    bgeot::parallelepiped_linear_geotrans(N) : 
                    bgeot::parallelepiped_geotrans(N, 1), tab3.begin());
      
      for (dim_type l = 0; l < N; l++) {
	tab[l]++; total++;
	if (l < N-1 && tab[l] >= iref[l]) { total -= tab[l]; tab[l] = 0; }
	else break;
      }
    }
  }

  /* deformation inside a unit square  -- ugly */
  static base_node shake_func(const base_node& x) {
    base_node z(x.size());
    scalar_type c1 = 1., c2 = 1.;
    for (size_type i=0; i < x.size(); ++i) {
      c1*=(x[i]*(1.-x[i]));
      c2*=(.5 - gmm::abs(x[i]-.5));
    }
    z[0] = x[0] + c1;
    for (size_type i=1; i < x.size(); ++i) {
      z[i] = x[i] + c2/10.;
    }
    return z;
  }

  void regular_unit_mesh(mesh& m, std::vector<size_type> nsubdiv, 
			 bgeot::pgeometric_trans pgt, bool noised) {
    mesh msh;
    size_type N = nsubdiv.size();
    base_node org(N); gmm::clear(org);
    std::vector<base_small_vector> vtab(N);
    for (dim_type i = 0; i < N; i++) { 
      vtab[i] = base_small_vector(N); gmm::clear(vtab[i]);
      (vtab[i])[i] = 1. / scalar_type(nsubdiv[i]) * 1.;
    }
    if (pgt->basic_structure() == bgeot::simplex_structure(N)) {
      getfem::parallelepiped_regular_simplex_mesh
	(msh, N, org, vtab.begin(), nsubdiv.begin());
    } else if (pgt->basic_structure() == bgeot::parallelepiped_structure(N)) {
      getfem::parallelepiped_regular_mesh
	(msh, N, org, vtab.begin(), nsubdiv.begin());
    } else if (pgt->basic_structure() == bgeot::prism_structure(N)) {
      getfem::parallelepiped_regular_prism_mesh
	(msh, N, org, vtab.begin(), nsubdiv.begin());
    } else {
      DAL_THROW(dal::failure_error, "cannot build a regular mesh for "
		<< bgeot::name_of_geometric_trans(pgt));
    }
    m.clear();
    /* build a mesh with a geotrans of degree K */
    for (dal::bv_visitor cv(msh.convex_index()); !cv.finished(); ++cv) {
      if (pgt == msh.trans_of_convex(cv)) {
	m.add_convex_by_points(msh.trans_of_convex(cv),
			       msh.points_of_convex(cv).begin()); 
      } else {
	std::vector<base_node> pts(pgt->nb_points());
	for (size_type i=0; i < pgt->nb_points(); ++i) {
	  pts[i] = msh.trans_of_convex(cv)->transform
	    (pgt->convex_ref()->points()[i], msh.points_of_convex(cv));
	}
	m.add_convex_by_points(pgt, pts.begin());
      }
    }
    m.optimize_structure();
    /* apply a continuous deformation + some noise */
    if (noised) {
      for (dal::bv_visitor ip(m.points().index()); !ip.finished(); ++ip) {
	bool is_border = false;
	base_node& P = m.points()[ip];
	for (size_type i=0; i < N; ++i) {
	  if (gmm::abs(P[i]) < 1e-10 || gmm::abs(P[i]-1.) < 1e-10)
	    is_border = true;
	}
	if (!is_border) { 
	  P = shake_func(P); 
	  for (size_type i=0; i < N; ++i)
	    P[i] += 0.20*(1./(nsubdiv[i]* pgt->complexity()))
			      * gmm::random(double());
	}
      }
    }
  }
  
    

  void regular_mesh(mesh& m, const std::string &st) {
    std::stringstream s(st);
    ftool::md_param PARAM;
    PARAM.read_param_file(s);
    
    std::string GT = PARAM.string_value("GT");
    if (GT.empty())
      DAL_THROW(failure_error, "regular mesh : you have at least to "
		"specify the geometric transformation");
    bgeot::pgeometric_trans pgt = 
      bgeot::geometric_trans_descriptor(GT);
    
    size_type N = pgt->dim();
    base_small_vector org(N);
    
    const std::vector<ftool::md_param::param_value> &o
      = PARAM.array_value("ORG");
    if (o.size() > 0) {
      if (o.size() != N)
	DAL_THROW(failure_error,
		  "ORG parameter should be an array of size " << N);
      for (size_type i = 0; i < N; ++i) {
	if (o[i].type_of_param() != ftool::md_param::REAL_VALUE)
	  DAL_THROW(failure_error, "ORG should be a real array.");
	org[i] = o[i].real();
      }
    }
    
    bool noised = (PARAM.int_value("NOISED") != 0);
    
    std::vector<size_type> nsubdiv(N);
    gmm::fill(nsubdiv, 2);
    const std::vector<ftool::md_param::param_value> &ns
      = PARAM.array_value("NSUBDIV");
    if (ns.size() > 0) {
      if (ns.size() != N)
	DAL_THROW(failure_error,
		  "NSUBDIV parameter should be an array of size " << N);
      for (size_type i = 0; i < N; ++i) {
	if (ns[i].type_of_param() != ftool::md_param::REAL_VALUE)
	  DAL_THROW(failure_error, "NSUBDIV should be an integer array.");
	nsubdiv[i] = size_type(ns[i].real()+0.5);
      }
    }
    
    base_small_vector sizes(N);
    gmm::fill(sizes, 1.0);
    
    const std::vector<ftool::md_param::param_value> &si
      = PARAM.array_value("SIZES");
    if (si.size() > 0) {
      if (si.size() != N)
	DAL_THROW(failure_error,
		  "SIZES parameter should be an array of size " << N);
      for (size_type i = 0; i < N; ++i) {
	if (si[i].type_of_param() != ftool::md_param::REAL_VALUE)
	  DAL_THROW(failure_error, "SIZES should be a real array.");
	sizes[i] = si[i].real();
      }
    }
    
    regular_unit_mesh(m, nsubdiv, pgt, noised);
        
    base_matrix M(N,N);
    for (size_type i=0; i < N; ++i) M(i,i) = sizes[i];
    m.transformation(M);
    m.translation(org);
    
  }




}  /* end of namespace getfem.                                             */
