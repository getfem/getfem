/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Elements Methods (getfem)            */
/* File    :  getfem_regular_meshes.C : basic refinement of meshes.        */
/*     									   */
/*                                                                         */
/* Date : December 20, 1999.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 1999-2002  Yves Renard.                                   */
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


#include <getfem_regular_meshes.h>

namespace getfem
{
  void parallelepiped_regular_simplex_mesh_(getfem_mesh &me, dim_type N,
    const base_node &org, const base_small_vector *ivect, const size_type *iref)
  {
    bgeot::mesh_structure cvt, sl;
    bgeot::convex<base_node>
      pararef = *(bgeot::parallelepiped_of_reference(N));
    base_node a = org;
    size_type i, nbpt = pararef.nb_points();

    // cout << pararef << endl;

    for (i = 0; i < nbpt; ++i)
    {
      a.fill(0.0);
      for (dim_type n = 0; n < N; ++n)
	gmm::add(gmm::scaled(ivect[n],pararef.points()[i][n]),a);
	//a.addmul(pararef.points()[i][n], ivect[n]);
      pararef.points()[i] = a;
    }

    // cout << pararef << endl;
    
    cvt.add_convex(pararef.structure(), dal::sequence_iterator<size_type>(0));
    cvt.to_edges();
    // cout << "cvt " << endl;
    // cvt.write_to_file(cout); getchar();
    bgeot::simplexify(cvt, sl, pararef.points(), N, me.eps());
    // cout << "sl " << endl;
    // sl.write_to_file(cout); getchar();

    size_type nbs = sl.nb_convex();
    std::vector<size_type> tab1(N+1), tab(N), tab3(nbpt);
    size_type total = 0;
    std::fill(tab.begin(), tab.end(), 0);
    while (tab[N-1] != iref[N-1])
    {
      for (a = org, i = 0; i < N; i++) 
	gmm::add(gmm::scaled(ivect[i],scalar_type(tab[i])),a);
        //a.addmul(scalar_type(tab[i]), ivect[i]);

      for (i = 0; i < nbpt; i++)
	tab3[i] = me.add_point(a + pararef.points()[i]);
	
      for (i = 0; i < nbs; i++)
      {
	bgeot::ref_mesh_point_ind_ct tab2 = sl.ind_points_of_convex(i);
	for (dim_type l = 0; l < N+1; l++)
	  tab1[l] = tab3[(tab2[l]
			  + (((total & 1) && N > 1) ? (nbpt/2) : 0)) % nbpt];
	me.add_simplex(N, tab1.begin());
      }

      for (dim_type l = 0; l < N; l++)
      {
	tab[l]++; total++;
	if (l < N-1 && tab[l] >= iref[l]) { total -= tab[l]; tab[l] = 0; }
	else break;
      }
    }
  }


  void parallelepiped_regular_prism_mesh_(getfem_mesh &me, dim_type N,
    const base_node &org, const base_small_vector *ivect, const size_type *iref)
  {
    getfem_mesh aux;
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



  void parallelepiped_regular_mesh_(getfem_mesh &me, dim_type N,
                                    const base_node &org, const base_small_vector *ivect, const size_type *iref, bool linear_gt)
  {
    bgeot::convex<base_node>
      pararef = *(bgeot::parallelepiped_of_reference(N));
    base_node a = org;
    size_type i, nbpt = pararef.nb_points();

    for (i = 0; i < nbpt; ++i)
    {
      a.fill(0.0);
      for (dim_type n = 0; n < N; ++n)
	gmm::add(gmm::scaled(ivect[n],pararef.points()[i][n]),a);
	//a.addmul(pararef.points()[i][n], ivect[n]);
      pararef.points()[i] = a;
    }

    std::vector<size_type> tab1(N+1), tab(N), tab3(nbpt);
    size_type total = 0;
    std::fill(tab.begin(), tab.end(), 0);
    while (tab[N-1] != iref[N-1])
    {
      for (a = org, i = 0; i < N; i++) 
	gmm::add(gmm::scaled(ivect[i], scalar_type(tab[i])),a);
	//a.addmul(scalar_type(tab[i]), ivect[i]);

      for (i = 0; i < nbpt; i++)
	tab3[i] = me.add_point(a + pararef.points()[i]);
      me.add_convex(linear_gt ? 
                    bgeot::parallelepiped_linear_geotrans(N) : 
                    bgeot::parallelepiped_geotrans(N, 1), tab3.begin());
     
      for (dim_type l = 0; l < N; l++)
      {
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
      c2*=(.5 - dal::abs(x[i]-.5));
    }
    z[0] = x[0] + c1;
    for (size_type i=1; i < x.size(); ++i) {
      z[i] = x[i] + c2/10.;
    }
    return z;
  }

  void regular_unit_mesh(getfem_mesh& m, std::vector<size_type> nsubdiv, 
			 bgeot::pgeometric_trans pgt, bool noised) {
    getfem_mesh mesh;
    size_type N = nsubdiv.size();
    base_node org(N); org.fill(0.0);
    std::vector<base_small_vector> vtab(N);
    for (dim_type i = 0; i < N; i++) { 
      vtab[i] = base_small_vector(N); vtab[i].fill(0.0);
      (vtab[i])[i] = 1. / scalar_type(nsubdiv[i]) * 1.;
    }
    if (pgt->basic_structure() == bgeot::simplex_structure(N)) {
      getfem::parallelepiped_regular_simplex_mesh
	(mesh, N, org, vtab.begin(), nsubdiv.begin());
    } else if (pgt->basic_structure() == bgeot::parallelepiped_structure(N)) {
      getfem::parallelepiped_regular_mesh
	(mesh, N, org, vtab.begin(), nsubdiv.begin());
    } else if (pgt->basic_structure() == bgeot::prism_structure(N)) {
      getfem::parallelepiped_regular_prism_mesh
	(mesh, N, org, vtab.begin(), nsubdiv.begin());
    } else {
      DAL_THROW(dal::failure_error, "cannot build a regular mesh for " << bgeot::name_of_geometric_trans(pgt));
    }
    m.clear();
    /* build a mesh with a geotrans of degree K */
    for (dal::bv_visitor cv(mesh.convex_index()); !cv.finished(); ++cv) {
      if (pgt == mesh.trans_of_convex(cv)) {
	m.add_convex_by_points(mesh.trans_of_convex(cv), mesh.points_of_convex(cv).begin()); 
      } else {
	std::vector<base_node> pts(pgt->nb_points());
	for (size_type i=0; i < pgt->nb_points(); ++i) {
	  pts[i] = mesh.trans_of_convex(cv)->transform(pgt->convex_ref()->points()[i], 
						       mesh.points_of_convex(cv));
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
	for (size_type i=0; i < N; ++i) { if (dal::abs(P[i]) < 1e-10 || dal::abs(P[i]-1.) < 1e-10) is_border = true; }
	if (!is_border) { 
	  P = shake_func(P); 
	  for (size_type i=0; i < N; ++i) P[i] += 0.40*(1./nsubdiv[i])*dal::random(double());
	}
      }
    }
  }
}  /* end of namespace getfem.                                             */
