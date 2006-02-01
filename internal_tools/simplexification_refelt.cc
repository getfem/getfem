// -*- c++ -*- (enables emacs c++ mode)
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2005 Yves Renard, Julien Pommier.                    */
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
  
/**
 * Linear Elastostatic problem with a crack.
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
*/

#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_mesher.h>
#include <gmm.h>

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */


size_type simplexify(const std::vector<base_node> &pts_d,
		     const std::vector<base_node> &pts, std::ostream &f) {

  size_type n = gmm::vect_size(pts[0]);
  gmm::dense_matrix<size_type> simplexes;
  getfem::mesh m;
  
  getfem::delaunay(pts_d, simplexes);

  for (size_type i = 0; i < pts.size(); ++i) m.add_point(pts[i]);
  
  for (size_type i = 0; i < gmm::mat_ncols(simplexes); ++i)
    m.add_simplex(n, gmm::vect_begin(gmm::mat_col(simplexes, i)));
  
  scalar_type qmin = 1.0;
  for (dal::bv_visitor i(m.convex_index()); !i.finished(); ++i) {
    scalar_type q = m.convex_quality_estimate(i);
    if (m.convex_quality_estimate(i) < 1e-5) m.sup_convex(i);
    else qmin = std::min(qmin, q);
  }
  
  cout << "quality min : " << qmin << " nbconvexes : "
       << m.convex_index().card() << endl;
  
  m.optimize_structure();

  f << "[" << m.convex_index().card() * (n+1) << "] = {\n  ";
  int nb_printed = 0;
  for (dal::bv_visitor i(m.convex_index()); !i.finished(); ++i) {
    for (size_type j = 0; j <= n; ++j) {
      if (nb_printed == 18) { f << "\n  "; nb_printed = 0; }
      if (m.ind_points_of_convex(i)[j] < 10) f << " ";
      f << " " << m.ind_points_of_convex(i)[j];
      if (j != n || i != m.convex_index().card()-1) f << ",";
      nb_printed ++;
    }
  }
  f << "\n  };\n";
  return m.convex_index().card();
}




int main(int argc, char *argv[]) {

  DAL_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  //getfem::getfem_mesh_level_set_noisy();


  try {

    bgeot::pconvex_ref pref;
    size_type nb;
    std::vector<base_node> pts;

    std::ofstream f("bgeot_convex_ref_simplexified.cc");

    f <<
      "// -*- c++ -*- (enables emacs c++ mode)\n"
      "//========================================================================\n"
      "//\n"
      "// Library : Basic GEOmetric Tool  (bgeot)\n"
      "// File    : bgeot_convex_ref_simplexified.cc : simplexification of\n"
      "//           convexes of reference\n"
      "//           \n"
      "// Date    : January 21, 2006.\n"
      "// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>\n"
      "//\n"
      "//========================================================================\n"
      "//\n"
      "// Copyright (C) 2006-2006 Yves Renard\n"
      "//\n"
      "// This file is a part of GETFEM++\n"
      "//\n"
      "// This program is free software; you can redistribute it and/or modify\n"
      "// it under the terms of the GNU General Public License as published by\n"
      "// the Free Software Foundation; version 2 of the License.\n"
      "//\n"
      "// This program is distributed in the hope that it will be useful,\n"
      "// but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
      "// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
      "// GNU General Public License for more details.\n"
      "// You should have received a copy of the GNU General Public License\n"
      "// along with this program; if not, write to the Free Software Foundation,\n"
      "// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.\n"
      "//\n"
      "//========================================================================\n\n\n";
      


    f << "#include <bgeot_convex_ref.h>\n\n";
    f << "\n namespace bgeot {\n\n";


    // parallelepipeds

//     f <<
//       "  static size_type simplexified_parallelepiped_2[6] = {\n"
//       "    3,  1,  0,  2,  3,  0\n"
//       "  };\n\n"
//       "  static size_type simplexified_parallelepiped_2_nb = 2;\n\n"
//       "  static size_type simplexified_parallelepiped_3[24] = {\n"
//       "    0,  4,  5,  7,  0,  3,  1,  7,  0,  5,  7,  1,  0,  4,  6,  7,  0,  2,\n"
//       "    3,  7,  0,  2,  6,  7\n"
//       "  };\n\n"
//       "  static size_type simplexified_parallelepiped_3_nb = 6;\n";

   for (size_type n = 2; n < 7; ++n) {
      
      pref = bgeot::parallelepiped_of_reference(n);
      cout << "simplexification of parallelepiped of dimension " << n << endl;
      pts = pref->points();
      
      // small decay in order to have matching meshes
      
      base_small_vector v(n); v.fill(0.5);
      if (n < 4) {
	for (size_type ip = 0; ip < pts.size(); ip += pts.size()-1)
	  { pts[ip] -= v; pts[ip] *= 0.9; pts[ip] += v; }
      }
      else {
	for (size_type ip = 0; ip < pts.size(); ++ip) {
	  size_type nb1 = 0;
	  for (size_type id = 0; id < n; ++id)
	    if (gmm::abs(pts[ip][id] - 1.0) < 1E-8) ++nb1;
	  if (nb1 == 2)
	    { pts[ip] -= v; pts[ip] *= 0.9; pts[ip] += v; }
	}
      }

      f << "\n  static size_type simplexified_parallelepiped_" << n;
      nb = simplexify(pts, pref->points(), f);
      f << "\n  static size_type simplexified_parallelepiped_" << n << "_nb = "
	<< nb << ";\n";
    }

    // prisms
    for (size_type n = 3; n < 7; ++n) {
      
      pref = bgeot::prism_of_reference(n);
      cout << "simplexification of prism of dimension " << n << endl;

      f << "\n  static size_type simplexified_prism_" << n;
      nb = simplexify(pref->points(), pref->points(), f);
      f << "\n  static size_type simplexified_prism_" << n << "_nb = "
	<< nb << ";\n";
    }

    f << "\n\n\n";
    f << "  size_type simplexified_tab(pconvex_structure cvs,\n"
      << "                             size_type **tab) {\n";
    for (size_type n = 2; n < 7; ++n) {
      f << "    if (cvs == parallelepiped_structure(" << n << ")) {\n";
      f << "      *tab = simplexified_parallelepiped_" << n << ";\n";
      f << "      return simplexified_parallelepiped_" << n << "_nb;\n";
      f << "    }\n\n";
    }
    for (size_type n = 3; n < 7; ++n) {
      f << "    if (cvs == prism_structure(" << n << ")) {\n";
      f << "      *tab = simplexified_prism_" << n << ";\n";
      f << "      return simplexified_prism_" << n << "_nb;\n";
      f << "    }\n\n";
    }
    f << "    DAL_THROW(failure_error, \"No simplexification "
      << " for this element\");\n";
    
    f << "  }\n\n";

    // refinement of simplexes
    
    for (size_type n = 1; n < 7; ++n) {
      
      cout << "refinement of simplex of dimension " << n << endl;

      pref = bgeot::equilateral_simplex_of_reference(n);
      bgeot::pconvex_ref pref2 = bgeot::simplex_of_reference(n, 2);
      pts = pref2->points();
      base_node barycentre = dal::mean_value(pref->points());

      bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(n, 1);
      for (size_type i = 0; i < pts.size(); ++i)
	pts[i] = pgt->transform(pts[i], pref->points());

      std::vector<base_node> pts2 = pts;

      for (size_type ip = 0; ip < pts.size(); ++ip) {
	size_type nb1 = 0;
	for (size_type id = 0; id < n; ++id)
	  if (gmm::abs(pref2->points()[ip][id] - 0.5) < 1E-8) ++nb1;
	if (nb1 >= 1) {
	  pts[ip] -= barycentre; pts[ip] *= 0.7; pts[ip] += barycentre;
	}
      }

      f << "\n  static size_type refinement_simplex_" << n;
      // nb = simplexify(pts, pref->points(), f);
      nb = simplexify(pts, pts2, f);
      f << "\n  static size_type refinement_simplex_" << n << "_nb = "
	<< nb << ";\n";
    }


    f << "\n\n\n";
    f << "  size_type refinement_simplexe_tab(size_type n,\n"
      << "                                    size_type **tab) {\n"
      << "    switch(n) {\n";
    for (size_type d = 1; d < 7; ++d)
      f  << "    case " << d << " : *tab = refinement_simplex_" << d << ";\n"
	 << "             return refinement_simplex_" << d << "_nb;\n";
    f << "    default : DAL_THROW(failure_error, \"No refinement for "
      << " this element\");\n    }\n";
    f << "  }\n\n";


    f << "}\n";




    f.close();

  }
  DAL_STANDARD_CATCH_ERROR;

  return 0; 
}
