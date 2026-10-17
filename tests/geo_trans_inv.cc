/*===========================================================================

 Copyright (C) 2002-2026 Yves Renard, Julien Pommier.

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

===========================================================================*/

/*
  various checks for the geometric transformations inversion
  (used by many other modules, such as interpolation from one mesh onto another, etc)
 */

#include "getfem/bgeot_geotrans_inv.h"
#include "getfem/getfem_regular_meshes.h"

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

using bgeot::size_type;
using bgeot::short_type;
using bgeot::dim_type;
using bgeot::scalar_type;
using bgeot::base_node;
using bgeot::base_vector;
using bgeot::base_small_vector;
using bgeot::base_matrix;

base_matrix random_base(size_type N) {
  base_matrix M(N,N);
  do {
    for (size_type i=0; i < N; ++i) {
      base_node n(N); gmm::fill_random(n); n /= gmm::vect_norm2(n);
      for (size_type j=0; j < N; ++j) M(i,j) = n[j];
    }
  } while (gmm::abs(gmm::lu_det(M)) < 0.7);
  return M;
}

void check_inversion(bgeot::pgeometric_trans pgt, const std::vector<base_node>& cvpts,
		     bgeot::geotrans_inv_convex& gic, 
		     const base_node& P, const base_node& expected_Pref, bool expected_in, bool verbose=true) {
  base_node Pref;
  bool is_in = gic.invert(P, Pref);
  base_node P2 = pgt->transform(Pref, cvpts.begin());
  scalar_type err = gmm::vect_dist2(P,P2);
  bool multiple_solutions = false;
  if (err < 1e-10 && gmm::vect_dist2(Pref,expected_Pref) > 1e-5) { multiple_solutions = true; }
  if ((is_in != expected_in && !multiple_solutions) || err > 1e-10) {
    cerr << "Error with inversion of " << bgeot::name_of_geometric_trans(pgt) << " on convex: " << gmm::vref(cvpts) << "\n";
    cerr << "  inversion of point " << P << " gave point " << Pref;
    if (expected_in) cerr << "  the point is known to be IN the convex.";
    else cerr << "  the point is known to be OUTSIDE the convex.";
    if (is_in) cerr << " but was reported as a point IN.\n";
    else cerr << " but was reported as a point OUTSIDE.\n";
    cerr << "  distance(transform(Pref)-P) = " << err;
    if (err >= 1e-10) cerr << ": TOO LARGE (>1e-10). Inversion failed miserabily.";
    cerr << "\n\n";
    GMM_ASSERT1(false, "geotrans_inv_convex failed\n");
  }
  if (verbose) { 
    cout << "SUCCESS: ";
    if (!multiple_solutions) {
      if (is_in) cout << "POINT IN "; else cout << "POINT OUT";
    } else cout << "PGT NON INVERTIBLE AT THIS POINT";
    cout << " inv(" << P << ")=" << Pref << ", (err=" << err << ")\n"; 
  }
}

void test_inversion(bgeot::pgeometric_trans pgt, bool verbose) {
  size_type N=pgt->dim();
  std::vector<base_node> cvpts(pgt->nb_points());
  cout << "Testing geotrans_inv with " << bgeot::name_of_geometric_trans(pgt) << "\n";
  base_matrix M = random_base(N);
  base_node translat(N); gmm::fill_random(translat);
  scalar_type scale = gmm::random()*10.;
  for (size_type i=0; i < pgt->nb_points(); ++i) {
    cvpts[i] = base_node(N);
    gmm::mult(M,pgt->convex_ref()->points()[i],cvpts[i]);
    cvpts[i] += translat;
    for (size_type j=0; j < cvpts[i].size(); ++j) { 
      cvpts[i][j] *= scale;
      cvpts[i][j] += gmm::random(double())*0.05*scale;
    }
  }
  bgeot::geotrans_inv_convex gic;
  gic.init(cvpts,pgt);
  for (size_type i=0; i < cvpts.size(); ++i) {
    check_inversion(pgt,cvpts,gic,cvpts[i],pgt->convex_ref()->points()[i],true,verbose);
  }
  for (size_type i=0; i < 100; ++i) {
    base_node Pref(pgt->dim());
    for (size_type j=0; j < Pref.size(); ++j) Pref[j] = (gmm::random() * 1.5 - 0.25);
    base_node P = pgt->transform(Pref, cvpts.begin());
    check_inversion(pgt,cvpts,gic,P,Pref,pgt->convex_ref()->is_in(Pref)<1e-10,verbose);
  }
}

/* problematic test-cases .. */
void test0() {
  bgeot::geotrans_inv_convex gic;
  std::vector<base_node> cvpts(4);
  cvpts[0] = base_node(2.5, 0.6);
  cvpts[1] = base_node(5., 0.);
  cvpts[2] = base_node(2.5, 1.8);
  cvpts[3] = base_node(3.2, 1.5);
  bgeot::pgeometric_trans pgt = bgeot::geometric_trans_descriptor("GT_QK(2,1)");
  gic.init(cvpts,pgt);
  base_node X(3.132, 1.617), Xref;
  gic.invert(X, Xref);
  cout << "Xref=" << Xref << "\n";
}

void test_inversion(bool verbose) {
  for (short_type N=1; N <= 4; ++N) {
    for (short_type K=1; K < ((N<4) ? 3 : 2); ++K) {
      test_inversion(bgeot::simplex_geotrans(N,K),verbose);
    }
  }
  for (short_type N=1; N <= 3; ++N) {
    for (short_type K=1; K < 3; ++K) {
      test_inversion(bgeot::parallelepiped_geotrans(N,K),verbose);
    }
  }
  for (short_type N=2; N <= 3; ++N) {
    for (short_type K=1; K < 3; ++K) {
      test_inversion(bgeot::prism_geotrans(N,K),verbose);
    }
  }
  test_inversion(bgeot::parallelepiped_linear_geotrans(5),verbose);
  test_inversion(bgeot::prism_linear_geotrans(3),verbose);
}

/* Inversion with projection into the element, for non-affine elements:
   interior points must be recovered, and exterior points must be
   projected into the reference element. */
void test_projection_into_element() {
  const char *names[5] = { "GT_QK(3,1)", "GT_QK(3,2)", "GT_PYRAMID(1)",
                           "GT_PYRAMID(2)", "GT_PRISM(3,1)" };
  const scalar_type interior[5][3] = { { 0.8, 0.8, 0.5}, { 0.9, 0.3, 0.7},
                                       {-0.5, 0.2, 0.3}, {-0.4,-0.3, 0.2},
                                       { 0.3, 0.4, 0.9} };
  for (size_type c = 0; c < 5; ++c) {
    bgeot::pgeometric_trans pgt = bgeot::geometric_trans_descriptor(names[c]);
    cout << "Testing projection into the element with " << names[c] << "\n";
    std::vector<base_node> cvpts; // distorted reference element
    for (size_type i = 0; i < pgt->nb_points(); ++i) {
      base_node p = pgt->convex_ref()->points()[i];
      p[0] += 0.05 * p[1] * p[2];
      cvpts.push_back(p);
    }
    bgeot::geotrans_inv_convex gic;
    gic.init(cvpts, pgt);
    bool converged;
    base_node xref(3), xr(3);
    for (size_type k = 0; k < 3; ++k) xref[k] = interior[c][k];
    gic.invert(pgt->transform(xref, cvpts), xr, converged, 1e-12, true);
    GMM_ASSERT1(gmm::vect_dist2(xr, xref) < 1e-9, "wrong inversion of the "
                "interior point " << xref << " for " << names[c] << ": " << xr);
    for (size_type i = 0; i < 6; ++i)
      for (size_type j = 0; j < 6; ++j)
        for (size_type k = 0; k < 6; ++k) {
          base_node xo(-1. + 0.6*double(i), -1. + 0.6*double(j),
                       -1. + 0.6*double(k));
          gic.invert(pgt->transform(xo, cvpts), xr, converged, 1e-12, true);
          GMM_ASSERT1(pgt->convex_ref()->is_in(xr) < 1e-12, "point " << xo
                      << " projected outside of " << names[c] << ": " << xr);
        }
  }
}

/* Inversion near the apex of pyramids, where Newton stalls at the level
   of rounding errors above the residual tolerance. */
void test_inversion_near_pyramid_apex() {
  for (const char *name : { "GT_PYRAMID(1)", "GT_PYRAMID(2)" }) {
    bgeot::pgeometric_trans pgt = bgeot::geometric_trans_descriptor(name);
    cout << "Testing inversion near the apex of " << name << "\n";
    std::vector<base_node> cvpts(pgt->convex_ref()->points().begin(),
                                 pgt->convex_ref()->points().end());
    bgeot::geotrans_inv_convex gic;
    gic.init(cvpts, pgt);
    for (scalar_type d : { 0.1, 0.03, 0.01 }) {
      base_node xref(0.2*d, -0.1*d, 1.-d), xr(3);
      bool converged;
      bool in = gic.invert(pgt->transform(xref, cvpts), xr, converged, 1e-12);
      GMM_ASSERT1(in && converged && gmm::vect_dist2(xr, xref) < 1e-10,
                  "inversion failed near the apex of " << name << " for "
                  << xref << ": " << xr);
    }
  }
}

/* Points off a curved element of lower dimension than the space must
   not be accepted, although the Newton correction vanishes for them. */
void test_inversion_off_lower_dimensional_element() {
  for (const char *name : { "GT_PK(1,2)", "GT_PK(2,2)", "GT_QK(2,2)" }) {
    bgeot::pgeometric_trans pgt = bgeot::geometric_trans_descriptor(name);
    size_type P = pgt->dim(), N = P + 1;
    cout << "Testing inversion off a curved " << name << " in dimension "
         << N << "\n";
    std::vector<base_node> cvpts;
    for (const base_node &pt : pgt->convex_ref()->points()) {
      base_node x(N);
      for (size_type k = 0; k < P; ++k) x[k] = pt[k];
      x[N-1] = 0.3 * pt[0] * pt[0];
      cvpts.push_back(x);
    }
    bgeot::geotrans_inv_convex gic;
    gic.init(cvpts, pgt);
    base_node xref(P), xr(P);
    for (size_type k = 0; k < P; ++k) xref[k] = 0.2 + 0.1 * double(k);
    base_node x = pgt->transform(xref, cvpts);
    bool converged;
    GMM_ASSERT1(gic.invert(x, xr, converged, 1e-12) && converged,
                "point on " << name << " not found");
    x[N-1] += 1e-6;
    GMM_ASSERT1(!gic.invert(x, xr, converged, 1e-12) && !converged,
                "point off " << name << " accepted");
  }
}

int main(int argc, char *argv[]) {
  dim_type N, MESH_TYPE;
  scalar_type LX, LY, LZ;
  size_type NX, NB_POINTS;
  bgeot::md_param PARAM;
  getfem::mesh mesh;

  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {
    test0();
    test_inversion(true);
    test_projection_into_element();
    test_inversion_near_pyramid_apex();
    test_inversion_off_lower_dimensional_element();
    PARAM.read_command_line(argc, argv);
    N = bgeot::dim_type(PARAM.int_value("N", "Domaine dimension"));
    NB_POINTS = PARAM.int_value("NB_POINTS", "Nb points");
    LX = PARAM.real_value("LX", "Size in X");
    LY = PARAM.real_value("LY", "Size in Y");
    LZ = PARAM.real_value("LZ", "Size in Y");
    NX = PARAM.int_value("NX", "Nomber of sace steps ");
    MESH_TYPE = bgeot::dim_type(PARAM.int_value("MESH_TYPE", "Mesh type "));
 
    cout << "Mesh generation\n";

    base_node org(N); gmm::fill(org,1);
    std::vector<base_small_vector> vtab(N);
    std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NX);
    for (dim_type i = 0; i < N; i++) { 
      vtab[i] = base_small_vector(N); gmm::clear(vtab[i]);
      (vtab[i])[i] = ((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NX) * 1.;
    }
    // if (N > 1) vtab[N-1][0] = incline * LX / scalar_type(NX);
    
    switch (MESH_TYPE) {
    case 0 : getfem::parallelepiped_regular_simplex_mesh
		     (mesh, N, org, vtab.begin(), ref.begin()); break;
    case 1 : getfem::parallelepiped_regular_mesh
		     (mesh, N, org, vtab.begin(), ref.begin()); break;
    case 2 : getfem::parallelepiped_regular_prism_mesh
		     (mesh, N, org, vtab.begin(), ref.begin()); break;
    default : GMM_ASSERT1(false, "Unknown type of mesh");
    }
    
    mesh.optimize_structure();


    scalar_type exectime = gmm::uclock_sec(), total_time = 0.0;

    bgeot::geotrans_inv gti;
    bgeot::base_node pt(N);
    dal::dynamic_array<base_node> ptab;
    dal::dynamic_array<size_type> itab;

    for (size_type i = 0; i < NB_POINTS; ++i) {
      for (dim_type k = 0; k < N; ++k) 
	pt[k] = gmm::random() + 1; //double());
      //cout << "point " << i << " : " << pt << "\n";
      gti.add_point(pt);
    }

    for (size_type i=0; i < 2; ++i) {
      total_time = 0;
      exectime = gmm::uclock_sec();
      if (i==0) { cout << " using points_in_box..\n"; }
      else { cout << " using brute force..\n"; }
      cout << "Time to sort points : " << gmm::uclock_sec() - exectime << endl;
      total_time += gmm::uclock_sec() - exectime;
      
      dal::bit_vector nn = mesh.convex_index();
      size_type nbtot = 0;
      for (size_type cv = nn.take_first(); cv != size_type(-1); cv << nn) {
	
	size_type nb = gti.points_in_convex(mesh.convex(cv),
					    mesh.trans_of_convex(cv),
					    ptab, itab, i==1);
	//if ((cv % 100) == 0) cout << "cv : " << cv << endl;
	nbtot += nb;
      }
      
      cout << "Time to invert geo trans : " << gmm::uclock_sec() - exectime
	   << endl;
      cout << "Total number : " << nbtot << endl;
      assert(nbtot == NB_POINTS);
    }
    total_time += gmm::uclock_sec() - exectime;
    
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0;
}
