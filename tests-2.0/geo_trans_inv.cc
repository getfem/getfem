/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2002-2012 Yves Renard, Julien Pommier.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/

/*
  various checks for the geometric transformations inversion
  (used by many other modules, such as interpolation from one mesh onto another, etc)
 */

#include <getfem/bgeot_geotrans_inv.h>
#include <getfem/getfem_regular_meshes.h>

using bgeot::size_type;
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
    cerr << "Error with inversion of " << bgeot::name_of_geometric_trans(pgt) << " on convex: " << cvpts << "\n";
    cerr << "  inversion of point " << P << " gave point " << Pref;
    if (expected_in) cerr << "  the point is known to be IN the convex.";
    else cerr << "  the point is known to be OUTSIDE the convex.";
    if (is_in) cerr << " but was reported as a point IN.\n";
    else cerr << " but was reported as a point OUTSIDE.\n";
    cerr << "  distance(transform(Pref)-P) = " << err;
    if (err >= 1e-10) cerr << ": TOO LARGE (>1e-10). Inversion failed miserabily.";
    cerr << "\n\n";
    GMM_THROW(dal::failure_error, "geotrans_inv_convex failed\n");
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
  for (dim_type N=1; N <= 4; ++N) {
    for (dim_type K=1; K < ((N<4) ? 3 : 2); ++K) {
      test_inversion(bgeot::simplex_geotrans(N,K),verbose);
    }
  }
  for (dim_type N=1; N <= 3; ++N) {
    for (dim_type K=1; K < 3; ++K) {
      test_inversion(bgeot::parallelepiped_geotrans(N,K),verbose);
    }
  }
  for (dim_type N=2; N <= 3; ++N) {
    for (dim_type K=1; K < 3; ++K) {
      test_inversion(bgeot::prism_geotrans(N,K),verbose);
    }
  }
  test_inversion(bgeot::parallelepiped_linear_geotrans(5),verbose);
  test_inversion(bgeot::prism_linear_geotrans(3),verbose);
}

int main(int argc, char *argv[]) {
  dim_type N, MESH_TYPE;
  scalar_type LX, LY, LZ;
  size_type NX, NB_POINTS;
  ftool::md_param PARAM;
  getfem::mesh mesh;

  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {
    test0();
    test_inversion(true);
    PARAM.read_command_line(argc, argv);
    N = dim_type(PARAM.int_value("N", "Domaine dimension"));
    NB_POINTS = PARAM.int_value("NB_POINTS", "Nb points");
    LX = PARAM.real_value("LX", "Size in X");
    LY = PARAM.real_value("LY", "Size in Y");
    LZ = PARAM.real_value("LZ", "Size in Y");
    NX = PARAM.int_value("NX", "Nomber of sace steps ");
    MESH_TYPE = dim_type(PARAM.int_value("MESH_TYPE", "Mesh type "));
 
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
    default : GMM_THROW(dal::internal_error, "Unknown type of mesh");
    }
    
    mesh.optimize_structure();


    scalar_type exectime = dal::uclock_sec(), total_time = 0.0;

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
      exectime = dal::uclock_sec();
      if (i==0) { cout << " using points_in_box..\n"; }
      else { cout << " using brute force..\n"; }
      cout << "Time to sort points : " << dal::uclock_sec() - exectime << endl;
      total_time += dal::uclock_sec() - exectime;
      
      dal::bit_vector nn = mesh.convex_index();
      size_type nbtot = 0;
      for (size_type cv = nn.take_first(); cv != size_type(-1); cv << nn) {
	
	size_type nb = gti.points_in_convex(mesh.convex(cv),
					    mesh.trans_of_convex(cv),
					    ptab, itab, i==1);
	//if ((cv % 100) == 0) cout << "cv : " << cv << endl;
	nbtot += nb;
      }
      
      cout << "Time to invert geo trans : " << dal::uclock_sec() - exectime
	   << endl;
      cout << "Total number : " << nbtot << endl;
      assert(nbtot == NB_POINTS);
    }
    total_time += dal::uclock_sec() - exectime;
    
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0;
}
