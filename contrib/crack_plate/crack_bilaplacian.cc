// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2006-2007 Yves Renard, Julien Pommier.
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

#include "crack_bilaplacian.h"

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::base_vector; /* dense vector. */
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* dense matrix. */

/* definition of some matrix/vector types. 
 * default types of getfem_model_solvers.h
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;


int main(int argc, char *argv[]) {

  try {
    bilaplacian_crack_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    plain_vector U;
    p.mesh.write_to_file("mesh.m") ;
    if (!p.solve(U)) GMM_ASSERT1(false, "Solve has failed");

    p.compute_error(U);
    p.compute_H2_error_field(U) ;

    // visualisation, export au format .vtk
    getfem::mesh mcut;
    p.mls.global_cut_mesh(mcut);
    unsigned Q = p.mf_u().get_qdim();
    assert( Q == 1 ) ;
    getfem::mesh_fem mf(mcut, Q);
    mf.set_classical_discontinuous_finite_element(2, 0.001);
    // mf.set_finite_element
    //	(getfem::fem_descriptor("FEM_PK_DISCONTINUOUS(2, 2, 0.0001)"));
    plain_vector V(mf.nb_dof());

    getfem::interpolation(p.mf_u(), mf, U, V);

    getfem::stored_mesh_slice sl;
    getfem::mesh mcut_refined;

    unsigned NX = p.PARAM.int_value("NX"), nn;
    if (NX < 6) nn = 12;
    else if (NX < 12) nn = 8;
    else if (NX < 30) nn = 3;
    else nn = 1;

    /* choose an adequate slice refinement based on the distance to the crack tip */
    std::vector<bgeot::short_type> nrefine(mcut.convex_index().last_true()+1);
    for (dal::bv_visitor cv(mcut.convex_index()); !cv.finished(); ++cv) {
      scalar_type dmin=0, d;
      base_node Pmin,P;
      for (unsigned i=0; i < mcut.nb_points_of_convex(cv); ++i) {
	P = mcut.points_of_convex(cv)[i];
	//d = gmm::vect_norm2(ls_function(P));
	d = gmm::vect_norm2(P) ; 
	if (d < dmin || i == 0) { dmin = d; Pmin = P; }
      }

      if (dmin < 1e-5)
	nrefine[cv] = nn*2;
      else if (dmin < .1) 
	nrefine[cv] = nn*2;
      else nrefine[cv] = nn;
      /*if (dmin < .01)
	cout << "cv: "<< cv << ", dmin = " << dmin << "Pmin=" << Pmin << " " << nrefine[cv] << "\n";*/
    }


    getfem::mesh_slicer slicer(mcut); 
    getfem::slicer_build_mesh bmesh(mcut_refined);
    slicer.push_back_action(bmesh);
    slicer.exec(nrefine, getfem::mesh_region::all_convexes());
     
    /*
      sl.build(mcut, 
      getfem::slicer_build_mesh(mcut_refined), nrefine);*/

    getfem::mesh_im mim_refined(mcut_refined); 
    mim_refined.set_integration_method(getfem::int_method_descriptor
				       ("IM_TRIANGLE(6)"));

    getfem::mesh_fem mf_refined(mcut_refined, Q);
    mf_refined.set_classical_discontinuous_finite_element(2, 0.001);
    plain_vector W(mf_refined.nb_dof());

    
    getfem::interpolation(p.mf_u(), mf_refined, U, W);


    int VTK_EXPORT = p.PARAM.int_value("VTK_EXPORT");
    if (VTK_EXPORT ) {
      cout << "exporting solution to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk", false);
      exp.exporting(mf_refined); 
      exp.write_point_data(mf_refined, W, "vertical_displacement");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d " << p.datafilename  << ".vtk -f "
	"WarpScalar -m BandedSurfaceMap -m Outline\n";
    }
    if ((VTK_EXPORT & (2|4))) {
      p.exact_sol.mf.set_qdim(Q);
      assert(p.exact_sol.mf.nb_dof() == p.exact_sol.U.size());   
      plain_vector EXACT(mf_refined.nb_dof());
      getfem::interpolation(p.exact_sol.mf, mf_refined, 
			    p.exact_sol.U, EXACT);
	
      plain_vector DIFF(EXACT); gmm::add(gmm::scaled(W,-1),DIFF);
      if ((VTK_EXPORT & 2)) {
	cout << "exporting exact solution to VTK\n";
	getfem::vtk_export exp2(p.datafilename + "_exact.vtk");
	exp2.exporting(mf_refined);
	exp2.write_point_data(mf_refined, EXACT, "reference solution");
	cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d " << p.datafilename << "_exact.vtk -f "
	"WarpScalar -m BandedSurfaceMap -m Outline\n";
      }
      if ((VTK_EXPORT & 4)) {
	cout << "exporting difference with exact solution to VTK\n";
	getfem::vtk_export exp2(p.datafilename + "_diff.vtk");
	exp2.exporting(mf_refined);
	exp2.write_point_data(mf_refined, DIFF, "difference solution");
	cout << "export done, you can view the data file with (for example)\n"
	  "mayavi -d " << p.datafilename  << "_diff.vtk -f "
	"WarpScalar -m BandedSurfaceMap -m Outline\n";
      }
    }

    int MATLAB_EXPORT = p.PARAM.int_value("MATLAB_EXPORT");
    if (MATLAB_EXPORT) {
      cout << "exporting solution to " << p.datafilename + ".mf" << " and " << p.datafilename << ".U\n";
      mf_refined.write_to_file(p.datafilename + ".mf", true);
      vecsave(p.datafilename + ".U", W);
      p.exact_sol.mf.set_qdim(Q);
      assert(p.exact_sol.mf.nb_dof() == p.exact_sol.U.size());   // ??
      plain_vector EXACT(mf_refined.nb_dof());
      getfem::interpolation(p.exact_sol.mf, mf_refined, 
			    p.exact_sol.U, EXACT);
      vecsave(p.datafilename + ".EXACT", EXACT);

      cout << "exporting original mesh to " << p.datafilename + "_base.mesh\n";
      p.mesh.write_to_file(p.datafilename + "_base.mesh");

      cout << "matlab export done, you can view the results with\n";
      cout << "mf=gfMeshFem('load', 'bilaplacian.mf'); U=load('bilaplacian.U')'; "
	"EXACT=load('bilaplacian.EXACT')'; m0=gfMesh('load','bilaplacian_base.mesh');\n";
      cout << "gf_plot(mf, U-EXACT, 'refine',1); hold on; gf_plot_mesh(m0); hold off; colorbar;\n";
    }

    
    //getchar(); 
    
  }
  GMM_STANDARD_CATCH_ERROR;
  return 0; 
}
























