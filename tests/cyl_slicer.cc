/*===========================================================================

 Copyright (C) 2010-2020 Roman Putanowicz.

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
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/


/**@file cyl_slicer.cc
   @brief Demonstrate function interpolation and cylinder slicing of 2D mesh.
 
   This program tests if slicer_cylinder class is working. In particular
   it check if it is possible to slice 2D mesh with 3D cylinder.

   The result of running this program are two VTK files:
   feminterpolation.vtk -- contains a mesh over rectangular domain
                           and function f(x,y,z) = x
   circleslice.vtk -- contains elliptic slice of the above mesh
*/

#include <getfem/getfem_mesh_slicers.h>
#include <getfem/getfem_mesh.h>
#include <getfem/bgeot_mesh_structure.h>
#include <getfem/getfem_export.h>
#include <getfem/getfem_regular_meshes.h>
#include <getfem/bgeot_config.h>

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;


bgeot::scalar_type func(const bgeot::base_node& x) {
  return x[0];
}

#ifdef GMM_USES_MPI
int main(int argc, char *argv[]) {
  GETFEM_MPI_INIT(argc, argv);
#else
int main(int, char **) {
#endif
  try {
    getfem::mesh mymesh;

    std::vector<getfem::size_type> nsubdiv(2);
    nsubdiv[0] = 8;  // number of mesh elements in X direction
    nsubdiv[1] = 4;  // number of mesh elements in Y direction

    // set pointer to geometric transformation
    bgeot::pgeometric_trans pgt = 
                             bgeot::geometric_trans_descriptor("GT_QK(2,1)");

    getfem::regular_unit_mesh(mymesh, nsubdiv, pgt);

    getfem::base_matrix M;
    int dim = pgt->dim();
    M.resize(dim, dim);
    gmm::clear(M); 
    M(0,0) = 2;
    M(1,1) = 1;
    mymesh.transformation(M);

    getfem::mesh_fem mf(mymesh, 1);
    mf.set_finite_element(getfem::QK_fem(2,1));
    std::vector<bgeot::scalar_type> U(mf.nb_dof());

    getfem::interpolation_function(mf, U, func);

    getfem::vtk_export exp("feminterpolation.vtk");
    exp.exporting(mymesh);
    exp.write_point_data(mf, U, "temperature");

    getfem::stored_mesh_slice sl;
    getfem::base_node x0(0.5, 0.5, -1.0);
    getfem::base_node x1(1.5, 0.5, 1.0);
    sl.build(mymesh, getfem::slicer_cylinder(x0, x1, 0.3, -1), 4);

    getfem::vtk_export expsl("circleslice.vtk");
    expsl.exporting(sl);
    expsl.write_point_data(mf, U, "temperature");

  } GMM_STANDARD_CATCH_ERROR;

  GETFEM_MPI_FINALIZE;

  return 0;
}
