/*===========================================================================

 Copyright (C) 2019-2020 Andriy Andreykiv, Konstantinos Poulios.

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
#include "getfem/getfem_generic_assembly.h"
#include "getfem/getfem_regular_meshes.h"

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  bgeot::md_param PARAM;
  PARAM.add_int_param("N", 10);
  PARAM.add_real_param("RESULT", 0.04478308527);
  PARAM.read_command_line(argc, argv);
  size_t N = PARAM.int_value("N", "Number of elements");
  double RESULT = PARAM.real_value("RESULT", "Expected result");

  getfem::mesh m;
  getfem::regular_unit_mesh(m, {N, N}, bgeot::geometric_trans_descriptor("GT_QK(2, 2)"));

  getfem::mesh_fem mf(m);
  mf.set_classical_finite_element(2);
  getfem::mesh_im mim(m);
  mim.set_integration_method(3);
  bgeot::multi_index scalarSize(0);
  getfem::im_data imd(mim, scalarSize, -1);
  std::vector<double> v(mf.nb_dof()), u(mf.nb_dof()), d(imd.nb_index());

  getfem::ga_workspace w;
  w.set_assembled_vector(v);
  w.add_fem_variable("u", mf, gmm::sub_interval(0, mf.nb_dof()), u);
  w.add_im_data("d", imd, d);
  w.add_assignment_expression("d", "Norm(X)", -1, 1, true);
  w.add_expression("d * Test_u", mim, -1);
  w.assembly(1);
  return (gmm::abs(gmm::vect_norm2(v) - RESULT) < 1e-10) ? 0 : 1;
}
