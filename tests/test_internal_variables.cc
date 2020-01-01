/*===========================================================================

 Copyright (C) 2019-2020 Konstantinos Poulios.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_generic_assembly.h"

using bgeot::dim_type;
using bgeot::size_type;
using bgeot::scalar_type;
using bgeot::base_node;

int main(int argc, char *argv[]) {

//  gmm::set_traces_level(1);

  bgeot::md_param PARAM;
  PARAM.add_int_param("NX", 3);
  PARAM.add_int_param("NY", 2);
  PARAM.add_int_param("FEM_ORDER", 1);
  PARAM.add_int_param("IM_ORDER", 1);
  PARAM.add_int_param("DIFFICULTY", 0);
  PARAM.read_command_line(argc, argv);
  size_type NX = PARAM.int_value("NX", "Number of elements in X direction");
  size_type NY = PARAM.int_value("NY", "Number of elements in Y direction");
  dim_type FEM_ORDER = dim_type(PARAM.int_value("FEM_ORDER", "Degree of finite element basis"));
  dim_type IM_ORDER = dim_type(PARAM.int_value("IM_ORDER", "Degree of integration method"));
  size_type DIFFICULTY = PARAM.int_value("DIFFICULTY", "Difficulty of test (0 or 1)");

  getfem::mesh m;
  getfem::regular_unit_mesh(m, {NX, NY}, bgeot::geometric_trans_descriptor("GT_QK(2, 2)"));

  getfem::mesh_region outer_faces;
  getfem::outer_faces_of_mesh(m, outer_faces);
  m.region(100) = getfem::select_faces_of_normal(m, outer_faces, base_node(-1, 0), 0.001);
  m.region(101) = getfem::select_faces_of_normal(m, outer_faces, base_node(1, 0), 0.001);
  m.region(102) = getfem::mesh_region::merge(m.region(100), m.region(101));

  dim_type N(2);
  getfem::mesh_fem mf(m, N), mf_intern(m);
  mf.set_classical_finite_element(FEM_ORDER);
  if (DIFFICULTY) mf_intern.set_qdim(3,4);
  mf_intern.set_classical_discontinuous_finite_element(IM_ORDER);

  getfem::mesh_im mim(m);
  mim.set_integration_method(dim_type(2*IM_ORDER+1));

  getfem::im_data mimd(mim);
  if (DIFFICULTY) mimd.set_tensor_size(bgeot::multi_index(3,4));

  getfem::model md1, md2;
  md1.add_fem_variable("u", mf);
  md2.add_fem_variable("u", mf);
  md1.add_im_variable("p", mimd);
  md2.add_fem_variable("p", mf_intern);

  md1.add_initialized_scalar_data("G", 1);
  md2.add_initialized_scalar_data("G", 1);
  md1.add_initialized_scalar_data("K", 1);
  md2.add_initialized_scalar_data("K", 1);

  std::string exprA, exprB;
  if (DIFFICULTY) {
    exprA = "(-1e3*asin(p(1,1))*Id(2)+2*G*(Sym(Grad_u)-Div_u*Id(2)/3)):Grad_Test_u";
    exprB = "(p+sin(0.001*K*Trace(Sym(Grad_u)))*[1,1,1,1;1,1,1,1;1,1,1,1]):Test_p";
  } else {
    exprA = "(-p*Id(2)+2*G*(Sym(Grad_u)-Div_u*Id(2)/3)):Grad_Test_u";
    exprB = "(p+K*Trace(Sym(Grad_u)))*Test_p";
  }
  getfem::add_nonlinear_generic_assembly_brick(md1, mim, exprA);
  getfem::add_nonlinear_generic_assembly_brick(md2, mim, exprA);
  getfem::add_nonlinear_generic_assembly_brick(md1, mim, exprB);
  getfem::add_nonlinear_generic_assembly_brick(md2, mim, exprB);

  md1.add_filtered_fem_variable("dirmult", mf, 102);
  md2.add_filtered_fem_variable("dirmult", mf, 102);
  getfem::add_linear_generic_assembly_brick(md1, mim, "(u-0.001*X(1)*[1;0]).dirmult", 102);
  getfem::add_linear_generic_assembly_brick(md2, mim, "(u-0.001*X(1)*[1;0]).dirmult", 102);

  gmm::iteration iter(1E-9, 1, 100);
  getfem::standard_solve(md1, iter);
  iter.init();
  getfem::standard_solve(md2, iter);

  for (const scalar_type &val : md1.real_variable("u"))
  std::cout<<val<<std::endl;

  std::cout<<std::endl;
  for (const scalar_type &val : md2.real_variable("u"))
  std::cout<<val<<std::endl;

  std::cout << "Displacement dofs: " << mf.nb_dof() << std::endl;
  std::cout << "Total dofs of model 1: " << md1.nb_dof() << std::endl;
  std::cout << "Total dofs of model 2: " << md2.nb_dof() << std::endl;

  return gmm::vect_dist2(md1.real_variable("u"), md2.real_variable("u")) < 1e-9 ? 0 : 1;
}
