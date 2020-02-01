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
#include "getfem/getfem_export.h"
#include "getfem/getfem_model_solvers.h"

using bgeot::dim_type;
using bgeot::size_type;
using bgeot::scalar_type;
using bgeot::base_node;

static bool debug=false;

int main(int argc, char *argv[]) {

  gmm::set_traces_level(1);

  bgeot::md_param PARAM;
  PARAM.add_int_param("NX", 1);
  PARAM.add_int_param("NY", 1);
  PARAM.add_int_param("FEM_ORDER", 1);
  PARAM.add_int_param("IM_ORDER", 1);
  PARAM.add_int_param("DIFFICULTY", 0);
  PARAM.read_command_line(argc, argv);
  size_type NX = PARAM.int_value("NX", "Number of elements in X direction");
  size_type NY = PARAM.int_value("NY", "Number of elements in Y direction");
  dim_type FEM_ORDER = dim_type(PARAM.int_value("FEM_ORDER", "Degree of finite element basis"));
  dim_type IM_ORDER = dim_type(PARAM.int_value("IM_ORDER", "Degree of integration method"));
  size_type DIFFICULTY = PARAM.int_value("DIFFICULTY", "Difficulty of test");

  getfem::mesh m;
  getfem::regular_unit_mesh(m, {NX, NY}, bgeot::geometric_trans_descriptor("GT_QK(2, 2)"));

  getfem::mesh_region outer_faces;
  getfem::outer_faces_of_mesh(m, outer_faces);
  m.region(98) = getfem::select_faces_of_normal(m, outer_faces, base_node(-1, 0), 0.001);
  m.region(99) = getfem::select_faces_of_normal(m, outer_faces, base_node(1, 0), 0.001);

  const dim_type N(m.dim());
  getfem::mesh_fem mf(m, N);
  mf.set_classical_finite_element(FEM_ORDER);
  if (DIFFICULTY % 1000 <= 99) {
    dal::bit_vector kept;
    kept.add(0, mf.nb_basic_dof());
    kept.setminus(mf.basic_dof_on_region(98));
    mf.reduce_to_basic_dof(kept);
  }

  getfem::mesh_im mim(m);
  mim.set_integration_method(dim_type(2*IM_ORDER+1));

  getfem::im_data mimd1(mim), mimd3(mim);
  {
    bgeot::multi_index vecsizes(1);
    vecsizes[0] = 3;
    mimd3.set_tensor_size(vecsizes);
  }

  getfem::model md1, md2;
  md1.add_fem_variable("u", mf);
  md2.add_fem_variable("u", mf);
  if (DIFFICULTY % 100 > 19) {
    md1.add_im_variable("eps", mimd3);
    md2.add_internal_im_variable("eps", mimd3);
  } else if (DIFFICULTY % 100 > 9) {
    md1.add_im_variable("eps", mimd3);
    md2.add_im_variable("eps", mimd3);
  }
  md1.add_im_variable("p", mimd1);
  if (DIFFICULTY % 100 > 19 && DIFFICULTY % 100 <= 29)
    md2.add_im_variable("p", mimd1);
  else
    md2.add_internal_im_variable("p", mimd1);


  md1.add_initialized_scalar_data("G", 1);
  md2.add_initialized_scalar_data("G", 1);
  md1.add_initialized_scalar_data("K", 1);
  md2.add_initialized_scalar_data("K", 1);

  std::string exprA, exprB, exprC;
  if (DIFFICULTY % 100 > 9) {
    exprA = "(-p*Id(2)+2*G*([[2/3,0;0,-1/3],[-1/3,0;0,2/3],[0,1;1,0]].eps)):Grad_Test_u";
    exprB = "(eps-Grad_u:[[1,0;0,0],[0,0;0,1],[0,0.5;0.5,0]]).Test_eps";
  } else
    exprA = "(-p*Id(2)+2*G*(Sym(Grad_u)-Div_u*Id(2)/3)):Grad_Test_u";
  if (DIFFICULTY % 10 == 0)
    exprC = "(p+K*Div_u)*Test_p";
  else
    exprC = "((p+1e3*pow(p,3))+K*Div_u)*Test_p";
  getfem::add_nonlinear_generic_assembly_brick(md1, mim, exprA);
  getfem::add_nonlinear_generic_assembly_brick(md2, mim, exprA);
  if (DIFFICULTY % 100 > 9) {
    getfem::add_nonlinear_generic_assembly_brick(md1, mim, exprB);
    getfem::add_nonlinear_generic_assembly_brick(md2, mim, exprB);
  }
  getfem::add_nonlinear_generic_assembly_brick(md1, mim, exprC);
  getfem::add_nonlinear_generic_assembly_brick(md2, mim, exprC);

  if (DIFFICULTY % 1000 > 99) {
    md1.add_filtered_fem_variable("dirmult", mf, 98);
    md2.add_filtered_fem_variable("dirmult", mf, 98);
    getfem::add_linear_generic_assembly_brick(md1, mim, "u.dirmult", 98);
    getfem::add_linear_generic_assembly_brick(md2, mim, "u.dirmult", 98);
  }

  // load
  getfem::add_linear_generic_assembly_brick(md1, mim, "1e-3*Test_u(2)", 99);
  getfem::add_linear_generic_assembly_brick(md2, mim, "1e-3*Test_u(2)", 99);

  std::cout << "Displacement dofs: " << mf.nb_dof() << std::endl;
  std::cout << "Total dofs of model 1: " << md1.nb_dof() << std::endl;
  std::cout << "Total dofs of model 2: " << md2.nb_dof() << std::endl;

  getfem::default_newton_line_search ls1A, ls2A;
  getfem::newton_search_with_step_control ls1B, ls2B;

  std::cout<<"SOLVING MODEL 1 (without internal variables)"<<std::endl;
  gmm::iteration iter(1E-9, 1, 30);
  if (DIFFICULTY % 10000 > 999)
    getfem::standard_solve(md1, iter, getfem::rselect_linear_solver(md1, "mumps"), ls1B);
  else
    getfem::standard_solve(md1, iter, getfem::rselect_linear_solver(md1, "mumps"), ls1A);

  std::cout<<"SOLVING MODEL 2 (with internal variables)"<<std::endl;
  iter.init();
  if (DIFFICULTY % 10000 > 999)
    getfem::standard_solve(md2, iter, getfem::rselect_linear_solver(md2, "mumps"), ls2B);
  else
    getfem::standard_solve(md2, iter, getfem::rselect_linear_solver(md2, "mumps"), ls2A);

  if (debug) {
    std::cout<<std::endl<<"u1:"<<std::endl;
    for (const scalar_type &val : md1.real_variable("u")) std::cout<<val<<std::endl;
    std::cout<<std::endl<<"u2:"<<std::endl;
    for (const scalar_type &val : md2.real_variable("u")) std::cout<<val<<std::endl;

    if (DIFFICULTY % 100 > 9) {
      std::cout<<std::endl<<"eps1:"<<std::endl;
      for (const scalar_type &val : md1.real_variable("eps")) std::cout<<val<<std::endl;
      std::cout<<std::endl<<"eps2:"<<std::endl;
      for (const scalar_type &val : md2.real_variable("eps")) std::cout<<val<<std::endl;
    }

    std::cout<<std::endl<<"p1:"<<std::endl;
    for (const scalar_type &val : md1.real_variable("p")) std::cout<<val<<std::endl;
    std::cout<<std::endl<<"p2:"<<std::endl;
    for (const scalar_type &val : md2.real_variable("p"))
    std::cout<<val<<std::endl;

    getfem::vtk_export exp("test_internal_variable_condensation.vtk", true);
    exp.exporting(mf);
    exp.write_point_data(mf, md1.real_variable("u"), "displ1");
    exp.write_point_data(mf, md2.real_variable("u"), "displ2");
    exp.write_point_data(mf, gmm::sub_vector(md2.real_rhs(true), md2.interval_of_variable("u")), "u residual");
  }

  int ret=0;
  ret += gmm::vect_dist1(md1.real_variable("u"), md2.real_variable("u")) < 1e-9 ? 0 : 1;
  if (DIFFICULTY % 100 > 9)
    ret += gmm::vect_dist1(md1.real_variable("eps"), md2.real_variable("eps")) < 1e-9 ? 0 : 2;
  ret += gmm::vect_dist1(md1.real_variable("p"), md2.real_variable("p")) < 1e-9 ? 0 : 4;

  std::cout<<"Test with difficulty "<<DIFFICULTY<<" returned "<<ret<<std::endl;
  return ret;
}
