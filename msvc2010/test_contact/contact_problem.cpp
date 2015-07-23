/*===========================================================================

 Copyright (C) 2011-2015 Andriy Andreykiv.

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

#include "contact_problem.h"
#include "contact_problem.h"
#include <getfem/getfem_regular_meshes.h>

contact_problem::contact_problem(int argc, char *argv[])
{
	std::cout<<"-- reading parameter file"<<std::endl;
	PARAM.read_command_line(argc,argv);
	tol_newton = PARAM.real_value("RESIDUAL"); 
	model_dim = PARAM.int_value("N","Dimension of the model");
	nstep = PARAM.int_value("NSTEP","Number of steps in the analysis");
	applied_disp = PARAM.real_value("APPLIED_DISP","Final value of the applied displacement");
	ls_offset = PARAM.real_value("LS_OFFSET","Adding a value to all LS nodes");

	//Master contact body
	size_type div_x_M = PARAM.int_value("DIVxM","Mesh division");
	size_type div_y_M = PARAM.int_value("DIVyM","Mesh division");
	size_type div_z_M = PARAM.int_value("DIVzM","Mesh di vision");
	scalar_type x_M = PARAM.real_value("xM","Origin x master");
	scalar_type y_M = PARAM.real_value("yM","Origin y master");
	scalar_type z_M = PARAM.real_value("zM","Origin z master");
	scalar_type L_x_M = PARAM.real_value("LxM","size X master");
	scalar_type L_y_M = PARAM.real_value("LyM","size Y master");
	scalar_type L_z_M = PARAM.real_value("LzM","size Z master");
	app_order_master = PARAM.int_value("APPROX_ORDER_MASTER","Aproximation order master");
	int_order_master = PARAM.int_value("INT_ORDER_MASTER","Intergration order master");
	std::string MESH_TYPE_PREFIX_MASTER = PARAM.string_value("MESH_TYPE_MASTER","mesh type");
	LM_INT_TYPE = PARAM.string_value("LM_INT_TYPE", "integration method for the contact surface");
	mu_master = PARAM.real_value("MU_MASTER", "First Elastic coefficient");
	lambda_master = PARAM.real_value("LAMBDA_MASTER", "Second Elastic coefficient");
	//
	std::cout<<"--generating the master mesh"<<std::endl;
	std::stringstream buf1;
	buf1<<"GT_"<<MESH_TYPE_PREFIX_MASTER<<"("<<model_dim<<","<<app_order_master<<")";
	std::string MESH_TYPE_MASTER=buf1.str();
	std::vector<size_type> nsubdiv_master(model_dim);
	nsubdiv_master[0] = div_x_M; if(model_dim>1) nsubdiv_master[1] = div_y_M; if(model_dim==3) nsubdiv_master[2] = div_z_M; 
	getfem::regular_unit_mesh(mesh_master, nsubdiv_master, bgeot::geometric_trans_descriptor(MESH_TYPE_MASTER));
	base_matrix M_m(model_dim,model_dim);
	bgeot::base_small_vector Origin_master(model_dim);
	Origin_master[0] = x_M;if(model_dim>1) Origin_master[1] = y_M; if(model_dim>2) Origin_master[2] = z_M;
	M_m(0,0) = L_x_M; if(model_dim>1) M_m(1,1) = L_y_M; if(model_dim==3) M_m(2,2) = L_z_M;
	mesh_master.transformation(M_m);
	mesh_master.translation(Origin_master);


	//Slave mesh
	size_type div_x_S = PARAM.int_value("DIVxS","Mesh division");
	size_type div_y_S = PARAM.int_value("DIVyS","Mesh division");
	size_type div_z_S = PARAM.int_value("DIVzS","Mesh di vision");
	scalar_type x_S = PARAM.real_value("xS","Origin x slave");
	scalar_type y_S = PARAM.real_value("yS","Origin y slave");
	scalar_type z_S = PARAM.real_value("zS","Origin z slave");
	scalar_type L_x_S = PARAM.real_value("LxS","size X slave");
	scalar_type L_y_S = PARAM.real_value("LyS","size Y slave");
	scalar_type L_z_S = PARAM.real_value("LzS","size Z slave");
	app_order_slave = PARAM.int_value("APPROX_ORDER_SLAVE","Aproximation order slave");
	int_order_slave = PARAM.int_value("INT_ORDER_SLAVE","Intergration order slave");
	std::string MESH_TYPE_PREFIX_SLAVE = PARAM.string_value("MESH_TYPE_SLAVE","mesh type");
	mu_slave = PARAM.real_value("MU_SLAVE", "First Elastic coefficient");
	lambda_slave = PARAM.real_value("LAMBDA_SLAVE", "Second Elastic coefficient");
	//
	std::cout<<"--generating the slave mesh"<<std::endl;
	std::stringstream buf2;
	buf2<<"GT_"<<MESH_TYPE_PREFIX_SLAVE<<"("<<model_dim<<","<<app_order_slave<<")";
	std::string MESH_TYPE_SLAVE=buf2.str();
	std::vector<size_type> nsubdiv_slave(model_dim);
	nsubdiv_slave[0] = div_x_S; if(model_dim>1) nsubdiv_slave[1] = div_y_S; if(model_dim==3) nsubdiv_slave[2] = div_z_S; 
	getfem::regular_unit_mesh(mesh_slave, nsubdiv_slave, bgeot::geometric_trans_descriptor(MESH_TYPE_SLAVE));
	base_matrix M_s(model_dim,model_dim);
	bgeot::base_small_vector Origin_slave(model_dim);
	Origin_slave[0] = x_S;if(model_dim>1) Origin_slave[1] = y_S; if(model_dim>2) Origin_slave[2] = z_S;
	M_s(0,0) = L_x_S; if(model_dim>1) M_s(1,1) = L_y_S; if(model_dim==3) M_s(2,2) = L_z_S;
	mesh_slave.transformation(M_s);
	mesh_slave.translation(Origin_slave);

	// define boundary regions
	std::cout<<"Building boundary regions for the master"<<std::endl;
	mark_boundary(mesh_master);
	std::cout<<"Building boundary regions for the slave"<<std::endl;
	mark_boundary(mesh_slave);
}

void mark_boundary(getfem::mesh& mesh)
{
	//	Create mesh regions for boundary 
	getfem::mesh_region border_faces;
	getfem::outer_faces_of_mesh(mesh, border_faces); 
	for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
		bgeot::base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
		un/= gmm::vect_norm2(un);
		if (gmm::abs(un[1] - 1.0) < 1.0E-7) { // Top face
			mesh.region(NORTH).add(i.cv(), i.f());
		} 
		if (gmm::abs(un[1] + 1.0) < 1.0E-7) {  // Bottom face
			mesh.region(SOUTH).add(i.cv(), i.f());}
		if (gmm::abs(un[0] - 1.0) < 1.0E-7) { // Right face
			mesh.region(EAST).add(i.cv(), i.f());} 
		if (gmm::abs(un[0] + 1.0) < 1.0E-7) { // Left face
			mesh.region(WEST).add(i.cv(), i.f());}
		if (mesh.dim()==3) 
		{
			if (gmm::abs(un[2] + 1.0) < 1.0E-7) { // front face
				mesh.region(FRONT).add(i.cv(), i.f());}
			if (gmm::abs(un[2] - 1.0) < 1.0E-7) { // back face
				mesh.region(BACK).add(i.cv(), i.f());}
		}
	}
	GMM_ASSERT1(mesh.region(NORTH).index().card()>0,"Region North is empty");
	GMM_ASSERT1(mesh.region(SOUTH).index().card()>0,"Region South is empty");
	GMM_ASSERT1(mesh.region(EAST).index().card()>0, "Region East is empty");
	GMM_ASSERT1(mesh.region(WEST).index().card()>0, "Region West is empty");
	if (mesh.dim()==3){
		GMM_ASSERT1(mesh.region(FRONT).index().card()>0,"Region Front is empty");
		GMM_ASSERT1(mesh.region(BACK).index().card()>0, "Region Back is empty");
	}
}
