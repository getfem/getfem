#pragma once
#include <getfem\deformable_mesh.h>

using getfem::size_type;
using getfem::scalar_type;
using bgeot::base_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

	enum  {NORTH = 1, EAST = 2, WEST = 3, SOUTH = 4, FRONT = 5, BACK = 6};
struct contact_problem{
	getfem::deformable_mesh mesh_master, mesh_slave;
	bgeot::md_param PARAM;
	scalar_type tol_newton,applied_disp;
	size_type model_dim, nstep;
	size_type app_order_master, app_order_slave;
	size_type int_order_master, int_order_slave;
	scalar_type ls_offset;
	scalar_type mu_master, lambda_master, mu_slave, lambda_slave;
	std::string LM_INT_TYPE; 
	contact_problem(int argc, char *argv[]);
	
};
void mark_boundary(getfem::mesh& m);