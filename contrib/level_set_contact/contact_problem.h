/* -*- c++ -*- (enables emacs c++ mode) */
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

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

#pragma once
#pragma once
#include <getfem/getfem_deformable_mesh.h>
#include <getfem/getfem_models.h>

using getfem::size_type;
using getfem::scalar_type;
using bgeot::base_matrix;
typedef getfem::model_real_plain_vector  plain_vector;

enum  { NORTH = 1, EAST = 2, WEST = 3, SOUTH = 4, FRONT = 5, BACK = 6};

struct contact_problem{
	getfem::mesh mesh_master, mesh_slave;
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
