/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2012-2012 Andriy Andreykiv
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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



#include <getfem/getfem_deformable_mesh.h>

getfem::deformable_mesh::deformable_mesh(bool _must_be_restored) : mesh(), must_be_restored(_must_be_restored){}

getfem::deformable_mesh::deformable_mesh(
	const getfem::deformable_mesh& _mesh) : 
mesh(), must_be_restored(_mesh.must_be_restored) {
	mesh::copy_from(_mesh);
}


getfem::deformable_mesh& getfem::make_deformable_mesh(const getfem::mesh& m){
    getfem::mesh* pmesh = &(const_cast<getfem::mesh&>(m));
	getfem::deformable_mesh* pm_deformable = dynamic_cast<getfem::deformable_mesh*>(pmesh);	
	GMM_ASSERT1(pm_deformable,"Cannot deform getfem::mesh. Use getfem::deformable_mesh !!!")
	return *pm_deformable;
}
