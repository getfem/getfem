#include <getfem\deformable_mesh.h>

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