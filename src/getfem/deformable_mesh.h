#pragma once
#include <getfem\getfem_mesh.h>
#include <getfem/getfem_mesh_fem.h>
#include <getfem\getfem_modeling.h>

namespace getfem
{

	template<class VECTOR> class temporary_mesh_deformator;


   /** This is a normal mesh, whith one extra method, allowing to displace points.       
	   The mesh can only be deformed by instance of class temporary_mesh_deformator, 
	   that restores the mesh on it's (deformator) destruction
  */
	class deformable_mesh : public mesh {
	public:
		mutable bool must_be_restored;
	private:
		template <class VECTOR> friend class temporary_mesh_deformator;

		/** says that if the mesh was deformed, it should be deformed back to the 
		    underformed state, as other bricks don't know they are dealing with a deformed mesh
			This mesh is used in Updated Lagrane based formulations, but the restore feature
			allows to use it with Total Lagrange as well*/
		inline bool to_be_restored() const {return must_be_restored;}

		/**displace the points by a given displacement vector
		@param U displacement vector as described by mf using dof index (NOT pts index)
		@param &mf mesh_fem object that corresponds to &U, should be compatible with the mesh
		*/
		template<typename VEC>
		void deform_mesh(VEC &dU, const mesh_fem& mf)
		{   
			PT_TAB& pts = points();
			size_type dim = pts.dim();

			GMM_ASSERT1((&mf.linked_mesh())==this,"in deform_mesh mf should be defined on the same mesh");

			GMM_ASSERT1(mf.get_qdim() == dim, 
				"input mesh_fem and the mesh dim are not compatible");
			GMM_ASSERT1(mf.nb_dof() == this->nb_points()*dim,
				"mesh_fem should be isoparametric to the mesh, with qdim == mesh dim");
			dal::bit_vector conv_indices = mf.convex_index(); 
			//this vector will track if a point can be deformed
			std::vector<bool> deform_pt_flag(pts.size(), true);
			size_type cv;
			for(cv << conv_indices; 
				cv!=bgeot::size_type(-1); cv << conv_indices) 
			{
				getfem::mesh::ind_cv_ct pt_index
					=  mf.linked_mesh().ind_points_of_convex(cv);
				getfem::mesh_fem::ind_dof_ct dof=mf.ind_dof_of_element(cv);
				bgeot::size_type num_points = 
					mf.linked_mesh().structure_of_convex(cv)->nb_points(); 
				for(size_type pt = 0; pt < num_points; ++pt) 
				{ 
					/** iterate through each components of point [pt]and deform the component*/
					if(deform_pt_flag[pt_index[pt]])
						for (size_type comp = 0; comp < dim; ++comp)
							//move pts by dU;
							pts[pt_index[pt]][comp] += dU[dof[pt*dim + comp]];

					//flag current [pt] to deformed
					deform_pt_flag[pt_index[pt]] = false;
				}
			}
		}

	public:

		deformable_mesh(bool _must_be_restored = true);
		deformable_mesh(const deformable_mesh&);
	};

	/**cast a conventional mesh into deformable one and remove the const*/
	deformable_mesh& make_deformable_mesh(const mesh&);


	/** a class that first deformes and then remembers to restore a deformable mesh 
	    if it has to be restored for other bricks*/
	template<class VECTOR = model_real_plain_vector> class temporary_mesh_deformator{
		VECTOR dU;
		const mesh_fem& mf;
		deformable_mesh& m;
	public:
		temporary_mesh_deformator(const mesh& _m, const mesh_fem &_mf, const VECTOR &_dU) : 
		  m(make_deformable_mesh(_m)), 
		  mf(_mf), 
		  dU(_dU) 
		  {m.deform_mesh(dU,mf);}

		  ~temporary_mesh_deformator(){
			  if (m.to_be_restored()){
				  m.deform_mesh(gmm::scaled(dU,scalar_type(-1.0)),mf); }
		  }
	};

}//end of getfem namespace