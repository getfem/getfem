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
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/

/** @file getfem_level_set_contact.h
    @author "Andriy Andreykiv" <andriy.andreykiv@gmail.com>
    @date July, 2012.
    @brief Non frictional level set based large sliding contact;
      for details see: 
      A. Andreykiv et al. A level set based large sliding contact
      algorithm for an easy analysis of implant positioning 
      2012 International Journal for Numerical Methods in Engineering, 
      89, pp. 1317-1336
      2D and 3D Examples of the usage: test_contact.cpp
 */


#pragma once
// #include <string>
// #include <memory>
// #include <map>
#include <getfem/getfem_models.h>
#include <getfem/getfem_model_solvers.h>
#include <getfem/getfem_deformable_mesh.h>
#include <gmm/gmm_except.h>

namespace level_set_contact {

	using getfem::mesh_fem;
	using getfem::mesh_im;
	using getfem::mesh;
	using getfem::model;
	using getfem::size_type;
	using getfem::scalar_type;
	using getfem::modeling_standard_plain_vector;
	typedef getfem::modeling_standard_plain_vector  plain_vector;
	typedef getfem::model_real_sparse_matrix sparse_matrix;


	/**build a level set function on mesh with zero on the boundary.
	Solves Laplace equation with zero Dirichlet on the boundary.
	Used to create simple level sets for contact testing*/
	template<class VECT> void boundary_level_set_field(
		const getfem::mesh& _mesh, 
		const getfem::mesh_fem& mf,
		const getfem::mesh_im& mim, 
		VECT& LS)
	{
		getfem::mesh& mesh = const_cast<getfem::mesh&>(_mesh);
		//model and vars
		getfem::model md;
		md.add_fem_variable("LS",mf);
		getfem::modeling_standard_plain_vector RHS(mf.nb_dof());

        //calculating the size of the LS RHS based on the size of the geometry
		getfem::base_node Pmin(mesh.dim()),Pmax(mesh.dim()),range(mesh.dim());
		mesh.bounding_box(Pmin,Pmax);
		gmm::add(Pmax,gmm::scaled(Pmin,-1.0),range);
		getfem::scalar_type mesh_size = *(std::max_element(range.begin(),range.end()));
		gmm::fill(RHS,mesh_size*5.0);
		md.add_initialized_fem_data("RHS",mf,RHS);

		//border region
		getfem::mesh_region border_faces;
		getfem::outer_faces_of_mesh(mesh, border_faces); 
    bgeot::size_type BORDER=getfem::mesh_region::free_region_id(mesh);
		for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) mesh.region(BORDER).add(i.cv(),i.f());

		//describing the PDE problem
		getfem::add_Laplacian_brick(md,mim,"LS");
		getfem::add_Dirichlet_condition_with_penalization(md,mim,"LS",1e9,BORDER);
		getfem::add_source_term_brick(md,mim,"LS","RHS");

		//solving
		gmm::iteration iter;
		GMM_TRACE2("building scalar level set with laplace equation..");
		getfem::standard_solve(md,iter);

		//extracting the result
		gmm::copy(md.real_variable("LS"),LS);

		//so, now the mesh is as it was, hence const is still valid
		mesh.sup_region(BORDER);
		GMM_TRACE2("..done")
	}


	/**base class for the master and the slave contact bodies.*/
	class contact_body{

		const std::string var_name;
		bool is_deformed;
		friend class contact_pair_update;

	protected:
		mesh& own_mesh; 
		const mesh_fem& own_mesh_fem;
		model& md;

	public:

		contact_body(model& _md, std::string _var_name);
		inline std::string get_var_name() const {return var_name;}
		inline mesh& get_mesh() {return own_mesh;}
		inline const mesh& get_mesh()const {return own_mesh;}
		inline const mesh_fem& get_mesh_fem() const {return own_mesh_fem;}
		inline const model& get_model() const {return md;}
		inline bool is_mesh_deformed() const {return is_deformed;}
	};


	/** Contact body that will be projected on the boundary 
	of the master. */
	class slave_contact_body: public contact_body {

		std::string ls_name;
		mesh_fem ls_mesh_fem;
		mesh_im* pmim;

	public:

		/**default constructor. Level set field will have zero value 
		right on the boundary of the contact body*/
		slave_contact_body(model& _md, const std::string& _var_name, 
			mesh_im* _pmim);

		/**Level set field is provided via the model variable name*/
		slave_contact_body(model& _md, std::string _var_name, 
			std::string _ls_name);
		inline std::string get_ls_name() const {return ls_name;}
		inline const plain_vector& ls_values() const 
		{return md.real_variable(ls_name);}
		inline plain_vector& ls_values() 
		{return md.set_real_variable(ls_name);}
		inline const mesh_fem& get_ls_mesh_fem() const {return md.mesh_fem_of_variable(ls_name);}
		template<class VECTOR> void set_level_set(const VECTOR& ls)
		{gmm::copy(ls,md.set_real_variable(ls_name));}

		/**adds a fixed value "off" to the level set field */
		void offset_level_set(scalar_type off);
	};


	class master_contact_body;

	/**Prepares the final information needed to pass to the contact 
	brick for every contact pair to assemble tangent terms*/
	class contact_pair_info {

		//obtain on construction
		master_contact_body& master_cb;
		slave_contact_body& slave_cb;
		const std::string mult_name;
		const size_type GIVEN_CONTACT_REGION;

		//to be built
		dal::bit_vector old_contact_elm_list;
		dal::bit_vector pre_old_ct_list;
		size_type ACTIVE_CONTACT_REGION;
		mutable dal::shared_ptr<mesh_im> pmim_contact;
		mutable getfem::pfem ifem_srf;
		mutable dal::shared_ptr<mesh_fem> pinterpolated_fem;
		mutable dal::shared_ptr<mesh_fem> pinterpolated_fem_U;
		mutable dal::shared_ptr<gmm::unsorted_sub_index> slave_ls_dofs;
		mutable dal::shared_ptr<gmm::unsorted_sub_index> slave_U_dofs;
		mutable size_type n_integrated_elems;

		// state of the object
		mutable bool members_are_computed;
		mutable bool init_cont_detect_done;
	public:

		// accessors
		inline const mesh_fem& slave_scalar_fem() const {
			if (dependecies_changed()) update();
			return *pinterpolated_fem;
		}
		inline const mesh_fem& slave_vector_fem() const  {
			if (dependecies_changed()) update();
			return *pinterpolated_fem_U;
		}
		inline const gmm::unsorted_sub_index& slave_scalar_dofs() const {
			if (dependecies_changed()) update();
			return *slave_ls_dofs;
		}
		inline const gmm::unsorted_sub_index& slave_vector_dofs() const {
			if (dependecies_changed()) update();				
			return *slave_U_dofs;
		}
		inline const mesh_im& contact_mesh_im() const {
			if (dependecies_changed()) update();
			return *pmim_contact;
		}

		inline size_type contact_region() const 
		{return ACTIVE_CONTACT_REGION;}

		inline const std::string& get_mult_name() const 
		{return mult_name;}

		inline size_type num_of_integr_elems() const {return n_integrated_elems;}
		// update
		inline bool dependecies_changed() const
		{return !members_are_computed;}
		inline void force_update() const 
		{members_are_computed=false;}

		/** Actual master/slave contact detection. Level set field is projected on the 
		boundary of the master and only the elements which nodes satisfy
		level_set + Multiplier > 0 
		become contact elements*/
		bool contact_changed();

		/**clearing contact element lists*/
		void clear_contact_history();

		/** updating contact information (mesh_fem's, mesh_im's)
		with the deformation. Contact detection is not performed*/
		void update(void) const;

		contact_pair_info(master_contact_body& underformed_mcb, 
			slave_contact_body& underformed_scb, const std::string& _mult_name, 
			size_type _GIVEN_CONTACT_REGION);

	private:
		/**prohibiting copying*/
		contact_pair_info(const contact_pair_info&);
		contact_pair_info& operator=(const contact_pair_info&);


	};

	struct face_type{
		size_type cv,f;
		face_type(size_type _cv=0, size_type _f=0):cv(_cv),f(_f){}
		face_type(const getfem::mr_visitor& i): cv(i.cv()),f(i.f()){}
	};

	/**Determines geometric transformation on the face of the element
	based on the geometric transformation of the element itself. Works
	only for PK and QK elements*/
	bgeot::pgeometric_trans face_trans_of_elem(bgeot::pgeometric_trans pelem_trans);


	/** Master contact body which surface will be used to project contact 
	stresses and stiffness terms. It contains and manages the slaves and 
	knows other masters.
	Master contact body must be created with mesh_fem that allows automatic 
	addition of mesh_fem description on new elements (use mesh_fem::set_auto_add or
	use set_classical_finite_element). This feature is used when new boundary 
	elements are created from faces. At the same time the mesh_im object that 
	is used to add for instance some structural bricks on the volume (elastostatic, 
	nonlinear_elastostatic, updated_lagrangian) should be either created before 
	master contact body, or set on master_contact_body::volume_region() if it's 
	created after. This is to avoid integration of the volume integrals on the 
	boundary elemenents of lower dimension. */
	class master_contact_body: public contact_body {


		const size_type mult_mim_order;
		const std::string mult_int_method;
		size_type BOUNDARY_ELEMENTS, VOLUME_ELEMENTS;
		std::vector<size_type> face_to_belem_ind;
		static std::vector<master_contact_body*> masters;
		std::map<std::string, dal::shared_ptr<contact_pair_info> > contact_table;
		std::map<size_type,face_type> border_faces;

	protected:

		/**contact detection for all slaves*/
		bool master_contact_changed(void);

		/** clearing previous contact elem lists*/
		void clear_contact_history(void);

	public:
	
		enum contact_integration{PER_ELEMENT=1,REGULARIZED_LEVEL_SET=2};

		/**approximation order for Lagrange multiplier on the contact surface*/
		const size_type mult_mf_order;

		/**integration approach for contact elements that are partially 
		crossed by level sets:
		PER_ELEMENT - a whole element is incuded into contact (default)
        REGULARIZED_LEVEL_SET - Gauss points where projected value of the level
		                        set is < zero are set to zero or small value 
								(with gradual transition)*/
		const contact_integration integration;

		/**width of transition for a regularazied Heaviside function in 
		case of REGULARIZED_LEVEL_SET*/
		const scalar_type regularized_tollerance;

		/**in case of REGULARIZED_LEVEL_SET this value
		scales weight of Gauss points that have negative level 
		set value*/
		const scalar_type small_weight_multiplier;

		/**if the angle (in degrees) between contact element and 
		level set contour exceed this value, this element is not included in 
		contact algorithm*/
		const scalar_type max_contact_angle;


		/** create master contact body with a model,
		name where masters displacements are defined, order for 
		Lagrange multiplier, order for the integration method*/
		master_contact_body(model& _md,
			const std::string& _var_name, 
			size_type _mult_order, size_type _mult_mim_order);

		/**the same as above, but specifically provide itegration 
		* method on the contact surface (_mult_int_method), additionally, 
		* specify if surface contact elements have to be cut by the level set.
		* The later ensures that contact surface is strictly a domain
		* that overlaps with the slave, hence this allows smooth growth of the contact
		* surface. The level set cutting is done using regularized Heaviside function
		*/
		master_contact_body(model& _md,
			const std::string& _var_name, 
			size_type _mult_order,
			const std::string& _mult_int_method,
			contact_integration _integration = PER_ELEMENT, 
			scalar_type _regularized_tollerance = 1e-6,
			scalar_type _small_weight_multiplier = 0.001,
			scalar_type _max_contact_angle = 45);

		/** associate a slave contact body with this master. \
		specify a region of expected contact interaction.  \
		(takes the whole master boundary if not specified)*/
		void add_slave(slave_contact_body& scb, 
			size_type slave_contact_region = -1);

		/** order of integration of boundary contact terms*/
		inline size_type contact_mim_order() const 
		{
		  GMM_ASSERT1(mult_mim_order!=size_type(-1),
			      "master body was not created with "					      "order of integration for contact area");
		  return mult_mim_order;
		}

		/** integration method on the contact surface, 
		* use it when the master is created with a specific
		* integration method and not the approx_order*/
		inline getfem::pintegration_method contact_int_method() const
	        {
		  GMM_ASSERT1(mult_mim_order==size_type(-1),
			      "master body was not created with integration "
			      "method for contact area");
		  return getfem::int_method_descriptor(mult_int_method);
		}

		/** region of all volume elements without the boundary*/
		inline size_type volume_region() const 
		{return VOLUME_ELEMENTS;}

		/**boundary elements, added after creation of 
		the master contact body */
		inline size_type boundary_region() const 
		{return BOUNDARY_ELEMENTS;}

		/**access to a structure that contains all the info 
		about contact pair between this master and a slave, defined 
		on @param slave_var_name*/
		const contact_pair_info& get_pair_info(
			const std::string& slave_var_name) const;

		/**the same as above, but non-const*/
		contact_pair_info& get_pair_info(
			const std::string& slave_var_name);


		/**contact detection for all masters/slave couples
		@return true if any of the contact areas where changed 
		(which requires new Newton method run)*/
		static bool any_contact_change();

		/** should be used in the beginning of a step
		to clean data structures that store previous
		contact element lists (used to verify if contact surface 
		is converged to one list)
		*/
		static void clear_all_contact_history();

		inline void update_for_slave(std::string slave_var_name)
		{contact_table[slave_var_name]->update();};

		/** return a pointer to mesh_im used for contact surface calculations
		*/
		dal::shared_ptr<mesh_im> build_mesh_im_on_boundary(
			size_type region);

		/**gives a face, corresponding to newly created 
		boundary element @param cv*/
		face_type ext_face_of_elem(size_type cv) const;

	private:
		/**prohibiting copying*/
		master_contact_body(const master_contact_body&);
		master_contact_body& operator=(const master_contact_body&);

	};

	enum update_depth{DEFORM_MESHES_ONLY,FULL_UPDATE};

	/**temporary object that updates contact pair, 
	deformes meshes and undeformes when it selfdestructs*/
	class contact_pair_update{
		dal::shared_ptr<getfem::temporary_mesh_deformator<> > def_master;
		dal::shared_ptr<getfem::temporary_mesh_deformator<> > def_slave;
		master_contact_body& mcb;
		slave_contact_body& scb;
	public:
		contact_pair_update(master_contact_body& _mcb,
				    slave_contact_body& _scb,
				    update_depth ud = FULL_UPDATE);

		~contact_pair_update();
	};


	/** adding level set based normal contact brick to the model.
	The contact is etablished between the 
	@param mcb - master contact body and 
	@param scb - slave contact body, defined on 
	@param md  - model object 
	@param rg  - optional assumed contact region 
	helping to narrow down contact search
	Note, this contact algorithm is note stabilized, hence,
	master contact body mesh should be coarser than slave's mesh.
	Otherwise this contact constraint will violate inf-sub condition
	and the solver will fail (or diverge, if it's iterative)
	*/
	size_type add_level_set_normal_contact_brick(model& md, 
		master_contact_body& mcb, 
		slave_contact_body& scb,
		size_type rg = -1);


	/** assembles normal contact terms on the boundary of 
	two contact bodies (master/slave)*/
	class level_set_contact_brick: public getfem::virtual_brick{

		model& md;
		master_contact_body& mcb;
		slave_contact_body& scb;

		/**id of the region of faces where contact has to be checked*/
		size_type given_contact_id;

		/**id of the region of boundary elements, 
		corresponding to the above faces*/
		size_type contact_region_id;

		/**actual region object,  with id = contact_region_id*/
		getfem::mesh_region contact_region;

	public:
		virtual void asm_real_tangent_terms(
			const model &md, size_type /* ib */,
			const model::varnamelist &vl,
			const model::varnamelist &dl,
			const model::mimlist &mims,
			model::real_matlist &matl,
			model::real_veclist &vecl,
			model::real_veclist &,
			size_type region,
			build_version version) const;			

		level_set_contact_brick(
			model& _md, 
			master_contact_body& _mcb, 
			slave_contact_body& scb, 
			size_type rg = -1);
	};


	/** A term, used in level set contact assemblies that 
	builds a surface projection matrix R = N^t X N 
	(where N is normal vector to the boundary)
	*/
	class NormalTerm : public getfem::nonlinear_elem_term
	{	private:
	const master_contact_body& mcb;
	bgeot::multi_index sizes_;
	bgeot::size_type version;
	bgeot::size_type dim;

	public:

		NormalTerm(const master_contact_body& _mcb, size_type version_ = 1) :
		  mcb(_mcb),
		  sizes_(version_),
		  version(version_),
		  dim(_mcb.get_mesh().dim()) {

				  GMM_ASSERT1(dim==2 || dim==3, "NormalTerm: wrong space dimension ");
				  GMM_ASSERT1(version==1 || version==2,"NormalTerm:: wrong version ");

				  if (version == 1)
					  if (dim == 2)
						  sizes_[0] = 2;
					  else
						  sizes_[0] = 3;
				  else 
					  if (dim == 2) {
						  sizes_[0] = 2;
						  sizes_[1] = 2; 
					  }
					  else {
						  sizes_[0] = 3;
						  sizes_[1] = 3;
					  } 
		  }
		  const bgeot::multi_index &sizes(size_type) const {return sizes_;};
		  void compute(getfem::fem_interpolation_context& ctx, bgeot::base_tensor &t);
	          void prepare(getfem::fem_interpolation_context& /* ctx */, size_type /* nl_part */) {}

	};

	/** Regularized Heaviside function.
	Can be used instead of mesh_im_level_set in assemblies.
	It's more stable, as it never fails in comparison to Delauney method
	(used inside mesh_im_level_set), but less accurate, as it has a 
	transition zone from 1 to 0 of epsilon width.
	The idea is taken from one of the articles of Ted Belytschko on XFem*/
	class HFunction : public getfem::nonlinear_elem_term
	{
	private:
		const mesh_fem &lsmf;
		const plain_vector &LS_U;
		scalar_type  m_Epsilon;
		scalar_type small_h;
		bgeot::multi_index sizes_;


	public:
		HFunction(
			const mesh_fem &lsmf_,
			const plain_vector &LS_U_,
			scalar_type epsilon=1e-9, 
			scalar_type small_h_=0);
		  const bgeot::multi_index &sizes(size_type) const;
		  void prepare(getfem::fem_interpolation_context& ctx, size_type nl_part);
		  void compute(getfem::fem_interpolation_context& ctx, bgeot::base_tensor &t);
		  scalar_type hRegularized(scalar_type x, scalar_type epsion, scalar_type small);
	};

	//A dummy nonlinear term, does nothing
	class Unity : public getfem::nonlinear_elem_term
	{
	private:
		const mesh_fem &mf;
		bgeot::multi_index sizes_;

	public:
		Unity(const mesh_fem &mf_);
		const bgeot::multi_index &sizes(size_type) const;
		void prepare(getfem::fem_interpolation_context& ctx, size_type nl_part);
		void compute(getfem::fem_interpolation_context& ctx, bgeot::base_tensor &t);
	};



	template<typename MAT, typename VECT> 
	void asm_level_set_contact_tangent_matrix(
		std::vector<MAT>& matl, 
		const master_contact_body& mcb, 
		const slave_contact_body& scb, 
		const VECT& LM,
		const getfem::mesh_region& contact_region)
	{
		//extract matrix references
		MAT& Kmm = matl[0];
		MAT& Kss = matl[1];
		//MAT& Kll = matl[2] remains zero
		MAT& Kms = matl[3];
		MAT& Kml = matl[4];
		MAT& Ksl = matl[5];

		const std::string& mult_name = 
			mcb.get_pair_info(scb.get_var_name()).get_mult_name();
	    const std::string ls_name = "ls_on"+mcb.get_var_name()+"_from_"+scb.get_var_name();

		//extract mfs, and mims
		const mesh_fem& mf_U_line = mcb.get_mesh_fem();
		const mesh_fem& mf_lambda = mcb.get_model().mesh_fem_of_variable(mult_name);
		const mesh_fem& mf_interpolate = 
			mcb.get_pair_info(scb.get_var_name()).slave_scalar_fem();
		const mesh_fem& mf_U_interpolate = 
			mcb.get_pair_info(scb.get_var_name()).slave_vector_fem();
		const mesh_fem& mf_master_ls = mcb.get_model().mesh_fem_of_variable(ls_name);
		const mesh_im&  mim_line         =  
			mcb.get_pair_info(scb.get_var_name()).contact_mesh_im();

		//build temp vectors for interpolated fems
		plain_vector LS_small(mf_interpolate.nb_dof());
		gmm::copy(gmm::sub_vector(scb.ls_values(),
			mcb.get_pair_info(scb.get_var_name()).slave_scalar_dofs()),LS_small);

		//nonlinear term to compute normal vector and R matrix
		NormalTerm R_matrix(mcb,2);

		//nonlinear term that describes regularized integration or dummy (unity) multiplier
		dal::shared_ptr<getfem::nonlinear_elem_term> integration(0);
		if (mcb.integration==master_contact_body::REGULARIZED_LEVEL_SET){
			integration.reset(new HFunction(mf_master_ls,mcb.get_model().real_variable(ls_name),
				mcb.regularized_tollerance,mcb.small_weight_multiplier));
		} else {integration.reset(new Unity(mf_master_ls));}


		//temp matrices due to different DOF indeces of the slave
		sparse_matrix Kms_small(mf_U_line.nb_dof(),mf_U_interpolate.nb_dof());
		sparse_matrix Kss_small(mf_U_interpolate.nb_dof(),mf_U_interpolate.nb_dof());
		sparse_matrix Ksl_small(mf_U_interpolate.nb_dof(),mf_lambda.nb_dof());

		//assembly
		getfem::generic_assembly assem_boundary;

		assem_boundary.set(
			"F=data$1(#3);"
			"L=data$2(#1);"
			"Kmm1 = comp(Base(#1).Grad(#3).vBase(#2).NonLin$1(#2).vGrad(#2).NonLin$2(#5))(i,j,k,:,k,m,n,:,n,m,1).L(i).F(j);"
			"Kmm2 = comp(Base(#1).NonLin$1(#2).vGrad(#2).Grad(#3).vBase(#2).NonLin$2(#5))(i,m,n,:,n,m,j,k,:,k,1).L(i).F(j);"
			"Kmm3 = comp(Base(#1).Base(#3).NonLin$1(#2).vGrad(#2).NonLin$1(#2).vGrad(#2).NonLin$2(#5))(i,j,k,l,:,l,k,m,n,:,n,m,1).L(i).F(j);"
			"Kmm4 = (-1.0)*comp(Base(#1).Base(#3).NonLin$1(#2).vGrad(#2).vGrad(#2).NonLin$2(#5))(i,j,m,n,:,n,l,:,l,m,1).L(i).F(j);"
			"M$1(#2,#2)+= sym(Kmm1+Kmm2+Kmm3+Kmm4);"
			"Ksm1=(-1.0)*comp(Base(#1).Grad(#3).vBase(#4).NonLin$1(#2).vGrad(#2).NonLin$2(#5))(i,j,k,:,k,m,n,:,n,m,1).L(i).F(j);"
			"Ksm2=(-1.0)*comp(Base(#1).Grad(#3).vGrad(#4).vBase(#2).NonLin$2(#5))(i,j,m,:,m,n,:,n,1).L(i).F(j);"
			"M$2(#4,#2)+= Ksm1+Ksm2;"
			"Kml1=comp(Base(#3).NonLin$1(#2).vGrad(#2).Base(#1).NonLin$2(#5))(i,m,n,:,n,m,:,1).F(i);"
			"Kml2=comp(Grad(#3).vBase(#2).Base(#1).NonLin$2(#5))(i,j,:,j,:,1).F(i);"
			"M$3(#2,#1)+= Kml1+Kml2;"
			"Kss_part = comp(Base(#1).Grad(#3).vGrad(#4).vBase(#4).NonLin$2(#5))(i,j,m,:,m,n,:,n,1).L(i).F(j);"
			"M$4(#4,#4)+=sym(Kss_part{1,2}+Kss_part{2,1});"
			"M$5(#4,#1)+=(-1.0)*comp(Grad(#3).vBase(#4).Base(#1).NonLin$2(#5))(i,k,:,k,:,1).F(i);"
			); /* Here we don't compute matrices that contain Hessian of 
			   the level set function, as Getfem does not compute Hessian 
			   for interpolated_fem class that we use for level set function */
		assem_boundary.push_mi(mim_line);        //mim on the contact surface
		assem_boundary.push_mf(mf_lambda);       //mf 1	Lambda
		assem_boundary.push_mf(mf_U_line);       //mf 2	Umaster 
		assem_boundary.push_mf(mf_interpolate);  //mf 3 LSslave
		assem_boundary.push_mf(mf_U_interpolate);//mf 4 Uslave
		assem_boundary.push_mf(mf_master_ls);    //mf 5 ls_on_master
		assem_boundary.push_nonlinear_term(&R_matrix); //matrix of the normal products
		assem_boundary.push_nonlinear_term(integration.get()); //term to limit integration domain
		assem_boundary.push_data(LS_small);      //   data Level set on interpolated
		assem_boundary.push_data(LM);            //   data Lagrange mult values
		assem_boundary.push_mat(Kmm);                        //   result mat 1
		assem_boundary.push_mat(gmm::transposed(Kms_small)); //   ..     mat 2
		assem_boundary.push_mat(Kml);                        //   ..     mat 3
		assem_boundary.push_mat(Kss_small);                  //   ..     mat 4
		assem_boundary.push_mat(Ksl_small);                  //   ..     mat 5
		assem_boundary.assembly(contact_region);	

		//transfering from interpolated mesh_fem into full slave mesh_fem mat's
		const gmm::sub_interval& Um_dof = gmm::sub_interval(0,mf_U_line.nb_dof());
		const gmm::unsorted_sub_index& Us_dof = 
			mcb.get_pair_info(scb.get_var_name()).slave_vector_dofs();
		const gmm::sub_interval& LM_dof = gmm::sub_interval(0,mf_lambda.nb_dof());
		gmm::copy(Kms_small,gmm::sub_matrix(Kms,Um_dof,Us_dof));
		gmm::copy(Kss_small,gmm::sub_matrix(Kss,Us_dof,Us_dof));
		gmm::copy(Ksl_small,gmm::sub_matrix(Ksl,Us_dof,LM_dof));

	}

	template<typename VECT0,typename VECT1> 
	void asm_level_set_contact_rhs(
		std::vector<VECT0>& vecl, 
		const master_contact_body& mcb, 
		const slave_contact_body& scb, 
		const VECT1& LM,
		const getfem::mesh_region& contact_region)
	{
		//extract vector references
		VECT0& RHS_Um = vecl[0];
		VECT0& RHS_Us = vecl[1];
		VECT0& RHS_LM = vecl[2];
		// vecl[3,  4 and 5] remain zero


		const std::string& mult_name = 
			mcb.get_pair_info(scb.get_var_name()).get_mult_name();
	    const std::string ls_name = "ls_on"+mcb.get_var_name()+"_from_"+scb.get_var_name();

		//extract mfs, and mims
		const mesh_fem& mf_U_line = mcb.get_mesh_fem();
		const mesh_fem& mf_lambda = 
			mcb.get_model().mesh_fem_of_variable(mult_name);
		const mesh_fem& mf_interpolate = 
			mcb.get_pair_info(scb.get_var_name()).slave_scalar_fem();
		const mesh_fem& mf_U_interpolate = 
			mcb.get_pair_info(scb.get_var_name()).slave_vector_fem();
		const mesh_fem& mf_master_ls = mcb.get_model().mesh_fem_of_variable(ls_name);
		const mesh_im&  mim_line         =  
			mcb.get_pair_info(scb.get_var_name()).contact_mesh_im();

		//build temp vectors for interpolated fems
		plain_vector LS_small(mf_interpolate.nb_dof());
		gmm::copy(gmm::sub_vector(scb.ls_values(),
			mcb.get_pair_info(scb.get_var_name()).slave_scalar_dofs()),LS_small);

		//nonlinear term to compute normal vector and R matrix
		NormalTerm R_matrix(mcb,2);

		//nonlinear term that describes regularized integration or dummy (unity) multiplier
		dal::shared_ptr<getfem::nonlinear_elem_term> integration(0);
		if (mcb.integration==master_contact_body::REGULARIZED_LEVEL_SET){
			integration.reset(new HFunction(mf_master_ls,mcb.get_model().real_variable(ls_name),
				mcb.regularized_tollerance,mcb.small_weight_multiplier));
		} else {integration.reset(new Unity(mf_master_ls));}

		// temp RHS vector due to diff DOF indeces for mesh_fem object of the slave
		plain_vector RHS_Us_small(mf_U_interpolate.nb_dof());

		getfem::generic_assembly assem_boundary;
		assem_boundary.set(
			"F=data$1(#3);"
			"L=data$2(#1);"
			"RHS_L_Us_1=comp(Base(#1).Base(#3).NonLin$1(#2).vGrad(#2).NonLin$2(#5))(i,j,m,n,:,n,m,1).L(i).F(j);"
			"RHS_L_Us_2=comp(Base(#1).Grad(#3).vBase(#2).NonLin$2(#5))(i,j,k,:,k,1).L(i).F(j);"
			"V$1(#2)+=RHS_L_Us_1+RHS_L_Us_2;"
			"V$2(#4)+=(-1.0)*comp(Base(#1).Grad(#3).vBase(#4).NonLin$2(#5))(i,j,k,:,k,1).L(i).F(j);"
			"V$3(#1)+=comp(Base(#1).Base(#3).NonLin$2(#5))(:,i,1).F(i);"
			);
		assem_boundary.push_mi(mim_line);        //mim on the contact surface
		assem_boundary.push_mf(mf_lambda);       //mf 1	Lambda
		assem_boundary.push_mf(mf_U_line);       //mf 2	Umaster 
		assem_boundary.push_mf(mf_interpolate);  //mf 3 LSslave
		assem_boundary.push_mf(mf_U_interpolate);//mf 4 Uslave
		assem_boundary.push_mf(mf_master_ls);    //mf 5 ls_on_master
		assem_boundary.push_nonlinear_term(&R_matrix); //matrix of the normal products
		assem_boundary.push_nonlinear_term(integration.get()); //term to limit integration domain
		assem_boundary.push_data(LS_small);      //   data Level set on interpolated
		assem_boundary.push_data(LM);            //   data Lagrange mult values
		assem_boundary.push_vec(RHS_Um);         //   result vec 1 
		assem_boundary.push_vec(RHS_Us_small);   //   ..     vec 2
		assem_boundary.push_vec(RHS_LM);         //   ..     vec 3
		assem_boundary.assembly(contact_region);

		//transfering from interpolated mesh_fem into full slave mesh_fem RHS 
		const gmm::unsorted_sub_index& Us_dof = 
			mcb.get_pair_info(scb.get_var_name()).slave_vector_dofs();
		gmm::copy(RHS_Us_small, gmm::sub_vector(RHS_Us,Us_dof));

	}



	typedef void(*SOLVE_FUNCTION)(
		getfem::model &md, 
		gmm::iteration &iter,
		getfem::rmodel_plsolver_type solver,
		getfem::abstract_newton_line_search &ls, 
		bool with_pseudo_potential);

	/** Solves a model that has contact in it.
	Function checks wheather the contact area has converged
	@param sf - a pointer to a newton solver function, 
	can be, for instance, getfem::standard_solve 
	@param it_newton - iteration object for newton method
	@param it_staggered - iteration object for staggered calculation
	between conact detection and newton method (only max
	num. of iterations should be provided)
	@param lsolver - solver for a linear system
	@param ls      - reference to line search method
	@param with_pseudo_potential - yes if the bricks have pseude potential*/

	void solve_with_contact(
		SOLVE_FUNCTION sf, 
		getfem::model& md, 
		gmm::iteration& it_newton,
		gmm::iteration& it_staggered,
		const std::string& lsolver,
		getfem::abstract_newton_line_search &ls,
		bool with_pseudo_potential = false);

} //end of the namespace level_set_contact
