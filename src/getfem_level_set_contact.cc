/*===========================================================================
 
 Copyright (C) 2012-2015 Andriy Andreykiv
 
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

#include <getfem/getfem_level_set_contact.h> 
#include <getfem/getfem_interpolated_fem.h> 
#include <getfem/getfem_derivatives.h> 
#include <algorithm>
#include <getfem/getfem_level_set.h>
#include <getfem/getfem_mesh_level_set.h>
#include <getfem/getfem_mesh_im_level_set.h>
#include <math.h>

level_set_contact::contact_body::contact_body(model& _md, std::string _var_name):
        var_name(_var_name),
	is_deformed(false),
	own_mesh(const_cast<mesh&>(_md.mesh_fem_of_variable(_var_name).
	linked_mesh())),
	own_mesh_fem(_md.mesh_fem_of_variable(_var_name)),
	md(_md)
{}



level_set_contact::slave_contact_body::slave_contact_body(
	getfem::model& _md, const std::string& _var_name,
	getfem::mesh_im* _pmim_ls) : contact_body(_md,_var_name),
	ls_name("ls_on_"+_var_name),
	ls_mesh_fem(_md.mesh_fem_of_variable(_var_name)),
	pmim(_pmim_ls)

{
	ls_mesh_fem.set_qdim(size_type(1));
	model_real_plain_vector LS(ls_mesh_fem.nb_dof());
	boundary_level_set_field(get_mesh(),ls_mesh_fem,
		*pmim,LS);
	md.add_initialized_fem_data(ls_name,ls_mesh_fem,LS);
}


level_set_contact::slave_contact_body::slave_contact_body(
	getfem::model& _md, std::string _var_name, std::string _ls_name): 
contact_body(_md,_var_name),
	ls_name(_ls_name),
	pmim(0)
{ }

void level_set_contact::slave_contact_body::offset_level_set(scalar_type off)
{
	for(size_type i=0;i<ls_values().size();i++) 
		ls_values()[i]+=off;
}

level_set_contact::master_contact_body::master_contact_body(
	model& _md, 
	const std::string& _var_name,
	size_type _mult_order, 
	size_type _mult_mim_order) : 

        contact_body(_md,_var_name),
	mult_mim_order(_mult_mim_order),
	mult_int_method(""),
	mult_mf_order(_mult_order),
	integration(PER_ELEMENT),
	regularized_tollerance(0),
	small_weight_multiplier(0),
	max_contact_angle(45)
{
	//store existing elements in VOLUME_ELEMENTS region
	//before boundary elements are created
	VOLUME_ELEMENTS = getfem::mesh_region::free_region_id(get_mesh()); 
	get_mesh().region(VOLUME_ELEMENTS).add(get_mesh().convex_index());

	//create boundary elements (current mesh_fem should be automatically
	// extended with these elements)
	BOUNDARY_ELEMENTS = getfem::mesh_region::free_region_id(get_mesh()); 
	get_mesh().region(BOUNDARY_ELEMENTS).clear();

	masters.push_back(this);
}

level_set_contact::master_contact_body::master_contact_body(
	model& _md, 
	const std::string& _var_name,
	size_type _mult_order, 
	const std::string& _mult_mim_method,
	contact_integration _integration,
	scalar_type _regularized_tollerance,
	scalar_type _small_weight_multiplier,
	scalar_type _max_contact_angle): 

		contact_body(_md,_var_name),
		mult_mim_order(size_type(-1)),
		mult_int_method(_mult_mim_method),
		mult_mf_order(_mult_order),
		integration(_integration),
		regularized_tollerance(_regularized_tollerance),
		small_weight_multiplier(_small_weight_multiplier),
		max_contact_angle(_max_contact_angle)
{
	//store existing elements in VOLUME_ELEMENTS region
	//before boundary elements are created
	VOLUME_ELEMENTS = getfem::mesh_region::free_region_id(get_mesh()); 
	get_mesh().region(VOLUME_ELEMENTS).add(get_mesh().convex_index());

	//create boundary elements (current mesh_fem should be automatically
	// extended with these elements)
	BOUNDARY_ELEMENTS = getfem::mesh_region::free_region_id(get_mesh()); 
	get_mesh().region(BOUNDARY_ELEMENTS).clear();
	masters.push_back(this);
}


const level_set_contact::contact_pair_info& 
	level_set_contact::master_contact_body::get_pair_info(
	const std::string& slave_var_name) const
{   
	std::map<std::string, dal::shared_ptr<contact_pair_info> >
		::const_iterator it = contact_table.find(slave_var_name);
	if (it!=contact_table.end()) return *(it->second);
	GMM_ASSERT1(false,"did not find info on slave contact body, \
					  defined on variable "+slave_var_name);
}

level_set_contact::contact_pair_info& 
	level_set_contact::master_contact_body::get_pair_info(
	const std::string& slave_var_name)
{   
	std::map<std::string, dal::shared_ptr<contact_pair_info> >
		::iterator it = contact_table.find(slave_var_name);
	if (it!=contact_table.end()) return *(it->second);
	GMM_ASSERT1(false,"did not find info on slave contact body, \
					  defined on variable "+slave_var_name);
}


level_set_contact::face_type level_set_contact::master_contact_body::
	ext_face_of_elem(size_type i) const
{
	std::map<size_type,face_type>::const_iterator it = border_faces.find(i);
	if(it!=border_faces.end()) return it->second;
	GMM_ASSERT1(false,"did not find a face, corresponding to element "<<i);
}


void level_set_contact::master_contact_body::
	add_slave(slave_contact_body& scb, size_type assumed_contact_region)
{
	//check input
	GMM_ASSERT1(&md==&scb.get_model(),
		"Model objects of master and slave are not the same");
	if (assumed_contact_region!=size_type(-1)) 
		GMM_ASSERT1(get_mesh().region(assumed_contact_region).is_boundary(),
		"Assumed_contact_region must be on the boundary");

	//add surface elements where contact will be computed
	size_type assumed_contact_elems = getfem::mesh_region::free_region_id(get_mesh());
	getfem::mesh_region& contact_elems = get_mesh().region(assumed_contact_elems);
	getfem::mesh_region& boundary_elems = get_mesh().region(BOUNDARY_ELEMENTS);
	dal::shared_ptr<getfem::mr_visitor> i;
	getfem::mesh_region outer_faces;
	outer_faces.clear();
	getfem::outer_faces_of_mesh(get_mesh(), outer_faces);

	if (assumed_contact_region==size_type(-1)){ //all faces will be searched for contact
		i.reset(new getfem::mr_visitor(outer_faces));
	}
	else // only specified faces will be searched
	{  
		getfem::mesh_region& assumed_region = get_mesh().
			region(assumed_contact_region);
		i.reset(new getfem::mr_visitor(assumed_region));
	}

	for (; !i->finished(); ++(*i)){
		getfem::size_type new_elem = 
			get_mesh().add_convex(
			level_set_contact::face_trans_of_elem(
			get_mesh().trans_of_convex(i->cv())),
			get_mesh().ind_points_of_face_of_convex(
			i->cv(), i->f()).begin());


		border_faces[new_elem] = face_type(*i);
		contact_elems.add(new_elem);
		boundary_elems.add(new_elem);
	}

	GMM_ASSERT1(get_mesh().region(BOUNDARY_ELEMENTS).index().card()!=0,
		"No boundary elements added !!!");

	GMM_ASSERT1(get_mesh().region(assumed_contact_elems).index().card()!=0,
		"No contact elements added !!!");

	//creating Lagrange multiplier
	std::string mult_name = md.new_name("mult_on_"+get_var_name()+
		"_and_"+scb.get_var_name());
	const mesh_fem &mf_mult = 
		getfem::classical_mesh_fem(get_mesh(),
		bgeot::dim_type(mult_mf_order), bgeot::dim_type(1));
	md.add_multiplier(mult_name,mf_mult,get_var_name());

	//adding variable to store level set, projected from the slave
	const mesh_fem& mf_ls = 
	  getfem::classical_mesh_fem(get_mesh(),bgeot::dim_type(mult_mf_order+1));
	plain_vector LS(mf_ls.nb_dof());
	md.add_initialized_fem_data("ls_on"+get_var_name()+
			"_from_"+scb.get_var_name(),mf_ls,LS);

	//register contact pair
	contact_table[scb.get_var_name()] = 
		dal::shared_ptr<contact_pair_info>
		(new contact_pair_info(*this,scb,mult_name,assumed_contact_elems)); 

}

std::vector<level_set_contact::master_contact_body*> 
	level_set_contact::master_contact_body::masters;


bool level_set_contact::master_contact_body::any_contact_change()
{
	if (masters.size()==0) GMM_WARNING3("Running contact detection, while no \
										contact bodies are registered");
	bool contact_surfaces_changed = false;
	for(size_type i=0;i<masters.size();i++)
		if (masters[i]->master_contact_changed()) 
			contact_surfaces_changed=true;

	return contact_surfaces_changed;
}

void level_set_contact::master_contact_body::clear_all_contact_history()
{
	if (masters.size()==0) GMM_WARNING3("Clearing contact lists, while no \
										contact bodies are registered");
	for(size_type i=0;i<masters.size();i++)
		masters[i]->clear_contact_history();
}

void level_set_contact::master_contact_body::clear_contact_history()
{
	std::map<std::string, dal::shared_ptr<contact_pair_info> >::
		iterator it = contact_table.begin();
	for(;it!=contact_table.end();it++)
		it->second->clear_contact_history();
}

bool level_set_contact::master_contact_body::master_contact_changed()
{
	bool contact_surfaces_changed = false;
	std::map<std::string, dal::shared_ptr<contact_pair_info> >::
		iterator it = contact_table.begin();
	for(;it!=contact_table.end();it++)
		if (it->second->contact_changed()) 
			contact_surfaces_changed=true;

	return contact_surfaces_changed;
}

dal::shared_ptr<getfem::mesh_im> level_set_contact::master_contact_body::
	build_mesh_im_on_boundary(size_type region_id)
{

	dal::shared_ptr<getfem::mesh_im> pmim_contact;

		pmim_contact.reset(new mesh_im(get_mesh()));
		if (mult_mim_order!=size_type(-1)){
			pmim_contact->set_integration_method
                          (get_mesh().region(region_id).index(),
                           bgeot::dim_type(contact_mim_order()));
		}
		else
		{
                  pmim_contact->set_integration_method
                    (get_mesh().region(region_id).index(), contact_int_method());
		}

	return pmim_contact;
}


level_set_contact::contact_pair_update::contact_pair_update(
	master_contact_body& _mcb,
	slave_contact_body& _scb,
	update_depth ud):
mcb(_mcb), scb(_scb)

{

	GMM_ASSERT1(!mcb.is_mesh_deformed(),"Trying to deform \
										already deformed Master Contact Body");
	GMM_ASSERT1(!scb.is_mesh_deformed(),"Trying to deform \
										already deformed Slave  Contact Body");

	const model_real_plain_vector& 
		Umaster=mcb.get_model().real_variable(mcb.get_var_name());
	// size_type dof_check = Umaster.size();
	// size_type node_check = mcb.get_mesh().nb_points();
	def_master.reset(new getfem::temporary_mesh_deformator<>
		(mcb.get_mesh_fem(),Umaster));
	mcb.is_deformed=true;
	if (&mcb.get_mesh()!=&scb.get_mesh()){ 
		//  not deforming the slave if the master and the slave are the same
		const model_real_plain_vector& 
			Uslave=scb.get_model().real_variable(scb.get_var_name());
		def_slave.reset(new getfem::temporary_mesh_deformator<>
			(scb.get_mesh_fem(),Uslave));
		scb.is_deformed=true;
	}
	if (ud == FULL_UPDATE) mcb.update_for_slave(scb.get_var_name());
}

level_set_contact::contact_pair_update::~contact_pair_update(){
	mcb.is_deformed=false;
	scb.is_deformed=false;
}



level_set_contact::contact_pair_info::contact_pair_info(
	master_contact_body& underformed_mcb, 
	slave_contact_body& underformed_scb, 
	const std::string& _mult_name,
	size_type _GIVEN_CONTACT_REGION) :

master_cb(underformed_mcb),
	slave_cb(underformed_scb),
	mult_name(_mult_name),
	GIVEN_CONTACT_REGION(_GIVEN_CONTACT_REGION),

	ACTIVE_CONTACT_REGION(getfem::mesh_region::free_region_id(master_cb.get_mesh())),
	pmim_contact(0),
	ifem_srf(0),
	pinterpolated_fem(0),
	pinterpolated_fem_U(0),
	members_are_computed(false),
	init_cont_detect_done(false)

{
	//input check (if mult_name is incorrect, exception will be generated)
	// const mesh_fem& mf_mult=
	//	master_cb.get_model().mesh_fem_of_variable(mult_name);
	GMM_ASSERT1(master_cb.get_mesh().
		region(GIVEN_CONTACT_REGION).index().card()!=0,
		"provided contact region for contact_pair_info class is empty!!!");
}

void level_set_contact::contact_pair_info::clear_contact_history()
{
	old_contact_elm_list.clear();
	pre_old_ct_list.clear();
}

bool level_set_contact::contact_pair_info::contact_changed()
{
  //deform master and slave meshes
  contact_pair_update 
    temp_mesh_deformation(master_cb,slave_cb,DEFORM_MESHES_ONLY);
  
  // create mf on the boundary of the master (copy from the master)
  mesh_fem mf_scalar(master_cb.get_mesh());
  for(size_type i=0;i<mf_scalar.linked_mesh().nb_convex();i++)
    mf_scalar.set_finite_element(i,master_cb.get_mesh_fem().fem_of_element(i));
  
  mf_scalar.set_qdim(1);
  getfem::partial_mesh_fem mf_boundary(mf_scalar);
  mf_boundary.adapt(mf_scalar.dof_on_region(GIVEN_CONTACT_REGION));
  
  // interpolate level set from the slave to the master
  model_real_plain_vector LS_on_contour(mf_boundary.nb_dof());
  getfem::interpolation(slave_cb.get_ls_mesh_fem(), mf_boundary, 
                        slave_cb.ls_values(), LS_on_contour);
  model_real_plain_vector LS(mf_scalar.nb_dof());
  mf_boundary.extend_vector(LS_on_contour,LS);
  gmm::copy(LS,master_cb.get_model().set_real_variable
            ("ls_on"+master_cb.get_var_name()+"_from_"+slave_cb.get_var_name()));
  
  // interpolate the gradient of the level set onto the master surfaces
  // (this is to obtain the normal direction of the level set)
  mesh_fem mf_gradient_ls(slave_cb.get_mesh());
  mf_gradient_ls.set_classical_discontinuous_finite_element(bgeot::dim_type(master_cb.mult_mf_order));
  mesh_fem mf_gradient_ls_vect(mf_gradient_ls);
  mf_gradient_ls_vect.set_qdim(slave_cb.get_mesh().dim());
  plain_vector GradLS(mf_gradient_ls.nb_dof()*slave_cb.get_mesh().dim());
  getfem::compute_gradient(slave_cb.get_ls_mesh_fem(), mf_gradient_ls, slave_cb.ls_values(), GradLS);
  getfem::partial_mesh_fem mf_boundary_vect(master_cb.get_mesh_fem());
  mf_boundary_vect.adapt(master_cb.get_mesh_fem().dof_on_region(GIVEN_CONTACT_REGION));
  plain_vector GradLS_boundary(mf_boundary.nb_dof()*slave_cb.get_mesh().dim());
  getfem::interpolation(mf_gradient_ls_vect,mf_boundary_vect,GradLS,GradLS_boundary);

  size_type dim = slave_cb.get_mesh().dim();
  const scalar_type TINY = 1e-15;
  for(size_type i=0;i<mf_boundary.nb_dof();i++){ //normalizing the projected ls field
    bgeot::base_node ls_node(dim);
    for(size_type j=0;j<dim;j++) ls_node[j]=GradLS_boundary[dim*i+j];
    ls_node/= (gmm::vect_norm2(ls_node)+TINY);
    for(size_type j=0;j<dim;j++) GradLS_boundary[dim*i+j]=ls_node[j];
  }
  plain_vector normLS_master(master_cb.get_mesh_fem().nb_dof());
  mf_boundary_vect.extend_vector(GradLS_boundary,normLS_master);
  
  
  // extend Lagrange Multiplier onto the whole boundary of the master
  const mesh_fem& mf_mult = 
    master_cb.get_model().mesh_fem_of_variable(mult_name);
  const model_real_plain_vector& lambda = 
    master_cb.get_model().real_variable(mult_name);
  model_real_plain_vector lambda_full(mf_mult.nb_basic_dof());
  if (lambda.size()>0) mf_mult.extend_vector(lambda,lambda_full);
  // update contact region
  dal::bit_vector cc = master_cb.get_mesh().
    region(GIVEN_CONTACT_REGION).index();
  master_cb.get_mesh().region(ACTIVE_CONTACT_REGION).clear();
  master_cb.get_mesh().region(ACTIVE_CONTACT_REGION).add(cc);
  bgeot::size_type o;
  for (o << cc; o != bgeot::size_type(-1); o << cc) {
    getfem::mesh_fem::ind_dof_ct dof_ls = mf_scalar.ind_basic_dof_of_element(o);
    getfem::mesh_fem::ind_dof_ct dof_lm = mf_mult.ind_basic_dof_of_element(o);
    
    //measure the angle between ls countour and the master face
    face_type face = master_cb.ext_face_of_elem(o);
    bgeot::base_node unit_face_normal = 
      master_cb.get_mesh().normal_of_face_of_convex(face.cv,face.f);
    unit_face_normal/=gmm::vect_norm2(unit_face_normal);
    scalar_type cosine_alpha = 0;
    for (size_type j = 0; j < dof_ls.size(); j++){ 
      bgeot::base_node ls_grad_node(dim);
      for(size_type k=0; k<dim; k++) 
        ls_grad_node[k]=normLS_master[dim*dof_ls[j]+k];
      cosine_alpha += gmm::vect_sp(ls_grad_node,unit_face_normal);
    }
    cosine_alpha /= scalar_type(dof_ls.size()); 
    scalar_type alpha = acos(cosine_alpha)*360/(2*M_PI); // now this is average angle
    // between master surface and ls zero contour
    
    scalar_type LS_extreeme = LS[dof_ls[0]];
    if (master_cb.integration==master_contact_body::PER_ELEMENT)
      for (size_type j = 0; j < dof_ls.size(); j++) 
        LS_extreeme=std::min(LS[dof_ls[j]],LS_extreeme);
    else
      for (size_type j = 0; j < dof_ls.size(); j++) 
        LS_extreeme=std::max(LS[dof_ls[j]],LS_extreeme);
    
    scalar_type LM_sum = 0;
    for (size_type j = 0; j < dof_lm.size(); j++) 
      LM_sum+=lambda_full[dof_lm[j]];
    
    const scalar_type TINY_2 = 1e-9;
    
    if (LS_extreeme+LM_sum < TINY_2 || alpha > master_cb.max_contact_angle) 
      master_cb.get_mesh().region(ACTIVE_CONTACT_REGION).sup(o);
  }
  
  // check whether contact areas have changed
  bool contact_surface_changed;
  const dal::bit_vector& current_contact_elm_list = 
    master_cb.get_mesh().region(ACTIVE_CONTACT_REGION).index();
  GMM_TRACE2("Current contact elements: "<< current_contact_elm_list);
  GMM_TRACE2("Old contact elements:     "<< old_contact_elm_list);
  GMM_TRACE2("Pre-old contact elements: "<< pre_old_ct_list);
  
  if (current_contact_elm_list == old_contact_elm_list && 
      current_contact_elm_list.card() == old_contact_elm_list.card()) {
    contact_surface_changed = false;
    GMM_TRACE2("   the contact area has not changed");
  } else {
    if (current_contact_elm_list == pre_old_ct_list &&
        current_contact_elm_list.card() == pre_old_ct_list.card()) {
      contact_surface_changed = false;
      GMM_TRACE2("   the contact area has changed, but cycling, \
						   so exiting active set search");
    } else {
      contact_surface_changed = true;
      GMM_TRACE2("   the contact area has changed");
      pre_old_ct_list = old_contact_elm_list;
      old_contact_elm_list = current_contact_elm_list;
    }
  }
  
  init_cont_detect_done = true;
  force_update();
  
  
  //building integration method
  pmim_contact = master_cb.build_mesh_im_on_boundary(ACTIVE_CONTACT_REGION);
  n_integrated_elems = pmim_contact->convex_index().card();
  GMM_ASSERT1(n_integrated_elems==current_contact_elm_list.card(),
              "Failure in integration method: The number of integrated elements "
              "does not correspond to the number of contact elements");
  
  return contact_surface_changed;  
  
}



void level_set_contact::contact_pair_info::update() const
{
	GMM_ASSERT1(master_cb.is_mesh_deformed(),"Master mesh is not deformed, \
											 cannot calucalte contact info");

	GMM_ASSERT1(slave_cb.is_mesh_deformed(),"Slave mesh is not deformed, \
											cannot calucalte contact info");

	GMM_ASSERT1(master_cb.get_mesh().region(ACTIVE_CONTACT_REGION).index().card()>0,
		"Internal error: Contact area is empty");

	//pinterpolated_fem for level set
	pinterpolated_fem.reset(new mesh_fem(master_cb.get_mesh()));
	if (ifem_srf.get()!=0) getfem::del_interpolated_fem(ifem_srf);
	ifem_srf=getfem::new_interpolated_fem(
		slave_cb.get_ls_mesh_fem(),*pmim_contact);
	pinterpolated_fem->set_finite_element(
		master_cb.get_mesh().region(ACTIVE_CONTACT_REGION).index(),ifem_srf);
	pinterpolated_fem->set_qdim(1);


	//pinterpolated_fem_U
	pinterpolated_fem_U.reset(new mesh_fem(master_cb.get_mesh()));
	pinterpolated_fem_U->set_finite_element(master_cb.get_mesh().
		region(ACTIVE_CONTACT_REGION).index(),ifem_srf);
	pinterpolated_fem_U->set_qdim(master_cb.get_mesh().dim());

	//slave_ls_dofs
	std::vector<size_type> index(pinterpolated_fem->nb_dof());
	dal::bit_vector cc = 
		master_cb.get_mesh().region(ACTIVE_CONTACT_REGION).index();
	for (dal::bv_visitor icv(cc); !icv.finished(); ++icv){
		for (size_type j = 0; j < pinterpolated_fem->nb_basic_dof_of_element(icv); 
			++j) 
		{index[pinterpolated_fem->ind_basic_dof_of_element(icv)[j]]
			= ifem_srf->index_of_global_dof(icv, j);}
	}

	slave_ls_dofs.reset(new gmm::unsorted_sub_index(index));

	//slave_U_dofs
	std::vector<size_type> indexU(pinterpolated_fem_U->nb_dof());
	size_type dim = pinterpolated_fem_U->get_qdim();
	for(size_type d=0;d<dim;d++)
		for(size_type i=0;i<pinterpolated_fem->nb_dof();i++)
			indexU[dim*i+d] = dim*index[i]+d;
	slave_U_dofs.reset(new gmm::unsorted_sub_index(indexU));

	members_are_computed=true;
}



getfem::size_type level_set_contact::add_level_set_normal_contact_brick(
	model& md, 
	master_contact_body& mcb, 
	slave_contact_body& scb,
	size_type rg)
{
	//level set contact class
	getfem::pbrick pbr = 
		new level_set_contact_brick(md,mcb,scb,rg);

	//term description
	const std::string& name_Um = mcb.get_var_name();
	const std::string& name_Us = scb.get_var_name();
	const std::string& name_LM = mcb.get_pair_info(name_Us).get_mult_name();
	model::termlist terms;
	terms.push_back(model::term_description(name_Um,name_Um,false));
	terms.push_back(model::term_description(name_Us,name_Us,false));
	terms.push_back(model::term_description(name_LM,name_LM,false));
	terms.push_back(model::term_description(name_Um,name_Us,true));
	terms.push_back(model::term_description(name_Um,name_LM,true));
	terms.push_back(model::term_description(name_Us,name_LM,true));

	//variables
	model::varnamelist variables;
	variables.push_back(name_Um);
	variables.push_back(name_Us);
	variables.push_back(name_LM);

	//empty data and integration method lists
	//(we don't have any properties or initial data,
	//while integration methods are created per iteration)
	model::varnamelist datalist;
	model::mimlist mimlist;

	//register the brick with the model and return its number
	return md.add_brick(pbr,variables,datalist,terms,mimlist,rg);
}


level_set_contact::level_set_contact_brick::
	level_set_contact_brick(
	model& _md,
	master_contact_body& _mcb, 
	slave_contact_body& _scb, 
	size_type rg) : 
md(_md),mcb(_mcb),scb(_scb), given_contact_id(rg)
{
	GMM_ASSERT1(&md == &mcb.get_model(),
		"Master body is defined on a different model then the input");

	//register master/slave pair
	mcb.add_slave(scb,given_contact_id);

	//Reduce computation to own MPI region
	contact_region_id = mcb.get_pair_info(scb.get_var_name()).contact_region();
	getfem::mesh_region& contact_region_ = mcb.get_mesh().region(contact_region_id);
	mcb.get_mesh().intersect_with_mpi_region(contact_region_);
	contact_region_id = contact_region_.id(); //probably not needed, but still


	set_flags("Level set contact brick", 
		false /* is linear*/,
		true  /* is symmetric */, 
		false /* is coercive */,
		true  /* is real */, 
		false /* is complex */);

}

void level_set_contact::level_set_contact_brick::asm_real_tangent_terms(
	const model &mdd, size_type /* ib */,
	const model::varnamelist &vl,
	const model::varnamelist &/* dl */,
	const model::mimlist &/* mims */,
	model::real_matlist &matl,
	model::real_veclist &vecl,
	model::real_veclist &,
	size_type region,
	build_version version) const {
  //input check
  GMM_ASSERT1(vl.size() == 3,
              "Level set contact  brick needs three variables");
  GMM_ASSERT1(matl.size() == 6,
              "Level set contact  brick needs six matrices");
  GMM_ASSERT1(vecl.size() == 6,
              "Level set contact  brick assembles size RHSs");
  GMM_ASSERT1(region==given_contact_id,
		"Assumed contact region has changed!!! \
		This implementation does not handle this \
		for efficiency reasons!!");
  
  if (version & model::BUILD_MATRIX ) 
    for(size_type i=0;i<matl.size();i++) gmm::clear(matl[i]);
  if (version & model::BUILD_RHS )
    for(size_type i=0;i<vecl.size();i++) gmm::clear(vecl[i]);
  
  const getfem::mesh_region& active_contact_region = 
    mcb.get_mesh().region(contact_region_id);
  if (active_contact_region.index().card()==0) return; //no contact -> no contact assembly
  
  //deform the meshes, update contact info
  contact_pair_update cp_update(mcb,scb,FULL_UPDATE);
  
  //extract DOF vectors
  const plain_vector &LM = mdd.real_variable(vl[2]);

  //Assemble Tangent Matrix
  if (version & model::BUILD_MATRIX ) {
    GMM_TRACE2("Level set contact brick stiffness matrix assembly on "
               << mcb.get_pair_info(scb.get_var_name()).num_of_integr_elems()
               << " elements");
    asm_level_set_contact_tangent_matrix(matl,mcb,scb,LM,active_contact_region);
  }

  //Assemble RHS
  if (version & model::BUILD_RHS ) {
    GMM_TRACE2("Level set contact brick RHS assembly on "
               << mcb.get_pair_info(scb.get_var_name()).num_of_integr_elems()
               << " elements");
    asm_level_set_contact_rhs(vecl,mcb,scb,LM,active_contact_region);
    for(size_type i=0;i<vecl.size();i++)	
      gmm::scale(vecl[i], scalar_type(-1));
  }
}


void level_set_contact::NormalTerm::compute(
	getfem::fem_interpolation_context& ctx, 
	bgeot::base_tensor &t) 
{
	size_type cv = ctx.convex_num();
	size_type cv_volume = mcb.ext_face_of_elem(cv).cv;
        bgeot::short_type f_volume  = mcb.ext_face_of_elem(cv).f;
	bgeot::base_node un = mcb.get_mesh().normal_of_face_of_convex
	  (cv_volume, bgeot::short_type(f_volume), ctx.xref());
	un /= gmm::vect_norm2(un);

	if (version == 1) {
		for (size_type i = 0; i < dim; i++) t[i] = un[i];
	} else {
		for (size_type i = 0; i < dim; i++)
			for (size_type j = 0; j < dim; j++)
				if (i == j) t(i, j) = 1.0 - un[i] * un[j];
				else t(i, j) =-un[i] * un[j];
	}
}


level_set_contact::HFunction::HFunction(
	const mesh_fem &lsmf_,
	const plain_vector &LS_U_,
	scalar_type epsilon, 
	scalar_type small_h_):

lsmf(lsmf_),
	LS_U(LS_U_),
	m_Epsilon(epsilon),
	small_h(small_h_),
	sizes_(1)
{sizes_[0]=1;}

const bgeot::multi_index& level_set_contact::HFunction::
	sizes(size_type) const { return sizes_;}

void level_set_contact::HFunction::
prepare(getfem::fem_interpolation_context& /*ctx*/, size_type /*nl_part*/) {}

void  level_set_contact::HFunction::compute(getfem::fem_interpolation_context& ctx, 
	bgeot::base_tensor &t)
{
	size_type cv = ctx.convex_num();
        plain_vector U;
        slice_vector_on_basic_dof_of_element(lsmf, LS_U, cv, U);
	// plain_vector U(lsmf.nb_basic_dof_of_element(cv));
	// gmm::copy(gmm::sub_vector(LS_U,gmm::sub_index(lsmf.ind_basic_dof_of_element(cv))),U);
	plain_vector ls_interpolated(1);
	ctx.pf()->interpolation(ctx,U,ls_interpolated,1);
	t[0] = hRegularized(ls_interpolated[0],m_Epsilon,small_h);
} 

bgeot::scalar_type level_set_contact::HFunction::
	hRegularized(scalar_type f, scalar_type epsilon, scalar_type small_h_)
{
	if (f>epsilon) return 1.0;
	if (f<(-epsilon)) return small_h_;
	return 0.5+0.125*(9.0*f/(epsilon)-5.0*pow(f/(epsilon),scalar_type(3)));
}


level_set_contact::Unity::Unity(const mesh_fem &mf_):mf(mf_),sizes_(1)
{sizes_[0]=1;}
const bgeot::multi_index& level_set_contact::Unity::sizes(size_type) const {return sizes_;}
void level_set_contact::Unity::
prepare(getfem::fem_interpolation_context& /*ctx*/, size_type /*nl_part*/) {}
void level_set_contact::Unity::
compute(getfem::fem_interpolation_context& /*ctx*/, bgeot::base_tensor &t){t[0]=1.0;}


void level_set_contact::
	solve_with_contact(
	SOLVE_FUNCTION sf, 
	getfem::model& md, 
	gmm::iteration& it_newton,
	gmm::iteration& it_staggered,
	const std::string& lsolver_name,
	getfem::abstract_newton_line_search &ls)
{
	bool active_set_converged = false;
	it_staggered.set_iteration(0);
	master_contact_body::clear_all_contact_history();

	do {
		active_set_converged = !master_contact_body::any_contact_change();
		it_newton.set_iteration(0);
		getfem::rmodel_plsolver_type plsolver=getfem::select_linear_solver<sparse_matrix,plain_vector>(md,lsolver_name);
		(*sf)(md,it_newton,plsolver,ls);
		GMM_TRACE2("Newton converged?  - "<<it_newton.converged());
		GMM_TRACE2("active set converged?  - "<<active_set_converged);
		GMM_ASSERT1(it_newton.converged(),"Newton method did not converge");
		it_staggered++;

	} while(!active_set_converged && 
		!it_staggered.finished(it_staggered.get_resmax()+1.0));

	if (active_set_converged && it_newton.converged()) it_staggered.enforce_converged();
}

bgeot::pgeometric_trans 
	level_set_contact::face_trans_of_elem(
	bgeot::pgeometric_trans pelem_trans)
{
	std::string name = bgeot::name_of_geometric_trans(pelem_trans);
	std::stringstream fname;
	fname<<name.substr(0,5);
	GMM_ASSERT1((fname.str()=="GT_QK" || fname.str()=="GT_PK"),
		"Cannot handle other transformations but QK or PK,\
		Sorry, to be implemented" );
	std::stringstream str1(name.substr(6,1));
	size_type dim;
	str1>>dim;

	std::istringstream str2(name.substr(8,1));
	size_type order;
	str2>>order;

	fname<<"("<<dim-1<<","<<order<<")";
	return bgeot::geometric_trans_descriptor(fname.str());
}
