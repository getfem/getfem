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

#include <vector>
#include <vector>
#include <gmm/gmm.h>
#include <getfem/getfem_interpolated_fem.h>
#include <getfem/bgeot_mesh.h> 
#include <getfem/getfem_import.h>
#include <getfem/getfem_assembling.h>
#include <getfem/getfem_interpolation.h>
#include <getfem/getfem_export.h>
#include <getfem/getfem_nonlinear_elasticity.h>
#include <getfem/getfem_level_set_contact.h>
#include <getfem/getfem_model_solvers.h>
#include <gmm/gmm_except.h>
#include "contact_problem.h"



int main(int argc, char *argv[])
{

	GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
	FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

        contact_problem p(argc,argv);

        //model
        getfem::model model;


        //mfs
        getfem::mesh_fem mf_master(p.mesh_master);
        //must use set_classical_finite_element
        // for the master, so that it has "auto_add" feature
        // contact algorithm will fail otherwise
        mf_master.set_classical_finite_element(bgeot::dim_type(p.app_order_master));
        mf_master.set_qdim(p.mesh_master.dim());
        getfem::mesh_fem mf_slave(p.mesh_slave);
        //set_classical_finite_element is not mandatory for slaves
        mf_slave.set_classical_finite_element(bgeot::dim_type(p.app_order_slave));
        mf_slave.set_qdim(p.mesh_slave.dim());
        
        //mims
        getfem::mesh_im mim_master(p.mesh_master);
        mim_master.set_integration_method(p.mesh_master.convex_index(),bgeot::dim_type(p.int_order_master));
        getfem::mesh_im mim_slave(p.mesh_slave);
        mim_slave.set_integration_method(p.mesh_slave.convex_index(),bgeot::dim_type(p.int_order_slave));
        
        //variables
        model.add_fem_variable("U_master",mf_master);
        model.add_fem_variable("U_slave" ,mf_slave);
        
        //materials
        getfem::phyperelastic_law mat_law = std::make_shared<getfem::SaintVenant_Kirchhoff_hyperelastic_law>();
        bgeot::base_vector mat_param_master(2); 
        mat_param_master[0] = p.lambda_master; mat_param_master[1] = p.mu_master;
        bgeot::base_vector mat_param_slave(2); 
        mat_param_slave[0] = p.lambda_slave; mat_param_slave[1] = p.mu_slave;
        
        
        //nonlinear elasticity bricks (can also use updated Lagrangian)
        model.add_initialized_fixed_size_data("params_master", mat_param_master);
        model.add_initialized_fixed_size_data("params_slave" , mat_param_slave );
        getfem::add_nonlinear_elasticity_brick(model, mim_master, "U_master", mat_law, "params_master");
        getfem::add_nonlinear_elasticity_brick(model, mim_slave,   "U_slave", mat_law, "params_slave");
        
        //Fixed Dirichlet on slaves bottom
        getfem::add_Dirichlet_condition_with_multipliers(model,mim_slave,
                                                         "U_slave",bgeot::dim_type(p.app_order_slave),SOUTH);
        
        //normal Dirichet's on all vertical wals of the master and the slave
        getfem::add_normal_Dirichlet_condition_with_multipliers
          (model,mim_master,"U_master",
           bgeot::dim_type(p.app_order_master),EAST);
        getfem::add_normal_Dirichlet_condition_with_multipliers
          (model,mim_master,"U_master",
           bgeot::dim_type(p.app_order_master),WEST);
        getfem::add_normal_Dirichlet_condition_with_multipliers
          (model,mim_slave,"U_slave",
           bgeot::dim_type(p.app_order_slave),EAST);
        getfem::add_normal_Dirichlet_condition_with_multipliers
          (model,mim_slave,"U_slave",
           bgeot::dim_type(p.app_order_slave),WEST);
        if (p.model_dim==3){
          getfem::add_normal_Dirichlet_condition_with_multipliers(model,mim_master,"U_master",
                                                                  bgeot::dim_type(p.app_order_master),FRONT);
          getfem::add_normal_Dirichlet_condition_with_multipliers(model,mim_master,"U_master",
                                                                  bgeot::dim_type(p.app_order_master),BACK);
          getfem::add_normal_Dirichlet_condition_with_multipliers
            (model,mim_slave,"U_slave",
             bgeot::dim_type(p.app_order_slave),FRONT);
          getfem::add_normal_Dirichlet_condition_with_multipliers
            (model,mim_slave,"U_slave",
             bgeot::dim_type(p.app_order_slave),BACK);
          
        }
        
        //Normal dirichlet condition assigned on the top side of the master
        // it gives the actual movement of the model
        bgeot::base_vector moving_dirichlet(1); moving_dirichlet[0]=0;
        model.add_initialized_fixed_size_data("moving_dirichlet",moving_dirichlet);
        getfem::add_normal_Dirichlet_condition_with_multipliers
          (model, mim_master,"U_master",bgeot::dim_type(p.app_order_master),NORTH,"moving_dirichlet");
        
        
        
        //CONTACT DEFINITION
        //   approximation order for Lagrange Mult
        size_type LM_approximation_order = p.app_order_master-1; //must be lower than app_order
        //   master contact body
        level_set_contact::master_contact_body mcb(
                                                   model,
                                                   "U_master",
                                                   LM_approximation_order,
                                                   p.LM_INT_TYPE,
                                                   level_set_contact::master_contact_body::/*::PER_ELEMENT*/REGULARIZED_LEVEL_SET, // integration approach
                                                   1e-7,                       // regularazied transition width
                                                   0.01,                       // negative ls gauss point weight
                                                   30                          // max allowed contact angle
                                                   );
        //   the slave
        level_set_contact::slave_contact_body scb(model,"U_slave",&mim_slave);
        //   can be used to move LS zero inside the mesh
        scb.offset_level_set(p.ls_offset); 
        //   contact brick
        level_set_contact::add_level_set_normal_contact_brick(model,mcb,scb);
        
        
        //solving
        gmm::iteration iter_newton(p.tol_newton,1,99999);
        gmm::iteration iter_contact(1e-4,1,100000);
        scalar_type applied_disp = p.applied_disp;
        size_type nstep=p.nstep;
        scalar_type disp_incr = applied_disp/scalar_type(nstep);
        
        for(size_type step=0; step<=nstep; step++)
          {
            std::stringstream s; s<<step;
            
            std::cout<<"step "<<s.str()<<std::endl;
            std::cout<<"Current displacement: "<<moving_dirichlet[0]<<std::endl;
            getfem::basic_newton_line_search  line_search;
            
            
            //actual step solving
            level_set_contact::solve_with_contact(getfem::standard_solve,model,
                                                  iter_newton,iter_contact,"superlu",line_search);
            
            GMM_ASSERT1(iter_contact.converged(),"ERROR: contact algorithm did not converge");
            std::cout << "update" << std::endl;
            
            //updating displacements
            moving_dirichlet[0]+=disp_incr;
            gmm::copy(moving_dirichlet,model.set_real_variable("moving_dirichlet"));
            
            std::cout << "post" << std::endl;
            //post-processing
            getfem::vtk_export exp_m("master_"+s.str()+".vtk",true);
            exp_m.exporting(mcb.get_mesh());
            exp_m.write_point_data(mf_master,model.real_variable("U_master"),"displacement");
            getfem::vtk_export exp_s("slave_"+s.str()+".vtk",true);
            exp_s.exporting(scb.get_mesh());
            exp_s.write_point_data(scb.get_mesh_fem(),model.real_variable("U_slave"),"displacement");
            exp_s.write_point_data(scb.get_ls_mesh_fem(),scb.ls_values(),"level set");
             std::cout << "end iter " << std::endl;
          }
        
	return 0;
        
}
