function [varargout]=subsref(obj, index)
% gfModel/subsref
  nout = max(nargout, 1); cnt=1;
  FGET = @gf_model_get;
  FSET = @gf_model_set;
  switch index(1).type
    case '{}'
      error('Cell array indexing not supported by gfModel objects')
    case '()'
      error('array indexing not supported by gfModel objects')
    case '.'
      switch index(1).subs
        case 'id'
          [varargout{1:nout}] = obj.id;
        case 'clear'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'clear', index(2).subs{:});
          else
            FSET(obj, 'clear', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_fem_variable'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_fem_variable', index(2).subs{:});
          else
            FSET(obj, 'add_fem_variable', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_filtered_fem_variable'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_filtered_fem_variable', index(2).subs{:});
          else
            FSET(obj, 'add_filtered_fem_variable', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_variable'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_variable', index(2).subs{:});
          else
            FSET(obj, 'add_variable', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'delete_variable'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'delete_variable', index(2).subs{:});
          else
            FSET(obj, 'delete_variable', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'resize_variable'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'resize_variable', index(2).subs{:});
          else
            FSET(obj, 'resize_variable', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_multiplier'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_multiplier', index(2).subs{:});
          else
            FSET(obj, 'add_multiplier', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_im_data'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_im_data', index(2).subs{:});
          else
            FSET(obj, 'add_im_data', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_fem_data'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_fem_data', index(2).subs{:});
          else
            FSET(obj, 'add_fem_data', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_initialized_fem_data'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_initialized_fem_data', index(2).subs{:});
          else
            FSET(obj, 'add_initialized_fem_data', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_data'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_data', index(2).subs{:});
          else
            FSET(obj, 'add_data', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_macro'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_macro', index(2).subs{:});
          else
            FSET(obj, 'add_macro', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'del_macro'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'del_macro', index(2).subs{:});
          else
            FSET(obj, 'del_macro', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_initialized_data'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_initialized_data', index(2).subs{:});
          else
            FSET(obj, 'add_initialized_data', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_variable'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'variable', index(2).subs{:});
          else
            FSET(obj, 'variable', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'to_variables'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'to_variables', index(2).subs{:});
          else
            FSET(obj, 'to_variables', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'delete_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'delete_brick', index(2).subs{:});
          else
            FSET(obj, 'delete_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'define_variable_group'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'define_variable_group', index(2).subs{:});
          else
            FSET(obj, 'define_variable_group', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_elementary_rotated_RT0_projection'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_elementary_rotated_RT0_projection', index(2).subs{:});
          else
            FSET(obj, 'add_elementary_rotated_RT0_projection', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_interpolate_transformation_from_expression'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_interpolate_transformation_from_expression', index(2).subs{:});
          else
            FSET(obj, 'add_interpolate_transformation_from_expression', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_element_extrapolation_transformation'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_element_extrapolation_transformation', index(2).subs{:});
          else
            FSET(obj, 'add_element_extrapolation_transformation', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_element_extrapolation_correspondance'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'set_element_extrapolation_correspondance', index(2).subs{:});
          else
            FSET(obj, 'set_element_extrapolation_correspondance', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_raytracing_transformation'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_raytracing_transformation', index(2).subs{:});
          else
            FSET(obj, 'add_raytracing_transformation', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_master_contact_boundary_to_raytracing_transformation'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_master_contact_boundary_to_raytracing_transformation', index(2).subs{:});
          else
            FSET(obj, 'add_master_contact_boundary_to_raytracing_transformation', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_slave_contact_boundary_to_raytracing_transformation'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_slave_contact_boundary_to_raytracing_transformation', index(2).subs{:});
          else
            FSET(obj, 'add_slave_contact_boundary_to_raytracing_transformation', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_rigid_obstacle_to_raytracing_transformation'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_rigid_obstacle_to_raytracing_transformation', index(2).subs{:});
          else
            FSET(obj, 'add_rigid_obstacle_to_raytracing_transformation', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_projection_transformation'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_projection_transformation', index(2).subs{:});
          else
            FSET(obj, 'add_projection_transformation', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_master_contact_boundary_to_projection_transformation'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_master_contact_boundary_to_projection_transformation', index(2).subs{:});
          else
            FSET(obj, 'add_master_contact_boundary_to_projection_transformation', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_slave_contact_boundary_to_projection_transformation'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_slave_contact_boundary_to_projection_transformation', index(2).subs{:});
          else
            FSET(obj, 'add_slave_contact_boundary_to_projection_transformation', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_rigid_obstacle_to_projection_transformation'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_rigid_obstacle_to_projection_transformation', index(2).subs{:});
          else
            FSET(obj, 'add_rigid_obstacle_to_projection_transformation', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_linear_term'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_linear_term', index(2).subs{:});
          else
            FSET(obj, 'add_linear_term', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_linear_generic_assembly_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_linear_generic_assembly_brick', index(2).subs{:});
          else
            FSET(obj, 'add_linear_generic_assembly_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_nonlinear_term'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_nonlinear_term', index(2).subs{:});
          else
            FSET(obj, 'add_nonlinear_term', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_nonlinear_generic_assembly_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_nonlinear_generic_assembly_brick', index(2).subs{:});
          else
            FSET(obj, 'add_nonlinear_generic_assembly_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_source_term'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_source_term', index(2).subs{:});
          else
            FSET(obj, 'add_source_term', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_source_term_generic_assembly_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_source_term_generic_assembly_brick', index(2).subs{:});
          else
            FSET(obj, 'add_source_term_generic_assembly_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_assembly_assignment'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_assembly_assignment', index(2).subs{:});
          else
            FSET(obj, 'add_assembly_assignment', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'clear_assembly_assignment'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'clear_assembly_assignment', index(2).subs{:});
          else
            FSET(obj, 'clear_assembly_assignment', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Laplacian_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Laplacian_brick', index(2).subs{:});
          else
            FSET(obj, 'add_Laplacian_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_generic_elliptic_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_generic_elliptic_brick', index(2).subs{:});
          else
            FSET(obj, 'add_generic_elliptic_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_source_term_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_source_term_brick', index(2).subs{:});
          else
            FSET(obj, 'add_source_term_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_normal_source_term_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_normal_source_term_brick', index(2).subs{:});
          else
            FSET(obj, 'add_normal_source_term_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Dirichlet_condition_with_simplification'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Dirichlet_condition_with_simplification', index(2).subs{:});
          else
            FSET(obj, 'add_Dirichlet_condition_with_simplification', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Dirichlet_condition_with_multipliers'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Dirichlet_condition_with_multipliers', index(2).subs{:});
          else
            FSET(obj, 'add_Dirichlet_condition_with_multipliers', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Dirichlet_condition_with_Nitsche_method'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Dirichlet_condition_with_Nitsche_method', index(2).subs{:});
          else
            FSET(obj, 'add_Dirichlet_condition_with_Nitsche_method', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Dirichlet_condition_with_penalization'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Dirichlet_condition_with_penalization', index(2).subs{:});
          else
            FSET(obj, 'add_Dirichlet_condition_with_penalization', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_normal_Dirichlet_condition_with_multipliers'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_normal_Dirichlet_condition_with_multipliers', index(2).subs{:});
          else
            FSET(obj, 'add_normal_Dirichlet_condition_with_multipliers', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_normal_Dirichlet_condition_with_penalization'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_normal_Dirichlet_condition_with_penalization', index(2).subs{:});
          else
            FSET(obj, 'add_normal_Dirichlet_condition_with_penalization', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_normal_Dirichlet_condition_with_Nitsche_method'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_normal_Dirichlet_condition_with_Nitsche_method', index(2).subs{:});
          else
            FSET(obj, 'add_normal_Dirichlet_condition_with_Nitsche_method', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_generalized_Dirichlet_condition_with_multipliers'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_generalized_Dirichlet_condition_with_multipliers', index(2).subs{:});
          else
            FSET(obj, 'add_generalized_Dirichlet_condition_with_multipliers', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_generalized_Dirichlet_condition_with_penalization'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_generalized_Dirichlet_condition_with_penalization', index(2).subs{:});
          else
            FSET(obj, 'add_generalized_Dirichlet_condition_with_penalization', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_generalized_Dirichlet_condition_with_Nitsche_method'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_generalized_Dirichlet_condition_with_Nitsche_method', index(2).subs{:});
          else
            FSET(obj, 'add_generalized_Dirichlet_condition_with_Nitsche_method', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_pointwise_constraints_with_multipliers'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_pointwise_constraints_with_multipliers', index(2).subs{:});
          else
            FSET(obj, 'add_pointwise_constraints_with_multipliers', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_pointwise_constraints_with_given_multipliers'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_pointwise_constraints_with_given_multipliers', index(2).subs{:});
          else
            FSET(obj, 'add_pointwise_constraints_with_given_multipliers', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_pointwise_constraints_with_penalization'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_pointwise_constraints_with_penalization', index(2).subs{:});
          else
            FSET(obj, 'add_pointwise_constraints_with_penalization', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'change_penalization_coeff'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'change_penalization_coeff', index(2).subs{:});
          else
            FSET(obj, 'change_penalization_coeff', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Helmholtz_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Helmholtz_brick', index(2).subs{:});
          else
            FSET(obj, 'add_Helmholtz_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Fourier_Robin_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Fourier_Robin_brick', index(2).subs{:});
          else
            FSET(obj, 'add_Fourier_Robin_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_constraint_with_multipliers'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_constraint_with_multipliers', index(2).subs{:});
          else
            FSET(obj, 'add_constraint_with_multipliers', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_constraint_with_penalization'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_constraint_with_penalization', index(2).subs{:});
          else
            FSET(obj, 'add_constraint_with_penalization', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_explicit_matrix'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_explicit_matrix', index(2).subs{:});
          else
            FSET(obj, 'add_explicit_matrix', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_explicit_rhs'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_explicit_rhs', index(2).subs{:});
          else
            FSET(obj, 'add_explicit_rhs', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_private_matrix'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'set_private_matrix', index(2).subs{:});
          else
            FSET(obj, 'set_private_matrix', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_private_rhs'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'set_private_rhs', index(2).subs{:});
          else
            FSET(obj, 'set_private_rhs', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_isotropic_linearized_elasticity_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_isotropic_linearized_elasticity_brick', index(2).subs{:});
          else
            FSET(obj, 'add_isotropic_linearized_elasticity_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_isotropic_linearized_elasticity_brick_pstrain'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_isotropic_linearized_elasticity_brick_pstrain', index(2).subs{:});
          else
            FSET(obj, 'add_isotropic_linearized_elasticity_brick_pstrain', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_isotropic_linearized_elasticity_brick_pstress'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_isotropic_linearized_elasticity_brick_pstress', index(2).subs{:});
          else
            FSET(obj, 'add_isotropic_linearized_elasticity_brick_pstress', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_linear_incompressibility_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_linear_incompressibility_brick', index(2).subs{:});
          else
            FSET(obj, 'add_linear_incompressibility_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_nonlinear_elasticity_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_nonlinear_elasticity_brick', index(2).subs{:});
          else
            FSET(obj, 'add_nonlinear_elasticity_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_finite_strain_elasticity_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_finite_strain_elasticity_brick', index(2).subs{:});
          else
            FSET(obj, 'add_finite_strain_elasticity_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_small_strain_elastoplasticity_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_small_strain_elastoplasticity_brick', index(2).subs{:});
          else
            FSET(obj, 'add_small_strain_elastoplasticity_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_elastoplasticity_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_elastoplasticity_brick', index(2).subs{:});
          else
            FSET(obj, 'add_elastoplasticity_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_finite_strain_elastoplasticity_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_finite_strain_elastoplasticity_brick', index(2).subs{:});
          else
            FSET(obj, 'add_finite_strain_elastoplasticity_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_nonlinear_incompressibility_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_nonlinear_incompressibility_brick', index(2).subs{:});
          else
            FSET(obj, 'add_nonlinear_incompressibility_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_finite_strain_incompressibility_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_finite_strain_incompressibility_brick', index(2).subs{:});
          else
            FSET(obj, 'add_finite_strain_incompressibility_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_bilaplacian_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_bilaplacian_brick', index(2).subs{:});
          else
            FSET(obj, 'add_bilaplacian_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Kirchhoff_Love_plate_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Kirchhoff_Love_plate_brick', index(2).subs{:});
          else
            FSET(obj, 'add_Kirchhoff_Love_plate_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_normal_derivative_source_term_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_normal_derivative_source_term_brick', index(2).subs{:});
          else
            FSET(obj, 'add_normal_derivative_source_term_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Kirchhoff_Love_Neumann_term_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Kirchhoff_Love_Neumann_term_brick', index(2).subs{:});
          else
            FSET(obj, 'add_Kirchhoff_Love_Neumann_term_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_normal_derivative_Dirichlet_condition_with_multipliers'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_normal_derivative_Dirichlet_condition_with_multipliers', index(2).subs{:});
          else
            FSET(obj, 'add_normal_derivative_Dirichlet_condition_with_multipliers', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_normal_derivative_Dirichlet_condition_with_penalization'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_normal_derivative_Dirichlet_condition_with_penalization', index(2).subs{:});
          else
            FSET(obj, 'add_normal_derivative_Dirichlet_condition_with_penalization', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Mindlin_Reissner_plate_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Mindlin_Reissner_plate_brick', index(2).subs{:});
          else
            FSET(obj, 'add_Mindlin_Reissner_plate_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_mass_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_mass_brick', index(2).subs{:});
          else
            FSET(obj, 'add_mass_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'shift_variables_for_time_integration'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'shift_variables_for_time_integration', index(2).subs{:});
          else
            FSET(obj, 'shift_variables_for_time_integration', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'perform_init_time_derivative'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'perform_init_time_derivative', index(2).subs{:});
          else
            FSET(obj, 'perform_init_time_derivative', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_time_step'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'set_time_step', index(2).subs{:});
          else
            FSET(obj, 'set_time_step', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set_time'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'set_time', index(2).subs{:});
          else
            FSET(obj, 'set_time', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_theta_method_for_first_order'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_theta_method_for_first_order', index(2).subs{:});
          else
            FSET(obj, 'add_theta_method_for_first_order', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_theta_method_for_second_order'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_theta_method_for_second_order', index(2).subs{:});
          else
            FSET(obj, 'add_theta_method_for_second_order', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Newmark_scheme'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Newmark_scheme', index(2).subs{:});
          else
            FSET(obj, 'add_Newmark_scheme', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'disable_bricks'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'disable_bricks', index(2).subs{:});
          else
            FSET(obj, 'disable_bricks', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'enable_bricks'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'enable_bricks', index(2).subs{:});
          else
            FSET(obj, 'enable_bricks', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'disable_variable'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'disable_variable', index(2).subs{:});
          else
            FSET(obj, 'disable_variable', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'enable_variable'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'enable_variable', index(2).subs{:});
          else
            FSET(obj, 'enable_variable', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'first_iter'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'first_iter', index(2).subs{:});
          else
            FSET(obj, 'first_iter', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'next_iter'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'next_iter', index(2).subs{:});
          else
            FSET(obj, 'next_iter', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_basic_contact_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_basic_contact_brick', index(2).subs{:});
          else
            FSET(obj, 'add_basic_contact_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_basic_contact_brick_two_deformable_bodies'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_basic_contact_brick_two_deformable_bodies', index(2).subs{:});
          else
            FSET(obj, 'add_basic_contact_brick_two_deformable_bodies', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'contact_brick_set_BN'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'contact_brick_set_BN', index(2).subs{:});
          else
            FSET(obj, 'contact_brick_set_BN', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'contact_brick_set_BT'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'contact_brick_set_BT', index(2).subs{:});
          else
            FSET(obj, 'contact_brick_set_BT', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_nodal_contact_with_rigid_obstacle_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_nodal_contact_with_rigid_obstacle_brick', index(2).subs{:});
          else
            FSET(obj, 'add_nodal_contact_with_rigid_obstacle_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_contact_with_rigid_obstacle_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_contact_with_rigid_obstacle_brick', index(2).subs{:});
          else
            FSET(obj, 'add_contact_with_rigid_obstacle_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_integral_contact_with_rigid_obstacle_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_integral_contact_with_rigid_obstacle_brick', index(2).subs{:});
          else
            FSET(obj, 'add_integral_contact_with_rigid_obstacle_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_penalized_contact_with_rigid_obstacle_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_penalized_contact_with_rigid_obstacle_brick', index(2).subs{:});
          else
            FSET(obj, 'add_penalized_contact_with_rigid_obstacle_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Nitsche_contact_with_rigid_obstacle_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Nitsche_contact_with_rigid_obstacle_brick', index(2).subs{:});
          else
            FSET(obj, 'add_Nitsche_contact_with_rigid_obstacle_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Nitsche_midpoint_contact_with_rigid_obstacle_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Nitsche_midpoint_contact_with_rigid_obstacle_brick', index(2).subs{:});
          else
            FSET(obj, 'add_Nitsche_midpoint_contact_with_rigid_obstacle_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Nitsche_fictitious_domain_contact_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Nitsche_fictitious_domain_contact_brick', index(2).subs{:});
          else
            FSET(obj, 'add_Nitsche_fictitious_domain_contact_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_nodal_contact_between_nonmatching_meshes_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_nodal_contact_between_nonmatching_meshes_brick', index(2).subs{:});
          else
            FSET(obj, 'add_nodal_contact_between_nonmatching_meshes_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_nonmatching_meshes_contact_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_nonmatching_meshes_contact_brick', index(2).subs{:});
          else
            FSET(obj, 'add_nonmatching_meshes_contact_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_integral_contact_between_nonmatching_meshes_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_integral_contact_between_nonmatching_meshes_brick', index(2).subs{:});
          else
            FSET(obj, 'add_integral_contact_between_nonmatching_meshes_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_penalized_contact_between_nonmatching_meshes_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_penalized_contact_between_nonmatching_meshes_brick', index(2).subs{:});
          else
            FSET(obj, 'add_penalized_contact_between_nonmatching_meshes_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_integral_large_sliding_contact_brick_raytracing'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_integral_large_sliding_contact_brick_raytracing', index(2).subs{:});
          else
            FSET(obj, 'add_integral_large_sliding_contact_brick_raytracing', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_rigid_obstacle_to_large_sliding_contact_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_rigid_obstacle_to_large_sliding_contact_brick', index(2).subs{:});
          else
            FSET(obj, 'add_rigid_obstacle_to_large_sliding_contact_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_master_contact_boundary_to_large_sliding_contact_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_master_contact_boundary_to_large_sliding_contact_brick', index(2).subs{:});
          else
            FSET(obj, 'add_master_contact_boundary_to_large_sliding_contact_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_slave_contact_boundary_to_large_sliding_contact_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_slave_contact_boundary_to_large_sliding_contact_brick', index(2).subs{:});
          else
            FSET(obj, 'add_slave_contact_boundary_to_large_sliding_contact_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_master_slave_contact_boundary_to_large_sliding_contact_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_master_slave_contact_boundary_to_large_sliding_contact_brick', index(2).subs{:});
          else
            FSET(obj, 'add_master_slave_contact_boundary_to_large_sliding_contact_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_Nitsche_large_sliding_contact_brick_raytracing'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_Nitsche_large_sliding_contact_brick_raytracing', index(2).subs{:});
          else
            FSET(obj, 'add_Nitsche_large_sliding_contact_brick_raytracing', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_rigid_obstacle_to_Nitsche_large_sliding_contact_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_rigid_obstacle_to_Nitsche_large_sliding_contact_brick', index(2).subs{:});
          else
            FSET(obj, 'add_rigid_obstacle_to_Nitsche_large_sliding_contact_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_master_contact_boundary_to_biased_Nitsche_large_sliding_contact_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_master_contact_boundary_to_biased_Nitsche_large_sliding_contact_brick', index(2).subs{:});
          else
            FSET(obj, 'add_master_contact_boundary_to_biased_Nitsche_large_sliding_contact_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_slave_contact_boundary_to_biased_Nitsche_large_sliding_contact_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_slave_contact_boundary_to_biased_Nitsche_large_sliding_contact_brick', index(2).subs{:});
          else
            FSET(obj, 'add_slave_contact_boundary_to_biased_Nitsche_large_sliding_contact_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'add_contact_boundary_to_unbiased_Nitsche_large_sliding_contact_brick'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, 'add_contact_boundary_to_unbiased_Nitsche_large_sliding_contact_brick', index(2).subs{:});
          else
            FSET(obj, 'add_contact_boundary_to_unbiased_Nitsche_large_sliding_contact_brick', index(2).subs{:});
            if (exist('ans', 'var') == 1)
              varargout{1}=ans;
            end;
          end;
          return;
        case 'set'
          if (nargout)
            [varargout{1:nargout}] = FSET(obj, index(2).subs{:});
          else
            FSET(obj,index(2).subs{:});
            if (exist('ans', 'var') == 1)
              h=ans;
              if (isstruct(h) & isfield(h,'id') & isfield(h,'cid')), h = gfObject(h); end;
              varargout{1}=h;
            end;
          end;
          return;
        case 'get'
          if (nargout) 
            h = FGET(obj, index(2).subs{:});
            if (isstruct(h) & isfield(h,'id') & isfield(h,'cid')), h = gfObject(h); end;
            [varargout{1:nargout}] = h;
          else
	     FGET(obj,index(2).subs{:});
            if (exist('ans', 'var') == 1)
              h=ans;
              if (isstruct(h) & isfield(h,'id') & isfield(h,'cid')), h = gfObject(h); end;
              varargout{1}=h;
            end;
          end;
          return;
        otherwise
          if ((numel(index) > 1) && (strcmp(index(2).type, '()')))
            h = FGET(obj,index(1).subs, index(2).subs{:});
            if (isstruct(h) & isfield(h,'id') & isfield(h,'cid')), h = gfObject(h); end;
            [varargout{1:nargout}] = h;
            cnt = cnt + 1;
          else
            h = FGET(obj, index(1).subs);
            if (isstruct(h) & isfield(h,'id') & isfield(h,'cid')), h = gfObject(h); end;
            [varargout{1:nargout}] = h;
          end
      end
  end
  % if there are others indexes, let matlab do its work
  if (numel(index) > cnt)
    for i=1:nout,
      varargout{i} = subsref(varargout{i}, index((cnt+1):end));
    end;
  end;
% autogenerated mfile;
