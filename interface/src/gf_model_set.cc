/*===========================================================================

 Copyright (C) 2009-2020 Yves Renard.

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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

/**\file gf_model_set.cc
   \brief getfemint_model setter.
*/

#include <getfem/getfem_im_data.h>
#include <getfem/getfem_contact_and_friction_nodal.h>
#include <getfem/getfem_contact_and_friction_integral.h>
#include <getfem/getfem_contact_and_friction_large_sliding.h>
#include <getfem/getfem_nonlinear_elasticity.h>
#include <getfem/getfem_plasticity.h>
#include <getfem/getfem_HHO.h>
#include <getfem/getfem_fourth_order.h>
#include <getfem/getfem_linearized_plates.h>
#include <getfemint_misc.h>
#include <getfemint_workspace.h>
#include <getfemint_gsparse.h>

using namespace getfemint;


/*@GFDOC
  Modifies a model object.
@*/


// Object for the declaration of a new sub-command.

struct sub_gf_md_set : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
                   getfemint::mexargs_out& out,
                   getfem::model *md) = 0;
};

typedef std::shared_ptr<sub_gf_md_set> psub_command;

typedef std::map<size_type, size_type> elt_corr_cont;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_md_set {                                \
      virtual void run(getfemint::mexargs_in& in,                       \
                       getfemint::mexargs_out& out,                     \
                       getfem::model *md)                               \
      { dummy_func(in); dummy_func(out); code }                         \
    };                                                                  \
    psub_command psubc = std::make_shared<subc>();                      \
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;         \
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;     \
    subc_tab[cmd_normalize(name)] = psubc;                              \
  }

static void filter_lawname(std::string &lawname) {
  for (auto &c : lawname)
    { if (c == ' ') c = '_'; if (c >= 'A' && c <= 'Z') c = char(c+'a'-'A'); }
}

void gf_model_set(getfemint::mexargs_in& m_in,
                  getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {

    /*@SET ('clear')
      Clear the model.@*/
    sub_command
      ("clear", 0, 0, 0, 1,
       md->clear();
       );


    /*@SET ('add fem variable', @str name, @tmf mf)
      Add a variable to the model linked to a @tmf. `name` is the variable
      name. @*/
    sub_command
      ("add fem variable", 2, 2, 0, 0,
       std::string name = in.pop().to_string();
       getfem::mesh_fem *mf = to_meshfem_object(in.pop());
       md->add_fem_variable(name, *mf);
       workspace().set_dependence(md, mf);
       );

    /*@SET ('add filtered fem variable', @str name, @tmf mf, @int region)
      Add a variable to the model linked to a @tmf. The variable is filtered
      in the sense that only the dof on the region are considered.
      `name` is the variable name. @*/
    sub_command
      ("add filtered fem variable", 3, 3, 0, 0,
       std::string name = in.pop().to_string();
       getfem::mesh_fem *mf = to_meshfem_object(in.pop());
       size_type region = in.pop().to_integer();
       md->add_filtered_fem_variable(name, *mf, region);
       workspace().set_dependence(md, mf);
       );


    /*@SET ('add im variable', @str name, @tmimd mimd)
      Add a variable to the model linked to a @tmimd. `name` is the variable
      name. @*/
    sub_command
      ("add im variable", 2, 2, 0, 0,
       std::string name = in.pop().to_string();
       getfem::im_data *mimd = to_meshimdata_object(in.pop());
       md->add_im_variable(name, *mimd);
       workspace().set_dependence(md, mimd);
       );

    /*@SET ('add internal im variable', @str name, @tmimd mimd)
      Add a variable to the model, which is linked to a @tmimd and will be
      condensed out during the assemblage of the tangent matrix. `name` is
      the variable name. @*/
    sub_command
      ("add internal im variable", 2, 2, 0, 0,
       std::string name = in.pop().to_string();
       getfem::im_data *mimd = to_meshimdata_object(in.pop());
       md->add_internal_im_variable(name, *mimd);
       workspace().set_dependence(md, mimd);
       );

    /*@SET ('add variable', @str name, sizes)
      Add a variable to the model of constant sizes. `sizes` is either a
      integer (for a scalar or vector variable) or a vector of dimensions
      for a tensor variable. `name` is the variable name. @*/
    sub_command
      ("add variable", 2, 2, 0, 0,
       std::string name = in.pop().to_string();
       mexarg_in argin = in.pop();
       bgeot::multi_index mi(1);
       if (argin.is_integer()) {
         mi[0] = argin.to_integer();
       } else {
         iarray v = argin.to_iarray();
         mi.resize(v.size());
         for (size_type i = 0; i < v.size(); ++i) mi[i] = v[i];
       }
       md->add_fixed_size_variable(name, mi);
       );

    /*@SET ('delete variable', @str name)
      Delete a variable or a data from the model. @*/
    sub_command
      ("delete variable", 1, 1, 0, 0,
       std::string name = in.pop().to_string();
       md->delete_variable(name);
       );


    /*@SET ('resize variable', @str name, sizes)
      Resize a  constant size variable of the model.  `sizes` is either a
      integer (for a scalar or vector variable) or a vector of dimensions
      for a tensor variable. `name` is the variable name. @*/
    sub_command
      ("resize variable", 2, 2, 0, 0,
       std::string name = in.pop().to_string();
       mexarg_in argin = in.pop();
       bgeot::multi_index mi(1);
       if (argin.is_integer()) {
         mi[0] = argin.to_integer();
       } else {
         iarray v = argin.to_iarray();
         mi.resize(v.size());
         for (size_type i = 0; i < v.size(); ++i) mi[i] = v[i];
       }
       md->resize_fixed_size_variable(name, mi);
       );


    /*@SET ('add multiplier', @str name, @tmf mf, @str primalname[, @tmim mim, @int region])
    Add a particular variable linked to a fem being a multiplier with
    respect to a primal variable. The dof will be filtered with the
    ``gmm::range_basis`` function applied on the terms of the model
    which link the multiplier and the primal variable. This in order to
    retain only linearly independent constraints on the primal variable.
    Optimized for boundary multipliers. @*/
    sub_command
      ("add multiplier", 3, 5, 0, 0,
       std::string name = in.pop().to_string();
       getfem::mesh_fem *mf = to_meshfem_object(in.pop());
       std::string primalname = in.pop().to_string();

       getfem::mesh_im *mim = 0;
       size_type region = size_type(-1);
       if (in.remaining()) {
         mexarg_in argin = in.pop();
         mim = to_meshim_object(argin);
         region = in.pop().to_integer();
       }
       if (mim)
         md->add_multiplier(name, *mf, primalname, *mim, region);
       else
         md->add_multiplier(name, *mf, primalname);
       workspace().set_dependence(md, mf);
       );


    /*@SET ('add im data', @str name, @tmimd mimd)
      Add a data set to the model linked to a @tmimd. `name` is the data
      name. @*/
    sub_command
      ("add im data", 2, 2, 0, 0,
       std::string name = in.pop().to_string();
       getfem::im_data *mimd = to_meshimdata_object(in.pop());
       md->add_im_data(name, *mimd);
       workspace().set_dependence(md, mimd);
       );


    /*@SET ('add fem data', @str name, @tmf mf[, sizes])
      Add a data to the model linked to a @tmf. `name` is the data name,
      `sizes` an optional parameter which is either an 
      integer  or a vector of suplementary dimensions with respect to `mf`. @*/
    sub_command
      ("add fem data", 2, 3, 0, 0,
       std::string name = in.pop().to_string();
       getfem::mesh_fem *mf = to_meshfem_object(in.pop());
       bgeot::multi_index mi(1); mi[0] = 1;
       if (in.remaining()) {
         mexarg_in argin = in.pop();
         if (argin.is_integer()) {
           mi[0] = argin.to_integer();
         } else {
           iarray v = argin.to_iarray();
           mi.resize(v.size());
           for (size_type i = 0; i < v.size(); ++i) mi[i] = v[i];
         }
       }
       md->add_fem_data(name, *mf, mi);
       workspace().set_dependence(md, mf);
       );


    /*@SET ('add initialized fem data', @str name, @tmf mf, @vec V[, sizes])
      Add a data to the model linked to a @tmf. `name` is the data name.
      The data is initiakized with `V`. The data can be a scalar or vector
      field. `sizes` an optional parameter which is either an 
      integer or a vector of suplementary dimensions with respect to `mf`.@*/
    sub_command
      ("add initialized fem data", 3, 4, 0, 0,
       std::string name = in.pop().to_string();
       getfem::mesh_fem *mf = to_meshfem_object(in.pop());
       if (!md->is_complex()) {
         darray st = in.pop().to_darray();
         std::vector<double> V(st.begin(), st.end());
         bgeot::multi_index mi(1);
         mi[0] = gmm::vect_size(V) / mf->nb_dof();
         if (in.remaining()) {
           mexarg_in argin = in.pop();
           if (argin.is_integer()) {
             mi[0] = argin.to_integer();
           } else {
             iarray v = argin.to_iarray();
             mi.resize(v.size());
             for (size_type i = 0; i < v.size(); ++i) mi[i] = v[i];
           }
         }
         md->add_initialized_fem_data(name, *mf, V, mi);
       } else {
         carray st = in.pop().to_carray();
         std::vector<std::complex<double> > V(st.begin(), st.end());
         bgeot::multi_index mi(1);
         mi[0] = gmm::vect_size(V) / mf->nb_dof();
         if (in.remaining()) {
           mexarg_in argin = in.pop();
           if (argin.is_integer()) {
             mi[0] = argin.to_integer();
           } else {
             iarray v = argin.to_iarray();
             mi.resize(v.size());
             for (size_type i = 0; i < v.size(); ++i) mi[i] = v[i];
           }
         }
         md->add_initialized_fem_data(name, *mf, V, mi);
       }
       workspace().set_dependence(md, mf);
       );


    /*@SET ('add data', @str name, @int size)
      Add a fixed size data to the model.  `sizes` is either a
      integer (for a scalar or vector data) or a vector of dimensions
      for a tensor data. `name` is the data name. @*/
    sub_command
      ("add data", 2, 2, 0, 0,
       std::string name = in.pop().to_string();
       mexarg_in argin = in.pop();
       bgeot::multi_index mi(1);
       if (argin.is_integer()) {
         mi[0] = argin.to_integer();
       } else {
         iarray v = argin.to_iarray();
         mi.resize(v.size());
         for (size_type i = 0; i < v.size(); ++i) mi[i] = v[i];
       }
       md->add_fixed_size_data(name, mi);
       );

    /*@SET ('add macro', @str name, @str expr)
      Define a new macro for the high generic assembly language.
      The name include the parameters. For instance name='sp(a,b)', expr='a.b'
      is a valid definition. Macro without parameter can also be defined.
      For instance name='x1', expr='X[1]' is valid. The form name='grad(u)',
      expr='Grad_u' is also allowed but in that case, the parameter 'u' will
      only be allowed to be a variable name when using the macro. Note that
      macros can be directly defined inside the assembly strings with the
      keyword 'Def'.
      @*/
    sub_command
      ("add macro", 2, 2, 0, 0,
       std::string name = in.pop().to_string();
       std::string expr = in.pop().to_string();
       md->add_macro(name, expr);
       );

    /*@SET ('del macro', @str name)
      Delete a previously defined macro for the high generic assembly language.
      @*/
    sub_command
      ("del macro", 1, 1, 0, 0,
       std::string name = in.pop().to_string();
       md->del_macro(name);
       );

    /*@SET ('add initialized data', @str name, @vec V[, sizes])
      Add an initialized fixed size data to the model. `sizes` an
      optional parameter which is either an 
      integer  or a vector dimensions that describes the format of the
      data. By default, the data is considered to b a vector field.
      `name` is the data name and `V` is the value of the data.@*/
    sub_command
      ("add initialized data", 2, 3, 0, 0,
       std::string name = in.pop().to_string();
       if (!md->is_complex()) {
         darray st = in.pop().to_darray();
         std::vector<double> V(st.begin(), st.end());
         bgeot::multi_index mi(1);
         mi[0] = gmm::vect_size(V);
         if (in.remaining()) {
           mexarg_in argin = in.pop();
           if (argin.is_integer()) {
             mi[0] = argin.to_integer();
           } else {
             iarray v = argin.to_iarray();
             mi.resize(v.size());
             for (size_type i = 0; i < v.size(); ++i) mi[i] = v[i];
           }
         }
         md->add_initialized_fixed_size_data(name, V, mi);
       } else {
         carray st = in.pop().to_carray();
         std::vector<std::complex<double> > V(st.begin(), st.end());
         bgeot::multi_index mi(1);
         mi[0] = gmm::vect_size(V);
         if (in.remaining()) {
           mexarg_in argin = in.pop();
           if (argin.is_integer()) {
             mi[0] = argin.to_integer();
           } else {
             iarray v = argin.to_iarray();
             mi.resize(v.size());
             for (size_type i = 0; i < v.size(); ++i) mi[i] = v[i];
           }
         }
         md->add_initialized_fixed_size_data(name, V, mi);
       }
       );


    /*@SET ('variable', @str name, @vec V)
      Set the value of a variable or data. `name` is the data name.@*/
    sub_command
      ("variable", 2, 2, 0, 0,
       std::string name = in.pop().to_string();
       if (!md->is_complex()) {
         darray st = in.pop().to_darray();
         GMM_ASSERT1(st.size()==md->real_variable(name).size(),
                     "Bad size in assignment");
         md->set_real_variable(name).assign(st.begin(),st.end());
       } else {
         carray st = in.pop().to_carray();
         GMM_ASSERT1(st.size() == md->complex_variable(name).size(),
                     "Bad size in assignment");
         md->set_complex_variable(name).assign(st.begin(),
                                                              st.end());
       }
       );


    /*@SET ('to variables', @vec V)
      Set the value of the variables of the model with the vector `V`.
      Typically, the vector `V` results of the solve of the tangent
      linear system (useful to solve your problem with you own solver).@*/
    sub_command
      ("to variables", 1, 1, 0, 0,
       if (!md->is_complex()) {
         darray st = in.pop().to_darray(-1);
         std::vector<double> V;
         V.assign(st.begin(), st.end());
         md->to_variables(V);
       } else {
         carray st = in.pop().to_carray(-1);
         std::vector<std::complex<double> > V;
         V.assign(st.begin(), st.end());
         md->to_variables(V);
       }
       );

    /*@SET ('delete brick', @int ind_brick)
      Delete a variable or a data from the model. @*/
    sub_command
      ("delete brick", 1, 1, 0, 0,
       size_type ib = in.pop().to_integer() - config::base_index();
       md->delete_brick(ib);
       );

    /*@SET ('define variable group', @str name[, @str varname, ...])
      Defines a group of variables for the interpolation (mainly for the
      raytracing interpolation transformation.@*/
    sub_command
      ("define variable group", 1, 1000, 0, 0,
       std::string name = in.pop().to_string();
       std::vector<std::string> nl;
       while (in.remaining()) nl.push_back(in.pop().to_string());
       md->define_variable_group(name, nl);
       );

    /*@SET ('add elementary rotated RT0 projection', @str transname)
      Add the elementary transformation corresponding to the projection
      on rotated RT0 element for two-dimensional elements to the model.
      The name is the name given to the elementary transformation. @*/
    sub_command
      ("add elementary rotated RT0 projection", 1, 1, 0, 0,
       std::string transname = in.pop().to_string();
       add_2D_rotated_RT0_projection(*md, transname);
       );

    /*@SET ('add elementary P0 projection', @str transname)
      Add the elementary transformation corresponding to the projection
      P0 element.
      The name is the name given to the elementary transformation. @*/
    sub_command
      ("add P0 projection", 1, 1, 0, 0,
       std::string transname = in.pop().to_string();
       add_P0_projection(*md, transname);
       );

    /*@SET ('add HHO reconstructed gradient', @str transname)
      Add to the model the elementary transformation corresponding to the
      reconstruction of a gradient for HHO methods.
      The name is the name given to the elementary transformation. @*/
    sub_command
      ("add HHO reconstructed gradient", 1, 1, 0, 0,
       std::string transname = in.pop().to_string();
       add_HHO_reconstructed_gradient(*md, transname);
       );

    /*@SET ('add HHO reconstructed symmetrized gradient', @str transname)
      Add to the model the elementary transformation corresponding to the
      reconstruction of a symmetrized gradient for HHO methods.
      The name is the name given to the elementary transformation. @*/
    sub_command
      ("add HHO reconstructed symmetrized gradient", 1, 1, 0, 0,
       std::string transname = in.pop().to_string();
       add_HHO_reconstructed_symmetrized_gradient(*md, transname);
       );

    /*@SET ('add HHO reconstructed value', @str transname)
      Add to the model the elementary transformation corresponding to the
      reconstruction of the variable for HHO methods.
      The name is the name given to the elementary transformation. @*/
    sub_command
      ("add HHO reconstructed value", 1, 1, 0, 0,
       std::string transname = in.pop().to_string();
       add_HHO_reconstructed_value(*md, transname);
       );

    /*@SET ('add HHO reconstructed symmetrized value', @str transname)
      Add to the model the elementary transformation corresponding to the
      reconstruction of the variable for HHO methods using a symmetrized
      gradient.
      The name is the name given to the elementary transformation. @*/
    sub_command
      ("add HHO reconstructed symmetrized value", 1, 1, 0, 0,
       std::string transname = in.pop().to_string();
       add_HHO_reconstructed_symmetrized_value(*md, transname);
       );

    /*@SET ('add HHO stabilization', @str transname)
      Add to the model the elementary transformation corresponding to the
      HHO stabilization operator.
      The name is the name given to the elementary transformation. @*/
    sub_command
      ("add HHO stabilization", 1, 1, 0, 0,
       std::string transname = in.pop().to_string();
       add_HHO_stabilization(*md, transname);
       );

    /*@SET ('add HHO symmetrized stabilization', @str transname)
      Add to the model the elementary transformation corresponding to the
      HHO stabilization operator using a symmetrized gradient.
      The name is the name given to the elementary transformation. @*/
    sub_command
      ("add HHO symmetrized stabilization", 1, 1, 0, 0,
       std::string transname = in.pop().to_string();
       add_HHO_symmetrized_stabilization(*md, transname);
       );


    /*@SET ('add interpolate transformation from expression', @str transname, @tmesh source_mesh, @tmesh target_mesh, @str expr)
      Add a transformation to the model from mesh `source_mesh` to mesh
      `target_mesh` given by the expression `expr` which corresponds to a
      high-level generic assembly expression which may contains some
      variable of the model. CAUTION: the derivative of the
      transformation with used variable is taken into account in the
      computation of the tangen system. However, order two derivative is not
      implemented, so such tranformation is not allowed in the definition
      of a potential. @*/
    sub_command
      ("add interpolate transformation from expression", 4, 4, 0, 0,
       std::string transname = in.pop().to_string();
       getfem::mesh *sm = extract_mesh_object(in.pop());
       getfem::mesh *tm = extract_mesh_object(in.pop());
       std::string expr = in.pop().to_string();
       add_interpolate_transformation_from_expression(*md, transname, *sm,
                                                      *tm, expr);
       );

    /*@SET ('add element extrapolation transformation', @str transname, @tmesh source_mesh, @mat elt_corr)
      Add a special interpolation transformation which represents the identity
      transformation but allows to evaluate the expression on another element
      than the current element by polynomial extrapolation. It is used for
      stabilization term in fictitious domain applications. the array elt_cor
      should be a two entry array whose first line contains the elements
      concerned by the transformation and the second line the respective
      elements on which the extrapolation has to be made. If an element
      is not listed in elt_cor the evaluation is just made on the current
      element. @*/
    sub_command
      ("add element extrapolation transformation", 3, 3, 0, 0,
       std::string transname = in.pop().to_string();
       getfem::mesh *sm = extract_mesh_object(in.pop());
       iarray v = in.pop().to_iarray();
       if (v.getm() != 2 || v.getp() != 1 || v.getq() != 1)
       	 THROW_BADARG("Invalid format for the convex correspondence list");
       elt_corr_cont elt_corr;
       for (size_type j=0; j < v.getn(); j++)
	 elt_corr[v(0,j)-config::base_index()] = v(1,j)-config::base_index();
       getfem::add_element_extrapolation_transformation(*md, transname, *sm,
							elt_corr);
       );

    /*@SET ('add standard secondary domain', @str name, @tmim mim, @int region = -1)
      Add a secondary domain to the model which can be used in a weak-form language expression for integration on the product of two domains. `name` is the name
      of the secondary domain, `mim` is an integration method on this domain
      and `region` the region on which the integration is to be performed. @*/
    sub_command
      ("add standard secondary domain", 3, 3, 0, 0,
       std::string name = in.pop().to_string();
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       add_standard_secondary_domain(*md, name, *mim, region);
       );


    /*@SET ('set element extrapolation correspondence', @str transname, @mat elt_corr)
      Change the correspondence map of an element extrapolation interpolate
     transformation. @*/
    sub_command
      ("set element extrapolation correspondence", 2, 2, 0, 0,
       std::string transname = in.pop().to_string();
       iarray v = in.pop().to_iarray();
       if (v.getm() != 2 || v.getp() != 1 || v.getq() != 1)
       	 THROW_BADARG("Invalid format for the convex correspondence list");
       elt_corr_cont elt_corr;
       for (size_type j=0; j < v.getn(); j++)
	 elt_corr[v(0,j)-config::base_index()] = v(1,j)-config::base_index();
       getfem::set_element_extrapolation_correspondence(*md, transname,
							elt_corr);
       );
    
    /*@SET ('add raytracing transformation', @str transname, @scalar release_distance)
      Add a raytracing interpolate transformation called `transname` to a model
      to be used by the generic assembly bricks.
      CAUTION: For the moment, the derivative of the
      transformation is not taken into account in the model solve. @*/
    sub_command
      ("add raytracing transformation", 2, 2, 0, 0,
       std::string transname = in.pop().to_string();
       scalar_type d = in.pop().to_scalar();
       add_raytracing_transformation(*md, transname, d);
       );

    /*@SET ('add master contact boundary to raytracing transformation', @str transname, @tmesh m, @str dispname, @int region)
      Add a master contact boundary with corresponding displacement variable
      `dispname` on a specific boundary `region` to an existing raytracing
      interpolate transformation called `transname`. @*/
    sub_command
      ("add master contact boundary to raytracing transformation", 4, 4, 0, 0,
       std::string transname = in.pop().to_string();
       getfem::mesh *sm = extract_mesh_object(in.pop());
       std::string dispname = in.pop().to_string();
       size_type region = in.pop().to_integer();
       add_master_contact_boundary_to_raytracing_transformation
       (*md, transname, *sm, dispname, region);
       );

    /*@SET ('add slave contact boundary to raytracing transformation', @str transname, @tmesh m, @str dispname, @int region)
      Add a slave contact boundary with corresponding displacement variable
      `dispname` on a specific boundary `region` to an existing raytracing
      interpolate transformation called `transname`. @*/
    sub_command
      ("add slave contact boundary to raytracing transformation", 4, 4, 0, 0,
       std::string transname = in.pop().to_string();
       getfem::mesh *sm = extract_mesh_object(in.pop());
       std::string dispname = in.pop().to_string();
       size_type region = in.pop().to_integer();
       add_slave_contact_boundary_to_raytracing_transformation
       (*md, transname, *sm, dispname, region);
       );

    /*@SET ('add rigid obstacle to raytracing transformation', @str transname, @str expr, @int N)
      Add a rigid obstacle whose geometry corresponds to the zero level-set
      of the high-level generic assembly expression `expr`
      to an existing raytracing interpolate transformation called `transname`.
      @*/
    sub_command
      ("add rigid obstacle to raytracing transformation", 3, 3, 0, 0,
       std::string transname = in.pop().to_string();
       std::string expr = in.pop().to_string();
       size_type N = in.pop().to_integer();
       add_rigid_obstacle_to_raytracing_transformation
       (*md, transname, expr, N);
       );
/*@SET ('add projection transformation', @str transname, @scalar release_distance)
      Add a projection interpolate transformation called `transname` to a model
      to be used by the generic assembly bricks.
      CAUTION: For the moment, the derivative of the
      transformation is not taken into account in the model solve. @*/  
      sub_command
      ("add projection transformation", 2, 2, 0, 0,
       std::string transname = in.pop().to_string();
       scalar_type d = in.pop().to_scalar();
       add_projection_transformation(*md, transname, d);
       );

    /*@SET ('add master contact boundary to projection transformation', @str transname, @tmesh m, @str dispname, @int region)
      Add a master contact boundary with corresponding displacement variable
      `dispname` on a specific boundary `region` to an existing projection
      interpolate transformation called `transname`. @*/
    sub_command
      ("add master contact boundary to projection transformation", 4, 4, 0, 0,
       std::string transname = in.pop().to_string();
       getfem::mesh *sm = extract_mesh_object(in.pop());
       std::string dispname = in.pop().to_string();
       size_type region = in.pop().to_integer();
       add_master_contact_boundary_to_projection_transformation
       (*md, transname, *sm, dispname, region);
       );

    /*@SET ('add slave contact boundary to projection transformation', @str transname, @tmesh m, @str dispname, @int region)
      Add a slave contact boundary with corresponding displacement variable
      `dispname` on a specific boundary `region` to an existing projection
      interpolate transformation called `transname`. @*/
    sub_command
      ("add slave contact boundary to projection transformation", 4, 4, 0, 0,
       std::string transname = in.pop().to_string();
       getfem::mesh *sm = extract_mesh_object(in.pop());
       std::string dispname = in.pop().to_string();
       size_type region = in.pop().to_integer();
       add_slave_contact_boundary_to_projection_transformation
       (*md, transname, *sm, dispname, region);
       );

    /*@SET ('add rigid obstacle to projection transformation', @str transname, @str expr, @int N)
      Add a rigid obstacle whose geometry corresponds to the zero level-set
      of the high-level generic assembly expression `expr`
      to an existing projection interpolate transformation called `transname`.
      @*/
    sub_command
      ("add rigid obstacle to projection transformation", 3, 3, 0, 0,
       std::string transname = in.pop().to_string();
       std::string expr = in.pop().to_string();
       size_type N = in.pop().to_integer();
       add_rigid_obstacle_to_projection_transformation
       (*md, transname, expr, N);
       );
      
    /*@SET ind = ('add linear term', @tmim mim, @str expression[, @int region[, @int is_symmetric[, @int is_coercive]]])
      Adds a matrix term given by the assembly string `expr` which will
      be assembled in region `region` and with the integration method `mim`.
      Only the matrix term will be taken into account, assuming that it is
      linear.
      The advantage of declaring a term linear instead of nonlinear is that
      it will be assembled only once and no assembly is necessary for the
      residual.
      Take care that if the expression contains some variables and if the
      expression is a potential or of first order (i.e. describe the weak
      form, not the derivative of the weak form), the expression will be
      derivated with respect to all variables.
      You can specify if the term is symmetric, coercive or not.
      If you are not sure, the better is to declare the term not symmetric
      and not coercive. But some solvers (conjugate gradient for instance)
      are not allowed for non-coercive problems.
      `brickname` is an optional name for the brick.@*/
    sub_command
      ("add linear term", 2, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string expr = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       int is_symmetric = 0;
       if (in.remaining()) is_symmetric = in.pop().to_integer();
       int is_coercive = 0;
       if (in.remaining()) is_coercive = in.pop().to_integer();
       
       size_type ind = getfem::add_linear_term
       (*md, *mim, expr, region, is_symmetric,
        is_coercive) + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add linear twodomain term', @tmim mim, @str expression, @int region, @str secondary_domain[, @int is_symmetric[, @int is_coercive]])
      Adds a linear term given by a weak form language expression like
      MODEL:SET('add linear term') but for an integration on a direct
      product of two domains, a first specfied by ``mim`` and ``region``
      and a second one by ``secondary_domain`` which has to be declared
      first into the model.@*/
    sub_command
      ("add linear twodomain term", 4, 6, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string expr = in.pop().to_string();
       size_type region = in.pop().to_integer();
       std::string secdom = in.pop().to_string();
       int is_symmetric = 0;
       if (in.remaining()) is_symmetric = in.pop().to_integer();
       int is_coercive = 0;
       if (in.remaining()) is_coercive = in.pop().to_integer();
       
       size_type ind = getfem::add_linear_twodomain_term
       (*md, *mim, expr, region, secdom, is_symmetric,
        is_coercive) + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add linear generic assembly brick', @tmim mim, @str expression[, @int region[, @int is_symmetric[, @int is_coercive]]])
      Deprecated. Use MODEL:SET('add linear term') instead. @*/
    sub_command
      ("add linear generic assembly brick", 2, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string expr = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       int is_symmetric = 0;
       if (in.remaining()) is_symmetric = in.pop().to_integer();
       int is_coercive = 0;
       if (in.remaining()) is_coercive = in.pop().to_integer();
       
       size_type ind
       = getfem::add_linear_term
       (*md, *mim, expr, region, is_symmetric,
	is_coercive) + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add nonlinear term', @tmim mim, @str expression[, @int region[, @int is_symmetric[, @int is_coercive]]])
      Adds a nonlinear term given by the assembly string `expr` which will
      be assembled in region `region` and with the integration method `mim`.
      The expression can describe a potential or a weak form. Second order
      terms (i.e. containing second order test functions, Test2) are not
      allowed.
      You can specify if the term is symmetric, coercive or not.
      If you are not sure, the better is to declare the term not symmetric
      and not coercive. But some solvers (conjugate gradient for instance)
      are not allowed for non-coercive problems.
      `brickname` is an optional name for the brick.@*/
    sub_command
      ("add nonlinear term", 2, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string expr = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       int is_symmetric = 0;
       if (in.remaining()) is_symmetric = in.pop().to_integer();
       int is_coercive = 0;
       if (in.remaining()) is_coercive = in.pop().to_integer();
       
       size_type ind
       = getfem::add_nonlinear_term
       (*md, *mim, expr, region, is_symmetric,
        is_coercive) + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add nonlinear twodomain term', @tmim mim, @str expression, @int region, @str secondary_domain[, @int is_symmetric[, @int is_coercive]])
      Adds a nonlinear term given by a weak form language expression like
      MODEL:SET('add nonlinear term') but for an integration on a direct
      product of two domains, a first specfied by ``mim`` and ``region``
      and a second one by ``secondary_domain`` which has to be declared
      first into the model.@*/
    sub_command
      ("add nonlinear twodomain term", 4, 6, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string expr = in.pop().to_string();
       size_type region = in.pop().to_integer();
       std::string secdom = in.pop().to_string();
       int is_symmetric = 0;
       if (in.remaining()) is_symmetric = in.pop().to_integer();
       int is_coercive = 0;
       if (in.remaining()) is_coercive = in.pop().to_integer();
       
       size_type ind
       = getfem::add_nonlinear_twodomain_term
       (*md, *mim, expr, region, secdom, is_symmetric,
        is_coercive) + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );
    

    /*@SET ind = ('add nonlinear generic assembly brick', @tmim mim, @str expression[, @int region[, @int is_symmetric[, @int is_coercive]]])
      Deprecated. Use MODEL:SET('add nonlinear term') instead.@*/
    sub_command
      ("add nonlinear generic assembly brick", 2, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string expr = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       int is_symmetric = 0;
       if (in.remaining()) is_symmetric = in.pop().to_integer();
       int is_coercive = 0;
       if (in.remaining()) is_coercive = in.pop().to_integer();
       
       size_type ind
       = getfem::add_nonlinear_term
       (*md, *mim, expr, region, is_symmetric,
        is_coercive) + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );
    
    /*@SET ind = ('add source term', @tmim mim, @str expression[, @int region])
      Adds a source term given by the assembly string `expr` which will
      be assembled in region `region` and with the integration method `mim`.
      Only the residual term will be taken into account.
      Take care that if the expression contains some variables and if the
      expression is a potential, the expression will be
      derivated with respect to all variables.
      `brickname` is an optional name for the brick.@*/
    sub_command
      ("add source term", 2, 3, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string expr = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       
       size_type ind
       = getfem::add_source_term
       (*md, *mim, expr, region) + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add twodomain source term', @tmim mim, @str expression, @int region, @str secondary_domain)
      Adds a source term given by a weak form language expression like
      MODEL:SET('add source term') but for an integration on a direct
      product of two domains, a first specfied by ``mim`` and ``region``
      and a second one by ``secondary_domain`` which has to be declared
      first into the model.@*/
    sub_command
      ("add twodomain source term", 4, 4, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string expr = in.pop().to_string();
       size_type region = in.pop().to_integer();
       std::string secdom = in.pop().to_string();
       
       size_type ind
       = getfem::add_twodomain_source_term
       (*md, *mim, expr, region, secdom) + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add source term generic assembly brick', @tmim mim, @str expression[, @int region])
      Deprecated. Use MODEL:SET('add source term') instead. @*/
    sub_command
      ("add source term generic assembly brick", 2, 3, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string expr = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       
       size_type ind
       = getfem::add_source_term_generic_assembly_brick
       (*md, *mim, expr, region) + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ('add assembly assignment', @str dataname, @str expression[, @int region[, @int order[, @int before]]]) 	 
      Adds expression `expr` to be evaluated at assembly time and being	 
      assigned to the data `dataname` which has to be of im_data type.
      This allows for instance to store a sub-expression of an assembly
      computation to be used on an other assembly. It can be used for instance
      to store the plastic strain in plasticity models.
      `order` represents the order of assembly where this assignement has to be
      done (potential(0), weak form(1) or tangent system(2) or at each
      order(-1)). The default value is 1.
      If before = 1, the the assignement is perfromed before the computation
      of the other assembly terms, such that the data can be used in the
      remaining of the assembly as an intermediary result (be careful that it is
      still considered as a data, no derivation of the expression is performed
      for the tangent system). 	 
      If before = 0 (default), the assignement is done after the assembly terms.
      @*/ 	 
    sub_command 	 
      ("add assembly assignment", 2, 5, 0, 0, 	 
       std::string dataname = in.pop().to_string(); 	 
       std::string expr = in.pop().to_string(); 	 
       size_type region = size_type(-1); 	 
       if (in.remaining()) region = in.pop().to_integer(); 	 
       size_type order = 1; 	 
       if (in.remaining()) order = in.pop().to_integer(); 	 
       bool before = false; 	 
       if (in.remaining()) before = (in.pop().to_integer() != 0); 	 
       
       md->add_assembly_assignments(dataname, expr, region, order, before);
       ); 	 
    
    /*@SET ('clear assembly assignment') 	 
      Delete all added assembly assignments 	 
      @*/ 	 
    sub_command 	 
      ("clear assembly assignment", 0, 0, 0, 0, 	 
       md->clear_assembly_assignments(); 	 
       );
    
    /*@SET ind = ('add Laplacian brick', @tmim mim, @str varname[, @int region])
    Add a Laplacian term to the model relatively to the variable `varname`
    (in fact with a minus : :math:`-\text{div}(\nabla u)`).
    If this is a vector valued variable, the Laplacian term is added
    componentwise. `region` is an optional mesh region on which the term
    is added. If it is not specified, it is added on the whole mesh. Return
    the brick index in the model.@*/
    sub_command
      ("add Laplacian brick", 2, 3, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind
       = getfem::add_Laplacian_brick(*md, *mim, varname, region)
       + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add generic elliptic brick', @tmim mim, @str varname, @str dataname[, @int region])
    Add a generic elliptic term to the model relatively to the variable `varname`.
    The shape of the elliptic term depends both on the variable and the data.
    This corresponds to a term
    :math:`-\text{div}(a\nabla u)`
    where :math:`a` is the data and :math:`u` the variable. The data can be
    a scalar,
    a matrix or an order four tensor. The variable can be vector valued or
    not. If the data is a scalar or a matrix and the variable is vector
    valued then the term is added componentwise. An order four tensor data
    is allowed for vector valued variable only. The data can be constant or
    describbed on a fem. Of course, when the data is a tensor describe on a
    finite element method (a tensor field) the data can be a huge vector.
    The components of the matrix/tensor have to be stored with the fortran
    order (columnwise) in the data vector (compatibility with blas). The
    symmetry of the given matrix/tensor is not verified (but assumed). If
    this is a vector valued variable, the elliptic term is added
    componentwise. `region` is an optional mesh region on which the term is
    added. If it is not specified, it is added on the whole mesh. Note that
    for the real
    version which uses the high-level generic assembly language, `dataname`
    can be any regular expression of the high-level generic assembly
    language (like "1", "sin(X(1))" or "Norm(u)" for instance) even
    depending on model variables. Return the
    brick index in the model.@*/
    sub_command
      ("add generic elliptic brick", 3, 4, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string dataname = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind
       = getfem::add_generic_elliptic_brick(*md, *mim,
                                            varname, dataname, region)
       + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add source term brick', @tmim mim, @str varname, @str dataexpr[, @int region[, @str directdataname]])
    Add a source term to the model relatively to the variable `varname`.
    The source term is
    represented by `dataexpr` which could be any regular expression of the
    high-level generic assembly language (except for the complex version
    where it has to be a declared data of the model).
    `region` is an optional mesh region
    on which the term is added. An additional optional data `directdataname`
    can be provided. The corresponding data vector will be directly added
    to the right hand side without assembly. Note that when region is a
    boundary, this brick allows to prescribe a nonzero Neumann boundary
    condition. Return the brick index in the model.@*/
    sub_command
      ("add source term brick", 3, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string dataname = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       std::string directdataname;
       if (in.remaining()) directdataname = in.pop().to_string();
       size_type ind
       = getfem::add_source_term_brick(*md, *mim,
                                 varname, dataname, region, directdataname)
       + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add normal source term brick', @tmim mim, @str varname, @str dataname, @int region)
      Add a source term on the variable `varname` on a boundary `region`.
      This region should be a boundary. The source term is
      represented by the data `dataepxpr` which could be any regular
      expression of the high-level generic assembly language (except
      for the complex version where it has to be a declared data of
      the model). A scalar
      product with the outward normal unit vector to the boundary is performed.
      The main aim of this brick is to represent a Neumann condition with a
      vector data without performing the scalar product with the normal as a
      pre-processing. Return the brick index in the model.@*/
    sub_command
      ("add normal source term brick", 4, 4, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string dataname = in.pop().to_string();
       size_type region = in.pop().to_integer();
       size_type ind
       = getfem::add_normal_source_term_brick(*md, *mim,
                                              varname, dataname, region)
       + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add Dirichlet condition with simplification', @str varname, @int region[, @str dataname])
      Adds a (simple) Dirichlet condition on the variable `varname` and
      the mesh region `region`. The Dirichlet condition is prescribed by
      a simple post-treatment of the final linear system (tangent system
      for nonlinear problems) consisting of modifying the lines corresponding
      to the degree of freedom of the variable on `region` (0 outside the
      diagonal, 1 on the diagonal of the matrix and the expected value on
      the right hand side).
      The symmetry of the linear system is kept if all other bricks are
      symmetric.
      This brick is to be reserved for simple Dirichlet conditions (only dof
      declared on the corresponding boundary are prescribed). The application
      of this brick on reduced dof may be problematic. Intrinsic vectorial
      finite element method are not supported. 
      `dataname` is the optional right hand side of  the Dirichlet condition.
      It could be constant (but in that case, it can only be applied to
      Lagrange f.e.m.) or (important) described on the same finite
      element method as `varname`.
      Returns the brick index in the model. @*/
    sub_command
      ("add Dirichlet condition with simplification", 2, 3, 0, 1,
       std::string varname = in.pop().to_string();
       size_type region = in.pop().to_integer();
       std::string dataname;
       if (in.remaining()) dataname = in.pop().to_string();

       size_type ind = config::base_index();
       ind += getfem::add_Dirichlet_condition_with_simplification
           (*md, varname, region, dataname);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add Dirichlet condition with multipliers', @tmim mim, @str varname, mult_description, @int region[, @str dataname])
      Add a Dirichlet condition on the variable `varname` and the mesh
      region `region`. This region should be a boundary. The Dirichlet
      condition is prescribed with a multiplier variable described by
      `mult_description`. If `mult_description` is a string this is assumed
      to be the variable name corresponding to the multiplier (which should be
      first declared as a multiplier variable on the mesh region in the model).
      If it is a finite element method (mesh_fem object) then a multiplier
      variable will be added to the model and build on this finite element
      method (it will be restricted to the mesh region `region` and eventually
      some conflicting dofs with some other multiplier variables will be
      suppressed). If it is an integer, then a  multiplier variable will be
      added to the model and build on a classical finite element of degree
      that integer. `dataname` is the optional right hand side of  the
      Dirichlet condition. It could be constant or described on a fem; scalar
      or vector valued, depending on the variable on which the Dirichlet
      condition is prescribed. Return the brick index in the model.@*/
    sub_command
      ("add Dirichlet condition with multipliers", 4, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       int version = 0;
       size_type degree = 0;
       std::string multname;
       getfem::mesh_fem *mf = 0;
       mexarg_in argin = in.pop();
       if (argin.is_integer()) {
         degree = argin.to_integer();
         version = 1;
       } else if (argin.is_string()) {
         multname = argin.to_string();
         version = 2;
       } else {
         mf = to_meshfem_object(argin);
         version = 3;
       }
       size_type region = in.pop().to_integer();
       std::string dataname;
       if (in.remaining()) dataname = in.pop().to_string();

       size_type ind = config::base_index();
       switch(version) {
       case 1:  ind += getfem::add_Dirichlet_condition_with_multipliers
           (*md, *mim, varname, dim_type(degree), region, dataname);
         break;
       case 2:  ind += getfem::add_Dirichlet_condition_with_multipliers
           (*md, *mim, varname, multname, region, dataname);
         break;
       case 3:  ind += getfem::add_Dirichlet_condition_with_multipliers
           (*md, *mim, varname, *mf, region, dataname);
         workspace().set_dependence(md, mf);
         break;
       }
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add Dirichlet condition with Nitsche method', @tmim mim, @str varname, @str Neumannterm, @str datagamma0, @int region[, @scalar theta][, @str dataname])
      Add a Dirichlet condition on the variable `varname` and the mesh
      region `region`. This region should be a boundary. `Neumannterm`
      is the expression of the Neumann term (obtained by the Green formula)
      described as an expression of the high-level
      generic assembly language. This term can be obtained by 
      MODEL:GET('Neumann term', varname, region) once all volumic bricks have
      been added to the model. The Dirichlet
      condition is prescribed with Nitsche's method. `datag` is the optional
      right hand side of the Dirichlet condition. `datagamma0` is the
      Nitsche's method parameter. `theta` is a scalar value which can be
      positive or negative. `theta = 1` corresponds to the standard symmetric
      method which is conditionally coercive for  `gamma0` small.
      `theta = -1` corresponds to the skew-symmetric method which is
      inconditionally coercive. `theta = 0` (default) is the simplest method
      for which the second derivative of the Neumann term is not necessary
      even for nonlinear problems. Return the brick index in the model.
    @*/
    sub_command
      ("add Dirichlet condition with Nitsche method", 5, 7, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string Neumannterm = in.pop().to_string();
       std::string gamma0name = in.pop().to_string();
       size_type region = in.pop().to_integer();
       scalar_type theta = scalar_type(0);
       std::string dataname;
       if (in.remaining()) {
         mexarg_in argin = in.pop();
         if (argin.is_string())
           dataname = argin.to_string();
         else
           theta = argin.to_scalar();
       }
       if (in.remaining()) dataname = in.pop().to_string();

       size_type ind = config::base_index();
       ind += getfem::add_Dirichlet_condition_with_Nitsche_method
       (*md, *mim, varname, Neumannterm,
        gamma0name, region, theta, dataname);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add Dirichlet condition with penalization', @tmim mim, @str varname, @scalar coeff, @int region[, @str dataname, @tmf mf_mult])
    Add a Dirichlet condition on the variable `varname` and the mesh
    region `region`. This region should be a boundary. The Dirichlet
    condition is prescribed with penalization. The penalization coefficient
    is initially `coeff` and will be added to the data of the model.
    `dataname` is the optional right hand side of the Dirichlet condition.
    It could be constant or described on a fem; scalar or vector valued,
    depending on the variable on which the Dirichlet condition is prescribed.
    `mf_mult` is an optional parameter which allows to weaken the
    Dirichlet condition specifying a multiplier space.
    Return the brick index in the model.@*/
    sub_command
      ("add Dirichlet condition with penalization", 4, 6, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       double coeff = in.pop().to_scalar();
       size_type region = in.pop().to_integer();
       std::string dataname;
       if (in.remaining()) dataname = in.pop().to_string();
       const getfem::mesh_fem *mf_mult = 0;
       if (in.remaining()) mf_mult = to_meshfem_object(in.pop());
       size_type ind = config::base_index();
       ind += getfem::add_Dirichlet_condition_with_penalization
       (*md, *mim, varname, coeff, region,
        dataname, mf_mult);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add normal Dirichlet condition with multipliers', @tmim mim, @str varname, mult_description, @int region[, @str dataname])
      Add a Dirichlet condition to the normal component of the vector
     (or tensor) valued variable `varname` and the mesh
      region `region`. This region should be a boundary. The Dirichlet
      condition is prescribed with a multiplier variable described by
      `mult_description`. If `mult_description` is a string this is assumed
      to be the variable name corresponding to the multiplier (which should be
      first declared as a multiplier variable on the mesh region in the model).
      If it is a finite element method (mesh_fem object) then a multiplier
      variable will be added to the model and build on this finite element
      method (it will be restricted to the mesh region `region` and eventually
      some conflicting dofs with some other multiplier variables will be
      suppressed). If it is an integer, then a  multiplier variable will be
      added to the model and build on a classical finite element of degree
      that integer. `dataname` is the optional right hand side of  the
      Dirichlet condition. It could be constant or described on a fem; scalar
      or vector valued, depending on the variable on which the Dirichlet
      condition is prescribed (scalar if the variable
      is vector valued, vector if the variable is tensor valued).
      Returns the brick index in the model.@*/
    sub_command
      ("add normal Dirichlet condition with multipliers", 4, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       int version = 0;
       size_type degree = 0;
       std::string multname;
       getfem::mesh_fem *mf = 0;
       mexarg_in argin = in.pop();
       if (argin.is_integer()) {
         degree = argin.to_integer();
         version = 1;
       } else if (argin.is_string()) {
         multname = argin.to_string();
         version = 2;
       } else {
         mf = to_meshfem_object(argin);
         version = 3;
       }
       size_type region = in.pop().to_integer();
       std::string dataname;
       if (in.remaining()) dataname = in.pop().to_string();

       size_type ind = config::base_index();
       switch(version) {
       case 1:  ind += getfem::add_normal_Dirichlet_condition_with_multipliers
           (*md, *mim, varname, dim_type(degree), region, dataname);
         break;
       case 2:  ind += getfem::add_normal_Dirichlet_condition_with_multipliers
           (*md, *mim, varname, multname, region, dataname);
         break;
       case 3:  ind += getfem::add_normal_Dirichlet_condition_with_multipliers
           (*md, *mim, varname, *mf, region, dataname);
         workspace().set_dependence(md, mf);
         break;
       }
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add normal Dirichlet condition with penalization', @tmim mim, @str varname, @scalar coeff, @int region[, @str dataname, @tmf mf_mult])
    Add a Dirichlet condition to the normal component of the vector
    (or tensor) valued variable `varname` and the mesh
    region `region`. This region should be a boundary. The Dirichlet
    condition is prescribed with penalization. The penalization coefficient
    is initially `coeff` and will be added to the data of the model.
    `dataname` is the optional right hand side of the Dirichlet condition.
    It could be constant or described on a fem; scalar or vector valued,
    depending on the variable on which the Dirichlet condition is prescribed
    (scalar if the variable
    is vector valued, vector if the variable is tensor valued).
    `mf_mult` is an optional parameter which allows to weaken the
    Dirichlet condition specifying a multiplier space.
    Returns the brick index in the model.@*/
    sub_command
      ("add normal Dirichlet condition with penalization", 4, 6, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       double coeff = in.pop().to_scalar();
       size_type region = in.pop().to_integer();
       std::string dataname;
       if (in.remaining()) dataname = in.pop().to_string();
       const getfem::mesh_fem *mf_mult = 0;
       if (in.remaining()) mf_mult = to_meshfem_object(in.pop());
       size_type ind = config::base_index();
       ind += getfem::add_normal_Dirichlet_condition_with_penalization
       (*md, *mim, varname, coeff, region,
        dataname, mf_mult);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add normal Dirichlet condition with Nitsche method', @tmim mim, @str varname, @str Neumannterm, @str gamma0name, @int region[, @scalar theta][, @str dataname])
      Add a Dirichlet condition to the normal component of the vector
      (or tensor) valued variable `varname` and the mesh region `region`.
      This region should be a boundary. `Neumannterm`
      is the expression of the Neumann term (obtained by the Green formula)
      described as an expression of the high-level
      generic assembly language. This term can be obtained by 
      MODEL:GET('Neumann term', varname, region) once all volumic bricks have
      been added to the model. The Dirichlet
      condition is prescribed with Nitsche's method. `dataname` is the optional
      right hand side of the Dirichlet condition. It could be constant or
      described on a fem. `gamma0name` is the
      Nitsche's method parameter. `theta` is a scalar value which can be
      positive or negative. `theta = 1` corresponds to the standard symmetric
      method which is conditionally coercive for  `gamma0` small.
      `theta = -1` corresponds to the skew-symmetric method which is
      inconditionally coercive. `theta = 0` is the simplest method
      for which the second derivative of the Neumann term is not necessary
      even for nonlinear problems. 
      Returns the brick index in the model.
      (This brick is not fully tested)
    @*/
    sub_command
      ("add normal Dirichlet condition with Nitsche method", 5, 7, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string Neumannterm = in.pop().to_string();
       std::string gamma0name = in.pop().to_string();
       size_type region = in.pop().to_integer();
       scalar_type theta = scalar_type(1);
       std::string dataname;
       if (in.remaining()) {
         mexarg_in argin = in.pop();
         if (argin.is_string())
           dataname = argin.to_string();
         else
           theta = argin.to_scalar();
       }
       if (in.remaining()) dataname = in.pop().to_string();

       size_type ind = config::base_index();
       ind += getfem::add_normal_Dirichlet_condition_with_Nitsche_method
       (*md, *mim, varname, Neumannterm,
        gamma0name, region, theta, dataname);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add generalized Dirichlet condition with multipliers', @tmim mim, @str varname, mult_description, @int region, @str dataname, @str Hname)
    Add a Dirichlet condition on the variable `varname` and the mesh
    region `region`.  This version is for vector field.
    It prescribes a condition :math:`Hu = r`
    where `H` is a matrix field. The region should be a boundary. The Dirichlet
    condition is prescribed with a multiplier variable described by
    `mult_description`. If `mult_description` is a string this is assumed
    to be the variable name corresponding to the multiplier (which should be
    first declared as a multiplier variable on the mesh region in the model).
    If it is a finite element method (mesh_fem object) then a multiplier
    variable will be added to the model and build on this finite element
    method (it will be restricted to the mesh region `region` and eventually
    some conflicting dofs with some other multiplier variables will be
    suppressed). If it is an integer, then a  multiplier variable will be
    added to the model and build on a classical finite element of degree
    that integer. `dataname` is the right hand side of  the
    Dirichlet condition. It could be constant or described on a fem; scalar
    or vector valued, depending on the variable on which the Dirichlet
    condition is prescribed. `Hname` is the data
    corresponding to the matrix field `H`.
    Returns the brick index in the model.@*/
    sub_command
      ("add generalized Dirichlet condition with multipliers", 6, 6, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       int version = 0;
       size_type degree = 0;
       std::string multname;
       getfem::mesh_fem *mf = 0;
       mexarg_in argin = in.pop();
       if (argin.is_integer()) {
         degree = argin.to_integer();
         version = 1;
       } else if (argin.is_string()) {
         multname = argin.to_string();
         version = 2;
       } else {
         mf = to_meshfem_object(argin);
         version = 3;
       }
       size_type region = in.pop().to_integer();
       std::string dataname = in.pop().to_string();
       std::string Hname = in.pop().to_string();
       size_type ind = config::base_index();
       switch(version) {
       case 1:  ind += getfem::add_generalized_Dirichlet_condition_with_multipliers
           (*md, *mim, varname, dim_type(degree), region, dataname, Hname);
         break;
       case 2:  ind += getfem::add_generalized_Dirichlet_condition_with_multipliers
           (*md, *mim, varname, multname, region, dataname, Hname);
         break;
       case 3:  ind += getfem::add_generalized_Dirichlet_condition_with_multipliers
           (*md, *mim, varname, *mf, region, dataname, Hname);
         workspace().set_dependence(md, mf);
         break;
       }
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add generalized Dirichlet condition with penalization', @tmim mim, @str varname, @scalar coeff, @int region, @str dataname, @str Hname[, @tmf mf_mult])
      Add a Dirichlet condition on the variable `varname` and the mesh
      region `region`. This version is for vector field.
      It prescribes a condition :math:`Hu = r`
      where `H` is a matrix field.
      The region should be a boundary. The Dirichlet
      condition is prescribed with penalization. The penalization coefficient
      is intially `coeff` and will be added to the data of the model.
      `dataname` is the right hand side of the Dirichlet condition.
      It could be constant or described on a fem; scalar or vector valued,
      depending on the variable on which the Dirichlet condition is prescribed.
      `Hname` is the data
      corresponding to the matrix field `H`. It has to be a constant matrix
      or described on a scalar fem.
      `mf_mult` is an optional parameter which allows to weaken the
      Dirichlet condition specifying a multiplier space.
      Return the brick index in the model.@*/
    sub_command
      ("add generalized Dirichlet condition with penalization", 6, 7, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       double coeff = in.pop().to_scalar();
       size_type region = in.pop().to_integer();
       std::string dataname= in.pop().to_string();
       std::string Hname= in.pop().to_string();
       const getfem::mesh_fem *mf_mult = 0;
       if (in.remaining()) mf_mult = to_meshfem_object(in.pop());
       size_type ind = config::base_index();
       ind += getfem::add_generalized_Dirichlet_condition_with_penalization
       (*md, *mim, varname, coeff, region,
        dataname, Hname, mf_mult);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add generalized Dirichlet condition with Nitsche method', @tmim mim, @str varname, @str Neumannterm, @str gamma0name, @int region[, @scalar theta], @str dataname, @str Hname)
      Add a Dirichlet condition on the variable `varname` and the mesh
      region `region`.
      This version is for vector field. It prescribes a condition
      @f$ Hu = r @f$ where `H` is a matrix field.
      CAUTION : the matrix H should have all eigenvalues equal to 1 or 0.
      The region should be a boundary.   `Neumannterm`
      is the expression of the Neumann term (obtained by the Green formula)
      described as an expression of the high-level
      generic assembly language. This term can be obtained by 
      MODEL:GET('Neumann term', varname, region) once all volumic bricks have
      been added to the model.  The Dirichlet
      condition is prescribed with Nitsche's method. `dataname` is the optional
      right hand side of the Dirichlet condition. It could be constant or
      described on a fem. `gamma0name` is the
      Nitsche's method parameter. `theta` is a scalar value which can be
      positive or negative. `theta = 1` corresponds to the standard symmetric
      method which is conditionally coercive for  `gamma0` small.
      `theta = -1` corresponds to the skew-symmetric method which is
      inconditionally coercive. `theta = 0` is the simplest method
      for which the second derivative of the Neumann term is not necessary
      even for nonlinear problems. `Hname` is the data
      corresponding to the matrix field `H`. It has to be a constant matrix
      or described on a scalar fem. Returns the brick index in the model.
      (This brick is not fully tested)
    @*/
    sub_command
      ("add generalized Dirichlet condition with Nitsche method", 7, 8, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string Neumannterm = in.pop().to_string();
       std::string gamma0name = in.pop().to_string();
       size_type region = in.pop().to_integer();
       scalar_type theta = scalar_type(1);
       std::string dataname;
       if (in.remaining()) {
         mexarg_in argin = in.pop();
         if (argin.is_string())
           dataname = argin.to_string();
         else
           theta = argin.to_scalar();
       }
       dataname = in.pop().to_string();
       std::string Hname= in.pop().to_string();

       size_type ind = config::base_index();
       ind += getfem::add_generalized_Dirichlet_condition_with_Nitsche_method
       (*md, *mim, varname, Neumannterm,
        gamma0name, region, theta, dataname, Hname);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add pointwise constraints with multipliers', @str varname, @str dataname_pt[, @str dataname_unitv] [, @str dataname_val])
    Add some pointwise constraints on the variable `varname` using
    multiplier. The multiplier variable is automatically added to the model.
    The conditions are prescribed on a set of points given in the data
    `dataname_pt` whose dimension is the number of points times the dimension
    of the mesh.
    If the variable represents a vector field, one has to give the data
    `dataname_unitv` which represents a vector of dimension the number of
    points times the dimension of the vector field which should store some
    unit vectors. In that case the prescribed constraint is the scalar
    product of the variable at the corresponding point with the corresponding
    unit vector.
    The optional data `dataname_val` is the vector of values to be prescribed
    at the different points.
    This brick is specifically designed to kill rigid displacement
    in a Neumann problem.
    Returns the brick index in the model.@*/
    sub_command
      ("add pointwise constraints with multipliers", 2, 4, 0, 1,
       std::string varname = in.pop().to_string();
       std::string dataname_pt = in.pop().to_string();
       const getfem::mesh_fem *mf_u
         = &(md->mesh_fem_of_variable(varname));
       GMM_ASSERT1(mf_u, "The variable should depend on a mesh_fem");
       std::string dataname_unitv;
       if (mf_u->get_qdim() > 1)
         dataname_unitv = in.pop().to_string();
       std::string dataname_val;
       if (in.remaining()) dataname_val = in.pop().to_string();

       size_type ind = config::base_index();
       ind += getfem::add_pointwise_constraints_with_multipliers
       (*md, varname, dataname_pt, dataname_unitv, dataname_val);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add pointwise constraints with given multipliers', @str varname, @str multname, @str dataname_pt[, @str dataname_unitv] [, @str dataname_val])
    Add some pointwise constraints on the variable `varname` using a given
    multiplier `multname`.
    The conditions are prescribed on a set of points given in the data
    `dataname_pt` whose dimension is the number of points times the dimension
    of the mesh.
    The multiplier variable should be a fixed size variable of size the
    number of points.
    If the variable represents a vector field, one has to give the data
    `dataname_unitv` which represents a vector of dimension the number of
    points times the dimension of the vector field which should store some
    unit vectors. In that case the prescribed constraint is the scalar
    product of the variable at the corresponding point with the corresponding
    unit vector.
    The optional data `dataname_val` is the vector of values to be prescribed
    at the different points.
    This brick is specifically designed to kill rigid displacement
    in a Neumann problem.
    Returns the brick index in the model.@*/
    sub_command
      ("add pointwise constraints with given multipliers", 3, 5, 0, 1,
       std::string varname = in.pop().to_string();
       std::string multname = in.pop().to_string();
       std::string dataname_pt = in.pop().to_string();
       const getfem::mesh_fem *mf_u
         = &(md->mesh_fem_of_variable(varname));
       GMM_ASSERT1(mf_u, "The variable should depend on a mesh_fem");
       std::string dataname_unitv;
       if (mf_u->get_qdim() > 1)
         dataname_unitv = in.pop().to_string();
       std::string dataname_val;
       if (in.remaining()) dataname_val = in.pop().to_string();

       size_type ind = config::base_index();
       ind += getfem::add_pointwise_constraints_with_given_multipliers
       (*md, varname, multname, dataname_pt, dataname_unitv,
        dataname_val);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add pointwise constraints with penalization', @str varname, @scalar coeff, @str dataname_pt[, @str dataname_unitv] [, @str dataname_val])
    Add some pointwise constraints on the variable `varname` thanks to
    a penalization. The penalization coefficient is initially
    `penalization_coeff` and will be added to the data of the model.
    The conditions are prescribed on a set of points given in the data
    `dataname_pt` whose dimension is the number of points times the dimension
    of the mesh.
    If the variable represents a vector field, one has to give the data
    `dataname_unitv` which represents a vector of dimension the number of
    points times the dimension of the vector field which should store some
    unit vectors. In that case the prescribed constraint is the scalar
    product of the variable at the corresponding point with the corresponding
    unit vector.
    The optional data `dataname_val` is the vector of values to be prescribed
    at the different points.
    This brick is specifically designed to kill rigid displacement
    in a Neumann problem.
    Returns the brick index in the model.@*/
    sub_command
      ("add pointwise constraints with penalization", 3, 5, 0, 1,
       std::string varname = in.pop().to_string();
       double coeff = in.pop().to_scalar();
       std::string dataname_pt = in.pop().to_string();
       const getfem::mesh_fem *mf_u
         = &(md->mesh_fem_of_variable(varname));
       GMM_ASSERT1(mf_u, "The variable should depend on a mesh_fem");
       std::string dataname_unitv;
       if (mf_u->get_qdim() > 1)
         dataname_unitv = in.pop().to_string();
       std::string dataname_val;
       if (in.remaining()) dataname_val = in.pop().to_string();

       size_type ind = config::base_index();
       ind += getfem::add_pointwise_constraints_with_penalization
       (*md, varname, coeff, dataname_pt, dataname_unitv,
        dataname_val);
       out.pop().from_integer(int(ind));
       );


    /*@SET ('change penalization coeff', @int ind_brick, @scalar coeff)
    Change the penalization coefficient of a Dirichlet condition with
    penalization brick. If the brick is not of this kind, this
    function has an undefined behavior.@*/
    sub_command
      ("change penalization coeff", 2, 2, 0, 0,
       size_type ind_brick = in.pop().to_integer() - config::base_index();
       double coeff = in.pop().to_scalar();
       getfem::change_penalization_coeff(*md, ind_brick, coeff);
       );


    /*@SET ind = ('add Helmholtz brick', @tmim mim, @str varname, @str dataexpr[, @int region])
      Add a Helmholtz term to the model relatively to the variable `varname`.
      `dataexpr` is the wave number. `region` is an optional mesh
      region on which the term is added. If it is not specified, it is added
      on the whole mesh. Return the brick index in the model.@*/
    sub_command
      ("add Helmholtz brick", 3, 4, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string dataname = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind
       = getfem::add_Helmholtz_brick(*md, *mim,
                                     varname, dataname, region)
       + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add Fourier Robin brick', @tmim mim, @str varname, @str dataexpr, @int region)
    Add a Fourier-Robin term to the model relatively to the variable
    `varname`. This corresponds to a weak term of the form
    :math:`\int (qu).v`. `dataexpr` is the parameter :math:`q` of
    the Fourier-Robin condition.  It can be an arbitrary valid expression
    of the high-level generic assembly language (except for the complex version
    for which it should be a data of the model). `region` is the mesh region
    on which the term is added. Return the brick index in the model.@*/
    sub_command
      ("add Fourier Robin brick", 4, 4, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string dataname = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind
       = getfem::add_Fourier_Robin_brick(*md, *mim,
                                         varname,dataname, region)
       + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add constraint with multipliers', @str varname, @str multname, @tspmat B, {@vec L | @str dataname})
    Add an additional explicit constraint on the variable `varname` thank to
    a multiplier `multname` peviously added to the model (should be a fixed
    size variable). The constraint is :math:`BU=L`
    with `B` being a rectangular sparse matrix. It is possible to change
    the constraint at any time with the methods MODEL:SET('set private matrix')
    and MODEL:SET('set private rhs'). If `dataname` is specified instead of `L`,
    the vector `L` is defined in the model as data with the given name.
    Return the brick index in the model.@*/
    sub_command
      ("add constraint with multipliers", 4, 4, 0, 1,
       std::string varname = in.pop().to_string();
       std::string multname = in.pop().to_string();
       std::shared_ptr<gsparse> B = in.pop().to_sparse();
       if (B->is_complex() && !md->is_complex())
         THROW_BADARG("Complex constraint for a real model");
       if (!B->is_complex() && md->is_complex())
         THROW_BADARG("Real constraint for a complex model");

       size_type ind
       = getfem::add_constraint_with_multipliers(*md,varname,multname);

       if (md->is_complex()) {
         if (B->storage()==gsparse::CSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->cplx_csc());
         else if (B->storage()==gsparse::WSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->cplx_wsc());
         else
           THROW_BADARG("Constraint matrix should be a sparse matrix");
       } else {
         if (B->storage()==gsparse::CSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->real_csc());
         else if (B->storage()==gsparse::WSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->real_wsc());
         else
           THROW_BADARG("Constraint matrix should be a sparse matrix");
       }

       if (in.front().is_string()) {
         std::string dataname = in.pop().to_string();
         getfem::set_private_data_rhs(*md, ind, dataname);
       } else if (!md->is_complex()) {
         darray st = in.pop().to_darray();
         std::vector<double> V(st.begin(), st.end());
         getfem::set_private_data_rhs(*md, ind, V);
       } else {
         carray st = in.pop().to_carray();
         std::vector<std::complex<double> > V(st.begin(), st.end());
         getfem::set_private_data_rhs(*md, ind, V);
       }

       out.pop().from_integer(int(ind + config::base_index()));
       );


    /*@SET ind = ('add constraint with penalization', @str varname, @scalar coeff, @tspmat B, {@vec L | @str dataname})
    Add an additional explicit penalized constraint on the variable `varname`.
    The constraint is :math`BU=L` with `B` being a rectangular sparse matrix.
    Be aware that `B` should not contain a plain row, otherwise the whole
    tangent matrix will be plain. It is possible to change the constraint
    at any time with the methods MODEL:SET('set private matrix')
    and MODEL:SET('set private rhs'). The method
    MODEL:SET('change penalization coeff') can be used.
    If `dataname` is specified instead of `L`, the vector `L` is defined
    in the model as data with the given name.
    Return the brick
    index in the model.@*/
    sub_command
      ("add constraint with penalization", 4, 4, 0, 1,
       std::string varname = in.pop().to_string();
       double coeff = in.pop().to_scalar();
       std::shared_ptr<gsparse> B = in.pop().to_sparse();
       if (B->is_complex() && !md->is_complex())
         THROW_BADARG("Complex constraint for a real model");
       if (!B->is_complex() && md->is_complex())
         THROW_BADARG("Real constraint for a complex model");

       size_type ind
       = getfem::add_constraint_with_penalization(*md, varname, coeff);

       if (md->is_complex()) {
         if (B->storage()==gsparse::CSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->cplx_csc());
         else if (B->storage()==gsparse::WSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->cplx_wsc());
         else
           THROW_BADARG("Constraint matrix should be a sparse matrix");
       } else {
         if (B->storage()==gsparse::CSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->real_csc());
         else if (B->storage()==gsparse::WSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->real_wsc());
         else
           THROW_BADARG("Constraint matrix should be a sparse matrix");
       }

       if (in.front().is_string()) {
         std::string dataname = in.pop().to_string();
         getfem::set_private_data_rhs(*md, ind, dataname);
       } else if (!md->is_complex()) {
         darray st = in.pop().to_darray();
         std::vector<double> V(st.begin(), st.end());
         getfem::set_private_data_rhs(*md, ind, V);
       } else {
         carray st = in.pop().to_carray();
         std::vector<std::complex<double> > V(st.begin(), st.end());
         getfem::set_private_data_rhs(*md, ind, V);
       }

       out.pop().from_integer(int(ind + config::base_index()));
       );


    /*@SET ind = ('add explicit matrix', @str varname1, @str varname2, @tspmat B[, @int issymmetric[, @int iscoercive]])
    Add a brick representing an explicit matrix to be added to the tangent
    linear system relatively to the variables `varname1` and `varname2`.
    The given matrix should have has many rows as the dimension of
    `varname1` and as many columns as the dimension of `varname2`.
    If the two variables are different and if `issymmetric` is set to 1
    then the transpose of the matrix is also added to the tangent system
    (default is 0). Set `iscoercive` to 1 if the term does not affect the
    coercivity of the tangent system (default is 0). The matrix can be
    changed by the command MODEL:SET('set private matrix'). Return the
    brick index in the model.@*/
    sub_command
      ("add explicit matrix", 3, 5, 0, 1,
       std::string varname1 = in.pop().to_string();
       std::string varname2 = in.pop().to_string();
       std::shared_ptr<gsparse> B = in.pop().to_sparse();
       bool issymmetric = false;
       bool iscoercive = false;
       if (in.remaining()) issymmetric = (in.pop().to_integer(0,1) != 0);
       if (!issymmetric && in.remaining())
         iscoercive = (in.pop().to_integer(0,1) != 0);

       size_type ind
       = getfem::add_explicit_matrix(*md, varname1, varname2,
                                     issymmetric, iscoercive);

       if (B->is_complex() && !md->is_complex())
         THROW_BADARG("Complex constraint for a real model");
       if (!B->is_complex() && md->is_complex())
         THROW_BADARG("Real constraint for a complex model");

       if (md->is_complex()) {
         if (B->storage()==gsparse::CSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->cplx_csc());
         else if (B->storage()==gsparse::WSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->cplx_wsc());
         else
           THROW_BADARG("Constraint matrix should be a sparse matrix");
       } else {
         if (B->storage()==gsparse::CSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->real_csc());
         else if (B->storage()==gsparse::WSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->real_wsc());
         else
           THROW_BADARG("Constraint matrix should be a sparse matrix");
       }

       out.pop().from_integer(int(ind + config::base_index()));
       );


    /*@SET ind = ('add explicit rhs', @str varname, @vec L)
      Add a brick representing an explicit right hand side to be added to
      the right hand side of the tangent linear system relatively to the
      variable `varname`. The given rhs should have the same size than the
      dimension of `varname`. The rhs can be changed by the command
      MODEL:SET('set private rhs'). If `dataname` is specified instead of
      `L`, the vector `L` is defined in the model as data with the given name.
      Return the brick index in the model.@*/
    sub_command
      ("add explicit rhs", 2, 2, 0, 1,
       std::string varname = in.pop().to_string();
       size_type ind
       = getfem::add_explicit_rhs(*md, varname);

       if (in.front().is_string()) {
         std::string dataname = in.pop().to_string();
         getfem::set_private_data_rhs(*md, ind, dataname);
       } else if (!md->is_complex()) {
         darray st = in.pop().to_darray();
         std::vector<double> V(st.begin(), st.end());
         getfem::set_private_data_rhs(*md, ind, V);
       } else {
         carray st = in.pop().to_carray();
         std::vector<std::complex<double> > V(st.begin(), st.end());
         getfem::set_private_data_rhs(*md, ind, V);
       }

       out.pop().from_integer(int(ind + config::base_index()));
       );


    /*@SET ('set private matrix', @int indbrick, @tspmat B)
      For some specific bricks having an internal sparse matrix
      (explicit bricks: 'constraint brick' and 'explicit matrix brick'),
      set this matrix. @*/
    sub_command
      ("set private matrix", 2, 2, 0, 0,
       size_type ind = in.pop().to_integer() - config::base_index();
       std::shared_ptr<gsparse> B = in.pop().to_sparse();

       if (B->is_complex() && !md->is_complex())
         THROW_BADARG("Complex constraint for a real model");
       if (!B->is_complex() && md->is_complex())
         THROW_BADARG("Real constraint for a complex model");

       if (md->is_complex()) {
         if (B->storage()==gsparse::CSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->cplx_csc());
         else if (B->storage()==gsparse::WSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->cplx_wsc());
         else
           THROW_BADARG("Constraint matrix should be a sparse matrix");
       } else {
         if (B->storage()==gsparse::CSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->real_csc());
         else if (B->storage()==gsparse::WSCMAT)
           getfem::set_private_data_matrix(*md, ind, B->real_wsc());
         else
           THROW_BADARG("Constraint matrix should be a sparse matrix");
       }
       );


    /*@SET ('set private rhs', @int indbrick, @vec B)
      For some specific bricks having an internal right hand side vector
      (explicit bricks: 'constraint brick' and 'explicit rhs brick'),
      set this rhs. @*/
    sub_command
      ("set private rhs", 2, 2, 0, 0,
       size_type ind = in.pop().to_integer() - config::base_index();

       if (!md->is_complex()) {
         darray st = in.pop().to_darray();
         std::vector<double> V(st.begin(), st.end());
         getfem::set_private_data_rhs(*md, ind, V);
       } else {
         carray st = in.pop().to_carray();
         std::vector<std::complex<double> > V(st.begin(), st.end());
         getfem::set_private_data_rhs(*md, ind, V);
       }
       );


    /*@SET ind = ('add isotropic linearized elasticity brick', @tmim mim, @str varname, @str dataname_lambda, @str dataname_mu[, @int region])
      Add an isotropic linearized elasticity term to the model relatively to
      the variable `varname`. `dataname_lambda` and `dataname_mu` should
      contain the Lame coefficients. `region` is an optional mesh region
      on which the term is added. If it is not specified, it is added
      on the whole mesh. Return the brick index in the model.@*/
    sub_command
      ("add isotropic linearized elasticity brick", 4, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string dataname_lambda = in.pop().to_string();
       std::string dataname_mu = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind
       = getfem::add_isotropic_linearized_elasticity_brick
       (*md, *mim, varname, dataname_lambda, dataname_mu, region)
       + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add isotropic linearized elasticity brick pstrain', @tmim mim, @str varname, @str data_E, @str data_nu[, @int region])
      Add an isotropic linearized elasticity term to the model relatively to
      the variable `varname`. `data_E` and `data_nu` should
      contain the Young modulus and Poisson ratio, respectively.
      `region` is an optional mesh region on which the term is added.
      If it is not specified, it is added
      on the whole mesh.
      On two-dimensional meshes, the term will correpsond to a plain strain
      approximation. On three-dimensional meshes, it will correspond to the
      standard model. 
      Return the brick index in the model.@*/
    sub_command
      ("add isotropic linearized elasticity brick pstrain", 4, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string data_E = in.pop().to_string();
       std::string data_nu = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind
       = getfem::add_isotropic_linearized_elasticity_brick_pstrain
       (*md, *mim, varname, data_E, data_nu, region)
       + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add isotropic linearized elasticity brick pstress', @tmim mim, @str varname, @str data_E, @str data_nu[, @int region])
      Add an isotropic linearized elasticity term to the model relatively to
      the variable `varname`. `data_E` and `data_nu` should
      contain the Young modulus and Poisson ratio, respectively.
      `region` is an optional mesh region on which the term is added.
      If it is not specified, it is added
      on the whole mesh.
      On two-dimensional meshes, the term will correpsond to a plain stress
      approximation. On three-dimensional meshes, it will correspond to the
      standard model. 
      Return the brick index in the model.@*/
    sub_command
      ("add isotropic linearized elasticity brick pstress", 4, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string data_E = in.pop().to_string();
       std::string data_nu = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind
       = getfem::add_isotropic_linearized_elasticity_brick_pstress
       (*md, *mim, varname, data_E, data_nu, region)
       + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add linear incompressibility brick', @tmim mim, @str varname, @str multname_pressure[, @int region[, @str dataexpr_coeff]])
    Add a linear incompressibility condition on `variable`. `multname_pressure`
    is a variable which represent the pressure. Be aware that an inf-sup
    condition between the finite element method describing the pressure and the
    primal variable has to be satisfied. `region` is an optional mesh region on
    which the term is added. If it is not specified, it is added on the whole
    mesh. `dataexpr_coeff` is an optional penalization coefficient for nearly
    incompressible elasticity for instance. In this case, it is the inverse
    of the Lame coefficient :math:`\lambda`. Return the brick index in the
    model.@*/
    sub_command
      ("add linear incompressibility brick", 3, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string multname = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       std::string dataname;
       if (in.remaining()) dataname = in.pop().to_string();
       size_type ind
       = getfem::add_linear_incompressibility
       (*md, *mim, varname, multname, region, dataname)
       + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add nonlinear elasticity brick', @tmim mim, @str varname, @str constitutive_law, @str dataname[, @int region])
    Add a nonlinear elasticity term to the model relatively to the
    variable `varname` (deprecated brick, use add_finite_strain_elaticity
    instead). `lawname` is the constitutive law which
    could be 'SaintVenant Kirchhoff', 'Mooney Rivlin', 'neo Hookean',
    'Ciarlet Geymonat' or 'generalized Blatz Ko'.
    'Mooney Rivlin' and 'neo Hookean' law names can be preceded with the word
    'compressible' or 'incompressible' to force using the corresponding version.
    The compressible version of these laws requires one additional material
    coefficient. By default, the incompressible version of 'Mooney Rivlin' law
    and the compressible one of the 'neo Hookean' law are considered. In
    general, 'neo Hookean' is a special case of the 'Mooney Rivlin' law that
    requires one coefficient less.
    IMPORTANT : if the variable is defined on a 2D mesh, the plane strain
    approximation is automatically used.
    `dataname` is a vector of parameters for the constitutive law. Its length
    depends on the law. It could be a short vector of constant values or a
    vector field described on a finite element method for variable
    coefficients. `region` is an optional mesh region on which the term
    is added. If it is not specified, it is added on the whole mesh.
    This brick use the low-level generic assembly.
    Returns the brick index in the model.@*/
    sub_command
      ("add nonlinear elasticity brick", 4, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       size_type N = mim->linked_mesh().dim();
       std::string varname = in.pop().to_string();
       std::string lawname = in.pop().to_string();
       std::string dataname = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind = config::base_index() +
       add_nonlinear_elasticity_brick
       (*md, *mim, varname,
        abstract_hyperelastic_law_from_name(lawname, N), dataname, region);

       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add finite strain elasticity brick', @tmim mim, @str constitutive_law, @str varname, @str params[, @int region])
    Add a nonlinear elasticity term to the model relatively to the
    variable `varname`. `lawname` is the constitutive law which
    could be 'SaintVenant Kirchhoff', 'Mooney Rivlin', 'Neo Hookean',
    'Ciarlet Geymonat' or 'Generalized Blatz Ko'.
    'Mooney Rivlin' and 'Neo Hookean' law names have to be preceeded with
    the word 'Compressible' or 'Incompressible' to force using the
    corresponding version.
    The compressible version of these laws requires one additional material
    coefficient.

    IMPORTANT : if the variable is defined on a 2D mesh, the plane strain
    approximation is automatically used.
    `params` is a vector of parameters for the constitutive law. Its length
    depends on the law. It could be a short vector of constant values or a
    vector field described on a finite element method for variable
    coefficients. `region` is an optional mesh region on which the term
    is added. If it is not specified, it is added on the whole mesh.
    This brick use the high-level generic assembly.
    Returns the brick index in the model.@*/
    sub_command
      ("add finite strain elasticity brick", 4, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       // size_type N = *mim.linked_mesh().dim();
       std::string lawname = in.pop().to_string();
       std::string varname = in.pop().to_string();
       std::string params = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();

       std::string ln = varname; // tolerance for the compatibility with 5.0
       filter_lawname(ln);
       if (ln.compare("saintvenant_kirchhoff") == 0 ||
           ln.compare("saint_venant_kirchhoff") == 0 ||
           ln.compare("generalized_blatz_ko") == 0 ||
           ln.compare("ciarlet_geymonat") == 0 ||
           ln.compare("incompressible_mooney_rivlin") == 0 ||
           ln.compare("compressible_mooney_rivlin") == 0 ||
           ln.compare("incompressible_neo_hookean") == 0 ||
           ln.compare("compressible_neo_hookean") == 0 ||
           ln.compare("compressible_neo_hookean_bonet") == 0 ||
           ln.compare("compressible_neo_hookean_ciarlet") == 0) {
         std::swap(lawname, varname);
       }

       size_type ind = config::base_index() +
       add_finite_strain_elasticity_brick
       (*md, *mim, lawname, varname, params, region);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add small strain elastoplasticity brick', @tmim mim,  @str lawname, @str unknowns_type [, @str varnames, ...] [, @str params, ...] [, @str theta = '1' [, @str dt = 'timestep']] [, @int region = -1])
      Adds a small strain plasticity term to the model `M`. This is the
      main GetFEM brick for small strain plasticity. `lawname` is the name
      of an implemented plastic law, `unknowns_type` indicates the choice
      between a discretization where the plastic multiplier is an unknown of
      the problem or (return mapping approach) just a data of the model
      stored for the next iteration. Remember that in both cases, a multiplier
      is stored anyway. `varnames` is a set of variable and data names with
      length which may depend on the plastic law (at least the displacement,
      the plastic multiplier and the plastic strain). `params` is a list of
      expressions for the parameters (at least elastic coefficients and the
      yield stress). These expressions can be some data names (or even
      variable names) of the model but can also be any scalar valid expression
      of the high level assembly language (such as '1/2', '2+sin(X[0])',
      '1+Norm(v)' ...). The last two parameters optionally provided in
      `params` are the `theta` parameter of the `theta`-scheme (generalized
      trapezoidal rule) used for the plastic strain integration and the
      time-step`dt`. The default value for `theta` if omitted is 1, which
      corresponds to the classical Backward Euler scheme which is first order
      consistent. `theta=1/2` corresponds to the Crank-Nicolson scheme
      (trapezoidal rule) which is second order consistent. Any value
      between 1/2 and 1 should be a valid value. The default value of `dt` is
      'timestep' which simply indicates the time step defined in the model
      (by md.set_time_step(dt)). Alternatively it can be any expression
      (data name, constant value ...). The time step can be altered from one
      iteration to the next one. `region` is a mesh region.

      The available plasticity laws are:

      - 'Prandtl Reuss' (or 'isotropic perfect plasticity').
        Isotropic elasto-plasticity with no hardening. The variables are the
        displacement, the plastic multiplier and the plastic strain.
        The displacement should be a variable and have a corresponding data
        having the same name preceded by 'Previous\_' corresponding to the
        displacement at the previous time step (typically 'u' and 'Previous_u').
        The plastic multiplier should also have two versions (typically 'xi'
        and 'Previous_xi') the first one being defined as data if
        `unknowns_type ` is 'DISPLACEMENT_ONLY' or the integer value 0, or as
        a variable if `unknowns_type` is DISPLACEMENT_AND_PLASTIC_MULTIPLIER
        or the integer value 1.
        The plastic strain should represent a n x n data tensor field stored
        on mesh_fem or (preferably) on an im_data (corresponding to `mim`).
        The data are the first Lame coefficient, the second one (shear modulus)
        and the uniaxial yield stress. A typical call is
        MODEL:GET('add small strain elastoplasticity brick', mim, 'Prandtl Reuss', 0, 'u', 'xi', 'Previous_Ep', 'lambda', 'mu', 'sigma_y', '1', 'timestep');
        IMPORTANT: Note that this law implements
        the 3D expressions. If it is used in 2D, the expressions are just
        transposed to the 2D. For the plane strain approximation, see below.
      - "plane strain Prandtl Reuss"
        (or "plane strain isotropic perfect plasticity")
        The same law as the previous one but adapted to the plane strain
        approximation. Can only be used in 2D.
      - "Prandtl Reuss linear hardening"
        (or "isotropic plasticity linear hardening").
        Isotropic elasto-plasticity with linear isotropic and kinematic
        hardening. An additional variable compared to "Prandtl Reuss" law:
        the accumulated plastic strain. Similarly to the plastic strain, it
        is only stored at the end of the time step, so a simple data is
        required (preferably on an im_data).
        Two additional parameters: the kinematic hardening modulus and the
        isotropic one. 3D expressions only. A typical call is
        MODEL:GET('add small strain elastoplasticity brick', mim, 'Prandtl Reuss linear hardening', 0, 'u', 'xi', 'Previous_Ep', 'Previous_alpha', 'lambda', 'mu', 'sigma_y', 'H_k', H_i', '1', 'timestep');
      - "plane strain Prandtl Reuss linear hardening"
        (or "plane strain isotropic plasticity linear hardening").
        The same law as the previous one but adapted to the plane strain
        approximation. Can only be used in 2D.

      See GetFEM user documentation for further explanations on the
      discretization of the plastic flow and on the implemented plastic laws.
      See also GetFEM user documentation on time integration strategy
      (integration of transient problems).

      IMPORTANT : remember that `small_strain_elastoplasticity_next_iter` has
      to be called at the end of each time step, before the next one
      (and before any post-treatment : this sets the value of the plastic
      strain and plastic multiplier).
      @*/
    sub_command
      ("add small strain elastoplasticity brick", 3, 15, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string lawname = in.pop().to_string();
       filter_lawname(lawname);
       size_type nb_var = 0; size_type nb_params = 0;
       if (lawname.compare("isotropic_perfect_plasticity") == 0 ||
           lawname.compare("prandtl_reuss") == 0 ||
           lawname.compare("plane_strain_isotropic_perfect_plasticity") == 0 ||
           lawname.compare("plane_strain_prandtl_reuss") == 0) {
         nb_var = nb_params = 3;
       } else if
         (lawname.compare("isotropic_plasticity_linear_hardening") == 0 ||
          lawname.compare("prandtl_reuss_linear_hardening") == 0 ||
          lawname.compare("plane_strain_isotropic_plasticity_linear_hardening") == 0 ||
          lawname.compare("plane_strain_prandtl_reuss_linear_hardening") == 0) {
         nb_var = 4; nb_params = 5;
       } else
         THROW_BADARG(lawname << " is not an implemented elastoplastic law");

       getfem::plasticity_unknowns_type unknowns_type(getfem::DISPLACEMENT_ONLY);
       mexarg_in argin = in.pop();
       if (argin.is_string()) {
         std::string opt = argin.to_string();
         filter_lawname(opt);
         if (opt.compare("displacement_only") == 0)
           unknowns_type = getfem::DISPLACEMENT_ONLY;
         else if (opt.compare("displacement_and_plastic_multiplier") == 0)
           unknowns_type = getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER;
         else
           THROW_BADARG("Wrong input");
       } else if (argin.is_integer())
         unknowns_type = static_cast<getfem::plasticity_unknowns_type>
                         (argin.to_integer(0,1));

       std::vector<std::string> varnames;
       for (size_type i = 0; i < nb_var; ++i)
         varnames.push_back(in.pop().to_string());

       std::vector<std::string> params;
       for (size_type i = 0; i < nb_params; ++i)
         params.push_back(in.pop().to_string());

       std::string theta = "1";
       std::string dt = "timestep";
       size_type region = size_type(-1);
       for (size_type i=0; i < 3 && in.remaining(); ++i) {
         argin = in.pop();
         if (argin.is_string()) {
           if (i==0)      theta = argin.to_string();
           else if (i==1) dt = argin.to_string();
           else           THROW_BADARG("Wrong input");
         } else if (argin.is_integer()) {
           region = argin.to_integer();
           GMM_ASSERT1(!in.remaining(), "Wrong input");
         }
       }
       params.push_back(theta);
       params.push_back(dt);

       size_type ind = config::base_index() +
         getfem::add_small_strain_elastoplasticity_brick
         (*md, *mim, lawname, unknowns_type, varnames, params, region);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add elastoplasticity brick', @tmim mim ,@str projname, @str varname, @str previous_dep_name, @str datalambda, @str datamu, @str datathreshold, @str datasigma[, @int region])
      Old (obsolete) brick which do not use the high level generic
      assembly. Add a nonlinear elastoplastic term to the model relatively
      to the variable `varname`, in small deformations, for an isotropic
      material and for a quasistatic model. `projname` is the type of
      projection that used: only the Von Mises projection is
      available with 'VM' or 'Von Mises'.
      `datasigma` is the variable representing the constraints on the material.
      `previous_dep_name` represents the displacement at the previous time step.
      Moreover, the finite element method on which `varname` is described
      is an K ordered mesh_fem, the `datasigma` one have to be at least
      an K-1 ordered mesh_fem.
      `datalambda` and `datamu` are the Lame coefficients of the studied
      material.
      `datathreshold` is the plasticity threshold of the material.
      The three last variables could be constants or described on the
      same finite element method.
      `region` is an optional mesh region on which the term is added.
      If it is not specified, it is added on the whole mesh.
      Return the brick index in the model.@*/
    sub_command
      ("add elastoplasticity brick", 8, 9, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string projname = in.pop().to_string();
       std::string varname = in.pop().to_string();
       std::string previous_dep = in.pop().to_string();
       std::string datalambda = in.pop().to_string();
       std::string datamu = in.pop().to_string();
       std::string datathreshold = in.pop().to_string();
       std::string datasigma = in.pop().to_string();
       size_type region = size_type(-1);

       // getfem::VM_projection proj(0);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind = config::base_index() +
       add_elastoplasticity_brick
       (*md, *mim,
        abstract_constraints_projection_from_name(projname),
        varname, previous_dep, datalambda, datamu,
        datathreshold, datasigma,
        region);

       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add finite strain elastoplasticity brick', @tmim mim , @str lawname, @str unknowns_type [, @str varnames, ...] [, @str params, ...] [, @int region = -1])
      Add a finite strain elastoplasticity brick to the model.
      For the moment there is only one supported law defined through 
      `lawname` as "Simo_Miehe".
      This law supports to possibilities of unknown variables to solve for
      defined by means of `unknowns_type` set to either
      'DISPLACEMENT_AND_PLASTIC_MULTIPLIER' (integer value 1) or
      'DISPLACEMENT_AND_PLASTIC_MULTIPLIER_AND_PRESSURE' (integer value 3).
      The  "Simo_Miehe" law expects as `varnames` a set of the
      following names that have to be defined as variables in the model:

      - the displacement variable which has to be defined as an unknown,
      - the plastic multiplier which has also defined as an unknown,
      - optionally the pressure variable for a mixed displacement-pressure
        formulation for 'DISPLACEMENT_AND_PLASTIC_MULTIPLIER_AND_PRESSURE'
        as `unknowns_type`,
      - the name of a (scalar) fem_data or im_data field that holds the
        plastic strain at the previous time step, and
      - the name of a fem_data or im_data field that holds all
        non-repeated components of the inverse of the plastic right
        Cauchy-Green tensor at the previous time step
        (it has to be a 4 element vector for plane strain 2D problems
        and a 6 element vector for 3D problems).

      The  "Simo_Miehe" law also expects as `params` a set of the
      following three parameters:

      - an expression for the initial bulk modulus K,
      - an expression for the initial shear modulus G,
      - the name of a user predefined function that decribes
        the yield limit as a function of the hardening variable
        (both the yield limit and the hardening variable values are
        assumed to be Frobenius norms of appropriate stress and strain
        tensors, respectively).

      As usual, `region` is an optional mesh region on which the term is added.
      If it is not specified, it is added on the whole mesh.
      Return the brick index in the model.@*/
    sub_command
      ("add finite strain elastoplasticity brick", 10, 11, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string lawname = in.pop().to_string();
       filter_lawname(lawname);
       size_type nb_var = 0; size_type nb_params = 0;
       if (lawname.compare("simo_miehe") == 0 ||
           lawname.compare("eterovic_bathe") == 0) {
         nb_var = 4;
         nb_params = 3;
       } else
         THROW_BADARG(lawname << " is not an implemented finite strain"
                              << " elastoplastic law");

       getfem::plasticity_unknowns_type unknowns_type(getfem::DISPLACEMENT_ONLY);
       mexarg_in argin = in.pop();
       if (argin.is_string()) {
         std::string opt = argin.to_string();
         filter_lawname(opt);
         if (opt.compare("displacement_and_plastic_multiplier") == 0)
           unknowns_type = getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER;
         else if (opt.compare("displacement_and_plastic_multiplier"
                              "_and_pressure") == 0)
           unknowns_type = getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER_AND_PRESSURE;
         else
           THROW_BADARG("Wrong input");
       } else if (argin.is_integer()) {
         unknowns_type = static_cast<getfem::plasticity_unknowns_type>
                         (argin.to_integer());
         GMM_ASSERT1
           (unknowns_type == getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER ||
            unknowns_type == getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER_AND_PRESSURE,
            "Not valid input for unknowns_type");
       }
       if (unknowns_type == getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER_AND_PRESSURE)
         nb_var += 1;

       std::vector<std::string> varnames;
       for (size_type i = 0; i < nb_var; ++i)
         varnames.push_back(in.pop().to_string());

       std::vector<std::string> params;
       for (size_type i = 0; i < nb_params; ++i)
         params.push_back(in.pop().to_string());

       size_type region(-1);
       if (in.remaining()) {
         argin = in.pop();
         if (!argin.is_integer())
           THROW_BADARG("Last optional argument must be an integer");
         region = argin.to_integer();
       }

       size_type ind = config::base_index() +
         add_finite_strain_elastoplasticity_brick
         (*md, *mim, lawname, unknowns_type, varnames, params, region);

       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add nonlinear incompressibility brick', @tmim mim, @str varname, @str multname_pressure[, @int region])
    Add a nonlinear incompressibility condition on `variable` (for large
    strain elasticity). `multname_pressure`
    is a variable which represent the pressure. Be aware that an inf-sup
    condition between the finite element method describing the pressure and the
    primal variable has to be satisfied. `region` is an optional mesh region on
    which the term is added. If it is not specified, it is added on the
    whole mesh. Return the brick index in the model.@*/
    sub_command
      ("add nonlinear incompressibility brick", 3, 4, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string multname = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind
       = getfem::add_nonlinear_incompressibility_brick
       (*md, *mim, varname, multname, region)
       + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add finite strain incompressibility brick', @tmim mim, @str varname, @str multname_pressure[, @int region])
    Add a finite strain incompressibility condition on `variable` (for large
    strain elasticity). `multname_pressure`
    is a variable which represent the pressure. Be aware that an inf-sup
    condition between the finite element method describing the pressure and the
    primal variable has to be satisfied. `region` is an optional mesh region on
    which the term is added. If it is not specified, it is added on the
    whole mesh. Return the brick index in the model.
    This brick is equivalent to the ``nonlinear incompressibility brick`` but
    uses the high-level generic assembly adding the term
    ``p*(1-Det(Id(meshdim)+Grad_u))`` if ``p`` is the multiplier and
    ``u`` the variable which represent the displacement.@*/
    sub_command
      ("add finite strain incompressibility brick", 3, 4, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string multname = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind
       = getfem::add_finite_strain_incompressibility_brick
       (*md, *mim, varname, multname, region)
       + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add bilaplacian brick', @tmim mim, @str varname, @str dataname [, @int region])
      Add a bilaplacian brick on the variable
      `varname` and on the mesh region `region`.
      This represent a term :math:`\Delta(D \Delta u)`.
      where :math:`D(x)` is a coefficient determined by `dataname` which
      could be constant or described on a f.e.m. The corresponding weak form
      is :math:`\int D(x)\Delta u(x) \Delta v(x) dx`.
      Return the brick index in the model.@*/
    sub_command
      ("add bilaplacian brick", 3, 4, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string dataname = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind = config::base_index() +
       add_bilaplacian_brick(*md, *mim,
                             varname, dataname, region);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );
    
    /*@SET ind = ('add Kirchhoff-Love plate brick', @tmim mim, @str varname, @str dataname_D, @str dataname_nu [, @int region])
      Add a bilaplacian brick on the variable
      `varname` and on the mesh region `region`.
      This represent a term :math:`\Delta(D \Delta u)` where :math:`D(x)`
      is a the flexion modulus determined by `dataname_D`. The term is
      integrated by part following a Kirchhoff-Love plate model
      with `dataname_nu` the poisson ratio.
      Return the brick index in the model.@*/
    sub_command
      ("add Kirchhoff-Love plate brick", 4, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string dataname_D = in.pop().to_string();
       std::string dataname_nu = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind = config::base_index() +
       add_bilaplacian_brick_KL(*md, *mim,
                                varname, dataname_D, dataname_nu, region);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add normal derivative source term brick', @tmim mim, @str varname, @str dataname, @int region)
      Add a normal derivative source term brick
      :math:`F = \int b.\partial_n v` on the variable `varname` and the
      mesh region `region`.

      Update the right hand side of the linear system.
      `dataname` represents `b` and `varname` represents `v`.
      Return the brick index in the model.@*/
    sub_command
      ("add normal derivative source term brick", 4, 4, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string dataname = in.pop().to_string();
       size_type region = in.pop().to_integer();
       size_type ind = config::base_index() +
       add_normal_derivative_source_term_brick(*md, *mim,
                                varname, dataname, region);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add Kirchhoff-Love Neumann term brick', @tmim mim, @str varname, @str dataname_M, @str dataname_divM, @int region)
       Add a Neumann term brick for Kirchhoff-Love model
      on the variable `varname` and the mesh region `region`.
      `dataname_M` represents the bending moment tensor and  `dataname_divM`
      its divergence.
      Return the brick index in the model.@*/
    sub_command
      ("add Kirchhoff-Love Neumann term brick", 5, 5, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string dataname_M = in.pop().to_string();
       std::string dataname_divM = in.pop().to_string();
       size_type region = in.pop().to_integer();
       size_type ind = config::base_index() +
       add_Kirchhoff_Love_Neumann_term_brick(*md, *mim,
                                varname, dataname_M, dataname_divM, region);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add normal derivative Dirichlet condition with multipliers', @tmim mim, @str varname, mult_description, @int region [, @str dataname, @int R_must_be_derivated])
       Add a Dirichlet condition on the normal derivative of the variable
      `varname` and on the mesh region `region` (which should be a boundary).
      The general form is
      :math:`\int \partial_n u(x)v(x) = \int r(x)v(x) \forall v`
      where :math:`r(x)` is
      the right hand side for the Dirichlet condition (0 for
      homogeneous conditions) and :math:`v` is in a space of multipliers
      defined by `mult_description`.
      If `mult_description` is a string this is assumed
      to be the variable name corresponding to the multiplier (which should be
      first declared as a multiplier variable on the mesh region in the model).
      If it is a finite element method (mesh_fem object) then a multiplier
      variable will be added to the model and build on this finite element
      method (it will be restricted to the mesh region `region` and eventually
      some conflicting dofs with some other multiplier variables will be
      suppressed). If it is an integer, then a  multiplier variable will be
      added to the model and build on a classical finite element of degree
      that integer. `dataname` is an optional parameter which represents
      the right hand side of the Dirichlet condition.
      If `R_must_be_derivated` is set to `true` then the normal
      derivative of `dataname` is considered.
      Return the brick index in the model.@*/
    sub_command
      ("add normal derivative Dirichlet condition with multipliers", 4, 6, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string multname;
       getfem::mesh_fem *mf = 0;
       size_type degree = 0;
       size_type version = 0;
       mexarg_in argin = in.pop();
       if (argin.is_integer()) {
         degree = argin.to_integer();
         version = 1;
       } else if (argin.is_string()) {
         multname = argin.to_string();
         version = 2;
       } else {
         mf = to_meshfem_object(argin);
         version = 3;
       }
       size_type region = in.pop().to_integer();
       std::string dataname;
       if (in.remaining()) dataname = in.pop().to_string();
       bool R_must_be_derivated = false;
       if (in.remaining())
         R_must_be_derivated = (in.pop().to_integer(0,1)) != 0;
       size_type ind = config::base_index();
       switch(version) {
       case 1:  ind +=
           add_normal_derivative_Dirichlet_condition_with_multipliers
           (*md, *mim, varname, dim_type(degree), region,
            dataname, R_must_be_derivated ); break;
       case 2:  ind +=
           add_normal_derivative_Dirichlet_condition_with_multipliers
           (*md, *mim, varname, multname, region,
            dataname, R_must_be_derivated ); break;
       case 3:  ind +=
           add_normal_derivative_Dirichlet_condition_with_multipliers
           (*md, *mim, varname,  *mf,
            region, dataname, R_must_be_derivated ); break;
       }
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add normal derivative Dirichlet condition with penalization', @tmim mim, @str varname, @scalar coeff, @int region [, @str dataname, @int R_must_be_derivated])
       Add a Dirichlet condition on the normal derivative of the variable
      `varname` and on the mesh region `region` (which should be a boundary).
      The general form is
      :math:`\int \partial_n u(x)v(x) = \int r(x)v(x) \forall v`
      where :math:`r(x)` is
      the right hand side for the Dirichlet condition (0 for
      homogeneous conditions).
      The penalization coefficient
      is initially `coeff` and will be added to the data of the model.
      It can be changed with the command MODEL:SET('change penalization coeff').
      `dataname` is an optional parameter which represents
      the right hand side of the Dirichlet condition.
      If `R_must_be_derivated` is set to `true` then the normal
      derivative of `dataname` is considered.
      Return the brick index in the model.@*/
    sub_command
      ("add normal derivative Dirichlet condition with penalization", 4, 6, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       double coeff = in.pop().to_scalar();
       size_type region = in.pop().to_integer();
       std::string dataname;
       if (in.remaining()) dataname = in.pop().to_string();
       bool R_must_be_derivated = false;
       if (in.remaining())
         R_must_be_derivated = (in.pop().to_integer(0,1)) != 0;
       size_type ind = config::base_index() +
       add_normal_derivative_Dirichlet_condition_with_penalization
           (*md, *mim, varname, coeff, region,
            dataname, R_must_be_derivated );
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add Mindlin Reissner plate brick', @tmim mim, @tmim mim_reduced, @str varname_u3, @str varname_theta , @str param_E, @str param_nu, @str param_epsilon, @str param_kappa [,@int variant [, @int region]])
      Add a term corresponding to the classical Reissner-Mindlin plate
      model for which `varname_u3` is the transverse displacement,
      `varname_theta` the rotation of
      fibers normal to the midplane, 'param_E' the Young Modulus,
      `param_nu` the poisson ratio,
      `param_epsilon` the plate thickness,
      `param_kappa` the shear correction factor. Note that since this brick
      uses the high level generic assembly language, the parameter can
      be regular expression of this language.
      There are three variants.
      `variant = 0` corresponds to the an
      unreduced formulation and in that case only the integration
      method `mim` is used. Practically this variant is not usable since
      it is subject to a strong locking phenomenon.
      `variant = 1` corresponds to a reduced integration where `mim` is
      used for the rotation term and `mim_reduced` for the transverse
      shear term. `variant = 2` (default) corresponds to the projection onto
      a rotated RT0 element of the transverse shear term. For the moment, this
      is adapted to quadrilateral only (because it is not sufficient to
      remove the locking phenomenon on triangle elements). Note also that if
      you use high order elements, the projection on RT0 will reduce the order
      of the approximation.
      Returns the brick index in the model.
      @*/
     sub_command
        ("add Mindlin Reissner plate brick", 7, 9, 0, 1,
         getfem::mesh_im *mim = to_meshim_object(in.pop());
         getfem::mesh_im *mim_reduced = to_meshim_object(in.pop());
         std::string varname_U = in.pop().to_string();
         std::string varname_theta = in.pop().to_string();
         std::string param_E = in.pop().to_string();
         std::string param_nu = in.pop().to_string();
         std::string param_epsilon = in.pop().to_string();
         std::string param_kapa = in.pop().to_string();
         size_type variant = size_type(2);
         if (in.remaining()) variant = in.pop().to_integer();
         size_type region = size_type(-1);
         if (in.remaining()) region = in.pop().to_integer();
         size_type ind = add_Mindlin_Reissner_plate_brick
         (*md, *mim, *mim_reduced,
          varname_U, varname_theta, param_E, param_nu, param_epsilon,
          param_kapa, variant, region);
         workspace().set_dependence(md, mim);
         out.pop().from_integer(int(ind));
         );
  
  
    /*@SET ind = ('add enriched Mindlin Reissner plate brick', @tmim mim, @tmim mim_reduced1, @tmim mim_reduced2, @str varname_ua, @str varname_theta,@str varname_u3, @str varname_theta3 , @str param_E, @str param_nu, @str param_epsilon [,@int variant [, @int region]])
    Add a term corresponding to the enriched Reissner-Mindlin plate
    model for which `varname_ua` is the membrane displacements,
    `varname_u3` is the transverse displacement,
    `varname_theta` the rotation of
    fibers normal to the midplane, 
    `varname_theta3` the pinching,     
    'param_E' the Young Modulus,
    `param_nu` the poisson ratio,
    `param_epsilon` the plate thickness. Note that since this brick
    uses the high level generic assembly language, the parameter can
    be regular expression of this language.
    There are four variants.
    `variant = 0` corresponds to the an
    unreduced formulation and in that case only the integration
    method `mim` is used. Practically this variant is not usable since
    it is subject to a strong locking phenomenon.
    `variant = 1` corresponds to a reduced integration where `mim` is
    used for the rotation term and `mim_reduced1` for the transverse
    shear term and `mim_reduced2` for the pinching term.
    `variant = 2` (default) corresponds to the projection onto
    a rotated RT0 element of the transverse shear term and a reduced integration for the pinching term.
    For the moment, this is adapted to quadrilateral only (because it is not sufficient to
    remove the locking phenomenon on triangle elements). Note also that if
    you use high order elements, the projection on RT0 will reduce the order
    of the approximation.
    `variant = 3` corresponds to the projection onto
    a rotated RT0 element of the transverse shear term and the projection onto P0 element of the pinching term.
    For the moment, this is adapted to quadrilateral only (because it is not sufficient to
    remove the locking phenomenon on triangle elements). Note also that if
    you use high order elements, the projection on RT0 will reduce the order
    of the approximation.   
    Returns the brick index in the model.
      @*/
     sub_command
        ("add enriched Mindlin Reissner plate brick", 10, 12, 0, 1,
         getfem::mesh_im *mim = to_meshim_object(in.pop());
         getfem::mesh_im *mim_reduced1 = to_meshim_object(in.pop());
         getfem::mesh_im *mim_reduced2 = to_meshim_object(in.pop());
         std::string varname_Ua = in.pop().to_string();
         std::string varname_theta = in.pop().to_string();
         std::string varname_U3 = in.pop().to_string();
         std::string varname_theta3 = in.pop().to_string();
         std::string param_E = in.pop().to_string();
         std::string param_nu = in.pop().to_string();
         std::string param_epsilon = in.pop().to_string();
         size_type variant = size_type(3);//2
         if (in.remaining()) variant = in.pop().to_integer();
         size_type region = size_type(-1);
         if (in.remaining()) region = in.pop().to_integer();
         size_type ind = add_enriched_Mindlin_Reissner_plate_brick
         (*md, *mim, *mim_reduced1,*mim_reduced2,
          varname_Ua, varname_theta, varname_U3, varname_theta3,
          param_E, param_nu, param_epsilon, variant, region);
         workspace().set_dependence(md, mim);
         out.pop().from_integer(int(ind));
         );
      

    /*@SET ind = ('add mass brick', @tmim mim, @str varname[, @str dataexpr_rho[, @int region]])
      Add mass term to the model relatively to the variable `varname`.
      If specified, the data `dataexpr_rho` is the
      density (1 if omitted). `region` is an optional mesh region on
      which the term is added. If it is not specified, it
      is added on the whole mesh. Return the brick index in the model.@*/
    sub_command
      ("add mass brick", 2, 4, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string dataname_rho;
       if (in.remaining()) dataname_rho = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind
       = getfem::add_mass_brick
       (*md, *mim, varname, dataname_rho, region)
       + config::base_index();
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ('shift variables for time integration')
      Function used to shift the variables of a model to the data
      corresponding of ther value on the previous time step for time
      integration schemes. For each variable for which a time integration
      scheme has been declared, the scheme is called to perform the shift.
      This function has to be called between two time steps. @*/
    sub_command
      ("shift variables for time integration", 0, 0, 0, 0,
       md->shift_variables_for_time_integration();
       );

    /*@SET ('perform init time derivative', @scalar ddt)
      By calling this function, indicates that the next solve will compute
      the solution for a (very) small time step `ddt` in order to initalize
      the data corresponding to the derivatives needed by time integration
      schemes (mainly the initial time derivative for order one in time
      problems  and the second order time derivative for second order in time
      problems). The next solve will not change the value of the variables. @*/
    sub_command
      ("perform init time derivative", 1, 1, 0, 0,
       double ddt = in.pop().to_scalar();
       md->perform_init_time_derivative(ddt);
       );

    /*@SET ('set time step', @scalar dt)
      Set the value of the time step to `dt`. This value can be change
      from a step to another for all one-step schemes (i.e. for the moment
      to all proposed time integration schemes). @*/
    sub_command
      ("set time step", 1, 1, 0, 0,
       double dt = in.pop().to_scalar();
       md->set_time_step(dt);
       );

    /*@SET ('set time', @scalar t)
      Set the value of the data `t` corresponding to the current time to `t`.
      @*/
    sub_command
      ("set time", 1, 1, 0, 0,
       double t = in.pop().to_scalar();
       md->set_time(t);
       );


    /*@SET ('add theta method for first order', @str varname, @scalar theta)
      Attach a theta method for the time discretization of the variable
      `varname`. Valid only if there is at most first order time derivative
      of the variable. @*/
    sub_command
      ("add theta method for first order", 2, 2, 0, 0,
       std::string varname = in.pop().to_string();
       double theta = in.pop().to_scalar();
       getfem::add_theta_method_for_first_order(*md, varname, theta);
       );

    /*@SET ('add theta method for second order', @str varname, @scalar theta)
      Attach a theta method for the time discretization of the variable
      `varname`. Valid only if there is at most second order time derivative
      of the variable. @*/
    sub_command
      ("add theta method for second order", 2, 2, 0, 0,
       std::string varname = in.pop().to_string();
       double theta = in.pop().to_scalar();
       getfem::add_theta_method_for_second_order(*md, varname, theta);
       );

    /*@SET ('add Newmark scheme', @str varname, @scalar beta, @scalar gamma)
      Attach a theta method for the time discretization of the variable
      `varname`. Valid only if there is at most second order time derivative
      of the variable. @*/
    sub_command
      ("add Newmark scheme", 3, 3, 0, 0,
       std::string varname = in.pop().to_string();
       double beta = in.pop().to_scalar();
       double gamma = in.pop().to_scalar();
       getfem::add_Newmark_scheme(*md, varname, beta, gamma);
       );

    /*@SET ('add_Houbolt_scheme', @str varname)
      Attach a Houbolt method for the time discretization of the variable
      `varname`. Valid only if there is at most second order time derivative
      of the variable  @*/
    sub_command
      ("add Houbolt scheme", 1, 1, 0, 0,
       std::string varname = in.pop().to_string();
       getfem::add_Houbolt_scheme(*md, varname);
       );

     /*@SET ('disable bricks', @ivec bricks_indices)
       Disable a brick (the brick will no longer participate to the
       building of the tangent linear system).@*/
     sub_command
       ("disable bricks", 1, 1, 0, 0,
        dal::bit_vector bv = in.pop().to_bit_vector();
        for (dal::bv_visitor ii(bv); !ii.finished(); ++ii)
          md->disable_brick(ii);
        );


     /*@SET ('enable bricks', @ivec bricks_indices)
       Enable a disabled brick. @*/
     sub_command
       ("enable bricks", 1, 1, 0, 0,
        dal::bit_vector bv = in.pop().to_bit_vector();
        for (dal::bv_visitor ii(bv); !ii.finished(); ++ii)
          md->enable_brick(ii);
        );

     /*@SET ('disable variable', @str varname)
       Disable a variable for a solve (and its attached multipliers).
       The next solve will operate only on
       the remaining variables. This allows to solve separately different
       parts of a model. If there is a strong coupling of the variables,
       a fixed point strategy can the be used. @*/
     sub_command
       ("disable variable", 1, 1, 0, 0,
        std::string varname = in.pop().to_string();
        md->disable_variable(varname);
        );


     /*@SET ('enable variable', @str varname)
       Enable a disabled variable (and its attached multipliers). @*/
     sub_command
       ("enable variable", 1, 1, 0, 0,
        std::string varname = in.pop().to_string();
        md->enable_variable(varname);
        );


     /*@SET ('first iter')
       To be executed before the first iteration of a time integration
       scheme. @*/
     sub_command
       ("first iter", 0, 0, 0, 0,
        md->first_iter();
        );


     /*@SET ('next iter')
       To be executed at the end of each iteration of a time
       integration scheme. @*/
     sub_command
       ("next iter", 0, 0, 0, 0,
        md->next_iter();
        );


     /*@SET ind = ('add basic contact brick', @str varname_u, @str multname_n[, @str multname_t], @str dataname_r, @tspmat BN[, @tspmat BT, @str dataname_friction_coeff][, @str dataname_gap[, @str dataname_alpha[, @int augmented_version[, @str dataname_gamma, @str dataname_wt]]])
       
     Add a contact with or without friction brick to the model.
     If U is the vector
     of degrees of freedom on which the unilateral constraint is applied,
     the matrix `BN` have to be such that this constraint is defined by
     :math:`B_N U \le 0`. A friction condition can be considered by adding
     the three parameters `multname_t`, `BT` and `dataname_friction_coeff`.
     In this case, the tangential displacement is :math:`B_T U` and
     the matrix `BT` should have as many rows as `BN` multiplied by
     :math:`d-1` where :math:`d` is the domain dimension.
     In this case also, `dataname_friction_coeff` is a data which represents
     the coefficient of friction. It can be a scalar or a vector representing a
     value on each contact condition.  The unilateral constraint is prescribed
     thank to a multiplier
     `multname_n` whose dimension should be equal to the number of rows of
     `BN`. If a friction condition is added, it is prescribed with a
     multiplier `multname_t` whose dimension should be equal to the number
     of rows of `BT`. The augmentation parameter `r` should be chosen in
     a range of
     acceptabe values (see Getfem user documentation). `dataname_gap` is an
     optional parameter representing the initial gap. It can be a single value
     or a vector of value. `dataname_alpha` is an optional homogenization
     parameter for the augmentation parameter
     (see Getfem user documentation).  The parameter `augmented_version`
     indicates the augmentation strategy : 1 for the non-symmetric
     Alart-Curnier augmented Lagrangian, 2 for the symmetric one (except for
     the coupling between contact and Coulomb friction), 3 for the
     unsymmetric method with augmented multipliers, 4 for the unsymmetric
     method with augmented multipliers and De Saxce projection. @*/
     sub_command
       ("add basic contact brick", 4, 12, 0, 1,

        bool friction = false;

        std::string varname_u = in.pop().to_string();
        std::string multname_n = in.pop().to_string();
        std::string dataname_r = in.pop().to_string();
        std::string multname_t; std::string friction_coeff;

        mexarg_in argin = in.pop();
        if (argin.is_string()) {
          friction = true;
          multname_t = dataname_r;
          dataname_r = argin.to_string();
          argin = in.pop();
        }

        std::shared_ptr<gsparse> BN = argin.to_sparse();
        if (BN->is_complex()) THROW_BADARG("Complex matrix not allowed");
        std::shared_ptr<gsparse> BT;
        if (friction) {
          BT = in.pop().to_sparse();
          if (BT->is_complex()) THROW_BADARG("Complex matrix not allowed");
          friction_coeff = in.pop().to_string();
        }

        std::string dataname_gap;
        dataname_gap = in.pop().to_string();
        std::string dataname_alpha;
        if (in.remaining()) dataname_alpha = in.pop().to_string();
        int augmented_version = 1;
        if (in.remaining()) augmented_version = in.pop().to_integer(1,4);

        std::string dataname_gamma;
        std::string dataname_wt;
        if (in.remaining()) {
          GMM_ASSERT1(friction,
                      "gamma and wt parameters are for the frictional brick only");
          dataname_gamma = in.pop().to_string();
          dataname_wt = in.pop().to_string();
        }

        getfem::CONTACT_B_MATRIX BBN;
        getfem::CONTACT_B_MATRIX BBT;
        if (BN->storage()==gsparse::CSCMAT) {
          gmm::resize(BBN, gmm::mat_nrows(BN->real_csc()),
                      gmm::mat_ncols(BN->real_csc()));
          gmm::copy(BN->real_csc(), BBN);
        }
        else if (BN->storage()==gsparse::WSCMAT) {
          gmm::resize(BBN, gmm::mat_nrows(BN->real_wsc()),
                      gmm::mat_ncols(BN->real_wsc()));
          gmm::copy(BN->real_wsc(), BBN);
        }
        else THROW_BADARG("Matrix BN should be a sparse matrix");

        if (friction) {
          if (BT->storage()==gsparse::CSCMAT) {
            gmm::resize(BBT, gmm::mat_nrows(BT->real_csc()),
                        gmm::mat_ncols(BT->real_csc()));
            gmm::copy(BT->real_csc(), BBT);
          }
          else if (BT->storage()==gsparse::WSCMAT) {
            gmm::resize(BBT, gmm::mat_nrows(BT->real_wsc()),
                        gmm::mat_ncols(BT->real_wsc()));
            gmm::copy(BT->real_wsc(), BBT);
          }
          else THROW_BADARG("Matrix BT should be a sparse matrix");
        }

        size_type ind;
        if (friction) {
          ind = getfem::add_basic_contact_brick
            (*md, varname_u, multname_n, multname_t, dataname_r, BBN, BBT,
             friction_coeff, dataname_gap, dataname_alpha, augmented_version,
             false, "", dataname_gamma, dataname_wt);
        } else {
          ind = getfem::add_basic_contact_brick
            (*md, varname_u, multname_n, dataname_r, BBN, dataname_gap,
             dataname_alpha, augmented_version);
        }

        out.pop().from_integer(int(ind + config::base_index()));
        );


     /*@SET ind = ('add basic contact brick two deformable bodies', @str varname_u1, @str varname_u2, @str multname_n, @str dataname_r, @tspmat BN1, @tspmat BN2[, @str dataname_gap[, @str dataname_alpha[, @int augmented_version]]])
       
     Add a frictionless contact condition to the model between two deformable
      bodies. If U1, U2 are the vector
      of degrees of freedom on which the unilateral constraint is applied,
      the matrices `BN1` and `BN2` have to be such that this condition
      is defined by
      $B_{N1} U_1 B_{N2} U_2 + \le gap$. The constraint is prescribed thank
      to a multiplier
      `multname_n` whose dimension should be equal to the number of lines of
      `BN`. The augmentation parameter `r` should be chosen in a range of
      acceptabe values (see Getfem user documentation). `dataname_gap` is an
      optional parameter representing the initial gap. It can be a single value
      or a vector of value. `dataname_alpha` is an optional homogenization
      parameter for the augmentation parameter
      (see Getfem user documentation). The parameter `aug_version` indicates
      the augmentation strategy : 1 for the non-symmetric Alart-Curnier
      augmented Lagrangian, 2 for the symmetric one, 3 for the unsymmetric
      method with augmented multiplier. @*/
     sub_command
       ("add basic contact brick two deformable bodies", 6, 9, 0, 1,

        std::string varname_u1 = in.pop().to_string();
        std::string varname_u2 = in.pop().to_string();
        std::string multname_n = in.pop().to_string();
        std::string dataname_r = in.pop().to_string();
       
        std::shared_ptr<gsparse> BN1 = in.pop().to_sparse();
        std::shared_ptr<gsparse> BN2 =  in.pop().to_sparse();
        if (BN1->is_complex()) THROW_BADARG("Complex matrix not allowed");
        if (BN2->is_complex()) THROW_BADARG("Complex matrix not allowed");

        std::string dataname_gap;
        if (in.remaining()) dataname_gap = in.pop().to_string();
        std::string dataname_alpha;
        if (in.remaining()) dataname_alpha = in.pop().to_string();
        int augmented_version = 1;
        if (in.remaining()) augmented_version = in.pop().to_integer(1,4);

        getfem::CONTACT_B_MATRIX BBN1; getfem::CONTACT_B_MATRIX BBN2;
        if (BN1->storage()==gsparse::CSCMAT) {
          gmm::resize(BBN1, gmm::mat_nrows(BN1->real_csc()),
                      gmm::mat_ncols(BN1->real_csc()));
          gmm::copy(BN1->real_csc(), BBN1);
        }
        else if (BN1->storage()==gsparse::WSCMAT) {
          gmm::resize(BBN1, gmm::mat_nrows(BN1->real_wsc()),
                      gmm::mat_ncols(BN1->real_wsc()));
          gmm::copy(BN1->real_wsc(), BBN1);
        }
        else THROW_BADARG("Matrix BN1 should be a sparse matrix");

        if (BN2->storage()==gsparse::CSCMAT) {
          gmm::resize(BBN2, gmm::mat_nrows(BN2->real_csc()),
                      gmm::mat_ncols(BN2->real_csc()));
          gmm::copy(BN2->real_csc(), BBN2);
        }
        else if (BN2->storage()==gsparse::WSCMAT) {
          gmm::resize(BBN2, gmm::mat_nrows(BN2->real_wsc()),
                      gmm::mat_ncols(BN2->real_wsc()));
          gmm::copy(BN2->real_wsc(), BBN2);
        }
        else THROW_BADARG("Matrix BN2 should be a sparse matrix");

        size_type ind;
        ind = getfem::add_basic_contact_brick_two_deformable_bodies
        (*md, varname_u1, varname_u2, multname_n, dataname_r, BBN1, BBN2,
         dataname_gap, dataname_alpha, augmented_version);

        out.pop().from_integer(int(ind + config::base_index()));
        );

    /*@SET ('contact brick set BN', @int indbrick, @tspmat BN)
    Can be used to set the BN matrix of a basic contact/friction brick. @*/
     sub_command
       ("contact brick set BN", 2, 2, 0, 0,
        size_type ind = in.pop().to_integer() - config::base_index();
        std::shared_ptr<gsparse> B = in.pop().to_sparse();

        if (B->is_complex())
          THROW_BADARG("BN should be a real matrix");

        if (B->storage()==gsparse::CSCMAT)
          gmm::copy(B->real_csc(),
                    getfem::contact_brick_set_BN(*md, ind));
        else if (B->storage()==gsparse::WSCMAT)
          gmm::copy(B->real_wsc(),
                    getfem::contact_brick_set_BN(*md, ind));
        else
          THROW_BADARG("BN should be a sparse matrix");
        );


    /*@SET ('contact brick set BT', @int indbrick, @tspmat BT)
      Can be used to set the BT matrix of a basic contact with
      friction brick. @*/
     sub_command
       ("contact brick set BT", 2, 2, 0, 0,
        size_type ind = in.pop().to_integer() - config::base_index();
        std::shared_ptr<gsparse> B = in.pop().to_sparse();

        if (B->is_complex())
          THROW_BADARG("BT should be a real matrix");

        if (B->storage()==gsparse::CSCMAT)
          gmm::copy(B->real_csc(), getfem::contact_brick_set_BT(*md, ind));
        else if (B->storage()==gsparse::WSCMAT)
          gmm::copy(B->real_wsc(), getfem::contact_brick_set_BT(*md, ind));
        else
          THROW_BADARG("BT should be a sparse matrix");
        );


    // CONTACT WITH RIGID OBSTACLE


    /*@SET ind = ('add nodal contact with rigid obstacle brick',  @tmim mim, @str varname_u, @str multname_n[, @str multname_t], @str dataname_r[, @str dataname_friction_coeff], @int region, @str obstacle[,  @int augmented_version])

    Add a contact with or without friction condition with a rigid obstacle
    to the model. The condition is applied on the variable `varname_u`
    on the boundary corresponding to `region`. The rigid obstacle should
    be described with the string `obstacle` being a signed distance to
    the obstacle. This string should be an expression where the coordinates
    are 'x', 'y' in 2D and 'x', 'y', 'z' in 3D. For instance, if the rigid
    obstacle correspond to :math:`z \le 0`, the corresponding signed distance
    will be simply "z". `multname_n` should be a fixed size variable whose size
    is the number of degrees of freedom on boundary `region`. It represents the
    contact equivalent nodal forces. In order to add a friction condition
    one has to add the `multname_t` and `dataname_friction_coeff` parameters.
    `multname_t` should be a fixed size variable whose size is
    the number of degrees of freedom on boundary `region` multiplied by
    :math:`d-1` where :math:`d` is the domain dimension. It represents
    the friction equivalent nodal forces.
    The augmentation parameter `r` should be chosen in a
    range of acceptabe values (close to the Young modulus of the elastic
    body, see Getfem user documentation).  `dataname_friction_coeff` is
    the friction coefficient. It could be a scalar or a vector of values
    representing the friction coefficient on each contact node. 
    The parameter `augmented_version`
    indicates the augmentation strategy : 1 for the non-symmetric
    Alart-Curnier augmented Lagrangian, 2 for the symmetric one (except for
    the coupling between contact and Coulomb friction),
    3 for the new unsymmetric method.
    Basically, this brick compute the matrix BN
    and the vectors gap and alpha and calls the basic contact brick. @*/
     sub_command
       ("add nodal contact with rigid obstacle brick", 6, 9, 0, 1,

        bool friction = false;

        getfem::mesh_im *mim = to_meshim_object(in.pop());
        std::string varname_u = in.pop().to_string();
        std::string multname_n = in.pop().to_string();
        std::string dataname_r = in.pop().to_string();
        std::string multname_t;  std::string dataname_fr;

        mexarg_in argin = in.pop();
        if (argin.is_string()) {
          friction = true;
          multname_t = dataname_r;
          dataname_r = argin.to_string();
          dataname_fr = in.pop().to_string();
          argin = in.pop();
        }

        size_type region = argin.to_integer();
        std::string obstacle = in.pop().to_string();
        int augmented_version = 1;
        if (in.remaining()) augmented_version = in.pop().to_integer(1,4);

        size_type ind;

        if (friction)
          ind = getfem::add_nodal_contact_with_rigid_obstacle_brick
            (*md, *mim, varname_u, multname_n,
             multname_t, dataname_r, dataname_fr, region, obstacle,
             augmented_version);
        else
          ind = getfem::add_nodal_contact_with_rigid_obstacle_brick
            (*md, *mim, varname_u, multname_n,
             dataname_r, region, obstacle, augmented_version);
        workspace().set_dependence(md, mim);
        out.pop().from_integer(int(ind + config::base_index()));
        );

    /*@SET ind = ('add contact with rigid obstacle brick',  @tmim mim, @str varname_u, @str multname_n[, @str multname_t], @str dataname_r[, @str dataname_friction_coeff], @int region, @str obstacle[,  @int augmented_version])
    DEPRECATED FUNCTION. Use 'add nodal contact with rigid obstacle brick' instead.@*/
     sub_command
       ("add contact with rigid obstacle brick", 6, 9, 0, 1,
        infomsg() << "WARNING : gf_mesh_fem_get('add contact with rigid obstacle "
        << "brick', ...) is a deprecated command.\n          Use gf_mesh_fem_get("
        << "'add nodal contact with rigid obstacle brick', ...) instead." << endl;
        SUBC_TAB::iterator it = subc_tab.find("add nodal contact with rigid obstacle brick");
        if (it != subc_tab.end())
            it->second->run(in, out, md);
        );

    /*@SET ind = ('add integral contact with rigid obstacle brick',  @tmim mim, @str varname_u, @str multname, @str dataname_obstacle, @str dataname_r [, @str dataname_friction_coeff], @int region [, @int option [, @str dataname_alpha [, @str dataname_wt [, @str dataname_gamma [, @str dataname_vt]]]]])

    Add a contact with or without friction condition with a rigid obstacle
    to the model. This brick adds a contact which is defined
    in an integral way. It is the direct approximation of an augmented
    Lagrangian formulation (see Getfem user documentation) defined at the
    continuous level. The advantage is a better scalability: the number of
    Newton iterations should be more or less independent of the mesh size.
    The contact condition is applied on the variable `varname_u`
    on the boundary corresponding to `region`. The rigid obstacle should
    be described with the data `dataname_obstacle` being a signed distance to
    the obstacle (interpolated on a finite element method).
    `multname` should be a fem variable representing the contact stress.
    An inf-sup condition beetween `multname` and `varname_u` is required.
    The augmentation parameter `dataname_r` should be chosen in a
    range of acceptabe values.
    The optional parameter `dataname_friction_coeff` is the friction
    coefficient which could be constant or defined on a finite element method.
    Possible values for `option` is 1 for the non-symmetric Alart-Curnier
    augmented Lagrangian method, 2 for the symmetric one, 3 for the
    non-symmetric Alart-Curnier method with an additional augmentation
    and 4 for a new unsymmetric method. The default value is 1.
    In case of contact with friction, `dataname_alpha` and `dataname_wt`
    are optional parameters to solve evolutionary friction problems.
    `dataname_gamma` and `dataname_vt` represent optional data for adding
    a parameter-dependent sliding velocity to the friction condition.
    @*/
     sub_command
       ("add integral contact with rigid obstacle brick", 6, 12, 0, 1,

        getfem::mesh_im *mim = to_meshim_object(in.pop());
        std::string varname_u = in.pop().to_string();
        std::string multname = in.pop().to_string();
        std::string dataname_obs = in.pop().to_string();
        std::string dataname_r = in.pop().to_string();

        size_type ind;
        int option = 1;
        mexarg_in argin = in.pop();
        if (argin.is_integer()) { // without friction
            size_type region = argin.to_integer();
            if (in.remaining()) option = in.pop().to_integer();

            ind = getfem::add_integral_contact_with_rigid_obstacle_brick
                    (*md, *mim, varname_u, multname,
                     dataname_obs, dataname_r, region, option);
        } else { // with friction
            std::string dataname_coeff = argin.to_string();
            size_type region = in.pop().to_integer();
            if (in.remaining()) option = in.pop().to_integer();
            std::string dataname_alpha = "";
            if (in.remaining()) dataname_alpha = in.pop().to_string();
            std::string dataname_wt = "";
            if (in.remaining()) dataname_wt = in.pop().to_string();
            std::string dataname_gamma = "";
            if (in.remaining()) dataname_gamma = in.pop().to_string();
            std::string dataname_vt = "";
            if (in.remaining()) dataname_vt = in.pop().to_string();

            ind = getfem::add_integral_contact_with_rigid_obstacle_brick
                    (*md, *mim, varname_u, multname,
                     dataname_obs, dataname_r, dataname_coeff, region, option,
                     dataname_alpha, dataname_wt, dataname_gamma, dataname_vt);
        }
        workspace().set_dependence(md, mim);
        out.pop().from_integer(int(ind + config::base_index()));
        );

    /*@SET ind = ('add penalized contact with rigid obstacle brick',  @tmim mim, @str varname_u, @str dataname_obstacle, @str dataname_r [, @str dataname_coeff], @int region [, @int option, @str dataname_lambda, [, @str dataname_alpha [, @str dataname_wt]]])

    Add a penalized contact with or without friction condition with a
    rigid obstacle to the model.
    The condition is applied on the variable `varname_u`
    on the boundary corresponding to `region`. The rigid obstacle should
    be described with the data `dataname_obstacle` being a signed distance to
    the obstacle (interpolated on a finite element method).
    The penalization parameter `dataname_r` should be chosen
    large enough to prescribe approximate non-penetration and friction
    conditions but not too large not to deteriorate too much the
    conditionning of the tangent system.
    `dataname_lambda` is an optional parameter used if option
    is 2. In that case, the penalization term is shifted by lambda (this
    allows the use of an Uzawa algorithm on the corresponding augmented
    Lagrangian formulation)
    @*/
     sub_command
       ("add penalized contact with rigid obstacle brick", 5, 10, 0, 1,

        getfem::mesh_im *mim = to_meshim_object(in.pop());
        std::string varname_u = in.pop().to_string();
        std::string dataname_obs = in.pop().to_string();
        std::string dataname_r = in.pop().to_string();

        size_type ind;
        int option = 1;
        mexarg_in argin = in.pop();
        if (argin.is_integer()) { // without friction
            size_type region = argin.to_integer();
            if  (in.remaining()) option = in.pop().to_integer();
            std::string dataname_n = "";
            if (in.remaining()) dataname_n = in.pop().to_string();

            ind = getfem::add_penalized_contact_with_rigid_obstacle_brick
                    (*md, *mim, varname_u,
                     dataname_obs, dataname_r, region, option, dataname_n);
        } else { // with friction
            std::string dataname_coeff = argin.to_string();
            size_type region = in.pop().to_integer();
            if (in.remaining()) option = in.pop().to_integer();
            std::string dataname_lambda = "";
            if (in.remaining()) dataname_lambda = in.pop().to_string();
            std::string dataname_alpha = "";
            if (in.remaining()) dataname_alpha = in.pop().to_string();
            std::string dataname_wt = "";
            if (in.remaining()) dataname_wt = in.pop().to_string();

            ind = getfem::add_penalized_contact_with_rigid_obstacle_brick
                    (*md, *mim, varname_u,
                     dataname_obs, dataname_r, dataname_coeff, region, option,
                     dataname_lambda, dataname_alpha, dataname_wt);
        }

        workspace().set_dependence(md, mim);
        out.pop().from_integer(int(ind + config::base_index()));
        );
     
         
     /*@SET ind = ('add Nitsche contact with rigid obstacle brick', @tmim mim, @str varname, @str Neumannterm, @str dataname_obstacle, @str gamma0name,  @int region[, @scalar theta[, @str dataname_friction_coeff[, @str dataname_alpha, @str dataname_wt]]])
      Adds a contact condition with or without Coulomb friction on the variable
      `varname` and the mesh boundary `region`. The contact condition
      is prescribed with Nitsche's method. The rigid obstacle should
      be described with the data `dataname_obstacle` being a signed distance to
      the obstacle (interpolated on a finite element method).
      `gamma0name` is the Nitsche's method parameter.
      `theta` is a scalar value which can be
      positive or negative. `theta = 1` corresponds to the standard symmetric
      method which is conditionally coercive for  `gamma0` small.
      `theta = -1` corresponds to the skew-symmetric method which is
      inconditionally coercive. `theta = 0` is the simplest method
      for which the second derivative of the Neumann term is not necessary.
      The optional parameter `dataname_friction_coeff` is the friction
      coefficient which could be constant or defined on a finite element
      method.
      CAUTION: This brick has to be added in the model after all the bricks
      corresponding to partial differential terms having a Neumann term.
      Moreover, This brick can only be applied to bricks declaring their
      Neumann terms. Returns the brick index in the model.
    @*/
    sub_command
      ("add Nitsche contact with rigid obstacle brick", 6, 10, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string Neumannterm = in.pop().to_string();
       std::string dataname_obs = in.pop().to_string();
       std::string gamma0name = in.pop().to_string();
       size_type region = in.pop().to_integer();

       scalar_type theta = scalar_type(1);
       std::string dataname_fr;
       if (in.remaining()) {
         mexarg_in argin = in.pop();
         if (argin.is_string())
           dataname_fr = argin.to_string();
         else
           theta = argin.to_scalar();
       }
       if (in.remaining()) dataname_fr = in.pop().to_string();
       std::string dataname_alpha;
       if (in.remaining()) dataname_alpha = in.pop().to_string();
       std::string dataname_wt;
       if (in.remaining()) dataname_wt = in.pop().to_string();

       size_type ind = config::base_index();
       ind += getfem::add_Nitsche_contact_with_rigid_obstacle_brick
       (*md, *mim, varname, Neumannterm, dataname_obs,
        gamma0name, theta,
        dataname_fr, dataname_alpha, dataname_wt, region);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

#ifdef EXPERIMENTAL_PURPOSE_ONLY

     /*@SET ind = ('add Nitsche midpoint contact with rigid obstacle brick', @tmim mim, @str varname, @str Neumannterm, @str Neumannterm_wt, @str dataname_obstacle, @str gamma0name,  @int region, @scalar theta, @str dataname_friction_coeff, @str dataname_alpha, @str dataname_wt)
      EXPERIMENTAL BRICK: for midpoint scheme only !!
      Adds a contact condition with or without Coulomb friction on the variable
      `varname` and the mesh boundary `region`. The contact condition
      is prescribed with Nitsche's method. The rigid obstacle should
      be described with the data `dataname_obstacle` being a signed distance to
      the obstacle (interpolated on a finite element method).
      `gamma0name` is the Nitsche's method parameter.
      `theta` is a scalar value which can be
      positive or negative. `theta = 1` corresponds to the standard symmetric
      method which is conditionally coercive for  `gamma0` small.
      `theta = -1` corresponds to the skew-symmetric method which is
      inconditionally coercive. `theta = 0` is the simplest method
      for which the second derivative of the Neumann term is not necessary.
      The optional parameter `dataname_friction_coeff` is the friction
      coefficient which could be constant or defined on a finite element
      method.
      Returns the brick index in the model.
      
    @*/
    sub_command
      ("add Nitsche midpoint contact with rigid obstacle brick", 11, 11, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string Neumannterm = in.pop().to_string();
       std::string Neumannterm_wt = in.pop().to_string();
       std::string dataname_obs = in.pop().to_string();
       std::string gamma0name = in.pop().to_string();
       size_type region = in.pop().to_integer();

       scalar_type theta = scalar_type(1);
       std::string dataname_fr;
       mexarg_in argin = in.pop();
       if (argin.is_string())
         dataname_fr = argin.to_string();
       else
         theta = argin.to_scalar();
       dataname_fr = in.pop().to_string();
       std::string dataname_alpha = in.pop().to_string();
       std::string dataname_wt = in.pop().to_string();

       size_type ind = config::base_index();
       ind += getfem::add_Nitsche_contact_with_rigid_obstacle_brick_modified_midpoint
       (*md, *mim, varname, Neumannterm, Neumannterm_wt,
        dataname_obs,
        gamma0name, theta,
        dataname_fr, dataname_alpha, dataname_wt, region);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

#endif

#if (0) // Deprecated brick : uses the old Neumann terms

    /*@SET ind = ('add Nitsche fictitious domain contact brick', @tmim mim, @str varname1, @str varname2, @str dataname_d1, @str dataname_d2, @str gamma0name [, @scalar theta[, @str dataname_friction_coeff[, @str dataname_alpha, @str dataname_wt1,@str dataname_wt2]]])
     Adds a contact condition with or without Coulomb friction between
     two bodies in a fictitious domain. The contact condition is applied on 
     the variable `varname_u1` corresponds with the first and slave body 
     with Nitsche's method and on the variable `varname_u2` corresponds 
     with the second and master body with Nitsche's method. 
     The contact condition is evaluated on the fictitious slave boundary.
     The first body should be described by the level-set `dataname_d1` 
     and the second body should be described by the level-set `dataname_d2`.
     `gamma0name` is the Nitsche's method parameter. 
     `theta` is a scalar value which can be positive or negative. 
     `theta = 1` corresponds to the standard symmetric method which is
     conditionally coercive for  `gamma0` small.
     `theta = -1` corresponds to the skew-symmetric method which is inconditionally coercive.
     `theta = 0` is the simplest method for which the second derivative of
     the Neumann term is not necessary. The optional parameter `dataname_friction_coeff`
     is the friction coefficient which could be constant or defined on a finite element method. 
     CAUTION: This brick has to be added in the model after all the bricks
     corresponding to partial differential terms having a Neumann term.
     Moreover, This brick can only be applied to bricks declaring their
     Neumann terms. Returns the brick index in the model. 
    @*/
    sub_command
      ("add Nitsche fictitious domain contact brick", 6, 11, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname1 = in.pop().to_string();
       std::string varname2 = in.pop().to_string();
       std::string dataname_d1 = in.pop().to_string();
       std::string dataname_d2 = in.pop().to_string();
       std::string gamma0name = in.pop().to_string();

       scalar_type theta = scalar_type(1);
       std::string dataname_fr;
       if (in.remaining()) {
         mexarg_in argin = in.pop();
         if (argin.is_string())
           dataname_fr = argin.to_string();
         else
           theta = argin.to_scalar();
       }
       if (in.remaining()) dataname_fr = in.pop().to_string();
       std::string dataname_alpha;
       if (in.remaining()) dataname_alpha = in.pop().to_string();
       std::string dataname_wt1;
       if (in.remaining()) dataname_wt1 = in.pop().to_string();
       std::string dataname_wt2;
       if (in.remaining()) dataname_wt2 = in.pop().to_string();

       size_type ind = config::base_index();
       ind += getfem::add_Nitsche_fictitious_domain_contact_brick
       (*md, *mim, varname1, varname2, dataname_d1,
        dataname_d2, gamma0name, theta,
        dataname_fr, dataname_alpha, dataname_wt1, dataname_wt2);
       workspace().set_dependence(md, mim);
       out.pop().from_integer(int(ind));
       );

#endif

    // CONTACT BETWEEN NON-MATCHING MESHES


    /*@SET ind = ('add nodal contact between nonmatching meshes brick',  @tmim mim1[, @tmim mim2], @str varname_u1[, @str varname_u2], @str multname_n[, @str multname_t], @str dataname_r[, @str dataname_fr], @int rg1, @int rg2[, @int slave1, @int slave2,  @int augmented_version])

    Add a contact with or without friction condition between two faces of
    one or two elastic bodies. The condition is applied on the variable
    `varname_u1` or the variables `varname_u1` and `varname_u2` depending
    if a single or two distinct displacement fields are given. Integers
    `rg1` and `rg2` represent the regions expected to come in contact with
    each other. In the single displacement variable case the regions defined
    in both `rg1` and `rg2` refer to the variable `varname_u1`. In the case
    of two displacement variables, `rg1` refers to `varname_u1` and `rg2`
    refers to `varname_u2`. `multname_n` should be a fixed size variable
    whose size is the number of degrees of freedom on those regions among
    the ones defined in `rg1` and `rg2` which are characterized as "slaves".
    It represents the contact equivalent nodal normal forces. `multname_t`
    should be a fixed size variable whose size corresponds to the size of
    `multname_n` multiplied by qdim - 1 . It represents the contact
    equivalent nodal tangent (frictional) forces. The augmentation parameter
    `r` should be chosen in a range of acceptabe values (close to the Young
    modulus of the elastic body, see Getfem user documentation). The
    friction coefficient stored in the parameter `fr` is either a single
    value or a vector of the same size as `multname_n`. The optional
    parameters `slave1` and `slave2` declare if the regions defined in `rg1`
    and `rg2` are correspondingly considered as "slaves". By default
    `slave1` is true and `slave2` is false, i.e. `rg1` contains the slave
    surfaces, while 'rg2' the master surfaces. Preferrably only one of
    `slave1` and `slave2` is set to true.  The parameter `augmented_version`
    indicates the augmentation strategy : 1 for the non-symmetric
    Alart-Curnier augmented Lagrangian, 2 for the symmetric one (except for
    the coupling between contact and Coulomb friction),
    3 for the new unsymmetric method.
    Basically, this brick computes the matrices BN and BT and the vectors
    gap and alpha and calls the basic contact brick. @*/
     sub_command
       ("add nodal contact between nonmatching meshes brick", 6, 13, 0, 1,

        bool two_variables = true;
        bool friction = false;

        getfem::mesh_im *mim1 = 0;
        getfem::mesh_im *mim2 = 0;
        std::string varname_u1;
        std::string varname_u2;
        bool slave1=true; bool slave2=false;
        int augmented_version = 1;

        mim1 = to_meshim_object(in.pop());
        mexarg_in argin = in.pop();
        if (argin.is_string()) {
          two_variables = false;
          mim2 = mim1;
          varname_u1 = argin.to_string();
          varname_u2 = varname_u1;
        } else {
          mim2 = to_meshim_object(argin);
          varname_u1 = in.pop().to_string();
          varname_u2 = in.pop().to_string();
	  cout << "ok here" << endl;
        }
        std::string multname_n = in.pop().to_string();
        std::string multname_t;
        std::string dataname_r = in.pop().to_string();
        std::string dataname_fr;
        argin = in.pop();
        if (!argin.is_integer()) {
          friction = true;
          multname_t = dataname_r;
          dataname_r = in.pop().to_string();
          dataname_fr = in.pop().to_string();
          argin = in.pop();
        }
        std::vector<size_type> vrg1(1,argin.to_integer());
        std::vector<size_type> vrg2(1,in.pop().to_integer());
        if (in.remaining()) slave1 = (in.pop().to_integer(0,1)) != 0;
        if (in.remaining()) slave2 = (in.pop().to_integer(0,1)) != 0;
        if (in.remaining()) augmented_version = in.pop().to_integer(1,4);

        size_type ind;
        if (!friction)
          ind = getfem::add_nodal_contact_between_nonmatching_meshes_brick
            (*md, *mim1, *mim2,
             varname_u1, varname_u2, multname_n, dataname_r,
             vrg1, vrg2, slave1, slave2, augmented_version);
        else
          ind = getfem::add_nodal_contact_between_nonmatching_meshes_brick
            (*md, *mim1, *mim2,
             varname_u1, varname_u2, multname_n, multname_t,
             dataname_r, dataname_fr,
             vrg1, vrg2, slave1, slave2, augmented_version);
        workspace().set_dependence(md, mim1);
        // if (two_variables)
        //   workspace().set_dependence(md, mim2);
        out.pop().from_integer(int(ind + config::base_index()));
        );

    /*@SET ind = ('add nonmatching meshes contact brick',  @tmim mim1[, @tmim mim2], @str varname_u1[, @str varname_u2], @str multname_n[, @str multname_t], @str dataname_r[, @str dataname_fr], @int rg1, @int rg2[, @int slave1, @int slave2,  @int augmented_version])
    DEPRECATED FUNCTION. Use 'add nodal contact between nonmatching meshes brick' instead.@*/
     sub_command
       ("add nonmatching meshes contact brick", 6, 13, 0, 1,
        infomsg() << "WARNING : gf_mesh_fem_get('add nonmatching meshes "
        << "contact brick', ...) is a deprecated command.\n          Use "
        << "gf_mesh_fem_get('add nodal contact between nonmatching meshes "
        << "brick', ...) instead." << endl;
        SUBC_TAB::iterator it = subc_tab.find("add nodal contact between nonmatching meshes brick");
        if (it != subc_tab.end())
            it->second->run(in, out, md);
        );

    /*@SET ind = ('add integral contact between nonmatching meshes brick',  @tmim mim, @str varname_u1, @str varname_u2, @str multname, @str dataname_r [, @str dataname_friction_coeff], @int region1, @int region2 [, @int option [, @str dataname_alpha [, @str dataname_wt1 , @str dataname_wt2]]])

    Add a contact with or without friction condition between nonmatching
    meshes to the model. This brick adds a contact which is defined
    in an integral way. It is the direct approximation of an augmented
    agrangian formulation (see Getfem user documentation) defined at the
    continuous level. The advantage should be a better scalability:
    the number of Newton iterations should be more or less independent
    of the mesh size.
    The condition is applied on the variables `varname_u1` and `varname_u2`
    on the boundaries corresponding to `region1` and `region2`.
    `multname` should be a fem variable representing the contact stress
    for the frictionless case and the contact and friction stress for the
    case with friction. An inf-sup condition between `multname` and
    `varname_u1` and `varname_u2` is required.
    The augmentation parameter `dataname_r` should be chosen in a
    range of acceptable values.
    The optional parameter `dataname_friction_coeff` is the friction
    coefficient which could be constant or defined on a finite element
    method on the same mesh as `varname_u1`.
    Possible values for `option` is 1 for the non-symmetric Alart-Curnier
    augmented Lagrangian method, 2 for the symmetric one, 3 for the
    non-symmetric Alart-Curnier method with an additional augmentation
    and 4 for a new unsymmetric method. The default value is 1.
    In case of contact with friction, `dataname_alpha`, `dataname_wt1` and
    `dataname_wt2` are optional parameters to solve evolutionary friction
    problems.
    @*/
     sub_command
       ("add integral contact between nonmatching meshes brick", 7, 12, 0, 1,

        getfem::mesh_im *mim = to_meshim_object(in.pop());
        std::string varname_u1 = in.pop().to_string();
        std::string varname_u2 = in.pop().to_string();
        std::string multname = in.pop().to_string();
        std::string dataname_r = in.pop().to_string();

        size_type ind;
        int option = 1;
        mexarg_in argin = in.pop();
        if (argin.is_integer()) { // without friction
            size_type region1 = argin.to_integer();
            size_type region2 = in.pop().to_integer();
            if (in.remaining()) option = in.pop().to_integer();

            ind = getfem::add_integral_contact_between_nonmatching_meshes_brick
                    (*md, *mim, varname_u1, varname_u2,
                     multname, dataname_r, region1, region2, option);
        } else { // with friction
            std::string dataname_coeff = argin.to_string();
            size_type region1 = in.pop().to_integer();
            size_type region2 = in.pop().to_integer();
            if (in.remaining()) option = in.pop().to_integer();
            std::string dataname_alpha = "";
            if (in.remaining()) dataname_alpha = in.pop().to_string();
            std::string dataname_wt1 = "";
            if (in.remaining()) dataname_wt1 = in.pop().to_string();
            std::string dataname_wt2 = "";
            if (in.remaining()) dataname_wt2 = in.pop().to_string();

            ind = getfem::add_integral_contact_between_nonmatching_meshes_brick
                    (*md, *mim, varname_u1, varname_u2,
                     multname, dataname_r, dataname_coeff, region1, region2,
                     option, dataname_alpha, dataname_wt1, dataname_wt2);
        }
        workspace().set_dependence(md, mim);
        out.pop().from_integer(int(ind + config::base_index()));
        );

    /*@SET ind = ('add penalized contact between nonmatching meshes brick',  @tmim mim, @str varname_u1, @str varname_u2, @str dataname_r [, @str dataname_coeff], @int region1, @int region2 [, @int option [, @str dataname_lambda, [, @str dataname_alpha [, @str dataname_wt1, @str dataname_wt2]]]])

    Add a penalized contact condition with or without friction between
    nonmatching meshes to the model.
    The condition is applied on the variables `varname_u1` and  `varname_u2`
    on the boundaries corresponding to `region1` and `region2`.
    The penalization parameter `dataname_r` should be chosen
    large enough to prescribe approximate non-penetration and friction
    conditions but not too large not to deteriorate too much the
    conditionning of the tangent system.
    The optional parameter `dataname_friction_coeff` is the friction
    coefficient which could be constant or defined on a finite element
    method on the same mesh as `varname_u1`.
    `dataname_lambda` is an optional parameter used if option
    is 2. In that case, the penalization term is shifted by lambda (this
    allows the use of an Uzawa algorithm on the corresponding augmented
    Lagrangian formulation)
    In case of contact with friction, `dataname_alpha`, `dataname_wt1` and
    `dataname_wt2` are optional parameters to solve evolutionary friction
    problems.
    @*/
     sub_command
       ("add penalized contact between nonmatching meshes brick", 6, 12, 0, 1,

        getfem::mesh_im *mim = to_meshim_object(in.pop());
        std::string varname_u1 = in.pop().to_string();
        std::string varname_u2 = in.pop().to_string();
        std::string dataname_r = in.pop().to_string();

        size_type ind;
        int option = 1;
        mexarg_in argin = in.pop();
        if (argin.is_integer()) { // without friction
            size_type region1 = argin.to_integer();
            size_type region2 = in.pop().to_integer();
            if  (in.remaining()) option = in.pop().to_integer();
            std::string dataname_n = "";
            if (in.remaining()) dataname_n = in.pop().to_string();

            ind = getfem::add_penalized_contact_between_nonmatching_meshes_brick
                    (*md, *mim, varname_u1, varname_u2,
                     dataname_r, region1, region2, option, dataname_n);
        } else { // with friction
            std::string dataname_coeff = argin.to_string();
            size_type region1 = in.pop().to_integer();
            size_type region2 = in.pop().to_integer();
            if (in.remaining()) option = in.pop().to_integer();
            std::string dataname_lambda = "";
            if (in.remaining()) dataname_lambda = in.pop().to_string();
            std::string dataname_alpha = "";
            if (in.remaining()) dataname_alpha = in.pop().to_string();
            std::string dataname_wt1 = "";
            if (in.remaining()) dataname_wt1 = in.pop().to_string();
            std::string dataname_wt2 = "";
            if (in.remaining()) dataname_wt2 = in.pop().to_string();

            ind = getfem::add_penalized_contact_between_nonmatching_meshes_brick
                    (*md, *mim, varname_u1, varname_u2,
                     dataname_r, dataname_coeff, region1, region2, option,
                     dataname_lambda, dataname_alpha, dataname_wt1, dataname_wt2);
        }

        workspace().set_dependence(md, mim);
        out.pop().from_integer(int(ind + config::base_index()));
        );

     /*@SET ind = ('add integral large sliding contact brick raytracing', @str dataname_r, @scalar release_distance, [, @str dataname_fr[, @str dataname_alpha[, @int version]]])
      Adds a large sliding contact with friction brick to the model.
      This brick is able to deal with self-contact, contact between
      several deformable bodies and contact with rigid obstacles.
      It uses the high-level generic assembly. It adds to the model
      a raytracing_interpolate_transformation object.
      For each slave boundary a multiplier variable should be defined.
      The release distance should be determined with care
      (generally a few times a mean element size, and less than the
      thickness of the body). Initially, the brick is added with no contact
      boundaries. The contact boundaries and rigid bodies are added with
      special functions. `version` is 0 (the default value) for the
      non-symmetric version and 1 for the more symmetric one
      (not fully symmetric even without friction). @*/

     sub_command
       ("add integral large sliding contact brick raytracing", 2, 6, 0, 1,
        
        std::string dataname_r = in.pop().to_string();
        scalar_type d = in.pop().to_scalar();
        std::string dataname_fr = "0";
        if (in.remaining()) dataname_fr = in.pop().to_string();
        if (dataname_fr.size() == 0) dataname_fr = "0";
        std::string dataname_alpha = "1";
        if (in.remaining()) dataname_alpha = in.pop().to_string();
        if (dataname_alpha.size() == 0) dataname_alpha = "1";
        bool sym_v = false;
        if (in.remaining()) sym_v = (in.pop().to_integer() != 0);
        bool frame_indifferent = false;
        if (in.remaining()) frame_indifferent = (in.pop().to_integer() != 0);

        size_type  ind
        = getfem::add_integral_large_sliding_contact_brick_raytracing
        (*md, dataname_r, d, dataname_fr, dataname_alpha, sym_v,
         frame_indifferent);
        out.pop().from_integer(int(ind + config::base_index()));
        );

     /*@SET ('add rigid obstacle to large sliding contact brick', @int indbrick, @str expr, @int N)
      Adds a rigid obstacle to an existing large sliding contact
      with friction brick. `expr` is an expression using the high-level
      generic assembly language (where `x` is the current point n the mesh)
      which should be a signed distance to the obstacle.
      `N` is the mesh dimension. @*/
    sub_command
      ("add rigid obstacle to large sliding contact brick", 3, 3, 0, 0,
       size_type ind = in.pop().to_integer() - config::base_index();
       std::string expr = in.pop().to_string();
       size_type N = in.pop().to_integer();

       add_rigid_obstacle_to_large_sliding_contact_brick
       (*md, ind, expr, N);
       );

     /*@SET ('add master contact boundary to large sliding contact brick', @int indbrick, @tmim mim, @int region, @str dispname[, @str wname])
      Adds a master contact boundary to an existing large sliding contact
      with friction brick. @*/
    sub_command
      ("add master contact boundary to large sliding contact brick", 4, 5, 0,0,
       size_type ind = in.pop().to_integer() - config::base_index();
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       size_type region = in.pop().to_integer();
       std::string dispname = in.pop().to_string();
       std::string wname;
       if (in.remaining()) wname = in.pop().to_string();

       add_contact_boundary_to_large_sliding_contact_brick
       (*md, ind, *mim, region, true, false,
        dispname, "", wname);
       );


    /*@SET ('add slave contact boundary to large sliding contact brick', @int indbrick, @tmim mim, @int region, @str dispname, @str lambdaname[, @str wname])
      Adds a slave contact boundary to an existing large sliding contact
      with friction brick. @*/
    sub_command
      ("add slave contact boundary to large sliding contact brick", 5, 6, 0,0,
       size_type ind = in.pop().to_integer() - config::base_index();
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       size_type region = in.pop().to_integer();
       std::string dispname = in.pop().to_string();
       std::string lambda = in.pop().to_string();
       std::string wname;
       if (in.remaining()) wname = in.pop().to_string();

       add_contact_boundary_to_large_sliding_contact_brick
       (*md, ind, *mim, region, false, true,
        dispname, lambda, wname);
       );

    /*@SET ('add master slave contact boundary to large sliding contact brick', @int indbrick, @tmim mim, @int region, @str dispname, @str lambdaname[, @str wname])
      Adds a contact boundary to an existing large sliding contact
      with friction brick which is both master and slave
      (allowing the self-contact). @*/
    sub_command
      ("add master slave contact boundary to large sliding contact brick",
       5, 6, 0, 0,
       size_type ind = in.pop().to_integer() - config::base_index();
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       size_type region = in.pop().to_integer();
       std::string dispname = in.pop().to_string();
       std::string lambda = in.pop().to_string();
       std::string wname;
       if (in.remaining()) wname = in.pop().to_string();

       add_contact_boundary_to_large_sliding_contact_brick
       (*md, ind, *mim, region, true, true,
        dispname, lambda, wname);
       );

       /*@SET ind = ('add Nitsche large sliding contact brick raytracing', @bool unbiased_version, @str dataname_r, @scalar release_distance[, @str dataname_fr[, @str dataname_alpha[, @int version]]])
      Adds a large sliding contact with friction brick to the model based on the Nitsche's method.
      This brick is able to deal with self-contact, contact between
      several deformable bodies and contact with rigid obstacles.
      It uses the high-level generic assembly. It adds to the model
      a raytracing_interpolate_transformation object. "unbiased_version" refers to the version of Nische's method to be used.
      (unbiased or biased one).
      For each slave boundary a  material law should be defined as a function of the dispacement variable on this boundary.
      The release distance should be determined with care
      (generally a few times a mean element size, and less than the
      thickness of the body). Initially, the brick is added with no contact
      boundaries. The contact boundaries and rigid bodies are added with
      special functions. `version` is 0 (the default value) for the
      non-symmetric version and 1 for the more symmetric one
      (not fully symmetric even without friction). @*/

     sub_command
       ("add Nitsche large sliding contact brick raytracing", 2, 6, 0, 1,
        
        bool unbiased=(in.pop().to_integer() != 0);
        std::string dataname_r = in.pop().to_string();
        scalar_type d = in.pop().to_scalar();
        std::string dataname_fr = "0";
        if (in.remaining()) dataname_fr = in.pop().to_string();
        if (dataname_fr.size() == 0) dataname_fr = "0";
        std::string dataname_alpha = "1";
        if (in.remaining()) dataname_alpha = in.pop().to_string();
        if (dataname_alpha.size() == 0) dataname_alpha = "1";
        bool sym_v = false;
        if (in.remaining()) sym_v = (in.pop().to_integer() != 0);
        bool frame_indifferent = false;
        if (in.remaining()) frame_indifferent = (in.pop().to_integer() != 0);

        size_type  ind
        = getfem::add_Nitsche_large_sliding_contact_brick_raytracing
        (*md, unbiased, dataname_r, d, dataname_fr, dataname_alpha, sym_v,
         frame_indifferent);
        out.pop().from_integer(int(ind + config::base_index()));
        );

     /*@SET ('add rigid obstacle to Nitsche large sliding contact brick', @int indbrick, @str expr, @int N)
      Adds a rigid obstacle to an existing large sliding contact
      with friction brick. `expr` is an expression using the high-level
      generic assembly language (where `x` is the current point n the mesh)
      which should be a signed distance to the obstacle.
      `N` is the mesh dimension. @*/
    sub_command
      ("add rigid obstacle to Nitsche large sliding contact brick", 3, 3, 0, 0,
       size_type ind = in.pop().to_integer() - config::base_index();
       std::string expr = in.pop().to_string();
       size_type N = in.pop().to_integer();

       add_rigid_obstacle_to_Nitsche_large_sliding_contact_brick
       (*md, ind, expr, N);
       );

     /*@SET ('add master contact boundary to biased Nitsche large sliding contact brick', @int indbrick, @tmim mim, @int region, @str dispname[, @str wname])
      Adds a master contact boundary to an existing biased Nitsche's large sliding contact
      with friction brick. @*/
    sub_command
      ("add master contact boundary to biased Nitsche large sliding contact brick", 4, 5, 0,0,
       size_type ind = in.pop().to_integer() - config::base_index();
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       size_type region = in.pop().to_integer();
       std::string dispname = in.pop().to_string();
       std::string wname;
       if (in.remaining()) wname = in.pop().to_string();

       add_contact_boundary_to_Nitsche_large_sliding_contact_brick
       (*md, ind, *mim, region, true, false, false,
        dispname, "", wname);
       );


    /*@SET ('add slave contact boundary to biased Nitsche large sliding contact brick', @int indbrick, @tmim mim, @int region, @str dispname, @str lambdaname[, @str wname])
      Adds a slave contact boundary to an existing biased Nitsche's large sliding contact
      with friction brick. @*/
    sub_command
      ("add slave contact boundary to biased Nitsche large sliding contact brick", 5, 6, 0,0,
       size_type ind = in.pop().to_integer() - config::base_index();
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       size_type region = in.pop().to_integer();
       std::string dispname = in.pop().to_string();
       std::string sigma = in.pop().to_string();
       std::string wname;
       if (in.remaining()) wname = in.pop().to_string();

       add_contact_boundary_to_Nitsche_large_sliding_contact_brick
       (*md, ind, *mim, region, false, true, false,
        dispname, sigma, wname);
       );

    /*@SET ('add contact boundary to unbiased Nitsche large sliding contact brick', @int indbrick, @tmim mim, @int region, @str dispname, @str lambdaname[, @str wname])
      Adds a contact boundary to an existing unbiased Nitschelarge sliding contact
      with friction brick which is both master and slave. @*/
    sub_command
      ("add contact boundary to unbiased Nitsche large sliding contact brick",
       5, 6, 0, 0,
       size_type ind = in.pop().to_integer() - config::base_index();
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       size_type region = in.pop().to_integer();
       std::string dispname = in.pop().to_string();
       std::string sigma = in.pop().to_string();
       std::string wname;
       if (in.remaining()) wname = in.pop().to_string();

       add_contact_boundary_to_Nitsche_large_sliding_contact_brick
       (*md, ind, *mim, region, true, true, true,
        dispname, sigma, wname);
       );
  }

  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  getfem::model *md  = to_model_object(m_in.pop());
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);


  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
              it->second->arg_in_max, it->second->arg_out_min,
              it->second->arg_out_max);
    it->second->run(m_in, m_out, md);
  }
  else bad_cmd(init_cmd);

}
