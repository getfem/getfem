/*===========================================================================

 Copyright (C) 2006-2020 Yves Renard, Julien Pommier.

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
// $Id$
#include <getfem_interface.h>
#include <getfemint.h>

using namespace getfemint;

void gf_workspace(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_delete(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_eltm(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_geotrans(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_geotrans_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_integ(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_integ_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_global_function(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_global_function_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_fem(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_fem_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_cont_struct(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_cont_struct_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_cvstruct_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesher_object(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesher_object_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_fem(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_fem_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_fem_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_im(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_im_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_im_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_im_data(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_im_data_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_im_data_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_model(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_model_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_model_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_slice(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_slice_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_slice_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_levelset(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_levelset_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_levelset_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_levelset(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_levelset_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_levelset_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_precond(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_precond_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_spmat(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_spmat_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_spmat_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_asm(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_compute(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_linsolve(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_util(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_exit(getfemint::mexargs_in&, getfemint::mexargs_out&) { exit(0); }

namespace getfemint {
  std::stringstream *global_pinfomsg = 0;
  std::ostream& infomsg() {
    return *global_pinfomsg;
  }

  config *config::cfg = 0;
  config::config(gfi_interface_type t) : current_function_(0) {
    switch (t) {
    case MATLAB_INTERFACE:
      base_index_ = 1;
      has_native_sparse_ = true;
      prefer_native_sparse_ = true;
      can_return_integer_ = false;
      has_1D_arrays_ = false;
      break;
    case PYTHON_INTERFACE:
      base_index_ = 0;
      has_native_sparse_ = false;
      prefer_native_sparse_ = false;
      can_return_integer_ = true;
      has_1D_arrays_ = true;
      break;
    case SCILAB_INTERFACE:
      base_index_ = 1;
      has_native_sparse_ = true;
      prefer_native_sparse_ = true;
      can_return_integer_ = false;
      has_1D_arrays_ = false;
      break;
    default:
      THROW_INTERNAL_ERROR;
    }
  }
}

typedef void (* psub_command)(getfemint::mexargs_in& in, getfemint::mexargs_out& out);


extern "C"
char* getfem_interface_main(int config_id, const char *function,
                            int nb_in_args, const gfi_array *in_args[],
                            int *nb_out_args, gfi_array ***pout_args,
			    char **pinfomsg, int scilab_flag) {

  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {
    subc_tab["workspace"] = gf_workspace;
    subc_tab["delete"] = gf_delete;
    subc_tab["eltm"] = gf_eltm;
    subc_tab["geotrans"] = gf_geotrans;
    subc_tab["geotrans_get"] = gf_geotrans_get;
    subc_tab["integ"] = gf_integ;
    subc_tab["integ_get"] = gf_integ_get;
    subc_tab["global_function"] = gf_global_function;
    subc_tab["global_function_get"] = gf_global_function_get;
    subc_tab["cont_struct"] = gf_cont_struct;
    subc_tab["cont_struct_get"] = gf_cont_struct_get;
    subc_tab["fem"] = gf_fem;
    subc_tab["fem_get"] = gf_fem_get;
    subc_tab["cvstruct_get"] = gf_cvstruct_get;
    subc_tab["mesher_object"] = gf_mesher_object;
    subc_tab["mesher_object_get"] = gf_mesher_object_get;
    subc_tab["mesh"] = gf_mesh;
    subc_tab["mesh_get"] = gf_mesh_get;
    subc_tab["mesh_set"] = gf_mesh_set;
    subc_tab["mesh_fem"] = gf_mesh_fem;
    subc_tab["mesh_fem_get"] = gf_mesh_fem_get;
    subc_tab["mesh_fem_set"] = gf_mesh_fem_set;
    subc_tab["mesh_im"] = gf_mesh_im;
    subc_tab["mesh_im_get"] = gf_mesh_im_get;
    subc_tab["mesh_im_set"] = gf_mesh_im_set;
    subc_tab["mesh_im_data"] = gf_mesh_im_data;
    subc_tab["mesh_im_data_get"] = gf_mesh_im_data_get;
    subc_tab["mesh_im_data_set"] = gf_mesh_im_data_set;
    subc_tab["model"] = gf_model;
    subc_tab["model_get"] = gf_model_get;
    subc_tab["model_set"] = gf_model_set;
    subc_tab["slice"] = gf_slice;
    subc_tab["slice_get"] = gf_slice_get;
    subc_tab["slice_set"] = gf_slice_set;
    subc_tab["levelset"] = gf_levelset;
    subc_tab["levelset_get"] = gf_levelset_get;
    subc_tab["levelset_set"] = gf_levelset_set;
    subc_tab["mesh_levelset"] = gf_mesh_levelset;
    subc_tab["mesh_levelset_get"] = gf_mesh_levelset_get;
    subc_tab["mesh_levelset_set"] = gf_mesh_levelset_set;
    subc_tab["asm"] = gf_asm;
    subc_tab["compute"] = gf_compute;
    subc_tab["precond"] = gf_precond;
    subc_tab["precond_get"] = gf_precond_get;
    subc_tab["spmat"] = gf_spmat;
    subc_tab["spmat_get"] = gf_spmat_get;
    subc_tab["spmat_set"] = gf_spmat_set;
    subc_tab["linsolve"] = gf_linsolve;
    subc_tab["util"] = gf_util;
    subc_tab["exit"] = gf_exit;
  }


  std::stringstream info;
  getfemint::global_pinfomsg = &info;
  *pinfomsg = 0; *pout_args = 0;
  // cout << "getfem_matlab_main(" << function << ", inarg=" << nb_in_args
  //      << ", outarg=" << *nb_out_args << ")\n";
  //  for (int i=0; i < nb_in_args; ++i) {
  //     cout << "  inarg[" << i << "]=";
  //     gfi_array_print((gfi_array*)in_args[i]); cout << "\n";
  //  }
  try {
    static getfemint::config *conf[3];
    if (!conf[config_id])
      conf[config_id] = new getfemint::config((gfi_interface_type)config_id);
    conf[config_id]->current_function_ = function;
    config::set_current_config(conf[config_id]);
    mexargs_in in(nb_in_args, in_args, false);
    mexargs_out out(*nb_out_args);
    out.set_scilab(bool(scilab_flag));

    SUBC_TAB::iterator it = subc_tab.find(function);
    if (it != subc_tab.end()) {
      it->second(in, out);
    }
    else {
      GMM_THROW(getfemint_bad_arg, "unknown function: " << function);
    }

    *pout_args = (gfi_array**)gfi_calloc(out.args().size(),sizeof(gfi_array*));
    if (!*pout_args) GMM_THROW(getfemint_error, "memory exhausted..");
    out.set_okay(1);
    *nb_out_args = int(out.args().size());
    std::copy(out.args().begin(), out.args().end(), *pout_args);
  }
  catch (const getfemint_bad_arg &e) {
    // cerr << "Bad argument!\n";
    return strdup(e.what());
  }
  catch (const getfemint_interrupted &) {
    cerr << "Ctrl-C catched!\n";
    return strdup("Interrupted (Ctrl-C)");
  }
  catch (const getfemint_error &e) {
    // cerr << "Error!\n";
    return strdup(e.what());
  }
  catch (const std::logic_error &e) {
    cerr << "logic_error exception caught\n";
    return strdup(e.what());
  }
  catch (const std::runtime_error &e) {
    cerr << "runtime_error exception caught\n";
    return strdup(e.what());
  }
  catch(const std::bad_alloc &) {
    return strdup("getfem caught a bad_alloc exception\n");
  }
  catch(const std::bad_cast &) {
    return strdup("getfem caught a bad_cast exception\n");
  }
  catch(const std::bad_typeid &) {
    return strdup("getfem caught a bad_typeid exception\n");
  }
  catch(const std::bad_exception &) {
    return strdup("getfem caught a bad_exception exception\n");
  }
  catch (...) {
    return strdup("getfem caught an unknown exception\n");
  }
  if (!info.str().empty()) {
    *pinfomsg = strdup(info.str().c_str());
  }
  //cout << "getfem_interface_main: exiting " << function << "\n";
  return 0;
}
