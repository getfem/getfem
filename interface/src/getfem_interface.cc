// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2008 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

#include "getfem_interface.h"
#include "getfemint.h"

using namespace getfemint;

void gf_workspace(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_delete(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_eltm(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_geotrans(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_geotrans_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_integ(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_integ_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_fem(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_fem_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_cvstruct_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_fem(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_fem_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_fem_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_im(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_im_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mesh_im_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mdbrick(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mdbrick_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mdbrick_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mdstate(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mdstate_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
void gf_mdstate_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out);
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
    default:
      THROW_INTERNAL_ERROR;
    }
  }
}


extern "C" 
char* getfem_interface_main(int config_id, const char *function, 
			    int nb_in_args,
			    const gfi_array *in_args[], 
			    int *nb_out_args,
			    gfi_array ***pout_args, char **pinfomsg)
{
  std::stringstream info;
  getfemint::global_pinfomsg = &info;
  *pinfomsg = 0; *pout_args = 0;
  //cout << "getfem_matlab_main(" << function << ", inarg=" << nb_in_args << ", outarg=" << *nb_out_args << ")\n";
  //for (int i=0; i < nb_in_args; ++i) { cout << "  inarg[" << i << "]=";gfi_array_print((gfi_array*)in_args[i]); cout << "\n"; }
  try {
    static getfemint::config *conf[2];  
    if (!conf[config_id]) conf[config_id] = new getfemint::config((gfi_interface_type)config_id);
    conf[config_id]->current_function_ = function;
    config::set_current_config(conf[config_id]);
    mexargs_in in(nb_in_args, in_args, false);
    mexargs_out out(*nb_out_args);    
    if (strcmp(function, "workspace")==0) gf_workspace(in,out);
    else if (strcmp(function, "delete")==0) gf_delete(in,out);
    else if (strcmp(function, "eltm")==0) gf_eltm(in,out);
    else if (strcmp(function, "geotrans")==0) gf_geotrans(in,out);
    else if (strcmp(function, "geotrans_get")==0) gf_geotrans_get(in,out);
    else if (strcmp(function, "integ")==0) gf_integ(in,out);
    else if (strcmp(function, "integ_get")==0) gf_integ_get(in,out);
    else if (strcmp(function, "fem")==0) gf_fem(in,out);
    else if (strcmp(function, "fem_get")==0) gf_fem_get(in,out);
    else if (strcmp(function, "cvstruct_get")==0) gf_cvstruct_get(in,out);
    else if (strcmp(function, "mesh")==0) gf_mesh(in,out);
    else if (strcmp(function, "mesh_get")==0) gf_mesh_get(in,out);
    else if (strcmp(function, "mesh_set")==0) gf_mesh_set(in,out);
    else if (strcmp(function, "mesh_fem")==0) gf_mesh_fem(in,out);
    else if (strcmp(function, "mesh_fem_get")==0) gf_mesh_fem_get(in,out);
    else if (strcmp(function, "mesh_fem_set")==0) gf_mesh_fem_set(in,out);
    else if (strcmp(function, "mesh_im")==0) gf_mesh_im(in,out);
    else if (strcmp(function, "mesh_im_get")==0) gf_mesh_im_get(in,out);
    else if (strcmp(function, "mesh_im_set")==0) gf_mesh_im_set(in,out);
    else if (strcmp(function, "mdbrick")==0) gf_mdbrick(in,out);
    else if (strcmp(function, "mdbrick_get")==0) gf_mdbrick_get(in,out);
    else if (strcmp(function, "mdbrick_set")==0) gf_mdbrick_set(in,out);
    else if (strcmp(function, "mdstate")==0) gf_mdstate(in,out);
    else if (strcmp(function, "mdstate_get")==0) gf_mdstate_get(in,out);
    else if (strcmp(function, "mdstate_set")==0) gf_mdstate_set(in,out);
    else if (strcmp(function, "slice")==0) gf_slice(in,out);
    else if (strcmp(function, "slice_get")==0) gf_slice_get(in,out);
    else if (strcmp(function, "slice_set")==0) gf_slice_set(in,out);
    else if (strcmp(function, "levelset")==0) gf_levelset(in,out);
    else if (strcmp(function, "levelset_get")==0) gf_levelset_get(in,out);
    else if (strcmp(function, "levelset_set")==0) gf_levelset_set(in,out);
    else if (strcmp(function, "mesh_levelset")==0) gf_mesh_levelset(in,out);
    else if (strcmp(function, "mesh_levelset_get")==0) gf_mesh_levelset_get(in,out);
    else if (strcmp(function, "mesh_levelset_set")==0) gf_mesh_levelset_set(in,out);
    else if (strcmp(function, "asm")==0) gf_asm(in,out);
    else if (strcmp(function, "compute")==0) gf_compute(in,out);
    else if (strcmp(function, "precond")==0) gf_precond(in,out);
    else if (strcmp(function, "precond_get")==0) gf_precond_get(in,out);
    else if (strcmp(function, "spmat")==0) gf_spmat(in,out);
    else if (strcmp(function, "spmat_get")==0) gf_spmat_get(in,out);
    else if (strcmp(function, "spmat_set")==0) gf_spmat_set(in,out);
    else if (strcmp(function, "linsolve")==0) gf_linsolve(in,out);
    else if (strcmp(function, "util")==0) gf_util(in,out);
    else if (strcmp(function, "exit")==0) exit(0);
    else {
      GMM_THROW(getfemint_bad_arg, "unknown function: " << function);
    }
    
    *pout_args = (gfi_array**)gfi_calloc(out.args().size(), sizeof(gfi_array*));
    if (!*pout_args) GMM_THROW(getfemint_error, "memory exhausted..");
    out.set_okay(1);
    *nb_out_args = int(out.args().size());
    std::copy(out.args().begin(), out.args().end(), *pout_args);
  }
  catch (getfemint_bad_arg e) {
    //cerr << "Bad argument!\n";
    return strdup(e.what());
  }
  catch (getfemint_interrupted) {
    cerr << "Ctrl-C catched!\n";
    return strdup("Interrupted (Ctrl-C)");
  }
  catch (getfemint_error e) {
    return strdup(e.what());
  }
  catch (std::logic_error e) {
    cerr << "logic_error exception caught\n";
    return strdup(e.what());
  }
  catch (std::runtime_error e) {
    cerr << "runtime_error exception caught\n";
    return strdup(e.what());
  }
  catch(std::bad_alloc) {
    return strdup("getfem caught a bad_alloc exception\n");
  }
  catch(std::bad_cast) {
    return strdup("getfem caught a bad_cast exception\n");
  }
  catch(std::bad_typeid) {
    return strdup("getfem caught a bad_typeid exception\n");
  }
  catch(std::bad_exception) {
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
