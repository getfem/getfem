/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2006-2012 Yves Renard, Julien Pommier.
 
 This file is a part of GetFEM++
 
 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
/* #include <unistd.h> */
#include "octave/mex.h"
#include "gfm_common.h"
#include "getfem_interface.h"

void set_superlu_callback(int (*cb)());
int handle_getfem_callback();

/* main file for the giant gf_matlab mex-file */
/*
char* getfem_interface_main(int config_id, const char *function, 
			    int nb_in_args,
			    const gfi_array *in_args[], 
			    int *nb_out_args,
			    gfi_array ***pout_args, char **pinfomsg, bool scilab_flag);
*/

#define FREE(p) { /*printf("%s@%d", __FILE__, __LINE__); */gfi_free(p); }
gfi_output *
call_getfem_interface(char *funname, gfi_array_list in, int nlhs)
{
  static gfi_output  result;
  /*  result.status = GFI_STATUS_OK;
  result.gfi_output_u.output.arg.arg_len = nlhs;
  result.gfi_output_u.output.arg.arg_val = gfi_calloc(nlhs, sizeof(gfi_array));
  return &result;
  */
  unsigned int i;
  gfi_array **pin;
  gfi_array **pout = NULL;
  char *errmsg=0, *infomsg=0;
  /*printf("** CALLED FUNCTION gf_%s(", funname);*/

  pin = gfi_calloc(in.arg.arg_len, sizeof(gfi_array*));
  for (i=0; i < in.arg.arg_len; ++i) {
    pin[i] = &in.arg.arg_val[i];
  }
  errmsg = getfem_interface_main(MATLAB_INTERFACE, funname, in.arg.arg_len, (const gfi_array **)pin, &nlhs, &pout, &infomsg,0);
  result.infomsg = infomsg;
  if (errmsg) {
    result.status = GFI_STATUS_ERROR;
    result.gfi_output_u.errmsg = errmsg;
  } else {
    result.status = GFI_STATUS_OK;
    result.gfi_output_u.output.arg.arg_len = nlhs;
    result.gfi_output_u.output.arg.arg_val = gfi_calloc(nlhs, sizeof(gfi_array));
    for (i=0; i < nlhs; ++i) {
      assert(pout[i]);
      result.gfi_output_u.output.arg.arg_val[i] = *pout[i];
      FREE(pout[i]);
    }
    if (pout) FREE(pout);
  }
  FREE(pin);
  return &result;
}

char *current_matlab_function = NULL;

void sigint_callback(int sig) {
  char *s = current_matlab_function; if (!s) s = "doh!!";
  fprintf(stderr, "*** CTRL-C hit during execution of the getfem_octave function: gf_%s...\n" \
	  "You will gain control as soon as the current operation is finished ***\n" \
	  "If you want to abort immediatly the current operation, hit CTRL-C again\n" \
	  "In that case, you will have to restart getfem_octave, using 'clear functions' for example:\n", s);
  set_cancel_flag(1);
  assert(handle_getfem_callback() == 1);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  gfi_output *out;
  gfi_array_list *in;

  set_cancel_flag(0);
  set_superlu_callback(is_cancel_flag_set);

  if (nrhs == 0 || !mxIsChar(prhs[0])) {
    mexErrMsgTxt("missing function name");
  }
  current_matlab_function = mxCalloc(mxGetNumberOfElements(prhs[0])+1,1);
  if (mxGetString(prhs[0], current_matlab_function, mxGetNumberOfElements(prhs[0])+1))
    mexErrMsgTxt("mxGetString failure");
  in = build_gfi_array_list(nrhs-1, prhs+1);
  if (!in) mexErrMsgTxt("a problem occured while reading arguments");


  install_custom_sigint(sigint_callback);
  out = call_getfem_interface(current_matlab_function, *in, nlhs);
  remove_custom_sigint(out->status == GFI_STATUS_OK);

  if (out == NULL) {
    mexPrintf("could not connect to getfem_octave server...");
  } else {
    if (out->infomsg) {
      mexPrintf("message from [gf_%s]:\n%s\n", current_matlab_function, out->infomsg);
    }

    if (out->status == GFI_STATUS_OK) {
      gfi_array_list *outl = &out->gfi_output_u.output;
      unsigned i;
      if (nlhs == 0 && outl->arg.arg_len > 0) { /* not very nice */
	mxArray *mx = gfi_array_to_mxarray(&outl->arg.arg_val[0]);
	mexPutVariable("caller", "ans", mx);
      }
      for (i=0; i < outl->arg.arg_len; ++i) {
        plhs[i] = gfi_array_to_mxarray(&outl->arg.arg_val[i]);
        gfi_array_destroy(&outl->arg.arg_val[i]);
      }
      FREE(outl->arg.arg_val);
    } else {
      /* duplicate the string into a octave-allocated buffer.. */
      char *s = mxCalloc(strlen(out->gfi_output_u.errmsg)+1,1);
      strcpy(s, out->gfi_output_u.errmsg);
      FREE(out->gfi_output_u.errmsg);
      mexErrMsgTxt(s);
    }
  }
}

