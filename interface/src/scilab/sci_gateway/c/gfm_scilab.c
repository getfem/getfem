/* -*- c++ -*- (enables emacs c++ mode) */
/*========================================================================

 Copyright (C) 2009 Yann Collette

 This file is a part of GETFEM++

 Getfem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 You should have received a copy of the GNU Lesser General Public
 License along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
 USA.

 ========================================================================*/

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <stdio.h>

#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <api_common.h>
#include <MALLOC.h>

#include "gfm_common.h"
#include "getfem_interface.h"

#define DEBUG

#define FREE(p) { /*printf("%s@%d", __FILE__, __LINE__); */gfi_free(p); }

gfi_output * call_getfem_interface(char *funname, gfi_array_list in, int nlhs)
{
  static gfi_output result;
  gfi_array **pin = NULL;
  gfi_array **pout = NULL;
  char *errmsg=0, *infomsg=0;
  unsigned int i;

#ifdef DEBUG
  sciprint("call_getfem_interface: len = %d \n",in.arg.arg_len);
  for(i=0;i<in.arg.arg_len;i++)
    {
      sciprint("element %d: \n",i);
      gfi_array_print(&in.arg.arg_val[i]);
    }
#endif

  pin = gfi_calloc(in.arg.arg_len, sizeof(gfi_array*));
  for (i=0; i < in.arg.arg_len; ++i) 
    {
      pin[i] = &in.arg.arg_val[i];
    }

#ifdef DEBUG
  sciprint("call_getfem_interface: funname = %s len = %d\n",funname, in.arg.arg_len);
#endif

  errmsg = getfem_interface_main(SCILAB_INTERFACE, funname, in.arg.arg_len, (const gfi_array **)pin, &nlhs, &pout, &infomsg);

#ifdef DEBUG
  sciprint("call_getfem_interface: end of getfem_interface_main, nlhs = %d infomsg = %s\n",nlhs,infomsg);
#endif

  result.infomsg = infomsg;
  if (errmsg) 
    {
      result.status = GFI_STATUS_ERROR;
      result.gfi_output_u.errmsg = errmsg;
    } 
  else 
    {
      result.status = GFI_STATUS_OK;
      result.gfi_output_u.output.arg.arg_len = nlhs;
      result.gfi_output_u.output.arg.arg_val = gfi_calloc(nlhs, sizeof(gfi_array));
      for(i=0; i<nlhs; ++i) 
	{
	  assert(pout[i]);
	  result.gfi_output_u.output.arg.arg_val[i] = *pout[i];
	  gfi_free(pout[i]);
	}
      if (pout) gfi_free(pout);
    }
  gfi_free(pin);

  return &result;
}

char *current_scilab_function = NULL;

void sigint_callback(int sig) 
{
  char *s = current_scilab_function; if (!s) s = "doh!!";
  fprintf(stderr, "*** CTRL-C hit during execution of the getfem_scilab function: gf_%s...\n" \
	  "You will gain control as soon as the current operation is finished ***\n" \
	  "If you want to abort immediatly the current operation, hit CTRL-C again\n" \
	  "In that case, you will have to restart getfem_matlab, using 'clear functions' for example:\n", s);
  set_cancel_flag(1);
  assert(handle_getfem_callback() == 1);
}

int sci_gf_scilab(char * fname) 
{
  gfi_output     * out  = NULL;
  gfi_array_list * in   = NULL;
  gfi_array_list * outl = NULL;
  int ** ptr_param = NULL;
  int *  sci_x     = NULL;
  int picol, pirow;
  unsigned int i;

  set_cancel_flag(0);
  set_superlu_callback(is_cancel_flag_set);
  
#ifdef DEBUG
  sciprint("sci_gf_scilab: Rhs = %d Lhs = %d\n", Rhs, Lhs);
#endif

  ptr_param = (int **)MALLOC((Rhs+1)*sizeof(int *));
  for(i=1;i<=Rhs;i++)
    {
#ifdef DEBUG
      sciprint("sci_gf_scilab: i = %d Rhs = %d\n", i, Rhs);
#endif
      getVarAddressFromPosition(i,ptr_param+i);
#ifdef DEBUG
      getVarDimension(ptr_param[i],&pirow,&picol);
      sciprint("sci_gf_scilab: position %d - address %d - type %d - dimension %d %d\n", i, ptr_param[i],getVarType(ptr_param[i]),pirow,picol);
#endif
    }

#ifdef DEBUG
      sciprint("sci_gf_scilab: list of parameters built\n", i, Rhs);
#endif

  in = build_gfi_array_list(Rhs, ptr_param);
  if (!in) 
    {
      Scierror(999,"%s: a problem occured while reading arguments.\n",fname);
      return 0;
    }
#ifdef DEBUG
  sciprint("sci_gf_scilab: end of build_gfi_array_list\n");
#endif

  install_custom_sigint(sigint_callback);

#ifdef DEBUG
  sciprint("sci_gf_scilab: fname = %s - calling call_getfem_interface with %s\n",fname, &fname[3]);
#endif

  out = call_getfem_interface(&fname[3], *in, Lhs); // Sans parametre de sorti ou avec un parametre de sortie, Lhs = 1 !!!

#ifdef DEBUG
  sciprint("sci_gf_scilab: end of call_getfem_interface\n");
#endif

  remove_custom_sigint(out->status == GFI_STATUS_OK);
  
  if (out == NULL) 
    {
      sciprint("%s: could not connect to getfem_scilab server...\n",fname);
      LhsVar(1) = 0;
    } 
  else 
    {
      if (out->infomsg) 
	{
	  sciprint("%s: message:\n%s\n",fname,out->infomsg);
	}
      
      if (out->status == GFI_STATUS_OK) 
	{
	  outl = &out->gfi_output_u.output;
	  for(i=0; i<outl->arg.arg_len; ++i) 
	    {
#ifdef DEBUG
	      sciprint("sci_gf_scilab: processing output argument %d / %d\n", i, outl->arg.arg_len);
	      sciprint("storage type = %d\n", outl->arg.arg_val[i].storage.type);
	      sciprint("output position: LhsVar(%d) = %d\n", i+1, Rhs+1+i);
#endif
	      sci_x = gfi_array_to_sci_array(&outl->arg.arg_val[i],Rhs+1+i);
	      LhsVar(i+1) = Rhs+1+i;
	      if (&outl->arg.arg_val[i]) gfi_array_destroy(&outl->arg.arg_val[i]);
	    }
	  
	  gfi_free(outl->arg.arg_val);
	} 
      else 
	{
	  Scierror(999,"%s: %s\n", fname,out->gfi_output_u.errmsg);
	  LhsVar(1) = 0;
	}
    }

  FREE(ptr_param);

  return 0;
}
