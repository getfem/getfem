/*===========================================================================

 Copyright (C) 2009-2020 Yann Collette

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

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <time.h>

#ifndef _MSC_VER
#include <unistd.h>
#endif
extern "C" {
#include <api_scilab.h> 
#include <sciprint.h>
#include <Scierror.h>
}

extern "C" {
#include "gfm_common.h"
}

#include "gfi_array.h"
#include "getfem_interface.h"

#ifndef _MSCVER
#include "stream_redirect.h"
#endif

//#define DEBUG_TIMER
//#define DEBUG
//#define DEBUG2

extern "C" int handle_getfem_callback();
extern "C" void set_superlu_callback(int (*cb)());

gfi_output * call_getfem_interface(char *funname, gfi_array_list in, int nlhs)
{
  static gfi_output result;
  gfi_array **pin  = NULL;
  gfi_array **pout = NULL;
  char *errmsg = 0, *infomsg = 0;
  int i;

#ifdef DEBUG2
  sciprint("call_getfem_interface: len = %d \n",in.arg.arg_len);
  for(i=0;i<in.arg.arg_len;i++)
    {
      sciprint("element %d: \n",i);
      gfi_array_print(&in.arg.arg_val[i]);
    }
#endif

  pin = (gfi_array **)gfi_calloc(in.arg.arg_len, sizeof(gfi_array*));
  for (i=0; i < in.arg.arg_len; ++i) 
    {
      pin[i] = &in.arg.arg_val[i];
    }

#ifdef DEBUG
  sciprint("call_getfem_interface: funname = %s len = %d\n",funname, in.arg.arg_len);
#endif

  errmsg = getfem_interface_main(SCILAB_INTERFACE, funname, in.arg.arg_len, (const gfi_array **)pin, &nlhs, &pout, &infomsg, 1);

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
      result.gfi_output_u.output.arg.arg_val = (gfi_array *)gfi_calloc(nlhs, sizeof(gfi_array));
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
  const char * s = current_scilab_function; if (!s) s = "doh!!";
  fprintf(stderr, "*** CTRL-C hit during execution of the getfem_scilab function: gf_%s...\n" \
	  "You will gain control as soon as the current operation is finished ***\n" \
	  "If you want to abort immediatly the current operation, hit CTRL-C again\n" \
	  "In that case, you will have to restart getfem_scilab:\n", s);
  set_cancel_flag(1);
  assert(handle_getfem_callback() == 1);
}

extern "C" int sci_gf_scilab(char * fname) 
{
  gfi_output     * out  = NULL;
  gfi_array_list * in   = NULL;
  gfi_array_list * outl = NULL;
  int ** ptr_param = NULL;
  int sci_x;
  unsigned int i;
  SciErr _SciErr;
#ifdef DEBUG_TIMER
  clock_t time_start, time_end;
#endif
#ifndef _MSCVER
  ScilabStream scicout(std::cout);
  ScilabStream scicerr(std::cerr);
#endif
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
      _SciErr = getVarAddressFromPosition(pvApiCtx,i,ptr_param+i);
#ifdef DEBUG
      _SciErr = getVarDimension(pvApiCtx,ptr_param[i],&pirow,&picol);
      _SciErr = getVarType(pvApiCtx,ptr_param[i],&var_type);
      sciprint("sci_gf_scilab: position %d - address %d - type %d - dimension %d %d\n", i, ptr_param[i],var_type,pirow,picol);
#endif
    }

#ifdef DEBUG
      sciprint("sci_gf_scilab: list of parameters built\n", i, Rhs);
#endif

#ifdef DEBUG_TIMER
  sciprint("DEBUG_TIMER: before build_gfi_array_list\n");
  time_start = clock()/CLOCKS_PER_SEC;
#endif
    
  in = build_gfi_array_list(Rhs, ptr_param);

#ifdef DEBUG_TIMER
  time_end = clock()/CLOCKS_PER_SEC;
  sciprint("DEBUG_TIMER: after build_gfi_array_list: %f\n", (double)(time_end - time_start));
  time_start = time_end;
#endif

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

  out = call_getfem_interface(&fname[3], *in, Lhs);

  // We now remove all the allocated memory for the input parameters
  // This memory is allocated in gfm_common.c
  if (in) 
    {
      for(i=0;i<in->arg.arg_len;i++)
	gfi_array_destroy(&in->arg.arg_val[i]);
      gfi_free(in->arg.arg_val);
      gfi_free(in); 
    }

#ifdef DEBUG_TIMER
  time_end = clock()/CLOCKS_PER_SEC;
  sciprint("DEBUG_TIMER: after call_getfem_interface: %f\n",(double)(time_end - time_start));
  time_start = time_end;
#endif

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

#ifdef DEBUG_TIMER
  time_end = clock()/CLOCKS_PER_SEC;
  sciprint("DEBUG_TIMER: after gfi_array_to_sci_array: %f\n", (double)(time_end - time_start));
#endif

  if (ptr_param) FREE(ptr_param);

  return 0;
}
