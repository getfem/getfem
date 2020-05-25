/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2006-2020 Julien Pommier.
 
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
#include <signal.h>
#include <string.h>
#include "gfm_common.h"

const char*
mxClassID2string(mxClassID id) {
  switch (id) {
    case mxUNKNOWN_CLASS: return "UNKNOWN";
    case mxCELL_CLASS:    return "CELL";
    case mxSTRUCT_CLASS:  return "STRUCT";
    case mxOBJECT_CLASS:  return "OBJECT";
    case mxCHAR_CLASS:    return "CHAR";
    case mxDOUBLE_CLASS:  return "DOUBLE";
    case mxSINGLE_CLASS:  return "SINGLE";
    case mxINT8_CLASS:    return "INT8";
    case mxUINT8_CLASS:   return "UINT8";
    case mxINT16_CLASS:   return "INT16";
    case mxUINT16_CLASS:  return "UINT16";
    case mxINT32_CLASS:   return "INT32";
    case mxUINT32_CLASS:  return "UINT32";
    case mxINT64_CLASS:   return "INT64";
    case mxUINT64_CLASS:  return "UINT64";
    case mxSPARSE_CLASS:  return "SPARSE";
#ifdef LOGICAL_IS_A_TYPE
    case mxLOGICAL_CLASS: return "LOGICAL";
#endif
    case mxOPAQUE_CLASS: return "OPAQUE_CLASS";
    default:
      if (id != (mxClassID)(-1))
	return "OBJECT";
      else
	return "unknown class...did you use the correct mex version ?";
  }
}

int
mxarray_to_gfi_array(const mxArray *mx, gfi_array *t)
{
  int n = mxGetNumberOfElements(mx);
  mwSize n2;
  n2=n;
  
  assert(t);
  switch (mxGetClassID(mx)) {    
    case mxCELL_CLASS: {
      mwSize i;
      t->storage.type = GFI_CELL;
      t->storage.gfi_storage_u.data_cell.data_cell_len = n;
      t->storage.gfi_storage_u.data_cell.data_cell_val = (gfi_array**)mxCalloc(n, sizeof(gfi_array*));
      for (i = 0; i < n; ++i) {
        if (mxGetCell(mx,i)) {
          t->storage.gfi_storage_u.data_cell.data_cell_val[i] = mxCalloc(1,sizeof(gfi_array));
          if (mxarray_to_gfi_array(mxGetCell(mx,i), t->storage.gfi_storage_u.data_cell.data_cell_val[i]) != 0) return 1;
        } else t->storage.gfi_storage_u.data_cell.data_cell_val[i] = NULL;
      }
    } break;
    case mxCHAR_CLASS: {
      t->storage.type = GFI_CHAR;
      t->storage.gfi_storage_u.data_char.data_char_len = n;    
      t->storage.gfi_storage_u.data_char.data_char_val = mxCalloc(n+1,sizeof(char));
      mxGetString(mx,t->storage.gfi_storage_u.data_char.data_char_val,n2+1);                           
      
    } break;
    case mxINT32_CLASS: {
      t->storage.type = GFI_INT32;
      t->storage.gfi_storage_u.data_int32.data_int32_len = n;
      t->storage.gfi_storage_u.data_int32.data_int32_val = mxGetData(mx);
    } break;
#ifdef LOGICAL_IS_A_TYPE
    case mxLOGICAL_CLASS: {
      int i;
      mxLogical *p = (mxLogical*)mxGetData(mx);
      t->storage.type = GFI_INT32;
      t->storage.gfi_storage_u.data_int32.data_int32_len = n;
      t->storage.gfi_storage_u.data_int32.data_int32_val = (int*)mxCalloc(n,sizeof(int));
      for (i = 0; i < n; ++i)
	t->storage.gfi_storage_u.data_int32.data_int32_val[i] = (p[i] ? 1 : 0);
    } break;
#endif
    case mxUINT32_CLASS: {
      t->storage.type = GFI_UINT32;
      t->storage.gfi_storage_u.data_uint32.data_uint32_len = n;
      t->storage.gfi_storage_u.data_uint32.data_uint32_val = mxGetData(mx);
    } break;
    case mxSPARSE_CLASS: /* for older (i.e. 6.1) matlab versions... */
    case mxDOUBLE_CLASS: {
      int is_complex = mxIsComplex(mx);
      
      if (!mxIsSparse(mx)) {
        t->storage.type = GFI_DOUBLE;
	t->storage.gfi_storage_u.data_double.is_complex = is_complex;
        t->storage.gfi_storage_u.data_double.data_double_len = n * (is_complex ? 2 : 1);
        if (!is_complex) {
          t->storage.gfi_storage_u.data_double.data_double_val = mxGetData(mx);
        } else {
          double *pr = mxGetData(mx), *pi = mxGetImagData(mx);
          int i;
          t->storage.gfi_storage_u.data_double.data_double_val = mxCalloc(n,sizeof(double)*2);
          for (i=0; i < n; ++i) { 
            t->storage.gfi_storage_u.data_double.data_double_val[i*2] = pr[i];
            t->storage.gfi_storage_u.data_double.data_double_val[i*2+1] = pi[i];
          }
        }
      } else {
        /*******************************************/
        int nnz = mxGetNzmax(mx);
        int i, lgs;
        mwIndex *ptr2;
        int *ptr;
        
        t->storage.type = GFI_SPARSE;
        t->storage.gfi_storage_u.sp.is_complex = is_complex;
        t->storage.gfi_storage_u.sp.ir.ir_len = nnz;
        ptr=mxCalloc(nnz,sizeof(int));
        ptr2=mxGetIr(mx);
        
        for (i=0; i < nnz; ++i)
        {
            ptr[i]=ptr2[i];
        }
        t->storage.gfi_storage_u.sp.ir.ir_val = ptr;
        /*t->storage.gfi_storage_u.sp.ir.ir_val = mxGetIr(mx);*/
        
        /*******************************************/
        lgs=mxGetN(mx)+1;
        t->storage.gfi_storage_u.sp.jc.jc_len = lgs;
        ptr=mxCalloc(lgs,sizeof(int));
        ptr2=mxGetJc(mx);
        
        for (i=0; i < lgs; ++i)
        {
            ptr[i]=ptr2[i];
        }
        
        t->storage.gfi_storage_u.sp.jc.jc_val = ptr;
        /*t->storage.gfi_storage_u.sp.jc.jc_val = mxGetJc(mx);*/  
        
        /********************************************/
        t->storage.gfi_storage_u.sp.pr.pr_len = nnz * (is_complex ? 2 : 1);
        if (!is_complex) {
          t->storage.gfi_storage_u.sp.pr.pr_val = mxGetPr(mx);        
        } else {
          double *pr = mxGetPr(mx), *pi = mxGetPi(mx);
          int i;
          t->storage.gfi_storage_u.sp.pr.pr_val = mxCalloc(nnz,sizeof(double)*2);
          for (i=0; i < nnz; ++i) {
            t->storage.gfi_storage_u.sp.pr.pr_val[i*2] = pr[i];
            t->storage.gfi_storage_u.sp.pr.pr_val[i*2+1] = pi[i];
          }
        }
      }
    } break;
    case mxOBJECT_CLASS:
    case mxSTRUCT_CLASS: 
    default: 
      {
      mxArray *fid = mxGetField(mx, 0, "id");
      mxArray *fcid = mxGetField(mx, 0, "cid");
      if (fid && mxGetClassID(fid) == mxUINT32_CLASS &&
          fcid && mxGetClassID(fcid) == mxUINT32_CLASS) {
        int n = mxGetNumberOfElements(fid),i;
        assert(n == mxGetNumberOfElements(fcid));
        if (mxGetNumberOfElements(fid) >= 1) {
          t->storage.type = GFI_OBJID;
          t->storage.gfi_storage_u.objid.objid_len = mxGetNumberOfElements(fid);
          t->storage.gfi_storage_u.objid.objid_val = mxCalloc(n, sizeof(gfi_object_id));
          for (i=0; i < n; ++i) {
            t->storage.gfi_storage_u.objid.objid_val[i].id = ((int*)mxGetData(fid))[i];
            t->storage.gfi_storage_u.objid.objid_val[i].cid = ((int*)mxGetData(fcid))[i];
          }
        } else {
          mexPrintf("empty arrays of getfem object ids not handled"); return 1;
        }
      } else { 
        mexPrintf("matlab structures (except getfem object ids) not handled"); return 1; 
      }
    } break;
  }
  {
    mwSize *pm;
    int k,kkk,dmlen;
    t->dim.dim_len = mxGetNumberOfDimensions(mx);
    dmlen=t->dim.dim_len;
    pm=(mwSize *)mxGetDimensions(mx);
    t->dim.dim_val = (u_int*)mxCalloc(dmlen,sizeof(int));
    
    for (k=0; k <dmlen; ++k) {   
      kkk=(int)pm[k];   
      t->dim.dim_val[k]=kkk;
    }
  }
  return 0;
}

mxArray*
gfi_array_to_mxarray(gfi_array *t) {
  mxArray *m;

  /* Matlab represent scalars as an array of size one */
  /* while gfi_array represents "scalar" values with 0-dimension array */
  int ndim = (t->dim.dim_len == 0 ? 1 : t->dim.dim_len);
  static const int one = 1;
  const int *dim = (t->dim.dim_len == 0 ? &one : (const int *)t->dim.dim_val);
  mwSize ndim2, *dim2;
  int number = 450;
  int ij;
  /****************** dimensions array in mwSize */  
  ndim2=ndim;
  dim2=mxCalloc(ndim,sizeof(mwSize));

  assert(t);
   
  for (ij=0;ij<ndim;ij++)
   {
        dim2[ij]=dim[ij];     
   }
  /*************** end dimensions array in mwSize */  

  switch (t->storage.type) {
    case GFI_UINT32: 
    case GFI_INT32: {
      int i;
      int *pr; 
      m = mxCreateNumericArray(ndim2, dim2, (t->storage.type == GFI_UINT32) ? mxUINT32_CLASS:mxINT32_CLASS, mxREAL);
      pr = mxGetData(m);
      for (i=0; i < t->storage.gfi_storage_u.data_int32.data_int32_len; ) 
      {   
        *pr++=(int)t->storage.gfi_storage_u.data_int32.data_int32_val[i++];      
      }
    } break;
    case GFI_DOUBLE: {
      int i;
      double *pr = 0, *pi = 0;
      if (!gfi_array_is_complex(t)) {
        m = mxCreateNumericArray(ndim2, dim2, mxDOUBLE_CLASS, mxREAL);
        pr=mxGetData(m);  
         for (i=0; i < t->storage.gfi_storage_u.data_double.data_double_len; ) 
         {
            *pr++=t->storage.gfi_storage_u.data_double.data_double_val[i++];
         }      
      } else {
        m = mxCreateNumericArray(ndim2, dim2, mxDOUBLE_CLASS, mxCOMPLEX);
        pr = mxGetData(m); pi = mxGetImagData(m);
        for (i=0; i < t->storage.gfi_storage_u.data_double.data_double_len; ) {
          *pr++ = t->storage.gfi_storage_u.data_double.data_double_val[i++];
          *pi++ = t->storage.gfi_storage_u.data_double.data_double_val[i++];
        }
      }
    } break;
    case GFI_CHAR: {
      char *s = calloc(t->storage.gfi_storage_u.data_char.data_char_len+1,1); 
      strncpy(s,t->storage.gfi_storage_u.data_char.data_char_val,t->storage.gfi_storage_u.data_char.data_char_len);
      m = mxCreateString(s); free(s);
    } break;
    case GFI_CELL: {
      unsigned i;
      m = mxCreateCellArray(ndim2, dim2);
      for (i=0; i < t->storage.gfi_storage_u.data_cell.data_cell_len; ++i)
        mxSetCell(m,i,gfi_array_to_mxarray(t->storage.gfi_storage_u.data_cell.data_cell_val[i]));
    } break;
    case GFI_OBJID: {
      unsigned i;
      mwSize j,i2;
      static const char *fields[] = {"id","cid"};
      mxArray *mxid, *mxcid;    
      j=t->storage.gfi_storage_u.objid.objid_len;
      m = mxCreateStructArray(1, (const mwSize*)&j, 2, fields);
      
      j=1;
      
      for (i=0; i < t->storage.gfi_storage_u.objid.objid_len; ++i) 
      {
        mxid = mxCreateNumericArray(1, (const mwSize*)&j, mxUINT32_CLASS, mxREAL);        
        *(int*)mxGetData(mxid) = t->storage.gfi_storage_u.objid.objid_val[i].id;       
        i2=i;       
        mxSetField(m,i2,fields[0], mxid);
        
        mxcid = mxCreateNumericArray(1, (const mwSize*)&j, mxUINT32_CLASS, mxREAL);                
        *(int*)mxGetData(mxcid) = t->storage.gfi_storage_u.objid.objid_val[i].cid;      
        mxSetField(m,i2,fields[1], mxcid);
      }
    } break;
    case GFI_SPARSE: {
        
     /****************** Ir treatement */
      int i;
      mwIndex d1,d2,nnz;
      mwIndex *pr2;
      assert(t);
      
      d1=t->dim.dim_val[0];
      d2=t->dim.dim_val[1];
      nnz=t->storage.gfi_storage_u.sp.pr.pr_len;
      m = mxCreateSparse(d1,d2,nnz,
                         gfi_array_is_complex(t) ? mxCOMPLEX : mxREAL);                 
      pr2=mxGetIr(m);
      
      for (i=0; i <  t->storage.gfi_storage_u.sp.ir.ir_len; i++) {
          *pr2 = t->storage.gfi_storage_u.sp.ir.ir_val[i];
          pr2++;
      }      
     /****************** end Ir treatement */
     /****************** Jc treatement */      
      pr2=mxGetJc(m);
      
        for (i=0; i <  t->storage.gfi_storage_u.sp.jc.jc_len; i++) {
          *pr2 =  t->storage.gfi_storage_u.sp.jc.jc_val[i];
          pr2++;
        }      
     /****************** end Jc treatement */
      
      if (!gfi_array_is_complex(t)) {
        memcpy(mxGetPr(m), t->storage.gfi_storage_u.sp.pr.pr_val, t->storage.gfi_storage_u.sp.pr.pr_len * sizeof(double));
      } else {
        double *pr, *pi; int i;
        pr = mxGetPr(m); pi = mxGetPi(m);
	for (i=0; i < t->storage.gfi_storage_u.sp.pr.pr_len; ) {
          *pr++ = t->storage.gfi_storage_u.sp.pr.pr_val[i++];
          *pi++ = t->storage.gfi_storage_u.sp.pr.pr_val[i++];
        }
      }
    } break;
    default:  {
      assert(0);
    } break;
  }
  mxFree(dim2); 
  return m;
}

/*******************************************************/

gfi_array_list *
build_gfi_array_list(int nrhs, const mxArray *prhs[]) {
  gfi_array_list *l;
  int i;
  gfi_array **pin;
  
  
  l = mxCalloc(1,sizeof(gfi_array_list));
  l->arg.arg_len = nrhs;
  l->arg.arg_val = mxCalloc(nrhs, sizeof(gfi_array));
  for (i=0; i < nrhs; ++i) {   
    if (mxarray_to_gfi_array(prhs[i], &l->arg.arg_val[i]) != 0) return NULL;
  }
  
  /*-------------*/
  pin = gfi_calloc(l->arg.arg_len, sizeof(gfi_array*));
  for (i=0; i < l->arg.arg_len; ++i) {
    pin[i] = &l->arg.arg_val[i];       
  }
  
  return l;
}


#ifndef WIN32
struct sigaction old_sigint;
#endif


static int sigint_hit = 0;
static getfem_sigint_handler_t sigint_callback;

static void sigint(int sig) {
  sigint_callback(sig);
  remove_custom_sigint(0);
  sigint_hit++;
}


void install_custom_sigint(getfem_sigint_handler_t h) {
#ifndef WIN32 /* matlab on win does not use signals so.. */
  struct sigaction new_sigint;
  new_sigint.sa_handler = sigint;
  sigint_callback = h;
  sigemptyset(&new_sigint.sa_mask);
  new_sigint.sa_flags = 0;
  sigaction (SIGINT, NULL, &old_sigint);
  if (old_sigint.sa_handler != SIG_IGN)
    sigaction(SIGINT, &new_sigint, NULL);
  sigint_hit = 0;
#endif
}

void remove_custom_sigint(int allow_rethrow) {
#ifndef WIN32
  struct sigaction act;
  sigaction (SIGINT, NULL, &act);
  if (act.sa_handler == sigint) {
    sigaction(SIGINT, &old_sigint, NULL);
  }
  if (allow_rethrow && sigint_hit) {
    fprintf(stderr, "ready, raising SIGINT now\n");
    raise(SIGINT); 
  }
  sigint_hit = 0; 
#endif
}
