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
#include <signal.h>
#include <string.h>
#include <stdio.h>

#include <stack-c.h>
#include <Scierror.h>
#include <sciprint.h>
#include <MALLOC.h>
#include <api_common.h>
#include <api_string.h>
#include <api_double.h>
#include <api_int.h>
#include <api_sparse.h>
#include <api_list.h>
#include <api_boolean.h>

#include "gfm_common.h"

//#define DEBUG

const char * sci_ClassID2string(sci_types id) 
{
  switch (id) 
    {
    case sci_matrix:             return "MATRIX";
    case sci_poly:               return "POLY";
    case sci_boolean:            return "BOOLEAN";
    case sci_boolean_sparse:     return "BOOLEAN_SPARSE";
    case sci_matlab_sparse:      return "MATLAB_SPARSE";
    case sci_ints:               return "INTS";
    case sci_handles:            return "HANDLES";
    case sci_strings:            return "STRINGS";
    case sci_u_function:         return "U_FUNCTION";
    case sci_c_function:         return "C_FUNCTION";
    case sci_lib:                return "LIB";
    case sci_list:               return "LIST";
    case sci_mlist:              return "MLIST";
    case sci_tlist:              return "TLIST";
    case sci_lufact_pointer:     return "LUFACT_POINTER";
    case sci_implicit_poly:      return "IMPLICIT_POLY";
    case sci_intrinsic_function: return "INTRINSIC_FUNCTION";
    default:
      return "unknown class: did you use the correct scilab version ?";
  }
}

int sci_array_to_gfi_array(int * sci_x, gfi_array *t)
{
  int i, n = 0;

#ifdef DEBUG
  int _picol, _pirow;
  getVarDimension(sci_x,&_pirow,&_picol);
  sciprint("sci_array_to_gfi_array: dimension of current variable %d %d - type = %d\n",_pirow,_picol,getVarType(sci_x));
#endif

  assert(t);

  switch (getVarType(sci_x)) 
    {
    case sci_list: 
      {
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: dealing with sci_list\n");
#endif
	int i;

	getListItemNumber(sci_x,&n);

	t->storage.type = GFI_CELL;
	t->storage.gfi_storage_u.data_cell.data_cell_len = n;
	t->storage.gfi_storage_u.data_cell.data_cell_val = (gfi_array**)MALLOC(n*sizeof(gfi_array*));

	for(i=0; i<n; ++i) 
	  {
	    t->storage.gfi_storage_u.data_cell.data_cell_val[i] = MALLOC(1*sizeof(gfi_array));
	    int * item_address = NULL;
	    getListItemAddress(sci_x,i+1,&item_address);
#ifdef DEBUG
	    sciprint("type of item %d: %d\n", i+1, getVarType(item_address));
#endif
	    if (sci_array_to_gfi_array(item_address, t->storage.gfi_storage_u.data_cell.data_cell_val[i]) != 0) return 1;
	  }

	t->dim.dim_len = 1;
	t->dim.dim_val = (u_int*)MALLOC(1*sizeof(u_int));
	t->dim.dim_val[0] = n;
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: end\n");
#endif
      } 
      break;
    case sci_mlist: 
      {
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: dealing with mlist\n");
#endif
	int pirow, picol, *pilen, *pilistaddress;
	char ** pstStrings;
	double * pdblDataID, * pdblDataCID;
	int * pintDims;
	double * pdblDims;

	getListItemNumber(sci_x,&n);
	getMatrixOfStringInList(sci_x,1,&pirow,&picol,NULL,NULL);
	pilen = (int *)MALLOC(pirow*picol*sizeof(int));
	getMatrixOfStringInList(sci_x,1,&pirow,&picol,pilen,NULL);
	pstStrings = (char **)MALLOC(pirow*picol*sizeof(char *));
	for(i=0;i<pirow*picol;i++)
	  {
	    pstStrings[i] = (char *)MALLOC((pilen[i]+1)*sizeof(char));
	  }
	getMatrixOfStringInList(sci_x,1,&pirow,&picol,pilen,pstStrings);

#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: pstStrings[0] = %s\n",pstStrings[0]);
#endif

	if (strcmp(pstStrings[0],"objid")==0)
	  {
#ifdef DEBUG
	    sciprint("sci_array_to_gfi_array: dealing with objid\n");
#endif
	    t->storage.type = GFI_OBJID;

	    getMatrixOfDoubleInList(sci_x, 2, &pirow, &picol, &pdblDataID);
	    getMatrixOfDoubleInList(sci_x, 3, &pirow, &picol, &pdblDataCID);

#ifdef DEBUG
	    sciprint("sci_array_to_gfi_array: pirow = %d picol = %d\n", pirow, picol);
#endif

	    n = pirow*picol;

	    t->storage.gfi_storage_u.objid.objid_len = n;
	    t->storage.gfi_storage_u.objid.objid_val = (struct gfi_object_id *)MALLOC(n*sizeof(struct gfi_object_id));

	    for(i=0;i<n;i++)
	      {
		t->storage.gfi_storage_u.objid.objid_val[i].id  = (int)pdblDataID[i];
		t->storage.gfi_storage_u.objid.objid_val[i].cid = (int)pdblDataCID[i];
#ifdef DEBUG
		sciprint("sci_array_to_gfi_array: objid[%d]: id = %d cid = %d\n", i, (int)pdblDataID[i], (int)pdblDataCID[i]);
#endif
	      }

	    t->dim.dim_len = 1;
	    t->dim.dim_val = (u_int*)MALLOC(1*sizeof(u_int));
	    t->dim.dim_val[0] = n;
	  }
	else if (strcmp(pstStrings[0],"hm")==0)
	  {
	    int piPrecision, ret, * p_item_address;
#ifdef DEBUG
	    sciprint("sci_array_to_gfi_array: dealing with hypermat\n");
#endif

	    // Get the dimensions
	    getListItemAddress(sci_x,2,&p_item_address);
	    switch(getVarType(p_item_address))
	      {
	      case sci_matrix:
		ret = getMatrixOfDoubleInList(sci_x, 2, &pirow, &picol, &pdblDims);
		t->dim.dim_len = pirow*picol;
		t->dim.dim_val = (u_int*)MALLOC(pirow*picol*sizeof(u_int));
		for(i=0;i<pirow*picol;i++) t->dim.dim_val[i] = (u_int)pdblDims[i];
		break;
	      case sci_ints:
		ret = getMatrixOfInteger32InList(sci_x, 2, &pirow, &picol, &pintDims);
		t->dim.dim_len = pirow*picol;
		t->dim.dim_val = (u_int*)MALLOC(pirow*picol*sizeof(u_int));
		for(i=0;i<pirow*picol;i++) t->dim.dim_val[i] = (u_int)pintDims[i];
		break;
	      default:
		Scierror(999,"wrong type for hypermatrix dimensions type: %d\n", getVarType(p_item_address));
		return 1;
	      }
	    // Get the matrixes (stored as a column vector of size prod(size(...)))
	    // We must detect if we have a INT UINT or DOUBLE
	    
	    getListItemAddress(sci_x, 3, &pilistaddress);

	    if (getVarType(pilistaddress)==sci_matrix)
	      {
		if (isVarComplex(pilistaddress))
		  {
		    t->storage.type = GFI_DOUBLE;
		    getComplexMatrixOfDoubleInList(sci_x, 3, &pirow, &picol, &pdblDataID, &pdblDataCID);

		    t->storage.gfi_storage_u.data_double.is_complex = GFI_COMPLEX;
		    
		    t->storage.gfi_storage_u.data_double.data_double_len = 2*pirow*picol;
		    t->storage.gfi_storage_u.data_double.data_double_val = (double *)MALLOC(2*pirow*picol*sizeof(double));
		    for(i=0;i<pirow*picol;++i) 
		      {
			t->storage.gfi_storage_u.data_double.data_double_val[2*i+0] = pdblDataID[i];
			t->storage.gfi_storage_u.data_double.data_double_val[2*i+1] = pdblDataCID[i];
		      }
		  }
		else
		  {
		    t->storage.type = GFI_DOUBLE;
		    getMatrixOfDoubleInList(sci_x, 3, &pirow, &picol, &pdblDataID);
		    
		    t->storage.gfi_storage_u.data_double.is_complex = GFI_REAL;

		    t->storage.gfi_storage_u.data_double.data_double_len = pirow*picol;
		    t->storage.gfi_storage_u.data_double.data_double_val = (double *)MALLOC(pirow*picol*sizeof(double));
		    for(i=0;i<pirow*picol;++i) t->storage.gfi_storage_u.data_double.data_double_val[i] = pdblDataID[i];
		  }
	      }
	    else if (getVarType(pilistaddress)==sci_ints)
	      {
		getMatrixOfIntegerPrecision(sci_x,&piPrecision);
		if ((piPrecision!=SCI_INT32)&&(piPrecision!=SCI_UINT32))
		  {
		    Scierror(999,"Can deal only with int32 or uint32\n");
		    return 1;
		  }
		
		if (piPrecision==SCI_INT32)
		  {
		    int * piData32;

		    t->storage.type = GFI_INT32;
		    getMatrixOfInteger32(sci_x,&pirow,&picol,&piData32);
		    
		    t->storage.gfi_storage_u.data_int32.data_int32_len = pirow*picol;
		    t->storage.gfi_storage_u.data_int32.data_int32_val = (int *)MALLOC(pirow*picol*sizeof(int));
		    for(i=0;i<pirow*picol;++i) t->storage.gfi_storage_u.data_int32.data_int32_val[i] = piData32[i];
		  }
		else if (piPrecision==SCI_UINT32)
		  {
		    unsigned int * puiData32;

		    t->storage.type = GFI_UINT32;
		    getMatrixOfUnsignedInteger32(sci_x,&pirow,&picol,&puiData32);
		    
		    t->storage.gfi_storage_u.data_uint32.data_uint32_len = pirow*picol;
		    t->storage.gfi_storage_u.data_uint32.data_uint32_val = (unsigned int *)MALLOC(pirow*picol*sizeof(unsigned int));
		    for(i=0;i<pirow*picol;++i) t->storage.gfi_storage_u.data_uint32.data_uint32_val[i] = puiData32[i];
		  }
		else
		  {
		    Scierror(999,"Can deal only with int32 or uint32\n");
		    return 1;
		  }
	      }
	    else
	      {
		Scierror(999,"Can deal only with double, int32 or uint32\n");
		return 1;
	      }
	  }
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: end\n");
#endif
      } 
      break;
    case sci_strings: 
      {
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: dealing with sci_strings\n");
#endif
	int picol, pirow, iret;
	int * pilen = NULL;
	char ** pstData = NULL;

	t->storage.type = GFI_CHAR;

	// First call to get picol and pirow
	getMatrixOfString(sci_x, &pirow, &picol, NULL, NULL);
	if ((pirow!=1)&&(picol!=1))
	  {
	    Scierror(999,"Can allocate only one string at a time\n");
	    return 1;
	  }
	pilen = (int *)MALLOC(pirow*picol*sizeof(int));

	// Second call to get pilen
	getMatrixOfString(sci_x, &pirow, &picol, pilen, NULL);
	pstData = (char **)MALLOC(pirow*picol*sizeof(char*));
	for(i=0; i<pirow*picol ; i++) pstData[i] = (char *)MALLOC((pilen[i] + 1)*sizeof(char));

	n = pilen[0] + 1;
	t->storage.gfi_storage_u.data_char.data_char_len = n;
	t->storage.gfi_storage_u.data_char.data_char_val = MALLOC((n)*sizeof(char));

	// Third call to retrieve data
	iret = getMatrixOfString(sci_x, &pirow, &picol, pilen, pstData);
	memcpy(t->storage.gfi_storage_u.data_char.data_char_val,pstData[0],(pilen[0] + 1)*sizeof(char));
#ifdef DEBUG
	sciprint("pirow = %d picol = %d pilen = %d\n", pirow, picol, pilen[0]);
	sciprint("storing |%s|\n",t->storage.gfi_storage_u.data_char.data_char_val);
#endif
	t->dim.dim_len = 1;
	t->dim.dim_val = (u_int*)MALLOC(1*sizeof(u_int));
	t->dim.dim_val[0] = pilen[0];

	FREE(pilen);
	for(i=0; i<pirow*picol ; i++) FREE(pstData[i]);
	FREE(pstData);

#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: end\n");
#endif
      } 
      break;
    case sci_ints: 
      {
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: dealing with sci_ints\n");
#endif
	int piPrecision;

	getMatrixOfIntegerPrecision(sci_x,&piPrecision);
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: precision = %d\n",piPrecision);
#endif
	if ((piPrecision!=SCI_INT32)&&(piPrecision!=SCI_UINT32))
	  {
	    Scierror(999,"Can deal only with int32 or uint32\n");
	    return 1;
	  }
	
	if (piPrecision==SCI_INT32)
	  {
	    t->storage.type = GFI_INT32;
	    int pirow, picol;
	    int * piData32;
	    getMatrixOfInteger32(sci_x,&pirow,&picol,&piData32);
#ifdef DEBUG
	    sciprint("DEBUG: %d dimensions - dim[0] = %d, dim[1] = %d\n", n, pirow, picol);
#endif
	    n = picol*pirow;
	    t->storage.gfi_storage_u.data_int32.data_int32_len = n;
	    t->storage.gfi_storage_u.data_int32.data_int32_val = (int *)MALLOC(n*sizeof(int));
	    for(i=0;i<n;++i)
	      {
		t->storage.gfi_storage_u.data_int32.data_int32_val[i] = piData32[i];
	      }

	    t->dim.dim_len = 2;
	    t->dim.dim_val = (u_int*)MALLOC(2*sizeof(u_int));
	    t->dim.dim_val[0]= pirow;
	    t->dim.dim_val[1]= picol;
	  }
	else
	  {
	    t->storage.type = GFI_UINT32;
	    int pirow, picol;
	    unsigned int * puiData32;
	    getMatrixOfUnsignedInteger32(sci_x,&pirow,&picol,&puiData32);

	    n = picol*pirow;
	    t->storage.gfi_storage_u.data_uint32.data_uint32_len = n;
	    t->storage.gfi_storage_u.data_uint32.data_uint32_val = (unsigned int *)MALLOC(n*sizeof(unsigned int));
	    for(i=0;i<n;++i)
	      {
		t->storage.gfi_storage_u.data_uint32.data_uint32_val[i] = puiData32[i];
	      }

	    t->dim.dim_len = 2;
	    t->dim.dim_val = (u_int*)MALLOC(2*sizeof(u_int));
	    t->dim.dim_val[0]= pirow;
	    t->dim.dim_val[1]= picol;
	  }
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: end\n");
#endif
      } 
      break;
    case sci_boolean: 
      {
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: dealing with sci_boolean\n");
#endif
	int pirow, picol;
	int * piBool;

	t->storage.type = GFI_INT32;

	getMatrixOfBoolean(sci_x,&pirow,&picol,&piBool);

	n = picol*pirow;

	t->storage.gfi_storage_u.data_int32.data_int32_len = n;
	t->storage.gfi_storage_u.data_int32.data_int32_val = (int *)MALLOC(n*sizeof(int));
	for(i=0;i<n;++i)
	  {
	    t->storage.gfi_storage_u.data_int32.data_int32_val[i] = piBool[i];
	  }

	t->dim.dim_len = 2;
	t->dim.dim_val = (u_int*)MALLOC(2*sizeof(u_int));
	t->dim.dim_val[0]= pirow;
	t->dim.dim_val[1]= picol;
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: end\n");
#endif
      } 
      break;
    case sci_matrix: 
      {
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: dealing with sci_matrix\n");
#endif
	int is_complex = isVarComplex(sci_x);
	int pirow, picol;
	double * pdblDataReal, * pdblDataImag;
	
	t->storage.type = GFI_DOUBLE;
	t->storage.gfi_storage_u.data_double.is_complex = is_complex;
	
	if (!is_complex) 
	  {
	    getMatrixOfDouble(sci_x,&pirow,&picol,&pdblDataReal);
	    
	    n = pirow*picol;

	    t->storage.gfi_storage_u.data_double.is_complex = GFI_REAL;
	    t->storage.gfi_storage_u.data_double.data_double_len = n;
	    t->storage.gfi_storage_u.data_double.data_double_val = (double *)MALLOC(n*sizeof(double));
	    for(i=0;i<n;++i)
	      {
		t->storage.gfi_storage_u.data_double.data_double_val[i] = pdblDataReal[i];
	      }
	  } 
	else 
	  {
	    getComplexMatrixOfDouble(sci_x,&pirow,&picol,&pdblDataReal,&pdblDataImag);

	    n = pirow*picol;

	    t->storage.gfi_storage_u.data_double.is_complex = GFI_COMPLEX;
	    t->storage.gfi_storage_u.data_double.data_double_len = 2*n;
	    t->storage.gfi_storage_u.data_double.data_double_val = (double *)MALLOC(2*n*sizeof(double));

	    for(i=0;i<n;++i) 
	      { 
		t->storage.gfi_storage_u.data_double.data_double_val[i*2]   = pdblDataReal[i];
		t->storage.gfi_storage_u.data_double.data_double_val[i*2+1] = pdblDataImag[i];
	      }
	  }

	t->dim.dim_len = 2;
	t->dim.dim_val = (u_int*)MALLOC(2*sizeof(u_int));
	t->dim.dim_val[0]= pirow;
	t->dim.dim_val[1]= picol;
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: end\n");
#endif
      }
      break;
    case sci_sparse:
      {
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: dealing with sci_sparse\n");
#endif
	int is_complex = isVarComplex(sci_x);
	int pirow, picol, nbitem;
	double * pdblDataReal, * pdblDataImag;
	int * nbitemrow, * picolpos, * offset;
	int j, k, Index;

        t->storage.type = GFI_SPARSE;
        t->storage.gfi_storage_u.sp.is_complex = is_complex;

	if (!is_complex) 
	  {
#ifdef DEBUG
	    sciprint("sci_array_to_gfi_array: not complex\n");
#endif
	    getSparseMatrix(sci_x,&pirow,&picol,&nbitem,&nbitemrow,&picolpos,&pdblDataReal);
	    
#ifdef DEBUG
	    for(i=0;i<pirow;i++) sciprint("nbitemrow[%d] = %d\n", i, nbitemrow[i]);
	    for(i=0;i<nbitem;i++) sciprint("picolpos[%d] = %d, pdblDataReal[%d] = %f\n",i,picolpos[i],i,pdblDataReal[i]);
#endif

	    t->storage.gfi_storage_u.sp.is_complex = GFI_REAL;

	    t->storage.gfi_storage_u.sp.ir.ir_len = nbitem;
	    t->storage.gfi_storage_u.sp.jc.jc_len = picol+1;
	    t->storage.gfi_storage_u.sp.pr.pr_len = nbitem;

	    t->storage.gfi_storage_u.sp.ir.ir_val = (int *)MALLOC(nbitem*sizeof(int));
	    t->storage.gfi_storage_u.sp.jc.jc_val = (int *)MALLOC((picol+1)*sizeof(int));
	    t->storage.gfi_storage_u.sp.pr.pr_val = (double *)MALLOC(nbitem*sizeof(double));

	    offset = (int *)MALLOC(pirow*sizeof(int));

	    // We compute position offset
	    offset[0] = 0;
	    for(i=1;i<pirow;i++)
	      {
		offset[i] = offset[i-1] + nbitemrow[i-1];
#ifdef DEBUG
		sciprint("offset[%d] = %d\n", i, offset[i]);
#endif
	      }

	    Index = 0;
	    for(i=0;i<picol;i++)
	      {
		t->storage.gfi_storage_u.sp.jc.jc_val[i] = Index;
		for(j=0;j<pirow;j++)
		  {
		    for(k=0;k<nbitemrow[j];k++)
		      {
			if (i==picolpos[offset[j]+k]-1)
			  {
			    t->storage.gfi_storage_u.sp.ir.ir_val[Index] = j;
			    t->storage.gfi_storage_u.sp.pr.pr_val[Index] = pdblDataReal[offset[j]+k];
#ifdef DEBUG
			    sciprint("Index = %d offset[%d] = %d, k = %d, pdblDataReal[%d] = %f\n", Index, j, offset[j], k, offset[j]+k, pdblDataReal[offset[j]+k]);
#endif
			    Index++;
			  }
		      }
		  }
	      }
	    t->storage.gfi_storage_u.sp.jc.jc_val[picol] = nbitem;
	  }
	else
	  {
#ifdef DEBUG
	    sciprint("sci_array_to_gfi_array: complex\n");
#endif
	    getComplexSparseMatrix(sci_x,&pirow,&picol,&nbitem,&nbitemrow,&picolpos,&pdblDataReal,&pdblDataImag);
	    
	    // We store the transposed matrix in t
	    t->storage.gfi_storage_u.sp.is_complex = GFI_COMPLEX;

	    t->storage.gfi_storage_u.sp.ir.ir_len = nbitem;
	    t->storage.gfi_storage_u.sp.jc.jc_len = picol+1;
	    t->storage.gfi_storage_u.sp.pr.pr_len = 2*nbitem;

	    t->storage.gfi_storage_u.sp.ir.ir_val = (int *)MALLOC(nbitem*sizeof(int));
	    t->storage.gfi_storage_u.sp.jc.jc_val = (int *)MALLOC((picol+1)*sizeof(int));
	    t->storage.gfi_storage_u.sp.pr.pr_val = (double *)MALLOC(2*nbitem*sizeof(double));

	    offset = (int *)MALLOC(pirow*sizeof(int));

	    // We compute position offset
	    offset[0] = 0;
	    for(i=1;i<pirow;i++)
	      offset[i] = offset[i-1] + nbitemrow[i-1];

	    Index = 0;
	    for(i=0;i<picol;i++)
	      {
		t->storage.gfi_storage_u.sp.jc.jc_val[i] = Index;
		for(j=0;j<pirow;j++)
		  {
		    for(k=0;k<nbitemrow[j];k++)
		      {
			if (i==picolpos[offset[j]+k]-1)
			  {
			    t->storage.gfi_storage_u.sp.ir.ir_val[Index]     = j;
			    t->storage.gfi_storage_u.sp.pr.pr_val[2*Index+0] = pdblDataReal[offset[j]+k];
			    t->storage.gfi_storage_u.sp.pr.pr_val[2*Index+1] = pdblDataImag[offset[j]+k];
			    Index++;
			  }
		      }
		  }
	      }
	    t->storage.gfi_storage_u.sp.jc.jc_val[picol] = nbitem;
	  }

#ifdef DEBUG
// 	for(i=0;i<t->storage.gfi_storage_u.sp.jc.jc_len;i++)
// 	  {
// 	    sciprint("The line %d starts at index %d. length = %d\n",i, t->storage.gfi_storage_u.sp.jc.jc_val[i], t->storage.gfi_storage_u.sp.jc.jc_len);
// 	    for(j=0;j<(t->storage.gfi_storage_u.sp.jc.jc_val[i+1]-t->storage.gfi_storage_u.sp.jc.jc_val[i]);j++)
// 	      {
// 		Index = t->storage.gfi_storage_u.sp.jc.jc_val[i];
// 		sciprint("[%f,%f] ",t->storage.gfi_storage_u.sp.pr.pr_val[2*(Index+j)+0],t->storage.gfi_storage_u.sp.pr.pr_val[2*(Index+j)+1]);
// 	      }
// 	    sciprint("\n");
// 	  }

	sciprint("sci_array_to_gfi_array: pirow = %d picol = %d\n",pirow, picol);
	sciprint("sci_array_to_gfi_array: end\n");
#endif
	t->dim.dim_len = 2;
	t->dim.dim_val = (u_int*)MALLOC(2*sizeof(u_int));
	t->dim.dim_val[0]= pirow;
	t->dim.dim_val[1]= picol;
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: pirow = %d picol = %d\n",pirow, picol);
	sciprint("sci_array_to_gfi_array: end\n");
#endif
	FREE(offset);
      }
      break;
    default: 
      {
	Scierror(999,"unhandled class type : %s\n", sci_ClassID2string(getVarType(sci_x)));
	return 1;
      } 
      break;
    }

  return 0;
}

int * gfi_array_to_sci_array(gfi_array *t, int ivar) 
{
  // Revoir cette partie. Surtout la partie GFI_CELL ...
  // Ajouter des fonctions
  // - exportString pour les variables simples
  // - addCellToCell, addStringToCell, etc.. pour gérer les listes imbriquées

  int *m_var;

  assert(t);

  /* Scilab represent scalars as an array of size one */
  /* while gfi_array represents "scalar" values with 0-dimension array */

  int ndim = (t->dim.dim_len == 0 ? 1 : t->dim.dim_len);
  static const int one = 1;
  const int *dim = (t->dim.dim_len == 0 ? &one : (const int *)t->dim.dim_val);
  
  switch (t->storage.type) 
    {
    case GFI_UINT32: 
      {
	int is_hypermat = 0;
	int nrow, ncol, i;
	int pirow, picol;

	if (ndim==1)
	  {
	    nrow = dim[0];
	    ncol = 1;
	  }
	else if (ndim==2)
	  {
	    nrow = dim[0];
	    ncol = dim[1];
	  }
	else
	  {
	    is_hypermat = 1;
	  }

#ifdef DEBUG
	sciprint("gfi_array_to_sci_array: create from a GFI_UINT32\n");
	sciprint("ndim = %d ivar = %d\n", ndim,ivar);
#endif
	if (~is_hypermat)
	  {
	    createMatrixOfUnsignedInteger32(ivar,dim[0],dim[1],t->storage.gfi_storage_u.data_uint32.data_uint32_val);
	    getVarAddressFromPosition(ivar, &m_var);
	  }
	else
	  {
	    char *fields[] = {"hm","dims","entries"};
	    int * dims;
	    unsigned int * entries;
	    int nb_elem = 1;

	    createMList(ivar,3,&m_var);
	    
	    createMatrixOfStringInList(ivar, m_var, 1, 1, 3, fields);
	    
	    dims = (int *)MALLOC(t->dim.dim_len*sizeof(int));
	    for(i=0;i<t->dim.dim_len;i++) 
	      {
		dims[i] = (int)t->dim.dim_val[i];
		nb_elem *= (int)t->dim.dim_val[i];
	      }

	    entries = (unsigned int *)MALLOC(nb_elem*sizeof(unsigned int));
	    for(i=0;i<nb_elem;i++) entries[i] = t->storage.gfi_storage_u.data_uint32.data_uint32_val[i];
	    
	    // Add a vector to the 'dims' field
	    createMatrixOfInteger32InList(ivar, m_var, 2, 1, t->dim.dim_len, dims);
	    // Add a vector to the 'entries' field
	    createMatrixOfUnsignedInteger32InList(ivar, m_var, 3, 1, nb_elem, entries);
		
	    FREE(dims);
	    FREE(entries);
	  }
      } 
      break;
    case GFI_INT32: 
      {
	int is_hypermat = 0;
	int nrow, ncol, i;
	int pirow, picol;

	if (ndim==1)
	  {
	    nrow = (int)dim[0];
	    ncol = 1;
	  }
	else if (ndim==2)
	  {
	    nrow = (int)dim[0];
	    ncol = (int)dim[1];
	  }
	else
	  {
	    is_hypermat = 1;
	  }

#ifdef DEBUG
	sciprint("gfi_array_to_sci_array: create from a GFI_INT32\n");
	sciprint("ndim = %d ivar = %d dim[0] = %d dim[1] = %d\n", ndim,ivar, (int)dim[0],(int)dim[1]);
#endif
	if (~is_hypermat)
	  {
	    createMatrixOfInteger32(ivar,dim[0],dim[1],t->storage.gfi_storage_u.data_int32.data_int32_val);
	    getVarAddressFromPosition(ivar, &m_var);
	  }
	else
	  {
	    char *fields[] = {"hm","dims","entries"};
	    int * dims;
	    double * entries;
	    int nb_elem = 1;

	    createMList(ivar,3,&m_var);
	    
	    createMatrixOfStringInList(ivar, m_var, 1, 1, 3, fields);
	    
	    dims = (int *)MALLOC(t->dim.dim_len*sizeof(int));
	    for(i=0;i<t->dim.dim_len;i++) 
	      {
		dims[i] = (int)t->dim.dim_val[i];
		nb_elem *= (int)t->dim.dim_val[i];
	      }

	    entries = (double *)MALLOC(nb_elem*sizeof(double));
	    for(i=0;i<nb_elem;i++) entries[i] = (double)t->storage.gfi_storage_u.data_int32.data_int32_val[i];
	    
	    // Add a vector to the 'dims' field
	    createMatrixOfInteger32InList(ivar, m_var, 2, 1, t->dim.dim_len, dims);
	    // Add a vector to the 'entries' field
	    createMatrixOfDoubleInList(ivar, m_var, 3, 1, nb_elem, entries);
		
	    FREE(entries);
	    FREE(dims);
	  }
      } 
      break;
    case GFI_DOUBLE: 
      {
#ifdef DEBUG
	sciprint("gfi_array_to_sci_array: create from a GFI_DOUBLE\n");
	sciprint("ndim = %d ivar = %d\n", ndim, ivar);
	sciprint("dim[0] = %d\n", dim[0]);
#endif
	int i, nrow, ncol;
	double *pr, *pi;
	int is_hypermat = 0;

	if (ndim==1)
	  {
	    nrow = dim[0];
	    ncol = 1;
	  }
	else if (ndim==2)
	  {
	    nrow = dim[0];
	    ncol = dim[1];
	  }
	else
	  {
	    is_hypermat = 1;
	  }

#ifdef DEBUG
	sciprint("DEBUG: hypermat = %d\n",is_hypermat);
#endif
	if (!is_hypermat)
	  {
	    if (!gfi_array_is_complex(t)) 
	      {
#ifdef DEBUG
		sciprint("DEBUG: array is not complex\n");
#endif
		createMatrixOfDouble(ivar, nrow, ncol, t->storage.gfi_storage_u.data_double.data_double_val);
	      } 
	    else 
	      {
#ifdef DEBUG
		sciprint("DEBUG: array is complex\n");
#endif
		double *pr, *pi; int i;
		
		pr = (double *)MALLOC(nrow*ncol*sizeof(double));
		pi = (double *)MALLOC(nrow*ncol*sizeof(double));
		
		for(i=0;i<nrow*ncol;i++)
		  {
		    pr[i] = t->storage.gfi_storage_u.data_double.data_double_val[2*i];
		    pi[i] = t->storage.gfi_storage_u.data_double.data_double_val[2*i+1];
		  }
		
		createComplexMatrixOfDouble(ivar, nrow, ncol, pr, pi);
		
		FREE(pr);
		FREE(pi);
	      }
	  }
	else
	  {
#ifdef DEBUG
	    sciprint("DEBUG: array is hypermat\n");
#endif
	    char *fields[] = {"hm","dims","entries"};
	    int * dims;
	    int nb_elem = 1;

	    createMList(ivar,3,&m_var);
	    
	    createMatrixOfStringInList(ivar, m_var, 1, 1, 3, fields);
	    
	    dims = (int *)MALLOC(t->dim.dim_len*sizeof(int));
	    for(i=0;i<t->dim.dim_len;i++) 
	      {
		dims[i] = (int)t->dim.dim_val[i];
		nb_elem *= (int)t->dim.dim_val[i];
	      }

	    if (!gfi_array_is_complex(t)) 
	      {
		double * entries;

		entries = (double *)MALLOC(nb_elem*sizeof(double));
		for(i=0;i<nb_elem;i++) entries[i] = t->storage.gfi_storage_u.data_double.data_double_val[i];
		
		// Add a vector to the 'dims' field
		createMatrixOfInteger32InList(ivar, m_var, 2, 1, t->dim.dim_len, dims);
		// Add a vector to the 'entries' field
		createMatrixOfDoubleInList(ivar, m_var, 3, 1, nb_elem, entries);
		
		FREE(entries);
#ifdef DEBUG
		sciprint("DEBUG: end array is hypermat\n");
#endif
	      }
	    else
	      {
		double * entries_pr, *entries_pi;

		entries_pr = (double *)MALLOC(nb_elem*sizeof(double));
		entries_pi = (double *)MALLOC(nb_elem*sizeof(double));
		for(i=0;i<nb_elem;i++) 
		  {
		    entries_pr[i] = t->storage.gfi_storage_u.data_double.data_double_val[2*i+0];
		    entries_pi[i] = t->storage.gfi_storage_u.data_double.data_double_val[2*i+1];
		  }
		
		// Add a vector to the 'dims' field
		createMatrixOfInteger32InList(ivar, m_var, 2, 1, t->dim.dim_len, dims);
		// Add a vector to the 'entries' field
		createComplexMatrixOfDoubleInList(ivar, m_var, 3, 1, nb_elem, entries_pr, entries_pi);
		
		FREE(entries_pr);
		FREE(entries_pi);
	      }
	    FREE(dims);
	  }

	getVarAddressFromPosition(ivar, &m_var);
      } 

      break;
    case GFI_CHAR: 
      {
#ifdef DEBUG
	sciprint("gfi_array_to_sci_array: create from a GFI_CHAR\n");
	sciprint("ndim = %d ivar = %d\n", dim[0],ivar);
#endif
	char * tmp_string = (char *)MALLOC((t->storage.gfi_storage_u.data_char.data_char_len+1)*sizeof(char));
	memcpy(tmp_string,t->storage.gfi_storage_u.data_char.data_char_val,t->storage.gfi_storage_u.data_char.data_char_len*sizeof(char));
	tmp_string[t->storage.gfi_storage_u.data_char.data_char_len] = '\0';
	createMatrixOfString(ivar, 1, 1, &tmp_string);
	getVarAddressFromPosition(ivar, &m_var);
#ifdef DEBUG
	sciprint("ivar = %d string = |%s| len = %d\n",ivar, tmp_string,t->storage.gfi_storage_u.data_char.data_char_len);
#endif
	if (tmp_string) FREE(tmp_string);
      } 
      break;
    case GFI_CELL: 
      {
	unsigned int i, nb_item;
	int * m_content;
#ifdef DEBUG
	sciprint("gfi_array_to_sci_array: create from a GFI_CELL\n");
	sciprint("ndim = %d ivar = %d\n", dim[0],ivar);

	for(i=0;i<ndim;i++)
	  sciprint("dim[%d] = %d\n",i,dim[i]);

	sciprint("now create list at pos %d, dimension %d\n", ivar, dim[0]);
#endif
	createList(ivar,dim[0],&m_var);

	for(i=0; i<t->storage.gfi_storage_u.data_cell.data_cell_len; ++i)
	  {
	    m_content = gfi_array_to_sci_array(t->storage.gfi_storage_u.data_cell.data_cell_val[i],ivar+i+1);

	    switch(getVarType(m_content))
	      {
	      case sci_list:
		{
#ifdef DEBUG
	sciprint("gfi_array_to_sci_array: create from a GFI_CELL - sci_list\n");
#endif
		  getListItemNumber(m_content,&nb_item);
		  createListInList(ivar, m_content, i+1, nb_item, &m_var);
		}
		break;
	      case sci_mlist:
		{
#ifdef DEBUG
	sciprint("gfi_array_to_sci_array: create from a GFI_CELL - sci_mlist\n");
#endif
		  getListItemNumber(m_content,&nb_item);
		  createMListInList(ivar, m_content, i+1, nb_item, &m_var);
		}
		break;
	      case sci_strings:
		{
#ifdef DEBUG
		  sciprint("gfi_array_to_sci_array: create from a GFI_CELL - sci_strings\n");
		  sciprint("add element %d to position %d\n", i, ivar);
#endif
		  int pirow, picol, *pilen, j;
		  char ** pstStrings;
		  
		  // Get the matrix of strings from the gfi_array
		  getMatrixOfString(m_content, &pirow, &picol, NULL, NULL);
		  pilen = (int *)MALLOC(pirow*picol*sizeof(int));
		  pstStrings = (char **)MALLOC(pirow*picol*sizeof(char *));
		  getMatrixOfString(m_content, &pirow, &picol, pilen, NULL);
		  for(j=0;j<pirow*picol;j++)
		    {
		      pstStrings[j] = (char *)MALLOC((pilen[j]+1)*sizeof(char));
		    }
		  getMatrixOfString(m_content, &pirow, &picol, pilen, pstStrings);
#ifdef DEBUG
		  sciprint("pirow = %d picol = %d, pilen[0] = %d\n", pirow, picol, pilen[0]);
#endif
		  // And now add it to the list
		  createMatrixOfStringInList(ivar, m_var, i+1, pirow, picol, pstStrings);
		  
		  // Desallocate
		  if (pilen) FREE(pilen);
		  for(j=0;j<pirow*picol;j++)
		    {
		      if (pstStrings[j]) FREE(pstStrings[j]);
		    }
		  if (pstStrings) FREE(pstStrings);
		}
		break;
	      case sci_ints:
#ifdef DEBUG
		sciprint("gfi_array_to_sci_array: create from a GFI_CELL - sci_ints\n");
#endif
		// UINT32 + INT32
		{
		  int pi_precision, pirow, picol;
		  int * piData32;
		  unsigned int * puiData32;
		  
		  getMatrixOfIntegerPrecision(m_content,&pi_precision);
		  if ((pi_precision!=SCI_INT32)&&(pi_precision!=SCI_UINT32))
		    {
		      Scierror(999,"Can deal only with int32 or uint32\n");
		    }
		  switch(pi_precision)
		    {
		    case SCI_INT32:
		      getMatrixOfInteger32(m_content, &pirow, &picol, &piData32);
		      createMatrixOfInteger32InList(ivar, m_var, i+1, pirow, picol, piData32);
		      break;
		    case SCI_UINT32:
		      getMatrixOfUnsignedInteger32(m_content, &pirow, &picol, &puiData32);
		      createMatrixOfUnsignedInteger32InList(ivar, m_var, i+1, pirow, picol, puiData32);
		      break;
		    default:
		      Scierror(999,"Can deal only with int32 or uint32\n");
		    }
		}
		break;
	      case sci_boolean:
		{
#ifdef DEBUG
		  sciprint("gfi_array_to_sci_array: create from a GFI_CELL - sci_boolean\n");
#endif
		  int pirow, picol;
		  int * piBool;
		  getMatrixOfBoolean(m_content, &pirow, &picol, &piBool);
		  
		  createMatrixOfBooleanInList(ivar, m_var, i+1, pirow, picol, piBool);
		}
		break;
	      case sci_matrix:
		{
#ifdef DEBUG
		  sciprint("gfi_array_to_sci_array: create from a GFI_CELL - sci_matrix\n");
#endif
		  int pirow, picol;
		  double * pdblReal, * pdblImag;
		  
		  if (isVarComplex(m_content))
		    {
		      getComplexMatrixOfDouble(m_content, &pirow, &picol, &pdblReal, &pdblImag);
		      createComplexMatrixOfDoubleInList(ivar, m_var, i+1, pirow, picol, pdblReal, pdblImag);
		    }
		  else
		    {
		      getMatrixOfDouble(m_content, &pirow, &picol, &pdblReal);
		      createMatrixOfDoubleInList(ivar, m_var, i+1, pirow, picol, pdblReal);
		    }
		}
		break;
	      case sci_sparse:
		{
#ifdef DEBUG
		  sciprint("gfi_array_to_sci_array: create from a GFI_CELL - sci_sparse\n");
#endif
		  int pirow, picol, nbitem, * nb_item_row, * pi_col_pos;
		  double * pdblReal, * pdblImag;
		  
		  if (isVarComplex(m_content))
		    {
		      getComplexSparseMatrix(m_content, &pirow, &picol, &nbitem, &nb_item_row, &pi_col_pos, &pdblReal, &pdblImag);
		      createComplexSparseMatrixInList(ivar, m_var, i+1, pirow, picol, nbitem, nb_item_row, pi_col_pos, pdblReal, pdblImag);
		    }
		  else
		    {
		      getSparseMatrix(m_content, &pirow, &picol, &nbitem, &nb_item_row, &pi_col_pos, &pdblReal);
		      createSparseMatrixInList(ivar, m_var, i+1, pirow, picol, nbitem, nb_item_row, pi_col_pos, pdblReal);
		    }
		}
		break;
	      }
	  }
      } 
      break;
    case GFI_OBJID: 
      {
	unsigned i, size_objid;
	double * pdblDataID, * pdblDataCID;
	char *fields[] = {"objid","id","cid"};

	createMList(ivar,3,&m_var);

	createMatrixOfStringInList(ivar, m_var, 1, 1, 3, fields);

	size_objid  = t->storage.gfi_storage_u.objid.objid_len;
	pdblDataID  = (double *)MALLOC(size_objid*sizeof(double));
	pdblDataCID = (double *)MALLOC(size_objid*sizeof(double));

	for(i=0; i<size_objid; ++i) 
	  {
	    pdblDataID[i]  = t->storage.gfi_storage_u.objid.objid_val[i].id;       
	    pdblDataCID[i] = t->storage.gfi_storage_u.objid.objid_val[i].cid;
	  }
	// Add a vector to the 'id' field
	createMatrixOfDoubleInList(ivar, m_var, 2, 1, size_objid, pdblDataID);

	// Add a vector to the 'cid' field
	createMatrixOfDoubleInList(ivar, m_var, 3, 1, size_objid, pdblDataCID);

	if (pdblDataID)  FREE(pdblDataID);
	if (pdblDataCID) FREE(pdblDataCID);
      } 
      break;
    case GFI_SPARSE: 
      {
	int picol, pirow, nbitem, * pi_col_pos, * nb_item_row;
	int i, j, k, nnz, Index;
	double * pdblDataReal, * pdblDataImag;
	int iscomplex = gfi_array_is_complex(t);

	nbitem = t->storage.gfi_storage_u.sp.ir.ir_len;
	pirow  = t->dim.dim_val[0];
	picol  = t->dim.dim_val[1];

	// Convert from Matlab to Scilab format
	nb_item_row  = (int *)MALLOC(pirow*sizeof(int));
	pi_col_pos   = (int *)MALLOC(nbitem*sizeof(int));
	pdblDataReal = (double *)MALLOC(nbitem*sizeof(double));
	if (iscomplex)
	  pdblDataImag = (double *)MALLOC(nbitem*sizeof(double));

#ifdef DEBUG
	sciprint("picol = %d pirow = %d nbitem = %d\n", picol, pirow, nbitem);

	for(i=0;i<picol+1;i++)
	  sciprint("jc_val[%d] = %d\n",i,t->storage.gfi_storage_u.sp.jc.jc_val[i]);
	for(i=0;i<nbitem;i++)
	  {
	    sciprint("ir_val[%d] = %d pr_val[%d] = %f\n",i,t->storage.gfi_storage_u.sp.ir.ir_val[i],i,t->storage.gfi_storage_u.sp.pr.pr_val[i]);
	  }
#endif

	Index = 0;
	if (iscomplex)
	  {
	    for(i=0;i<pirow;i++)
	      {
		nb_item_row[i] = 0;
		for(j=0;j<nbitem;j++)
		  {
		    if (t->storage.gfi_storage_u.sp.ir.ir_val[j]==i)
		      {
			// We found the row number, we need now to find the corresponding col number.
			for(k=0;k<pirow;k++)
			  {
			    if (t->storage.gfi_storage_u.sp.jc.jc_val[k]>j)
			      {
				pi_col_pos[Index] = k;
				break;
			      }
			  }
			
			pdblDataReal[Index] = t->storage.gfi_storage_u.sp.pr.pr_val[2*j+0];
			pdblDataImag[Index] = t->storage.gfi_storage_u.sp.pr.pr_val[2*j+1];
			nb_item_row[i]++;
			Index++;
		      }
		  }
	      }
	  }
	else
	  {
	    for(i=0;i<pirow;i++)
	      {
		nb_item_row[i] = 0;
		for(j=0;j<nbitem;j++)
		  {
		    if (t->storage.gfi_storage_u.sp.ir.ir_val[j]==i)
		      {
			// We found the row number, we need now to find the corresponding col number.
			for(k=0;k<pirow;k++)
			  {
			    if (t->storage.gfi_storage_u.sp.jc.jc_val[k]>j)
			      {
				pi_col_pos[Index] = k;
				break;
			      }
			  }
			
			pdblDataReal[Index] = t->storage.gfi_storage_u.sp.pr.pr_val[j];
			nb_item_row[i]++;
			Index++;
		      }
		  }
	      }
	  }

#ifdef DEBUG
	for(i=0;i<pirow;i++)
	  {
	    sciprint("nb_item_row[%d] = %d\n", i, nb_item_row[i]);
	  }
	for(i=0;i<nbitem;i++)
	  {
	    sciprint("pi_col_pos[%d] = %d val = %f\n", i, pi_col_pos[i],pdblDataReal[i]);
	  }
#endif

	if (iscomplex)
	  {
	    createComplexSparseMatrix(ivar, pirow, picol, nbitem, nb_item_row, pi_col_pos, pdblDataReal, pdblDataImag);
	  }
	else
	  {
	    createSparseMatrix(ivar, pirow, picol, nbitem, nb_item_row, pi_col_pos, pdblDataReal);
	  }

	// Free allocated memory
	if (nb_item_row) FREE(nb_item_row);
	if (pi_col_pos)  FREE(pi_col_pos);
	if (pdblDataReal) FREE(pdblDataReal);
	if (iscomplex)
	  if (pdblDataImag) FREE(pdblDataImag);
      } 
      break;
    default:  
      {
	assert(0);
      } break;
    }

  return m_var;
}

/*******************************************************/

gfi_array_list * build_gfi_array_list(int nrhs, int  ** prhs) 
{
  gfi_array_list *l;
  int i;
  
  l = MALLOC(1*sizeof(gfi_array_list));
  l->arg.arg_len = nrhs;
  l->arg.arg_val = MALLOC(nrhs*sizeof(gfi_array));

  for(i=1; i<=nrhs; i++) 
    {
#ifdef DEBUG
      sciprint("build_gfi_array_list: processing parameter %d\n", i);
#endif
      if (sci_array_to_gfi_array(prhs[i], &l->arg.arg_val[i-1]) != 0) return NULL;
    }
  
#ifdef DEBUG
  sciprint("build_gfi_array_list: end of processing\n");
#endif
  
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
