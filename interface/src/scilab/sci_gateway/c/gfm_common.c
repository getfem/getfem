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
#include <signal.h>
#include <string.h>
#include <stdio.h>

#include <api_scilab.h> 
#include <api_stack_common.h> 
#include <Scierror.h>
#include <sciprint.h>
#include <localization.h>

#include "gfm_common.h"

//#define DEBUG

// The spt function is a fortran Scilab function which transpose a sparse matrix
// This function is used to convert a matlab sparse format into a scilab sparse format

#define USE_SPT

// a function to transpose a sparse matrix. From modules/sparse/src/fortran
extern void C2F(spt)(int * m, int * n, int * nel, int * it, int * ptr, double * A_R,  double * A_I,  
	                 int * A_mnel,  int * A_icol, double * At_R, double * At_I, int * At_mnel, int * At_icol);
	    

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

#define CHECK_ERROR_API_SCILAB if(_SciErr.iErr)			    \
                             {					    \
			       printError(&_SciErr, 0);		    \
			       return 0;			    \
			     }					    

int sci_array_to_gfi_array(int * sci_x, gfi_array *t)
{
  SciErr _SciErr;
  int i, n = 0, var_type = 0;
  int * item_address = NULL;
  int pirow, picol, *pilen = NULL, *pilistaddress, size_pistring = 0, nbitem;
  int * nbitemrow = NULL, * picolpos = NULL;
  char ** pstStrings = NULL, ** pstData = NULL;
  double * pdblDataID = NULL, * pdblDataCID = NULL;
  int * pintDims = NULL, * ptr = NULL;
  double * pdblDims = NULL, * pdblDataReal = NULL, * pdblDataImag = NULL;
  int piPrecision, * p_item_address = NULL;
  int * piData32 = NULL, * piBool = NULL;
  unsigned int * puiData32 = NULL;
  int is_complex;
  int tmp_cnt, tmp_cnt2;
  double * tmp_dblDataReal = NULL, * tmp_dblDataImag = NULL;
  
#ifdef DEBUG
  int _picol, _pirow;
  _SciErr = getVarDimension(sci_x,&_pirow,&_picol); CHECK_ERROR_API_SCILAB;
  _SciErr = getVarType(sci_x,&var_type); CHECK_ERROR_API_SCILAB;

  sciprint("sci_array_to_gfi_array: dimension of current variable %d %d - type = %d\n",_pirow,_picol,var_type);
#endif

  assert(t);

  _SciErr = getVarType(pvApiCtx,sci_x,&var_type); CHECK_ERROR_API_SCILAB;

  switch (var_type) 
    {
    case sci_list: 
      {
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: dealing with sci_list\n");
#endif

	_SciErr = getListItemNumber(pvApiCtx,sci_x,&n); CHECK_ERROR_API_SCILAB;

	t->storage.type = GFI_CELL;
	t->storage.gfi_storage_u.data_cell.data_cell_len = n;
	t->storage.gfi_storage_u.data_cell.data_cell_val = (gfi_array**)MALLOC(n*sizeof(gfi_array*));

	for(i=0; i<n; ++i) 
	  {
	    t->storage.gfi_storage_u.data_cell.data_cell_val[i] = MALLOC(1*sizeof(gfi_array));
	    _SciErr = getListItemAddress(pvApiCtx,sci_x,i+1,&item_address); CHECK_ERROR_API_SCILAB;
#ifdef DEBUG
	    _SciErr = getVarType(pvApiCtx,item_address,&var_type); CHECK_ERROR_API_SCILAB;
	    sciprint("type of item %d: %d\n", i+1, var_type);
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

	_SciErr = getListItemNumber(pvApiCtx,sci_x,&n); CHECK_ERROR_API_SCILAB;
	_SciErr = getMatrixOfStringInList(pvApiCtx,sci_x,1,&pirow,&picol,NULL,NULL); CHECK_ERROR_API_SCILAB;
	pilen = (int *)MALLOC(pirow*picol*sizeof(int));
	_SciErr= getMatrixOfStringInList(pvApiCtx,sci_x,1,&pirow,&picol,pilen,NULL); CHECK_ERROR_API_SCILAB;
	pstStrings = (char **)MALLOC(pirow*picol*sizeof(char *));
	for(i=0;i<pirow*picol;i++)
	  {
	    pstStrings[i] = (char *)MALLOC((pilen[i]+1)*sizeof(char));
	  }
	_SciErr = getMatrixOfStringInList(pvApiCtx,sci_x,1,&pirow,&picol,pilen,pstStrings); CHECK_ERROR_API_SCILAB;
	size_pistring = pirow*picol;

#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: pstStrings[0] = %s\n",pstStrings[0]);
#endif

	if (strcmp(pstStrings[0],"objid")==0)
	  {
#ifdef DEBUG
	    sciprint("sci_array_to_gfi_array: dealing with objid\n");
#endif
	    t->storage.type = GFI_OBJID;

	    _SciErr = getMatrixOfDoubleInList(pvApiCtx, sci_x, 2, &pirow, &picol, &pdblDataID); CHECK_ERROR_API_SCILAB;
	    _SciErr = getMatrixOfDoubleInList(pvApiCtx, sci_x, 3, &pirow, &picol, &pdblDataCID); CHECK_ERROR_API_SCILAB;

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
#ifdef DEBUG
	    sciprint("sci_array_to_gfi_array: dealing with hypermat\n");
#endif

	    // Get the dimensions
	    _SciErr = getListItemAddress(pvApiCtx,sci_x,2,&p_item_address); CHECK_ERROR_API_SCILAB;
	    _SciErr = getVarType(pvApiCtx,p_item_address,&var_type); CHECK_ERROR_API_SCILAB

	    switch(var_type)
	      {
	      case sci_matrix:
		_SciErr = getMatrixOfDoubleInList(pvApiCtx,sci_x, 2, &pirow, &picol, &pdblDims); CHECK_ERROR_API_SCILAB;
		t->dim.dim_len = pirow*picol;
		t->dim.dim_val = (u_int*)MALLOC(pirow*picol*sizeof(u_int));
		for(i=0;i<pirow*picol;i++) t->dim.dim_val[i] = (u_int)pdblDims[i];
		break;
	      case sci_ints:
		_SciErr = getMatrixOfInteger32InList(pvApiCtx,sci_x, 2, &pirow, &picol, &pintDims); CHECK_ERROR_API_SCILAB;
		t->dim.dim_len = pirow*picol;
		t->dim.dim_val = (u_int*)MALLOC(pirow*picol*sizeof(u_int));
		for(i=0;i<pirow*picol;i++) t->dim.dim_val[i] = (u_int)pintDims[i];
		break;
	      default:
		_SciErr = getVarType(pvApiCtx,p_item_address,&var_type); CHECK_ERROR_API_SCILAB;
		Scierror(999,"wrong type for hypermatrix dimensions type: %d\n", var_type);
		return 1;
	      }
	    // Get the matrixes (stored as a column vector of size prod(size(...)))
	    // We must detect if we have a INT UINT or DOUBLE
	    
	    _SciErr = getListItemAddress(pvApiCtx, sci_x, 3, &pilistaddress); CHECK_ERROR_API_SCILAB;
	    _SciErr = getVarType(pvApiCtx,pilistaddress,&var_type); CHECK_ERROR_API_SCILAB;

	    if (var_type==sci_matrix)
	      {
		if (isVarComplex(pvApiCtx,pilistaddress))
		  {
		    t->storage.type = GFI_DOUBLE;
		    _SciErr = getComplexMatrixOfDoubleInList(pvApiCtx, sci_x, 3, &pirow, &picol, &pdblDataID, &pdblDataCID); CHECK_ERROR_API_SCILAB;
	    
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
		    _SciErr = getMatrixOfDoubleInList(pvApiCtx, sci_x, 3, &pirow, &picol, &pdblDataID); CHECK_ERROR_API_SCILAB;
		    
		    t->storage.gfi_storage_u.data_double.is_complex = GFI_REAL;

		    t->storage.gfi_storage_u.data_double.data_double_len = pirow*picol;
		    t->storage.gfi_storage_u.data_double.data_double_val = (double *)MALLOC(pirow*picol*sizeof(double));
		    for(i=0;i<pirow*picol;++i) t->storage.gfi_storage_u.data_double.data_double_val[i] = pdblDataID[i];
		  }
	      }
	    else if (var_type==sci_ints)
	      {
		_SciErr = getMatrixOfIntegerPrecision(pvApiCtx,sci_x,&piPrecision); CHECK_ERROR_API_SCILAB;
		if ((piPrecision!=SCI_INT32)&&(piPrecision!=SCI_UINT32))
		  {
		    Scierror(999,"Can deal only with int32 or uint32\n");
		    return 1;
		  }
		
		if (piPrecision==SCI_INT32)
		  {
		    t->storage.type = GFI_INT32;
		    _SciErr = getMatrixOfInteger32(pvApiCtx,sci_x,&pirow,&picol,&piData32); CHECK_ERROR_API_SCILAB;
		    
		    t->storage.gfi_storage_u.data_int32.data_int32_len = pirow*picol;
		    t->storage.gfi_storage_u.data_int32.data_int32_val = (int *)MALLOC(pirow*picol*sizeof(int));
		    for(i=0;i<pirow*picol;++i) t->storage.gfi_storage_u.data_int32.data_int32_val[i] = piData32[i];
		  }
		else if (piPrecision==SCI_UINT32)
		  {
		    t->storage.type = GFI_UINT32;
		    _SciErr = getMatrixOfUnsignedInteger32(pvApiCtx,sci_x,&pirow,&picol,&puiData32); CHECK_ERROR_API_SCILAB;
		    
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
	// Free the allocated memory
	for(i=0;i<size_pistring;i++)
	  {
	    if (pstStrings[i]) FREE(pstStrings[i]);
	  }
	if (pstStrings) FREE(pstStrings);
	if (pilen)      FREE(pilen);
      } 
      break;
    case sci_strings: 
      {
#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: dealing with sci_strings\n");
#endif

	t->storage.type = GFI_CHAR;

	// First call to get picol and pirow
	_SciErr = getMatrixOfString(pvApiCtx,sci_x, &pirow, &picol, NULL, NULL); CHECK_ERROR_API_SCILAB;
	if ((pirow!=1)&&(picol!=1))
	  {
	    Scierror(999,"Can allocate only one string at a time\n");
	    return 1;
	  }
	pilen = (int *)MALLOC(pirow*picol*sizeof(int));

	// Second call to get pilen
	_SciErr = getMatrixOfString(pvApiCtx,sci_x, &pirow, &picol, pilen, NULL); CHECK_ERROR_API_SCILAB;
	pstData = (char **)MALLOC(pirow*picol*sizeof(char*));
	for(i=0; i<pirow*picol ; i++) pstData[i] = (char *)MALLOC((pilen[i] + 1)*sizeof(char));

	n = pilen[0] + 1;
	t->storage.gfi_storage_u.data_char.data_char_len = n;
	t->storage.gfi_storage_u.data_char.data_char_val = MALLOC((n)*sizeof(char));

	// Third call to retrieve data
	_SciErr = getMatrixOfString(pvApiCtx,sci_x, &pirow, &picol, pilen, pstData); CHECK_ERROR_API_SCILAB;
	memcpy(t->storage.gfi_storage_u.data_char.data_char_val,pstData[0],(pilen[0] + 1)*sizeof(char));
#ifdef DEBUG
	sciprint("pirow = %d picol = %d pilen = %d\n", pirow, picol, pilen[0]);
	sciprint("storing |%s|\n",t->storage.gfi_storage_u.data_char.data_char_val);
#endif
	t->dim.dim_len = 1;
	t->dim.dim_val = (u_int*)MALLOC(1*sizeof(u_int));
	t->dim.dim_val[0] = pilen[0];

	for(i=0; i<pirow*picol ; i++) 
	  if (pstData[i]) FREE(pstData[i]);
	if (pstData) FREE(pstData);
	if (pilen)   FREE(pilen);

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

	_SciErr = getMatrixOfIntegerPrecision(pvApiCtx,sci_x,&piPrecision); CHECK_ERROR_API_SCILAB;
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

	    _SciErr = getMatrixOfInteger32(pvApiCtx,sci_x,&pirow,&picol,&piData32); CHECK_ERROR_API_SCILAB;
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

	    _SciErr = getMatrixOfUnsignedInteger32(pvApiCtx,sci_x,&pirow,&picol,&puiData32); CHECK_ERROR_API_SCILAB;

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

	t->storage.type = GFI_INT32;

	_SciErr = getMatrixOfBoolean(pvApiCtx,sci_x,&pirow,&picol,&piBool); CHECK_ERROR_API_SCILAB;

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
	is_complex = isVarComplex(pvApiCtx,sci_x);
	
	t->storage.type = GFI_DOUBLE;
	t->storage.gfi_storage_u.data_double.is_complex = is_complex;
	
	if (!is_complex) 
	  {
	    _SciErr = getMatrixOfDouble(pvApiCtx,sci_x,&pirow,&picol,&pdblDataReal); CHECK_ERROR_API_SCILAB;
	    
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
	    _SciErr = getComplexMatrixOfDouble(pvApiCtx,sci_x,&pirow,&picol,&pdblDataReal,&pdblDataImag); CHECK_ERROR_API_SCILAB;

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
	is_complex = isVarComplex(pvApiCtx,sci_x);

        t->storage.type = GFI_SPARSE;
        t->storage.gfi_storage_u.sp.is_complex = is_complex;

	if (!is_complex) 
	  {
#ifdef DEBUG
	    sciprint("sci_array_to_gfi_array: not complex\n");
#endif
	    _SciErr = getSparseMatrix(pvApiCtx,sci_x,&pirow,&picol,&nbitem,&nbitemrow,&picolpos,&pdblDataReal); CHECK_ERROR_API_SCILAB;
	    
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

#ifdef USE_SPT
	    // We use the spt function to transpose the matlab sparse matrix and so to ease the conversion
	    // to the scilab sparse model.
	    // There is a memory overhead: 
	    // picol integers allocated in the real case
	    // picol integers + 2*nbitem doubles allocated in the complex case

	    ptr = (int *)MALLOC(picol*sizeof(int));

#ifdef DEBUG
	    sciprint("DEBUG: pirow = %d picol = %d nbitem = %d - real case\n", pirow, picol, nbitem);
#endif

	    C2F(spt)(&pirow, &picol, &nbitem, &is_complex, ptr, pdblDataReal,  NULL,  
		     nbitemrow,  picolpos, 
		     t->storage.gfi_storage_u.sp.pr.pr_val, 
		     NULL, 
		     t->storage.gfi_storage_u.sp.jc.jc_val, 
		     t->storage.gfi_storage_u.sp.ir.ir_val);
	    
	    // We compute position offset
	    tmp_cnt = t->storage.gfi_storage_u.sp.jc.jc_val[0];
	    t->storage.gfi_storage_u.sp.jc.jc_val[0] = 0;
	    for(i=1;i<picol;i++)
	      {
		tmp_cnt2 = t->storage.gfi_storage_u.sp.jc.jc_val[i];
		t->storage.gfi_storage_u.sp.jc.jc_val[i] = t->storage.gfi_storage_u.sp.jc.jc_val[i-1] + tmp_cnt;
		tmp_cnt = tmp_cnt2;
	      }
	    t->storage.gfi_storage_u.sp.jc.jc_val[picol] = nbitem;

	    for(i=0;i<nbitem;i++) 
	      {
		t->storage.gfi_storage_u.sp.ir.ir_val[i]--;
	      }

	    if (ptr) FREE(ptr);
#else
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

	    if (offset) FREE(offset);
#endif
	  }
	else
	  {
#ifdef DEBUG
	    sciprint("sci_array_to_gfi_array: complex\n");
#endif
	    _SciErr = getComplexSparseMatrix(pvApiCtx,sci_x,&pirow,&picol,&nbitem,&nbitemrow,&picolpos,&pdblDataReal,&pdblDataImag); CHECK_ERROR_API_SCILAB;
	    
	    // We store the transposed matrix in t
	    t->storage.gfi_storage_u.sp.is_complex = GFI_COMPLEX;

	    t->storage.gfi_storage_u.sp.ir.ir_len = nbitem;
	    t->storage.gfi_storage_u.sp.jc.jc_len = picol+1;
	    t->storage.gfi_storage_u.sp.pr.pr_len = 2*nbitem;

	    t->storage.gfi_storage_u.sp.ir.ir_val = (int *)MALLOC(nbitem*sizeof(int));
	    t->storage.gfi_storage_u.sp.jc.jc_val = (int *)MALLOC((picol+1)*sizeof(int));
	    t->storage.gfi_storage_u.sp.pr.pr_val = (double *)MALLOC(2*nbitem*sizeof(double));

#ifdef USE_SPT
	    sciprint("DEBUG: testing the new method - complex\n");

	    ptr = (int *)MALLOC(picol*sizeof(int));
	    tmp_dblDataReal = (double *)MALLOC(nbitem*sizeof(double));
	    tmp_dblDataImag = (double *)MALLOC(nbitem*sizeof(double));

	    sciprint("pirow = %d picol = %d nbitem = %d\n", pirow, picol, nbitem);

	    C2F(spt)(&pirow, &picol, &nbitem, &is_complex, ptr, pdblDataReal,  pdblDataImag,  
		     nbitemrow,  picolpos, 
		     tmp_dblDataReal, 
		     tmp_dblDataImag, 
		     t->storage.gfi_storage_u.sp.jc.jc_val, 
		     t->storage.gfi_storage_u.sp.ir.ir_val);
	    
	    // We compute position offset
	    tmp_cnt = t->storage.gfi_storage_u.sp.jc.jc_val[0];
	    t->storage.gfi_storage_u.sp.jc.jc_val[0] = 0;
	    for(i=1;i<picol;i++)
	      {
		tmp_cnt2 = t->storage.gfi_storage_u.sp.jc.jc_val[i];
		t->storage.gfi_storage_u.sp.jc.jc_val[i] = t->storage.gfi_storage_u.sp.jc.jc_val[i-1] + tmp_cnt;
		tmp_cnt = tmp_cnt2;
	      }
	    t->storage.gfi_storage_u.sp.jc.jc_val[picol] = nbitem;

	    for(i=0;i<nbitem;i++) 
	      {
		t->storage.gfi_storage_u.sp.ir.ir_val[i]--;
	      }

	    // Now store the real + imag data in getfem
	    for(i=0;i<nbitem;i++)
	      {
		t->storage.gfi_storage_u.sp.pr.pr_val[2*i+0] = tmp_dblDataReal[i];
		t->storage.gfi_storage_u.sp.pr.pr_val[2*i+1] = tmp_dblDataImag[i];
	      }

	    if (tmp_dblDataReal) FREE(tmp_dblDataReal);
	    if (tmp_dblDataImag) FREE(tmp_dblDataImag);
	    if (ptr)             FREE(ptr);
#else
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

	    if (offset) FREE(offset);
#endif
	  }

#ifdef DEBUG
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
      }
      break;
    default: 
      {
	_SciErr = getVarType(pvApiCtx,sci_x,&var_type); CHECK_ERROR_API_SCILAB;
	Scierror(999,"unhandled class type : %s\n", sci_ClassID2string(var_type));
	return 1;
      } 
      break;
    }

  return 0;
}

int gfi_array_to_sci_array(gfi_array *t, int ivar) 
{
  // Revoir cette partie. Surtout la partie GFI_CELL ...
  // Ajouter des fonctions
  // - exportString pour les variables simples
  // - addCellToCell, addStringToCell, etc.. pour gérer les listes imbriquées

  SciErr _SciErr;
  int * m_var = NULL, var_type;

  /* Scilab represent scalars as an array of size one */
  /* while gfi_array represents "scalar" values with 0-dimension array */

  int ndim = (t->dim.dim_len == 0 ? 1 : t->dim.dim_len);
  static const int one = 1;
  const int * dim = (t->dim.dim_len == 0 ? &one : (const int *)t->dim.dim_val);
  int is_hypermat = 0;
  int nrow, ncol, nb_item, pi_precision;
  int pirow, picol, nbitem, * pi_col_pos = NULL, * nb_item_row = NULL;
  int * nb_item_row_tmp = NULL, * pi_col_pos_tmp = NULL, * pilen = NULL;
  int * dims = NULL;
  unsigned int * entries = NULL;
  int nb_elem = 1;
  double * pr = NULL, * pi = NULL;
  int i, j;
  double * pdblDataReal = NULL, * pdblDataImag = NULL, * pdblDataReal_tmp = NULL, * pdblDataImag_tmp = NULL;
  double * dbl_entries = NULL;
  double * entries_pr = NULL, * entries_pi = NULL;
  int iscomplex;
  int * ptr = NULL, * m_content = NULL, * piData32 = NULL;
  char ** pstStrings = NULL;
  double * pdblReal = NULL, * pdblImag = NULL;
  unsigned int * puiData32 = NULL;
  int * piBool = NULL;

  assert(t);
  
  switch (t->storage.type) 
    {
    case GFI_UINT32: 
      {
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
	    _SciErr = createMatrixOfUnsignedInteger32(pvApiCtx,ivar,dim[0],dim[1],t->storage.gfi_storage_u.data_uint32.data_uint32_val); CHECK_ERROR_API_SCILAB;
	  }
	else
	  {
	    const char * const fields[] = {"hm","dims","entries"};
	    nb_elem = 1;

	    _SciErr = createMList(pvApiCtx,ivar,3,&m_var); CHECK_ERROR_API_SCILAB;
	    _SciErr = createMatrixOfStringInList(pvApiCtx,ivar, m_var, 1, 1, 3, fields); CHECK_ERROR_API_SCILAB;
	    
	    dims = (int *)MALLOC(t->dim.dim_len*sizeof(int));
	    for(i=0;i<t->dim.dim_len;i++) 
	      {
		dims[i] = (int)t->dim.dim_val[i];
		nb_elem *= (int)t->dim.dim_val[i];
	      }

	    entries = (unsigned int *)MALLOC(nb_elem*sizeof(unsigned int));
	    for(i=0;i<nb_elem;i++) entries[i] = t->storage.gfi_storage_u.data_uint32.data_uint32_val[i];
	    
	    // Add a vector to the 'dims' field -> a row vector
	    _SciErr = createMatrixOfInteger32InList(pvApiCtx,ivar, m_var, 2, 1, t->dim.dim_len, dims); CHECK_ERROR_API_SCILAB;
	    // Add a vector to the 'entries' field -> a column vector
	    _SciErr = createMatrixOfUnsignedInteger32InList(pvApiCtx,ivar, m_var, 3, nb_elem, 1, entries); CHECK_ERROR_API_SCILAB;
		
	    if (dims)    FREE(dims);
	    if (entries) FREE(entries);
	  }
      } 
      break;
    case GFI_INT32: 
      {
	is_hypermat = 0;

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
	    _SciErr = createMatrixOfInteger32(pvApiCtx,ivar,dim[0],dim[1],t->storage.gfi_storage_u.data_int32.data_int32_val); CHECK_ERROR_API_SCILAB;
	  }
	else
	  {
	    const char * const fields[] = {"hm","dims","entries"};
	    nb_elem = 1;

	    _SciErr = createMList(pvApiCtx,ivar,3,&m_var); CHECK_ERROR_API_SCILAB;
	    _SciErr = createMatrixOfStringInList(pvApiCtx,ivar, m_var, 1, 1, 3, fields); CHECK_ERROR_API_SCILAB;
	    
	    dims = (int *)MALLOC(t->dim.dim_len*sizeof(int));
	    for(i=0;i<t->dim.dim_len;i++) 
	      {
		dims[i] = (int)t->dim.dim_val[i];
		nb_elem *= (int)t->dim.dim_val[i];
	      }

	    dbl_entries = (double *)MALLOC(nb_elem*sizeof(double));
	    for(i=0;i<nb_elem;i++) dbl_entries[i] = (double)t->storage.gfi_storage_u.data_int32.data_int32_val[i];
	    
	    // Add a vector to the 'dims' field -> a row vector
	    _SciErr = createMatrixOfInteger32InList(pvApiCtx,ivar, m_var, 2, 1, t->dim.dim_len, dims); CHECK_ERROR_API_SCILAB;
	    // Add a vector to the 'entries' field -> a column vector
	    _SciErr = createMatrixOfDoubleInList(pvApiCtx,ivar, m_var, 3, nb_elem, 1, dbl_entries); CHECK_ERROR_API_SCILAB;
		
	    if (entries) FREE(dbl_entries);
	    if (dims)    FREE(dims);
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
	is_hypermat = 0;

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
		_SciErr = createMatrixOfDouble(pvApiCtx,ivar, nrow, ncol, t->storage.gfi_storage_u.data_double.data_double_val); CHECK_ERROR_API_SCILAB;
	      } 
	    else 
	      {
#ifdef DEBUG
		sciprint("DEBUG: array is complex\n");
#endif
		
		pr = (double *)MALLOC(nrow*ncol*sizeof(double));
		pi = (double *)MALLOC(nrow*ncol*sizeof(double));
		
		for(i=0;i<nrow*ncol;i++)
		  {
		    pr[i] = t->storage.gfi_storage_u.data_double.data_double_val[2*i];
		    pi[i] = t->storage.gfi_storage_u.data_double.data_double_val[2*i+1];
		  }
		
		_SciErr = createComplexMatrixOfDouble(pvApiCtx,ivar, nrow, ncol, pr, pi); CHECK_ERROR_API_SCILAB;
		
		if (pr) FREE(pr);
		if (pi) FREE(pi);
	      }
	  }
	else
	  {
#ifdef DEBUG
	    sciprint("DEBUG: array is hypermat\n");
#endif
	    const char * const fields[] = {"hm","dims","entries"};
	    nb_elem = 1;

	    _SciErr = createMList(pvApiCtx,ivar,3,&m_var); CHECK_ERROR_API_SCILAB;
	    _SciErr = createMatrixOfStringInList(pvApiCtx,ivar, m_var, 1, 1, 3, fields); CHECK_ERROR_API_SCILAB;
	    
	    dims = (int *)MALLOC(t->dim.dim_len*sizeof(int));
	    for(i=0;i<t->dim.dim_len;i++) 
	      {
		dims[i] = (int)t->dim.dim_val[i];
		nb_elem *= (int)t->dim.dim_val[i];
	      }

	    if (!gfi_array_is_complex(t)) 
	      {

		dbl_entries = (double *)MALLOC(nb_elem*sizeof(double));
		for(i=0;i<nb_elem;i++) dbl_entries[i] = t->storage.gfi_storage_u.data_double.data_double_val[i];
		
		// Add a vector to the 'dims' field -> a row vector
		_SciErr = createMatrixOfInteger32InList(pvApiCtx,ivar, m_var, 2, 1, t->dim.dim_len, dims); CHECK_ERROR_API_SCILAB;
		// Add a vector to the 'entries' field -> a column vector
		_SciErr = createMatrixOfDoubleInList(pvApiCtx,ivar, m_var, 3, nb_elem, 1, dbl_entries); CHECK_ERROR_API_SCILAB;
		
		if (entries) FREE(entries);
#ifdef DEBUG
		sciprint("DEBUG: end array is hypermat\n");
#endif
	      }
	    else
	      {
		entries_pr = (double *)MALLOC(nb_elem*sizeof(double));
		entries_pi = (double *)MALLOC(nb_elem*sizeof(double));
		for(i=0;i<nb_elem;i++) 
		  {
		    entries_pr[i] = t->storage.gfi_storage_u.data_double.data_double_val[2*i+0];
		    entries_pi[i] = t->storage.gfi_storage_u.data_double.data_double_val[2*i+1];
		  }
		
		// Add a vector to the 'dims' field
		_SciErr = createMatrixOfInteger32InList(pvApiCtx,ivar, m_var, 2, 1, t->dim.dim_len, dims); CHECK_ERROR_API_SCILAB;
		// Add a vector to the 'entries' field
		_SciErr = createComplexMatrixOfDoubleInList(pvApiCtx,ivar, m_var, 3, 1, nb_elem, entries_pr, entries_pi); CHECK_ERROR_API_SCILAB;
		
		if (entries_pr) FREE(entries_pr);
		if (entries_pi) FREE(entries_pi);
	      }
	    if (dims) FREE(dims);
	  }
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

	_SciErr = createMatrixOfString(pvApiCtx,ivar, 1, 1, &tmp_string); CHECK_ERROR_API_SCILAB;
#ifdef DEBUG
	sciprint("ivar = %d string = |%s| len = %d\n",ivar, tmp_string,t->storage.gfi_storage_u.data_char.data_char_len);
#endif
	if (tmp_string) FREE(tmp_string);
      } 
      break;
    case GFI_CELL: 
      {
#ifdef DEBUG
	sciprint("gfi_array_to_sci_array: create from a GFI_CELL\n");
	sciprint("ndim = %d ivar = %d\n", dim[0],ivar);

	for(i=0;i<ndim;i++)
	  sciprint("dim[%d] = %d\n",i,dim[i]);

	sciprint("now create list at pos %d, dimension %d\n", ivar, dim[0]);
#endif
	_SciErr = createList(pvApiCtx,ivar,dim[0],&m_var); CHECK_ERROR_API_SCILAB;
	
	for(i=0; i<t->storage.gfi_storage_u.data_cell.data_cell_len; ++i)
	  {
	    m_content = gfi_array_to_sci_array(t->storage.gfi_storage_u.data_cell.data_cell_val[i],ivar+i+1);

	    _SciErr = getVarType(pvApiCtx,m_content,&var_type); CHECK_ERROR_API_SCILAB;
	    switch(var_type)
	      {
	      case sci_list:
		{
#ifdef DEBUG
		  sciprint("gfi_array_to_sci_array: create from a GFI_CELL - sci_list\n");
#endif
		  _SciErr = getListItemNumber(pvApiCtx,m_content,&nb_item); CHECK_ERROR_API_SCILAB;
		  _SciErr = createListInList(pvApiCtx,ivar, m_content, i+1, nb_item, &m_var); CHECK_ERROR_API_SCILAB;
		}
		break;
	      case sci_mlist:
		{
#ifdef DEBUG
		  sciprint("gfi_array_to_sci_array: create from a GFI_CELL - sci_mlist\n");
#endif
		  _SciErr = getListItemNumber(pvApiCtx,m_content,&nb_item); CHECK_ERROR_API_SCILAB;
		  _SciErr = createMListInList(pvApiCtx,ivar, m_content, i+1, nb_item, &m_var); CHECK_ERROR_API_SCILAB;
		}
		break;
	      case sci_strings:
		{
#ifdef DEBUG
		  sciprint("gfi_array_to_sci_array: create from a GFI_CELL - sci_strings\n");
		  sciprint("add element %d to position %d\n", i, ivar);
#endif
		  
		  // Get the matrix of strings from the gfi_array
		  _SciErr = getMatrixOfString(pvApiCtx,m_content, &pirow, &picol, NULL, NULL); CHECK_ERROR_API_SCILAB;
		  pilen = (int *)MALLOC(pirow*picol*sizeof(int));
		  pstStrings = (char **)MALLOC(pirow*picol*sizeof(char *));
		  _SciErr = getMatrixOfString(pvApiCtx,m_content, &pirow, &picol, pilen, NULL); CHECK_ERROR_API_SCILAB;

		  for(j=0;j<pirow*picol;j++)
		    {
		      pstStrings[j] = (char *)MALLOC((pilen[j]+1)*sizeof(char));
		    }
		  _SciErr = getMatrixOfString(pvApiCtx,m_content, &pirow, &picol, pilen, pstStrings); CHECK_ERROR_API_SCILAB;
#ifdef DEBUG
		  sciprint("pirow = %d picol = %d, pilen[0] = %d\n", pirow, picol, pilen[0]);
#endif
		  // And now add it to the list
		  _SciErr = createMatrixOfStringInList(pvApiCtx,ivar, m_var, i+1, pirow, picol, pstStrings); CHECK_ERROR_API_SCILAB;
		  
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
		  _SciErr = getMatrixOfIntegerPrecision(pvApiCtx,m_content,&pi_precision); CHECK_ERROR_API_SCILAB;
		  if ((pi_precision!=SCI_INT32)&&(pi_precision!=SCI_UINT32))
		    {
		      Scierror(999,"Can deal only with int32 or uint32\n");
		    }
		  switch(pi_precision)
		    {
		    case SCI_INT32:
		      _SciErr = getMatrixOfInteger32(pvApiCtx,m_content, &pirow, &picol, &piData32); CHECK_ERROR_API_SCILAB;
		      _SciErr = createMatrixOfInteger32InList(pvApiCtx,ivar, m_var, i+1, pirow, picol, piData32); CHECK_ERROR_API_SCILAB;
		      break;
		    case SCI_UINT32:
		      _SciErr = getMatrixOfUnsignedInteger32(pvApiCtx,m_content, &pirow, &picol, &puiData32); CHECK_ERROR_API_SCILAB;
		      _SciErr = createMatrixOfUnsignedInteger32InList(pvApiCtx,ivar, m_var, i+1, pirow, picol, puiData32); CHECK_ERROR_API_SCILAB;
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

		  _SciErr = getMatrixOfBoolean(pvApiCtx,m_content, &pirow, &picol, &piBool); CHECK_ERROR_API_SCILAB;
		  _SciErr = createMatrixOfBooleanInList(pvApiCtx,ivar, m_var, i+1, pirow, picol, piBool); CHECK_ERROR_API_SCILAB;
		}
		break;
	      case sci_matrix:
		{
#ifdef DEBUG
		  sciprint("gfi_array_to_sci_array: create from a GFI_CELL - sci_matrix\n");
#endif
		  if (isVarComplex(pvApiCtx,m_content))
		    {
		      _SciErr = getComplexMatrixOfDouble(pvApiCtx,m_content, &pirow, &picol, &pdblReal, &pdblImag); CHECK_ERROR_API_SCILAB;
		      _SciErr = createComplexMatrixOfDoubleInList(pvApiCtx,ivar, m_var, i+1, pirow, picol, pdblReal, pdblImag); CHECK_ERROR_API_SCILAB;
		    }
		  else
		    {
		      _SciErr = getMatrixOfDouble(pvApiCtx,m_content, &pirow, &picol, &pdblReal); CHECK_ERROR_API_SCILAB;
		      _SciErr = createMatrixOfDoubleInList(pvApiCtx,ivar, m_var, i+1, pirow, picol, pdblReal); CHECK_ERROR_API_SCILAB;
		    }
		}
		break;
	      case sci_sparse:
		{
#ifdef DEBUG
		  sciprint("gfi_array_to_sci_array: create from a GFI_CELL - sci_sparse\n");
#endif
		  if (isVarComplex(pvApiCtx,m_content))
		    {
		      _SciErr = getComplexSparseMatrix(pvApiCtx,m_content, &pirow, &picol, &nbitem, &nb_item_row, &pi_col_pos, &pdblReal, &pdblImag); CHECK_ERROR_API_SCILAB;
		      _SciErr = createComplexSparseMatrixInList(pvApiCtx,ivar, m_var, i+1, pirow, picol, nbitem, nb_item_row, pi_col_pos, pdblReal, pdblImag); CHECK_ERROR_API_SCILAB;
		    }
		  else
		    {
		      _SciErr = getSparseMatrix(pvApiCtx,m_content, &pirow, &picol, &nbitem, &nb_item_row, &pi_col_pos, &pdblReal); CHECK_ERROR_API_SCILAB;
		      _SciErr = createSparseMatrixInList(pvApiCtx,ivar, m_var, i+1, pirow, picol, nbitem, nb_item_row, pi_col_pos, pdblReal); CHECK_ERROR_API_SCILAB;
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

	_SciErr = createMList(pvApiCtx,ivar,3,&m_var); CHECK_ERROR_API_SCILAB;
	_SciErr = createMatrixOfStringInList(pvApiCtx,ivar, m_var, 1, 1, 3, fields); CHECK_ERROR_API_SCILAB;

	size_objid  = t->storage.gfi_storage_u.objid.objid_len;
	pdblDataID  = (double *)MALLOC(size_objid*sizeof(double));
	pdblDataCID = (double *)MALLOC(size_objid*sizeof(double));

	for(i=0; i<size_objid; ++i) 
	  {
	    pdblDataID[i]  = t->storage.gfi_storage_u.objid.objid_val[i].id;       
	    pdblDataCID[i] = t->storage.gfi_storage_u.objid.objid_val[i].cid;
	  }
	// Add a vector to the 'id' field
	_SciErr = createMatrixOfDoubleInList(pvApiCtx,ivar, m_var, 2, 1, size_objid, pdblDataID); CHECK_ERROR_API_SCILAB;
	// Add a vector to the 'cid' field
	_SciErr = createMatrixOfDoubleInList(pvApiCtx,ivar, m_var, 3, 1, size_objid, pdblDataCID); CHECK_ERROR_API_SCILAB;

	if (pdblDataID)  FREE(pdblDataID);
	if (pdblDataCID) FREE(pdblDataCID);
      } 
      break;
    case GFI_SPARSE: 
      {
	iscomplex = gfi_array_is_complex(t);

	nbitem = t->storage.gfi_storage_u.sp.ir.ir_len;
	pirow  = t->dim.dim_val[0];
	picol  = t->dim.dim_val[1];

	// Convert from Matlab to Scilab format
	nb_item_row  = (int *)MALLOC(pirow*sizeof(int));
	pi_col_pos   = (int *)MALLOC(nbitem*sizeof(int));
	pdblDataReal = (double *)MALLOC(nbitem*sizeof(double));
	if (iscomplex)
	  pdblDataImag = (double *)MALLOC(nbitem*sizeof(double));
	
#ifdef USE_SPT
	// We use the spt function to transpose the matlab sparse matrix and so to ease the conversion
	// to the scilab sparse model.
	// There is a memory overhead: 
	// (picol + pirow + nbitem) integers allocated in the real case
	// (picol + pirow + nbitem) integers + 2*nbitem doubles allocated in the complex case
	nb_item_row_tmp = (int *)MALLOC(picol*sizeof(int));
	pi_col_pos_tmp  = (int *)MALLOC(nbitem*sizeof(int));
	ptr             = (int *)MALLOC(pirow*sizeof(int));


	if (iscomplex)
	  {
#ifdef DEBUG
	    sciprint("DEBUG: pirow = %d picol = %d nbitem = %d - complex\n", pirow, picol, nbitem);
#endif
	    pdblDataImag_tmp = (double *)MALLOC(nbitem*sizeof(double));
	    pdblDataReal_tmp = (double *)MALLOC(nbitem*sizeof(double));
	    
	    for(i=0;i<nbitem;i++)
	      {
		pdblDataReal_tmp[i] = t->storage.gfi_storage_u.sp.pr.pr_val[2*i+0];
		pdblDataImag_tmp[i] = t->storage.gfi_storage_u.sp.pr.pr_val[2*i+1];
		pi_col_pos_tmp[i]   = t->storage.gfi_storage_u.sp.ir.ir_val[i] + 1;
	      }
	    
	    for(i=0;i<picol;i++)
	      {
		nb_item_row_tmp[i] = t->storage.gfi_storage_u.sp.jc.jc_val[i+1] - t->storage.gfi_storage_u.sp.jc.jc_val[i];
	      }
	    
	    // Transpose the matrix
	    C2F(spt)(&picol, &pirow, &nbitem, &iscomplex, ptr, 
		     pdblDataReal_tmp,  
		     pdblDataImag_tmp,  
		     nb_item_row_tmp, 
		     pi_col_pos_tmp,
		     pdblDataReal, 
		     pdblDataImag, 
		     nb_item_row, 
		     pi_col_pos);
	    
	    if (pdblDataImag_tmp) FREE(pdblDataImag_tmp);
	    if (pdblDataReal_tmp) FREE(pdblDataReal_tmp);
	  }
	else
	  {
#ifdef DEBUG
	    sciprint("DEBUG: pirow = %d picol = %d nbitem = %d - real\n", pirow, picol, nbitem);
#endif
	    for(i=0;i<nbitem;i++)
	      {
		pi_col_pos_tmp[i] = t->storage.gfi_storage_u.sp.ir.ir_val[i] + 1;
	      }
	    
	    for(i=0;i<picol;i++)
	      {
		nb_item_row_tmp[i] = t->storage.gfi_storage_u.sp.jc.jc_val[i+1] - t->storage.gfi_storage_u.sp.jc.jc_val[i];
	      }
	    
	    // Transpose the matrix
	    C2F(spt)(&picol, &pirow, &nbitem, &iscomplex, ptr, 
		     t->storage.gfi_storage_u.sp.pr.pr_val,  
		     NULL,  
		     nb_item_row_tmp, 
		     pi_col_pos_tmp,
		     pdblDataReal, 
		     NULL, 
		     nb_item_row, 
		     pi_col_pos);
	  }

	if (ptr)             FREE(ptr);
	if (nb_item_row_tmp) FREE(nb_item_row_tmp);
	if (pi_col_pos_tmp)  FREE(pi_col_pos_tmp);
#else
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
			    if ((j>=t->storage.gfi_storage_u.sp.jc.jc_val[k])&&(j<=t->storage.gfi_storage_u.sp.jc.jc_val[k+1]-1))
			      {
				pi_col_pos[Index] = k+1;
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
			    if ((j>=t->storage.gfi_storage_u.sp.jc.jc_val[k])&&(j<=t->storage.gfi_storage_u.sp.jc.jc_val[k+1]-1))
			      {
				pi_col_pos[Index] = k+1;
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
#endif

	if (iscomplex)
	  {
	    _SciErr = createComplexSparseMatrix(pvApiCtx,ivar, pirow, picol, nbitem, nb_item_row, pi_col_pos, pdblDataReal, pdblDataImag); CHECK_ERROR_API_SCILAB;
	  }
	else
	  {
	    _SciErr = createSparseMatrix(pvApiCtx,ivar, pirow, picol, nbitem, nb_item_row, pi_col_pos, pdblDataReal); CHECK_ERROR_API_SCILAB;
	  }

	// Free allocated memory
	if (nb_item_row)  FREE(nb_item_row);
	if (pi_col_pos)   FREE(pi_col_pos);
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
  return 0;
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

static void sigint(int sig) 
{
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

void remove_custom_sigint(int allow_rethrow) 
{
#ifndef WIN32
  struct sigaction act;
  sigaction (SIGINT, NULL, &act);
  if (act.sa_handler == sigint) 
    {
      sigaction(SIGINT, &old_sigint, NULL);
    }
  if (allow_rethrow && sigint_hit) 
    {
      fprintf(stderr, "ready, raising SIGINT now\n");
      raise(SIGINT); 
    }
  sigint_hit = 0; 
#endif
}
