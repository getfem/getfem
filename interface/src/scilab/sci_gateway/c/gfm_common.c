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

#define DEBUG

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

int sci_array_to_gfi_array(int * sci_x, int ivar, gfi_array *t)
{
  int i, n = 0;

#ifdef DEBUG
  int _picol, _pirow;
  getVarDimension(sci_x,&_pirow,&_picol);
  sciprint("sci_array_to_gfi_array: dimension of current variable %d %d\n",_pirow,_picol);
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
	    getListItemAddress(sci_x,i, &item_address);
	    if (sci_array_to_gfi_array(item_address, ivar, t->storage.gfi_storage_u.data_cell.data_cell_val[i]) != 0) return 1;
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
	sciprint("sci_array_to_gfi_array: dealing with hypermat\n");
#endif
	int pirow, picol, *pilen;
	char ** pstStrings;
	double * pdblDataID, * pdblDataCID;

	getListItemNumber(sci_x,&n);
	getMatrixOfStringInList(ivar,sci_x,1,&pirow,&picol,NULL,NULL);
	pilen = (int *)MALLOC(pirow*picol*sizeof(int));
	getMatrixOfStringInList(ivar,sci_x,1,&pirow,&picol,pilen,NULL);
	pstStrings = (char **)MALLOC(pirow*picol*sizeof(char *));
	for(i=0;i<pirow*picol;i++)
	  {
	    pstStrings[i] = (char *)MALLOC((pilen[i]+1)*sizeof(char));
	  }
	getMatrixOfStringInList(ivar,sci_x,1,&pirow,&picol,pilen,pstStrings);

	// On doit traiter les hypermatrice
	// Et les GFI_OBJID
	// On récupère les labels

#ifdef DEBUG
	sciprint("sci_array_to_gfi_array: pstStrings[0] = %s\n",pstStrings[0]);
#endif

	if (strcmp(pstStrings[0],"objid")==0)
	  {
#ifdef DEBUG
	    sciprint("sci_array_to_gfi_array: dealing with objid\n");
#endif
	    t->storage.type = GFI_OBJID;

	    getMatrixOfDoubleInList(ivar, sci_x, 2, &pirow, &picol, &pdblDataID);
	    getMatrixOfDoubleInList(ivar, sci_x, 3, &pirow, &picol, &pdblDataCID);

#ifdef DEBUG
	    sciprint("sci_array_to_gfi_array: pirow = %d picol = %d\n", pirow, picol);
#endif

	    n = pirow*picol;

	    t->storage.gfi_storage_u.objid.objid_len = n;
	    t->storage.gfi_storage_u.objid.objid_val = (struct gfi_object_id *)MALLOC(n*sizeof(struct gfi_object_id));

	    for(i=0;i<n;i++)
	      {
		t->storage.gfi_storage_u.objid.objid_val[i].id  = pdblDataID[i];
		t->storage.gfi_storage_u.objid.objid_val[i].cid = pdblDataCID[i];
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
	    // On vérifie que c'est une hm
	    // On récupère les dimensions dans dims (2nd label)
	    // On récupère le vecteur de double dans entries (3ieme label)
	    exit(1);
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
	for(i=0; i<pirow*picol ; i++) pstData[i] = (char *)MALLOC((pilen[i] + 1)*sizeof(char));//+ 1 for null termination

	n = pilen[0] + 1; // YC:
	t->storage.gfi_storage_u.data_char.data_char_len = n;
	t->storage.gfi_storage_u.data_char.data_char_val = MALLOC((n)*sizeof(char));

	// Third call to retrieve data
	iret = getMatrixOfString(sci_x, &pirow, &picol, pilen, pstData);
	memcpy(t->storage.gfi_storage_u.data_char.data_char_val,pstData[0],(pilen[0] + 1)*sizeof(char));
#ifdef DEBUG
	sciprint("pirow = %d picol = %d pilen = %d\n", pirow, picol, pilen[0]);
	sciprint("storing |%s|\n",pstData[0]);
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

	    n = picol*pirow;
	    t->storage.gfi_storage_u.data_int32.data_int32_len = n;
	    t->storage.gfi_storage_u.data_int32.data_int32_val = (int *)MALLOC(n*sizeof(in));
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
	int * nbitemrow, * picolpos;
	int j, Index;

        t->storage.type = GFI_SPARSE;
        t->storage.gfi_storage_u.sp.is_complex = is_complex;

	if (!is_complex) 
	  {
	    getSparseMatrix(sci_x,&pirow,&picol,&nbitem,&nbitemrow,&picolpos,&pdblDataReal);
	    
	    t->storage.gfi_storage_u.sp.ir.ir_len = nbitem;
	    t->storage.gfi_storage_u.sp.jc.jc_len = nbitem;
	    t->storage.gfi_storage_u.sp.pr.pr_len = nbitem;

	    t->storage.gfi_storage_u.sp.ir.ir_val = (int *)MALLOC(nbitem*sizeof(int));
	    t->storage.gfi_storage_u.sp.jc.jc_val = (int *)MALLOC(nbitem*sizeof(int));
	    t->storage.gfi_storage_u.sp.pr.pr_val = (double *)MALLOC(nbitem*sizeof(double));

	    Index = 0;
	    for(i=0;i<pirow;++i)
	      {
		for(j=0;j<nbitemrow[i];++j)
		  {
		    t->storage.gfi_storage_u.sp.ir.ir_val[Index] = i;
		    t->storage.gfi_storage_u.sp.jc.jc_val[Index] = picolpos[Index]-1;
		    t->storage.gfi_storage_u.sp.pr.pr_val[Index] = pdblDataReal[Index];
		  }
	      }
	  }
	else
	  {
	    getComplexSparseMatrix(sci_x,&pirow,&picol,&nbitem,&nbitemrow,&picolpos,&pdblDataReal,&pdblDataImag);
	    
	    t->storage.gfi_storage_u.sp.ir.ir_len = nbitem;
	    t->storage.gfi_storage_u.sp.jc.jc_len = nbitem;
	    t->storage.gfi_storage_u.sp.pr.pr_len = 2*nbitem;

	    t->storage.gfi_storage_u.sp.ir.ir_val = (int *)MALLOC(nbitem*sizeof(int));
	    t->storage.gfi_storage_u.sp.jc.jc_val = (int *)MALLOC(nbitem*sizeof(int));
	    t->storage.gfi_storage_u.sp.pr.pr_val = (double *)MALLOC(2*nbitem*sizeof(double));

	    Index = 0;
	    for(i=0;i<pirow;++i)
	      {
		for(j=0;j<nbitemrow[i];++j)
		  {
		    t->storage.gfi_storage_u.sp.ir.ir_val[Index] = i;
		    t->storage.gfi_storage_u.sp.jc.jc_val[Index] = picolpos[Index]-1;
		    t->storage.gfi_storage_u.sp.pr.pr_val[2*Index]   = pdblDataReal[Index];
		    t->storage.gfi_storage_u.sp.pr.pr_val[2*Index+1] = pdblDataImag[Index];
		  }
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

  int *m;

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
#ifdef DEBUG
	sciprint("gfi_array_to_sci_array: create from a GFI_UINT32\n");
	sciprint("ndim = %d ivar = %d\n", ndim,ivar);
#endif
	createMatrixOfUnsignedInteger32(ivar,ndim,dim[0],t->storage.gfi_storage_u.data_uint32.data_uint32_val);
	getVarAddressFromPosition(ivar, &m);
      } 
      break;
    case GFI_INT32: 
      {
#ifdef DEBUG
	sciprint("gfi_array_to_sci_array: create from a GFI_INT32\n");
	sciprint("ndim = %d ivar = %d\n", ndim,ivar);
#endif
	createMatrixOfInteger32(ivar,ndim,dim[0],t->storage.gfi_storage_u.data_int32.data_int32_val);
	getVarAddressFromPosition(ivar, &m);
      } 
      break;
    case GFI_DOUBLE: 
      {
#ifdef DEBUG
	sciprint("gfi_array_to_sci_array: create from a GFI_DOUBLE\n");
	sciprint("ndim = %d ivar = %d\n", ndim,ivar);
#endif
	int i;
	double *pr, *pi;
	if (!gfi_array_is_complex(t)) 
	  {
	    allocMatrixOfDouble(ivar, dim[0], dim[1], &t->storage.gfi_storage_u.data_double.data_double_val);
	  } 
	else 
	  {
	    double *pr, *pi; int i;

	    pr = (double *)MALLOC(dim[0]*dim[1]*sizeof(double));
	    pi = (double *)MALLOC(dim[0]*dim[1]*sizeof(double));

	    for(i=0;i<dim[0]*dim[1];i++)
	      {
		pr[i] = t->storage.gfi_storage_u.data_double.data_double_val[2*i];
		pi[i] = t->storage.gfi_storage_u.data_double.data_double_val[2*i+1];
	      }

	    createComplexMatrixOfDouble(ivar, dim[0], dim[1], pr, pi);

	    FREE(pr);
	    FREE(pi);
	  }
	getVarAddressFromPosition(ivar, &m);
      } 

      break;
    case GFI_CHAR: 
      {
#ifdef DEBUG
	sciprint("gfi_array_to_sci_array: create from a GFI_CHAR\n");
	sciprint("ndim = %d ivar = %d\n", ndim,ivar);
#endif
	char * tmp_string = (char *)MALLOC((t->storage.gfi_storage_u.data_char.data_char_len+1)*sizeof(char));
	memcpy(tmp_string,t->storage.gfi_storage_u.data_char.data_char_val,t->storage.gfi_storage_u.data_char.data_char_len*sizeof(char));
	tmp_string[t->storage.gfi_storage_u.data_char.data_char_len] = '\0';
	createMatrixOfString(ivar, 1, 1, &tmp_string);
	getVarAddressFromPosition(ivar, &m);
#ifdef DEBUG
	sciprint("string = |%s| len = %d\n",tmp_string,t->storage.gfi_storage_u.data_char.data_char_len);
#endif
	if (tmp_string) FREE(tmp_string);
      } 
      break;
    case GFI_CELL: 
      {
	unsigned int i, nb_item;
	int * m_list;
#ifdef DEBUG
	sciprint("gfi_array_to_sci_array: create from a GFI_CELL\n");
	sciprint("ndim = %d ivar = %d\n", ndim,ivar);

	for(i=0;i<ndim;i++)
	  sciprint("dim[%d] = %d\n",i,dim[i]);
#endif
	createList(ivar,ndim,&m);

	for(i=0; i<t->storage.gfi_storage_u.data_cell.data_cell_len; ++i)
	  {
	    // YC: il vaudrait mieux ne pas appeler de façon récursive gfi_array_to_sci_array
	    // on ne peut pas créer 2 variables au même endroit sur la stack avec des adresses differentes.
	    // Peut-on avoir des GFI_CELL dans des GFI_CELL ??? oui ...

	    m_list = gfi_array_to_sci_array(t->storage.gfi_storage_u.data_cell.data_cell_val[i],ivar);

	    switch(getVarType(m_list))
	      {
	      case sci_list:
		{
#ifdef DEBUG
	sciprint("gfi_array_to_sci_array: create from a GFI_CELL - sci_list\n");
#endif
		  getListItemNumber(m,&nb_item);
		  createListInList(ivar, m, i, nb_item, &m);
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
		  getMatrixOfString(m_list, &pirow, &picol, NULL, NULL);
		  pilen = (int *)MALLOC(pirow*picol*sizeof(int));
		  pstStrings = (char **)MALLOC(pirow*picol*sizeof(char *));
		  getMatrixOfString(m_list, &pirow, &picol, pilen, NULL);
		  for(j=0;j<pirow*picol;j++)
		    {
		      pstStrings[j] = (char *)MALLOC((pilen[j]+1)*sizeof(char));
		    }
		  getMatrixOfString(m_list, &pirow, &picol, pilen, pstStrings);
		  // And now add it to the list
		  createMatrixOfStringInList(ivar, m, i, pirow, picol, pstStrings);
		  
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
		  
		  getMatrixOfIntegerPrecision(m,&pi_precision);
		  if ((pi_precision!=SCI_INT32)&&(pi_precision!=SCI_UINT32))
		    {
		      Scierror(999,"Can deal only with int32 or uint32\n");
		      return 1;
	  }
		  switch(pi_precision)
		    {
		    case SCI_INT32:
		      getMatrixOfInteger32(m_list, &pirow, &picol, &piData32);
		      createMatrixOfInteger32InList(ivar, m, i, pirow, picol, piData32);
		      break;
		    case SCI_UINT32:
		      getMatrixOfUnsignedInteger32(m_list, &pirow, &picol, &puiData32);
		      createMatrixOfUnsignedInteger32InList(ivar, m, i, pirow, picol, puiData32);
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
		  getMatrixOfBoolean(m_list, &pirow, &picol, &piBool);
		  
		  createMatrixOfBooleanInList(ivar, m, i, pirow, picol, piBool);
		}
		break;
	      case sci_matrix:
		{
#ifdef DEBUG
		  sciprint("gfi_array_to_sci_array: create from a GFI_CELL - sci_matrix\n");
#endif
		  int pirow, picol;
		  double * pdblReal, * pdblImag;
		  
		  if (isVarComplex(m_list))
		    {
		      getComplexMatrixOfDouble(m_list, &pirow, &picol, &pdblReal, &pdblImag);
		      createComplexMatrixOfDoubleInList(ivar, m, i, pirow, picol, pdblReal, pdblImag);
		    }
		  else
		    {
		      getMatrixOfDouble(m_list, &pirow, &picol, &pdblReal);
		      createMatrixOfDoubleInList(ivar, m, i, pirow, picol, pdblReal);
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
		  
		  if (isVarComplex(m_list))
		    {
		      getComplexSparseMatrix(m_list, &pirow, &picol, &nbitem, &nb_item_row, &pi_col_pos, &pdblReal, &pdblImag);
		      createComplexSparseMatrix(ivar, pirow, picol, nbitem, nb_item_row, pi_col_pos, pdblReal, pdblImag);
		    }
		  else
		    {
		      getSparseMatrix(m_list, &pirow, &picol, &nbitem, &nb_item_row, &pi_col_pos, &pdblReal);
		      createSparseMatrixInList(ivar, m, 1, pirow, picol, nbitem, nb_item_row, pi_col_pos, pdblReal);
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
	static const char *fields[] = {"objid","id","cid"};

	createMList(ivar,3,&m);

	createMatrixOfStringInList(ivar, m, 1, 1, 3, fields);

	size_objid  = t->storage.gfi_storage_u.objid.objid_len;
	pdblDataID  = (double *)MALLOC(size_objid*sizeof(double));
	pdblDataCID = (double *)MALLOC(size_objid*sizeof(double));

	for(i=0; i<size_objid; ++i) 
	  {
	    pdblDataID[i]  = t->storage.gfi_storage_u.objid.objid_val[i].id;       
	    pdblDataCID[i] = t->storage.gfi_storage_u.objid.objid_val[i].cid;
	  }
	// Add a vector to the 'id' field
	createMatrixOfDoubleInList(ivar, m, 2, 1, size_objid, pdblDataID);

	// Add a vector to the 'cid' field
	createMatrixOfDoubleInList(ivar, m, 3, 1, size_objid, pdblDataCID);

	if (pdblDataID)  FREE(pdblDataID);
	if (pdblDataCID) FREE(pdblDataCID);
      } 
      break;
    case GFI_SPARSE: 
      {
	int picol, pirow, nbitem, * pi_col_pos, * nb_item_row;
	int i, j, nnz, Index;
	double * pdblDataReal, * pdblDataImag;
	int iscomplex = gfi_array_is_complex(t);

	int *sp_icol, *sp_jrow;
	double * sp_real, *sp_imag;

	nbitem = t->storage.gfi_storage_u.sp.pr.pr_len;
	pirow  = t->dim.dim_val[0];
	picol  = t->dim.dim_val[1];

	sp_icol = (int *)MALLOC(nbitem*sizeof(int));
	sp_jrow = (int *)MALLOC(nbitem*sizeof(int));
	sp_real = (double *)MALLOC(nbitem*sizeof(double));
	if (iscomplex)
	  sp_imag = (double *)MALLOC(nbitem*sizeof(double));

	// Copy the sparse matrix content
	for (i=0; i<nnz; i++) 
	  {
	    sp_icol[i] = t->storage.gfi_storage_u.sp.ir.ir_val[i];
	    sp_jrow[i] = t->storage.gfi_storage_u.sp.jc.jc_val[i];
	    if (iscomplex)
	      {
		sp_real[i] = t->storage.gfi_storage_u.sp.pr.pr_val[2*i];
		sp_imag[i] = t->storage.gfi_storage_u.sp.pr.pr_val[2*i+1];
	      }
	    else
	      {
		sp_real[i] = t->storage.gfi_storage_u.sp.pr.pr_val[i];
	      }
	  }

	// Convert to Scilab format
	nb_item_row = (int *)MALLOC(pirow*sizeof(int));
	pi_col_pos  = (int *)MALLOC(nbitem*sizeof(int));
	pdblDataReal = (double *)MALLOC(nbitem*sizeof(double));
	if (iscomplex)
	  pdblDataImag = (double *)MALLOC(nbitem*sizeof(double));

	Index = 0;
	for(i=0; i<pirow; ++i)
	  {
	    for(j=0; j<nbitem; ++j)
	      {
		if (sp_jrow[j]==i) 
		  {
		    nb_item_row[i]++;
		    pi_col_pos[Index] = sp_icol[j];
		    pdblDataReal[Index] = sp_real[j];
		    if (iscomplex)
		      pdblDataImag[Index] = sp_imag[j];
		    Index++;
		  }
	      }
	  }

	if (iscomplex)
	  {
	    createComplexSparseMatrix(ivar, pirow, picol, nbitem, nb_item_row, pi_col_pos, pdblDataReal, pdblDataImag);
	  }
	else
	  {
	    createSparseMatrix(ivar, pirow, picol, nbitem, nb_item_row, pi_col_pos, pdblDataReal);
	  }

	// Free allocated memory
	if (sp_icol) FREE(sp_icol);
	if (sp_jrow) FREE(sp_jrow);
	if (sp_real) FREE(sp_real);
	if (iscomplex)
	  if (sp_imag) FREE(sp_imag);

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

  return m;
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
      if (sci_array_to_gfi_array(prhs[i], i, &l->arg.arg_val[i-1]) != 0) return NULL;
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
