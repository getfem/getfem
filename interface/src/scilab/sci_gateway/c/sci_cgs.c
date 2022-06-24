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

#include <stdio.h>

#include <gfm_common.h>
#include <Scierror.h>
#include <sciprint.h>

#include <sparse2.h>
#include <iter.h>
#include <err.h>

//#define DEBUG

// x = cgs(A,b)
// cgs(A,b,tol)
// cgs(A,b,tol,maxit)
// cgs(A,b,tol,maxit,M)
// cgs(A,b,tol,maxit,M1,M2)
// cgs(A,b,tol,maxit,M1,M2,x0)
// [x,flag] = cgs(A,b,...)
// [x,flag,relres] = cgs(A,b,...)
// [x,flag,relres,iter] = cgs(A,b,...)
// [x,flag,relres,iter,resvec] = cgs(A,b,...)

// k : no. of direction (search) vectors; =0 - none
// maxit: upper bound on the no. of iter. steps
// steps: no. of iter. steps done 
// tol: accuracy required

// iter_spcgs -- simple interface for SPMAT data structures 
// use always as follows:
//    x = iter_spcgs(A,B,b,r0,tol,x,limit,steps);
// or 
//    x = iter_spcgs(A,B,b,r0,tol,VNULL,limit,steps);
// In the second case the solution vector is created.  
// If B is not NULL then it is a preconditioner. 
// VEC * iter_spcgs(SPMAT * A, SPMAT * B, VEC * b, VEC * r0, double tol, VEC * x, int limit, int * steps)

int sci_spcgs(char * fname)
{
  // [x,flag,[relres,[iter,[resvec]]]] = cgs(A,b,tol,[maxit,[M,[x0]]])
  int * A_pi_address = NULL, A_pi_nb_rows, A_pi_nb_cols, A_pi_nb_items, * A_pi_nb_items_row = NULL, * A_pi_col_pos = NULL;
  double * A_pdbl_real = NULL;
  int * b_pi_address = NULL, b_pi_nb_rows, b_pi_nb_cols;
  double * b_pdbl_real = NULL;
  int * tol_pi_address = NULL, tol_pi_nb_rows, tol_pi_nb_cols;
  double * tol_pdbl_real = NULL;
  int * maxit_pi_address = NULL, maxit_pi_nb_rows, maxit_pi_nb_cols;
  double * maxit_pdbl_real = NULL;
  int * M_pi_address = NULL, M_pi_nb_rows, M_pi_nb_cols, M_pi_nb_items, * M_pi_nb_items_row = NULL, * M_pi_col_pos = NULL;
  double * M_pdbl_real = NULL;
  int * x0_pi_address = NULL, x0_pi_nb_rows, x0_pi_nb_cols;
  double * x0_pdbl_real = NULL;
  int xsol_pi_nb_rows, xsol_pi_nb_cols;
  double * xsol_pdbl_real = NULL;
  int iter_pi_nb_rows, iter_pi_nb_cols;
  double * iter_pdbl_real = NULL;
  int resvec_pi_nb_rows, resvec_pi_nb_cols;
  double * resvec_pdbl_real = NULL;
  SciErr _SciErr;
  int var_type;
  SPMAT  * A = NULL;
  VEC * b = NULL, * x0 = NULL, * r0 = NULL, * xsol = NULL;
  SPMAT * M = NULL;
  int Index, steps, i, j;

  CheckRhs(3,7);
  CheckLhs(1,5);

  // Get A
  _SciErr = getVarAddressFromPosition(pvApiCtx,1,&A_pi_address);

  _SciErr = getVarType(pvApiCtx,A_pi_address,&var_type);
  if (var_type!=sci_sparse)
    {
      Scierror(999,"%s: wrong parameter, a sparse matrix is needed\n",fname);
      return 0;
    }

  if (isVarComplex(pvApiCtx,A_pi_address))
    {
      Scierror(999,"%s: wrong parameter, a real sparse matrix is needed\n",fname);
      return 0;
    }

  _SciErr = getSparseMatrix(pvApiCtx,A_pi_address, &A_pi_nb_rows, &A_pi_nb_cols, 
			    &A_pi_nb_items, &A_pi_nb_items_row, &A_pi_col_pos, &A_pdbl_real);

  // Convert Scilab sparse into SPMAT
  A = sp_get(A_pi_nb_rows, A_pi_nb_cols, 5);
  Index = 0;
  for(i=0;i<A_pi_nb_rows;i++)
    {
      for(j=0;j<A_pi_nb_items_row[i];j++)
	{
	  sp_set_val(A,i,A_pi_col_pos[Index]-1, A_pdbl_real[Index]);
	  Index++;
	}
    }

  // Get b
  _SciErr = getVarAddressFromPosition(pvApiCtx,2,&b_pi_address);
  _SciErr = getMatrixOfDouble(pvApiCtx,b_pi_address, &b_pi_nb_rows, &b_pi_nb_cols, &b_pdbl_real);

  // Convert Scilab vector into VEC
  b  = v_get(b_pi_nb_rows);
  r0 = v_get(b_pi_nb_rows);
  for(i=0;i<b_pi_nb_rows;i++)
    {
      v_set_val(b,i,b_pdbl_real[i]);
      v_set_val(r0,i,1.0);
    }

  // Get tol
  _SciErr = getVarAddressFromPosition(pvApiCtx,3,&tol_pi_address);
  _SciErr = getMatrixOfDouble(pvApiCtx,tol_pi_address, &tol_pi_nb_rows, &tol_pi_nb_cols, &tol_pdbl_real);

  // Get optional maxit
  if (Rhs>=4)
    {
      _SciErr = getVarAddressFromPosition(pvApiCtx,4,&maxit_pi_address);
      _SciErr = getMatrixOfDouble(pvApiCtx,maxit_pi_address, &maxit_pi_nb_rows, &maxit_pi_nb_cols, &maxit_pdbl_real);
    }

  // Get optional M
  if (Rhs>=5)
    {
      _SciErr = getVarAddressFromPosition(pvApiCtx,5,&M_pi_address);
      _SciErr = getVarType(pvApiCtx,M_pi_address,&var_type);
      if (var_type!=sci_sparse)
	{
	  Scierror(999,"%s: wrong parameter, a sparse matrix is needed\n",fname);
	  return 0;
	}
      
      if (isVarComplex(pvApiCtx,M_pi_address))
	{
	  Scierror(999,"%s: wrong parameter, a real sparse matrix is needed\n",fname);
	  return 0;
	}

      _SciErr = getSparseMatrix(pvApiCtx,M_pi_address, &M_pi_nb_rows, &M_pi_nb_cols, 
				&M_pi_nb_items, &M_pi_nb_items_row, &M_pi_col_pos, &M_pdbl_real);

      // Convert SPMAT into Scilab sparse
      M = sp_get(M_pi_nb_rows, M_pi_nb_cols, 5);
      Index = 0;
      for(i=0;i<M_pi_nb_rows;i++)
	{
	  for(j=0;j<M_pi_nb_items_row[i];j++)
	    {
	      sp_set_val(M,i,M_pi_col_pos[Index]-1, M_pdbl_real[Index]);
	      Index++;
	    }
	}
    }

  // Get optional x0
  if (Rhs>=6)
    {
      _SciErr = getVarAddressFromPosition(pvApiCtx,6,&x0_pi_address);
      _SciErr = getMatrixOfDouble(pvApiCtx,x0_pi_address, &x0_pi_nb_rows, &x0_pi_nb_cols, &x0_pdbl_real);

      // Convert Scilab vector into VEC
      x0 = v_get(x0_pi_nb_rows);
      for(i=0;i<x0_pi_nb_rows;i++)
	{
	  v_set_val(x0,i,x0_pdbl_real[i]);
	}
    }
  else
    {
      x0 = v_get(b_pi_nb_rows);
      for(i=0;i<b_pi_nb_rows;i++)
	{
	  v_set_val(x0,i,0.0);
	}
    }

  // call iter_spcgs method.

  catchall(xsol = iter_spcgs(A, M, b, r0, *tol_pdbl_real, x0, (int)*maxit_pdbl_real, &steps),
  	   Scierror(999,"%s: an error (%d) occured.\n",fname,_err_num); return 0);

  // Transfert xsol to Scilab
  xsol_pdbl_real = (double *)MALLOC(b_pi_nb_rows*sizeof(double));
  memcpy(xsol_pdbl_real,xsol->ve,b_pi_nb_rows*sizeof(double));
  xsol_pi_nb_rows = b_pi_nb_rows;
  xsol_pi_nb_cols = 1;
  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+1, xsol_pi_nb_rows, xsol_pi_nb_cols, xsol_pdbl_real);
  if (xsol_pdbl_real) FREE(xsol_pdbl_real);

  LhsVar(1) = Rhs+1;

  if (Lhs>=2)
    {
      iter_pdbl_real  = (double *)MALLOC(1*sizeof(double));
      *iter_pdbl_real = (double)steps;
      iter_pi_nb_rows = 1;
      iter_pi_nb_cols = 1;
      _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+2, iter_pi_nb_rows, iter_pi_nb_cols, iter_pdbl_real);
      if (iter_pdbl_real) FREE(iter_pdbl_real);

      LhsVar(2) = Rhs+2;
    }

  if (Lhs>=3)
    {
      resvec_pdbl_real = (double *)MALLOC(b_pi_nb_rows*sizeof(double));
      memcpy(resvec_pdbl_real,r0->ve,b_pi_nb_rows*sizeof(double));
      resvec_pi_nb_rows = b_pi_nb_rows;
      resvec_pi_nb_cols = 1;
      _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+3, resvec_pi_nb_rows, resvec_pi_nb_cols, resvec_pdbl_real);
      if (resvec_pdbl_real) FREE(resvec_pdbl_real);

      LhsVar(3) = Rhs+3;
    }

  if (A)    sp_free(A);
  if (b)    v_free(b);
  if (x0)   v_free(x0);
  if (r0)   v_free(r0);
  //if (xsol) v_free(xsol);
  if (M)    sp_free(M);

  return 0;
}
