#include <stdio.h>

#include <api_common.h>
#include <api_sparse.h>
#include <api_double.h>
#include <MALLOC.h>
#include <stack-c.h>

#include <sparse2.h>
#include <iter.h>
#include <err.h>

#define DEBUG

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

// iter_spmgcr - a simple interface to iter_mgcr 
// VEC * iter_spmgcr(SPMAT * A, SPMAT * B, VEC * b, double tol, VEC * x, int k, int limit, int * steps)

int sci_spmgcr(char * fname)
{
  // [x,[iter]] = pmgcr(A,b,tol,[maxit,[k,[B,[x0]]]])
  int * A_pi_address, A_pi_nb_rows, A_pi_nb_cols, A_pi_nb_items, * A_pi_nb_items_row, * A_pi_col_pos;
  double * A_pdbl_real;
  int * B_pi_address, B_pi_nb_rows, B_pi_nb_cols, B_pi_nb_items, * B_pi_nb_items_row, * B_pi_col_pos;
  double * B_pdbl_real;
  int * b_pi_address, b_pi_nb_rows, b_pi_nb_cols;
  double * b_pdbl_real;
  int * tol_pi_address, tol_pi_nb_rows, tol_pi_nb_cols;
  double * tol_pdbl_real;
  int * maxit_pi_address, maxit_pi_nb_rows, maxit_pi_nb_cols;
  double * maxit_pdbl_real;
  int * k_pi_address, k_pi_nb_rows, k_pi_nb_cols;
  double * k_pdbl_real;
  int * x0_pi_address, x0_pi_nb_rows, x0_pi_nb_cols;
  double * x0_pdbl_real;
  int * xsol_pi_address, xsol_pi_nb_rows, xsol_pi_nb_cols;
  double * xsol_pdbl_real;
  int * iter_pi_address, iter_pi_nb_rows, iter_pi_nb_cols;
  double * iter_pdbl_real;
  SciErr _SciErr;
  StrCtx _StrCtx;
  int var_type;
  SPMAT  * A = NULL, * B = NULL;
  VEC * b = NULL, * x0 = NULL, * xsol = NULL;
  int Index, steps, i, j, k;

  CheckRhs(3,7);
  CheckLhs(1,2);

  // Get A
  _SciErr = getVarAddressFromPosition(&_StrCtx,1,&A_pi_address);

  _SciErr = getVarType(&_StrCtx,A_pi_address,&var_type);
  if (var_type!=sci_sparse)
    {
      Scierror(999,"%s: wrong parameter, a sparse matrix is needed\n",fname);
      return 0;
    }

  if (isVarComplex(&_StrCtx,A_pi_address))
    {
      Scierror(999,"%s: wrong parameter, a real sparse matrix is needed\n",fname);
      return 0;
    }

  _SciErr = getSparseMatrix(&_StrCtx,A_pi_address, &A_pi_nb_rows, &A_pi_nb_cols, 
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
  _SciErr = getVarAddressFromPosition(&_StrCtx,2,&b_pi_address);
  _SciErr = getMatrixOfDouble(&_StrCtx,b_pi_address, &b_pi_nb_rows, &b_pi_nb_cols, &b_pdbl_real);

  // Convert Scilab vector into VEC
  b  = v_get(b_pi_nb_rows);
  for(i=0;i<b_pi_nb_rows;i++)
    {
      v_set_val(b,i,b_pdbl_real[i]);
    }

  // Get tol
  _SciErr = getVarAddressFromPosition(&_StrCtx,3,&tol_pi_address);
  _SciErr = getMatrixOfDouble(&_StrCtx,tol_pi_address, &tol_pi_nb_rows, &tol_pi_nb_cols, &tol_pdbl_real);

   // Get optional maxit
  if (Rhs>=4)
    {
      _SciErr = getVarAddressFromPosition(&_StrCtx,4,&maxit_pi_address);
      _SciErr = getMatrixOfDouble(&_StrCtx,maxit_pi_address, &maxit_pi_nb_rows, &maxit_pi_nb_cols, &maxit_pdbl_real);
    }

   // Get optional k
  if (Rhs>=5)
    {
      _SciErr = getVarAddressFromPosition(&_StrCtx,5,&k_pi_address);
      _SciErr = getMatrixOfDouble(&_StrCtx,k_pi_address, &k_pi_nb_rows, &k_pi_nb_cols, &k_pdbl_real);
    }

  // Get optional B
  if (Rhs>=6)
    {
      _SciErr = getVarAddressFromPosition(&_StrCtx,6,&B_pi_address);
      _SciErr = getVarType(&_StrCtx,B_pi_address,&var_type);
      if (var_type!=sci_sparse)
	{
	  Scierror(999,"%s: wrong parameter, a sparse matrix is needed\n",fname);
	  return 0;
	}
      
      if (isVarComplex(&_StrCtx,B_pi_address))
	{
	  Scierror(999,"%s: wrong parameter, a real sparse matrix is needed\n",fname);
	  return 0;
	}

      _SciErr = getSparseMatrix(&_StrCtx,B_pi_address, &B_pi_nb_rows, &B_pi_nb_cols, 
				&B_pi_nb_items, &B_pi_nb_items_row, &B_pi_col_pos, &B_pdbl_real);

      // Convert SPMAT into Scilab sparse
      B = sp_get(B_pi_nb_rows, B_pi_nb_cols, 5);
      Index = 0;
      for(i=0;i<B_pi_nb_rows;i++)
	{
	  for(j=0;j<B_pi_nb_items_row[i];j++)
	    {
	      sp_set_val(B,i,B_pi_col_pos[Index]-1, B_pdbl_real[Index]);
	      Index++;
	    }
	}
    }

  // Get optional x0
  if (Rhs>=7)
    {
      _SciErr = getVarAddressFromPosition(&_StrCtx,7,&x0_pi_address);
      _SciErr = getMatrixOfDouble(&_StrCtx,x0_pi_address, &x0_pi_nb_rows, &x0_pi_nb_cols, &x0_pdbl_real);

      // Convert Scilab vector into VEC
      x0 = v_get(x0_pi_nb_rows);
      for(i=0;i<x0_pi_nb_rows;i++)
	{
	  v_set_val(x0,i,x0_pdbl_real[i]);
	}
    }
  
  // call iter_spmgcr method.
  catchall(xsol = iter_spmgcr(A, B, b, *tol_pdbl_real, x0, k, (int)*maxit_pdbl_real, &steps),
	   Scierror(999,"%s: an error occured.\n",fname); return 0);

  // Transfert xsol to Scilab
  xsol_pdbl_real = (double *)MALLOC(b_pi_nb_rows*sizeof(double));
  memcpy(xsol_pdbl_real,xsol->ve,b_pi_nb_rows*sizeof(double));
  xsol_pi_nb_rows = b_pi_nb_rows;
  xsol_pi_nb_cols = 1;
  _SciErr = createMatrixOfDouble(&_StrCtx, Rhs+1, xsol_pi_nb_rows, xsol_pi_nb_cols, xsol_pdbl_real);

  LhsVar(1) = Rhs+1;

  if (Lhs>=2)
    {
      iter_pdbl_real  = (double *)MALLOC(1*sizeof(double));
      *iter_pdbl_real = (double)steps;
      iter_pi_nb_rows = 1;
      iter_pi_nb_cols = 1;
      _SciErr = createMatrixOfDouble(&_StrCtx, Rhs+2, iter_pi_nb_rows, iter_pi_nb_cols, iter_pdbl_real);

      LhsVar(2) = Rhs+2;
    }

  if (A)    sp_free(A);
  if (B)    sp_free(B);
  if (b)    v_free(b);
  if (x0)   v_free(x0);
  if (xsol) v_free(xsol);

  if (xsol_pdbl_real)   FREE(xsol_pdbl_real);
  if (iter_pdbl_real)   FREE(iter_pdbl_real);

  return 0;
}