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
  int * A_pi_address, A_pi_nb_rows, A_pi_nb_cols, A_pi_nb_items, * A_pi_nb_items_row, * A_pi_col_pos;
  double * A_pdbl_real;
  int * b_pi_address, b_pi_nb_rows, b_pi_nb_cols;
  double * b_pdbl_real;
  int * tol_pi_address, tol_pi_nb_rows, tol_pi_nb_cols;
  double * tol_pdbl_real;
  int * maxit_pi_address, maxit_pi_nb_rows, maxit_pi_nb_cols;
  double * maxit_pdbl_real;
  int * M_pi_address, M_pi_nb_rows, M_pi_nb_cols, M_pi_nb_items, * M_pi_nb_items_row, * M_pi_col_pos;
  double * M_pdbl_real;
  int * x0_pi_address, x0_pi_nb_rows, x0_pi_nb_cols;
  double * x0_pdbl_real;
  int * xsol_pi_address, xsol_pi_nb_rows, xsol_pi_nb_cols;
  double * xsol_pdbl_real;
  int * iter_pi_address, iter_pi_nb_rows, iter_pi_nb_cols;
  double * iter_pdbl_real;
  int * resvec_pi_address, resvec_pi_nb_rows, resvec_pi_nb_cols;
  double * resvec_pdbl_real;
  SciErr _SciErr;
  StrCtx _StrCtx;
  int var_type;
  SPMAT  * A = NULL;
  VEC * b = NULL, * x0 = NULL, * r0 = NULL, * xsol = NULL;
  SPMAT * M = NULL;
  int Index, steps, i, j;

  CheckRhs(3,7);
  CheckLhs(1,5);

  // Get A
#ifdef DEBUG
  sciprint("get A\n");
#endif
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
#ifdef DEBUG
  sciprint("get b\n");
#endif
  _SciErr = getVarAddressFromPosition(&_StrCtx,2,&b_pi_address);
  _SciErr = getMatrixOfDouble(&_StrCtx,b_pi_address, &b_pi_nb_rows, &b_pi_nb_cols, &b_pdbl_real);

  // Convert Scilab vector into VEC
  b  = v_get(b_pi_nb_rows);
  r0 = v_get(b_pi_nb_rows);
  for(i=0;i<b_pi_nb_rows;i++)
    {
      v_set_val(b,i,b_pdbl_real[i]);
      v_set_val(r0,i,1.0);
    }

  // Get tol
#ifdef DEBUG
  sciprint("get tol\n");
#endif
  _SciErr = getVarAddressFromPosition(&_StrCtx,3,&tol_pi_address);
  _SciErr = getMatrixOfDouble(&_StrCtx,tol_pi_address, &tol_pi_nb_rows, &tol_pi_nb_cols, &tol_pdbl_real);

  // Get optional maxit
  if (Rhs>=4)
    {
#ifdef DEBUG
      sciprint("get maxit\n");
#endif
      _SciErr = getVarAddressFromPosition(&_StrCtx,4,&maxit_pi_address);
      _SciErr = getMatrixOfDouble(&_StrCtx,maxit_pi_address, &maxit_pi_nb_rows, &maxit_pi_nb_cols, &maxit_pdbl_real);
    }

  // Get optional M
  if (Rhs>=5)
    {
#ifdef DEBUG
      sciprint("get M\n");
#endif
      _SciErr = getVarAddressFromPosition(&_StrCtx,5,&M_pi_address);
      _SciErr = getVarType(&_StrCtx,M_pi_address,&var_type);
      if (var_type!=sci_sparse)
	{
	  Scierror(999,"%s: wrong parameter, a sparse matrix is needed\n",fname);
	  return 0;
	}
      
      if (isVarComplex(&_StrCtx,M_pi_address))
	{
	  Scierror(999,"%s: wrong parameter, a real sparse matrix is needed\n",fname);
	  return 0;
	}

      _SciErr = getSparseMatrix(&_StrCtx,M_pi_address, &M_pi_nb_rows, &M_pi_nb_cols, 
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
#ifdef DEBUG
  sciprint("get x0\n");
#endif
  if (Rhs>=6)
    {
      _SciErr = getVarAddressFromPosition(&_StrCtx,6,&x0_pi_address);
      _SciErr = getMatrixOfDouble(&_StrCtx,x0_pi_address, &x0_pi_nb_rows, &x0_pi_nb_cols, &x0_pdbl_real);

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
  
#ifdef DEBUG
  sciprint("starting cgs done\n");
  sciprint("A = %d - row = %d, col = %d\n", A, A_pi_nb_rows, A_pi_nb_cols);
  sciprint("M = %d - row = %d, col = %d\n", M, M_pi_nb_rows, M_pi_nb_cols);
  sciprint("b = %d - len = %d\n", b,b_pi_nb_cols*b_pi_nb_rows);
  sciprint("r0 = %d\n", r0);
  sciprint("x0 = %d\n", x0);
  sciprint("tol = %f - len = %d\n",*tol_pdbl_real,tol_pi_nb_cols*tol_pi_nb_rows);
  sciprint("maxit = %d - len = %d\n", (int)*maxit_pdbl_real,maxit_pi_nb_cols*maxit_pi_nb_rows);
#endif

  // call iter_spcgs method.
  catchall(xsol = iter_spcgs(A, M, b, r0, *tol_pdbl_real, x0, (int)*maxit_pdbl_real, &steps),
  	   Scierror(999,"%s: an error (%d) occured.\n",fname,_err_num); return 0);

#ifdef DEBUG
  sciprint("cgs done - xsol = %d\n",xsol->dim);
#endif

  // Transfert xsol to Scilab
  xsol_pdbl_real = (double *)MALLOC(b_pi_nb_rows*sizeof(double));
  memcpy(xsol_pdbl_real,xsol->ve,b_pi_nb_rows*sizeof(double));
  xsol_pi_nb_rows = b_pi_nb_rows;
  xsol_pi_nb_cols = 1;
  _SciErr = createMatrixOfDouble(&_StrCtx, Rhs+1, xsol_pi_nb_rows, xsol_pi_nb_cols, xsol_pdbl_real);
  if (xsol_pdbl_real) FREE(xsol_pdbl_real);

  LhsVar(1) = Rhs+1;

  if (Lhs>=2)
    {
      iter_pdbl_real  = (double *)MALLOC(1*sizeof(double));
      *iter_pdbl_real = (double)steps;
      iter_pi_nb_rows = 1;
      iter_pi_nb_cols = 1;
      _SciErr = createMatrixOfDouble(&_StrCtx, Rhs+2, iter_pi_nb_rows, iter_pi_nb_cols, iter_pdbl_real);
      if (iter_pdbl_real) FREE(iter_pdbl_real);

      LhsVar(2) = Rhs+2;
    }

  if (Lhs>=3)
    {
      resvec_pdbl_real = (double *)MALLOC(b_pi_nb_rows*sizeof(double));
      memcpy(resvec_pdbl_real,r0->ve,b_pi_nb_rows*sizeof(double));
      resvec_pi_nb_rows = b_pi_nb_rows;
      resvec_pi_nb_cols = 1;
      _SciErr = createMatrixOfDouble(&_StrCtx, Rhs+3, resvec_pi_nb_rows, resvec_pi_nb_cols, resvec_pdbl_real);
      if (resvec_pdbl_real) FREE(resvec_pdbl_real);

      LhsVar(3) = Rhs+3;
    }

#ifdef DEBUG
  sciprint("free memory\n");
#endif

  if (A)    sp_free(A);
  if (b)    v_free(b);
  if (x0)   v_free(x0);
  if (r0)   v_free(r0);
  //if (xsol) v_free(xsol);
  if (M)    sp_free(M);

  return 0;
}
