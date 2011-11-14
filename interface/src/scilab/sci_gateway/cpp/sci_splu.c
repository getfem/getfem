#include <stdio.h>

#include <api_common.h>
#include <api_sparse.h>
#include <api_double.h>
#include <MALLOC.h>
#include <stack-c.h>

#include "sparse2.h"
#include "err.h"

#define DEBUG

int sci_splu(char * fname)
{
  int      p_in_spmat_nb_rows, p_in_spmat_nb_cols, p_in_spmat_nb_items;
  int    * p_in_spmat_address;
  int    * p_in_spmat_items_row = NULL;
  int    * p_in_spmat_col_pos   = NULL;
  double * p_in_spmat_val       = NULL;
  int      p_in_dbl_nb_rows, p_in_dbl_nb_cols;
  double * p_in_dbl_matrix  = NULL;
  int    * p_in_dbl_address = NULL;
  SPMAT  * A = NULL;
  PERM   * pivot = NULL;
  int      Index, i, j, res;
  int    * p_out_spmat_item_row = NULL;
  int    * p_out_spmat_col_pos  = NULL;
  double * p_out_spmat_val      = NULL;
  double   value, alpha = 1.0;
  int      nnz = 0, var_type;
  SciErr _SciErr;
  StrCtx _StrCtx;

  CheckRhs(1,2);
  CheckLhs(1,3);

  // First, access to the input variable (a matrix of strings)
  _SciErr = getVarAddressFromPosition(&_StrCtx,1,&p_in_spmat_address);

  _SciErr = getVarType(&_StrCtx,p_in_spmat_address,&var_type);
  if (var_type!=sci_sparse)
    {
      Scierror(999,"%s: wrong parameter, a sparse matrix is needed\n",fname);
      return 0;
    }

  if (isVarComplex(&_StrCtx,p_in_spmat_address))
    {
      Scierror(999,"%s: wrong parameter, a real sparse matrix is needed\n",fname);
      return 0;
    }

  _SciErr = getSparseMatrix(&_StrCtx,p_in_spmat_address, &p_in_spmat_nb_rows, &p_in_spmat_nb_cols, 
			    &p_in_spmat_nb_items, &p_in_spmat_items_row, &p_in_spmat_col_pos, &p_in_spmat_val);

  if (Rhs==2)
    {
      // Second, get the alpha parameter
      // First, access to the input variable (a matrix of doubles)
      _SciErr = getVarAddressFromPosition(&_StrCtx,2,&p_in_dbl_address);
      
      _SciErr = getMatrixOfDouble(&_StrCtx,p_in_dbl_address, &p_in_dbl_nb_rows, &p_in_dbl_nb_cols, &p_in_dbl_matrix);
      alpha = *p_in_dbl_matrix;
    }

  // Proceed the factorization
  A = sp_get(p_in_spmat_nb_rows, p_in_spmat_nb_cols, 5);
  Index = 0;
  for(i=0;i<p_in_spmat_nb_rows;i++)
    {
      for(j=0;j<p_in_spmat_items_row[i];j++)
	{
	  sp_set_val(A,i,p_in_spmat_col_pos[Index]-1, p_in_spmat_val[Index]);
	  Index++;
	}
    }

  pivot = px_get(A->m);

  catchall(spLUfactor(A,pivot,alpha),Scierror(999,"%s: an error occured.\n",fname); return 0);

  // Now, create the result
  for(i=0;i<A->m; i++) nnz += A->row[i].len;

  p_out_spmat_item_row = (int *)MALLOC(p_in_spmat_nb_rows*sizeof(int));
  p_out_spmat_col_pos  = (int *)MALLOC(nnz*sizeof(int));
  p_out_spmat_val      = (double *)MALLOC(nnz*sizeof(double));

  // Get the L matrix
  if (Lhs>=1)
    {
      Index = 0;
      for(i=0;i<p_in_spmat_nb_rows;i++)
	{
	  p_out_spmat_item_row[i] = 0;
	  for(j=0;j<A->row[i].len;j++)
	    {
	      if (A->row[i].elt[j].col<i) // <= before
		{
		  p_out_spmat_item_row[i]++;
		  p_out_spmat_col_pos[Index] = A->row[i].elt[j].col+1;
		  p_out_spmat_val[Index]     = A->row[i].elt[j].val;
		  Index++;
		}
	      else if (A->row[i].elt[j].col==i) // <= before
		{
		  p_out_spmat_item_row[i]++;
		  p_out_spmat_col_pos[Index] = i+1;
		  p_out_spmat_val[Index]     = 1;
		  Index++;
		}
	    }
	}
      
      _SciErr = createSparseMatrix(&_StrCtx,Rhs+1, p_in_spmat_nb_rows, p_in_spmat_nb_cols, Index, 
				   p_out_spmat_item_row, p_out_spmat_col_pos, p_out_spmat_val);

      LhsVar(1) = Rhs+1;
    }

  // Get the U matrix
  if (Lhs>=2)
    {
      Index = 0;
      for(i=0;i<p_in_spmat_nb_rows;i++)
	{
	  p_out_spmat_item_row[i] = 0;
	  for(j=0;j<A->row[i].len;j++)
	    {
	      if (A->row[i].elt[j].col>=i)
		{
		  p_out_spmat_item_row[i]++;
		  p_out_spmat_col_pos[Index] = A->row[i].elt[j].col+1;
		  p_out_spmat_val[Index]     = A->row[i].elt[j].val;
		  Index++;
		}
	    }
	}
      
      _SciErr = createSparseMatrix(&_StrCtx,Rhs+2, p_in_spmat_nb_rows, p_in_spmat_nb_cols, Index, 
				   p_out_spmat_item_row, p_out_spmat_col_pos, p_out_spmat_val);

      LhsVar(2) = Rhs+2;
    }

  // Get the permutation matrix
  if (Lhs==3)
    {
      Index = 0;
      for(i=0;i<p_in_spmat_nb_rows;i++)
	{
	  p_out_spmat_item_row[i] = 1;
	  p_out_spmat_col_pos[i]  = pivot->pe[i]+1;
	  p_out_spmat_val[i]      = 1.0;
	}

      _SciErr = createSparseMatrix(&_StrCtx,Rhs+3, p_in_spmat_nb_rows, p_in_spmat_nb_cols, p_in_spmat_nb_rows, 
				   p_out_spmat_item_row, p_out_spmat_col_pos, p_out_spmat_val);

      LhsVar(3) = Rhs+3;
    }

  if (A) sp_free(A);

  if (p_out_spmat_item_row) FREE(p_out_spmat_item_row);
  if (p_out_spmat_col_pos)  FREE(p_out_spmat_col_pos);
  if (p_out_spmat_val)      FREE(p_out_spmat_val);

  return 0;
}
