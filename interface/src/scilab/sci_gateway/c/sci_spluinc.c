#include <api_common.h>
#include <api_sparse.h>
#include <api_double.h>
#include <MALLOC.h>
#include <stack-c.h>

#include <sparse2.h>
#include <err.h>

#define DEBUG

int sci_spluinc(char * fname)
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
  int      Index, i, j, res;
  int    * p_out_spmat_item_row = NULL;
  int    * p_out_spmat_col_pos  = NULL;
  double * p_out_spmat_val      = NULL;
  double   value, alpha = 1.0;
  int      nnz = 0, var_type;
  StrErr  _StrErr;
  StrCtx  _StrCtx;

  CheckRhs(1,2);
  CheckLhs(1,2);

#ifdef DEBUG
  sciprint("Lhs = %d Rhs = %d\n", Lhs, Rhs);
#endif

  // First, access to the input variable (a matrix of strings)

  _StrErr = getVarAddressFromPosition(&_StrCtx,1,&p_in_spmat_address);

  _StrErr = getVarType(&_StrCtx,p_in_spmat_address,&var_type);

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

  _StrErr = getSparseMatrix(&_StrCtx,p_in_spmat_address, &p_in_spmat_nb_rows, &p_in_spmat_nb_cols, 
			    &p_in_spmat_nb_items, &p_in_spmat_items_row, &p_in_spmat_col_pos, &p_in_spmat_val);

  if (Rhs==2)
    {
      // Second, get the alpha parameter
      // First, access to the input variable (a matrix of doubles)
      _StrErr = getVarAddressFromPosition(&_StrCtx,2,&p_in_dbl_address);
      
      _StrErr = getMatrixOfDouble(&_StrCtx,p_in_dbl_address, &p_in_dbl_nb_rows, &p_in_dbl_nb_cols, &p_in_dbl_matrix);
      alpha = *p_in_dbl_matrix;
    }

  ///////////////////////////////
  // Proceed the factorization //
  ///////////////////////////////

  // Fill the Meschash matrix
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

  // Factorization
  catchall(spILUfactor(A,alpha),Scierror(999,"%s: an error occured.\n",fname); return 0);

  // Now, create the result
  A = sp_col_access(A);
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
	      if (A->row[i].elt[j].col<i)
		{
		  p_out_spmat_item_row[i]++;
		  p_out_spmat_col_pos[Index] = A->row[i].elt[j].col+1;
		  p_out_spmat_val[Index]     = A->row[i].elt[j].val;
		  Index++;
		}
	      else if (A->row[i].elt[j].col==i)
		{
		  p_out_spmat_item_row[i]++;
		  p_out_spmat_col_pos[Index] = i+1;
		  p_out_spmat_val[Index]     = 1;
		  Index++;
		}
	    }
	}
      
      _StrErr = createSparseMatrix(&_StrCtx,Rhs+1, p_in_spmat_nb_rows, p_in_spmat_nb_cols, Index, 
				   p_out_spmat_item_row, p_out_spmat_col_pos, p_out_spmat_val);

      LhsVar(1) = Rhs+1;
    }

  // Get the U matrix
  if (Lhs==2)
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
      
      _StrErr = createSparseMatrix(&_StrCtx,Rhs+2, p_in_spmat_nb_rows, p_in_spmat_nb_cols, Index, 
				   p_out_spmat_item_row, p_out_spmat_col_pos, p_out_spmat_val);

      LhsVar(2) = Rhs+2;
    }

  if (A) sp_free(A);

  if (p_out_spmat_item_row) FREE(p_out_spmat_item_row);
  if (p_out_spmat_col_pos)  FREE(p_out_spmat_col_pos);
  if (p_out_spmat_val)      FREE(p_out_spmat_val);

  return 0;
}
