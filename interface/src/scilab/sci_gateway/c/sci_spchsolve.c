#include <api_common.h>
#include <api_sparse.h>
#include <api_double.h>
#include <MALLOC.h>
#include <stack-c.h>

#include <string.h>

#include <sparse2.h>
#include <err.h>

int sci_spchsolve(char * fname)
{
  int      p_in_spmat_nb_rows, p_in_spmat_nb_cols, p_in_spmat_nb_items;
  int    * p_in_spmat_address;
  int    * p_in_spmat_items_row = NULL;
  int    * p_in_spmat_col_pos   = NULL;
  double * p_in_spmat_val       = NULL;
  int      p_in_b_nb_rows, p_in_b_nb_cols;
  double * p_in_b_dbl_matrix  = NULL;
  int    * p_in_b_dbl_address = NULL;
  int      p_out_x_nb_rows, p_out_x_nb_cols;
  double * p_out_x_dbl_matrix  = NULL;
  int    * p_out_x_dbl_address = NULL;
  SPMAT  * A  = NULL;
  VEC    * vB = NULL, * vOut = NULL;
  int      Index, i, j, res;
  double   value, alpha = 1.0;
  StrErr   _StrErr;
  StrCtx   _StrCtx;
  int      nnz = 0, var_type;

  CheckRhs(1,2);
  CheckLhs(1,1);

  // First, access to the input variable (a matrix of strings)
  _StrErr = getVarAddressFromPosition(&_StrCtx, 1,&p_in_spmat_address);

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

  _StrErr = getSparseMatrix(&_StrCtx, p_in_spmat_address, &p_in_spmat_nb_rows, &p_in_spmat_nb_cols, 
			&p_in_spmat_nb_items, &p_in_spmat_items_row, &p_in_spmat_col_pos, &p_in_spmat_val);

  // Second, get b
  _StrErr = getVarAddressFromPosition(&_StrCtx, 2, &p_in_b_dbl_address);

  _StrErr = getMatrixOfDouble(&_StrCtx, p_in_b_dbl_address, &p_in_b_nb_rows, &p_in_b_nb_cols, &p_in_b_dbl_matrix);

  ////////////////////////////
  // Proceed the resolution //
  ////////////////////////////

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

  // Fill the Meschash vector
  vB   = v_get(p_in_b_nb_rows);
  vOut = v_get(p_in_b_nb_rows);
  for(i=0;i<p_in_b_nb_rows;i++)
    {
      v_set_val(vB,i,p_in_b_dbl_matrix[i]);
    }

  // Resolution
  catchall(spCHsolve(A,vB,vOut),Scierror(999,"%s: an error (%d) occured.\n",fname,_err_num); return 0);

  // Now, create the result
  p_out_x_dbl_matrix = (double *)MALLOC(p_in_b_nb_rows*sizeof(double));
  memcpy(p_out_x_dbl_matrix,vOut->ve,p_in_b_nb_rows*sizeof(double));

  _StrErr = createMatrixOfDouble(&_StrCtx, Rhs+1, p_in_b_nb_rows, p_in_b_nb_cols, p_out_x_dbl_matrix);

  LhsVar(1) = Rhs+1;

  if (A) sp_free(A);
  if (p_out_x_dbl_matrix) FREE(p_out_x_dbl_matrix);

  return 0;
}
