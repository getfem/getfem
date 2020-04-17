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

#include <gfm_common.h>
#include <sciprint.h>
#include <Scierror.h>

#include <sparse2.h>
#include <err.h>

//#define DEBUG

int sci_spcholinc(char * fname)
{
  int      p_in_spmat_nb_rows, p_in_spmat_nb_cols, p_in_spmat_nb_items;
  int    * p_in_spmat_address;
  int    * p_in_spmat_items_row = NULL;
  int    * p_in_spmat_col_pos   = NULL;
  double * p_in_spmat_val       = NULL;
  SPMAT  * A = NULL;
  int      Index, i, j;
  int    * p_out_spmat_item_row = NULL;
  int    * p_out_spmat_col_pos  = NULL;
  double * p_out_spmat_val      = NULL;
  int      nnz = 0, var_type;
  SciErr _SciErr;

  CheckRhs(1,1);
  CheckLhs(1,1);

  // First, access to the input variable (a matrix of strings)
    
  _SciErr = getVarAddressFromPosition(pvApiCtx,1,&p_in_spmat_address);

  _SciErr = getVarType(pvApiCtx,p_in_spmat_address,&var_type);
  if (var_type!=sci_sparse)
    {
      Scierror(999,"%s: wrong parameter, a sparse matrix is needed\n",fname);
      return 0;
    }

  if (isVarComplex(pvApiCtx,p_in_spmat_address))
    {
      Scierror(999,"%s: wrong parameter, a real sparse matrix is needed\n",fname);
      return 0;
    }

  _SciErr = getSparseMatrix(pvApiCtx,p_in_spmat_address, &p_in_spmat_nb_rows, &p_in_spmat_nb_cols, 
			    &p_in_spmat_nb_items, &p_in_spmat_items_row, &p_in_spmat_col_pos, &p_in_spmat_val);

#ifdef DEBUG
  sciprint("DEBUG: %d, %d\n",p_in_spmat_nb_rows, p_in_spmat_nb_cols);
#endif

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
  catchall(spICHfactor(A),Scierror(999,"%s: an error occured.\n",fname); return 0);

  // Now, create the result
  A = sp_col_access(A);
  for(i=0;i<A->m; i++) nnz += A->row[i].len;

  p_out_spmat_item_row = (int *)MALLOC(p_in_spmat_nb_rows*sizeof(int));
  p_out_spmat_col_pos  = (int *)MALLOC(nnz*sizeof(int));
  p_out_spmat_val      = (double *)MALLOC(nnz*sizeof(double));

  // Get the L matrix
  Index = 0;
  for(i=0;i<p_in_spmat_nb_rows;i++)
    {
      p_out_spmat_item_row[i] = 0;
      for(j=0;j<A->row[i].len;j++)
	{
	  if (A->row[i].elt[j].col<=i)
	    {
	      p_out_spmat_item_row[i]++;
	      p_out_spmat_col_pos[Index] = A->row[i].elt[j].col+1;
	      p_out_spmat_val[Index]     = A->row[i].elt[j].val;
	      Index++;
	    }
	}
    }
  
  _SciErr = createSparseMatrix(pvApiCtx,Rhs+1, p_in_spmat_nb_rows, p_in_spmat_nb_cols, Index, 
			       p_out_spmat_item_row, p_out_spmat_col_pos, p_out_spmat_val);
  
  LhsVar(1) = Rhs+1;

  if (A) sp_free(A);

  if (p_out_spmat_item_row) FREE(p_out_spmat_item_row);
  if (p_out_spmat_col_pos)  FREE(p_out_spmat_col_pos);
  if (p_out_spmat_val)      FREE(p_out_spmat_val);

  return 0;
}
