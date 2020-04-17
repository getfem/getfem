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
#include <Scierror.h>
#include <sciprint.h>

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
  double * p_out_x_dbl_matrix  = NULL;
  SPMAT  * A  = NULL;
  VEC    * vB = NULL, * vOut = NULL;
  int      Index, i, j;
  SciErr   _SciErr;
  int      var_type;

  CheckRhs(1,2);
  CheckLhs(1,1);

  // First, access to the input variable (a matrix of strings)
  _SciErr = getVarAddressFromPosition(pvApiCtx, 1,&p_in_spmat_address);

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

  _SciErr = getSparseMatrix(pvApiCtx, p_in_spmat_address, &p_in_spmat_nb_rows, &p_in_spmat_nb_cols, 
			&p_in_spmat_nb_items, &p_in_spmat_items_row, &p_in_spmat_col_pos, &p_in_spmat_val);

  // Second, get b
  _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_in_b_dbl_address);

  _SciErr = getMatrixOfDouble(pvApiCtx, p_in_b_dbl_address, &p_in_b_nb_rows, &p_in_b_nb_cols, &p_in_b_dbl_matrix);

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

  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+1, p_in_b_nb_rows, p_in_b_nb_cols, p_out_x_dbl_matrix);

  LhsVar(1) = Rhs+1;

  if (A) sp_free(A);
  if (p_out_x_dbl_matrix) FREE(p_out_x_dbl_matrix);

  return 0;
}
