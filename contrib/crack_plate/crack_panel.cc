/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2002-2012 Yves Renard, Julien Pommier.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/

/**
 * Problem dealing with a problem related to an industrial situation.
 * 
*/
#include "crack_bilaplacian.h"

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_linearized_plates.h"
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "gmm/gmm.h"
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_mesh_fem_product.h"
#include "getfem/getfem_mesh_fem_global_function.h"
#include "getfem/getfem_mesh_fem_sum.h"

//#include "../tests/crack.cc"




/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::short_type;  
using bgeot::dim_type;
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector; 




/************************************************************
 * main program
 ************************************************************/
 
 
int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.  

  try {
    bilaplacian_crack_problem flex_pb ;
    flex_pb.PARAM.read_command_line(argc, argv);
    flex_pb.init() ;
    plain_vector U;
    scalar_type ring_radius = flex_pb.PARAM.real_value("RING_RADIUS");
    if (flex_pb.PARAM.int_value("SOL_REF") == 2) {
       if (!flex_pb.solve(U)) GMM_ASSERT1(false, "Solve has failed");
       cout.precision(16);
       flex_pb.compute_sif(U, ring_radius);
    }
    if (p.PARAM.int_value("ENRICHMENT_OPTION") > 2){
        p.sif_direct_estimation(U) ;
    }
    
    
    // Export solutions for visualisation (bending pb only, for the moment)
    int VTK_EXPORT = int(p.PARAM.int_value("VTK_EXPORT"));
    int MATLAB_EXPORT = int(p.PARAM.int_value("MATLAB_EXPORT"));
    int DX_EXPORT = int(p.PARAM.int_value("DX_EXPORT"));
    if (VTK_EXPORT || MATLAB_EXPORT || DX_EXPORT){
       flex_pb.export_solution(U) ;
    }
    //crack_problem memb_pb ;
    cout << "fin du programme atteinte \n" ;
  }

      GMM_STANDARD_CATCH_ERROR;
  
  return 0; 
}



















