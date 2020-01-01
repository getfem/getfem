/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2006-2020 Yves Renard, Julien Pommier.
 
 This file is a part of GetFEM++
 
 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/

#ifndef GFM_COMMON_H
#define GFM_COMMON_H

#include "mex.h"
#include "gfi_array.h"

const char* mxClassID2string(mxClassID id);
int mxarray_to_gfi_array(const mxArray *mx, gfi_array *t);
mxArray* gfi_array_to_mxarray(gfi_array *t);
gfi_array_list *build_gfi_array_list(int nrhs, const mxArray *prhs[]);

typedef void (*getfem_sigint_handler_t)(int);
void install_custom_sigint(getfem_sigint_handler_t h);
void remove_custom_sigint(int allow_rethrow);

#endif
