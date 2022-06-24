/* -*- c++ -*- (enables emacs c++ mode) */
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

#include <api_scilab.h> 
#include "gfi_array.h"

extern StrCtx* pvApiCtx; // valid for Scilab 6.0 ? 

const char* sci_ClassID2string(sci_types id);
int sci_array_to_gfi_array(int * sci_x, gfi_array *t);
int gfi_array_to_sci_array(gfi_array *t, int i);
gfi_array_list *build_gfi_array_list(int nrhs, int ** prhs);

typedef void (*getfem_sigint_handler_t)(int);
void install_custom_sigint(getfem_sigint_handler_t h);
void remove_custom_sigint(int allow_rethrow);

#endif
