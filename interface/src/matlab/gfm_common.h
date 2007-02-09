/* -*- c++ -*- (enables emacs c++ mode) */
/*========================================================================

 Copyright (C) 2006-2006 Yves Renard, Julien Pommier.

 This file is a part of GETFEM++

 Getfem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; version 2.1 of the License.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 You should have received a copy of the GNU Lesser General Public
 License along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
 USA.

 ========================================================================*/

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
