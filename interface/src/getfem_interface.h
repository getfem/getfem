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

#ifndef GETFEM_INTERFACE_H
#define GETFEM_INTERFACE_H

#include "gfi_array.h"

#ifdef __cplusplus
extern "C" {
#endif

char* getfem_interface_main(int config_id, const char *function, 
			    int nb_in_args,
			    const gfi_array *in_args[], 
			    int *nb_out_args,
			    gfi_array ***pout_args, char **pinfomsg);

#ifdef __cplusplus
}
#endif

#endif
