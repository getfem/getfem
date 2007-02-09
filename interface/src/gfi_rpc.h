/* -*- c++ -*- (enables emacs c++ mode)
   ========================================================================
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

========================================================================
*/


#ifndef GFI_RPC_H
#define GFI_RPC_H

#include <rpc/rpc.h>
#include <gfi_array.h>

#ifdef __cplusplus
extern "C" {
#endif


struct gfmrpc_call_1_argument {
  int config_id;
  char *fname;
  gfi_array_list in;
  int nlhs;
};
typedef struct gfmrpc_call_1_argument gfmrpc_call_1_argument;

#define GFMRPC 400000
#define GFMRPC_VERS_1 1

#define GFMRPC_NULL 0
extern  void * gfmrpc_null_1(CLIENT *);
extern  void * gfmrpc_null_1_svc(struct svc_req *);
#define GFMRPC_CHDIR 1
extern  void * gfmrpc_chdir_1(char *, CLIENT *);
extern  void * gfmrpc_chdir_1_svc(char *, struct svc_req *);
#define GFMRPC_CALL 2
extern  gfi_output * gfmrpc_call_1(int, char *, gfi_array_list , int , CLIENT *);
extern  gfi_output * gfmrpc_call_1_svc(int, char *, gfi_array_list , int , struct svc_req *);
extern int gfmrpc_1_freeresult (SVCXPRT *, xdrproc_t, caddr_t);

/* the xdr functions */

extern  bool_t xdr_gfi_type_id (XDR *, gfi_type_id*);
extern  bool_t xdr_gfi_object_id (XDR *, gfi_object_id*);
extern  bool_t xdr_gfi_sparse (XDR *, gfi_sparse*);
extern  bool_t xdr_pgfi_array (XDR *, pgfi_array*);
extern  bool_t xdr_gfi_storage (XDR *, gfi_storage*);
extern  bool_t xdr_gfi_array (XDR *, gfi_array*);
extern  bool_t xdr_gfi_array_list (XDR *, gfi_array_list*);
extern  bool_t xdr_gfi_status (XDR *, gfi_status*);
extern  bool_t xdr_gfi_output (XDR *, gfi_output*);
extern  bool_t xdr_gfmrpc_call_1_argument (XDR *, gfmrpc_call_1_argument*);

#ifdef __cplusplus
}
#endif

#endif /* !GFI_RPC_H */
