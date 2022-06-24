/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2006-2020 Yves Renard, Julien Pommier.

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
