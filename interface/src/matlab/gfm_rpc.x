enum gfi_type_id {GFI_INT32,GFI_UINT32,GFI_DOUBLE,GFI_CHAR,GFI_CELL,GFI_OBJID,GFI_SPARSE};


struct gfi_object_id {
        int id;
        int cid;
};

struct gfi_sparse {
        int ir<>;
        int jc<>;
        double pr<>;
};

typedef struct gfi_array* pgfi_array;

union gfi_storage switch (gfi_type_id type) {
  case GFI_INT32:
        int data_int32<>;
  case GFI_UINT32:
        unsigned data_uint32<>;
  case GFI_DOUBLE:
        double data_double<>;
  case GFI_CHAR:
        char data_char<>;
  case GFI_CELL:
        pgfi_array data_cell<>;
  case GFI_OBJID:
        struct gfi_object_id objid<>;
  case GFI_SPARSE:
        struct gfi_sparse sp;
};

struct gfi_array {
        unsigned dim<>;
        gfi_storage storage;
};

struct gfi_array_list {
        gfi_array arg<>;
};

enum gfi_status {GFI_STATUS_OK, GFI_STATUS_ERROR};
union gfi_output switch (gfi_status status) {
  case GFI_STATUS_OK:
        gfi_array_list output;
  case GFI_STATUS_ERROR:
        string errmsg<>;
};

program GFMRPC {
      version GFMRPC_VERS_1 {
         void GFMRPC_NULL(void) = 0;
        void GFMRPC_CHDIR(string dir) = 1;        
        gfi_output GFMRPC_CALL(string fname, gfi_array_list in, int nlhs) = 2;
      } = 1;
   } = 400000;
