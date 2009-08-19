#include <mex.h> 
#include <sci_gateway.h>
static int direct_gateway(char *fname,void F(void)) { F();return 0;};
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
extern Gatefunc sci_gf_scilab;
static GenericTable Tab[]={
  {(Myinterfun)sci_gateway,sci_gf_scilab,"util"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"cvstruct_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"geotrans"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"geotrans_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"compute"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mesh_fem"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mesh_fem_set"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mesh_fem_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mesh_im"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mesh_im_set"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mesh_im_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"eltm"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mesh"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mesh_set"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mesh_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mdbrick"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mdbrick_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mdbrick_set"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mdstate"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mdstate_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mdstate_set"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"model"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"model_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"model_set"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"slice"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"slice_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"slice_set"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"levelset"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"levelset_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"levelset_set"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mesh_levelset"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mesh_levelset_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"mesh_levelset_set"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"precond"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"precond_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"linsolve"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"spmat"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"spmat_set"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"spmat_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"asm"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"fem"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"fem_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"integ"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"integ_get"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"workspace"},
  {(Myinterfun)sci_gateway,sci_gf_scilab,"delete"},
};
 
int C2F(libscigetfem_c)()
{
  Rhs = Max(0, Rhs);
  (*(Tab[Fin-1].f))(Tab[Fin-1].name,Tab[Fin-1].F);
  return 0;
}
