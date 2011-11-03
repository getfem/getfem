function varargout = %objid_e(varargin)
  gf_obj = varargin(2);
  other_param = varargin(1);
  varargout(1) = [];
  
  // The CID value for each class of object is defined in
  // gfi_array.h
  
  select gf_obj('cid')
  case 0 then
    // gfCvStruct
    varargout = gf_cvstruct_get(gf_obj,other_param);
  case 1 then
    // gfEltm
    // No gf_eltm_get function
  case 2 then
    // gfFem
    varargout = gf_fem_get(gf_obj,other_param);
  case 3 then
    // gfGeoTrans
    varargout = gf_geotrans_get(gf_obj,other_param);
  case 4 then
    // gfGlobalFunction
    varargout = gf_global_function_get(gf_obj,other_param);
  case 5 then
    // gfInteg
    varargout = gf_integ_get(gf_obj,other_param);
  case 6 then
    // gfLevelSet
    varargout = gf_levelset_get(gf_obj,other_param);
  case 7 then
    // gfMdBrick
    varargout = gf_mdbrick_get(gf_obj,other_param);
  case 8 then
    // gfMdState
    varargout = gf_mdstate_get(gf_obj,other_param);
  case 9 then
    // gfMesh
    varargout = gf_mesh_get(gf_obj,other_param);
  case 10 then
    // gfMeshFem
    varargout = gf_mesh_fem_get(gf_obj,other_param);
  case 11 then
    // gfMeshIm
    varargout = gf_mesh_im_get(gf_obj,other_param);
  case 12 then
    // gfMeshLevelSet
    varargout = gf_mesh_levelset_get(gf_obj,other_param);
  case 13 then
    // gfMesherObject
    varargout = gf_mesher_object_get(gf_obj,other_param);
  case 14 then
    // gfModel
    varargout = gf_model_get(gf_obj,other_param);
  case 15 then
    // gfPrecond
    varargout = gf_precond_get(gf_obj,other_param);
  case 16 then
    // gfSlice
    varargout = gf_slice_get(gf_obj,other_param);
  case 17 then
    // gfSpmat
    varargout = gf_spmat_get(gf_obj,other_param);
  case 18 then
    // gfPoly
    // No gf_poly_get function
  else
    error('wrong object ID');
  end
  varargout = list(varargout);
endfunction
