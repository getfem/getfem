function %objid_set(varargin)
  gf_obj = varargin(1);
  other_param = list(varargin(2:$));
  
  // The CID value for each class of object is defined in
  // gfi_array.h

  select gf_obj('cid')
  case 0 then
    // gfContStruct
    // No gf_cont_struct_set function
  case 1 then
    // gfCvStruct
    // No gf_cvstruct_set function
  case 2 then
    // gfEltm
    // No gf_eltm_set function
  case 3 then
    // gfFem
    // No gf_fem_set function
  case 4 then
    // gfGeoTrans
    // No gf_geotrans_set function
  case 5 then
    // gfGlobalFunction
    // No gf_global_function_set function
  case 6 then
    // gfInteg
    // No gf_integ_set function
  case 7 then
    // gfLevelSet
    gf_levelset_set(gf_obj,other_param(:));
  case 8 then
    // gfMesh
    gf_mesh_set(gf_obj,other_param(:));
  case 9 then
    // gfMeshFem
    gf_mesh_fem_set(gf_obj,other_param(:));
  case 10 then
    // gfMeshIm
    gf_mesh_im_set(gf_obj,other_param(:));
  case 11 then
    // gfMeshIm
    gf_mesh_im_data_set(gf_obj,other_param(:));
  case 12 then
    // gfMeshLevelSet
    gf_mesh_levelset_set(gf_obj,other_param(:));
  case 13 then
    // gfMesherObject
    // No gf_mesher_object_set function
  case 14 then
    // gfModel
    gf_model_set(gf_obj,other_param(:));
  case 15 then
    // gfPrecond
    // No gf_precond_set function
  case 16 then
    // gfSlice
    gf_slice_set(gf_obj,other_param(:));
  case 17 then
    // gfSpmat
    gf_spmat_set(gf_obj,other_param(:));
  case 18 then
    // gfPoly
    // No gf_poly_set function
  else
    error('wrong object ID');
  end
endfunction
