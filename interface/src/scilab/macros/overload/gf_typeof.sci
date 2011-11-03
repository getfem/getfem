function res = gf_typeof(gf_var)
  res = '';

  // The CID value for each class of object is defined in
  // gfi_array.h

  if (typeof(gf_var)~='objid') then
    error('gf_typeof: only objid structures accepted');
  end

  select gf_var('cid')
  case 0 then
    res = 'gfCvStruct';
  case 1 then
    res =  'gfEltm';
  case 2 then
    res = 'gfFem';
  case 3 then
    res = 'gfGeoTrans';
  case 4 then
    res = 'gfGlobalFunction';
  case 5 then
    res = 'gfInteg';
  case 6 then
    res = 'gfLevelSet';
  case 7 then
    res = 'gfMdBrick';
  case 8 then
    res = 'gfMdState';
  case 9 then
    res = 'gfMesh';
  case 10 then
    res = 'gfMeshFem';
  case 11 then
    res = 'gfMeshIm';
  case 12 then
    res = 'gfMeshLevelSet';
  case 13 then
    res = 'gfMesherObject';
  case 14 then
    res = 'gfModel';
  case 15 then
    res = 'gfPrecond';
  case 16 then
    res = 'gfSlice';
  case 17 then
    res = 'gfSpmat';
  case 18 then
    res = 'gfPoly';
  else
    error('wrong object ID');
  end
endfunction
