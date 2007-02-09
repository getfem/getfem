function a=loadobj(b)
% gfMeshFem/loadobj
% this function is automatically called by matlab when objects of class
% gfMeshFem are loaded from a MAT-file 
  a=gfMeshFem('from string',b.txt);
