function a=loadobj(b)
% gfMeshIm/loadobj
% this function is automatically called by matlab when objects of class
% gfMeshIm are loaded from a MAT-file 
  a=gfMeshIm('from string',b.txt);
