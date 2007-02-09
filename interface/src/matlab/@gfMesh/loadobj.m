function a=loadobj(b)
% gfMesh/loadobj
% this function is automatically called by matlab when objects of class
% gfMesh are loaded from a MAT-file 

  a=gfMesh('from string',b.txt);
