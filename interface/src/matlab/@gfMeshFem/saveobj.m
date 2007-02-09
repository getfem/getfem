function b=saveobj(a)
% gfMeshFem/saveobj
% this function is automatically called by matlab when objects of class
% gfMeshFem are saved in a MAT-file 
  disp('saving gfMeshFem object..');
  b=a; b.txt=char(a);