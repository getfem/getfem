function b=saveobj(a)
% gfMeshIm/saveobj
% this function is automatically called by matlab when objects of class
% gfMeshIm are saved in a MAT-file 
  disp('saving gfMeshIm object..');
  b=a; b.txt=char(a);
