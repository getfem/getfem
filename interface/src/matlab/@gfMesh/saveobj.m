function b=saveobj(a)
  % gfMesh/saveobj
  % this function is automatically called by matlab when objects of class
  % gfMesh are saved in a MAT-file 
  disp('saving gfMesh object..');
  b=a; b.txt=char(a);
