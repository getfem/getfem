function [Pu,Pt,xyP,xyT]=gf_mesh_permuts(p,t)
//[Pu,Pt]=gf_mesh_permuts(p,t)
//[Pu,Pt,xyP,xyT]=gf_mesh_permuts(p,t)

[nargout,nargin] = argn();
    
nodelst=new_intset;
add_to_intset(nodelst,1,size(t,2));
m=new_mesh;
mesh_from_pdetool(m,p,t);
mf_u=new_mesh_fem(m);
fem_u=PK_fem(2,1);
// YC: function ID too long
set_finite_element(mf_u,nodelst,fem_u, ...
		  simplex_exact_integration(2));

mf_p=new_mesh_fem(m);
fem_p=PK_fem(2,0);
set_finite_element(mf_p,nodelst,fem_p, ...
		  simplex_exact_integration(2));

xypts = gf_get_interpolation_pts(mf_u);
Pu=spzeros(size(p,2)*2,size(p,2)*2);
for i=1:size(p,2)
  ii=find((xypts(1,i) == p(1,:)) & (xypts(2,i) == p(2,:)));
  if (length(ii) ~= 1) then pause; end;
  Pu((i-1)*2+1,ii)=1;
  Pu((i-1)*2+2,ii+size(p,2))=1;
end

if (nargout >= 3) then
  xyP = xypts;
end;

xypts = gf_get_interpolation_pts(mf_p);
Pt=spzeros(size(t,2),size(t,2));
for i=1:size(t,2)
  c(:,i)=(p(:,t(1,i))+p(:,t(2,i))+p(:,t(3,i)))/3;
end

for i=1:size(t,2)
//    for j=1,size(t,2)
//      M=[p(1,t(1,j)) p(1,t(2,j)) p(1,t(3,j));...
//	 p(2,t(1,j)) p(2,t(2,j)) p(2,t(3,j));...
//	 1 1 1];
//      MB=[xypts(1,i);xypts(2,i);1];
 
  d=c; 
  d(1,:) = d(1,:)-xypts(1,i);
  d(2,:) = d(2,:)-xypts(2,i);
  d = sum(d.^2);
  [y,j]=min(d);
  if (length(j) ~= 1) then pause; end;
  Pt(i,j)=1;
end

if (nargout >= 4) then
  xyT = xypts;
end
  
del_mesh_fem(mf_u);
del_mesh_fem(mf_p);
del_mesh(m);
del_intset(nodelst);
endfunction

