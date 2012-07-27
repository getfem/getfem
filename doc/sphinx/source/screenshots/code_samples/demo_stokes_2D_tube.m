% this example uses the "old" gf_solve instead of the bricks
% framework..

gf_workspace('clear all');
disp('2D stokes demonstration on a quadratic mesh');
clear pde; 
pde.type = 'stokes';
pde.viscos=1.0;
pde.F = {0,0};
pde.bound{1}.R  = {'-y.*(y-5)',0};
pde.bound{2}.R  = {0,'+(x-20).*(x-25)'};
pde.bound{3}.R  = {0,0};
pde.bound{1}.type = 'Dirichlet';
pde.bound{2}.type = 'Dirichlet';
pde.bound{3}.type = 'Dirichlet';
m=gf_mesh('import','GiD','../meshes/tube_2D_spline.GiD.msh');
pde.mf_u=gf_mesh_fem(m,2);
mfulag=gf_mesh_fem(m,2);
pde.mf_p=gf_mesh_fem(m,1);
pde.mf_d=gf_mesh_fem(m,1);
pde.mim=gf_mesh_im(m,gf_integ('IM_TRIANGLE(5)'));
% this is a good example of the usefulness of the cubic bubble
% -> if not used, the pressure has strange values
gf_mesh_fem_set(pde.mf_u,'fem',gf_fem('FEM_PK_WITH_CUBIC_BUBBLE(2,2)'));
gf_mesh_fem_set(pde.mf_d,'fem',gf_fem('FEM_PK(2,2)'));
gf_mesh_fem_set(pde.mf_p,'fem',gf_fem('FEM_PK_DISCONTINUOUS(2,1)'));

% we use a P3 mesh fem for interpolation of the U field, since
% because of its cubic bubble function, the pde.mf_u is not lagrangian 
gf_mesh_fem_set(mfulag,'fem',gf_fem('FEM_PK(2,3)'));

all_faces = gf_mesh_get(m, 'outer faces', gf_mesh_get(m, 'cvid'));
P=gf_mesh_get(m,'pts');
INpid=find(abs(P(1,:)) < 1e-4);
OUTpid=find(abs(P(2,:)+20) < 1e-4);
INfaces=gf_mesh_get(m, 'faces from pid', INpid);
OUTfaces=gf_mesh_get(m, 'faces from pid', OUTpid);
gf_mesh_set(m, 'boundary', 1, INfaces);
gf_mesh_set(m, 'boundary', 2, OUTfaces);
gf_mesh_set(m, 'boundary', 3, setdiff(all_faces',union(INfaces',OUTfaces','rows'),'rows')');

tic; [U,P]=gf_solve(pde); disp(sprintf('solve done in %.2f sec', toc));

Ul=gf_compute(pde.mf_u,U,'interpolate on',mfulag);
subplot(2,2,1); 
gf_plot(mfulag,Ul,'norm','on','deformation',Ul,'deformation_scale','10%',...
	'deformed_mesh','on');
colorbar; title('|U| plotted on the deformed mesh');

subplot(2,2,2); 
gf_plot(pde.mf_p,P(:)','deformation',U,'deformation_mf',pde.mf_u); 
colorbar; title('Pression on the deformed mesh');

subplot(2,2,3); 
gf_plot(mfulag,Ul(:)','mesh','on','meshopts',{});
hold on; gf_plot(pde.mf_p,P(:)','refine',1); hold off; 
colorbar; title('Quiver plot of U, with color plot of the pression');

subplot(2,2,4); 
gf_plot(mfulag,Ul(:)','mesh','on','meshopts',{}, ...
	'quiver_density',100,'quiver_scale',0.4); 
hold on; gf_plot(pde.mf_p,P(:)'); 
axis([27 33 3 9]); title('Quiver plot zoomed');
