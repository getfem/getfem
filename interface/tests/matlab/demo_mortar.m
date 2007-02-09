% do a partition of the mesh into two disjoint areas, and then
% solve the linear elasticity problem with a mortar join on 
% the interface between the two areas

gf_workspace('clear all'); 
NX=9;
radius = 0.3; xc = .5; yc = .5;
m=gfMesh('cartesian', 0:1/NX:1, 0:1/NX:1);
[pid,idx] = get(m, 'pid_from_cvid');
P=get(m,'pts');

areap=[];
for cv=1:(numel(idx)-1),
  areap(cv) = 1;
  for i=idx(cv):(idx(cv+1)-1),
    if (norm(P(:,pid(i)) - [xc;yc]) > radius),
      areap(cv)=0;
      break;
    end;
  end;
end;

mfu=gfMeshFem(m, 2); set(mfu, 'fem', gfFem('FEM_QK(2,2)'));
mfd=gfMeshFem(m, 1); set(mfd, 'fem', gfFem('FEM_QK(2,1)'));
mfm=gfMeshFem(m, 2); set(mfm, 'fem', gfFem('FEM_QK(2,2)'));
mfdu=gfMeshFem(m);  set(mfdu,'fem',gfFem('FEM_QK_DISCONTINUOUS(2,2)'));
mim=gfMeshIm(m, 5);




set(mfu, 'dof_partition', areap);

b_in  = get(m, 'outer faces', find(areap==1));
b_out = get(m, 'outer faces', find(areap==0));
b_border = get(m, 'outer faces');
b_out = setdiff(b_out', b_border', 'rows')';

fleft =gf_mesh_get(m,'faces from pid',find(abs(P(1,:))<1e-6));
fright=gf_mesh_get(m,'faces from pid',find(abs(P(1,:) - 1)<1e-6));
% assign boundary numbers
set(m,'region',1,fleft);
set(m,'region',2,fright);

MORTAR_BOUNDARY_IN = 40;
MORTAR_BOUNDARY_OUT = 41;
set(m, 'region', MORTAR_BOUNDARY_IN, b_in);
set(m, 'region', MORTAR_BOUNDARY_OUT, b_out);


gf_plot_mesh(m,'boundaries',40);
disp('This is the mortar interface (press a key to continue)'); pause;

indm=get(mfm, 'dof on region', MORTAR_BOUNDARY_OUT);
expr = 'M(#1,#2)+=comp(vBase(#1).vBase(#2))(:,i,:,i)';
M =   gf_asm('boundary', MORTAR_BOUNDARY_IN , expr, mim, mfm, mfu);
M = M-gf_asm('boundary', MORTAR_BOUNDARY_OUT, expr, mim, mfm, mfu);
M=M(indm, :);

b0=gfMdBrick('isotropic_linearized_elasticity',mim, mfu);

%b0=gfMdBrick('generic elliptic',mim, mfu);
b1=gfMdBrick('dirichlet',b0, 1, mfu,'augmented');
b2=gfMdBrick('source term',b1, 2);
F=get(mfd, 'eval', {0; 'y+5'});
set(b2, 'param','source_term', mfd, F);

b3=gfMdBrick('constraint', b2, 'augmented');
set(b3, 'constraints', M, zeros(numel(indm),1));

bn = b3;

mds=gfMdState(bn);
get(bn, 'solve', mds, 'very noisy');
U=get(mds, 'state');U=U(1:get(mfu, 'nbdof'));
VM = get(b0, 'von mises', mds, mfdu);
gf_plot(mfdu,VM,'deformed_mesh','on', 'deformation',U,...
	'deformation_mf',mfu,'refine', 4, 'deformation_scale',1); 
caxis([0 500]);