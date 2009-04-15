% clears every getfem object 
gf_workspace('clear all');
clear pde;
disp('getfem_matlab basic example -- 1D laplacian with high order FEM');
disp(['this demo illustrates some current drawbacks of this getfem ' ...
      'release. Hopefully, these numerical problems with high order elements ' ...
      'will be solved soon!']);
% creation of a simple cartesian mesh
m = gf_mesh('cartesian',[0:1:100]);

% visualisation of the mesh, with the numbering
% of vertices and convexes.
gf_plot_mesh(m, 'vertices','convexes');
%disp('press enter to continue..'); pause


mfu = gfMeshFem(m,1);% set(mfu, 'classical fem', 0);
mfd = gfMeshFem(m,1);% set(mfd, 'classical fem', 0);
mim = gfMeshIm(m);
expr_u = 'sin(x/2)+x/100'; %sin(x).*cos(y)';
expr_f = 'sin(x/2)/4'; %'-2*sin(x).*sin(y)';

mferr=gfMeshFem(m,1); 
gf_mesh_fem_set(mferr,'fem',gf_fem('FEM_PK(1,4)'));
x = gf_mesh_fem_get(mferr, 'basic dof nodes');
eU = eval(expr_u);

bord = gf_mesh_get(m,'outer faces');
gf_mesh_set(m, 'boundary', 1, bord);

    
b0=gfMdBrick('generic elliptic',mim,mfu);
b1=gfMdBrick('dirichlet', b0, 1, mfu,'augmented');
b2=gfMdBrick('source term',b1);


gf_util('trace level', 0); % do not clutter the output with
                           % boring messages..
gf_util('warning level', 0);

for j=1:3
  if (j==1),
    disp(['Results for the basic PK FEM.. as you can see the ' ...
          'numerical noise of PK(K>4) has an impact on our ' ...
          'evaluation of Dirichlet conditions..']);
  elseif (j==2),
    disp(['Results for the basic PK FEM. Quite better but it is well-known that ' ...
          'polynomial basis of high order are quite bad for ' ...
          'interpolation: this has an impact on the volumic source']);
  else
    disp(['Results with the PK_HIERARCHICAL. Actually the best, since the ' ...
	  'hierarchical basis is much better conditionned than the ' ...
	  'PK. Note that the condition number of the matrix ' ...
	  'increases slower']);
  end;
  for iK=1:6,
    if (j==3), 
      K = 2*iK; % avoid prime numbers: on primes numbers, the hierarchical basic is the PK basis
    else K=iK; end;
    if (j<=2),
      fem_u=gf_fem(sprintf('FEM_PK(1,%d)',K));
      fem_d=gf_fem(sprintf('FEM_PK(1,%d)',K));
    else
      fem_u=gf_fem(sprintf('FEM_PK_HIERARCHICAL(1,%d)',K));
      fem_d=gf_fem(sprintf('FEM_PK(1,%d)',K));
    end
    
    im = gf_integ('IM_EXACT_SIMPLEX(1)');
    gf_mesh_fem_set(mfu,'fem',fem_u);
    gf_mesh_fem_set(mfd,'fem',fem_d);
    gf_mesh_im_set(mim, 'integ', im);
 
    set(b1, 'param', 'R', mfd, get(mfd, 'eval', { expr_u }));
    set(b2, 'param', 'source_term', mfd, get(mfd, 'eval', { expr_f }));
    
    mds = gfMdState(b2);
    
    gf_mdbrick_get(b2, 'solve', mds);
    Uu=gf_mdstate_get(mds, 'state'); Uu=Uu(1:get(mfu, 'nbdof'));

    %gf_plot(pde.mf_u, U,'mesh','regions'); colorbar;
    U = gf_compute(mfu, Uu, 'interpolate on', mferr);
    plot(x,U,'r+-',x,eU,'bx:'); legend('approx','exact');
    [mx,pos]=max(abs(eU-U));
    disp(sprintf('K=%d .. max_rel(err)=%1.5g at x=%3.2f [condition number=%e]',K,mx/max(abs(eU)), x(pos),condest(get(mds, 'reduced_tangent_matrix'))));
  end;
end;
