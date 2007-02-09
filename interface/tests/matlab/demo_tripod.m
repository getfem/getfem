disp('This demo is an adaption of the original tripod demo')
disp('which uses the new "brick" framework of getfem')
disp('The code is shorter, faster and much more powerful')
disp('You can easily switch between linear/non linear')
disp('compressible/incompressible elasticity!')

linear = 1
incompressible = 0


gf_workspace('clear all');
% import the mesh
m=gfMesh('import','gid','../meshes/tripod.GiD.msh');
mfu=gfMeshFem(m,3);     % mesh-fem supporting a 3D-vector field
mfd=gfMeshFem(m,1);     % scalar mesh_fem, for data fields.
% the mesh_im stores the integration methods for each tetrahedron
mim=gfMeshIm(m,gf_integ('IM_TETRAHEDRON(5)'));
% we choose a P2 fem for the main unknown
gf_mesh_fem_set(mfu,'fem',gf_fem('FEM_PK(3,2)'));
% the material is homogeneous, hence we use a P0 fem for the data
gf_mesh_fem_set(mfd,'fem',gf_fem('FEM_PK(3,0)'));
% display some informations about the mesh
disp(sprintf('nbcvs=%d, nbpts=%d, nbdof=%d',gf_mesh_get(m,'nbcvs'),...
             gf_mesh_get(m,'nbpts'),gf_mesh_fem_get(mfu,'nbdof')));
P=gf_mesh_get(m,'pts'); % get list of mesh points coordinates
pidtop=find(abs(P(2,:)-13)<1e-6); % find those on top of the object
pidbot=find(abs(P(2,:)+10)<1e-6); % find those on the bottom
% build the list of faces from the list of points
ftop=gf_mesh_get(m,'faces from pid',pidtop); 
fbot=gf_mesh_get(m,'faces from pid',pidbot);
% assign boundary numbers
gf_mesh_set(m,'boundary',1,ftop);
gf_mesh_set(m,'boundary',2,fbot);

E = 1e3; Nu = 0.3;
% set the Lame coefficients
lambda = E*Nu/((1+Nu)*(1-2*Nu));
mu = E/(2*(1+Nu));

% create a meshfem for the pressure field (used if incompressible ~= 0)
mfp=gfMeshFem(m); set(mfp, 'fem',gfFem('FEM_PK_DISCONTINUOUS(3,0)'));
if (linear)
  % the linearized elasticity , for small displacements
  b0 = gfMdBrick('isotropic_linearized_elasticity',mim,mfu)
  set(b0, 'param','lambda', lambda);
  set(b0, 'param','mu', mu);
  if (incompressible)
    b1 = gfMdBrick('linear incompressibility term', b0, mfp);
  else
    b1 = b0;
  end;
else
  % See also demo_nonlinear_elasticity for a better example
  if (incompressible)
    b0 = gfMdBrick('nonlinear elasticity',mim, mfu, 'Mooney Rivlin');
    b1 = gfMdBrick('nonlinear elasticity incompressibility term',b0,mfp);
    set(b0, 'param','params',[lambda;mu]);
  else
    % large deformation with a linearized material law.. not
    % a very good choice!
    b0 = gfMdBrick('nonlinear elasticity',mim, mfu, 'SaintVenant Kirchhoff');
    set(b0, 'param','params',[lambda;mu]);
    %b0 = gfMdBrick('nonlinear elasticity',mim, mfu, 'Ciarlet Geymonat');
    b1 = b0;
  end;
end

% set a vertical force on the top of the tripod
b2 = gfMdBrick('source term', b1, 1);
set(b2, 'param', 'source_term', mfd, get(mfd, 'eval', {0;-10;0}));

% attach the tripod to the ground
b3 = gfMdBrick('dirichlet', b2, 2, mfu, 'penalized');

mds=gfMdState(b3)

disp('running solve...')

t0=cputime; 

get(b3, 'solve', mds, 'noisy', 'max_iter', 1000, 'max_res', 1e-6, 'lsolver', 'superlu');
disp(sprintf('solve done in %.2f sec', cputime-t0));

mfdu=gf_mesh_fem(m,1);
% the P2 fem is not derivable across elements, hence we use a discontinuous
% fem for the derivative of U.
gf_mesh_fem_set(mfdu,'fem',gf_fem('FEM_PK_DISCONTINUOUS(3,1)'));
VM=get(b0, 'von mises',mds,mfdu);

U=get(mds, 'state'); U=U(1:get(mfu, 'nbdof'));

disp('plotting ... can also take some minutes!');

% we plot the von mises on the deformed object, in superposition
% with the initial mesh.
if (linear),
  gf_plot(mfdu,VM,'mesh','on', 'cvlst', get(m, 'outer faces'),...
	  'deformation',U,'deformation_mf',mfu);
else
  gf_plot(mfdu,VM,'mesh','on', 'cvlst', get(m, 'outer faces'),...
	  'deformation',U,'deformation_mf',mfu,'deformation_scale',1);
end;

caxis([0 100]);
colorbar; view(180,-50); camlight;
gf_colormap('tripod');

% the von mises stress is exported into a VTK file
% (which can be viewed with 'mayavi -d tripod.vtk -m BandedSurfaceMap')
% see http://mayavi.sourceforge.net/
gf_mesh_fem_get(mfdu,'export to vtk','tripod.vtk','ascii',VM,'vm')

