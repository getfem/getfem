function [varargout]=gf_solve(varargin)
% function varargout=gf_solve(what, varargin)
% OBSOLETE FUNCTION used in some old scripts. Kept for compatibility reason.
% It solve a few set of pde. DO NOT USE ANYMORE.
  
  if (nargin==0), error('not enough input arguments'); end;
  pde = varargin{1};
  if (~isstruct(pde)),
    error('the first argument should be a pde description structure');
  end;
  
  
  pde = getopt(pde,varargin{2:end});
  
  if (~isfield(pde, 'verbosity')),
    pde.verbosity = 0;
  end;
  if (~isfield(pde, 'mim')),
    error(['since v2.0, the pde structure for gf_solve should contain' ...
	   ' a mesh_im object in its "mim" field']);
  end;
  if (~isfield(pde,'type')),
    error('the pde structure should have a ''type'' field');
  end;
  if (~isfield(pde,'asm')),
    pde.asm = rmfield(struct('dummy',1),'dummy'); % pde.asm = struct; does not work with matlab 6.1.0.450...    
  end;  
  if (~isfield(pde,'solver')),
    pde.solver = 'default';
  end;
  nout=max(nargout,1);
  switch pde.type
   case 'laplacian'
    [varargout{1:nout}]=do_laplacian(pde);
   case 'linear elasticity'
    [varargout{1:nout}]=do_linear_elasticity(pde);
   case 'stokes'
    [varargout{1:nout}]=do_stokes(pde);
   otherwise
    error(['unhandled PDE.type : ' pde.type]);
  end;
  
function assert_field(pde,varargin)
  for i=1:numel(varargin),
    if (~isfield(pde,varargin{i}))
      error(['no member ' varargin{i} ' in struct pde!']); 
    end;
  end;
function ok=has_field(pde,varargin)
  ok = 0;
  for i=1:numel(varargin),
    if (~isfield(pde,varargin{i}))
      return;
    end;
  end;
  ok=1;

  
function pde=eval_asm_data(in_pde,dname,default_value,mf)
  pde = in_pde;
  if (nargin == 3) 
    mf = pde.mf_d;
  end;
  if (~isfield(pde.asm, dname))
    if (isfield(pde,dname)),
      z = getfield(pde,dname);      
    else
      warning(['you did not define the "' dname '" data for the ' pde.type ' pde struct']);
      disp(['setting "' dname '" to its default value of ']); disp(default_value);
      z = default_value;
    end;
    if (iscell(z)) z = reshape(z,numel(z),1); end;
    pde.asm = setfield(pde.asm, dname, gf_mesh_fem_get_eval(mf, z));
  end;

% solves the scalar laplacian
function [U,pde]=do_laplacian(in_pde)
  pde = in_pde; U=[];
  assert_field(pde, 'mf_u','mf_d');
  pde=eval_asm_data(pde,'lambda', {1});
  if (~isfield(pde.asm, 'K')),
    pde.asm.K = gf_asm('laplacian',pde.mim, pde.mf_u, pde.mf_d, pde.asm.lambda);
  end;
  pde = do_classical_bc(pde);
  [U,pde] = do_classical_solve(pde);%.solver,pde.asm.K, pde.asm.Q, pde.asm.G(:), pde.asm.H, pde.asm.R(:), pde.asm.F(:));
% END OF do_laplacian

% solves linear elasticity
function [U,pde]=do_linear_elasticity(in_pde)
  pde = in_pde; U=[];
  assert_field(pde, 'mf_u','mf_d');
  if (~has_field(pde.asm,'lambda','mu'))
    if (has_field(pde,'lambda', 'mu')),
      pde=eval_asm_data(pde,'lambda', {1});
      pde=eval_asm_data(pde,'mu', {1});
    elseif (isfield(pde,'E') & isfield(pde, 'PR')) % young modulus and poisson ratio
      tmpE = gf_mesh_fem_get_eval(pde.mf_d, pde.E);
      tmpnu = gf_mesh_fem_get_eval(pde.mf_d, pde.PR);
      pde.asm.lambda = tmpE.*tmpnu./((1+tmpnu).*(1-2*tmpnu)); 
      pde.asm.mu     = tmpE./(2*(1+tmpnu)); % shear modulus
      % if (is_plane_stress) then
      %   lambda = 2*lambda.*mu./(lambda+2*mu);
      % end;      
    else
      error(['no description of either (young modulus E and poisson '...
	     'ratio nu) or (mu and lambda) in pde structure']);
    end;
  end;
  if (~isfield(pde.asm, 'K')),
    pde.asm.K = gf_asm('linear elasticity',pde.mim,pde.mf_u, pde.mf_d, pde.asm.lambda,pde.asm.mu);
  end;
  pde = do_classical_bc(pde);
  
  %at this poiint, the boundary conditions and volumic source term should have been assembled
  [U,pde] = do_classical_solve(pde);
% END OF do_linear_elasticity

function [U,P,pde]=do_stokes(in_pde)
  pde = in_pde; U=[]; P=[];
  assert_field(pde, 'mf_u','mf_d');
  pde=eval_asm_data(pde,'viscos', {1});
  if (~isfield(pde.asm, 'K')),
    [pde.asm.K,pde.asm.B] = gf_asm('stokes',pde.mim,pde.mf_u, pde.mf_p, pde.mf_d, pde.asm.viscos);
    if (nnz(pde.asm.K-pde.asm.K')),
      error('K not symetric, you found a bug!');
    end;
  end;
  pde = do_classical_bc(pde);
  [U,P,pde] = do_stokes_solve(pde);
% END OF do_stokes
  
function [U,P,pde]=do_stokes_solve(in_pde)
  pde = in_pde; U=[]; P=[];
  assert_field(pde.asm, 'H','R','K','Q','F','G');  
  [null,ud]=gf_spmat_get(pde.asm.H,'dirichlet nullspace', pde.asm.R);
  K=pde.asm.K+pde.asm.Q;
  if nnz(K-K')
    sym=0; disp('non symmetric matrix, aborting; keyboard mode'); keyboard;
  else
    sym=1;
  end    
  Fu=null'*((pde.asm.F(:)+pde.asm.G(:))-K*ud(:));
  Fp=-pde.asm.B'*ud(:);
  K=null'*K*null;
  B=null'*pde.asm.B;
  K=(K+K')/2; % make sure that the matrix is absolutely symetric
	      %  pde.solver.type = 'cg';
	      %  pde.solver = set_default_values(pde.solver,'type','cg','maxiter',1000,'residu',1e-6);
  if (strcmpi(pde.solver,'brute_stokes')),
    [U,P]=do_solve_stokes_cg2(K,B,Fu(:),Fp(:));
  else
    [U,P]=do_solve_stokes_cg(K,B,Fu(:),Fp(:));
  end;
  U=null*U+ud(:);
  U = U(:)';
  P = -P(:)';
% END OF do_stokes_solve

function pde=do_classical_bc(pde)
  q_dim = gf_mesh_fem_get(pde.mf_u, 'qdim');
  do_F = ~isfield(pde.asm,'F');
  do_H = ~isfield(pde.asm,'H');
  do_R = ~isfield(pde.asm,'R');
  do_Q = ~isfield(pde.asm,'Q');
  do_G = ~isfield(pde.asm,'G');
  if (do_F),
    pde=eval_asm_data(pde,'F', num2cell(zeros(q_dim,1)));
    pde.asm.F = gf_asm('volumic source', pde.mim, pde.mf_u, pde.mf_d, ...
                       pde.asm.F);
  end;
  if (isfield(pde, 'pdetool')),
    [pde.asm.Q,pde.asm.G,pde.asm.H,pde.asm.R]=...
	gf_asm('pdetool boundary conditions',...
	       pde.mim,pde.mf_u,pde.mf_d,pde.pdetool.b,pde.pdetool.e);
  else
    assert_field(pde,'bound');
    q_dim = gf_mesh_fem_get(pde.mf_u, 'qdim');
    u_nbdof = gf_mesh_fem_get(pde.mf_u, 'nbdof');
    d_nbdof = gf_mesh_fem_get(pde.mf_d, 'nbdof');
    if (do_H), pde.asm.H = sparse(u_nbdof, u_nbdof); end;
    if (do_Q), pde.asm.Q = sparse(u_nbdof, u_nbdof); end;
    if (do_R), pde.asm.R = zeros(u_nbdof,1); end;
    if (do_G), pde.asm.G = zeros(u_nbdof,1); end;
    for bnum=1:numel(pde.bound),
      assert_field(pde.bound{bnum},'type');
      is_dirichlet = 0; is_neumann = 0;
      switch (pde.bound{bnum}.type)
       case 'None'
       case 'Dirichlet'
	is_dirichlet=1;
       case 'Neumann'
	is_neumann=1;
       case 'Mixed'
	is_dirichlet=1; is_neumann=1;
       otherwise
	disp(['bc type ' pde.bound{bnum}.type 'unhandled']);
      end;
      if (is_dirichlet),
	assert_field(pde.bound{bnum},'R');
	if (do_R | do_H),
          vR = gf_mesh_fem_get_eval(pde.mf_d, {pde.bound{bnum}.R{:}}');
          if (isfield(pde.bound{bnum},'H')),
            vH = gf_mesh_fem_get_eval(pde.mf_d, {pde.bound{bnum}.H{:}}');
          else 
            h = num2cell(eye(q_dim)); vH = gf_mesh_fem_get_eval(pde.mf_d, {h{:}}'); 
          end;
          [bH,bR] = gf_asm('dirichlet', bnum, pde.mim,pde.mf_u, pde.mf_d, reshape(vH,q_dim*q_dim,d_nbdof), vR);
        end;
        if (do_R), pde.asm.R = pde.asm.R + bR; end;
	if (do_H), pde.asm.H = pde.asm.H + bH; end;
      end;
      if (is_neumann),
	assert_field(pde.bound{bnum},'G');
	if (do_G), 
          vG = gf_mesh_fem_get_eval(pde.mf_d, {pde.bound{bnum}.G{:}}');
          vG = gf_asm('boundary source', bnum, pde.mim,pde.mf_u, pde.mf_d, vG);	
          pde.asm.G = pde.asm.G + vG; 
        end;
        if (do_Q),
          if (isfield(pde.bound{bnum},'Q')),
            vQ = gf_mesh_fem_get_eval(pde.mf_d, {pde.bound{bnum}.Q{:}}');
          else 
            q = num2cell(eye(q_dim)); vQ = gf_mesh_fem_get_eval(pde.mf_d, {q{:}}'); 	  
          end;
          bQ = gf_asm('boundary qu term',bnum,pde.mim,pde.mf_u,pde.mf_d, ...
                      reshape(vQ,q_dim*q_dim,d_nbdof));
          pde.asm.Q = pde.asm.Q + bQ;
        end;
      end;
    end;
  end;

%END OF do_classical_bc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   solvers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% solves (K+Q)U=F+G
% under constraint HU=R
function [U,pde]=do_classical_solve(in_pde)
  pde = in_pde;
  assert_field(pde.asm,'K','Q','G','H','R','F');
  [null,ud]=gf_spmat_get(pde.asm.H,'dirichlet nullspace', pde.asm.R);
  RK=pde.asm.K+pde.asm.Q;
  if nnz(RK-RK')
    sym=0; disp('non symmetric matrix');
  else
    sym=1; 
  end
  RF=null'*((pde.asm.F(:)+pde.asm.G(:))-RK*ud(:));
  RK=null'*RK*null;
  if sym
    RK=(RK+RK')/2;
  end
  pde.asm.RK = RK;
  RB=null;
%  gf_matlab('util','save matrix','hb','solve.hb',RK);
  U=RB*(RK\RF)+ud(:);
  U=U(:)'; % row vector
%END OF do_solve


% solves [K  B][U] = [Fu]
%        [B' 0][P]   [Fp]
% with K *positive* definite
function [U,P]=do_solve_stokes_cg(K,B,Fu,Fp)
  verbos_disp_start(sprintf('factorizing K (n=%d,nnz=%d)',size(K,1),nnz(K)));
  R=chol(K);
  verbos_disp_end;
  verbos_disp(sprintf('K factored, nnz(R)=%d',nnz(R)));
% we have to avoid transpositions on sparse matrix since this
% operation has high cost of n*(nnz/n)*log((nnz/n)) ESPECIALY
% for triangular matrices from factorisations: the cost of the
% transposition is greater than the cost of a triangular 
% solve which is n*(nnz/n).
  F=((R\(Fu'/R)')'*B)' - Fp;
  tol=1e-8;
  verbos_disp_start('running Conjugate gradientG');
%  P = cg(F,R,B,10000,1e-6);
  %[P,flag,relres,iter,resvec] = pcg(@multA, F, tol, 500, @multM, @multM, [], R, B);

  [P,flag,relres,iter,resvec] = gmres(@multA, F, 100, tol, 50, @multM, @multM, [], R, B);
  %  figure(5); plot(resvec);
  %  disp(sprintf('    .. flag = %d, relres=%g, iter=%d',flag, relres, ...
  %	       iter));
  
  verbos_disp_end;
  if (flag), 
    warning(sprintf('conjugate gradient did not converge! flag=%d, res=%g, iter=%d',flag,relres,iter));
  else
    verbos_disp(sprintf('pcg: flag=%d, res=%g, iter=%d', flag, relres, iter));
  end;
  U=R\(((Fu-B*P)'/R)');
  verbos_disp('do_solve_stokes_cg all done');


% solves [K  B][U] = [Fu]
%        [B' 0][P]   [Fp]
% with K *positive* definite
function [U,P]=do_solve_stokes_cg3(K,B,Fu,Fp)
  nu = size(K,2); np = size(B,2);

  disp('solve stokes uzawa cholinc');
  [pcB]=cholinc(K,'0');
  pcBt = pcB';
  disp('solve stokes uzawa first pcg');
  P = zeros(np,1);
  U = pcg(K,Fu - B*P,1e-6,100,pcBt,pcB);
  disp('solve stokes uzawa : got U');
  for k=1:10000,
    r = Fp - B'*U;
    res = norm(r);
    if (res < 1e-10) break; end;
    disp(sprintf('solve stokes : iter=%d res=%g',k, res));
    z = pcg(K, B*r, 1e-6, 100, pcBt, pcB);
    rho = res*res/dot(r,(B'*z));
    P = P - rho*r;
    U = U + rho*z;
  end;
						
% try to apply gmres to the global system
function [U,P]=do_solve_stokes_cg2(K,B,Fu,Fp)
  tic;
  nu = size(K,2); np = size(B,2);
  Z = [K B; B' sparse(np,np)];
  Z2=Z + [sparse(nu,nu) sparse(nu,np); sparse(np,nu) spdiags(0.001* ...
						  ones(np,1),0,np, ...
						  np)];
  disp(sprintf('begin luinc [nu=%d,np=%d, nnz=%d]', nu, np, nnz(Z2)));
  [L,U] = luinc(Z2,'0');
  disp('begin gmres');
  [UP,FLAG,RELRES,ITER,RESVEC]=gmres(Z,[Fu;Fp],50,1e-9,1000,L,U);
  U=UP(1:nu); P=UP((nu+1):(nu+np));
  disp(sprintf('do_solve_stokes_cg2 done in %g sec (%d iter, flag=%d)',toc,ITER,FLAG));
  resU=norm(K*U+B*P-Fu,2);
  resP=norm(B'*U-Fp,2);
  disp(sprintf('resU=%g, resP=%g',resU,resP));  

% try to apply gmres to the global system
function [U,P]=do_solve_stokes_cg2_test(K,B,Fu,Fp)
  tic;
  nu = size(K,2); np = size(B,2);
  Z = [K B; B' sparse(np,np)];
  Z2=Z + [sparse(nu,nu) sparse(nu,np); sparse(np,nu) spdiags(0.001* ...
						  ones(np,1),0,np, ...
						  np)];
  disp(sprintf('begin luinc [nu=%d,np=%d, nnz=%d]', nu, np, nnz(Z2)));

  keyboard
  
  [L,U] = luinc(Z2,'0');
  disp('begin gmres');
  [UP,FLAG,RELRES,ITER,RESVEC]=gmres(Z,[Fu;Fp],50,1e-9,1000,L,U);
  U=UP(1:nu); P=UP((nu+1):(nu+np));
  disp(sprintf('do_solve_stokes_cg2 done in %g sec (%d iter, flag=%d)',toc,ITER,FLAG));
  resU=norm(K*U+B*P-Fu,2);
  resP=norm(B'*U-Fp,2);
  disp(sprintf('resU=%g, resP=%g',resU,resP));  
   
function [U,P]=do_solve_stokes_cg2_old(K,B,Fu,Fp)
  alpha=1e-6;
  tic;
  if (0),
    R=chol(K);
    RB=full(R'\B);
    T=(alpha*speye(size(B,2))-RB'*RB); 
    P=T\(Fp-B'*(K\Fu));
    U=R\(((Fu-B*P)'/R)');
  else
    % unfortunately, the basic stokes solver is very slow...
    % on small 3D problems, the fastest way is to reduce to a (full) linear system on the pression...
    % drawback: it eats a lot of memory..
    disp('using the "brute force" solver for stokes..');
    R=chol(K);
    RB=full(R'\B);
    T=(-RB'*RB);
    F=(Fp-B'*(K\Fu));
    T(1,:)=0; T(1,1)=1;F(1)=0;
    P=T\F;
    U=R\(((Fu-B*P)'/R)');
  end;
  disp(sprintf('do_solve_stokes_cg2 done in %g sec',toc));
  resU=norm(K*U+B*P-Fu,2);
  resP=norm(B'*U-Fp,2);
  disp(sprintf('resU=%g, resP=%g',resU,resP));

function Y=multlup(X,L,U,P)
  Y=U\(L\(P*X));
  
% DO NOT USE THIS ONE... BROKEN
function X=cg(F,R,B,maxit,tol)
  X = rand(size(F));
  r = F-multA(X,R,B);
  nr0 = norm(r,2);
  nr = nr0;
  d = r;
  it = 1;
  while (nr/nr0 > tol & it < maxit),
    Ad = multA(d,R,B);
    lambda = (nr^2)/(dot(d, Ad));
    X = X + lambda*d;
    r = r - lambda*Ad;
    nrp = nr;
    nr = norm(r,2);
    beta = (nr*nr)/(nrp*nrp);
    d = r + beta*d;
    it = it+1;
  end;
  disp(sprintf('iterations: %d , res=%g', it, nr/nr0));

function AX=multA(X,R,B)
  tic;
%  AX=B'*(R\(R'\(B*X)));
  BX=(B*X)';
  BXR=(BX/R)';
  RBXR=R\BXR;
  AX=(RBXR'*B)';
%  AX=((R\((B*X)'/R)')'*B)';
  t = toc; verbos_disp(sprintf('iter : %f sec r=%g',t,norm(AX,2)));

function MX=multM(X,R,B)
  MX=X;

function verbos_disp(what)
  global verbosity
  if (verbosity > 0),
    disp(what);
  end;
function verbos_disp_start(what)
  global verbosity
  if (verbosity > 0),
    disp([what '...']); tic;
  end;
function verbos_disp_end
  global verbosity
  if (verbosity > 0)
    disp(sprintf('done (%2.3f sec)', toc));
  end;
    
function r=ison(v)
  r = strcmpi(v,'on');
