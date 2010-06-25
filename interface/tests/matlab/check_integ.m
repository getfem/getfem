function check_integ(iverbose,idebug)
  global gverbose;
  global gdebug;  
  if (nargin >= 1),
    gverbose = iverbose;
    if (nargin == 2),
      gdebug = idebug;
    else 
      gdebug = 0; end;
  else 
    gverbose = 0;
  end;
  gf_workspace('clear all');
  im=gf_integ('IM_TRIANGLE(3)');
  dim = gf_integ_get(im,'dim');  
  gfassert('dim==2');
  ise = gf_integ_get(im,'is_exact');
  gfassert('~ise');
  npt = gf_integ_get(im,'nbpts');
  gfassert('npt==[4 2 2 2]');
  pts = gf_integ_get(im,'pts');
  c = gf_integ_get(im,'coeffs');
  gfassert('size(pts)==[2 10]');
  gfassert('size(c)==[1 10]');
  im=gf_integ('IM_TRIANGLE(7)');  
  c = gf_integ_get(im,'coeffs');
  C=[0.0267    0.0267    0.0267    0.0386    0.0386    0.0386 ...
     0.0386    0.0386    0.0386    0.0878    0.0878    0.0878 ...
    -0.0748    0.2460    0.2460    0.4611    0.4611    0.1739 ...
     0.1739    0.3261    0.3261    0.1739    0.1739    0.3261    0.3261];
%   C=[0.0386    0.0386    0.0267    0.0267    0.0878    0.0878 ...
%      0.0386    0.0386   -0.0748    0.0878    0.0386    0.0386 ...
%      0.0267    0.2460    0.4611    0.4611    0.2460    0.1739 ...
%      0.3261    0.3261    0.1739    0.1739    0.3261    0.3261    0.1739];
  gfassert('norm(c(:)-C(:))<1e-3');
  for i=-1:4,
    if (i >= 1 & i <= 3)
      gf_integ_get(im,'face_pts',i);
      gf_integ_get(im,'face_coeffs',i);
    else
      asserterr('gf_integ_get(im,''face_pts'',i)');
      asserterr('gf_integ_get(im,''face_coeffs'',i)');
    end;
  end;
  gf_integ_get(im,'char');
  asserterr('gf_integ(''IM_TRIANGLE(0)'')');
  
  im=gf_integ('IM_EXACT_SIMPLEX(3)');
  dim = gf_integ_get(im,'dim');  
  gfassert('dim==3');
  ise = gf_integ_get(im,'is_exact');
  gfassert('ise');
  asserterr('gf_integ_get(im,''nbpts'')');
  gf_integ_get(im,'char');
  
  asserterr('gf_integ(''IM_EXACT_SIMPLEX(0)'')');
