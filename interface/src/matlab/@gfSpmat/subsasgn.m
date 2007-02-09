function [A] = subsasgn(A,S,B)
  nout = max(nargout,1);  
  cnt=1;
  switch S(1).type    
   case '{}'
    error('Cell array indexing not supported by gfSpmat objects')
   case {'()'}
    ri=S(1).subs{1}; rj = S(1).subs{2}; sz = gf_spmat_get(A, 'size');
    if (ri == ':') ri = 1:double(sz(1)); end;
    if (rj == ':') rj = 1:double(sz(2)); end;
    gf_spmat_set(A, 'set', ri, rj, B);
   case '.'
    switch index(1).subs
     case 'diag'
      gf_spmat_set(A, 'diag', B);
     otherwise
      gf_spmat_get(A,S(1).subs,B);
    end
  end


