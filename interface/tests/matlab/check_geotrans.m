% Copyright (C) 2005-2020 Julien Pommier.
%
% This file is a part of GetFEM
%
% GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
% under  the  terms  of the  GNU  Lesser General Public License as published
% by  the  Free Software Foundation;  either version 3 of the License,  or
% (at your option) any later version along with the GCC Runtime Library
% Exception either version 3.1 or (at your option) any later version.
% This program  is  distributed  in  the  hope  that it will be useful,  but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
% License and GCC Runtime Library Exception for more details.
% You  should  have received a copy of the GNU Lesser General Public License
% along  with  this program;  if not, write to the Free Software Foundation,
% Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.


function check_geotrans(iverbose,idebug)
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
  gt=gf_geotrans('GT_PK(3,1)');
  dim = gf_geotrans_get(gt,'dim');  
  gfassert('dim==3');
  islin = gf_geotrans_get(gt,'is_linear');
  gfassert('islin==1');
  npt=gf_geotrans_get(gt,'nbpts');
  gfassert('npt==4');  
  p = gf_geotrans_get(gt,'pts');
  gfassert('size(p)==[3 4]');
  n = gf_geotrans_get(gt, 'normals');
  gfassert('norm(n(:,1)-0.5774)<0.001');
  s = gf_geotrans_get(gt, 'char');
  gfassert('s==''GT_PK(3,1)''');
  
  
  s='GT_PRODUCT(GT_PRODUCT(GT_PK(2,3),GT_PK(1,1)),GT_QK(2,3))';
  gt=gf_geotrans(s);
  islin = gf_geotrans_get(gt,'is_linear');
  gfassert('islin==0');
  gf_geotrans_get(gt, 'char');
  
  s='GT_LINEAR_PRODUCT(GT_PK(1,1),GT_PK(1,1))';
  gt=gf_geotrans(s);
  islin = gf_geotrans_get(gt,'is_linear');
  gfassert('islin==1');
