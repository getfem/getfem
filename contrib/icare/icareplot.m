% Copyright (C) 2012-2012 Yves Renard, Michel Fournie.
%
% This file is a part of GETFEM++
%
% Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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


function icareplot(nn)
  global mfu U mfdu DU rot
  gf_workspace('clear all');
  mfu=gfMeshFem('load','icare.mf_u');
  mfdu=gfMeshFem(get(mfu,'linked mesh')); 
  set(mfdu,'fem',gfFem('FEM_PK_DISCONTINUOUS(2,1)'));

  for n=nn,
    U=load(sprintf('icare.U%d',n))';
    DU=gf_compute(mfu,U,'gradient',mfdu);
    rot=DU(1,2,:)-DU(2,1,:); 
  
    subplot(2,1,1); gf_plot(mfu,U,'refine',2,'norm','on'); 
    %caxis([0 1.5]); 
    colorbar;
    subplot(2,1,2); gf_plot(mfdu,rot(:)','refine',1); 
    %caxis([-2 2]); 
    colorbar;
    disp('press any key'); pause;
  end;
