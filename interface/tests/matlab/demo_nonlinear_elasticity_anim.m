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


% replay all the computations of demo_nonlinear_elasticity.m

load demo_nonlinear_elasticity_U.mat
nbstep = size(UU,1);
m=gfMesh('from string', m_char);
mfu=gfMeshFem('from string',mfu_char,m);
mfdu=gfMeshFem('from string',mfdu_char,m);



sl=gfSlice({'boundary'}, m, 16, gf_mesh_get(m,'outer faces'));
P0=gf_slice_get(sl,'pts');
for step=1:1:nbstep,
  U=UU(step,:);
  VM=VVM(step,:);
  
  slU=gf_compute(mfu,U,'interpolate on',sl);
  slVM=gf_compute(mfdu,VM,'interpolate on',sl);
  
  gf_slice_set(sl,'pts', P0+slU);
  
  gf_plot_slice(sl, 'data', slVM, 'mesh_edges','on', 'mesh','on'); 
  
  %gf_plot(mfdu,VM,'mesh','on', 'cvlst',gf_mesh_get(mfdu,'outer faces'), 'deformation',U,'deformation_mf',mfu,'deformation_scale', 1, 'refine', 16); 
  axis([-3     6     0    20    -2     2]);
  caxis([0 .15]);
  view(30+20*w, 23+30*w);  
  campos([50 -30 80]);
  camva(8);
  camup;
  camlight; 
  axis off;
  print(gcf, '-dpng', '-r150', sprintf('torsion%03d',step));
end;
