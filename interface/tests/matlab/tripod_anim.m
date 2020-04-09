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


% You should run demo_tripod first ...
m=gf_mesh('import','gid','../meshes/tripod.GiD.msh');
% mfu=gfMeshFem('from string', smfu, m);
% mfdu=gfMeshFem('from string', smfdu, m);

%gf_plot_mesh(m,'cvlst',gf_mesh_get(m,'outer
%faces'),'curved','on','edges_color',[1 0 0]); hold on;
pr=1; haut=0;
for r=[6 14 26],
  disp('slicing...'); tic;
  %sl=gfSlice({'cylinder', 0, [0;0;0], [0;1;0], r}, m, 4);
  %sl=gf_slice(m,{'boundary',{'cylinder', [0;0;0], [0;1;0], 15}},4);

  %sl=gf_slice(m,{'union',{'cylinderb', [0;0;0], [0;1;0], 15},{'cylinderb', [0;0;0], [0;1;0], 7}},4);
  sl=gfSlice({'boundary',{'diff',{'cylinder', -1, [0;0;0], [0;1;0], r},{'cylinder', -1, [0;0;0], [0;1;0], r-4}}},m, 6);
  %sl=gf_slice(m,{'none'},4, gf_mesh_get(m,'outer faces'));
  %sl=gf_slice(m,{'boundary', {'none'}},4);
  %sl=gf_slice(m,{'ballb', [0;0;0], 10},2);
  %sl=gf_slice(m,{'planarb',[0;0;0],[0;0;1]},1);
  disp(sprintf('..........done in %3.2f sec',toc));

  
  P=gf_slice_get(sl,'pts'); 
  P(2,:)=P(2,:)-haut;
  sl
  D=gf_compute(mfdu,VM,'interpolate on',sl);
  gf_plot_slice(sl, 'mesh','on','data',D,'pcolor','on','mesh_edges_color',[1 1 .7]);
  hold on;
  pr=r; haut=haut+24;
end;
gf_workspace('stats');
hold off; caxis([0 5]); view(180,-50); 
gf_colormap('earth'); caxis([0 100]); camlight;
