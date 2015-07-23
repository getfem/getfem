% Copyright (C) 2012-2015 Yves Renard, Julien Pommier.
%
% This file is a part of GetFEM++
%
% GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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


% addpath ~/source++/getfem++/contrib/xfem_stab_unilat_contact/

gf_workspace('clear all');
mf = gf_mesh_fem('load', 'xfem_stab_unilat_contact_ls.mf');
%lsU = -load('xfem_stab_unilat_contact_ls.U')';




clf

disp('plotting the lagrange multipliers on the contact interface');
  sll=gf_Slice('load','xfem_stab_unilat_contact.sll');
  slL=load('xfem_stab_unilat_contact.slL')';
  P0=gf_slice_get(sll, 'pts');
  [h1,h2,h3,h4]=gf_plot_slice(sll, 'tube','off','mesh_slice_edges_color',[.3 .3 .3]);
  hold on;
  gf_slice_set(sll,'pts',[P0 ; -max(slL,-100)]);
  [hh1,hh2,hh3,hh4]=gf_plot_slice(sll, 'tube','off','mesh_slice_edges_color','black','mesh_slice_edges_width',1.5);
  sl=gf_Slice('load','xfem_stab_unilat_contact.sl');
  gf_plot_slice(sl,'mesh','on');
  
  npt = size(P0, 2);
  P0 = [P0;zeros(1,npt)];
  P1 = gf_slice_get(sll,'pts'); 
  lseg = gf_slice_get(sll,'splxs', 1);
  F=[lseg(1,:) lseg(2,:); lseg(2,:) npt+lseg(2,:); npt+lseg(1,:) npt+lseg(1,:)];
  %F=[lseg; npt+lseg(2,:)];
  h=patch('Vertices',[P0 P1]', 'Faces', F');
  hold off;
  set(h,'FaceAlpha',0.3);
  set(h,'LineStyle','none');
  set(gcf,'renderer','opengl');
  set(gcf,'color','white');
  axis on;
  view(3);
  camzoom(1.2);
  axis([0    1   -0.500    0.5000 -1 1]);
  %print(gcf,'-dpng','-r300', 'lagrange_multipliers.png');
