function [hmesh,hbound,hfill,hvert,hconv,hdof]=gf_plot_mesh(M, varargin)
% function [hmesh,hbound,hfill,hvert,hconv,hdof]=gf_plot_mesh(M, [,properties])
%                     [,'cvlst',CVLST] ['boundaries'[BLST]])
%   General mesh plotting function.
%  
%   H=gf_plot_mesh(M) displays a mesh.
%
%   properties are:
%    'vertices', {'off'|'on'}     displays also vertices numbers. 
%    'convexes', {'off'|'on'}     displays also convexes numbers. 
%    'dof',{'off'|'on'}           displays also finite element nodes.
%    'regions',BLST               displays the boundaries listed in BLST.
%    'cvlst',CVLST                display only the listed convexes. If
%   CVLST has two rows, display only the faces listed in the second row.
%    'edges', {'on' | 'off'}      display edges ?
%    'faces', {'off'|'on'}        fills each 2D-face of the mesh
%    'curved', {'off'|'on'}       displays curved edges
%    'refine',N                   refine curved edges and filled faces N times  
%    'deformation', Udef          optionnal deformation applied to the mesh (M must be a mesh_fem object)
%    'edges_color',[.6 .6 1]      RGB values for the color of edges
%    'edges_width',1              
%    'faces_color',[.75 .75 .75]) RGB values for the color of faces
%    'quality',{ 'off' | 'on' }   Display the quality of the mesh.
%
%   CAUTION:
%     For 'dof', M should be a mesh_fem identifier, 
%   not a simple mesh object.
%  
%   $Id$
%  Copyright (C) 1999-2017 A. Huard, Y. Renard, J. Pommier
%
%  This file is a part of GetFEM++
%
%  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
%  under  the  terms  of the  GNU  Lesser General Public License as published
%  by  the  Free Software Foundation;  either version 3 of the License,  or
%  (at your option) any later version along with the GCC Runtime Library
%  Exception either version 3.1 or (at your option) any later version.
%  This program  is  distributed  in  the  hope  that it will be useful,  but
%  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
%  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
%  License and GCC Runtime Library Exception for more details.
%  You  should  have received a copy of the GNU Lesser General Public License
%  along  with  this program;  if not, write to the Free Software Foundation,
%  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

if nargin<1,
  error('Too few input arguments')
end
hmesh=[];
hbound=[];
hfill=[];
hvert=[];
hconv=[];
hdof=[];

mdim = gf_mesh_get(M,'dim');
if (mdim <= 2),
  defaultref=8;
else 
  defaultref=4;
end;
opt = struct('vertices','off',...
	     'convexes','off',...
	     'dof','off',...
	     'regions',[],...
	     'boundaries',[],...
	     'cvlst',[],...
	     'edges','on',...
	     'faces','off',...
	     'explode',0,...
	     'quality','off',...
	     'curved','off',...
	     'refine',defaultref,...
	     'deformation',[],...
	     'edges_color',[.6 .6 1],...
	     'edges_width',.7,...
	     'faces_color',[.75 .75 .75]);

% parse argument list
opt = getopt(opt,varargin);

if (numel(opt.boundaries) == 0),
  opt.boundaries = opt.regions;
end;

if (strcmpi(opt.boundaries, 'all')),
  opt.boundaries = gf_mesh_get(M, 'boundaries');
end;

% init cvlst and cvflst
if (isempty(opt.cvlst)) 
  cvlst = gf_mesh_get(M,'cvid'); 
  cvflst = [cvlst; zeros(1,numel(cvlst))];
else 
  cvlst = opt.cvlst;
  cvflst = cvlst;  
  if (size(cvflst,1)==2),
    cvlst = unique(cvlst(1,:));
  end;
end;

PXY = gf_mesh_get(M, 'pts');

E = gf_mesh_get(M, 'edges', cvflst);
if (~ison(opt.curved) && isempty(opt.deformation)),
  X = [PXY(1,E(1,:)); PXY(1,E(2,:))];
  if (mdim == 1),
    Y = zeros(size(X));
    PXY = [PXY; zeros(size(PXY))];
  elseif (mdim >= 2),
    Y = [PXY(2,E(1,:)); PXY(2,E(2,:))];
    if (mdim == 3),
      Z = [PXY(3,E(1,:)); PXY(3,E(2,:))];
    end;
  end;
else
  % here, mdim is always >= 2
  if (~isempty(opt.deformation)),
    vE = gf_compute(M,opt.deformation,'mesh edges deformation',opt.refine,cvflst);
  else
    vE = gf_mesh_get(M, 'curved edges', opt.refine, cvflst);
  end;
  ni = size(vE,2);
  ne = size(vE,3);
  X = reshape(vE(1,:,:), [ni ne]);
  Y = reshape(vE(2,:,:), [ni ne]);
  if (mdim == 3),
    Z = reshape(vE(3,:,:), [ni ne]);
  end;
end;

% get the viewable pts id
PID = union(E(1,:),E(2,:));
if (mdim > 3), error('sorry, only mesh of dimension <= 3 allowed'); end;
nbpts = size(PXY,2);

Bmax  = max(PXY(:,PID)')';
Bmin  = min(PXY(:,PID)')';
Bdiff = Bmax - Bmin;
Bdiff = Bdiff + (Bdiff == 0); % remplace 0.0 par 1.0
ecart=Bdiff/150;



if (ison(opt.convexes)),
  cv_center = zeros(max(mdim,2),numel(cvlst));
  % find convexes centers
  [cv_pid, cv_idx] = gf_mesh_get(M, 'pid from cvid',cvlst);
  for i=1:length(cvlst),
    cv_center(:,i) = mean(PXY(:, cv_pid(double(cv_idx(i)):double(cv_idx(i+1))-1)),2);
  end;
end;

if (ison(opt.dof)),
  Q = gf_mesh_fem_get(M, 'qdim');
  dofid = gf_mesh_fem_get(M, 'dof from cv', cvlst);
  [dofpos] = gf_mesh_fem_get(M, 'dof nodes', dofid);
  [keep] = find([1 any(dofpos(:,2:end) ~= dofpos(:,1:end-1),1)]);
  dofmult = [keep(2:end)-keep(1:end-1) size(dofpos,2)+1-keep(end)];
  dofpos = dofpos(:, keep); dofid = dofid(keep);
  if (mdim == 1) dofpos = [dofpos; zeros(size(dofpos))]; end;
end;

for bnum=1:length(opt.boundaries),
  cvf = gf_mesh_get(M, 'boundary', opt.boundaries(bnum));
  
  bid = gf_mesh_get(M, 'edges', cvf, 'merge convex');
  if (bnum == 8) disp(bid); end;
  
  bedge{bnum} = zeros(2, size(bid,2), mdim);
  for i=1:max(mdim,2),
    bedge{bnum}(:,:,i) = [PXY(i,bid(1,:)); PXY(i,bid(2,:))];
  end;
end;

% save graphical context
cax = newplot;
cfig = get(cax,'Parent');
hold_state = ishold;
ax_nextplot = lower(get(cax,'NextPlot'));
fig_nextplot = lower(get(cfig,'NextPlot'));

disp('plotting mesh...');
if (mdim <= 2),
  if (ison(opt.edges)) 
    hmesh = line(X, Y, 'Color',opt.edges_color,'LineWidth',opt.edges_width); 
  end;
  for bnum=1:length(opt.boundaries),
    hbound{bnum} = line(bedge{bnum}(:,:,1), bedge{bnum}(:,:,2), 'Color','red','LineWidth',2);
  end
  if (ison(opt.vertices)),
    hvert = text(PXY(1,PID)+ecart(1), PXY(2,PID)+ecart(2), num2str(double(PID')),...
		 'HorizontalAlignment','center','VerticalAlignment','middle');
  end;
  if (ison(opt.convexes)),
    hconv = text(cv_center(1,:), cv_center(2,:), num2str(double(cvlst')), ...
		 'HorizontalAlignment','center','VerticalAlignment','middle',...
		 'Color', [.7 0 0]);
  end;
  if (ison(opt.dof)),
    hdof = zeros(numel(dofid),1);
    for i=1:numel(dofid),
      if (dofmult(i)==1) s=int2str(dofid(i)); else s=sprintf('%d*%d',dofid(i),dofmult(i)); end;
      hdof(i) = text(dofpos(1,i)-ecart(1), dofpos(2,i)-ecart(2), s,...
		      'HorizontalAlignment','center','VerticalAlignment','middle',...
		      'Color', [0 .4 0]);
    end;
  end;
else
  if (ison(opt.edges)) 
    hmesh = line(X, Y, Z, 'Color',opt.edges_color,'LineWidth',opt.edges_width); 
  end;
  for bnum=1:length(opt.boundaries),
    hbound{bnum} = line(bedge{bnum}(:,:,1), bedge{bnum}(:,:,2), bedge{bnum}(:,:,3), 'Color','red','LineWidth',2);
  end
  if (ison(opt.vertices)),
    hvert = text(PXY(1,PID)+ecart(1), PXY(2,PID)+ecart(2), PXY(3,PID)+ecart(3), num2str(PID'),...
		 'HorizontalAlignment','center','VerticalAlignment','middle','Color', [.0 0 0]);
  end;
  if (ison(opt.convexes)),
    hconv = text(cv_center(1,:), cv_center(2,:), cv_center(3,:), num2str(cvlst'), ...
		 'HorizontalAlignment','center','VerticalAlignment','middle',...
		 'Color', [.7 0 0]);
  end;
  if (ison(opt.dof)),
    hdof = text(dofpos(1,:)-ecart(1), dofpos(2,:)-ecart(2), dofpos(3,:)-ecart(3), num2str(dofid'),...
		 'HorizontalAlignment','center','VerticalAlignment','middle',...
		 'Color', [0 .4 0]);
  end;
end

if (ison(opt.quality)),
  q=gf_mesh_get(M,'quality', cvflst(1,:));
  qmf=gf_mesh_fem(M);
  gf_mesh_fem_set(qmf, 'classical fem', 0);
  [a,b] = gf_mesh_fem_get(qmf, 'dof from cvid', cvflst(1,:));
  Q=zeros(1, gf_mesh_fem_get(qmf, 'nbdof'));
  for k=1:numel(b)-1,
    Q(a(b(k))) = q(k);
  end
end;

if (opt.explode ~= 0),
  sl=gf_slice({'explode',opt.explode},M,opt.refine,cvflst);
  data={};
  if (ison(opt.quality)),
    sQ = gf_compute(qmf,Q,'interpolate on',sl);
    data={'data',sQ};
  end;
  gf_plot_slice(sl, data{:},'mesh_faces',opt.faces,...
		'mesh_edges_color', opt.edges_color, ...
		'mesh_edges_width',opt.edges_width, ...
		'mesh_faces_color',opt.faces_color);
  gf_delete(sl); light;
elseif (ison(opt.quality)),
  hold on; gf_plot(qmf, Q, 'cvlst', cvflst); hold off;
elseif (ison(opt.faces)),
  
  % should be replaced by a gf_plot_slice ..
  T = gf_mesh_get(M, 'triangulated surface', opt.refine, cvflst);
  if (mdim == 2),
    hfill=patch(T(1:mdim:(mdim*3),:),T(2:mdim:(mdim*3),:), ...
		opt.faces_color, 'Erasemode','normal','Edgecolor','none');
  else
    hfill=patch(T(1:mdim:(mdim*3),:),T(2:mdim:(mdim*3),:), T(3:mdim:(mdim*3),:), ...
		opt.faces_color, 'Erasemode','normal','Edgecolor','none');
  end;
  light;
end;

if (ison(opt.quality))
  gf_delete(qmf);
end;

if (~hold_state),
  set(cax,'DataAspectRatio', [1 1 1]);
  set(cax,'XLimMode','auto');
  set(cax,'YLimMode','auto');
  set(cax,'ZLimMode','auto');
end;
% restore graphical context
set(cax,'NextPlot',ax_nextplot)
set(cfig,'NextPlot',fig_nextplot)

 
function r=ison(v)
  r = strcmpi(v,'on');
