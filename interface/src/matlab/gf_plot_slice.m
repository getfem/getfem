function [hfaces, htube, hquiver, hmesh]=gf_plot_slice(sl,varargin)
% function [hfaces, htube, hquiver, hmesh]=gf_plot_slice(sl,varargin)
% this function is used to plot a slice of mesh/mesh_fem (see gf_slice)
%
% The options are specified as pairs of "option name"/"option value"
%
%           OPTION NAME       DEFAULT VALUE         ACTION
%                    data    []                  data to be plotted (one value per slice node)
%             convex_data    []                  data to be plotted (one value per mesh convex)
%                    mesh    'auto'              'on' -> show the mesh (faces of edges), 
%                                                'off' -> ignore mesh
%              mesh_edges    'on'                show mesh edges ?
%        mesh_edges_color    [0.60 0.60 1]       color of mesh edges
%        mesh_edges_width    0.70                width of mesh edges
%        mesh_slice_edges    'on'                show edges of the slice ?
%  mesh_slice_edges_color    [0.70 0 0]
%  mesh_slice_edges_width    0.50
%              mesh_faces    'off'               'on' -> fill mesh faces (otherwise they are transparent)
%        mesh_faces_color    [0.75 0.75 0.75]
%                  pcolor    'on'                if the field is scalar, a color plot of its values is plotted
%                  quiver    'on'                if the field is vector, represent arrows
%          quiver_density    50                  density of arrows in quiver plot
%            quiver_scale    1                   density of arrows in quiver plot 
%                    tube    'on'                use tube plot for 'filar' (1D) parts of the slice
%              tube_color    'red'               color of tubes (ignored if 'data' is not empty and 'pcolor' is on)
%             tube_radius    '0.5%'              tube radius; you can use a constant, or a percentage 
%                                                (of the mesh size) or a vector of nodal values
%             showoptions    'on'                display the list of options
%  
% the 'data' and 'convex_data' are mutually exclusive.
%  
% RETURNS: handles to the various graphical objects created.  
%  Copyright (C) 1999-2020 Yves Renard
%
%  This file is a part of GetFEM
%
%  GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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

  %mf=struct(mf);
  hfaces=[];
  hquiver=[];
  hmesh=[];
  htube=[];
  mdim = gf_slice_get(sl, 'dim');
  
  if (gf_slice_get(sl, 'nbsplxs', 3)),
    warning('won''t plot 3D slices, extract the slice boundary first');
  end;
  
  if (mdim ~= 2 & mdim ~= 3),
    error('only 2D and 3D mesh are handled by this function');
  end;
  
  opt = struct('data',[],... % data to be plotted on the slice (on slice nodes)
	       'convex_data',[],...% data to be plotted (given on the mesh convexes)
	       'mesh','auto',...  % show the mesh ?
	       'mesh_edges','on',... % show mesh edges ?
               'mesh_edges_color',[.6 .6 1],...
               'mesh_edges_width',.7,...
               'mesh_slice_edges','on',...
	       'mesh_slice_edges_color',[.7 .0 0],...
               'mesh_slice_edges_width',.5,...
               'mesh_faces','off',... % fill mesh faces (otherwise they are transparent)
               'mesh_faces_color',[.75 .75 .75],...
	       'pcolor','on',... % if the field is scalar, a color plot of its values is plotted
	       'quiver','on',... % if the field is vector, represent arrows 	       
	       'quiver_density',50,... % density of arrows in quiver plot
	       'quiver_scale',1,... % scaling of arrows (0=>no scaling)
	       'tube','on',... % use tube plot for linear parts of the slice
	       'tube_color','red',... % color of tubes (ignored if 'data' is not empty)
               'tube_radius','0.5%',...  % tube radius; you can use a constant, or a percentage (of the mesh size) or a vector of nodal values
               'showoptions','off'); % list options used

  opt = getopt(opt,varargin);
  
  
  
  %qdim = numel(U) / gf_slice_get(sl, 'nbpts');
  %if (fix(qdim) ~= qdim), error('wrong number of elements for U'); end;

  if (ison(opt.showoptions)), disp(opt); end;
  
  if (~isempty(opt.convex_data)), 
    if (~isempty(opt.data)),
      error('"data" and "convex_data" are mutually exclusive');
    end;
    opt.data = gf_slice_get(sl, 'interpolate_convex_data', opt.convex_data);
  end;
  
  if (isauto(opt.mesh)),
    if (isempty(opt.data)), opt.mesh = 'on'; 
    else opt.mesh = 'off'; end; 
  end;
  
  Pm = gf_slice_get(sl,'pts'); 
  if (numel(Pm) == 0) return; end;
  if (~isempty(opt.data) & size(opt.data,2) ~= size(Pm,2)),
    error(sprintf('wrong dimensions for the data (has %d columns, should have %d columns)',...
		  size(opt.data,2),size(Pm,2)));
  end;

  P = cell(mdim,1);
  for i=1:mdim,
    P{i} = Pm(i,:);
    T{i} = gf_slice_get(sl,'splxs', i);
    box(i,:) = [min(P{i}) max(P{i})];
  end;
  
  if (ischar(opt.tube_radius) & numel(opt.tube_radius) & opt.tube_radius(end)=='%')
    opt.tube_radius = max(abs(box(:,2)-box(:,1))) * 0.01 * str2num(opt.tube_radius(1:end-1));
  end;
  
  % save graphical context
  cax = newplot;
  cfig = get(cax,'Parent');
  hold_state = ishold;
  ax_nextplot = lower(get(cax,'NextPlot'));
  fig_nextplot = lower(get(cfig,'NextPlot'));

 % handle simplexes of dimension 1

  if (~isempty(T{1})),
    [htube,hmesh]=do_plot_1D(P,T{1},opt);
  end;
  [hfaces,h,hquiver]=do_plot_2D(sl,P,T{2},opt); hmesh=[hmesh(:)' h(:)'];
  
  if (mdim == 3), view(3); else view(2); end;
  if (strcmpi(get(gcf,'renderer'),'opengl')),
    warning('OpenGL renderer does not work well with getfem, changing to zbuffer');
  end;
  set(gcf,'renderer','zbuffer');
 
  if (~hold_state),
    set(cax,'DataAspectRatio', [1 1 1]);
    set(cax,'XLimMode','auto');
    set(cax,'YLimMode','auto');
    set(cax,'ZLimMode','auto');
  end;

  % restore graphical context
  set(cax,'NextPlot',ax_nextplot);
  set(cfig,'NextPlot',fig_nextplot);


function [htube,hmesh]=do_plot_1D(P,T,opt)
  htube=[]; hmesh=[];
  if (isempty(T)) return; end;
  
  if (~ison(opt.tube)),
    for j=1:numel(P)
      C{j}=[P{j}(T(1,:));P{j}(T(2,:))];
    end;
    if (numel(P)==1) C{2}=zeros(size(C{1})); end;
    hmesh=line(C{:},'Color',opt.mesh_edges_color,'LineWidth',opt.mesh_edges_width);
  else
    if (~isempty(opt.data) & ison(opt.pcolor)),
      qdim = size(opt.data,1);
      if (qdim == 1),
	if (iscell(opt.data))
	  htube=plot_tube(P,T,opt.data{1},opt.tube_radius,opt.tube_color);
	else
	  htube=plot_tube(P,T,opt.data,opt.tube_radius,opt.tube_color);
	end;
      else
	warning('1D slices not supported for vector data..');
      end;
    else 
      htube=plot_tube(P,T,[],opt.tube_radius,opt.tube_color);
    end;
  end;

% cell2mat not available in matlab R12
function M=mycell2mat(C)
  M=cat(1,C{:});

% plots a "tube" along edges, with color given by D, and a possibly varying radius
% radius: constant or equal to nb points
% D(ata): empty or equal to nb points or nb segments
function h=plot_tube(P, T, D, radius, tubecolor)
  h = [];
  P = mycell2mat(P);
  if (isempty(T)) return; end;
  T=double(T); % matlab 6.5 is not able to handle operator '+' on
               % int32 ...
  it0 = T(1,1); nT = size(T,2); nP = size(P,2); mdim=size(P,1);
  if (mdim == 2) P = [P; zeros(1,nP)]; mdim = 3; end; % handle 2D slices
  % convert radius to node data
  if (numel(radius)==1), 
    radius = radius*ones(1,nP); %radius(1)=0.5; radius(end)=0.5;
  elseif (numel(radius)==nT),
    radius = ([radius(1) radius(:)']+[radius(:)' radius(end)])/2;
  end;
  if (size(D,1) > 1), error('only scalar data can be represented on a tube_plot'); end;
  if (size(D,2)==nP), 
    point_data=1; 
  else 
    point_data=0; 
  end;
  nsubdiv = 20; ct = cos((0:nsubdiv)*2*pi/nsubdiv); ct(end)=ct(1); st = sin((0:nsubdiv)*2*pi/nsubdiv); st(end)=st(1);
  cnt=0;
  while (1),
    % search for consecutive edge points
    it1 = it0;
    while (it1 < nT & T(1,it1+1) == T(2,it1)) it1 = it1+1; end;
    %disp(sprintf('sequence: %d - %d -- [%d-%d] - [%d-%d]',it0,it1,T(1,it0),T(2,it0),T(1,it1),T(2,it1)))
    % extract the sequence of points
    ip = [T(1,it0) T(2,it0:it1)];
    p = P(:,ip);     
    if (numel(D)),
      if (point_data) d = D(ip); else d = D(it0:it1); end;
    end;
    nseg = it1-it0+1;
    % compute the normals of edges
    normals = zeros(3, 2, nseg);
    tang = p(:,2:end) - p(:,1:end-1); tang = tang ./ repmat(sqrt(sum(tang.^2,1)),size(tang,1),1);
    for i=1:nseg
      normals(:,:,i) = null(tang(:,i)'); % won't be ok if normals have an
                                         % important rotation from a segment
                                         % to another
                                         % VERY PROBABLE BUG!!!      
    end;    
    X=zeros(mdim,nsubdiv+1,numel(ip));
    for i=1:numel(ip),
      if (i == 1),
        n=normals(:,:,i); 
      elseif (i == numel(ip)) 
        n=normals(:,:,end);
      else
        n = (normals(:,:,i-1)+normals(:,:,i))/2;
      end
      for k=1:nsubdiv+1,
        X(:,k,i) = p(:,i) + radius(ip(i))*(n(:,1)*ct(k) + n(:,2)*st(k));
      end;
    end;
    if (numel(D)),
      C=repmat(d,nsubdiv+1,1);
      h = [h; surface(squeeze(X(1,:,:)), squeeze(X(2,:,:)), squeeze(X(3,:,:)),C,...
		      'linestyle','none','FaceColor','interp')];
    else
      h = [h; surface(squeeze(X(1,:,:)), squeeze(X(2,:,:)), squeeze(X(3,:,:)),...
		      'linestyle','none','facecolor',tubecolor)];
    end;
    cnt=cnt+1;
    it0 = it1+1;
    if (it0 > nT) return; end;
  end;

% draw faces
function [hfaces,hmesh,hquiver]=do_plot_2D(sl,P,T,opt)
  hfaces=[]; hmesh=[]; hquiver=[];
  mdim=numel(P);
  if (numel(T)),
    d={};  
    if (ison(opt.pcolor) & size(opt.data,1)==1 & ~isempty(opt.data)),
      d={'FaceVertexCData',opt.data(:),'FaceColor','interp'};
    elseif (isempty(opt.data) & ison(opt.mesh_faces))
      d={'FaceVertexCData',opt.mesh_faces_color, 'FaceColor','flat'};
    end;    
    if (~isempty(d)),
      hfaces=patch('Vertices',mycell2mat(P)','Faces',T',d{:}, ...
		   'EdgeColor','none');
    end;
    if (ison(opt.quiver)),
      if (size(opt.data,1)>1),
	hquiver=do_quiver_plot(P,opt.data,opt);
      end;
    end;  
  end;
  if (ison(opt.mesh) & (ison(opt.mesh_edges) | ison(opt.mesh_slice_edges))),
    [p,t1,t2] = gf_slice_get(sl,'edges');
    if (ison(opt.mesh_edges))
      hmesh=patch('Vertices',p','Faces',t1','EdgeColor',opt.mesh_edges_color,'LineWidth',opt.mesh_edges_width);
    end;
    if (ison(opt.mesh_slice_edges)),
      hmesh=[hmesh patch('Vertices',p','Faces',t2','EdgeColor',opt.mesh_slice_edges_color,'LineWidth',opt.mesh_slice_edges_width)];
    end;
  end;  
  
  
% arrow plot
function hquiver=do_quiver_plot(P,U,opt)
  hquiver=[];
  P=mycell2mat(P); mdim=size(P,1);qdim=size(U,1);
  nP = size(P,2);
  ptlst=1:nP;
  bmin = min(P); bmax = max(P);
  xyscale = max(bmax-bmin);
  qradius2 = (xyscale/opt.quiver_density)^2;
  vscale  = max(max(abs(U)));
  qlst = [];
  rm=[];
  while (numel(ptlst)>0)
    ii = ptlst(1); 
    qlst = [qlst ii];
    x = P(1,ii); y = P(2, ii); 
    if (mdim == 2),
      rm = (find((P(1,:)-x).^2 + (P(2,:)-y).^2 < qradius2));
    elseif (mdim == 3),
      z = P(3,ii);
      rm = (find((P(1,:)-x).^2 + (P(2,:)-y).^2 + (P(3,:)-z).^2 < qradius2));
    end;
    if (numel(rm)==0) error('internal error in gf_plot'); end;
    ptlst = setdiff(ptlst, rm);
  end;
  if (qdim == 2),
    hquiver=quiver(P(1,qlst),P(2,qlst),U(1,qlst),U(2,qlst),opt.quiver_scale);
    set(hquiver,'Color', [0 .4 0]);
  else
    hquiver=quiver3(P(1,qlst),P(2,qlst),P(3,qlst),U(1,qlst),U(2,qlst),U(3,qlst),...
                    opt.quiver_scale);
    set(hquiver,'Color', [0 .4 0]);
  end;
  
function r=ison(v)
  r = strcmpi(v,'on');
  
function r=isauto(v)
  r = strcmpi(v,'auto');
