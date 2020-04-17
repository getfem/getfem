// Copyright (C) 2009-2020 Yann Colette
// 
//  This file is a part of GetFEM++
// 
//  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
//  under  the  terms  of the  GNU  Lesser General Public License as published
//  by  the  Free Software Foundation;  either version 3 of the License,  or
//  (at your option) any later version along with the GCC Runtime Library
//  Exception either version 3.1 or (at your option) any later version.
//  This program  is  distributed  in  the  hope  that it will be useful,  but
//  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
//  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
//  License and GCC Runtime Library Exception for more details.
//  You  should  have received a copy of the GNU Lesser General Public License
//  along  with  this program;  if not, write to the Free Software Foundation,
//  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

function [hfaces, htube, hquiver, hmesh]=gf_plot_slice(sl,varargin)
// function [hfaces, htube, hquiver, hmesh]=gf_plot_slice(sl,varargin)
// this function is used to plot a slice of mesh/mesh_fem (see gf_slice)
//
// The options are specified as pairs of 'option name'/'option value'
//
//           OPTION NAME       DEFAULT VALUE         ACTION
//                    data    []                  data to be plotted (one value per slice node)
//             convex_data    []                  data to be plotted (one value per mesh convex)
//                    mesh    'auto'              'on' -> show the mesh (faces of edges), 
//                                                'off' -> ignore mesh
//              mesh_edges    'on'                show mesh edges ?
//        mesh_edges_color    [0.60 0.60 1]       color of mesh edges
//        mesh_edges_width    0.70                width of mesh edges
//        mesh_slice_edges    'on'                show edges of the slice ?
//  mesh_slice_edges_color    [0.70 0 0]
//  mesh_slice_edges_width    0.50
//              mesh_faces    'off'               'on' -> fill mesh faces (otherwise they are transparent)
//        mesh_faces_color    [0.75 0.75 0.75]
//                  pcolor    'on'                if the field is scalar, a color plot of its values is plotted
//                  quiver    'on'                if the field is vector, represent arrows
//          quiver_density    50                  density of arrows in quiver plot
//            quiver_scale    1                   density of arrows in quiver plot 
//                    tube    'on'                use tube plot for 'filar' (1D) parts of the slice
//              tube_color    'red'               color of tubes (ignored if 'data' is not empty and 'pcolor' is on)
//             tube_radius    0.05                tube radius; you can use a constant or a vector of nodal values
//             showoptions    'on'                display the list of options
//  
// the 'data' and 'convex_data' are mutually exclusive.
//  
// RETURNS: handles to the various graphical objects created.  
////////////////////////

[nargout,nargin] = argn();

if nargin<1 then
  error('Too few input arguments')
end

opts = build_options_list(varargin(:));

hfaces  = [];
hquiver = [];
hmesh   = [];
htube   = [];
mdim = gf_slice_get(sl, 'dim');

if (gf_slice_get(sl, 'nbsplxs', 3)) then
  warning('won''t plot 3D slices, extract the slice boundary first');
end

if (mdim ~= 2 & mdim ~= 3) then
  error('only 2D and 3D mesh are handled by this function');
end

[o_data,err]                  = get_param(opts,'data',[]); // data to be plotted on the slice (on slice nodes)
[o_convex_data,err]           = get_param(opts,'convex_data',[]); // data to be plotted (given on the mesh convexes)
[o_msh,err]                   = get_param(opts,'mesh','auto'); // show the mesh ?
[o_msh_edges,err]             = get_param(opts,'mesh_edges','on'); // show mesh edges ?
[o_msh_edges_color,err]       = get_param(opts,'mesh_edges_color',[.6 .6 1]);
[o_msh_edges_width,err]       = get_param(opts,'mesh_edges_width',.7);
[o_msh_slice_edges,err]       = get_param(opts,'mesh_slice_edges','on');
[o_msh_slice_edges_color,err] = get_param(opts,'mesh_slice_edges_color',[.7 0 0]);
[o_msh_slice_edges_width,err] = get_param(opts,'mesh_slice_edges_width',.5);
[o_msh_faces,err]             = get_param(opts,'mesh_faces','off'); // fill mesh faces (otherwise they are transparent)
[o_msh_faces_color,err]       = get_param(opts,'mesh_faces_color',[.75 .75 .75]);
[o_pcolor,err]                = get_param(opts,'pcolor','on'); // if the field is scalar, a color plot of its values is plotted
[o_quiver,err]                = get_param(opts,'quiver','on'); // if the field is vector, represent arrows 	 
[o_quiver_density,err]        = get_param(opts,'quiver_density',50); // density of arrows in quiver plot
[o_quiver_scale,err]          = get_param(opts,'quiver_scale',1); // scaling of arrows (0=>no scaling)
[o_tube,err]                  = get_param(opts,'tube','on'); // use tube plot for linear parts of the slice
[o_tube_color,err]            = get_param(opts,'tube_color','red'); // color of tubes (ignored if 'data' is not empty)
[o_tube_radius,err]           = get_param(opts,'tube_radius',0.05); // tube radius; you can use a constant, or a percentage (of the mesh size) or a vector of nodal values
[o_showoptions,err]           = get_param(opts,'showoptions','off'); // list options used

if (ison(o_showoptions)) then disp(opts); end;

if (~isempty(o_convex_data)) then 
  if (~isempty(o_data)) then
    error('''data'' and ''convex_data'' are mutually exclusive');
  end
  o_data = gf_slice_get(sl, 'interpolate_convex_data', o_convex_data);
end

if (isauto(o_msh)) then
  if (isempty(o_data)) then o_msh = 'on'; 
  else o_msh = 'off'; end; 
end

Pm = gf_slice_get(sl,'pts'); 
if (length(Pm) == 0) then return; end;
if (~isempty(o_data) & size(o_data,2) ~= size(Pm,2)) then
  error(sprintf('wrong dimensions for the data (has %d columns, should have %d columns)', size(o_data,2),size(Pm,2)));
end

P = list();
T = list();
for i=1:mdim
  P(i) = Pm(i,:);
  T(i) = gf_slice_get(sl,'splxs', i);
  box(i,:) = [min(P(i),'r') max(P(i),'r')];
end

// handle simplexes of dimension 1

if (~isempty(T(1))) then
  [htube,hmesh]=do_plot_1D(P,T(1),opts);
end

[hfaces,h,hquiver] = do_plot_2D(sl,P,T(2),opts); hmesh=[hmesh(:)' h(:)'];

h_current = gca();

if (mdim == 3) then 
  h_current.view = '3d'; 
else 
  h_current.view = '2d'; 
end
endfunction

////////////////
// do_plot_1D //
////////////////

function [htube,hmesh]=do_plot_1D(P,T,opt)

htube=[]; hmesh=[];
if (isempty(T)) then
  return; 
end

[o_data,err]                  = get_param(opt,'data',[]); // data to be plotted on the slice (on slice nodes)
[o_convex_data,err]           = get_param(opt,'convex_data',[]); // data to be plotted (given on the mesh convexes)
[o_msh,err]                   = get_param(opt,'mesh','auto'); // show the mesh ?
[o_msh_edges,err]             = get_param(opt,'mesh_edges','on'); // show mesh edges ?
[o_msh_edges_color,err]       = get_param(opt,'mesh_edges_color',[.6 .6 1]);
[o_msh_edges_width,err]       = get_param(opt,'mesh_edges_width',.7);
[o_msh_slice_edges,err]       = get_param(opt,'mesh_slice_edges','on');
[o_msh_slice_edges_color,err] = get_param(opt,'mesh_slice_edges_color',[.7 0 0]);
[o_msh_slice_edges_width,err] = get_param(opt,'mesh_slice_edges_width',.5);
[o_msh_faces,err]             = get_param(opt,'mesh_faces','off'); // fill mesh faces (otherwise they are transparent)
[o_msh_faces_color,err]       = get_param(opt,'mesh_faces_color',[.75 .75 .75]);
[o_pcolor,err]                = get_param(opt,'pcolor','on'); // if the field is scalar, a color plot of its values is plotted
[o_quiver,err]                = get_param(opt,'quiver','on'); // if the field is vector, represent arrows 	 
[o_quiver_density,err]        = get_param(opt,'quiver_density',50); // density of arrows in quiver plot
[o_quiver_scale,err]          = get_param(opt,'quiver_scale',1); // scaling of arrows (0=>no scaling)
[o_tube,err]                  = get_param(opt,'tube','on'); // use tube plot for linear parts of the slice
[o_tube_color,err]            = get_param(opt,'tube_color','red'); // color of tubes (ignored if 'data' is not empty)
[o_tube_radius,err]           = get_param(opt,'tube_radius',0.05); // tube radius; you can use a constant, or a percentage (of the mesh size) or a vector of nodal values
[o_showoptions,err]           = get_param(opt,'showoptions','off'); // list options used

if (~ison(o_tube)) then
  C = list();
  for j=1:length(P)
    C(j)=[P(j)(T(1,:));P(j)(T(2,:))];
  end
  if (length(P)==1) C(2)=zeros(size(C(1))); end;
  // hmesh = line(C(:),'Color',o_msh_edges_color);
  if length(C)==2 then
    plot2d(C(:));
    hmesh = gce();
    hmesh.children.thickness  = o_msh_edges_width;
    hmesh.children.foreground = o_msh_edges_color;
  else
    plot3d(C(:));
    hmesh = gce();
    hmesh.thickness  = o_msh_edges_width;
    hmesh.foreground = o_msh_edges_color;
  end
else
  if (~isempty(o_data) & ison(o_pcolor)) then
    qdim = size(o_data,1);
    if (qdim == 1) then
      if (typeof(o_data)=='list')
        plot_tube(P,T,o_data(1),o_tube_radius,o_tube_color);
      else
        plot_tube(P,T,o_data,o_tube_radius,o_tube_color);
      end
    else
      warning('1D slices not supported for vector data..');
    end
  else 
    plot_tube(P,T,[],o_tube_radius,o_tube_color);
  end
end
endfunction

////////////////
// mycell2mat //
////////////////

function M=mycell2mat(C)
// M=cat(1,C{:});
M = [];
for i=1:length(C)
  M = [M C(i)'];
end
endfunction

///////////////
// plot_tube //
///////////////

// plots a 'tube' along edges, with color given by D, and a possibly varying radius
// radius: constant or equal to nb points
// D(ata): empty or equal to nb points or nb segments
function h=plot_tube(P, T, D, radius, tubecolor)

h = [];
P = mycell2mat(P);

if (isempty(T)) then return; end;
it0  = T(1,1);  nT = size(T,2); nP = size(P,1); mdim=size(P,2);

if (mdim == 2) then 
  P = [P'; zeros(1,nP)]'; 
  mdim = 3; 
end // handle 2D slices

// convert radius to node data
if (length(radius)==1) then
  radius = radius*ones(1,nP); //radius(1)=0.5; radius($)=0.5;
elseif (length(radius)==nT) then
  radius = ([radius(1) radius(:)']+[radius(:)' radius($)])/2;
end

if (size(D,1) > 1) then error('only scalar data can be represented on a tube_plot'); end;
if (size(D,2)==nP) then
  point_data=1; 
else 
  point_data=0; 
end

nsubdiv = 20; 
ct   = cos((0:nsubdiv)*2*%pi/nsubdiv); 
ct($)= ct(1); 
st   = sin((0:nsubdiv)*2*%pi/nsubdiv);
st($)= st(1);
cnt  = 0;

h = [];

// Size P:  158.    2.  
// Size T:   2.    120. 
while (1)
  // search for consecutive edge points
  it1 = it0;
  //while (it1 < nT & T(1,it1+1) == T(2,it1)) it1 = it1+1; end;
  while (it1 < nT & T(1,it1+1) == T(2,it1)), it1 = it1+1; end;
  //disp(sprintf('sequence: %d - %d -- [%d-%d] - [%d-%d]',it0,it1,T(1,it0),T(2,it0),T(1,it1),T(2,it1)))
  // extract the sequence of points
  ip = [T(1,it0) T(2,it0:it1)];
  p = P(ip,:);      // P(:,ip)
  if (length(D)) then
    if (point_data) then 
      d = D(ip); 
    else 
      d = D(it0:it1); 
    end
  end
  nseg = it1-it0+1;
  // compute the normals of edges
  // normals = zeros(3, 2, nseg); // produce a hypermat
  normals = [];
  tang = p(2:$,:) - p(1:$-1,:); 
  for i=1:size(tang,1)
    tang(i,:) = tang(i,:) / max(%eps,sqrt(sum(tang(i,:).^2))); 
  end
  for i=1:nseg
    normals(:,:,i) = null_space(tang(i,:)); // won't be ok if normals have an
                                            // important rotation from a segment
                                            // to another - VERY PROBABLE BUG!!!      
  end
  X = zeros(mdim,nsubdiv+1,length(ip));
  for i=1:length(ip),
    if (i == 1) then
      n = normals(:,:,i)'; 
    elseif (i == length(ip)) then
      n = normals(:,:,$)';
    else
      n = ((normals(:,:,i-1)+normals(:,:,i))/2)';
    end
    for k=1:nsubdiv+1
      X(:,k,i) = (p(i,:) + radius(ip(i))*(n(1,:)*ct(k) + n(2,:)*st(k)))';
    end;
  end;
  if (length(D)) then
    C = repmat(d,nsubdiv+1,1);
    surf(squeeze(X(1,:,:)), squeeze(X(2,:,:)), squeeze(X(3,:,:)),C); // 'linestyle','none','FaceColor','interp')];
    h($+1) = gce();
    h($).thickness  = 0; // corresponds to linestyle none
    h($).color_flag = 3; // 2 -> flat shadding 3 -> interpolated shadding
  else
    //surf(squeeze(X(1,:,:)), squeeze(X(2,:,:)), squeeze(X(3,:,:))); // 'linestyle','none','facecolor',tubecolor)];
    // Workaround for bug 4042
    MyX1 = squeeze(X(1,:,:));
    MyX1 = matrix(MyX1.entries,double(MyX1.dims));
    MyX2 = squeeze(X(2,:,:));
    MyX2 = matrix(MyX2.entries,double(MyX2.dims));
    MyX3 = squeeze(X(3,:,:));
    MyX3 = matrix(MyX3.entries,double(MyX3.dims));

    //surf(squeeze(X(1,:,:)), squeeze(X(2,:,:)), squeeze(X(3,:,:)),'edgeco','cya'); // 'linestyle','none','facecolor',tubecolor)];
    surf(MyX1, MyX2, MyX3,'edgeco','cya'); // 'linestyle','none','facecolor',tubecolor)];
    h($+1) = gce();
    h($).thickness  = 0; // corresponds to linestyle none
    h($).color_mode = color(tubecolor);
    h($).color_flag = 0;
  end
  cnt = cnt+1;
  it0 = it1+1;
  if (it0 > nT) then return; end;
end
endfunction

////////////////
// do_plot_2D //
////////////////

// draw faces
function [hfaces,hmesh,hquiver] = do_plot_2D(sl,P,T,opt)

hfaces  = []; 
hmesh   = []; 
hquiver = [];
mdim    = length(P);

[o_data,err]                  = get_param(opt,'data',[]); // data to be plotted on the slice (on slice nodes)
[o_convex_data,err]           = get_param(opt,'convex_data',[]); // data to be plotted (given on the mesh convexes)
[o_msh,err]                   = get_param(opt,'mesh','auto'); // show the mesh ?
[o_msh_edges,err]             = get_param(opt,'mesh_edges','on'); // show mesh edges ?
[o_msh_edges_color,err]       = get_param(opt,'mesh_edges_color',[.6 .6 1]);
[o_msh_edges_width,err]       = get_param(opt,'mesh_edges_width',.7);
[o_msh_slice_edges,err]       = get_param(opt,'mesh_slice_edges','on');
[o_msh_slice_edges_color,err] = get_param(opt,'mesh_slice_edges_color',[.7 0 0]);
[o_msh_slice_edges_width,err] = get_param(opt,'mesh_slice_edges_width',.5);
[o_msh_faces,err]             = get_param(opt,'mesh_faces','off'); // fill mesh faces (otherwise they are transparent)
[o_msh_faces_color,err]       = get_param(opt,'mesh_faces_color',[.75 .75 .75]);
[o_pcolor,err]                = get_param(opt,'pcolor','on'); // if the field is scalar, a color plot of its values is plotted
[o_quiver,err]                = get_param(opt,'quiver','on'); // if the field is vector, represent arrows 	 
[o_quiver_density,err]        = get_param(opt,'quiver_density',50); // density of arrows in quiver plot
[o_quiver_scale,err]          = get_param(opt,'quiver_scale',1); // scaling of arrows (0=>no scaling)
[o_tube,err]                  = get_param(opt,'tube','on'); // use tube plot for linear parts of the slice
[o_tube_color,err]            = get_param(opt,'tube_color','red'); // color of tubes (ignored if 'data' is not empty)
[o_tube_radius,err]           = get_param(opt,'tube_radius',0.05); // tube radius; you can use a constant, or a percentage (of the mesh size) or a vector of nodal values
[o_showoptions,err]           = get_param(opt,'showoptions','off'); // list options used

d = mlist(['dlist','FaceVertexCData','FaceColor'],[],[]);  
d_is_set = %F;

if (length(T)) then
  if (ison(o_pcolor) & size(o_data,1)==1 & ~isempty(o_data)) then
    d('FaceVertexCData') = o_data(:);
    d('FaceColor') = 'interp';
    d_is_set = %T;
  elseif (isempty(o_data) & ison(o_msh_faces)) then
    d('FaceVertexCData') = o_msh_faces_color;
    d('FaceColor') = 'flat';
    d_is_set = %T;
  end  
  if (d_is_set) then
    //hfaces = patch('Vertices',mycell2mat(P)','Faces',T',d(:), 'EdgeColor','none'); // YC:
    p_tmp = mycell2mat(P);

    h = gcf();
    ctmp = d('FaceVertexCData');
    if (size(ctmp,1)==1) then // just one color
      ctmp = ones(size(T,2),1) * color(round(255*ctmp(1)), ...
                                       round(255*ctmp(2)), ...
                                       round(255*ctmp(3)));
    else
      ctmp = ceil((ctmp - min(ctmp)) / max(%eps,(max(ctmp) - min(ctmp))) * size(h.color_map,1));
      ctmp = matrix(ctmp(T,1),size(T,1),length(p_tmp(T,2))/size(T,1))';
    end
  
    if (size(p_tmp,2)==2) then
      xtmp = matrix(p_tmp(T,1),size(T,1),length(p_tmp(T,2))/size(T,1))';
      ytmp = matrix(p_tmp(T,2),size(T,1),length(p_tmp(T,2))/size(T,1))';
      ztmp = matrix(ones(p_tmp(T,2)),size(T,1),length(p_tmp(T,2))/size(T,1))';
      plot3d(xtmp', ytmp', list(ztmp',ctmp'));
    else
      xtmp = matrix(p_tmp(T,1),size(T,1),length(p_tmp(T,1))/size(T,1))';
      ytmp = matrix(p_tmp(T,2),size(T,1),length(p_tmp(T,1))/size(T,1))';
      ztmp = matrix(p_tmp(T,3),size(T,1),length(p_tmp(T,1))/size(T,1))';
      plot3d(xtmp', ytmp', list(ztmp',ctmp'));
    end
    hfaces = gce();
    select d('FaceColor')
      case 'interp' then  hfaces.color_flag = 3;
      case 'flat'   then  hfaces.color_flag = 2;
    end
    hfaces.thickness  = 0; ///o_msh_edges_width;
    hfaces.line_style = 1;
    hfaces.foreground = color(round(255*o_msh_edges_color(1)), ...
                              round(255*o_msh_edges_color(2)), ...
                              round(255*o_msh_edges_color(3)));
    hfaces.hiddencolor = -1;
  end
  if (ison(o_quiver)) then
    if (size(o_data,1)>1) then
      hquiver = do_quiver_plot(P,o_data,opt);
    end
  end 
end
if (ison(o_msh) & (ison(o_msh_edges) | ison(o_msh_slice_edges))) then
  [p,t1,t2] = gf_slice_get(sl,'edges');
  if (ison(o_msh_edges)) then
    // p: 2 x 1661
    // t1: 2 x 1760
    // t2: 0

    p = p';
    if (size(p,2)==2) & size(t1,1)~=0 then // 2D plot
      xtmp = matrix(p(t1,1),size(t1,1),length(p(t1,1))/size(t1,1))';
      ytmp = matrix(p(t1,2),size(t1,1),length(p(t1,1))/size(t1,1))';
      ztmp = matrix(ones(p(t1,2)),size(t1,1),length(p(t1,1))/size(t1,1))';
      plot3d(xtmp', ytmp', ztmp');
    elseif size(t1,1)~=0 then
      xtmp = matrix(p(t1,1),size(t1,1),length(p(t1,1))/size(t1,1))';
      ytmp = matrix(p(t1,2),size(t1,1),length(p(t1,1))/size(t1,1))';
      ztmp = matrix(p(t1,3),size(t1,1),length(p(t1,1))/size(t1,1))';
      plot3d(xtmp', ytmp', ztmp');
    end
    hmesh = gce();
    hmesh.thickness = o_msh_edges_width;
    hmesh.line_style = 1;
    hmesh.foreground = color(round(255*o_msh_edges_color(1)), ...
                             round(255*o_msh_edges_color(2)), ...
                             round(255*o_msh_edges_color(3)));
    hmesh.hiddencolor = -1;
    p = p'; //t1 = t1';
  end
  if (ison(o_msh_slice_edges)) then
    //hmesh = [hmesh patch('Vertices',p','Faces',t2','EdgeColor',o_msh_slice_edges_color,'LineWidth',o_msh_slice_edges_width)]; // YC
    p = p'; //t2 = t2';
    if (size(p,2)==2) & size(t2,1)~=0 then
      xtmp = matrix(p(t2,1),size(t2,1),length(p(t2,1))/size(t2,1))';
      ytmp = matrix(p(t2,2),size(t2,1),length(p(t2,1))/size(t2,1))';
      ztmp = matrix(ones(p(t2,2)),size(t2,1),length(p(t2,1))/size(t2,1))';
      plot3d(xtmp', ytmp', ztmp');
    elseif size(t2,1)~=0 then
      xtmp = matrix(p(t2,1),size(t2,1),length(p(t2,1))/size(t2,1))';
      ytmp = matrix(p(t2,2),size(t2,1),length(p(t2,1))/size(t2,1))';
      ztmp = matrix(p(t2,3),size(t2,1),length(p(t2,1))/size(t2,1))';
      plot3d(xtmp', ytmp', ztmp');
    end
    
    hmesh_tmp = gce();
    hmesh_tmp.thickness = o_msh_slice_edges_width;
    hmesh_tmp.line_style = 1;
    hmesh_tmp.foreground = color(round(255*o_msh_slice_edges_color(1)), ...
                                 round(255*o_msh_slice_edges_color(2)), ...
                                 round(255*o_msh_slice_edges_color(3)));
    hmesh.hiddencolor = -1;
    p = p'; //t2 = t2';
    hmesh = [hmesh hmesh_tmp];
  end
end 
endfunction

////////////////////
// do_quiver_plot //
////////////////////

// arrow plot
function hquiver = do_quiver_plot(P,U,opt)

[o_data,err]                  = get_param(opt,'data',[]); // data to be plotted on the slice (on slice nodes)
[o_convex_data,err]           = get_param(opt,'convex_data',[]); // data to be plotted (given on the mesh convexes)
[o_msh,err]                   = get_param(opt,'mesh','auto'); // show the mesh ?
[o_msh_edges,err]             = get_param(opt,'mesh_edges','on'); // show mesh edges ?
[o_msh_edges_color,err]       = get_param(opt,'mesh_edges_color',[.6 .6 1]);
[o_msh_edges_width,err]       = get_param(opt,'mesh_edges_width',.7);
[o_msh_slice_edges,err]       = get_param(opt,'mesh_slice_edges','on');
[o_msh_slice_edges_color,err] = get_param(opt,'mesh_slice_edges_color',[.7 0 0]);
[o_msh_slice_edges_width,err] = get_param(opt,'mesh_slice_edges_width',.5);
[o_msh_faces,err]             = get_param(opt,'mesh_faces','off'); // fill mesh faces (otherwise they are transparent)
[o_msh_faces_color,err]       = get_param(opt,'mesh_faces_color',[.75 .75 .75]);
[o_pcolor,err]                = get_param(opt,'pcolor','on'); // if the field is scalar, a color plot of its values is plotted
[o_quiver,err]                = get_param(opt,'quiver','on'); // if the field is vector, represent arrows 	 
[o_quiver_density,err]        = get_param(opt,'quiver_density',50); // density of arrows in quiver plot
[o_quiver_scale,err]          = get_param(opt,'quiver_scale',1); // scaling of arrows (0=>no scaling)
[o_tube,err]                  = get_param(opt,'tube','on'); // use tube plot for linear parts of the slice
[o_tube_color,err]            = get_param(opt,'tube_color','red'); // color of tubes (ignored if 'data' is not empty)
[o_tube_radius,err]           = get_param(opt,'tube_radius',0.05); // tube radius; you can use a constant, or a percentage (of the mesh size) or a vector of nodal values
[o_showoptions,err]           = get_param(opt,'showoptions','off'); // list options used

hquiver  = [];
P        = mycell2mat(P)'; 
mdim     = size(P,1);
qdim     = size(U,1);
nP       = size(P,2);
ptlst    = 1:nP;
bmin     = min(P);
bmax     = max(P);
xyscale  = max(bmax-bmin);
qradius2 = (xyscale/o_quiver_density)^2;
vscale   = max(max(abs(U)));
qlst     = [];
rm       = [];

while (length(ptlst)>0)
  ii   = ptlst(1); 
  qlst = [qlst ii];
  x    = P(1,ii);
  y    = P(2, ii); 
  if (mdim == 2) then
    rm = (find((P(1,:)-x).^2 + (P(2,:)-y).^2 < qradius2));
  elseif (mdim == 3) then
    z = P(3,ii);
    rm = (find((P(1,:)-x).^2 + (P(2,:)-y).^2 + (P(3,:)-z).^2 < qradius2));
  end
  if (length(rm)==0) then error('internal error in gf_plot'); end;
  ptlst = _setdiff(ptlst, rm);
end
if (qdim == 2) then
  nx = ones(1,2*length(P(1,qlst)));
  ny = ones(1,2*length(P(1,qlst)));
  deltaUmax = max(U(:,qlst),'c') - min(U(:,qlst),'c');
  deltaP    = max(P(:,qlst),'c') - min(P(:,qlst),'c');
  nx(1:2:$) = P(1,qlst);
  nx(2:2:$) = P(1,qlst) + 0.025 * deltaP(1) .* U(1,qlst) / max(%eps,norm(deltaUmax));
  ny(1:2:$) = P(2,qlst);
  ny(2:2:$) = P(2,qlst) + 0.025 * deltaP(2) .* U(2,qlst) / max(%eps,norm(deltaUmax));
  xarrows(nx,ny);
  a = gca();
  a.data_bounds = [min(P,'c')';max(P,'c')'];
  hquiver = gce();
  hquiver.arrow_size = o_quiver_scale;
else
  champ3(P(1,qlst),P(2,qlst),P(3,qlst),U(1,qlst),U(2,qlst),U(3,qlst)); // green
  hquiver = gce();
  hquiver.arrow_size = o_quiver_scale;
end
endfunction

