function [hmesh,hbound,hfill,hvert,hconv,hdof]=gf_plot_mesh(M, varargin)
// function [hmesh,hbound,hfill,hvert,hconv,hdof]=gf_plot_mesh(M, [,properties])
//                     [,'cvlst',CVLST] ['boundaries'[BLST]])
//   General mesh plotting function.
//  
//   H=gf_plot_mesh(M) displays a mesh.
//
//   properties are:
//    'vertices', {'off'|'on'}     displays also vertices numbers. 
//    'convexes', {'off'|'on'}     displays also convexes numbers. 
//    'dof',{'off'|'on'}           displays also finite element nodes.
//    'regions',BLST               displays the boundaries listed in BLST.
//    'cvlst',CVLST                display only the listed convexes. If
//   CVLST has two rows, display only the faces listed in the second row.
//    'edges', {'on' | 'off'}      display edges ?
//    'faces', {'off'|'on'}        fills each 2D-face of the mesh
//    'curved', {'off'|'on'}       displays curved edges
//    'refine',N                   refine curved edges and filled faces N times  
//    'deformation', Udef          optionnal deformation applied to the mesh (M must be a mesh_fem object)
//    'edges_color',[.6 .6 1]      RGB values for the color of edges
//    'edges_width',1              
//    'faces_color',[.75 .75 .75]) RGB values for the color of faces
//    'quality',{ 'off' | 'on' }   Display the quality of the mesh.
//
//   CAUTION:
//     For 'dof', M should be a mesh_fem identifier, 
//   not a simple mesh object.
//  
//   $Id: gf_plot_mesh.m 2282 2006-02-23 16:24:13Z pommier $

//   A. Huard, Y. Renard, J. Pommier 

printf('DEBUG: in gf_plot_mesh\n');

[nargout,nargin] = argn();

if nargin<1 then
  error('Too few input arguments')
end

opts = build_options_list(varargin(:));

hmesh  = [];
hbound = [];
hfill  = [];
hvert  = [];
hconv  = [];
hdof   = [];

mdim = gf_mesh_get(M,'dim');
if (mdim <= 2) then
  defaultref = 8;
else 
  defaultref = 4;
end

[o_vertices,err]    = get_param(opts,'vertices','off');
[o_convexes,err]    = get_param(opts,'convexes','off');
[o_dof,err]         = get_param(opts,'dof','off');
[o_regions,err]     = get_param(opts,'regions','');
[o_boundaries,err]  = get_param(opts,'boundaries','');
[o_cvlst,err]       = get_param(opts,'cvlst',[]);
[o_edges,err]       = get_param(opts,'edges','on');
[o_faces,err]       = get_param(opts,'faces','off');
[o_explode,err]     = get_param(opts,'explode',0);
[o_quality,err]     = get_param(opts,'quality','off');
[o_curved,err]      = get_param(opts,'curved','off');
[o_refine,err]      = get_param(opts,'refine',defaultref);
[o_deformation,err] = get_param(opts,'deformation',[]);
[o_edges_color,err] = get_param(opts,'edges_color',[.6 .6 1]);
[o_edges_width,err] = get_param(opts,'edges_width',.7);
[o_faces_color,err] = get_param(opts,'faces_color',[.75 .75 .75]);

if (length(o_boundaries) == 0) then
  o_boundaries = o_regions;
end

if (typeof(o_boundaries)=='string') then
  if (convstr(o_boundaries,'l')=='all') then
    o_boundaries = gf_mesh_get(M, 'boundaries');
  end
end

// init cvlst and cvflst
if (isempty(o_cvlst)) then
  cvlst  = gf_mesh_get(M,'cvid'); 
  cvflst = [cvlst; int32(zeros(1,length(cvlst)))]; // int32 is the type of cvlst
else 
  cvlst  = o_cvlst;
  cvflst = cvlst;  
  if (size(cvflst,1)==2) then
    cvlst = unique(cvlst(1,:));
  end
end

PXY = gf_mesh_get(M, 'pts');

E = gf_mesh_get(M, 'edges', cvflst);

if (~ison(o_curved) & isempty(o_deformation)) then
  X = [PXY(1,E(1,:)); PXY(1,E(2,:))];
  if (mdim == 1) then
    Y = zeros(size(X));
    PXY = [PXY; zeros(PXY)];
  elseif (mdim >= 2) then
    Y = [PXY(2,E(1,:)); PXY(2,E(2,:))];
    if (mdim == 3) then
      Z = [PXY(3,E(1,:)); PXY(3,E(2,:))];
    end
  end
else
  // here, mdim is always >= 2
  if (~isempty(o_deformation)) then
    vE = gf_compute(M,o_deformation,'mesh edges deformation',o_refine,cvflst);
  else
    vE = gf_mesh_get(M, 'curved edges', o_refine, cvflst);
  end
  ni = size(vE,2);
  ne = size(vE,3);
  X  = matrix(vE(1,:,:), [ni ne]);
  Y  = matrix(vE(2,:,:), [ni ne]);
  if (mdim == 3) then
    Z = matrix(vE(3,:,:), [ni ne]);
  end
end

// get the viewable pts id
PID = union(E(1,:),E(2,:));
if (mdim > 3) then error('sorry, only mesh of dimension <= 3 allowed'); end;
nbpts = size(PXY,2);

Bmax  = max(PXY(:,PID)','r')';
Bmin  = min(PXY(:,PID)','r')';
Bdiff = Bmax - Bmin;
Bdiff = Bdiff + (Bdiff == 0); // remplace 0.0 par 1.0
ecart = Bdiff/150;

if (ison(o_convexes)) then
  cv_center = zeros(max(mdim,2),length(cvlst));
  // find convexes centers
  [cv_pid, cv_idx] = gf_mesh_get(M, 'pid from cvid',cvlst);
  for i=1:length(cvlst)
    cv_center(:,i) = mean(PXY(:, cv_pid(double(cv_idx(i)):double(cv_idx(i+1))-1)),2);
  end
end

if (ison(o_dof)) then
  Q = gf_mesh_fem_get(M, 'qdim');
  dofid    = gf_mesh_fem_get(M, 'dof from cv', cvlst);
  [dofpos] = gf_mesh_fem_get(M, 'dof nodes', dofid);
  [keep]   = find([1 or(dofpos(:,2:$) ~= dofpos(:,1:$-1),1)]);
  dofmult  = [keep(2:$)-keep(1:$-1) size(dofpos,2)+1-keep($)];
  dofpos   = dofpos(:, keep); dofid = dofid(keep);
  if (mdim == 1) then dofpos = [dofpos; zeros(size(dofpos))]; end;
end

bedge = list();
for bnum=1:length(o_boundaries)
  cvf = gf_mesh_get(M, 'boundary', o_boundaries(bnum));
  
  bid = gf_mesh_get(M, 'edges', cvf, 'merge convex');
  if (bnum == 8) then disp(bid); end;
  
  bedge(bnum) = zeros(2, size(bid,2), mdim);
  
  for i=1:max(mdim,2)
    bedge(bnum)(:,:,i) = [PXY(i,bid(1,:)); PXY(i,bid(2,:))];
  end
end

// save graphical context
cax = gcf();

disp('plotting mesh...');

if (mdim <= 2) then
  if (ison(o_edges)) then
    plot(X, Y);
    hmesh = gce();
    hmesh.children(:).thickness = o_edges_width;
    hmesh.children(:).line_style = 0; // Continous lines
    hmesh.children(:).foreground = o_edges_color;
  end
  for bnum=1:length(o_boundaries),
    plot(bedge(bnum)(:,:,1), bedge(bnum)(:,:,2));
    hbound(bnum) = gce();
    hbound(bnum).children(:).thickness  = 2;
    hbound(bnum).children(:).line_style = 0; // Continous lines
    hbound(bnum).children(:).foreground = 5; // red
  end
  if (ison(o_vertices)) then
    for i=1:length(PID)
      xstring(PXY(1,PID(i))+ecart(1), PXY(2,PID(i))+ecart(2), string(double(PID(i)))); // 'HorizontalAlignment','center','VerticalAlignment','middle'
      hvert = gce();
      hvert.alignment = 'center';
    end
  end
  if (ison(o_convexes)) then
    for i=1:size(cv_center,2)
      xstring(cv_center(1,i), cv_center(2,i), string(double(cvlst(i)))); // 'HorizontalAlignment','center','VerticalAlignment','middle', 'Color', [.7 0 0]
      hconv = gce();
      hconv.alignment = 'center';
      hconv.font_foreground = 5; // Red
    end
  end
  if (ison(o_dof)) then
    hdof = zeros(length(dofid),1);
    for i=1:length(dofid),
      if (dofmult(i)==1) then 
        s=string(dofid(i)); 
      else 
        s=sprintf('%d*%d',dofid(i),dofmult(i)); 
      end
      xstring(dofpos(1,i)-ecart(1), dofpos(2,i)-ecart(2), s); /// 'HorizontalAlignment','center','VerticalAlignment','middle', 'Color', [0 .4 0]);
      hdof(i) = gce();
      hdof(i).alignment = 'center';
      hdof(i).font_foreground = 2; // Blue
    end
  end
else
  if (ison(o_edges)) then
    plot3d(X, Y, Z); // 'Color',o_edges_color,'LineWidth',o_edges_width
    hmesh = gce();
    hmesh.children(:).thickness  = o_edges_width;
    hmesh.children(:).line_style = 0; // Continuous line
    hmesh.children(:).foreground = o_edges_color; 
  end
  for bnum=1:length(o_boundaries),
    plot3d(bedge(bnum)(:,:,1), bedge(bnum)(:,:,2), bedge(bnum)(:,:,3)); // 'Color','red','LineWidth',2);
    hbound(bnum) = gce();
    hbound(bnum).children(:).thickness  = 2;
    hbound(bnum).children(:).line_style = 0; // Continuous line
    hbound(bnum).children(:).foreground = 5; // Red
  end
  if (ison(o_vertices)) then
    for i=1:length(PID)
      xstring(PXY(1,PID(i))+ecart(1), PXY(2,PID(i))+ecart(2), string(PID(i))); // 'HorizontalAlignment','center','VerticalAlignment','middle','Color', [.0 0 0]
      hvert = gce();
      hvert.data = [hvert.data PXY(3,PID(i))+ecart(3)]; // We add the 3rd component
      hvert.alignment = 'center';
      hvert.font_foreground = 1; // Black
    end
  end
  if (ison(o_convexes)) then
    for i=1:size(cv_center,2)
      xstring(cv_center(1,i), cv_center(2,i), string(cvlst(i))); // 'HorizontalAlignment','center','VerticalAlignment','middle','Color', [.7 0 0]
      hconv = gce();
      hconv.data = [hconv.data cv_center(3,i)]; // We add the 3rd component
      hconv.alignment = 'center';
      hconv.font_foreground = 5; // Red
    end
  end
  if (ison(o_dof)) then
    for i=1:size(dofpos,2)
      xstring(dofpos(1,i)-ecart(1), dofpos(2,i)-ecart(2), string(dofid(i))); // 'HorizontalAlignment','center','VerticalAlignment','middle' 'Color', [0 .4 0]);
      hdof = gce();
      hdof.data = [hdof.data dofpos(3,i)-ecart(3)]; // We add the 3rd component
      hdof.alignment = 'center';
      hdof.font_foreground = 3; // Green
    end
  end
end

if (ison(o_quality)) then
  q = gf_mesh_get(M,'quality', cvflst(1,:));
  qmf = gf_mesh_fem(M);
  gf_mesh_fem_set(qmf, 'classical fem', 0);
  [a,b] = gf_mesh_fem_get(qmf, 'dof from cvid', cvflst(1,:));
  Q = zeros(1, gf_mesh_fem_get(qmf, 'nbdof'));
  for k=1:length(b)-1
    Q(a(b(k))) = q(k);
  end
end

if (o_explode ~= 0) then
  sl = gf_slice(list('explode',o_explode),M,o_refine,cvflst);
  data = list();
  if (ison(o_quality)) then
    sQ = gf_compute(qmf,Q,'interpolate on',sl);
    data = list('data',sQ);
  end
  gf_plot_slice(sl, data(:),'mesh_faces',       o_faces, ...
                            'mesh_edges_color', o_edges_color, ...
                            'mesh_edges_width', o_edges_width, ...
                            'mesh_faces_color', o_faces_color);
  gf_delete(sl); // light;
elseif (ison(o_quality)) then
  gf_plot(qmf, Q, 'cvlst', cvflst); 
elseif (ison(o_faces)) then
  // should be replaced by a gf_plot_slice ..
  T = gf_mesh_get(M, 'triangulated surface', o_refine, cvflst);
  if (mdim == 2) then
    // YC: trouver un equivalent à patch
    //hfill=plot2d([T(1:mdim:(mdim*3),:) T(1,:)],[T(2:mdim:(mdim*3),:) T(2,:)], o_faces_color, flag = [-1 0 4]); //, 'Erasemode','normal','Edgecolor','none');
    plot3d(T(1:mdim:(mdim*3),:),T(2:mdim:(mdim*3),:), zeros(T(2:mdim:(mdim*3),:)), o_faces_color, flag = [-1 0 4]); //, 'Erasemode','normal','Edgecolor','none');
    hfill= gca();
    hfill.view='2d';
  else
    plot3d(T(1:mdim:(mdim*3),:),T(2:mdim:(mdim*3),:), T(3:mdim:(mdim*3),:), o_faces_color, flag = [-1 0 4]); //'Erasemode','normal','Edgecolor','none');
  end
end

if (ison(o_quality)) then
  gf_delete(qmf);
end
endfunction

