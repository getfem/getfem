function [hline, hdof] = gf_plot_1D(mf,U, varargin)
// function h=gf_plot_1D(mf,U,...)
// this function plots a 1D finite elements field.
//
// The options are specified as pairs of 'option name'/'option value'
//  'style', 'bo-'       : the line style and dof marker style (same
//                         syntax as in the matlab command 'plot').
//  'color', []          : override the line color.
//  'dof_color', [1,0,0] : color of the markers for the degrees of freedom.
//  'width', 2           : line width.

opts = build_options_list(varargin);

try 
  gf_workspace('push', 'gf_plot_1D');
  [hline, hdof] = gf_plot_1D_aux(mf,U, opts);
catch
  disp('error in gf_plot_1D : ' + lasterror());
end;
gf_workspace('pop');
endfunction

  
function [hline, hdof] = gf_plot_1D_aux(mf, U, opts)
[opt_style,err]      = get_param(opts,'style','bo-');
[opt_dof_marker,err] = get_param(opts,'dof_marker','');
[opt_width,err]      = get_param(opts,'width',2);
[opt_dof_color,err]  = get_param(opts,'dof_color',[1 0 0]);

// try to remove markers from the line style
s              = opt_style; 
opt_style      = ''; 
opt_dof_marker = '';

for i=s
  if (isempty(strindex('ox+*sdv^<>ph',i))) then 
    opt_style = [opt_style i];
  else
    opt_dof_marker = i;
  end
end

// save graphical context
cax = gcf();

nbd = gf_mesh_fem_get(mf, 'nbdof');
if (nbd < 100) then
  REFINE = 32;
elseif (nbd < 1000) then
  REFINE = 6;
else
  REFINE = 2;
end

m   = gf_mesh_fem_get(mf, 'linked_mesh');
sl  = gf_slice(list('none'),m, REFINE); 
Usl = gf_compute(mf,U,'interpolate on', sl);
D   = unique(gf_mesh_fem_get(mf, 'dof nodes'));
slD = gf_slice('points', m, D);
UD  = gf_compute(mf,U,'interpolate on',slD);

X = gf_slice_get(sl, 'pts');
Y = Usl;
plot(X, Y, opt_style); 
hline = gce();
hline.children.thickness = opt_width;
if (~isempty(opt_color)) then 
  hline.children.line_style = opt_color;
end

hdof = [];
if (~isempty(opt_dof_marker)) then
  hdof = plot2d(gf_slice_get(slD, 'pts'), UD, opt_dof_marker); // opt_dof_marker must be < 0 
end
endfunction

