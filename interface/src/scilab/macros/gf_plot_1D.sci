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

function [hline, hdof] = gf_plot_1D(mf,U, varargin)
// function h=gf_plot_1D(mf,U,...)
// this function plots a 1D finite element field.
//
// Available options are specified as pairs of 'option name'/'option value'
//  'style', 'bo-'       : line style and dof marker style (same
//                         syntax as in the Scilab command 'plot');
//  'color', ''          : override line color (by a given color name);
//  'dof_color', ''      : override color of dof markers;
//  'width', 2           : line width.

opts = build_options_list(varargin(:));

try 
  gf_workspace('push', 'gf_plot_1D');
  [hline, hdof] = gf_plot_1D_aux(mf,U, opts);
catch
  [str,n,line,func]=lasterror();
  disp('error in gf_plot_1D: ' + str);
  disp(sprintf('error %d in %s at line %d\n', n, func, line));
  error('');
end
gf_workspace('pop');
endfunction

  
function [hline, hdof] = gf_plot_1D_aux(mf, U, opts)

[opt_style,err]      = get_param(opts,'style','bo-');
[opt_color,err]      = get_param(opts,'color','');
[opt_dof_color,err]  = get_param(opts,'dof_color','');
[opt_width,err]      = get_param(opts,'width',2);

// remove eventual markers from the line style
s              = opt_style; 
opt_style      = ''; 
opt_dof_marker = '';

for i = 1:length(s)
  if (isempty(strindex('ox+*.sdv^<>p', part(s, i)))) then
    opt_style = opt_style + part(s, i);
  elseif i == 1 then
    opt_dof_marker = part(s, i);
  elseif '-.' == part(s, [i-1,i]) then
    opt_style = opt_style + '.';
  else opt_dof_marker = part(s, i);
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

m = gf_mesh_fem_get(mf, 'linked_mesh');
sl = gf_slice(list('none'),m, REFINE); 
Usl = gf_compute(mf,U,'interpolate on', sl);
D = unique(gf_mesh_fem_get(mf, 'basic dof nodes'));
slD = gf_slice('points', m, D);
UD = gf_compute(mf,U,'interpolate on',slD);

X = gf_slice_get(sl, 'pts');
Y = Usl;
plot(X, Y, opt_style); 
hline = gce();
hline.children.thickness = opt_width;
if (~isempty(opt_color)) then
  hline.children.foreground = color(opt_color);
end

hdof = [];
if (~isempty(opt_dof_marker)) then
  // add color to the marker if it is given in opt_style
  for i = 1:length(opt_style)
    if (~isempty(strindex('rgbcmykw', part(opt_style, i)))) then
      opt_dof_marker = part(s, i) + opt_dof_marker;
    end
  end
  
  plot(gf_slice_get(slD, 'pts'), UD, opt_dof_marker);
  hdof = gce();
  if (~isempty(opt_color)) then
    hdof.children.mark_foreground = color(opt_dof_color);
  end
end
endfunction

