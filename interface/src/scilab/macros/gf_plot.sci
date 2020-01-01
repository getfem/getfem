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

function [hsurf, hcontour, hquiver, hmesh, hdefmesh]=gf_plot(mf,U,varargin)
// function h=gf_plot(mf,U,varargin)
// this function plots a 2D or 3D finite elements field.
//
// The options are specified as pairs of 'option name'/'option value'
//  
//  'zplot',{'off'|'on'}           : only for qdim=1, mdim=2
//  'norm', {'off'|'on'}           : if qdim >= 2, color-plot the norm of the field
//  'dir',[]	                      : or the scalar product of the field with 'dir' 
//                                   (can be a vector, or 'x', 'y' etc..)
//  'refine',8		      : nb of refinments for curved edges and surface plots
//  'interpolated',{'off'|'on'}    : if triangular patch are interpolated
//  'pcolor',{'on'|'off'}          : if the field is scalar, a color plot of its values is plotted
//  'quiver',{'on'|'off'}          : if the field is vector, represent arrows 	       
//  'quiver_density',50            : density of arrows in quiver plot
//  'quiver_scale',1               : scaling of arrows (0=>no scaling)
//  'mesh',{'off'|'on'}	           : show the mesh ?
//  'meshopts',{cell(0)}	          : cell array of options passed to gf_plot_slice for the mesh 
//  'deformed_mesh', {'off'|'on'}  : shows the deformed mesh (only when qdim == mdim)
//  'deformed_meshopts', {cell(0)} : cell array of options passed to gf_plot_slice 
//                                   for the deformed mesh 
//  'deformation',[]	              : plots on the deformed object (only when qdim == mdim)
//  'deformation_mf',[]            : plots on the deformed object (only when qdim == mdim)
//  'deformation_scale',0.1'       : indicate the amplitude of the deformation. Can be 
//                                   an absolute value if given as a number
//  'cvlst',[]		                   : list of convexes to plot (empty=>all convexes)
//  'title',[]                     : set the title
//  'contour',[]                   : list of contour values
//  'disp_options','off'           : displays the options passed to gf_plot

hsurf    = [];
hcontour = list();
hquiver  = [];
hmesh    = [];
hdefmesh = [];

opts = build_options_list(varargin(:));

if has_field(opts,'disp_options') then
  if ison(opts('disp_options')) then
    disp(opts);
  end
end

try 
  gf_workspace('push');
  [hsurf, hcontour, hquiver, hmesh, hdefmesh] = gf_plot_aux(mf,U,opts);
catch
  [str,n,line,func]=lasterror();
  disp('error in gf_plot: ' + str);
  disp(sprintf('error %d in %s at line %d\n', n, func, line));
  error('');
end
gf_workspace('pop');
endfunction

function [hsurf, hcontour, hquiver, hmesh, hdefmesh] = gf_plot_aux(mf,U,opts)

[nargout,nargin] = argn();

if nargin<2 then
  error('Too few input arguments')
end

hsurf    = [];
hcontour = list();
hquiver  = [];
hmesh    = [];
hdefmesh = [];
  
try
  qdim = gf_mesh_fem_get(mf, 'qdim');
  mdim = gf_mesh_get(mf, 'dim'); mdim3=mdim*3;
catch
  [str,n,line,func]=lasterror();
  disp('invalid mesh_fem: ' + str);
  disp(sprintf('error %d in %s at line %d\n', n, func, line));
  error('');
end

if (mdim == 1) then
  hsurf = gf_plot_1D(mf,U,opts);
  return;
end

if (mdim ~= 2 & mdim ~= 3) then
  error('only 2D and 3D mesh are handled by this function');
end

[o_zplot,err]             = get_param(opts,'zplot',             'off');  // only for qdim=1, mdim=2
[o_norm,err]              = get_param(opts,'norm',              'off');  // if qdim >= 2, color-plot the norm of the field
[o_dir,err]               = get_param(opts,'dir',               []);     // or the scalar product of the field with 'dir' (can be a vector, or 'x', 'y' etc..)
[o_refine,err]            = get_param(opts,'refine',            8);      // nb of refinments for curved edges and surface plots
[o_interpolated,err]      = get_param(opts,'interpolated',      'off');  //if triangular patch are interpolated
[o_pcolor,err]            = get_param(opts,'pcolor',            'on');   // if the field is scalar, a color plot of its values is plotted
[o_quiver,err]            = get_param(opts,'quiver',            'on');   // if the field is vector, represent arrows 	       
[o_quiver_density,err]    = get_param(opts,'quiver_density',    50);     // density of arrows in quiver plot
[o_quiver_scale,err]      = get_param(opts,'quiver_scale',      1);      // scaling of arrows (0=>no scaling)
[o_mesh,err]              = get_param(opts,'mesh',              'off');  // show the mesh ?
[o_meshopts,err]          = get_param(opts,'meshopts',          list()); // cell array of options passed to gf_plot_slice for the mesh 
[o_deformed_mesh,err]     = get_param(opts,'deformed_mesh',     'off');  // shows the deformed mesh (only when qdim == mdim)
[o_deformed_meshopts,err] = get_param(opts,'deformed_meshopts', list()); // cell array of options passed to gf_plot_slice for the deformed mesh 
[o_deformation,err]       = get_param(opts,'deformation',       []);     // plots on the deformed object (only when qdim == mdim)
[o_deformation_mf,err]    = get_param(opts,'deformation_mf',    []);     // plots on the deformed object (only when qdim == mdim)
[o_deformation_scale,err] = get_param(opts,'deformation_scale', '10%');    // indicate the amplitude of the deformation. Can be a percentage of the mesh width if given as a string, or an absolute value if given as a number
[o_cvlst,err]             = get_param(opts,'cvlst',             []);     // list of convexes to plot
[o_title,err]             = get_param(opts,'title',             []);
[o_contour,err]           = get_param(opts,'contour',           []);
[o_mesh_level_set,err]    = get_param(opts,'mesh_level_set',    []);     // list of contour values

if (ison(o_zplot)) then
  if (mdim ~= 2) then
    error('zplot allowed only on 2D scalar mesh_fems');
  else 
    o_interpolated = 'on'; // or patch won't work
  end
end

is_scalarplot = (ison(o_norm) + ~isempty(o_dir));
if (is_scalarplot > 1) then
  error('only one occurence of the options ''norm'' and ''dir'' is allowed');
end

if (typeof(o_dir)=='string') then
  v = zeros(1,qdim);
  if (o_dir=='x') then
    v(1)=1;
  elseif (o_dir=='y') then
    v(2)=1;
  elseif (o_dir=='z') then
    v(3)=1;
  else error('wrong direction');
  end
  o_dir=v;
end

scalarplot_dir=o_dir(:);
if (qdim == 1) then 
  is_scalarplot  = 1; 
  scalarplot_dir = 1;
end
if (~isempty(o_contour) & ~is_scalarplot) then
  error('contour plot has no meaning for a vector field');
end

mfdef = mf;
if (~isempty(o_deformation_mf)) then mfdef = o_deformation_mf; end;
dqdim = gf_mesh_fem_get(mfdef,'qdim');
if (~isempty(o_deformation') | ison(o_deformed_mesh)) then
  if (mdim ~= dqdim & ~ison(o_zplot)) then
    error(sprintf('can''t plot the deformation of an %dD-object by a %dD-field',mdim,dqdim));
  end
end

if (isempty(o_cvlst)) then o_cvlst = gf_mesh_get(mf,'cvid'); end;

nbdof = gf_mesh_fem_get(mf,'nbdof');
if (nbdof <= 0) then error('invalid finite element mf argument'); end

if (length(U) ~= nbdof) then
  error('wrong dimensions for U, should have ' + string(nbdof) + ' columns'); 
end

if (ison(o_zplot) | mdim == 3) then
  h = gca();
  h.view = '3d';
else 
  h = gca();
  h.view = '2d';
end

// build the slice object
try
  if (~isempty(o_mesh_level_set)) then
    sl = gf_slice(list('none'),o_mesh_level_set,o_refine,o_cvlst);
  elseif (gf_mesh_fem_get(mf, 'has_linked_mesh_levelset')) then
    sl = gf_slice(list('none'),gf_mesh_fem_get(mf,'linked_mesh_levelset'),o_refine,o_cvlst);
  else
    sl = gf_slice(list('none'),gf_mesh_fem_get(mf,'linked mesh'),o_refine,o_cvlst);
  end
catch
  [str,n,line,func]=lasterror();
  disp('can''t build slice: ' + str);
  disp(sprintf('error %d in %s at line %d\n', n, func, line));
  error('');
end

try
  Usl = gf_compute(mf,U,'interpolate on',sl);
catch
  [str,n,line,func]=lasterror();
  disp('can''t interpolate on slice: ' + str);
  disp(sprintf('error %d in %s at line %d\n', n, func, line));
  error('');
end
Psl = gf_slice_get(sl,'pts');

// plot the original mesh
if (ison(o_mesh)) then
  hmesh = list(gf_plot_slice(sl,'mesh','on',o_meshopts(:)));
end

// apply the optional deformation
if (~isempty(o_deformation) | mfdef('id') ~= mf('id')) then
  ida = gf_mesh_fem_get(mfdef,'linked mesh');
  idb = gf_mesh_fem_get(mf,'linked mesh');
  if (ida.id ~= idb.id) then
    error('the deformation mesh_fem and the data mesh_fem do not seem to share the same mesh');
  end
  if (~isempty(o_deformation)) then 
    Udef = o_deformation;
  else 
    Udef = U;
  end
  Pdef  = gf_compute(mfdef, Udef, 'interpolate on', sl);
  if (isnumeric(o_deformation_scale)) then 
    dscale = o_deformation_scale;
  elseif (typeof(o_deformation_scale)=='string' & length(o_deformation_scale) & part(o_deformation_scale,length(o_deformation_scale))=='%') then
    dscale = eval(part(o_deformation_scale,1:length(o_deformation_scale)-1));
    mwidth = max([max(Psl,'c') - min(Psl,'c')],'r');
    defmax = max(abs(Pdef(:)));
    if (defmax) then
      dscale = dscale * 0.01 * mwidth / defmax;
    end
  else 
    error('wrong value for deformation_scale: should be a number or a percentage in a string');
  end
  Psl = Psl + Pdef*dscale;
  gf_slice_set(sl,'pts', Psl);    
  clear Pdef Udef dscale mwidth defmax;
end

if (is_scalarplot) then
  // compute scalar values if necessary
  if (ison(o_norm)) then
    sV = sqrt(sum(Usl.*conj(Usl),1));
  else
    sV = scalarplot_dir(:)'*Usl;
  end
  // and optionally apply the zplot deformation
  if (ison(o_pcolor) & ison(o_zplot) & is_scalarplot) then
    Psl = [Psl;sV]; gf_slice_set(sl,'pts',Psl);
  end
end

// plot the deformed mesh
if (ison(o_deformed_mesh)) then
  hdefmesh = list(gf_plot_slice(sl,'mesh','on', o_deformed_meshopts(:))); // YC: ??
end

if (is_scalarplot) then
  // plot the 'surfacic' part
  if (ison(o_pcolor)) then
    gf_plot_slice(sl,'mesh','off','data',sV);
  end
  
  // basic contour plot (should work also on 3D surfaces)
  contour_colors = [0.9 0.0 0.0; 
                    0.0 0.8 0.0;
                    0.0 0.0 0.6;
                    0.6 0.6 0.0;
                    0.7 0.0 0.7;
                    0.0 0.7 0.9]; 
                    
  hcontour = list();
  for cnum=1:length(o_contour)
    c = o_contour(cnum);
    slC = gf_slice(list('isovalues',0,mf,U,c),sl);
    [a,b,c,hcontour(cnum)] = gf_plot_slice(slC,'tube','off','mesh','off');
    
    if (~isempty(hcontour(cnum))) then
      //printf('here isovalue');
      //disp(hcontour(cnum).parent)
      //hcontour(cnum).children.color_mode = color(round(255*contour_colors(modulo(cnum,size(contour_colors,1))+1,1)), ...
      //                                           round(255*contour_colors(modulo(cnum,size(contour_colors,1))+1,2)), ...
      //                                           round(255*contour_colors(modulo(cnum,size(contour_colors,1))+1,3)));
      //hcontour(cnum).children.color_flag = 0;
      hcontour(cnum).children.thickness = 1;
//      set(hcontour(cnum),...
//          'Color',contour_colors(modulo(cnum,size(contour_colors,1))+1,:),...
//          'LineStyle',contour_linestyle(modulo(cnum,length(contour_linestyle))+1),...
//          'LineWidth',1);
    end
    gf_delete(slC);
  end
else
  [a,b,hquiver,c] = gf_plot_slice(sl,'data',Usl,'quiver','on','quiver_density',o_quiver_density,'quiver_scale',o_quiver_scale);
end
endfunction

