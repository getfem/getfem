function [hline, hdof] = gf_plot_1D(mf,U, varargin)
% function h=gf_plot_1D(mf,U,...)
% this function plots a 1D finite elements field.
%
% The options are specified as pairs of "option name"/"option value"
%  'style', 'bo-'       : the line style and dof marker style (same
%                         syntax as in the matlab command 'plot').
%  'color', []          : override the line color.
%  'dof_color', [1,0,0] : color of the markers for the degrees of freedom.
%  'width', 2           : line width.
%  Copyright (C) 1999-2017 Yves Renard
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

  try 
    gf_workspace('push', 'gf_plot_1D');
    [hline, hdof] = gf_plot_1D_aux(mf,U, varargin{:});
  catch
    disp(['error in gf_plot_1D : ' lasterr]);
  end;
  gf_workspace('pop');

  
function [hline, hdof] = gf_plot_1D_aux(mf, U, varargin)
  opt = struct('style', 'bo-',...
	       'color', [], ...
	       'dof_color', [1, 0, 0],...
	       'width', 2);
  
  % parse argument list
  opt = getopt(opt,varargin);
  
  % try to remove markers from the line style
  s = opt.style; opt.style = ''; opt.dof_marker = '';
  for i=s,
    if (isempty(strfind('ox+*sdv^<>ph',i))),
      opt.style = [opt.style i];
    else
      opt.dof_marker = i;
    end;
  end;
  
  % save graphical context
  cax = newplot;
  cfig = get(cax,'Parent');
  hold_state = ishold;
  ax_nextplot = lower(get(cax,'NextPlot'));
  fig_nextplot = lower(get(cfig,'NextPlot'));

  nbd=gf_mesh_fem_get(mf, 'nbdof');
  if (nbd < 100)
    REFINE = 32;
  elseif (nbd < 1000)
    REFINE = 6;
  else
    REFINE = 2;
  end;
  m=gf_mesh_fem_get(mf, 'linked_mesh');
  sl=gf_slice({'none'},m, REFINE); 
  Usl = gf_compute(mf,U,'interpolate on', sl);
  D = unique(gf_mesh_fem_get(mf, 'basic dof nodes'));
  slD = gf_slice('points', m, D);
  UD  = gf_compute(mf,U,'interpolate on',slD);
  
  X = gf_slice_get(sl, 'pts');
  Y = Usl;
  hline = plot(X, Y, opt.style); 
  set(hline, 'LineWidth', opt.width);
  if (~isempty(opt.color)), set(hline, 'Color', opt.color); end;
  %get(hline, 'Color')
  hdof = [];
  if (~isempty(opt.dof_marker)),
    hold on;
    hdof = scatter(gf_slice_get(slD, 'pts'), UD, opt.dof_marker);
    set(hdof, 'CData', opt.dof_color);
    hold off;
  end;
  
  if (~hold_state),
    set(cax,'XLimMode','auto');
    set(cax,'YLimMode','auto');
    set(cax,'ZLimMode','auto');
  end;
  
  % restore graphical context
  set(cax,'NextPlot',ax_nextplot);
  set(cfig,'NextPlot',fig_nextplot);
  
