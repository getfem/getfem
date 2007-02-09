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
  D = unique(gf_mesh_fem_get(mf, 'dof nodes'));
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
  