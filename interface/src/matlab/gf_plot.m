function [hsurf, hcontour, hquiver, hmesh, hdefmesh]=gf_plot(mf,U,varargin)
% function h=gf_plot(mf,U,varargin)
% this function plots a 2D or 3D finite elements field.
%
% The options are specified as pairs of "option name"/"option value"
%  
%  'zplot',{'off'|'on'}       : only for qdim=1, mdim=2
%  'norm', {'off'|'on'}       : if qdim >= 2, color-plot the norm of the field
%  'dir',[]	              : or the scalar product of the field with 'dir' 
%                               (can be a vector, or 'x', 'y' etc..)
%  'refine',8		      : nb of refinments for curved edges and surface plots
%  'interpolated',{'off'|'on'}: if triangular patch are interpolated
%  'pcolor',{'on'|'off'}      : if the field is scalar, a color plot of its values is plotted
%  'quiver',{'on'|'off'}      : if the field is vector, represent arrows 	       
%  'quiver_density',50        : density of arrows in quiver plot
%  'quiver_scale',1           : scaling of arrows (0=>no scaling)
%  'mesh',{'off'|'on'}	      : show the mesh ?
%  'meshopts',{cell(0)}	         : cell array of options passed to gf_plot_slice for the mesh 
%  'deformed_mesh', {'off'|'on'} : shows the deformed mesh (only when qdim == mdim)
%  'deformed_meshopts', {cell(0)}: cell array of options passed to gf_plot_slice 
%                                  for the deformed mesh 
%  'deformation',[]	      : plots on the deformed object (only when qdim == mdim)
%  'deformation_mf',[]        : plots on the deformed object (only when qdim == mdim)
%  'deformation_scale','10%'  : indicate the amplitude of the deformation. Can be 
%                               a percentage of the mesh width if given as a string, 
%                               or an absolute value if given as a number
%  'cvlst',[]		      : list of convexes to plot (empty=>all convexes)
%  'title',[]                 : set the title
%  'contour',[]               : list of contour values
%  'disp_options', {'off'|'on'} : shows the option or not.
%
%
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

  try 
    gf_workspace('push');
    [hsurf, hcontour, hquiver, hmesh, hdefmesh]=gf_plot_aux(mf,U,varargin{:});
  catch
    disp(['error in gf_plot : ' lasterr]);
  end;
  gf_workspace('pop');
  

function [hsurf, hcontour, hquiver, hmesh, hdefmesh]=gf_plot_aux(mf,U,varargin)
  if nargin<2,
    error('Too few input arguments')
  end

  hsurf=[];
  hcontour={};
  hquiver=[];
  hmesh=[];
  hdefmesh=[];
    
  try
    mf=struct(mf);
    qdim = gf_mesh_fem_get(mf, 'qdim');
    mdim = gf_mesh_get(mf, 'dim'); mdim3=mdim*3;
  catch
    error(['invalid mesh_fem ? -- ' lasterr]);
  end;

  if (mdim == 1),
    hsurf = gf_plot_1D(mf,U,varargin{:});
    return;
  end;
  
  if (mdim ~= 2 & mdim ~= 3),
    error('only 2D and 3D mesh are handled by this function');
  end;
  
  opt = struct('zplot','off',...   % only for qdim=1, mdim=2
	       'norm','off',...     % if qdim >= 2, color-plot the norm of the field
	       'dir',[],...       % or the scalar product of the field with 'dir' (can be a vector, or 'x', 'y' etc..)
	       'refine',8,...         % nb of refinments for curved edges and surface plots
	       'interpolated','off',... %if triangular patch are interpolated
	       'pcolor','on',... % if the field is scalar, a color plot of its values is plotted
	       'quiver','on',... % if the field is vector, represent arrows 	       
	       'quiver_density',50,... % density of arrows in quiver plot
	       'quiver_scale',1,... % scaling of arrows (0=>no scaling)
	       'mesh','off',...  % show the mesh ?
	       'meshopts',{cell(0)},... % cell array of options passed to gf_plot_slice for the mesh 
	       'deformed_mesh', 'off',... % shows the deformed mesh (only when qdim == mdim)
	       'deformed_meshopts', {cell(0)},... % cell array of options passed to gf_plot_slice for the deformed mesh 
	       'deformation',[],... % plots on the deformed object (only when qdim == mdim)
	       'deformation_mf',[],... % plots on the deformed object (only when qdim == mdim)
	       'deformation_scale','10%',... % indicate the amplitude of the deformation. Can be a percentage of the mesh width if given as a string, or an absolute value if given as a number
	       'cvlst',[],... % list of convexes to plot
	       'title',[],...
	       'contour',[],...
               'mesh_level_set',[], ... % list of contour values
	       'disp_options', 'on');

  opt = getopt(opt,varargin);

  if (ison(opt.disp_options))
    disp(opt);
  end;
  if (ison(opt.zplot))
    if (mdim ~= 2),
      error('zplot allowed only on 2D scalar mesh_fems');
    else opt.interpolated = 'on'; % or patch won't work
    end;
  end;
  is_scalarplot = (ison(opt.norm) + ~isempty(opt.dir));
  if (is_scalarplot > 1),
    error('only one occurence of the options ''norm'' and ''dir'' is allowed');
  end;
  if (ischar(opt.dir)),
    v = zeros(1,qdim);
    if (strcmp(opt.dir,'x'))
      v(1)=1;
    elseif (strcmp(opt.dir,'y'))
      v(2)=1;
    elseif (strcmp(opt.dir,'z'))
      v(3)=1;
    else error('wrong direction');
    end;
    opt.dir=v;
  end;
  scalarplot_dir=opt.dir(:);
  if (qdim == 1) is_scalarplot = 1; scalarplot_dir=1; end;
  if (~isempty(opt.contour) & ~is_scalarplot),
    error('contour plot has no meaning for a vector field');
  end;
  mfdef = mf;
  if (~isempty(opt.deformation_mf)) mfdef = struct(opt.deformation_mf); end;
  dqdim = gf_mesh_fem_get(mfdef,'qdim');
  if (~isempty(opt.deformation) | ison(opt.deformed_mesh)),
    if (mdim ~= dqdim & ~ison(opt.zplot)),
      error(sprintf('can''t plot the deformation of an %dD-object by a %dD-field',mdim,dqdim));
    end;
  end;

  if (isempty(opt.cvlst)), opt.cvlst = gf_mesh_get(mf,'cvid'); end;

  nbdof=gf_mesh_fem_get(mf,'nbdof');
  if (nbdof <= 0), error('invalid finite element mf argument'); end

  if (length(U) ~= nbdof),
    error(['wrong dimensions for U, should' ...
	   ' have ' int2str(nbdof) 'columns']); 
  end

  % save graphical context
  cax = newplot;
  cfig = get(cax,'Parent');
  hold_state = ishold;
  ax_nextplot = lower(get(cax,'NextPlot'));
  fig_nextplot = lower(get(cfig,'NextPlot'));

  if (ison(opt.zplot) | mdim == 3), view(3); else view(2); end;
  
  % build the slice object
  try
    if (~isempty(opt.mesh_level_set)),
      sl = gf_slice({'none'},opt.mesh_level_set,opt.refine,opt.cvlst);
    elseif (gf_mesh_fem_get(mf, 'has_linked_mesh_levelset')),
      sl  = gf_slice({'none'},gf_mesh_fem_get(mf,'linked_mesh_levelset'),opt.refine,opt.cvlst);
    else
      sl  = gf_slice({'none'},gf_mesh_fem_get(mf,'linked mesh'),opt.refine,opt.cvlst);
    end;
  catch
    error(['can''t build slice : ' lasterr]);
  end;
  try
    Usl = gf_compute(mf,U,'interpolate on',sl);
  catch
    error(['can''t interpolate on slice : ' lasterr]);
  end;  
  Psl = gf_slice_get(sl,'pts');

  % plot the original mesh
  if (ison(opt.mesh)), 
    hmesh={gf_plot_slice(sl,'mesh','on',opt.meshopts{:})};
    hold on;
  end;

  % apply the optional deformation
  if (~isempty(opt.deformation) | mfdef.id ~= mf.id),
    ida = gf_mesh_fem_get(mfdef,'linked mesh'); idb = gf_mesh_fem_get(mf,'linked mesh');
    if (ida.id ~= idb.id),
      error('the deformation mesh_fem and the data mesh_fem do not seem to share the same mesh');
    end;
    if (~isempty(opt.deformation)), Udef = opt.deformation;
    else Udef = U;
    end;
    Pdef  = gf_compute(mfdef, Udef, 'interpolate on', sl);
    if (isnumeric(opt.deformation_scale)), dscale = opt.deformation_scale;
    elseif (ischar(opt.deformation_scale) & ...
	    numel(opt.deformation_scale) & opt.deformation_scale(end)=='%'),
      dscale = str2num(opt.deformation_scale(1:end-1));
      mwidth = max([max(Psl,[],2) - min(Psl,[],2)]);
      defmax = max(abs(Pdef(:)));
      if (defmax),	
	dscale = dscale * 0.01 * mwidth / defmax;
      end;
    else error('wrong value for deformation_scale: should be a number or a percentage in a string');
    end;
    Psl = Psl + Pdef*dscale;
    gf_slice_set(sl,'pts', Psl);    
    clear Pdef Udef dscale mwidth defmax;
  end;

  if (is_scalarplot),
    % compute scalar values if necessary
    if (ison(opt.norm)),
      sV = sqrt(sum(Usl.*conj(Usl),1));
    else
      sV = scalarplot_dir(:)'*Usl;
    end;
    % and optionally apply the zplot deformation
    if (ison(opt.pcolor) & ison(opt.zplot) & is_scalarplot),
      Psl = [Psl;sV]; gf_slice_set(sl,'pts',Psl);
    end;
  end;
  
  % plot the deformed mesh
  if (ison(opt.deformed_mesh)),
    hdefmesh={gf_plot_slice(sl,'mesh','on', opt.deformed_meshopts{:})};
    hold on; 
  end;

  if (is_scalarplot),
    % plot the 'surfacic' part
    if (ison(opt.pcolor)),
      gf_plot_slice(sl,'mesh','off','data',sV);
    end;
    
    % basic contour plot (should work also on 3D surfaces)
    contour_colors = [.9 0 0; 0 .8 0; .0 0 .6; .6 .6 0; .7 0 .7; 0 .7 .9]; 
    contour_linestyle = get(cax,'LineStyleOrder');
    hcontour=cell(numel(opt.contour),1);
    for cnum=1:numel(opt.contour),
      c=opt.contour(cnum);
      slC = gf_slice({'isovalues',0,mf,U,c},sl);
      [a,b,c,hcontour{cnum}]=gf_plot_slice(slC,'tube','off','mesh','off');
      if (~isempty(hcontour{cnum})),
        set(hcontour{cnum},...
            'Color',contour_colors(mod(cnum,size(contour_colors,1))+1,:),...
            'LineStyle',contour_linestyle(mod(cnum,numel(contour_linestyle))+1),...
            'LineWidth',1);
      end;
      gf_delete(slC);
    end;
  else
    [a,b,hquiver,c]=gf_plot_slice(sl,'data',Usl,'quiver','on','quiver_density',opt.quiver_density,'quiver_scale',opt.quiver_scale);
  end
  
  if (~hold_state),
    set(cax,'DataAspectRatio', [1 1 1]);
    set(cax,'XLimMode','auto');
    set(cax,'YLimMode','auto');
    set(cax,'ZLimMode','auto');
  end;
  
  % restore graphical context
  set(cax,'NextPlot',ax_nextplot);
  set(cfig,'NextPlot',fig_nextplot);


function r=ison(v)
  r = strcmpi(v,'on');
  
