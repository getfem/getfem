gf_workspace('clear all');
mf = gfMeshFem('load', 'xfem_dirichlet_ls.mf');
lsU = -load('xfem_dirichlet_ls.U')';

nn=4 % choose what is plotted


clf
if nn==0,
  disp('plot the exported mf');
  [hsur, hcont] = gf_plot(mf,lsU,'refine',16,'contour',0,'mesh','on', 'pcolor','off');
  set(hcont{1}, 'LineWidth', 2);
  set(hcont{1}, 'Color', 'black');
elseif nn==1
  disp('plot the cut mesh');
  mc = gfMesh('load','cut.mesh');
  mfc = gfMeshFem(mc); 
  set(mfc,'classical_fem',2);
  
  lsUc = gf_compute(mf, lsU, 'interpolate on', mfc);
  
  %[hsur, hcont] = gf_plot(mf, lsU, 'refine', 24, 'zplot', 'on');
  %hold on;
  [hsur, hcont] = gf_plot(mfc, lsUc.*(lsUc>0), 'refine', 8, 'zplot', 'on', 'mesh','on', 'pcolor','on'); hold on;
  colormap([.8 1 .8]);
  %[hsur, hcont] = gf_plot(mf,lsU,'refine',4,'contour',0,'pcolor','off');
  
  %set(hcont{1}, 'LineWidth', 2);
  %set(hcont{1}, 'Color', 'black');
  %axis('tight'); axis off;
elseif nn==2 || nn==3,
  disp('plot the solution, with the 0 isovalue');
  sl=gfSlice('load','xfem_dirichlet.sl');
  slU=load('xfem_dirichlet.slU')';
  P=gf_slice_get(sl,'pts'); P=[P(1:2,:);slU];
  gf_slice_set(sl,'pts',P);
  gf_plot_slice(sl, 'data', slU, 'mesh','on','mesh_edges','on');
  
  

  m=get(mf, 'linked_mesh');
  slc=gfSlice({'isovalues', 0, mf, lsU, 0}, m, 16);
  hold on;
  P2=gf_slice_get(slc, 'pts');
  gf_slice_set(slc, 'pts', [P2;0.1 * ones(1,size(P2,2))]);
  [h1,h2,h3,h4]=gf_plot_slice(slc, 'tube','off','mesh_slice_edges_color','black');
  set(h4, 'LineWidth', 2);
  
  if (nn == 2),
    %set(hcont{1}, 'Color', 'black');
    view(3); camlight; view(2);
  else
    slc2=gfSlice('load', 'xfem_dirichlet.sl0');
    hold on;
    [h1,h2,h3,h4]=gf_plot_slice(slc2, 'tube','off','mesh_slice_edges_color','white');
    set(h4, 'LineWidth', 4);
    caxis([-.2 .3]); gf_colormap('froid');
  end;
elseif nn==4,
  disp('plotting the lagrange multipliers on the dirichlet boundary');
  sll=gfSlice('load','xfem_dirichlet.sll');
  slL=load('xfem_dirichlet.slL')';
  P0=gf_slice_get(sll, 'pts');
  [h1,h2,h3,h4]=gf_plot_slice(sll, 'tube','off','mesh_slice_edges_color',[.3 .3 .3]);
  hold on;
  gf_slice_set(sll,'pts',[P0 ; slL*0.05]);
  [hh1,hh2,hh3,hh4]=gf_plot_slice(sll, 'tube','off','mesh_slice_edges_color','black','mesh_slice_edges_width',1.5);
  sl=gfSlice('load','xfem_dirichlet.sl');
  gf_plot_slice(sl,'mesh','on');
  
  npt = size(P0, 2);
  P0 = [P0;zeros(1,npt)];
  P1 = gf_slice_get(sll,'pts'); 
  lseg = gf_slice_get(sll,'splxs', 1);
  F=[lseg(1,:) lseg(2,:); lseg(2,:) npt+lseg(2,:); npt+lseg(1,:) npt+lseg(1,:)];
  %F=[lseg; npt+lseg(2,:)];
  h=patch('Vertices',[P0 P1]', 'Faces', F');
  hold off;
  set(h,'FaceAlpha',0.3);
  set(h,'LineStyle','none');
  set(gcf,'renderer','opengl');
  set(gcf,'color','white');
  axis off;
  view(3);
  camzoom(1.7);
  print(gcf,'-dpng','-r300', 'lagrange_multipliers.png');
end;

