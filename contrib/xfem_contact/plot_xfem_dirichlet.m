gf_workspace('clear all');
mf = gfMeshFem('load', 'xfem_dirichlet_ls.mf');
lsU = -load('xfem_dirichlet_ls.U')';

nn=3 % choose what is plotted


clf
if nn==0,
  [hsur, hcont] = gf_plot(mf,lsU,'refine',16,'contour',0,'mesh','on', 'pcolor','off');
  set(hcont{1}, 'LineWidth', 2);
  set(hcont{1}, 'Color', 'black');
elseif nn==1
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
  
  
end;

