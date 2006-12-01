function icareplot(nn)
  global mfu U mfdu DU rot
  gf_workspace('clear all');
  mfu=gfMeshFem('load','icare.mf_u');
  mfdu=gfMeshFem(get(mfu,'linked mesh')); 
  set(mfdu,'fem',gfFem('FEM_PK_DISCONTINUOUS(2,1)'));

  for n=nn,
    U=load(sprintf('icare.U%d',n))';
    DU=gf_compute(mfu,U,'gradient',mfdu);
    rot=DU(1,2,:)-DU(2,1,:); 
  
    subplot(2,1,1); gf_plot(mfu,U,'refine',2,'norm','on'); 
    %caxis([0 1.5]); 
    colorbar;
    subplot(2,1,2); gf_plot(mfdu,rot(:)','refine',1); 
    %caxis([-2 2]); 
    colorbar;
    disp('press any key'); pause;
  end;