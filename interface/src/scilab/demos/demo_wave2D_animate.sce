dt = 2*%pi/20;
t  = 0:dt:2*pi-dt/2;

for i=1:length(t),  
  disp(sprintf('theta=%1.3f', t(i)));
  gf_plot(mfu,imag(U(:)'*exp(1*%i*t(i))),'refine',28,'contour',0); 
  //axis([-11 11 -11 11]); caxis([-1 1]);
  xs2png(gcf,sprintf('wave%02d.png',i));
  //F = getframe(gca);
  //mov = addframe(mov,F);
end
//mov = close(mov);

