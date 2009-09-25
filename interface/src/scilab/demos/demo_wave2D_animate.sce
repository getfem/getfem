dt = 2*%pi/20;
t  = 0:dt:2*%pi-dt/2;

h = scf();;
h.color_map = jetcolormap(255);

for i=1:length(t),  
  disp(sprintf('theta=%1.3f', t(i)));
  drawlater;
  gf_plot(mfu,imag(U(:)'*exp(1*%i*t(i))),'refine',28,'contour',0); 
  h.color_map = jetcolormap(255);
  drawnow;
  
  // use:
  // convert -delay 50 -loop 0 wave*.png animatewave.gif
  // To produce the animated gif image.
  // Convert is an ImageMagick tool.
  xs2png(h.figure_id,sprintf('wave%02d.png',i));
  clf;
end

