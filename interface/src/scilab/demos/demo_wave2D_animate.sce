// Copyright (C) 2010-2020 Yann COLLETTE.
//
// This file is a part of GetFEM
//
// GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 3 of the License,  or
// (at your option) any later version along with the GCC Runtime Library
// Exception either version 3.1 or (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License and GCC Runtime Library Exception for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

lines(0);
stacksize('max');

path = get_absolute_file_path('demo_wave2D_animate.sce');

printf('demo wave2D_animate started\n');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

dt = 2*%pi/20;
t  = 0:dt:2*%pi-dt/2;

h = scf();
h.color_map = jetcolormap(255);

for i=1:length(t),  
  disp(sprintf('theta=%1.3f', t(i)));
  drawlater;
  clf;
  gf_plot(mfu,imag(U(:)'*exp(1*%i*t(i))),'refine',28,'contour',0); 
  h.color_map = jetcolormap(255);
  drawnow;
  
  // use:
  // convert -delay 50 -loop 0 wave*.png animatewave.gif
  // To produce the animated gif image.
  // Convert is an ImageMagick tool.
  xs2png(h.figure_id, path + sprintf('/wave%02d.png',i));
end

printf('demo wave2D_animate terminated\n');
