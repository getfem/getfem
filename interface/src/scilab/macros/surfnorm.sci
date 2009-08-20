// Copyright (C) 2007 David Bateman
//
// This file is part of Octave.
//
// Octave is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// Octave is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Octave; see the file COPYING.  If not, see
// <http://www.gnu.org/licenses/>.

// -*- texinfo -*-
// @deftypefn {Function File} {} surfnorm (@var{x}, @var{y}, @var{z})
// @deftypefnx {Function File} {} surfnorm (@var{z})
// @deftypefnx {Function File} {[@var{nx}, @var{ny}, @var{nz}] =} surfnorm (@dots{})
// @deftypefnx {Function File} {} surfnorm (@var{h}, @dots{})
// Find the vectors normal to a meshgridded surface. The meshed gridded 
// surface is defined by @var{x}, @var{y}, and @var{z}. If @var{x} and 
// @var{y} are not defined, then it is assumed that they are given by
//
// @example
// [@var{x}, @var{y}] = meshgrid (1:size(@var{z}, 1), 
//                      1:size(@var{z}, 2));
// @end example
//
// If no return arguments are requested, a surface plot with the normal 
// vectors to the surface is plotted. Otherwise the componets of the normal
// vectors at the mesh gridded points are returned in @var{nx}, @var{ny},
// and @var{nz}.
//
// The normal vectors are calculated by taking the cross product of the 
// diagonals of eash of teh quadrilaterals in the meshgrid to find the 
// normal vectors of the centers of these quadrilaterals. The four nearest
// normal vectors to the meshgrid points are then averaged to obtain the 
// normal to the surface at the meshgridded points.
//
// An example of the use of @code{surfnorm} is
//
// @example
// surfnorm (peaks (25));
// @end example
// @seealso{surf, quiver3}
// @end deftypefn

function [Nx, Ny, Nz] = surfnorm (varargin)

[nargout,nargin] = argn();

if (nargin == 1) then
  z = varargin(1);
  [x, y] = meshgrid (1:size(z,1), 1:size(z,2));
  ioff = 2;
else
  x = varargin(1);
  y = varargin(2);
  z = varargin(3);
  ioff = 4;
end

// Make life easier, and avoid having to do the extrapolation later, do
// a simpler linear extrapolation here. This is approximative, and works
// badly for closed surfaces like spheres.
xx = [2 .* x(:,1) - x(:,2), x, 2 .* x(:,$) - x(:,$-1)];
xx = [2 .* xx(1,:) - xx(2,:); xx; 2 .* xx($,:) - xx($-1,:)];
yy = [2 .* y(:,1) - y(:,2), y, 2 .* y(:,$) - y(:,$-1)];
yy = [2 .* yy(1,:) - yy(2,:); yy; 2 .* yy($,:) - yy($-1,:)];
zz = [2 .* z(:,1) - z(:,2), z, 2 .* z(:,$) - z(:,$-1)];
zz = [2 .* zz(1,:) - zz(2,:); zz; 2 .* zz($,:) - zz($-1,:)];

u_x = xx(1:$-1,1:$-1) - xx(2:$,2:$);
u_y = yy(1:$-1,1:$-1) - yy(2:$,2:$);
u_z = zz(1:$-1,1:$-1) - zz(2:$,2:$);
v_x = xx(1:$-1,2:$) - xx(2:$,1:$-1);
v_y = yy(1:$-1,2:$) - yy(2:$,1:$-1);
v_z = zz(1:$-1,2:$) - zz(2:$,1:$-1);

c = cross ([u_x(:), u_y(:), u_z(:)], [v_x(:), v_y(:), v_z(:)]);
w_x = matrix(c(:,1), size(u_x));
w_y = matrix(c(:,2), size(u_y));
w_z = matrix(c(:,3), size(u_z));

// Create normal vectors as mesh vectices from normals at mesh centers
Nx = (w_x(1:$-1,1:$-1) + w_x(1:$-1,2:$) + w_x(2:$,1:$-1) + w_x(2:$,2:$)) ./ 4; 
Ny = (w_y(1:$-1,1:$-1) + w_y(1:$-1,2:$) + w_y(2:$,1:$-1) + w_y(2:$,2:$)) ./ 4; 
Nz = (w_z(1:$-1,1:$-1) + w_z(1:$-1,2:$) + w_z(2:$,1:$-1) + w_z(2:$,2:$)) ./ 4; 

// Normalize the normal vectors
len = sqrt (Nx.^2 + Ny.^2 + Nz.^2);
Nx  = Nx ./ len;
Ny  = Ny ./ len;
Nz  = Nz ./ len;
endfunction
