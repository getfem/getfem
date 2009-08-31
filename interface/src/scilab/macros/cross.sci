// Copyright (C) 1995, 1996, 1997, 1999, 2000, 2002, 2004, 2005, 2006,
//               2007 Kurt Hornik
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
// @deftypefn {Function File} {} cross (@var{x}, @var{y}, @var{dim})
// Computes the vector cross product of the two 3-dimensional vectors
// @var{x} and @var{y}.
//
// @example
// @group
// cross ([1,1,0], [0,1,1])
//      @result{} [ 1; -1; 1 ]
// @end group
// @end example
//
// If @var{x} and @var{y} are matrices, the cross product is applied 
// along the first dimension with 3 elements. The optional argument 
// @var{dim} is used to force the cross product to be calculated along
// the dimension defined by @var{dim}.
// @end deftypefn

// Author: Kurt Hornik <Kurt.Hornik@wu-wien.ac.at>
// Created: 15 October 1994
// Adapted-By: jwe

function z = cross (x, y, dim)
	
[nargout,nargin] = argn();

if (ndims (x) < 3 & ndims (y) < 3 & nargin < 3) then
  // COMPATIBILITY -- opposite behaviour for cross(row,col)
  // Swap x and y in the assignments below to get the matlab behaviour.
  // Better yet, fix the calling code so that it uses conformant vectors.
  if (size(x,2) == 1 & size(y,1) == 1) then
    warning('cross: taking cross product of column by row');
    y = y.';
  elseif (size(x,1) == 1 & size(y,2) == 1) then
    warning('cross: taking cross product of row by column');
    x = x.';
  end
end

if (nargin == 2) then
   dim = find (size (x) == 3, 1);
   if (isempty (dim)) then 
     error('cross: must have at least one dimension with 3 elements');
   end
 else
   if (size (x) ~= 3) then
     error('cross: dimension dim must have 3 elements');
   end
end

nd = ndims (x);
sz = size (x);
idx1 = list();
for i = 1:nd
  idx1(i) = 1:sz(i);
end
idx2 = idx1;
idx3 = idx1;
idx1(dim) = 1;
idx2(dim) = 2;
idx3(dim) = 3;

if (and(size(x)==size(y))) then
  z = cat(dim, ...
          (x(idx2(:)) .* y(idx3(:)) - x(idx3(:)) .* y(idx2(:))), ...
          (x(idx3(:)) .* y(idx1(:)) - x(idx1(:)) .* y(idx3(:))), ...
          (x(idx1(:)) .* y(idx2(:)) - x(idx2(:)) .* y(idx1(:))));
else
  error('cross: x and y must have the same dimensions');
end
endfunction
