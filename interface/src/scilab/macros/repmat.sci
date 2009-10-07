// Copyright (C) 2000 Paul Kienzle
//
// This file is part of Octave.
//
// Octave is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// Octave is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Octave; see the file COPYING.  If not, write to the Free
// Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.

// -*- texinfo -*-
// @deftypefn {Function File} {} repmat (@var{A}, @var{m}, @var{n})
// @deftypefnx {Function File} {} repmat (@var{A}, [@var{m} @var{n}])
// Form a block matrix of size @var{m} by @var{n}, with a copy of matrix
// @var{A} as each element.  If @var{n} is not specified, form an 
// @var{m} by @var{m} block matrix.
// @end deftypefn

// Author: Paul Kienzle <pkienzle@kienzle.powernet.co.uk>
// Created: July 2000

function x = repmat (a, m, n)

[nargout, nargin] = argn();

if (nargin < 2 | nargin > 3) then
  error('repmat (a, m, n)');
end

if (nargin == 3) then
  if (~(isscalar (m) & isscalar (n))) then
    error('repmat: with 3 arguments m and n must be scalar');
  end
  idx = [m, n];
else 
  if (isscalar (m)) then
    idx = [m, m];
    n   = m;
  elseif (isvector (m) & length (m) > 1) then
    // Ensure that we have a row vector
    idx = m(:).';
  else
    error('repmat: invalid dimensional argument');
  end
end

x = [];

if (length(a) == 1) then
  if (type(a)==10) then
    x = char (ascii(a) * ones (idx));
  else
    x = a*ones(idx(1),idx(2));
  end
elseif (ndims (a) == 2 & length (idx) < 3) then
  if (type(a)==10)
    x    = char (kron (ones (idx), ascii (a)));
    aidx = size(a);
    x    = a (kron (ones (1, idx(1)), 1:aidx(1)), kron (ones (1, idx(2)), 1:aidx(2)));
  else
    aidx = size(a);
    if (length(aidx) > length(idx)) then
      idx = [idx, ones(1,length(aidx)-length(idx))];
    elseif (length(aidx) < length(idx)) then
      aidx = [aidx, ones(1,length(idx)-length(aidx))];
    end
    cidx = list();
    for i=1:length(aidx)
      cidx(i) = kron (ones (1, idx(i)), 1:aidx(i));
    end
    x = a(cidx(:));
  end
end
endfunction
