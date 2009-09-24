// Copyright (C) 2000, 2005, 2006, 2007 Paul Kienzle
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
// @deftypefn {Function File} {} setdiff (@var{a}, @var{b})
// @deftypefnx {Function File} {} setdiff (@var{a}, @var{b}, "rows")
// Return the elements in @var{a} that are not in @var{b}, sorted in
// ascending order.  If @var{a} and @var{b} are both column vectors
// return a column vector, otherwise return a row vector.
//
// Given the optional third argument @samp{"rows"}, return the rows in
// @var{a} that are not in @var{b}, sorted in ascending order by rows.
// @seealso{unique, union, intersect, setxor, ismember}
// @end deftypefn

// Author: Paul Kienzle
// Adapted-by: jwe

function c = _setdiff (a, b, byrows_arg)

[nargout,nargin] = argn();

if (nargin < 2 | nargin > 3) then
  error('setdiff: 2 or 3 arguments allowed');
end

byrows = %F;

if (nargin == 3) then
  if (byrows_arg~="rows") then
    error('expecting third argument to be ''rows''');
  elseif (typeof(a)=='list' | typeof(b)=='list') then
    warning('setdiff: ''rows'' not valid for cell arrays');
  else
    byrows = %T;
  end
end

if (byrows) then
  c = unique (a, 'r');
  if (~isempty (c) & ~isempty (b)) then
    // Form a and b into combined set.
    b = unique (b, 'r');
    [dummy, idx] = gsort([c; b],'lr','i');
    // Eliminate those elements of a that are the same as in b.
    dups = find (and(dummy(1:$-1,:) == dummy(2:$,:),2));
    c = [c;b];
    dummy(unique([dups dups+1]),:) = [];
    c = dummy;
  end
else
  c = unique(a);
  if (~isempty (c) & ~isempty (b)) then
    // Form a and b into combined set.
    b = unique(b);
    // Doesn't work with string
    if (typeof([c(:); b(:)])=='string') then
      [dummy, idx] = gsort(-[ascii(c(:)); ascii(b(:))]);
    else
      [dummy, idx] = gsort(-double([c(:); b(:)]));
    end
    
    // Eliminate those elements of a that are the same as in b.
    dups = find (dummy(1:$-1)==dummy(2:$)); 
    c(idx(dups)) = [];
    // Reshape if necessary.
    if (size (c, 1) ~= 1 & size (b, 1) == 1) then
      c = c.';
    end
  end
end  
endfunction
  
//!assert(setdiff(["bb";"zz";"bb";"zz"],["bb";"cc";"bb"],"rows"), "zz") // OK ??
//!assert(setdiff(["b";"z";"b";"z"],["b";"c";"b"],"rows"), "z") // OK ??
//!assert(setdiff(["b";"z";"b";"z"],["b";"c";"b"]), "z") // NOK
//!assert(setdiff([1, 1; 2, 2; 3, 3; 4, 4], [1, 1; 2, 2; 4, 4], "rows"), [3 3]) // OK
//!assert(setdiff([1; 2; 3; 4], [1; 2; 4], "rows"), 3) // OK
//!assert(setdiff([1, 2; 3, 4], [1, 2; 3, 6], "rows"), [3, 4]) // OK ??
//!assert(setdiff(list("one","two";"three","four"),list("one","two";"three","six")), list("four")) // NOK
