%  Copyright (C) 1999-2020 Yves Renard
%
%  This file is a part of GetFEM++
%
%  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
%  under  the  terms  of the  GNU  Lesser General Public License as published
%  by  the  Free Software Foundation;  either version 3 of the License,  or
%  (at your option) any later version along with the GCC Runtime Library
%  Exception either version 3.1 or (at your option) any later version.
%  This program  is  distributed  in  the  hope  that it will be useful,  but
%  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
%  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
%  License and GCC Runtime Library Exception for more details.
%  You  should  have received a copy of the GNU Lesser General Public License
%  along  with  this program;  if not, write to the Free Software Foundation,
%  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

function o=getopt(opt,v)
  o = opt;
  if (nargin==1) return; end;
  if (mod(length(v),2) ~= 0) error('odd number of property/value pairs'); end;
  for i=1:2:length(v),
    optname=v{i};
    optval =v{i+1};
    if (~ischar(optname)) error(['expecting a property name, found a ' class(optname)]); end;
    if (~isfield(opt,optname)) error(['unknown property :"' optname '"']); end;
    o = setfield(o, optname, optval);
  end;
