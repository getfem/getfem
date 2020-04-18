% Copyright (C) 2004-2020 Julien Pommier.
%
% This file is a part of GetFEM
%
% GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
% under  the  terms  of the  GNU  Lesser General Public License as published
% by  the  Free Software Foundation;  either version 3 of the License,  or
% (at your option) any later version along with the GCC Runtime Library
% Exception either version 3.1 or (at your option) any later version.
% This program  is  distributed  in  the  hope  that it will be useful,  but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
% License and GCC Runtime Library Exception for more details.
% You  should  have received a copy of the GNU Lesser General Public License
% along  with  this program;  if not, write to the Free Software Foundation,
% Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

function R=axrot_matrix(A, B, theta)
  n=(B-A); n = n/norm(n);
  a=n(1); b=n(2); c=n(3); 
  d=sqrt(b^2+c^2);
  T=eye(4); T(1:3,4)=-A(:);
  Rx=eye(4); 
  if (norm(n(2:3))>1e-6)
    Rx(2:3,2:3)=[c -b; b c]/d;
  end;
  Ry=eye(4); Ry([1 3],[1 3])=[d -a; a d];
  Rz=eye(4); Rz(1:2,1:2)=[cos(theta) sin(theta); -sin(theta) cos(theta)];
  R = inv(T)*inv(Rx)*inv(Ry)*Rz*Ry*Rx*T;
  
