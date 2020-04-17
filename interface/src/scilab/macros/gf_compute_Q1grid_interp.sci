// Copyright (C) 2009-2020 Yann Colette
// 
//  This file is a part of GetFEM++
// 
//  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
//  under  the  terms  of the  GNU  Lesser General Public License as published
//  by  the  Free Software Foundation;  either version 3 of the License,  or
//  (at your option) any later version along with the GCC Runtime Library
//  Exception either version 3.1 or (at your option) any later version.
//  This program  is  distributed  in  the  hope  that it will be useful,  but
//  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
//  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
//  License and GCC Runtime Library Exception for more details.
//  You  should  have received a copy of the GNU Lesser General Public License
//  along  with  this program;  if not, write to the Free Software Foundation,
//  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

function [U2,Iq,MF2]=gf_compute_Q1grid_interp(MF1,U1,varargin)
// See help on gf_compute
//  $Id: gf_compute_Q1grid_interp.m 1937 2005-03-08 16:50:13Z pommier $

[nargout,nargin] = argn();

if (nargin < 3) then
  error('not enough input arguments');
end;
  
gf_workspace('push', 'gf_compute_Q1grid_interp');

meshpts = gf_mesh_get(MF1, 'pts');
zmin = min(meshpts,'c');
zmax = max(meshpts,'c');
ndim = length(zmin);
if (ndim > 3) then
  error('this number of dim is not supported (patch me)'); 
end;

X = list();

try
  select varargin(1)
    case 'regular h' then
      if (length(varargin) ~= 2) then error('wrong number of arguments'); end;
      if (length(varargin(2)) ~= ndim) then error('invalid dimension'); end;
      for i=1:ndim
        if (varargin(2)(i) <= 0) then error('invalid step value'); end;
        X(i) = zmin(i):(varargin(2)(i)):zmax(i);
      end
    case 'regular N' then
      if (length(varargin) ~= 2) then error('wrong number of arguments'); end;
      if (length(varargin(2)) ~= ndim) then error('invalid dimension'); end;
      for i=1:ndim,
        if (varargin(2)(i) <= 0) then error('invalid number of cells'); end;
        h = (zmax(i) - zmin(i))/(varargin(2)(i));
        X(i) = zmin(i):h:zmax(i);
      end
   else
    X = varargin(1);
    if (typeof(X)~='list') then error('grid points should be stored in a list array of size nbdim'); end;
    if (length(X) ~= ndim) then error('wrong number of dimension in the grid points argument'); end;
  end

  Q=gf_mesh_fem_get(MF1,'qdim');
  M = gf_mesh('cartesian', X(:));
  MF2 = gf_mesh_fem(M,Q);
  gf_mesh_fem_set(MF2, 'classical fem', 1); // Q1 fem
  mfU2 = gf_compute(MF1,U1, 'interpolate on', MF2);

  PTS = gf_mesh_fem_get(MF2, 'basic dof nodes');

  PTS = PTS($:-1:1,1:Q:$);   // (x,y,z)->(z,y,x) and remove duplicate dof
  [PTS,I] = gsort(PTS','lr','i'); // sort points, by z then by y then by x etc..
  I = Q*(I-1) + 1;
  sz = Q;
  for i=1:length(X) sz = [sz length(X(i))]; end;
  Iq=zeros(Q,length(I));
  for q=1:Q,
    Iq(q,:) = I'+(q-1);
  end;
  Iq = Iq(:);
  U2 = matrix(mfU2(Iq),sz);
  if (nargout == 3) then
    gf_workspace('keep', MF2);
  end;
catch
  gf_workspace('pop');
  error(lasterror());
end
gf_workspace('pop');
endfunction

