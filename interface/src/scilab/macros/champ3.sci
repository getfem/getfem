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

function champ3(x,y,z,fx,fy,fz,c)

X = zeros(2*length(x),1);
Y = zeros(2*length(y),1);
Z = zeros(2*length(z),1);

X(1:2:$) = matrix(x,length(x),1);
X(2:2:$) = matrix(x+fx,length(x),1);
Y(1:2:$) = matrix(y,length(y),1);
Y(2:2:$) = matrix(y+fy,length(y),1);
Z(1:2:$) = matrix(z,length(x),1);
Z(2:2:$) = matrix(z+fz,length(z),1);

xsegs(X,Y,c);
e = gce();
e.arrow_size = 1;
e.data(:,3) = Z;
a = gca();
a.view = '3d';
a.data_bounds(:,3) = [min(e.data(:,3)); max(e.data(:,3))];
endfunction

