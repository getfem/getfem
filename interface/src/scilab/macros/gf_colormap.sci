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

function out = gf_colormap(name),
// function  c=gf_colormap(name)
//   return a colormap, or change the current colormap.
//   name can be: 'tripod', 'chouette', 'froid', 'tank'
//   or 'earth'.

[nargout,nargin] = argn();

if (name=='tripod') then
  r  = [0.7 0.7 0.7]; 
  l  = r($,:); 
  s  = 63; 
  s1 = 20;
  s2 = 25;
  s3 = 48;
  s4 = 55; 
  for i=1:s, 
    c1 = max(min((i-s1)/(s2-s1),1),0);
    c2 = max(min((i-s3)/(s4-s3),1),0); 
    r($+1,:)=(1-c2)*((1-c1)*l + c1*[1 0 0]) + c2*[1 .8 .2]; 
  end
elseif (name=='chouette') then
  gg = [0.8 1.0  0.8; 
        0.7 0.9  0.4;
        0.3 0.8  0.2;
        0.1 0.7  0.4;
        0.2 0.7  1.0;
        0.3 0.3  1.0;
        1.0 0.8  0.1;
        1.0 0.6  0.1;
        1.0 0.45 0.1;
        1.0 0.3  0.1];
  r = matrix(repmat(gg',6,1),3,60)'; 
elseif (name=='froid') then
  gg = [0.8 1.0 0.8; 
        0.7 0.9 0.4;
        0.3 0.8 0.2;
        0.1 0.7 0.4;
        0.2 0.7 1.0;
        0.3 0.3 1.0];
  r = matrix(repmat(gg',10,1),3,60)';     
elseif (name=='tank') then
  r=[0.0 0.0 1.0; 
     0.0 0.5 1.0; 
     0.0 1.0 0.5; 
     0.0 1.0 0.0; 
     0.5 1.0 0.0; 
     1.0 0.5 0.0; 
     1.0 0.4 0.0; 
     1.0 0.0 0.0;
     1.0 0.2 0.0; 
     1.0 0.4 0.0;
     1.0 0.6 0.0;
     1.0 0.8 0.0];
    r = matrix(repmat(r',5,1),3,60)';     
elseif (name=='earth') then
  r=[252 233  79; //	Butter 1
     247 222  30;
     237 212   0; //	Butter 2
     216 180   0;
     196 160   0; //	Butter 3
     138 226  52; //	Chameleon 1
     115 210  22; //	Chameleon 2
      78 154   6];
  r = matrix(r'/255, 8,1);
  r = matrix(r,3,length(r)/3)';     
elseif (name=='getfem') then
  r = [252 233  79; //	Butter 1
       237 212   0; //	Butter 2
       196 160   0; //	Butter 3
       138 226  52; //	Chameleon 1
       115 210  22; //	Chameleon 2
        78 154   6; //	Chameleon 3
       252 175  62; //	Orange 1
       245 121   0; //	Orange 2
       206  92   0; //	Orange 3
       114 159 207; //	Sky Blue 1
       114 159 207; //	Sky Blue 1
        52 101 164; //	Sky Blue 2
        52 101 164; //	Sky Blue 2
        32  74 135; //	Sky Blue 3
        32  74 135; //	Sky Blue 3
       173 127 168; //	Plum 1
       173 127 168; //	Plum 1
       173 127 168; //	Plum 1
       117  80 123; //	Plum 2
       117  80 123; //	Plum 2
       117  80 123; //	Plum 2
        92  53 102; //	Plum 3
        92  53 102; //	Plum 3
        92  53 102; //	Plum 3
       233 185 110; //	Chocolate 1
       233 185 110; //	Chocolate 1
       233 185 110; //	Chocolate 1
       233 185 110; //	Chocolate 1
       193 125  17; //	Chocolate 2
       193 125  17; //	Chocolate 2
       193 125  17; //	Chocolate 2
       193 125  17; //	Chocolate 2
       143  89   2; //	Chocolate 3
       143  89   2; //	Chocolate 3
       143  89   2; //	Chocolate 3
       143  89   2; //	Chocolate 3
       239  41  41; //	Scarlet Red 1
       239  41  41; //	Scarlet Red 1
       239  41  41; //	Scarlet Red 1
       239  41  41; //	Scarlet Red 1
       239  41  41; //	Scarlet Red 1
       204   0   0; //	Scarlet Red 2
       204   0   0; //	Scarlet Red 2
       204   0   0; //	Scarlet Red 2
       204   0   0; //	Scarlet Red 2
       204   0   0; //	Scarlet Red 2
       164   0   0; //	Scarlet Red 3
       164   0   0; //	Scarlet Red 3
       164   0   0; //	Scarlet Red 3
       164   0   0; //	Scarlet Red 3
       164   0   0]; //	Scarlet Red 3

  r = r/255;
else
  error('wrong colormap');
end

if (nargout) then
  out = r;
else
  f = gcf();
  f.color_map = r;
end
endfunction

