% Copyright (C) 2013-2020 Franz Chouly.
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
%


function [ u , gu ] = uAnalytic ( x , t )

    % The solution is periodic of order 3
    % Shift the time 't' with t=0 : beginning of the period
    
    tp = rem(t,3);
    
    % The solution has 3 phases
    % Shift the time 'tp' with t=0 : beginning of a phase
    % and get also the phase number
    
    tf = rem(tp,1);
    nf = floor(tp);
    
    % Get the index of the zone in each phase : I,II,III,IV
    % (zones are given according to characteristics of the wave equation)
    
    if (tf<=x)
        if (tf<=(1-x))
            zone = 1;
        else
            zone = 3;
        end
    else
        if (tf<=(1-x))
            zone = 2;
        else
            zone = 4;
        end
    end
 
    % Return the solution according to the Phase (1,2,3) and the 
    % zone (I,II,III,IV)
    
    switch nf
        
        case 0
            
            switch zone
                case 1
                    u = 1/2-x/2;
                    gu = -1/2;
                case 2
                    u = 1/2-tf/2;
                    gu = 0;
                case 3
                    u = 1/2-x/2;
                    gu = -1/2;
                case 4
                    u = 1/2-tf/2;
                    gu = 0;
            end
            
        case 1 
            
            switch zone
                case 1
                    u = -tf/2;
                    gu = 0;
                case 2
                    u = -x/2;
                    gu = -1/2;
                case 3
                    u = -1/2+x/2;
                    gu = 1/2;
                case 4
                    u = -1/2+tf/2;
                    gu = 0;
            end
            
        case 2
            
            switch zone
                case 1
                    u = tf/2;
                    gu = 0;
                case 2
                    u = tf/2;
                    gu = 0;
                case 3
                    u = 1/2-x/2;
                    gu = -1/2;
                case 4
                    u = 1/2-x/2;
                    gu = -1/2;
            end
        
    end
 
end

