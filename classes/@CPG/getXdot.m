function x_dot = getXdot(object,t,x)
% getXdot    Computes the derivative of the state
%
%  x_dot = getXdot(object,t,x)
%   compute the time derivative of the model at time instant t, in state x
%   and with I_{ext} = 0;
%
% Contributors:
%
% Matteo Lodi (matteo.lodi@edu.unige.it)
%
% Copyright (C) 2016 University of Genoa, Italy.

% Legal note:
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc.,
% 59 Temple Place, Suite 330,
% Boston, MA  02111-1307  USA



N = object.N;
g_in = object.g_in;
g_ex = object.g_ex;
g_el = object.g_el;

inhActivation = object.inhActivation;
excActivation = object.excActivation;

neurons = object.neurons;

%The neurons must be all the same name of state variable
nx = neurons{1}.getnx();

x_dot = zeros(nx*N,1);

%Valori sinapsi
EsynIn = object.EsynIn;
EsynEx = object.EsynEx;



for i=1:N
    Isyn = 0;
    beginI = (i-1)*nx+1;
    endI = (i-1)*nx+nx;
    %Compute exc syn
    for j=1:N
        if g_ex(i,j) ~= 0
            otherI = (j-1)*nx+1;
            Isyn = Isyn + g_ex(i,j)*(-x(beginI)+EsynEx)*excActivation{i,j}.getActivation(x(otherI));
        end
    end
    %compute inh syn
    for j=1:N
        if g_in(i,j) ~= 0
            otherI = (j-1)*nx+1;
            Isyn = Isyn + g_in(i,j)*(-x(beginI)+EsynIn)*inhActivation{i,j}.getActivation(x(otherI));
        end
    end
    %compute el syn
    for j=1:N
        if g_el(i,j) ~= 0
            otherI = (j-1)*nx+1;
            Isyn = Isyn + g_el(i,j)*(+x(otherI)-x(beginI));
        end
    end
    
    x_dot(beginI:endI) = neurons{i}.getXdot(t,x(beginI:endI),Isyn);
end

end
