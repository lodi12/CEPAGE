function J = getJacobian(object,t,x)

% getJacobian    Computes the Jacobian of the vector field in point x 
%
%  J = getJacobian(object,x)
%   compute theJacobian of the model in state x
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

nx = object.nx;

if size(x,1) ~= nx
    error('input vector must be a nx x npoints vector');
end

np = size(x,2);

J = cell(np,1);

gK = 8;

VCa = 120;
Vk = -80;

V1 = -1.2;
V2 = 18;

CM = object.CM;
I = object.I;
gCa = object.gCa;
phi = object.phi;
V3 = object.V3;
V4 = object.V4;
gL = object.gL;
Vl = object.EL;
for i=1:np

xx = x(:,i);
    
V = xx(1);
N = xx(2);

% J1 = -(gL + N*gK + (gCa*(tanh((V - V1)/V2) + 1))/2 - (gCa*(V - VCa)*(tanh((V - V1)/V2)^2 - 1))/(2*V2))/CM;
J1 = 1/CM*(-gL-1/V2*0.5*gCa*(V-VCa)/(cosh((V-V1)/V2).^2) -0.5*gCa*(1+tanh((V-V1)/V2)) -gK*N);
J2 = -(gK*(V - Vk))/CM;
J3 = 0.5*phi/2/V4*(sinh((V-V3)/2/V4)+sinh((V-V3)/2/V4)*tanh((V-V3)/V4)+2*cosh((V-V3)/2/V4)./(cosh((V-V3)/V4).^2))-phi*N/2/V4*sinh((V-V3)/2/V4);
% J3 = (phi*cosh((V - V3)/(2*V4))*(tanh((V - V3)/V)^2 - 1)*((V - V3)/V^2 - 1/V))/2 + (phi*sinh((V - V3)/(2*V4))*(tanh((V - V3)/V)/2 - N + 1/2))/(2*V4);
J4 = -phi*cosh((V - V3)/(2*V4));

J{i} = [J1 J2; J3 J4];

if size(x,2) == 1
    JJ = J{1};
    J = [];
    J = JJ;
end

end
