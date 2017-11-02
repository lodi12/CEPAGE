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

eps = object.eps;
I = object.I;

for i=1:np

xx = x(:,i);
    
J1 = 1-3*xx(1)^2;
J2 = -1;
J3 = eps*10*exp(-10*xx(1))/((1+exp(-10*xx(1)))^2);
J4 = -eps;

J{i} = [J1 J2; J3 J4];

if size(x,2) == 1
    JJ = J{1};
    J = [];
    J = JJ;
end

end