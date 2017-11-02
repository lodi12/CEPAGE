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

N = object.getN;

g_in = object.get_g_in;
g_ex = object.get_g_ex;
g_el = object.get_g_el; 

tetaSyn = object.tetaSyn;
uSyn = object.uSyn;

Esyn_Ex = object.EsynEx;
Esyn_In = object.EsynIn;


ndim_vector = zeros(N,1);

for i=1:N
    ndim_vector(i) = object.neurons{i}.getnx;
end

ndim = sum(ndim_vector);

if size(x,1) ~= ndim
    error('input vector must be a nx x npoints vector');
end

np = size(x,2);

J = cell(np,1);


for i=1:np
    tmpJ = zeros(ndim,ndim);
    xx = x(:,i);
    currentPointer = 1;
    for k=1:N
        neurNx = ndim_vector(k);
        xsingolo = xx(currentPointer:currentPointer+neurNx-1);
        neuronJac = object.neurons{k}.getJacobian(0,xsingolo);
        
        %uncoupled neuron
        tmpJ(currentPointer:currentPointer+neurNx-1,currentPointer:currentPointer+neurNx-1) = neuronJac;
       
        % Add synaptic contribute
        % Add dVidot/dVi
        otherI = 1;
        for j=1:N
            tmpJ(currentPointer,currentPointer) =  tmpJ(currentPointer,currentPointer) - ...
                g_in(k,j)/(1+exp(-uSyn*(xx(otherI)-tetaSyn))) - ...
                g_ex(k,j)/(1+exp(-uSyn*(xx(otherI)-tetaSyn))) - g_el(k,j);
            otherI = otherI + ndim_vector(j);
        end
        
        % Add dVidot/dVj
        otherI = 1;
        for j=1:N
            tmpJ(currentPointer,otherI) =  tmpJ(currentPointer,otherI) + ...
                g_in(k,j)*(Esyn_In-xx(currentPointer))*(uSyn*exp(-uSyn*(xx(otherI)-tetaSyn)))/(1+exp(-uSyn*(xx(otherI)-tetaSyn)))^2 + ...
                g_ex(k,j)*(Esyn_Ex-xx(currentPointer))*(uSyn*exp(-uSyn*(xx(otherI)-tetaSyn)))/(1+exp(-uSyn*(xx(otherI)-tetaSyn)))^2 + ...
                g_el(k,j);
            otherI = otherI + ndim_vector(j);
        end
        
        
        currentPointer = currentPointer+ndim_vector(k);
        
    end
    
    J{i} = tmpJ;
    
end

if size(x,2) == 1
    JJ = J{1};
    J = [];
    J = JJ;
end
