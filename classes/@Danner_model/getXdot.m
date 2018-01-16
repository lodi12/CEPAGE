function x_dot = getXdot(object,t,x,varargin)
% getXdot    Computes the derivative of the state
%
%  x_dot = getXdot(object,t,x)
%   compute the time derivative of the model at time instant t, in state x
%   and with I_{ext} = 0;
%
%  x_dot = getXdot(object,t,x,I_{ext})
%   compute the time derivative of the model at time instant t, in state x
%   and with I_{ext} specified by the user;

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

% synapsis input
Isyn = 0;
if nargin == 4
    if isscalar(varargin{end})
        Isyn = varargin{end};
    end
end


% Define param
C = object.C;

gNaP = object.gNaP;
ENa = object.ENa;

gl = object.gl;
El = object.El;
VhalfM = object.VhalfM;
km = object.km;
Vhalfh = object.Vhalfh;
kh = object.kh;
tau0 = object.tau0;
tauMax = object.tauMax;
VhalfTau = object.VhalfTau;
kTau = object.kTau;
Di = object.Di;
gSynE = object.gSynE;
EsynE = object.EsynE;
I = object.I;
V = x(1);
h = x(2);

m = 1/(1+exp((V-VhalfM)/km));

tauH = tau0 +((tauMax-tau0)/cosh((V-VhalfTau)/kTau));

hInf = 1/(1+exp((V-Vhalfh)/kh));
dh = 1/tauH*(hInf-h);


Idrive = gSynE*(V-EsynE)*Di;

Vdot = 1/C*(-gl*(V-El)-gNaP*m*h*(V-ENa)-Idrive+Isyn+I);

x_dot = [Vdot;dh];


end