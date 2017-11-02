function disp(object)
% disp   Displays some information about the ML_model object
%
% disp(OBJ)
% OBJ is the HR_model object.

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
disp('Morris-Lecar like relaxation oscillators neuron model.');
disp('The Neuron is modeld as follow');
disp('   _');
disp('  |');
disp('  | dV/dt = (-gl*(V - Vl) - gCa*Minf*(V - VCa) - gK*N*(V - Vk ) + I - Isyn)/CM');
disp(' <');
disp('  | dN/dt = lamda_N*(Ninf-N)');
disp('  |_');
disp('   ');
disp('  Minf = 0.5*(1 + tanh((V - V1 )/V2 ));');
disp('  Ninf = 0.5*(1 + tanh((V - V3 )/V4 ));');
disp('  lamda_N = phi*cosh((V - V3 )/(2*V4) );');
disp('    ');
disp('  CM = 5;');
disp('  gK = 8;');
disp('  gl = 2;');
disp('  VCa = 120;');
disp('  Vk = -80;');
disp('  Vl = -60;');
disp('  V1 = -1.2;');
disp('  V2 = 18;');
disp(' ');
disp('Parameters value:');
disp(['     gCa = ',num2str(object.gCa)]);
disp(['     V3 = ',num2str(object.V3)]);
disp(['     V4 = ',num2str(object.V4)]);
disp(['     phi = ',num2str(object.phi)]);
disp(['     I = ',num2str(object.I)]);
end