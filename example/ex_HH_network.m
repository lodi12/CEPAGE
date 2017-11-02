%__________________________________________________________________________
%                                                                         
% ex_HH_Network.m
%__________________________________________________________________________
%
% In this example a network composed of Hodgkin-Huxley neurons is analyzed 
% in time and phase state space
%__________________________________________________________________________

clear variables
close all
clc

%% Create Hodgkin–Huxley neuron network object and view its properties

% nueron parameters
gna = 160;
ENa = 0.045;
gk2 = 30;
Ek = -0.07;
gl = 8;
El = -0.046;
tNa = 0.0405;
tk2 = 0.9;
C = 0.5;
Iapp = 0.006;
VshiftK2 = -0.024;

% Create neuron object
neuron = HH_model(gna,ENa,gk2,Ek,gl,El,tNa,tk2,C,Iapp,VshiftK2);

theta = -0.03;
nu = 1000;

synapse = FTM_synapse_model(nu,theta);

% network parameters
g_in =  1e-3*[0 1 1;1 0 1;1 1 0];
g_ex =  zeros(3);
g_el =  zeros(3);


EsynIn = -0.0625;
EsynEx = 0;
N = 3;

% Create network object
netw = CPG(N,neuron,g_in,g_ex,g_el,EsynIn,EsynEx,synapse,synapse);

% Plot network structure
netw.plot

%% Simulate mixed network with odeint integrator with different starting condition
% And get the phase evolution from the time evolution

% Find initial condition of the networks
Nphi = 5;
phi = linspace(0.05,0.95,Nphi);
[phi1,phi2] = meshgrid(phi,phi);
options.Vth = -0.04;
CI = netw.getCIfromPhi([phi1(:),phi2(:)],options);

% Set integrator options
opt.integrator = 'eulero';
opt.integratorOptions.dt = 1e-3;

% Compute state evolution
[T,X] = netw.sim([0,500],CI,opt);



%%
% Get phase evolution

Vth = -0.04;
phiRes = cell(numel(T),1);

for k=1:size(CI,1)
    phiRes{k} = netw.getPhaseRepresentationFromTrack(T{k},X{k},Vth);
end

% Plot phase evolution
plotPhaseSpace(phiRes);

%% Get phase evolution of a network with different network configuration

% Set new network configuration
g_in =  1e-3*[0 0 1;0 0 1;1 1 0];
g_ex = 1e-3*[0 1 0;1 0 0;0 0 0];

netw = netw.set_g_in(g_in);
netw = netw.set_g_ex(g_ex);

% And plot new configuration
netw.plot


% Set integrator option
options.Vth = Vth;
options.integrator = 'odeint';

% Simulate and get phase representation
phi = netw.getPhaseRepresentation(1000,CI,options);

plotPhaseSpace(phi);
