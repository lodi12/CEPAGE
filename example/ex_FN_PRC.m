%__________________________________________________________________________
%                                                                         
% ex_FN_PRC.m
%__________________________________________________________________________
%
% In this example a network composed of Fitzhugh-Nagumo neurons is analyzed 
% trought both exact and approximate phase representation.
%__________________________________________________________________________

clear variables
close all
clc

%% Create neuron and network objects and view its properties

% nueron parameters
eps = 0.2;
I = 0.5;

% Create neuron object
neuron = FN_relaxation_model(I,eps);

% Display information about neuron
neuron.disp

% network parameters
N = 3;

g_in = 1e-3*[0 1 1 ; 1 0 1 ; 1 1 0];
g_ex = zeros(N);
g_el = zeros(N);

EsynIn = -1.5;
EsynEx = 0;


nu = 100;
theta = 0;

synapse = FTM_synapse_model(nu,theta);

% Create network object
netw = CPG(N,neuron,g_in,g_ex,g_el,EsynIn,EsynEx,synapse,synapse);

% Plot CPG structure
netw.plot;




%% Compute exact phase evolution of the network

% Compute initial condition of the networks
Nphi = 11;
phi = linspace(0,1,Nphi);
[phi1,phi2] = meshgrid(phi,phi);
options.Vth = 0;
CI = netw.getCIfromPhi([phi1(:) phi2(:)],options);

% Set integrator options
options.Vth = 0;
options.integrator = 'odeint';

% Compute phase evolution
phi = netw.getPhaseRepresentation(40000,CI,options);


plotPhaseSpace(phi);

%% Compute exact phase evolution of the network

% Compute initial condition of the networks
Nphi = 11;
phi = linspace(0,1,Nphi);
[phi1,phi2] = meshgrid(phi,phi);
options.Vth = 0;
CI = netw.getCIfromPhi([phi1(:) phi2(:)],options);

% Change synapses
g_in = 1e-3*[0 0 1 ; 0 0 1 ; 1 1 0];
g_ex = 1e-3*[0 1 0 ; 1 0 0 ; 0 0 0];

netw = netw.set_g_in(g_in);
netw = netw.set_g_ex(g_ex);

netw.plot();

% Set integrator options
options.Vth = 0;
options.integrator = 'eulero';
options.integratorOptions.dt = 1e-2;

% Compute phase evolution
phi = netw.getPhaseRepresentation(10000,CI,options);

h = plotPhaseSpace(phi);
%% Compute approximate phase evolution

% Compute PRC
PRCopt.Vth = 0.5;
[PRC,limitCycle] = neuron.computePRC(800,5000,PRCopt);

% plot PRC
plot(PRC.phi,PRC.PRC(:,1))

% Define approximate phase evolution parameters
opts.nStepPhase = 40;
opts.nStepIntegral = 100;

% Get Vectorial Field of the Phase evolution
[deltaPhiMatrix, deltaPhiDot] = netw.computeApproxVectorField(PRC,limitCycle,opts);

%% Plot and compare results obtained with the two methods

figure(h{2});
hold on

quiver(deltaPhiMatrix(:,1),deltaPhiMatrix(:,2),deltaPhiDot(:,1),deltaPhiDot(:,2));