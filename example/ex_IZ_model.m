%__________________________________________________________________________
%                                                                         
% ex_IZ_model.m
%__________________________________________________________________________
%
% In this example a neuron described by Izhikevich model is simulated
% throught the tool contained in SFNN toolbox                 
%__________________________________________________________________________


clear variables
close all
clc

%% Create Izhikevich neuron object and view its properties
a = 0.02;
b = 0.2;
    
I = 10;
    
c = -50;
d = 2;
    
gL = 0;
El = -60;
    
neuron = IZ_model(a,b,c,d,I,gL,El);

neuron.disp();
    
%% Simulate and plot neuron behavior - Bursting    

% Set starting condition
x0 = [0 0];

% Set simulation time
Tspan = [0 1000];

% Choose integrator and set options
opt.integrator = 'eulero';
opt.integratorOptions.dt = 0.001;

% Run simulation
[T1,X1] = neuron.simplot(Tspan,x0,opt);

%% Create CPG

% network parameters
g_in = 1e-2*[0 1;1 0];
g_ex = zeros(2);
g_el = zeros(2);

EsynIn = -80;
EsynEx = 0;
N = 2;

theta = -50;
nu = 100;

synapse = FTM_synapse_model(nu,theta);

% Create network object
netw = CPG(N,neuron,g_in,g_ex,g_el,EsynIn,EsynEx,synapse,synapse);

% Set starting condition
x0 = [0 0 -50 0];

% Set simulation time
Tspan = [0 10000];

% Choose integrator and set options
opt.integrator = 'odeint';
opt.integratorOptions.dt = 0.001;
opt.Vth = -60;
% Run simulation

[T2,X2] = netw.simplot(Tspan,x0,opt);
% phi = netw.getPhaseRepresentation(Tspan,x0,opt);
