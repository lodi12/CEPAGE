%__________________________________________________________________________
%                                                                         
% ex_HR_model.m
%__________________________________________________________________________
%
% In this example a neuron described by Hindmarsh Rose model is simulated
% throught the tool contained in CEPAGE toolbox                 
%__________________________________________________________________________

clear variables
close all
clc

%% Create Hindmarsh Rose neuron object and view its properties

% nueron parameters
b = 2;
I = 2;
mu = 0.01;
s = 4;
x_rest = -1.6;

% Create neuron object

neuron = HR_model(b,mu,s,I,x_rest);

% Display information about neuron
neuron.disp

%% Simulate and plot neuron behavior - Bursting

% Set starting condition
x0 = [0 0 0];

% Set simulation time
Tspan = [0 1000];

% Run simulation
[T,X] = neuron.sim(Tspan,x0);

%% Simulate and plot neuron behavior with eulero - Quiescence

% Change neuron parameter I
neuron = neuron.set_I(1);

% Set starting condition
x0 = [0 0 0];

% Set simulation time
Tspan = [0 1000];

% Choose integrator and set options
opt.integrator = 'eulero';
opt.integratorOptions.dt = 0.001;

% Run simulation
[T,X] = neuron.simplot(Tspan,x0,opt);

%% Change parameters and plot new simulation with odeint - Spiking

% Change neuron parameters I and b
neuron = neuron.set_I(5);
neuron = neuron.set_b(3);

% Set starting condition
x0 = [0 0 0];

% Set simulation time
Tspan = [0 1000];

% Choose integrator and set options
opt.integrator = 'odeint';
opt.integratorOptions.dt = 0.0001;

% Run simulation
[T,X] = neuron.simplot([0,1000],[0 0 0],opt);
