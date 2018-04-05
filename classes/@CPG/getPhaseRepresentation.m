function phi = getPhaseRepresentation(object,Tspan,CI,varargin )
% getPhaseRepresentation   Gets the phase evolution of the network
%
% phi = getPhaseRepresentation(object,T,CI )
% phi is a cell array with N-1 elements describing the evolution of the
% phase difference between neuron 1 and neuron i+1
% The network is simulated for T seconds.
% CI is the initial conditions of the networks and must be a 1xN array. If
% CI is a Ns x N matrix the function compute the phase evolutions of the
% networks starting from the Ns starting conditions
%
% phi = getPhaseRepresentation(object,T,CI,OPTS)
% A structure OPTS can be provided with % the following fields:
% - integrator: string indicating the solver used to integrate the
%              differential equations. It can be either 'ode45','ode23',
%              'ode113','ode15s','ode23s','ode23t' , 'ode23tb', 'odeint' or 'eulero.
%               If you choose eulero or odeint the simulation is made in C throught mex
%               file and a supported mex compiler is required.
%              Default: 'ode45'.
% - integratorOptions: options to provide to the ODE solver. Type doc odeset to
%               get help. If you choose eulero the variable must have a
%               field dt that describe the integration step size.
%               Default: odeset.
% - Vth :  is the value of the membrane potential that correspond to an
%        event (it must be between max and minimum value of the membrane
%        potential)
%        Default : 0
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


if ~object.is_delayed
    
    possibleSolver = {'ode45','ode23','ode113','ode15s','ode23s','ode23t','ode23tb','eulero','odeint'};
    
    % nx = object.neurons{1}.getnx;
    N  = object.N;
    Nstati = object.totState;
    
    if size(CI,2) ~= Nstati
        error(['initial conditions must be a vector with ',num2str(Nstati),' columns']);
    end
    
    if numel(Tspan) == 1
        Tspan(2) = Tspan(1);
        Tspan(1) = 0;
    elseif numel(Tspan) > 2
        error('Tspan must be a vector with two elements [Tstart, Tend]');
    end
    
    integrator = 'ode45';
    integratorOptions = odeset;
    
    if nargin == 4
        if isfield(varargin{1},'integratorOptions')
            integratorOptions = varargin{1}.integratorOptions;
        end
        if isfield(varargin{1},'integrator')
            integrator = varargin{1}.integrator;
        end
    end
    
    intOk = -1;
    for i=1:numel(possibleSolver)
        if strcmp(integrator,possibleSolver{i})
            intOk = 1;
            break;
        end
    end
    
    if intOk == -1
        error('Choosen integrator not allowed');
    end
    
    
    Vth = 0;
    stopThreshold = eps;
    if nargin == 4
        if isfield(varargin{1},'Vth')
            Vth = varargin{1}.Vth;
        end
        if isfield(varargin{1},'stopThreshold')
            stopThreshold = varargin{1}.stopThreshold;
        end
    end
    
    
    
    
    
    if strcmp(integrator,'eulero')
        if ~isfield(integratorOptions,'dt')
            error('Integrator step dt must be provided when using eulero integrator');
        end
        oldFolder = cd;
        
        ii = 0;
        nameFolder = 'tmp';
        
        while(exist(nameFolder,'dir') == 7)
            ii = ii+1;
            nameFolder = ['tmp',num2str(ii)];
        end
        
        mkdir(nameFolder);
        dt = integratorOptions.dt;
        cd(nameFolder);
        
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,[object.getCbuilder,';\n}']);
        fclose(fout);
        
        cpth = getcpath();
        
        if ispc
            copyfile([cpth,'euleroEvents.obj']);
            eval(['mex -silent euleroEvents.obj vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output euleroEvents']);
        elseif isunix
            copyfile([cpth,'euleroEvents.o']);
            eval(['mex -silent euleroEvents.o vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output euleroEvents']);
        else
            error('Unsopported operative system');
        end
        
        
        nStep = (Tspan(2)-Tspan(1))/dt;
        
        if size(CI,1) > 1
            phi = cell(size(CI,1),1);
            parfor kk = 1:size(CI,1)
                phi{kk} = mod(euleroEvents(Nstati,nStep,dt,CI(kk,:),Vth,stopThreshold,N),1);
            end
        else
            phi = mod(euleroEvents(Nstati,nStep,dt,CI,Vth,stopThreshold,N),1);
            
        end
        
        clear euleroEvents;
        
        cd(oldFolder);
        a = 0;
        tic
        while (a == 0) && (toc < 5)
            a = rmdir(nameFolder,'s');
        end
        
        if a ~= 1
            warning(['Cannot delete ',nameFolder,' folder, remove manually!']);
        end
        
        
    elseif strcmp(integrator,'odeint')
        
        useBoost = getCEPAGEPar();
        useBoost = useBoost.useBoost;
        
        if ~useBoost
            error('Boost c++ integrator is not installed');
        end
        
        if ~isfield(integratorOptions,'dt')
            dt = 1e-3;
        else
            dt = integratorOptions.dt;
        end
        oldFolder = cd;
        
        ii = 0;
        nameFolder = 'tmp';
        
        while(exist(nameFolder,'dir') == 7)
            ii = ii+1;
            nameFolder = ['tmp',num2str(ii)];
        end
        
        mkdir(nameFolder);
        cd(nameFolder);
        
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,[object.getCbuilder,';\n}']);
        fclose(fout);
        
        boostDir = getCEPAGEPar();
        boostDir = boostDir.boostDir;
        
        cpth = getcpath();
        
        if ispc
            copyfile([cpth,'odeintEvents.obj']);
            eval(['mex -silent odeintEvents.obj vectorField.cpp "-I',cpth,'" "-I',boostDir,'/include" "-L',boostDir,'/lib" -L"',cpth,'" -lCEPAGE -output odeintEvents']);
        elseif isunix
            copyfile([cpth,'odeintEvents.o']);
            eval(['mex -silent odeintEvents.o vectorField.cpp "-I',cpth,'" "-I',boostDir,'/include" "-L',boostDir,'/lib" -L"',cpth,'" -lCEPAGE -output odeintEvents']);
        else
            error('Unsopported operative system');
        end
        
        if size(CI,1) > 1
            phi = cell(size(CI,1),1);
            parfor kk = 1:size(CI,1)
                phi{kk} = mod(odeintEvents(Nstati,Tspan(2),dt,CI(kk,:),Vth,N),1);
            end
        else
            phi = mod(odeintEvents(Nstati,Tspan(2),dt,CI,Vth,N),1);
        end
        
        clear odeintEvents;
        
        cd(oldFolder);
        a = 0;
        tic
        while (a == 0) && (toc < 5)
            a = rmdir(nameFolder,'s');
        end
        
        if a ~= 1
            warning(['Cannot delete ',nameFolder,' folder, remove manually!']);
        end
        
    else
        
        if ~object.is_continuous
            error(['Pahse representation employing MATLAB integrator '...
                'can be used only for continuous system']);
        end
        
        Te = cell(size(CI,1),1);
        ie = cell(size(CI,1),1);
        eventFun = @(T, Y) object.eventsTh(T, Y, Vth);
        integratorOptions.Events = eventFun;
        
        
        if size(CI,1) > 1
            phi = cell(size(CI,1),1);
        end
        
        for i=1:size(CI,1)
            x0  = CI(i,:);
            
            command = [integrator,'(@object.getXdot,[Tspan(1) Tspan(2)],x0,integratorOptions);'];
            [~,~,Te,~,ie] = eval(command);
            
            
            len = Inf;
            
            for n=1:N
                len = min([len,sum(ie == n)]);
            end
            
            phiTmp = zeros(len-1,N-1);
            
            T1 = Te(ie == 1);
            
            for n=2:N
                Tn = Te(ie == n);
                phiTmp(:,n-1) = mod((Tn(2:len)-T1(2:len))./(T1(2:len)-T1(1:len-1)),1);
            end
            
            
            if size(CI,1) > 1
                phi{i} = phiTmp;
            else
                phi = phiTmp;
            end
            
        end
        
        
        
    end
    
    
else
    
    
    possibleSolver = {'dde23','eulero','odeint'};
    
    % nx = object.neurons{1}.getnx;
    N  = object.N;
    Nstati = object.totState;
    
    if size(CI,2) ~= Nstati
        error(['initial conditions must be a vector with ',num2str(Nstati),' columns']);
    end
    
    if numel(Tspan) == 1
        Tspan(2) = Tspan(1);
        Tspan(1) = 0;
    elseif numel(Tspan) > 2
        error('Tspan must be a vector with two elements [Tstart, Tend]');
    end
    
    integrator = 'dde23';
    integratorOptions = ddeset;
    
    if nargin == 4
        if isfield(varargin{1},'integratorOptions')
            integratorOptions = varargin{1}.integratorOptions;
        end
        if isfield(varargin{1},'integrator')
            integrator = varargin{1}.integrator;
        end
    end
    
    intOk = -1;
    for i=1:numel(possibleSolver)
        if strcmp(integrator,possibleSolver{i})
            intOk = 1;
            break;
        end
    end
    
    if intOk == -1
        error('Choosen integrator not allowed');
    end
    
    
    Vth = 0;
    stopThreshold = eps;
    if nargin == 4
        if isfield(varargin{1},'Vth')
            Vth = varargin{1}.Vth;
        end
        if isfield(varargin{1},'stopThreshold')
            stopThreshold = varargin{1}.stopThreshold;
        end
    end
    
    nDelay = numel(object.delays);
    if nargin == 4
        if isfield(varargin{1},'x0_delayed')
            x0_del = varargin{1}.x0_delayed;
            
            if(any(size(x0_del) ~= [nDelay,object.totState]))
                if (all(size(x0_del') == [nDelay,object.totState]))
                    x0_del = x0_del';
                else
                    x0_del = zeros(object.totState,nDelay);
                end
            end
            
        else
            x0_del = zeros(object.totState,nDelay);
        end
        
    else
        x0_del = zeros(object.totState,nDelay);
    end
    
    
    
    if strcmp(integrator,'eulero')
        if ~isfield(integratorOptions,'dt')
            error('Integrator step dt must be provided when using eulero integrator');
        end
        oldFolder = cd;
        
        ii = 0;
        nameFolder = 'tmp';
        
        while(exist(nameFolder,'dir') == 7)
            ii = ii+1;
            nameFolder = ['tmp',num2str(ii)];
        end
        
        mkdir(nameFolder);
        dt = integratorOptions.dt;
        cd(nameFolder);
        
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,[object.getCbuilder,';\n}']);
        fclose(fout);
        
        cpth = getcpath();
        
        if ispc
            copyfile([cpth,'euleroEvents_delayed.obj']);
            eval(['mex -silent euleroEvents_delayed.obj vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output euleroEvents_delayed']);
        elseif isunix
            copyfile([cpth,'euleroEvents_delayed.o']);
            eval(['mex -silent euleroEvents_delayed.o vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output euleroEvents_delayed']);
        else
            error('Unsopported operative system');
        end
        
        nStep = (Tspan(2)-Tspan(1))/dt;
        
        if size(CI,1) > 1
            phi = cell(size(CI,1),1);
            parfor kk = 1:size(CI,1)
                phi{kk} = mod(euleroEvents_delayed(Nstati,nStep,dt,CI(kk,:),Vth,stopThreshold,N,x0_del),1);
            end
        else
            phi = mod(euleroEvents_delayed(Nstati,nStep,dt,CI,Vth,stopThreshold,N,x0_del),1);
            
        end
        
        clear euleroEvents_delayed;
        
        cd(oldFolder);
        a = 0;
        tic
        while (a == 0) && (toc < 5)
            a = rmdir(nameFolder,'s');
        end
        
        if a ~= 1
            warning(['Cannot delete ',nameFolder,' folder, remove manually!']);
        end
        
    elseif strcmp(integrator,'odeint')
        
        useBoost = getCEPAGEPar();
        useBoost = useBoost.useBoost;
        
        if ~useBoost
            error('Boost c++ integrator is not installed');
        end
        
        if ~isfield(integratorOptions,'dt')
            dt = 1e-3;
        else
            dt = integratorOptions.dt;
        end
        oldFolder = cd;
        
        ii = 0;
        nameFolder = 'tmp';
        
        while(exist(nameFolder,'dir') == 7)
            ii = ii+1;
            nameFolder = ['tmp',num2str(ii)];
        end
        
        mkdir(nameFolder);
        cd(nameFolder);
        
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,[object.getCbuilder,';\n}']);
        fclose(fout);
        
        boostDir = getCEPAGEPar();
        boostDir = boostDir.boostDir;
        
        cpth = getcpath();
        
        if ispc
            copyfile([cpth,'odeintEvents_delayed.obj']);
            eval(['mex -silent odeintEvents_delayed.obj vectorField.cpp "-I',cpth,'" "-I',boostDir,'/include" "-L',boostDir,'/lib" -L"',cpth,'" -lCEPAGE -output odeintEvents_delayed']);
        elseif isunix
            copyfile([cpth,'odeintEvents_delayed.o']);
            eval(['mex -silent odeintEvents_delayed.o vectorField.cpp "-I',cpth,'" "-I',boostDir,'/include" "-L',boostDir,'/lib" -L"',cpth,'" -lCEPAGE -output odeintEvents_delayed']);
        else
            error('Unsopported operative system');
        end
        
        if size(CI,1) > 1
            phi = cell(size(CI,1),1);
            parfor kk = 1:size(CI,1)
                phi{kk} = mod(odeintEvents_delayed(Nstati,Tspan(2),dt,CI(kk,:),Vth,N),1);
            end
        else
            phi = mod(odeintEvents_delayed(Nstati,Tspan(2),dt,CI,Vth,N),1);
        end
        
        clear odeintEvents;
        
        cd(oldFolder);
        a = 0;
        tic
        while (a == 0) && (toc < 5)
            a = rmdir(nameFolder,'s');
        end
        
        if a ~= 1
            warning(['Cannot delete ',nameFolder,' folder, remove manually!']);
        end
        
    else
        
        if ~object.is_continuous
            error(['Pahse representation employing MATLAB integrator '...
                'can be used only for continuous system']);
        end
        
        Te = cell(size(CI,1),1);
        ie = cell(size(CI,1),1);
        eventFun = @(T, Y) object.eventsTh(T, Y, Vth);
        integratorOptions.Events = eventFun;
        
        
        if size(CI,1) > 1
            phi = cell(size(CI,1),1);
        end
        
        for i=1:size(CI,1)
            x0  = CI(i,:);
            
            command = [integrator,'(@object.getXdot,[Tspan(1) Tspan(2)],x0,integratorOptions);'];
            [~,~,Te,~,ie] = eval(command);
            
            
            len = Inf;
            
            for n=1:N
                len = min([len,sum(ie == n)]);
            end
            
            phiTmp = zeros(len-1,N-1);
            
            T1 = Te(ie == 1);
            
            for n=2:N
                Tn = Te(ie == n);
                phiTmp(:,n-1) = mod((Tn(2:len)-T1(2:len))./(T1(2:len)-T1(1:len-1)),1);
            end
            
            
            if size(CI,1) > 1
                phi{i} = phiTmp;
            else
                phi = phiTmp;
            end
            
        end
        
        
        
    end
    
    
end

