function [T,X] = sim(object,Tspan,x0,varargin)
% sim    Simulates the neuron
%
%
% [T,X] = sim(OBJ,T,X0)
% Simulates the neuron model OBJ, starting from initial condition X0 between
% Tspan(1) and Tspan(2). Output T is time and output X is the state evolution
%
% SIGNALS = sim(OBJ,T,X0,OPTS)
% A structure OPTS can be provided with % the following fields:
% - integrator: string indicating the solver used to integrate the
%              differential equations. It can be either 'ode45','ode23',
%              'ode113','ode15s','ode23s','ode23t' , 'ode23tb','dde23','odeint' or 'eulero.
%               If you choose eulero the simulation is made in C throught mex
%               file and a supported mex compiler is required. See the
%              documentation of ode45 to get help (doc ode45).
%              Default: 'ode45' for non delayed systems and 'dde23' for delayed systems.
% - integrProp: options to provide to the ODE solver. Type doc odeset to
%               get help. If you choose eulero the variable must have a
%               field dt that describe the integration step size.
%               Default: odeset.
% - x0_delayed: delayed initial condition from lower to higher delays. x0_delayed must be a Ndelay x
%               Nstates matrix. default : 0.
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
    
    nx = object.getnx;
    
    if numel(x0) ~= nx
        error(['x0 must be a vector with ',num2str(nx),' elements']);
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
    
    
    if strcmp(integrator,'eulero')
        if ~isfield(integratorOptions,'dt')
            error('Integrator step dt must be provided when using eulero integrator');
        end
        
        oldFolder = cd;
        mkdir('tmp');
        dt = integratorOptions.dt;
        cd('tmp');
        %     object.generateC('vectorField');
        
        % create vectorField.cpp
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,['\t*vf = new ',object.getCbuilder,';\n}']);
        fclose(fout);
        
        cpth = getcpath();
        copyfile([cpth,'eulero.o']);
        
        
        eval(['mex -silent -c vectorField.cpp -I"',cpth,'" -L"',cpth,'" -lCEPAGE']);
        eval(['mex -silent eulero.o vectorField.o -L"',cpth,'" -lCEPAGE']);
        nStep = (Tspan(2)-Tspan(1))/dt;
        X = eulero(nx,nStep,dt,x0);
        cd(oldFolder);
        rmdir('tmp','s');
        X = X';
        T = linspace(Tspan(1),Tspan(2),nStep);
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
        mkdir('tmp');
        cd('tmp');
        %     object.generateC('vectorField');
        %     movefile('vectorField.c','vectorField.cpp');
        
        % create vectorField.cpp
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,['\t*vf = new ',object.getCbuilder,';\n}']);
        fclose(fout);
        
        cpth = getcpath();
        copyfile([cpth,'odeint.o']);
        
        boostDir = getCEPAGEPar();
        boostDir = boostDir.boostDir;
        
        eval(['mex -silent -c vectorField.cpp -I"',cpth,'" -L"',cpth,'" -lCEPAGE']);
        eval(['mex -silent odeint.o vectorField.o  -I"',boostDir,'/include" -L"',boostDir,'/lib" -L"',cpth,'" -lCEPAGE']);
        [T,X] = odeint(nx,Tspan(2)-Tspan(1),dt,x0);
        T = T+Tspan(1);
        cd(oldFolder);
        rmdir('tmp','s');
        
    else
        %     integratorOptions.Jacobian = @object.getJacobian;
        integratorOptions.Events = @object.getResetConditions;
        
        currentT = Tspan(1);
        
        T = [0];
        X = [zeros(1,object.nx)];
        
        while(currentT < Tspan(2))
            
            command = [integrator,'(@object.getXdot,[currentT Tspan(2)],x0,integratorOptions);'];
            [Ttmp,Xtmp] = eval(command);
            [x0,object] = object.resetStates(Ttmp(end),Xtmp(end,:));
            
            T(end) = [];
            X(end,:) = [];
            
            T = [T;Ttmp];
            X = [X;Xtmp];
            
            currentT = T(end);
            
        end
        
    end
    
else
    
    %% TIME DELAYED
    
    possibleSolver = {'dde23','eulero','odeint'};
    
    nx = object.getnx;
    
    if numel(x0) ~= nx
        error(['x0 must be a vector with ',num2str(nx),' elements']);
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
    
    
    nDelay = numel(object.delays);
    if nargin == 4
        if isfield(varargin{1},'x0_delayed')
            x0_del = varargin{1}.x0_delayed;
            
            if(any(size(x0_del) ~= [nDelay,nx]))
                if (all(size(x0_del') == [nDelay,nx]))
                    x0_del = x0_del';
                else
                    x0_del = zeros(nx,nDelay);
                end
            end
            
        else
            x0_del = zeros(nx,nDelay);
        end
        
    else
        x0_del = zeros(nx,nDelay);
    end
    
    
    if strcmp(integrator,'eulero')
        if ~isfield(integratorOptions,'dt')
            error('Integrator step dt must be provided when using eulero integrator');
        end
        
        oldFolder = cd;
        mkdir('tmp');
        dt = integratorOptions.dt;
        cd('tmp');
        %     object.generateC('vectorField');
        
        % create vectorField.cpp
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,['\t*vf = new ',object.getCbuilder,';\n}']);
        fclose(fout);
        
        cpth = getcpath();
        copyfile([cpth,'eulero_delayed.o']);
        
        
        eval(['mex -silent -c vectorField.cpp -I"',cpth,'" -L"',cpth,'" -lCEPAGE']);
        eval(['mex -silent eulero_delayed.o vectorField.o -L"',cpth,'" -lCEPAGE']);
        nStep = (Tspan(2)-Tspan(1))/dt;
        X = eulero_delayed(nx,nStep,dt,x0,x0_del');
        cd(oldFolder);
        rmdir('tmp','s');
        X = X';
        T = linspace(Tspan(1),Tspan(2),nStep);
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
        mkdir('tmp');
        cd('tmp');
        %     object.generateC('vectorField');
        %     movefile('vectorField.c','vectorField.cpp');
        
        % create vectorField.cpp
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,['\t*vf = new ',object.getCbuilder,';\n}']);
        fclose(fout);
        
        cpth = getcpath();
        copyfile([cpth,'odeint_delayed.o']);
        
        boostDir = getCEPAGEPar();
        boostDir = boostDir.boostDir;
        
        eval(['mex -silent -c vectorField.cpp -I"',cpth,'" -L"',cpth,'" -lCEPAGE']);
        eval(['mex -silent odeint_delayed.o vectorField.o  -I"',boostDir,'/include" -L"',boostDir,'/lib" -L"',cpth,'" -lCEPAGE']);
        [T,X] = odeint_delayed(nx,Tspan(2)-Tspan(1),dt,x0,x0_del);
        T = T+Tspan(1);
        cd(oldFolder);
        rmdir('tmp','s');
        
    else
        %     integratorOptions.Jacobian = @object.getJacobian;
        integratorOptions.Events = @object.getResetConditions;
        
        currentT = Tspan(1);
        
        delays = object.delays(:);
        
        
        T = [-delays(end:-1:1);0]';
        X = [x0_del(:,end:-1:1),x0(:)]';
        
        
        jumps = [];
        inity0 = [];
        
        while(currentT < Tspan(2))
            
            hyst = @(t)object.hystoryFun(t,T,X);
            integratorOptions.Jumps = jumps;
            integratorOptions.InitialY = inity0;
            command = [integrator,'(@object.getXdot,[delays],hyst, [currentT Tspan(2)],integratorOptions);'];
            sol = eval(command);
            
            Ttmp = sol.x';
            Xtmp = sol.y';
            
            T(end) = [];
            X(end,:) = [];
            
            T = [T(:);Ttmp(:)];
            X = [X;Xtmp];
            
            Zpast = object.hystoryFun(T(end)-delays,T,X);
            
            [x0,object] = object.resetStates(Ttmp(end),Xtmp(end,:),Zpast);
            
            
            
            currentT = T(end);
            
            jumps = [jumps;currentT];
            inity0 = x0;
        end
        
    end
    
    
end

end

