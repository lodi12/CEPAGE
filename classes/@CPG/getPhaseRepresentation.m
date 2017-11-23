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
    
    %     totalMexString = '';
    %
    %     for i=1:N
    %         nm = ['neuron',num2str(i)];
    %         object.neurons{i}.generateC(nm);
    %         eval(['mex -c ',nm,'.c'])
    %         totalMexString = [totalMexString,' ',nm,'.o'];
    %     end
    %     object.generateC('vectorField');
    %     mex -c vectorField.c;
    
    fout = fopen('vectorField.cpp','w');
    fprintf(fout,'#include "vectorField.h"\n');
    fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
    fprintf(fout,[object.getCbuilder,';\n}']);
    fclose(fout);
    
    cpth = getcpath();
    copyfile([cpth,'euleroEvents.o']);
    
    eval(['mex -silent -c vectorField.cpp -I',cpth,' -L',cpth,' -lCEPAGE']);
    str = ['mex euleroEvents.o vectorField.o -L',cpth,' -lCEPAGE'];
    eval(str);
    nStep = (Tspan(2)-Tspan(1))/dt;
    
    if size(CI,1) > 1
        phi = cell(size(CI,1),1);
        parfor kk = 1:size(CI,1)
            phi{kk} = mod(euleroEvents(Nstati,nStep,dt,CI(kk,:),Vth,stopThreshold,N),1);
        end
    else
        phi = mod(euleroEvents(Nstati,nStep,dt,CI,Vth,stopThreshold,N),1);
        
    end
    cd(oldFolder);
    rmdir(nameFolder,'s');
elseif strcmp(integrator,'odeint')
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
    
    %     totalMexString = '';
    %
    %     for i=1:N
    %         nm = ['neuron',num2str(i)];
    %         object.neurons{i}.generateC(nm);
    %         movefile([nm,'.c'],[nm,'.cpp']);
    %         eval(['mex -c ',nm,'.cpp'])
    %         totalMexString = [totalMexString,' ',nm,'.o'];
    %     end
    %     object.generateC('vectorField');
    %     movefile('vectorField.c','vectorField.cpp');
    %     mex -c vectorField.cpp;
    
    fout = fopen('vectorField.cpp','w');
    fprintf(fout,'#include "vectorField.h"\n');
    fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
    fprintf(fout,[object.getCbuilder,';\n}']);
    fclose(fout);
    
    boostDir = getCEPAGEPar();
    boostDir = boostDir.boostDir;
    
    cpth = getcpath();
    copyfile([cpth,'odeintEvents.o']);
    
    eval(['mex -silent -c vectorField.cpp -I',cpth,' -L',cpth,' -lCEPAGE']);
    str = ['mex odeintEvents.o vectorField.o -I',boostDir,'/include -L',boostDir,'/lib -L',cpth,' -lCEPAGE'];
    eval(str);
    
    if size(CI,1) > 1
        phi = cell(size(CI,1),1);
        parfor kk = 1:size(CI,1)
            phi{kk} = mod(odeintEvents(Nstati,Tspan(2),dt,CI(kk,:),Vth,N),1);
        end
    else
        phi = mod(odeintEvents(Nstati,Tspan(2),dt,CI,Vth,N),1);
    end
    cd(oldFolder);
    rmdir(nameFolder,'s');
    
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



% N = object.N;
% nx = object.neurons{1}.getnx;
%
% nc = size(X,2);
%
% if nc ~= N*nx
%     error(['X must have ',num2str(N*nx),'columns']);
% end
%
% ii = (0:N-1)*nx+1;
%
% X = X(:,ii);
%
% Te = cell(N,1);
%
% low = X < Vth;
% up = X >= Vth;
%
% len = Inf;
%
% for i=1:N
%     evDif = low(1:end-1,i)-up(2:end,i);
%     evSum = low(1:end-1,i)+up(2:end,i);
%     Te{i} = T(evDif == 0 & evSum == 2);
%     len = min(len,numel(Te{i}));
% end
%
% T1 = Te{1}(end)-Te{1}(end-1);
%
% phi = cell(N-1,1);
%
% for i=2:N
%     phi{i-1} = mod((Te{i}(1:len)-Te{1}(1:len))/T1,1);
% end


