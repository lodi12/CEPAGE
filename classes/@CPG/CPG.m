classdef CPG
    
    % CPG   Generic neuron network model
    %
    % This object represents a generic Central Pattern Generator model.
    % The network is characterized by a set of states (x)
    %
    % Each neuron i-th of the network is influenced by neuron j-th by mean
    % the synapsis current:
    %
    %   $$ Isyn = g_{in}(i,j)\frac{V_i-EsynIn}{1+e^{-uSyn(V_j-tetaSyn)}} +
    %   g_{ex}(i,j)\frac{V_i-EsynEx}{1+e^{-uSyn(V_j-tetaSyn)}} +
    %   g_{el}(V_j-V_i) $$
    %
    % 	OBJ = CPG()
    % 	Builds a CPG object OBJ with all parameters equal to 0.
    %
    % 	OBJ = CPG(g_in,g_ex,g_el,Esyn_In,Esyn_Ex,inhActivation,excActivation)
    % 	Builds a CPG object OBJ with user assigned parameters value
    %
    %   CPG methods:
    %   getPhaseRepresentation - Gets the phase evolution of the CPG%
    %   getCIfromPhi - Get the initial condition of the CPG
    %   getApproxPhaseRepresentation - Get the phase difference evolution
    %                                  of the network using standard
    %                                  phase reduction method
    %   getXdot - computes the derivative of the state
    %   sim - simulates the neuron.
    %   simplot - simulates the neuron system and plots time evolution of
    %             states.
    %   disp - displays some information about the HR_model object
    
    %   getN - Gets the number of neurons in the CPG
    %   get_b - gets the value of parameter b
    %   get_s - gets the value of parameter s
    %   get_I - gets the value of parameter I
    %   get_x_rest - gets the value of parameter x_{rest}
    %   get_mu - gets the value of parameter \mu
    %   set_b - sets the value of parameter b
    %   set_s - sets the value of parameter s
    %   set_I - sets the value of parameter I
    %   set_x_rest - sets the value of parameter x_{rest}
    %   set_mu - sets the value of parameter \mu
    %
    % See also HH_model, neuron_model
    %
    % Contributors:
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
    
    %Properties
    
    properties (Access = protected)
        N = 0; % number of neurons in the network
        neurons = []; % Array containing neuron object
        g_in = []; % matrix describing chimical inibhitory synapsis link
        g_ex = [];  % matrix describing chimical excitatory synapsis link
        g_el = [];  % matrix describing electrical synapsis link
        EsynIn = 0; % Chimical inhibitory synapsis reverse potential
        EsynEx = 0;  % Chimical excitatory synapsis reverse potential
        inhActivation = {};
        excActivation = {};
        
        totState = 0;
        incrementalIndexState = [];
        
    end
    
    %methods
    methods
        
        %Constructor
        function object = CPG(N,varargin)
            
            object.N = N;
            object.neurons = cell(N,1);
            object.g_in = zeros(N,N);
            object.g_ex = zeros(N,N);
            object.g_el = zeros(N,N);
            object.inhActivation = cell(N,N);
            object.excActivation = cell(N,N);
            
            object.incrementalIndexState = zeros(N+2*N*N+1,1);
            object.incrementalIndexState(1) = 1;
            if nargin >= 2
                %set neurons model
                nModel = varargin{1};
                if isscalar(nModel) && isa(nModel,'neuron_model')
                    
                    for i=1:N
                        object.neurons{i} = nModel;
                        object.incrementalIndexState(i+1) = object.incrementalIndexState(i)+nModel.getnx;
                    end
                elseif numel(nModel) == N && iscell(nModel)
                    for i=1:N
                        if ~isa(nModel{i},'neuron_model')
                            error('Error in object input');
                        end
                        object.neurons{i} = nModel{i};
                        object.incrementalIndexState(i+1) = object.incrementalIndexState(i)+nModel{i}.getnx;
                    end
                else
                    error('Error in object input');
                end
            end
            
            if nargin == 9
                %set synaps strenght
                object.g_in = varargin{2};
                object.g_ex = varargin{3};
                object.g_el = varargin{4};
                object.EsynIn = varargin{5};
                object.EsynEx = varargin{6};
                
                act = varargin{7};
                
                ii = 2;
                
                if isscalar(act) && isa(act,'synapse_model')
                    for i=1:N
                        for j=1:N
                            object.inhActivation{i,j} = act;
                            object.incrementalIndexState(N+ii) = object.incrementalIndexState(N+ii-1)+act.getnx;
                            ii = ii+1;
                        end
                    end
                elseif all(size(act) == [N,N]) && iscell(act)
                    for i=1:N
                        for j=1:N
                            if ~isa(act{i,j},'synapse_model')
                                error('Error in object input');
                            end
                            
                            object.inhActivation{i,j} = act{i,j};
                            object.incrementalIndexState(N+ii) = object.incrementalIndexState(N+ii-1)+act.getnx;
                            ii = ii+1;
                        end
                    end
                else
                    error('Error in object input');
                end
                
                
                act = varargin{8};
                
                if isscalar(act) && isa(act,'synapse_model')
                    for i=1:N
                        for j=1:N
                            object.excActivation{i,j} = act;
                            object.incrementalIndexState(N+ii) = object.incrementalIndexState(N+ii-1)+act.getnx;
                            ii = ii+1;
                        end
                    end
                elseif all(size(act) == [N,N]) && iscell(act)
                    for i=1:N
                        for j=1:N
                            if ~isa(act{i,j},'synapse_model')
                                error('Error in object input');
                            end
                            
                            object.excActivation{i,j} = act{i,j};
                            object.incrementalIndexState(N+ii) = object.incrementalIndexState(N+ii-1)+act.getnx;
                            ii = ii+1;
                        end
                    end
                else
                    error('Error in object input');
                end
                
                object.totState = object.incrementalIndexState(end)-1;
                
            end
        end
        
        x_dot = getXdot(object,t,x,Z,Isynglobale);
        J = getJacobian(object,t,x);
        str = getCbuilder(object);
        plot(object);
        disp(object);
        CI = getCIfromPhi(object,deltaPhi,Ttrans);
        phi = getPhaseRepresentationFromTrack(object,T,X,Vth)
        [Te,ie] = getPhaseRepresentation(object,Tspan,CI,varargin )
        [deltaPhiMatrix, deltaPhiDot] = computeApproxVectorField(object,PRC,orbit,varargin)
        writeApproxVectorField(w_in,w_ex,w_el,M_in,M_ex,M_el,varargin);
        % get methods
        N = getN(object);
        g_in = get_g_in(object);
        g_ex = get_g_ex(object);
        g_el = get_g_el(object);
        % set methods
        object = set_g_in(object,g_in);
        object = set_g_ex(object,g_ex);
        object = set_g_el(object,g_el);
        [w_in,w_ex,w_el,M_in,M_ex,M_el] = computeMw(object,PRC,limitCycle,varargin);
        
    end
    
    
    
    %     methods (Access = private)
    %         [value,isterminal,direction] = eventsTh(object,t,y,varargin);
    %     end
end
