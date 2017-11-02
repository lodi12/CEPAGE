classdef neuron_model
    
    % neuron_model   Generic neuron model
    %
    % This object represents a generic neuron model.
    % The neuron is characterized by a set of states (x)
    %
    % The continuous time formulation is:
    %       dx/dt = f(x,Iext)
    % where Iext is an external current that act on the neuron.
    % The first state variable is the membrane potential.
    % This is an abstract class which cannot be istantiated.
    %
    % neuron_model methods:
    %   computePRC - Compute the Phase resetting curve for the neuron
    %   getStateNames - gets the mnemonical name of the neuron states
    %   getnx - Gets the number of state varible nx
    %
    % See also HH_model, neuron_network
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
        nx = 0 % number of state variable
        xnames = cell(0); % States name
        modelName = '';
        isContinuous = false;
    end
    
    methods (Abstract)
        % Other methods
        str = getCbuilder(object);    
        x_dot = getXdot(object,t,x,varargin);
        J = getJacobian(object,x);
    end
    

    
    methods
        %constructor
        function object = neuron_model(nx,varargin)
            % abstract class empty cconstructor
        end
        
        
        nx = getnx(object);
        names = getStateNames(object)
        [phi,PRC] = computePRC(object,nPoints,Ttrans,varargin);
        [T,X] = sim(object,Tspan,x0,varargin);
        [T,X] = simplot(object,Tspan,x0,varargin);
        
       
        
         [position,isterminal,direction] = getResetConditions(object,t,y);
        [xreset,object] = resetStates(object,t,x);
        
    end
    
end
