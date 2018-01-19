classdef synapse_model
    
    % v   Generic synapse model
    %
    % This object represents a generic synapse model.
    % The synapse is characterized by an activation function f(x)
    %
    %
    % synapse_model methods:
    %   getActivation - Compute the activation function
    %   getXdot       - Compute the synapse states derivative
    %
    % See also FTM_synapse_model
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
        delays = [];
    end
    
    methods (Abstract)
        % Other methods
        act = getActivation(object,Vpre,varargin);    
        x_dot = getXdot(object,t,x,varargin);
        str = getCbuilder(object);
    end
    

    
    methods
        %constructor
        function object = synapse_model(varargin)
            % abstract class empty cconstructor
        end
        
        
        nx = getnx(object);
        names = getStateNames(object);
        [position,isterminal,direction] = getResetConditions(object,t,y,varargin);
        [xreset,object] = resetStates(object,t,x,varargin);
        cont = is_continuous(object);
        del = is_delayed(object);
        del = getDelays(object)
    end
    
end
