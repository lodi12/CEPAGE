classdef synapse_model
    
    % v   Generic synapse model
    %
    % This object represents a generic synapse model.
    % The synapse is characterized by an activation function f(x)
    %
    %
    % synapse_model methods:
    %   getActivation - Compute the activation function
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
        modelName = '';
    end
    
    methods (Abstract)
        % Other methods
        act = getActivation(object,Vpre);    
        str = getCbuilder(object);
    end
    

    
    methods
        %constructor
        function object = synapse_model(varargin)
            % abstract class empty cconstructor
        end
        
    end
    
end
