classdef Heaviside_synapse_model < synapse_model
    
    % Heaviside synapse model
    %
    % This object represents a synapse described by the Fast threshold 
    % modulation model.
    % The synapse is characterized by an activation function 
    %
    %   f(x) = H(V-theta)
    %
    % where H is the Heaviside function
    %
    % OBJ = FTM_synapse_model()
    % Builds an FTM_synapse_model object OBJ with all parameters equal to 0.
    %
    % OBJ = Heaviside_synapse_model(theta)
    % Builds an Heaviside_synapse_model object OBJ with user assigned parameters value
    %
    % Heaviside_synapse_model methods:
    %   getActivation - Compute the activation function
    %
    % See also synapse_model
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
        theta = 0;
    end
       

    
    methods
        %constructor
        function object = Heaviside_synapse_model(varargin)
            if nargin == 0
                object.theta = 0;
                object.nx = 0;
                object.modelName = 'Heaviside_synapse';
                object.isContinuous = true;
            elseif nargin == 1
                object.theta = varargin{1};
                object.nx = 0;
                object.modelName = 'Heaviside_synapse';
                object.isContinuous = true;
            else
                error('Wrong number of input arguments');
            end
        end
        
        act = getActivation(object,Vpre,varargin);    
        str = getCbuilder(object);
        disp(object);
        x_dot = getXdot(object,t,x,varargin);
    end
    
end
