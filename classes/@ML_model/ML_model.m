classdef ML_model < neuron_model
    
    % ML_model   Neuron represented by the Morris-Lecar model
    %
    % The Neuron is modeld as follow:
    %   _
    %  |
    %  | dV/dt = (-gl*(V - Vl) - gCa*Minf*(V - VCa) - gK*N*(V - Vk ) + I - I_{ext})/CM
    % <
    %  | dN/dt = lamda_N*(Ninf-N)
    %  |_-
    %
    %   Minf = 0.5*(1 + tanh((V - V1 )/V2 ));
    %   Ninf = 0.5*(1 + tanh((V - V3 )/V4 ));
    %   lamda_N = phi*cosh((V - V3 )/(2*V4) );
    %
    %   CM = 5;
    %   gK = 8;
    %   gl = 2;
    %   VCa = 120;
    %   Vk = -80;
    %   Vl = -60;
    %   V1 = -1.2;
    %   V2 = 18;
    %
    % Ther model has two state variable; V represent the membrane potential.
    % gCa,V3,V4,phi and I are parameters that allows to change neuron behaviours.
    %
    % OBJ = ML_model()
    % Builds an ML_model object OBJ with all parameters equal to 0.
    %
    % OBJ = ML_model(CM,gCa,V3,V4,phi,gl,Vl, I)
    % Builds an ML_model object OBJ with user assigned parameters value
    %
    %
    %   ML_model methods:
    %   getXdot - omputes the derivative of the state
    %   sim - simulates the neuron.
    %   simplot - simulates the neuron system and plots time evolution of
    %             states.
    %   disp - displays some information about the FN_relaxation_model object
    %   get_gCa - gets the value of parameter gCa
    %   get_V3 - gets the value of parameter V3
    %   get_V4 - gets the value of parameter V4
    %   get_I - gets the value of parameter I
    %   set_gCa - sets the value of parameter gCa
    %   set_V3 - sets the value of parameter V3
    %   set_V4 - sets the value of parameter V4
    %   set_I - sets the value of parameter I
    %   getJacobian -Computes the Jacobian of the vector field in point x
    %   generateC - Generates C files for the computation of model vector field
    %   getJacobian -Computes the Jacobian of the vector field in point x
    %   generateC - Generates C files for the computation of model vector field
    %
    % The ML_model object is derived from neuron_model and
    % inherits all its methods.
    %
    % See also FN_relaxation_model, neuron_model
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
    
    %properties
    properties(Access = protected)
        gCa = 0;
        V3 = 0;
        V4 = 0;
        phi = 0;
        I = 0;
        CM = 0;
        gL = 0;
        EL = 0;
    end
    
    methods
        function obj = ML_model(varargin)
            obj.nx = 2;
            obj.xnames = {'V','N'};
            obj.modelName = 'Morris-Lecar';
            obj.isContinuous = true;
            if nargin == 8
                obj.CM = varargin{1};
                obj.gCa = varargin{2};
                obj.V3 = varargin{3};
                obj.V4 = varargin{4};
                obj.phi = varargin{5};
                obj.gL = varargin{6};
                obj.EL = varargin{7};
                obj.I = varargin{8};
            end
        end
        
        % Other methods
        x_dot = getXdot(object,t,x,varargin);
        disp(object);
        str = getCbuilder(object);
        
        % get methods
        CM = get_CM(object);
        gCa = get_gCa(object);
        V3 = get_V3(object);
        V4 = get_V4(object);
        phi = get_phi(object);
        gL = get_gL(object);
        EL = get_EL(object);
        I = get_I(object);
        
        % set methods
        object = set_CM(object,CM);
        object = set_gCa(object,gCa);
        object = set_V3(object,V3);
        object = set_V4(object,V4);
        object = set_phi(object,phi);
        object = set_gL(object,gL);
        object = set_EL(object,EL);
        object = set_I(object,I);
    end
    
end
