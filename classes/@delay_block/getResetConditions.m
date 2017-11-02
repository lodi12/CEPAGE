function [position,isterminal,direction] = getResetConditions(object,t,y)
position = [0;-1]; % The value that we want to be zero
isterminal = [1;0];  % Halt integration 
direction = [0;0];   % The zero can be approached from either direction

