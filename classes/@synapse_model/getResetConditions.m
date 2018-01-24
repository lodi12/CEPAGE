%%

function [position,isterminal,direction] = getResetConditions(object,t,x,Vpre)
position = -1*ones(object.nx,1); % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction

