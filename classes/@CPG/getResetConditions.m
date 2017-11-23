function [position,isterminal,direction] = getResetConditions(object,t,y,ie)

N = object.N;

totState = object.totState;
incrementalIndexState = object.incrementalIndexState;

position = zeros(totState,1); % The value that we want to be zero
isterminal = zeros(totState,1);  % Halt integration
direction = zeros(totState,1);   % The zero can be approached from either direction


inhActivation = object.inhActivation;
excActivation = object.excActivation;

neurons = object.neurons;

for i=1:N
    beginI = incrementalIndexState(i);
    endI = incrementalIndexState(i+1)-1;
    
    
        [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = neurons{i}.getResetConditions(t,y(beginI:endI));
end


% synapses reset
ii = 2;
for i=1:N
    for j=1:N
        beginI = incrementalIndexState(N+ii);
        endI = incrementalIndexState(N+ii)-1;
        
        [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = inhActivation{i}.getResetConditions(t,y(beginI:endI));
        ii = ii+1;
    end
end

for i=1:N
    for j=1:N
        beginI = incrementalIndexState(N+ii);
        endI = incrementalIndexState(N+ii)-1;
        
        [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = excActivation{i}.getResetConditions(t,y(beginI:endI));
        ii = ii+1;
    end
end

