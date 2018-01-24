function [position,isterminal,direction] = getResetConditions(object,t,y,varargin)

N = object.N;

totState = object.totState;
incrementalIndexState = object.incrementalIndexState;

position = zeros(totState,1); % The value that we want to be zero
isterminal = zeros(totState,1);  % Halt integration
direction = zeros(totState,1);   % The zero can be approached from either direction


inhActivation = object.inhActivation;
excActivation = object.excActivation;

neurons = object.neurons;


if ~object.is_delayed
    
    for i=1:N
        beginI = incrementalIndexState(i);
        endI = incrementalIndexState(i+1)-1;
        

        
        [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = neurons{i}.getResetConditions(t,y(beginI:endI));
    end
    
    
    % synapses reset
    ii = 1;
    for i=1:N
        for j=1:N
            beginI = incrementalIndexState(N+ii);
            endI = incrementalIndexState(N+ii+1)-1;
            
                    otherI = incrementalIndexState(j);
            
            [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = inhActivation{i,j}.getResetConditions(t,y(beginI:endI),y(otherI));
            ii = ii+1;
        end
    end
    
    for i=1:N
        for j=1:N
            beginI = incrementalIndexState(N+ii);
            endI = incrementalIndexState(N+ii+1)-1;
            
            otherI = incrementalIndexState(j);
            
            [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = excActivation{i,j}.getResetConditions(t,y(beginI:endI),y(otherI));
            ii = ii+1;
        end
    end
    
    
else
    
    
    if nargin == 4
        Z = varargin{1};
    else
        Z = zeros(object.totState,numel(object.delays));
    end
    
    
    delayIndexNeur = object.delayIndexNeur;
    delayInhSyn = object.delayInhSyn;
    delayExcSyn = object.delayExcSyn;
    
    
    for i=1:N
        beginI = incrementalIndexState(i);
        endI = incrementalIndexState(i+1)-1;
        
        Ztmp = Z(:,delayIndexNeur{i});
        
        
        [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = neurons{i}.getResetConditions(t,y(beginI:endI),Ztmp);
    end
    
    
    % synapses reset
    ii = 1;
    for i=1:N
        for j=1:N
            beginI = incrementalIndexState(N+ii);
            endI = incrementalIndexState(N+ii+1)-1;
            
            otherI = incrementalIndexState(j);
            
            Ztmp = Z(otherI,delayInhSyn{i,j});
            
            
            
            [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = inhActivation{i,j}.getResetConditions(t,y(beginI:endI),x(otherI),Ztmp);
            ii = ii+1;
        end
    end
    
    for i=1:N
        for j=1:N
            beginI = incrementalIndexState(N+ii);
            endI = incrementalIndexState(N+ii+1)-1;
            
            otherI = incrementalIndexState(j);
            
            Ztmp = Z(:,delayExcSyn{i,j});
            
            [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = excActivation{i,j}.getResetConditions(t,y(beginI:endI),x(otherI),Ztmp);
            ii = ii+1;
        end
    end
    
end

