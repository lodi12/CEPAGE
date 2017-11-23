function [xreset,object] = resetStates(object,t,x,ie)

xreset = zeros(size(x));

N = object.N;

incrementalIndexState = object.incrementalIndexState;

inhActivation = object.inhActivation;
excActivation = object.excActivation;

neurons = object.neurons;

for i=1:N
    beginI = incrementalIndexState(i);
    endI = incrementalIndexState(i+1)-1;
    
    if any(ie >= beginI & ie <= endI)
        xreset(beginI:endI) = neurons{i}.resetStates(t,x(beginI:endI));
    else
            xreset(beginI:endI) = x(beginI:endI);
        end
end


% synapses reset
ii = 2;
for i=1:N
    for j=1:N
        beginI = incrementalIndexState(N+ii);
        endI = incrementalIndexState(N+ii)-1;
        
        if any(ie >= beginI & ie <= endI)
            xreset(beginI:endI) = inhActivation{i}.resetStates(t,x(beginI:endI));
        else
            xreset(beginI:endI) = x(beginI:endI);
        end
        ii = ii+1;
    end
end

for i=1:N
    for j=1:N
        beginI = incrementalIndexState(N+ii);
        endI = incrementalIndexState(N+ii)-1;
        if any(ie >= beginI & ie <= endI)
            xreset(beginI:endI) = excActivation{i}.resetStates(t,x(beginI:endI));
        else
            xreset(beginI:endI) = x(beginI:endI);
        end
        ii = ii+1;
    end
end

end

