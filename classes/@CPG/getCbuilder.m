function str = getCbuilder(object)

    N = object.N;
    neurons = object.neurons;

    inhActivation = object.inhActivation;
    excActivation = object.excActivation;
    
    g_in = object.g_in;
    g_ex = object.g_ex;
    g_el = object.g_el;
    
    str = [];
    
    Ninh = sum(g_in(:) ~= 0);
    Nexc = sum(g_ex(:) ~= 0);
    Nel = sum(g_el(:) ~= 0);
    

    
    
    str = [str,'neuron_model **neuron;\n\n'];
    str = [str,'CPG::synStruct_t **inhSyn;\n'];
    str = [str,'CPG::synStruct_t **excSyn;\n\n'];
    str = [str,'CPG::synStruct_t **elSyn;\n\n'];
    
    str = [str,sprintf('neuron = new neuron_model*[%d];\n',N)];
    str = [str,sprintf('inhSyn = new CPG::synStruct*[%d];\n',Ninh)];
    str = [str,sprintf('excSyn = new CPG::synStruct*[%d];\n',Nexc)];
    str = [str,sprintf('elSyn = new CPG::synStruct*[%d];\n',Nel)];

  
    
    for i=1:N
    str = [str,'neuron[',num2str(i-1),'] = new ',neurons{i}.getCbuilder(),';\n'];
    end
    
    str = [str,'\n'];
    
    
    ii = 0;
    
    for i=1:N
        for j=1:N
            if g_in(i,j) ~= 0
                str = [str,'inhSyn[',num2str(ii),'] = new CPG::synStruct(',num2str(i-1),',',num2str(j-1),',',num2str(g_in(i,j)),',',num2str(object.EsynIn),', new ',inhActivation{i,j}.getCbuilder(),');\n'];
                ii = ii + 1;
            end
        end
    end
    
    str = [str,'\n'];
    
    ii = 0;
    
    for i=1:N
        for j=1:N
            if g_ex(i,j) ~= 0
                str = [str,'excSyn[',num2str(ii),'] = new CPG::synStruct(',num2str(i-1),',',num2str(j-1),',',num2str(g_ex(i,j)),',',num2str(object.EsynEx),', new ',excActivation{i,j}.getCbuilder(),');\n'];
                ii = ii + 1;
            end
        end
    end
    
    str = [str,'\n'];
    
    ii = 0;
    
    for i=1:N
        for j=1:N
            if g_el(i,j) ~= 0
                str = [str,'elSyn[',num2str(ii),'] = new CPG::synStruct(',num2str(i-1),',',num2str(j-1),',',num2str(g_el(i,j)),',0, new FTM_synapse_model());\n'];
                ii = ii + 1;
            end
        end
    end
    
    str = [str,'\n\n'];
    
    
    str2 = 'double delays[] ={';
    
    if numel(object.delays > 0)
    for i=1:numel(object.delays)-1
        str2 = [str2,num2str(object.delays(i)),','];
    end
    str2 = [str2,num2str(object.delays(end))];
    end
    
    
    
    
    
    str = [str,str2,'};\n\n'];
    str = [str,sprintf('*vf = new CPG(%d,neuron,%d,%d,%d,inhSyn,excSyn,elSyn,%d,delays)', ...
        object.N,Ninh,Nexc,Nel,numel(object.delays))];
    
    
    
    
end