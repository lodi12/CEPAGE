function str = getCbuilder(object)

    N = object.N;
    neurons = object.neurons;

    inhActivation = object.inhActivation;
excActivation = object.excActivation;
    
    g_in = object.g_in;
    g_ex = object.g_ex;
    g_el = object.g_el;
    
    str = [];
    
    % g_in
    str = [str,'double g_in[] ={'];
    for i=1:N-1
        for j=1:N-1
        str = [str,sprintf('%e,',g_in(i,j))];
        end
        str = [str,sprintf('%e,\n',g_in(i,end))]; 
    end   
    for j=1:N-1
        str = [str,sprintf('%e,',g_in(end,j))];
    end  
    str = [str,sprintf('%e};\n\n',g_in(end,end))];
    
    
    % g_ex
    str = [str,'double g_ex[] ={'];
    for i=1:N-1
        for j=1:N-1
        str = [str,sprintf('%e,',g_ex(i,j))];
        end
        str = [str,sprintf('%e,\n',g_ex(i,end))]; 
    end   
    for j=1:N-1
        str = [str,sprintf('%e,',g_ex(end,j))];
    end  
    str = [str,sprintf('%e};\n\n',g_ex(end,end))];
    
    
    
    % g_el
    str = [str,'double g_el[] ={'];
    for i=1:N-1
        for j=1:N-1
        str = [str,sprintf('%e,',g_el(i,j))];
        end
        str = [str,sprintf('%e,\n',g_el(i,end))]; 
    end   
    for j=1:N-1
        str = [str,sprintf('%e,',g_el(end,j))];
    end  
    str = [str,sprintf('%e};\n\n',g_el(end,end))];
    
    
    
    
    str = [str,'neuron_model **neuron;\n\n'];
    str = [str,'synapse_model **inhSynapses;\n'];
    str = [str,'synapse_model **excSynapses;\n\n'];
    
    str = [str,sprintf('neuron = (neuron_model **)malloc(%d*sizeof(neuron_model *));\n',N)];
    str = [str,sprintf('inhSynapses = (synapse_model **)malloc(%d*%d*sizeof(synapse_model *));\n',N,N)];
    str = [str,sprintf('excSynapses = (synapse_model **)malloc(%d*%d*sizeof(synapse_model *));\n\n',N,N)];
    
    for i=1:N
    str = [str,'neuron[',num2str(i-1),'] = new ',neurons{i}.getCbuilder(),';\n'];
    end
    
    str = [str,'\n'];
    
    for i=1:N
        for j=1:N
            indexSyn = (i-1)*N+(j-1)+1;
            str = [str,'inhSynapses[',num2str(indexSyn-1),'] = new ',inhActivation{i,j}.getCbuilder(),';\n'];
        end
    end
    
    str = [str,'\n'];
    
    for i=1:N
        for j=1:N
            indexSyn = (i-1)*N+(j-1)+1;
            str = [str,'excSynapses[',num2str(indexSyn-1),'] = new ',excActivation{i,j}.getCbuilder(),';\n'];
        end
    end
    
    str = [str,'\n'];
    
    
    
    str = [str,'\n\n'];
    
    
    str2 = 'double delays[] ={';
    
    if numel(object.delays > 0)
    for i=1:numel(object.delays)-1
        str2 = [str2,num2str(object.delays(i)),','];
    end
    str2 = [str2,num2str(object.delays(end))];
    end
    
    
    
    
    
    str = [str,str2,'};\n\n'];
    str = [str,sprintf('*vf = new CPG(%d,neuron,g_in, g_ex, g_el,%e,%e,inhSynapses,excSynapses,delays)', ...
        object.N,object.EsynIn,object.EsynEx)];
    
    
    
    
end