function str = getCbuilder(object)

N = object.N;
neurons = object.neurons;

inhActivation = object.inhActivation;
excActivation = object.excActivation;

g_in = object.g_in;
g_ex = object.g_ex;
g_el = object.g_el;

str = [];

%     % g_in
%     str = [str,'double g_in[] ={'];
%     for i=1:N-1
%         for j=1:N-1
%         str = [str,sprintf('%e,',g_in(i,j))];
%         end
%         str = [str,sprintf('%e,\n',g_in(i,end))];
%     end
%     for j=1:N-1
%         str = [str,sprintf('%e,',g_in(end,j))];
%     end
%     str = [str,sprintf('%e};\n\n',g_in(end,end))];
%
%
%     % g_ex
%     str = [str,'double g_ex[] ={'];
%     for i=1:N-1
%         for j=1:N-1
%         str = [str,sprintf('%e,',g_ex(i,j))];
%         end
%         str = [str,sprintf('%e,\n',g_ex(i,end))];
%     end
%     for j=1:N-1
%         str = [str,sprintf('%e,',g_ex(end,j))];
%     end
%     str = [str,sprintf('%e};\n\n',g_ex(end,end))];



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
str = [str,'t_SynStruct **inhSynapses;\n'];
str = [str,'t_SynStruct **excSynapses;\n\n'];

str = [str,sprintf('neuron = (neuron_model **)malloc(%d*sizeof(neuron_model *));\n',N)];
%     str = [str,sprintf('inhSynapses = (synapse_model **)malloc(%d*%d*sizeof(synapse_model *));\n',N,N)];
%     str = [str,sprintf('excSynapses = (synapse_model **)malloc(%d*%d*sizeof(synapse_model *));\n\n',N,N)];

str = [str,sprintf('inhSynapses = new SynStruct*[%d];\n',sum(sum(g_in ~= 0)))];
str = [str,sprintf('excSynapses = new SynStruct*[%d];\n\n',sum(sum(g_ex ~= 0)))];

for i=1:N
    str = [str,'neuron[',num2str(i-1),'] = new ',neurons{i}.getCbuilder(),';\n'];
end

str = [str,'\n'];


indexSyn = 0;

for i=1:N
    for j=1:N
        if g_in(i,j) ~= 0
            indexSyn = indexSyn+1;
            str = [str,'inhSynapses[',num2str(indexSyn-1),'] = new SynStruct(',num2str(i-1),',',num2str(j-1),', new ',inhActivation{i,j}.getCbuilder(),',',num2str(g_in(i,j)),');\n'];
        end
    end
end

str = [str,'\n'];

for i=1:N
    for j=1:N
        if g_ex(i,j) ~= 0
            indexSyn = indexSyn+1;
            str = [str,'excSynapses[',num2str(indexSyn-1),'] = new SynStruct(',num2str(i-1),',',num2str(j-1),', new ',excActivation{i,j}.getCbuilder(),',',num2str(g_ex(i,j)),');\n'];
        end
    end
end

str = [str,'\n'];



str = [str,'\n\n'];

str = [str,sprintf('*vf = new CPG(%d,neuron, g_el,%e,%e,inhSynapses,excSynapses )', ...
    object.N,object.EsynIn,object.EsynEx)];




end