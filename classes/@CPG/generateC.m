function generateC(object,fileName,varargin)
% generateC Generates C files for the computation of network vector field
%
% generateC(OBJ,filename)
% Generates the C filename.c and filename.h files for the computation of
% model vector field
%
% generateC(OBJ,filename,folder)
% Generates the C filename.c and filename.h files for the computation of
% model vector field in the directory folder
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
N = object.N;

if nargin == 3
    folder = [varargin{1},'/'];
    if(exist(varargin{1},'dir') ~= 7)
        mkdir(folder);
    end
else
    folder = '';
end

% First generate neuron Cfile
for i=1:N
    nm = ['neuron',num2str(i)];
    if isempty(folder)
    object.neurons{i}.generateC(nm);
    else
       object.neurons{i}.generateC(nm,folder); 
    end
end

nx = object.neurons{1}.getnx();
g_in = object.g_in;
g_ex = object.g_ex;
g_el = object.g_el;

fout = fopen([folder,fileName,'.h'],'w+');
fprintf(fout,['#ifndef ',fileName,'_H\n']);
fprintf(fout,['#define ',fileName,'_H\n\n']);
fprintf(fout,'#include <math.h>\n');
fprintf(fout,'#include <stdlib.h>\n');
for i=1:N
    fprintf(fout,'#include "neuron%d.h"\n',i);
end
fprintf(fout,'#define N %d\n',object.N);
fprintf(fout,'#define nx %d\n',nx);
fprintf(fout,['void ',fileName,'(double *x,double *xdot,double Iext);\n\n']);
fprintf(fout,['void ',fileName,'_jac(double *x,double **jac,double Iext);\n\n']);
fprintf(fout,'#endif');
fclose(fout);

fout = fopen([folder,fileName,'.c'],'w+');
fprintf(fout,['#include "',[fileName,'.h'],'"\n\n']);
fprintf(fout,['void ',fileName,'(double *x,double *xdot,double Iext)\n']);
fprintf(fout,'{\n');
fprintf(fout,'\tconst double EsynIn = %f;\n',object.EsynIn);
fprintf(fout,'\tconst double EsynEx = %f;\n',object.EsynEx);
fprintf(fout,'\tconst double tetaSyn = %f;\n',object.tetaSyn);
fprintf(fout,'\tconst double uSyn = %f;\n',object.uSyn);

fprintf(fout,'\tconst double g_in[N][N] = {');
for i=1:N
    fprintf(fout,'\t\t{');
    for j=1:N-1
        fprintf(fout,'%f,',g_in(i,j));
    end
    fprintf(fout,'%f}',g_in(i,end));
    if i ~= N
        fprintf(fout,',\n');
    end
end
fprintf(fout,'};\n');

fprintf(fout,'\tconst double g_ex[N][N] = {');
for i=1:N
    fprintf(fout,'\t\t{');
    for j=1:N-1
        fprintf(fout,'%f,',g_ex(i,j));
    end
    fprintf(fout,'%f}',g_ex(i,end));
    if i ~= N
        fprintf(fout,',\n');
    end
end
fprintf(fout,'};\n');

fprintf(fout,'\tconst double g_el[N][N] = {');
for i=1:N
    fprintf(fout,'\t\t{');
    for j=1:N-1
        fprintf(fout,'%f,',g_el(i,j));
    end
    fprintf(fout,'%f}',g_el(i,end));
    if i ~= N
        fprintf(fout,',\n');
    end
end
fprintf(fout,'};\n');

fprintf(fout,'\tint i,j;\n');
fprintf(fout,'\tdouble Isyn[N];\n');

fprintf(fout,'\n');

fprintf(fout,'\tfor(i=0;i<N;i++)\n');
fprintf(fout,'\t{\n');
fprintf(fout,'\t\tIsyn[i] = 0;\n');
fprintf(fout,'\t\tfor(j=0;j<N;j++)\n');
fprintf(fout,['\t\t\tIsyn[i] = Isyn[i] + g_ex[i][j]*(EsynEx-x[i*nx])/(1+exp(-uSyn*(x[j*nx]-tetaSyn)))+'...
    'g_in[i][j]*(EsynIn-x[i*nx])/(1+exp(-uSyn*(x[j*nx]-tetaSyn))) +'...
    'g_el[i][j]*(x[j*nx]-x[i*nx]);\n']);
fprintf(fout,'\t}\n');

for i=1:N
    fprintf(fout,'neuron%d(x+%d,xdot+%d,Isyn[%d]);\n',i,(i-1)*nx,(i-1)*nx,i-1);
end
fprintf(fout,'}\n\n');

fprintf(fout,['void ',fileName,'_jac(double *x,double **jac,double Iext)']);

fprintf(fout,'{\n');
fprintf(fout,'\tconst double Esyn_In = %f;\n',object.EsynIn);
fprintf(fout,'\tconst double Esyn_Ex = %f;\n',object.EsynEx);
fprintf(fout,'\tconst double tetaSyn = %f;\n',object.tetaSyn);
fprintf(fout,'\tconst double uSyn = %f;\n',object.uSyn);

fprintf(fout,'\tconst double g_in[N][N] = {');
for i=1:N
    fprintf(fout,'\t\t{');
    for j=1:N-1
        fprintf(fout,'%f,',g_in(i,j));
    end
    fprintf(fout,'%f}',g_in(i,end));
    if i ~= N
        fprintf(fout,',\n');
    end
end
fprintf(fout,'};\n');

fprintf(fout,'\tconst double g_ex[N][N] = {');
for i=1:N
    fprintf(fout,'\t\t{');
    for j=1:N-1
        fprintf(fout,'%f,',g_ex(i,j));
    end
    fprintf(fout,'%f}',g_ex(i,end));
    if i ~= N
        fprintf(fout,',\n');
    end
end
fprintf(fout,'};\n');

fprintf(fout,'\tconst double g_el[N][N] = {');
for i=1:N
    fprintf(fout,'\t\t{');
    for j=1:N-1
        fprintf(fout,'%f,',g_el(i,j));
    end
    fprintf(fout,'%f}',g_el(i,end));
    if i ~= N
        fprintf(fout,',\n');
    end
end
fprintf(fout,'};\n\n');

ndim_vector = zeros(N,1);

for i=1:N
    ndim_vector(i) = object.neurons{i}.getnx;
end

fprintf(fout,'\tconst int ndim_vector[N] = {');
for i=1:N-1
    fprintf(fout,'%d,',ndim_vector(i));
end
fprintf(fout,'%d};\n',ndim_vector(end));

fprintf(fout,'\tconst int ndim = %d;\n\nint i,j,k,currentPointer,tmpi;\n\n',sum(ndim_vector));

fprintf(fout,'double *xsingolo,**neuronJac;\n\n');

fprintf(fout,'// Jacobian function vector\n');
fprintf(fout,'void (*JacobianFunction[N])(double *, double **, double );\n');

% MATLAB
for i=1:N
    fprintf(fout,'JacobianFunction[%d] = neuron%d_jac;\n',i-1,i);
end
% C

%fprintf(fout,'\n// init jac\n');
%fprintf(fout,'jac = (double **)malloc(sizeof(double *)*ndim);\n');
%fprintf(fout,'for(i=0;i<ndim;i++)\n');
%fprintf(fout,'\tjac[i] = (double *)malloc(sizeof(double)*ndim);\n\n');

fprintf(fout,'currentPointer = 0;\n');
fprintf(fout,'for(k=0;k<N;k++)\n');
fprintf(fout,'{\n');
fprintf(fout,'\tint neurNx = ndim_vector[k];\n\n');
    
fprintf(fout,'\t// Fill k-th neuron state\n');
fprintf(fout,'\txsingolo = (double *)malloc(sizeof(double)*neurNx);\n');
fprintf(fout,'\tfor(i=0;i<neurNx;i++)\n');
fprintf(fout,'\t\txsingolo[i] = x[currentPointer+i];\n\n');
    
fprintf(fout,'\t// Fill k-th neuron Jacobian\n');
fprintf(fout,'\tneuronJac =(double **) malloc(sizeof(double *)*neurNx);\n');
fprintf(fout,'\tfor(i=0;i<neurNx;i++)\n');
fprintf(fout,'\t\tneuronJac[i] = (double *)malloc(sizeof(double)*neurNx);\n\n');
    
fprintf(fout,'\t// Compute Jac of single neuron\n');
fprintf(fout,'\t(*JacobianFunction[k])(xsingolo,neuronJac,0);\n\n');
    
fprintf(fout,'\t// And set in big jac matrix \n');
fprintf(fout,'\tfor(i=0;i<neurNx;i++)\n');
fprintf(fout,'\t\tfor(j=0;j<neurNx;j++)\n');
fprintf(fout,'\t\t\tjac[currentPointer+i][currentPointer+j] = neuronJac[i][j];\n\n');
    
fprintf(fout,'\t// Add dfi/dVi contribute\n');
fprintf(fout,'\ttmpi = 0;\n');
fprintf(fout,'\tfor(j=0;j<N;j++)\n');
fprintf(fout,'\t{\n');
fprintf(fout,'\t\tjac[currentPointer][currentPointer] -= g_in[k][j]/(1.0+exp(-uSyn*(x[tmpi]-tetaSyn)));\n');
fprintf(fout,'\t\tjac[currentPointer][currentPointer] -= g_ex[k][j]/(1.0+exp(-uSyn*(x[tmpi]-tetaSyn)));\n');
fprintf(fout,'\t\tjac[currentPointer][currentPointer] -= g_el[k][j];\n');
fprintf(fout,'\t\ttmpi += ndim_vector[j];\n');
fprintf(fout,'\t}\n\n');
        

fprintf(fout,'\t// Add dfi/dVj contribute\n');
fprintf(fout,'\ttmpi = 0;\n');
fprintf(fout,'\tfor(j=0;j<N;j++)\n');
fprintf(fout,'\t{\n');
fprintf(fout,'\t\tjac[currentPointer][tmpi] += g_in[k][j]*(Esyn_In-x[currentPointer])*(uSyn*exp(-uSyn*(x[tmpi]-tetaSyn)))/pow(1+exp(-uSyn*(x[tmpi]-tetaSyn)),2);\n');
fprintf(fout,'\t\tjac[currentPointer][tmpi] += g_ex[k][j]*(Esyn_Ex-x[currentPointer])*(uSyn*exp(-uSyn*(x[tmpi]-tetaSyn)))/pow(1+exp(-uSyn*(x[tmpi]-tetaSyn)),2);\n');
fprintf(fout,'\t\tjac[currentPointer][tmpi] += g_el[k][j];\n');
fprintf(fout,'\t\ttmpi += ndim_vector[j];\n');
fprintf(fout,'\t}\n\n');
           
    
fprintf(fout,'\tfree(xsingolo);\n');
fprintf(fout,'\tfor(i=0;i<neurNx;i++)\n');
fprintf(fout,'\t\tfree(neuronJac[i]);\n');
fprintf(fout,'\tfree(neuronJac);\n\n');
    
fprintf(fout,'\tcurrentPointer += neurNx;\n');
        
fprintf(fout,'}\n');




fprintf(fout,'}');


fclose(fout);


end