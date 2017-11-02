function writeContinuationInterface(object,limitCycle,varargin)

% writeContinuationInterface Generates files for the continuation of the
% limit cycle described by the network
%
% writeContinuationInterface(OBJ,limitCycle)
% Generates the  filename.m file for the continuation of the
% limit cycle limitCycle. limitCycle must be a struct with fields T and X.
% writeContinuationInterface(OBJ,limitCycle,opt )
% Generates files for the computation of the continuation of the
% limit cycle limitCycle. limitCycle must be a struct with fields T and X.
% opt  is a struct with the following fields:
% - code: a string describing the code language to be generated. It can be
%  'auto' (c+auto07 code) or 'matcont'  ('m'+matcont code).
% - filename : name of the function and of the file that describe the
% vectorial field
%
% Contributors:
%
% Matteo Lodi (matteo.lodi@edu.unige.it)
%
% Copyright (C) 2016 University of Genoa, Italy.

% Legal note:
% This program is free software; you can redistribute it and\or
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

filename = 'continuationInterface';
code = 'matcont';

if nargin == 3
    opt = varargin{1};
    if ~isstruct(opt)
        error('opt must be a struct');
    else
        if isfield(opt,'code')
            code = opt.code;
        end
        if isfield(opt,'filename')
            filename = opt.filename;
        end
    end
end


if strcmp(code,'matcont')
    warning('Only auto code can be written');
elseif strcmp(code,'auto')
    
    %%%%%%% AUTO INTERFACE FILE  %%%%%%%%
    
    % save limitCycle
    fin = fopen('limitCycle.dat','w');
    for i=1:numel(limitCycle.T)-1
        fprintf(fin,'%f\t',limitCycle.T(i));
        for j=1:size(limitCycle.X,2)-1
            fprintf(fin,'%f\t',limitCycle.X(i,j));
        end
        fprintf(fin,'%f\n',limitCycle.X(i,end));
    end
    fprintf(fin,'%f\t',limitCycle.T(end));
    for j=1:size(limitCycle.X,2)-1
        fprintf(fin,'%f\t',limitCycle.X(end,j));
    end
    fprintf(fin,'%f',limitCycle.X(end,end));
    fclose(fin);
    
    % Write .c file
    N = object.N;
    % First generate neuron Cfile
    for i=1:N
        nm = ['neuron',num2str(i)];
%         if isempty(folder)
            object.neurons{i}.generateC(nm);
%         else
%             object.neurons{i}.generateC(nm,folder);
%         end
    end
    
    fin = fopen([filename,'.c'],'w');
    %fprintf(fin,'#include "CPG.h"\n');
    fprintf(fin,'#include "auto_f2c.h"\n\n');
    for i=1:N
        fprintf(fin,'#include "neuron%d.h"\n',i);
    end
    fprintf(fin,'#define N %d\n',object.N);
    fprintf(fin,'#define nx %d\n',object.neurons{1}.getnx);
    
    fprintf(fin,'int func(integer ndim, const doublereal *u, const integer *icp,\n');
    fprintf(fin,'\tconst doublereal *par, integer ijac, doublereal *f, doublereal *dfdu, \n');
    fprintf(fin,'\tdoublereal *dfdp)\n');
    fprintf(fin,'{\n');
    
    
    nx = object.neurons{1}.getnx();
    g_in = object.g_in;
    g_ex = object.g_ex;
    g_el = object.g_el;
    
    
    
    fprintf(fin,'\t/**** Write here your continuation problem ****/\n\n');
    fprintf(fin,'\tconst double EsynIn = %f;\n',object.EsynIn);
    fprintf(fin,'\tconst double EsynEx = %f;\n',object.EsynEx);
    fprintf(fin,'\tconst double tetaSyn = %f;\n',object.tetaSyn);
    fprintf(fin,'\tconst double uSyn = %f;\n',object.uSyn);
    
    fprintf(fin,'\tconst double g_in[N][N] = {');
    for i=1:N
        fprintf(fin,'\t\t{');
        for j=1:N-1
            fprintf(fin,'%f,',g_in(i,j));
        end
        fprintf(fin,'%f}',g_in(i,end));
        if i ~= N
            fprintf(fin,',\n');
        end
    end
    fprintf(fin,'};\n');
    
    fprintf(fin,'\tconst double g_ex[N][N] = {');
    for i=1:N
        fprintf(fin,'\t\t{');
        for j=1:N-1
            fprintf(fin,'%f,',g_ex(i,j));
        end
        fprintf(fin,'%f}',g_ex(i,end));
        if i ~= N
            fprintf(fin,',\n');
        end
    end
    fprintf(fin,'};\n');
    
    fprintf(fin,'\tconst double g_el[N][N] = {');
    for i=1:N
        fprintf(fin,'\t\t{');
        for j=1:N-1
            fprintf(fin,'%f,',g_el(i,j));
        end
        fprintf(fin,'%f}',g_el(i,end));
        if i ~= N
            fprintf(fin,',\n');
        end
    end
    fprintf(fin,'};\n');
    
    fprintf(fin,'\tint i,j;\n');
    fprintf(fin,'\tdouble *xdot,*x;\n');
    fprintf(fin,'\tdouble Isyn[N];\n');
    
    fprintf(fin,'\tx = malloc(N*nx*sizeof(double));\n');
    fprintf(fin,'\txdot = malloc(N*nx*sizeof(double));\n');
    fprintf(fin,'\tfor(i=0;i<N*nx;i++)\n\t\tx[i] = u[i];\n');
    fprintf(fin,'\n');
    
    fprintf(fin,'\tfor(i=0;i<N;i++)\n');
    fprintf(fin,'\t{\n');
    fprintf(fin,'\t\tIsyn[i] = 0;\n');
    fprintf(fin,'\t\tfor(j=0;j<N;j++)\n');
    fprintf(fin,['\t\t\tIsyn[i] = Isyn[i] + g_ex[i][j]*(EsynEx-x[i*nx])/(1+exp(-uSyn*(x[j*nx]-tetaSyn)))+'...
        'g_in[i][j]*(EsynIn-x[i*nx])/(1+exp(-uSyn*(x[j*nx]-tetaSyn))) +'...
        'g_el[i][j]*(x[j*nx]-x[i*nx]);\n']);
    fprintf(fin,'\t}\n');
    
    for i=1:N
        fprintf(fin,'neuron%d(x+%d,xdot+%d,Isyn[%d]);\n',i,(i-1)*nx,(i-1)*nx,i-1);
    end    

    fprintf(fin,'\tfor(i=0;i<N*nx;i++)\n\t\tf[i] = xdot[i];\n');
    fprintf(fin,'\treturn 0;\n}\n\n');
    
    
    
    fprintf(fin,'int stpnt(integer ndim, doublereal t, doublereal *u, doublereal *par)\n');
    fprintf(fin,'{\n');
    fprintf(fin,'return 0;\n}\n');
    
    
    fprintf(fin,'int pvls(integer ndim, const doublereal *u, doublereal *par)\n');
    fprintf(fin,'{\n');
    fprintf(fin,'\treturn 0;\n');
    fprintf(fin,'}\n\n');
    
    
    fprintf(fin,'int bcnd(integer ndim, const doublereal *par, const integer *icp,\n');
    fprintf(fin,'\tinteger nbc, const doublereal *u0, const doublereal *u1, integer ijac,\n');
    fprintf(fin,'\tdoublereal *fb, doublereal *dbc)\n');
    fprintf(fin,'{\n');
    fprintf(fin,'\treturn 0;\n');
    fprintf(fin,'}\n\n');
    
    
    fprintf(fin,'int icnd(integer ndim, const doublereal *par, const integer *icp,\n');
    fprintf(fin,'\tinteger nint, const doublereal *u, const doublereal *uold,\n');
    fprintf(fin,'\tconst doublereal *udot, const doublereal *upold, integer ijac,\n');
    fprintf(fin,'\tdoublereal *fi, doublereal *dint)\n');
    fprintf(fin,'{\n');
    fprintf(fin,'\treturn 0;\n');
    fprintf(fin,'}\n\n\n');
    
    
    fprintf(fin,'int fopt(integer ndim, const doublereal *u, const integer *icp,\n');
    fprintf(fin,'\tconst doublereal *par, integer ijac,\n');
    fprintf(fin,'\tdoublereal *fs, doublereal *dfdu, doublereal *dfdp)\n');
    fprintf(fin,'{\n');
    fprintf(fin,'\treturn 0;\n');
    fprintf(fin,'}\n');
    
    totDim = 0;
    for i=1:N
        totDim = totDim + object.neurons{i}.getnx();
    end
    
    fin2 = fopen(['c.',filename],'w');
    fprintf(fin2,'unames = {');
    ii = 1;
    for i=1:N-1
        names = object.neurons{i}.getStateNames;
        for j=1:numel(names)
            fprintf(fin2,[num2str(ii),' : ''',names{j},num2str(i),''',']);
            ii = ii+1;
        end
    end
    names = object.neurons{end}.getStateNames;
    for j=1:numel(names)-1
        fprintf(fin2,[num2str(ii),' : ''',names{j},num2str(N),''',']);
        ii = ii+1;
    end
    fprintf(fin2,[num2str(ii),' : ''',names{end},num2str(N),'''}\n']);
    
    fprintf(fin2,'parnames = {}\n');
    fprintf(fin2,'dat = ''limitCycle''\n');
    fprintf(fin2,'NDIM = %d, IPS = 2 , IRS = 0, ILP = 0\n',totDim);
    fprintf(fin2,'ICP =  [1 11] \n');
    fprintf(fin2,'NTST=  700, NCOL=   5, IAD =   3, ISP =   2, ISW = 1, IPLT= 0, NBC= 0, NINT= 0\n');
    fprintf(fin2,'NMX=  7100, NPR=  0, MXBF=  1, IID =   2, ITMX= 8, ITNW= 7, NWTN= 5, JAC= 0\n');
    fprintf(fin2,'EPSL= 1e-05, EPSU = 1e-04, EPSS = 1e-04\n');
    fprintf(fin2,'DS  =   0.01, DSMIN= 0.00001, DSMAX=   0.002, IADS=   1\n');
    fprintf(fin2,'NPAR=   2, THL =  {11: 0.0}, THU =  {}\n');
    
    fclose(fin2);
    
    fin2 = fopen([filename,'.auto'],'w');
    fprintf(fin2,'import numpy as np\n');
    fprintf(fin2,'import scipy.io as sio\n');
    
    fprintf(fin2,'#***Compute a solution family***\n');
    fprintf(fin2,['cycleCont=run(e=''',filename,''',c=''',filename,''')\n\n']);
    
    fclose(fin2);
    
    msgbox(['Warning! In order to compute continuation you must set continuation problem and ',...
        'initial condition in file ',filename,'.c and set the continuation constant in ',...
        'c.',filename,' file.']);
    
    
    fclose(fin);
    
    
else
    error('Choise for code must be matcont or auto')
end

end
