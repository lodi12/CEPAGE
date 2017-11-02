function generateC(object,fileName,varargin)
% generateC Generates C files for the computation of model vector field
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

if nargin == 3
    folder = [varargin{1},'/'];
    if(exist(varargin{1},'dir') ~= 7)
        mkdir(folder);
    end
else
    folder = '';
end
    
fout = fopen([folder,fileName,'.h'],'w+');
fprintf(fout,['#ifndef ',fileName,'_H\n']);
fprintf(fout,['#define ',fileName,'_H\n\n']);
fprintf(fout,'#include <math.h>\n');
fprintf(fout,['void ',fileName,'(double *x,double *xdot,double Iext);\n\n']);
fprintf(fout,['void ',fileName,'_jac(double *x,double **jac,double Iext);\n\n']);
fprintf(fout,'#endif');
fclose(fout);

fout = fopen([folder,fileName,'.c'],'w+');
fprintf(fout,['#include "',[fileName,'.h'],'"\n\n']);
fprintf(fout,['void ',fileName,'(double *x,double *xdot,double Iext)\n']);
fprintf(fout,'{\n');
fprintf(fout,'\tconst double eps = %f;\n',object.get_eps);
fprintf(fout,'\tconst double I = %f;\n',object.get_I);
fprintf(fout,'\n');
fprintf(fout,'\txdot[0] = x[0]-x[0]*x[0]*x[0]+I-x[1]+Iext;\n');
fprintf(fout,'\txdot[1] = eps*(1/(1+exp(-10*x[0]))-x[1]);\n');
fprintf(fout,'}\n\n');

fprintf(fout,['void ',fileName,'_jac(double *x,double **jac,double Iext)\n']);
fprintf(fout,'{\n');
fprintf(fout,'\tconst double eps = %f;\n',object.get_eps);
fprintf(fout,'\tconst double I = %f;\n',object.get_I);
fprintf(fout,'\n');

fprintf(fout,'\tjac[0][0] = 1.0-3.0*x[0]*x[0];\n');
fprintf(fout,'\tjac[0][1] = -1.0;\n');
fprintf(fout,'\tjac[1][0] = eps*10.0*exp(-10.0*x[0])/(pow((1+exp(-10.0*x[0])),2));\n');
fprintf(fout,'\tjac[1][1] = -eps;\n');
fprintf(fout,'}');

fclose(fout);


end