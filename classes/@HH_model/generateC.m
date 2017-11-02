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

% Define param
fprintf(fout,'\tconst double gna = %f;\n',object.gna);
fprintf(fout,'\tconst double ENa = %f;\n',object.ENa);
fprintf(fout,'\tconst double gk2 = %f;\n',object.gk2);
fprintf(fout,'\tconst double Ek = %f;\n',object.Ek);
fprintf(fout,'\tconst double gl = %f;\n',object.gl);
fprintf(fout,'\tconst double El = %f;\n',object.El);
fprintf(fout,'\tconst double tNa = %f;\n',object.tNa);
fprintf(fout,'\tconst double tk2 = %f;\n',object.tk2);
fprintf(fout,'\tconst double C = %f;\n',object.C);
fprintf(fout,'\tconst double Iapp = %f;\n',object.Iapp);
fprintf(fout,'\tconst double VshiftK2 = %f;\n',object.VshiftK2);
fprintf(fout,'\n');

fprintf(fout,'\t double hInf = 1/(1+exp(500*(x[0]+0.0333)));\n');
fprintf(fout,'\t double nInf = 1/(1+exp(-150*(x[0]+0.0305)));\n');
fprintf(fout,'\t double mInf = 1/(1+exp(-83*(x[0]+0.018+VshiftK2)));\n');

fprintf(fout,'\t double INa = gna*nInf*nInf*nInf*x[1]*(x[0]-ENa);\n');
fprintf(fout,'\t double Ik2 = gk2*x[2]*x[2]*(x[0]-Ek);\n');
fprintf(fout,'\t double Il = gl*(x[0]-El);\n');

fprintf(fout,'\txdot[0] = (-INa-Ik2-Il-Iapp+Iext)/C;\n');
fprintf(fout,'\txdot[1] = (hInf- x[1])/tNa;\n');
fprintf(fout,'\txdot[2] = (mInf- x[2])/tk2;\n');
fprintf(fout,'}\n\n');

fprintf(fout,['void ',fileName,'_jac(double *x,double **jac,double Iext)\n']);
fprintf(fout,'{\n');
fprintf(fout,'\tconst double gna = %f;\n',object.gna);
fprintf(fout,'\tconst double ENa = %f;\n',object.ENa);
fprintf(fout,'\tconst double gk2 = %f;\n',object.gk2);
fprintf(fout,'\tconst double Ek = %f;\n',object.Ek);
fprintf(fout,'\tconst double gl = %f;\n',object.gl);
fprintf(fout,'\tconst double El = %f;\n',object.El);
fprintf(fout,'\tconst double tNa = %f;\n',object.tNa);
fprintf(fout,'\tconst double tk2 = %f;\n',object.tk2);
fprintf(fout,'\tconst double C = %f;\n',object.C);
fprintf(fout,'\tconst double Iapp = %f;\n',object.Iapp);
fprintf(fout,'\tconst double VshiftK2 = %f;\n',object.VshiftK2);
fprintf(fout,'\n');
fprintf(fout,'\n');
fprintf(fout,'\tdouble V = x[0];\n');
fprintf(fout,'\tdouble mK2 = x[1];\n');
fprintf(fout,'\tdouble hNa = x[2];\n');
fprintf(fout,'\n');

fprintf(fout,'\tjac[0][0] = -(gl + gk2*mK2*mK2 + (gna*hNa)/pow((exp(- 150.0*V - 183.0/40.0) + 1.0),3) - (450.0*gna*hNa*exp(- 150.0*V - 183.0/40.0)*(ENa - V))/pow((exp(- 150.0*V - 183.0/40.0) + 1.0),4))/C;\n');
fprintf(fout,'\tjac[0][1] = (2.0*gk2*mK2*(Ek - V))/C;\n');
fprintf(fout,'\tjac[0][2] = (gna*(ENa - V))/(C*pow((exp(- 150.0*V - 183.0/40.0) + 1.0),3));\n');

fprintf(fout,'\tjac[1][0] =(83.0*exp(- 83.0*V - 83.0*VshiftK2 - 747.0/500.0))/(tk2*pow((exp(- 83.0*V - 83.0*VshiftK2 - 747.0/500.0) + 1.0),2));\n');
fprintf(fout,'\tjac[1][1] =-1.0/tk2;\n');
fprintf(fout,'\tjac[1][2] =0;\n');

fprintf(fout,'\tjac[2][0] =-(500.0*exp(500.0*V + 65.0/4.0))/(tNa*pow((exp(500.0*V + 65.0/4.0) + 1.0),2));\n');
fprintf(fout,'\tjac[2][1] =0;\n');
fprintf(fout,'\tjac[2][2] =-1.0/tNa;\n');

                                                                              
fprintf(fout,'}');

fclose(fout);


end