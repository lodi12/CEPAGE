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
if nargin >= 3
    folder = [varargin{1},'/'];
    if(exist(varargin{1},'dir') ~= 7)
        mkdir(folder);
    end
else
    folder = '';
end
    
gpustring = '';
if nargin == 4
    if(varargin{2})
        gpustring = ' __device__ ';
    end
end

fout = fopen([folder,fileName,'.h'],'w+');
fprintf(fout,['#ifndef ',fileName,'_H\n']);
fprintf(fout,['#define ',fileName,'_H\n\n']);
fprintf(fout,'#include <math.h>\n');
fprintf(fout,['void ',gpustring,fileName,'(double *x,double *xdot,double Iext);\n\n']);
fprintf(fout,['void ',gpustring,fileName,'_jac(double *x,double **jac,double Iext);\n\n']);
fprintf(fout,'#endif');
fclose(fout);

fout = fopen([folder,fileName,'.c'],'w+');
fprintf(fout,['#include "',[fileName,'.h'],'"\n\n']);
fprintf(fout,['void ',gpustring,fileName,'(double *x,double *xdot,double Iext)\n']);
fprintf(fout,'{\n');
fprintf(fout,'\tconst double C = %f;\n',object.C);
fprintf(fout,'\tconst double gNaP = %f;\n',object.gNaP);
fprintf(fout,'\tconst double ENa =  %f;\n',object.ENa);
fprintf(fout,'\tconst double gl = %f;\n',object.gl);
fprintf(fout,'\tconst double El = %f;\n',object.El);
fprintf(fout,'\tconst double VhalfM = %f;\n',object.VhalfM);
fprintf(fout,'\tconst double km = %f;\n',object.km);
fprintf(fout,'\tconst double Vhalfh = %f;\n',object.Vhalfh);
fprintf(fout,'\tconst double kh = %f;\n',object.kh);
fprintf(fout,'\tconst double tau0 = %f;\n',object.tau0);
fprintf(fout,'\tconst double tauMax = %f;\n',object.tauMax);
fprintf(fout,'\tconst double VhalfTau = %f;\n',object.VhalfTau);
fprintf(fout,'\tconst double kTau = %f;\n',object.kTau);
fprintf(fout,'\tconst double Di = %f;\n',object.Di);
fprintf(fout,'\tconst double gSynE = %f;\n',object.gSynE);
fprintf(fout,'\tconst double EsynE = %f;\n',object.EsynE);


fprintf(fout,'\n');

fprintf(fout,'\tdouble V = x[0];\n');
fprintf(fout,'\tdouble h = x[1];\n');

fprintf(fout,'\tdouble m = 1.0/(1.0+exp((V-VhalfM)/km));\n');

fprintf(fout,'\tdouble tauH = tau0 +((tauMax-tau0)/cosh((V-VhalfTau)/kTau));\n');

fprintf(fout,'\tdouble hInf = 1.0/(1+exp((V-Vhalfh)/kh));\n');
fprintf(fout,'\tdouble Idrive = gSynE*(V-EsynE)*Di;\n');
fprintf(fout,'\txdot[1] = 1.0/tauH*(hInf-h);\n');
fprintf(fout,'\txdot[0] = 1.0/C*(-gl*(V-El)-gNaP*m*h*(V-ENa)-Idrive+Iext);\n\n}\n\n');


fprintf(fout,['void ',gpustring,fileName,'_jac(double *x,double **jac,double Iext)\n']);
fprintf(fout,'{\n');
% fprintf(fout,'\tconst double b = %f;\n',object.get_b);
% fprintf(fout,'\tconst double I = %f;\n',object.get_I);
% fprintf(fout,'\tconst double s = %f;\n',object.get_s);
% fprintf(fout,'\tconst double mu = %f;\n',object.get_mu);
% fprintf(fout,'\tconst double x_rest = %f;\n',object.get_x_rest);
% fprintf(fout,'\n');
% 
% fprintf(fout,'\tjac[0][0] = -3.0*x[0]*x[0]+2.0*b*x[0];\n');
% fprintf(fout,'\tjac[0][1] = 1.0;\n');
% fprintf(fout,'\tjac[0][2] = -1.0;\n');
% fprintf(fout,'\tjac[1][0] = -10.0*x[0];\n');
% fprintf(fout,'\tjac[1][1] = -1.0;\n');
% fprintf(fout,'\tjac[1][2] = 0;\n');
% fprintf(fout,'\tjac[2][0] = mu*s;\n');
% fprintf(fout,'\tjac[2][1] = 0;\n');
% fprintf(fout,'\tjac[2][2] = -mu;\n');
fprintf(fout,'}');
fclose(fout);


end