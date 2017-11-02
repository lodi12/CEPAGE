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
fprintf(fout,'\tconst double gCa = %f;\n',object.get_gCa);
fprintf(fout,'\tconst double I = %f;\n',object.get_I);
fprintf(fout,'\tconst double V3 = %f;\n',object.get_V3);
fprintf(fout,'\tconst double V4 = %f;\n',object.get_V4);
fprintf(fout,'\tconst double phi = %f;\n',object.get_phi);

fprintf(fout,'\tconst double CM = %f;\n',object.get_CM);
fprintf(fout,'\tconst double gK = 8;\n');
fprintf(fout,'\tconst double gl = %f;\n',object.get_gL);
fprintf(fout,'\tconst double VCa = 120;\n');
fprintf(fout,'\tconst double Vk = -80;\n');
fprintf(fout,'\tconst double Vl = %f;\n',object.get_EL);
fprintf(fout,'\tconst double V1 = -1.2;\n');
fprintf(fout,'\tconst double V2 = 18;\n');
fprintf(fout,'\n');
fprintf(fout,'\tfloat Minf = 0.5*(1 + tanh((x[0] - V1 )/V2 ));\n');
fprintf(fout,'\tfloat Ninf = 0.5*(1 + tanh((x[0] - V3 )/V4 ));\n');
fprintf(fout,'\tfloat lamda_N = phi*cosh((x[0] - V3 )/(2*V4) );\n');
fprintf(fout,'\n');
fprintf(fout,'\txdot[1] = lamda_N*(Ninf-x[1]);\n');
fprintf(fout,'\txdot[0] = (-gl*(x[0] - Vl) - gCa*Minf*(x[0] - VCa) - gK*x[1]*(x[0] - Vk ) + I + Iext)/CM;\n');
fprintf(fout,'}\n\n');

fprintf(fout,['void ',gpustring,fileName,'_jac(double *x,double **jac,double Iext)\n']);
fprintf(fout,'{\n');
fprintf(fout,'\tconst double gCa = %f;\n',object.get_gCa);
fprintf(fout,'\tconst double I = %f;\n',object.get_I);
fprintf(fout,'\tconst double V3 = %f;\n',object.get_V3);
fprintf(fout,'\tconst double V4 = %f;\n',object.get_V4);
fprintf(fout,'\tconst double phi = %f;\n',object.get_phi);

fprintf(fout,'\tconst double CM = %f;\n',object.get_CM);
fprintf(fout,'\tconst double gK = 8;\n');
fprintf(fout,'\tconst double gl = %f;\n',object.get_gL);
fprintf(fout,'\tconst double VCa = 120;\n');
fprintf(fout,'\tconst double Vk = -80;\n');
fprintf(fout,'\tconst double Vl = %f;\n',object.get_EL);
fprintf(fout,'\tconst double V1 = -1.2;\n');
fprintf(fout,'\tconst double V2 = 18;\n');
fprintf(fout,'\n');
fprintf(fout,'\tdouble V = x[0];\n');
fprintf(fout,'\tdouble N = x[1];\n');
fprintf(fout,'\n');
fprintf(fout,'\tjac[0][0] = -(gl + N*gK + (gCa*(tanh((V - V1)/V2) + 1.0))/2.0 - (gCa*(V - VCa)*(pow(tanh((V - V1)/V2),2) - 1))/(2.0*V2))/CM;\n');
fprintf(fout,'\tjac[0][1] = -(gK*(V - Vk))/CM;\n');
fprintf(fout,'\tjac[1][0] = (phi*cosh((V - V3)/(2.0*V4))*(pow(tanh((V - V3)/V),2) - 1.0)*((V - V3)/(V*V) - 1.0/V))/2.0 + (phi*sinh((V - V3)/(2.0*V4))*(tanh((V - V3)/V)/2.0 - N + 1/2.0))/(2.0*V4);\n');
fprintf(fout,'\tjac[1][1] = -phi*cosh((V - V3)/(2.0*V4));\n');
fprintf(fout,'}');

fclose(fout);


end