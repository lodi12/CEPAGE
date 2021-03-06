disp(' ')
disp('********************************')
disp('* Welcome to CEPAGE Toolbox! *')
disp('********************************')
disp(' ')
disp('Checking mex compiler installation:')

warning off

myCCompiler = mex.getCompilerConfigurations('C','Selected');

if isempty(myCCompiler) || isempty(myCCompiler.Name)
    error('Unable to found C compiler for mex file')
else
    disp(['C compiler found: ',myCCompiler.Name,' compiler']);
end

myCppCompiler = mex.getCompilerConfigurations('C++','Selected');

if isempty(myCppCompiler) || isempty(myCppCompiler.Name)
    error('Unable to found C++ compiler for mex file')
else
    disp(['C++ compiler found: ',myCppCompiler.Name,' compiler']);
end
disp(' ')
boostDir = input('Insert boost installation directory\n','s');

if boostDir(end) == '/' && numel(boostDir) ~= 1
    boostDir = boostDir(1:end-1);
end

fin = fopen('function/getCEPAGEPar.m','w');
fprintf(fin,'function CEPAGEPar = getCEPAGEPar()\n');
fprintf(fin,['CEPAGEPar.boostDir = ''',boostDir,''';']);
fclose(fin);

disp(' ')
disp('Checking odeint installation...')


fout = fopen('try.cpp','w');
fprintf(fout,['#include <iostream>\n',...
    '#include <boost/serialization/array_wrapper.hpp>',...
    '\n',...
    '#include <boost/array.hpp>\n',...
    '\n',...
    '#include <boost/numeric/odeint.hpp>\n',...
    '\n',...
    'using namespace std;\n',...
    'using namespace boost::numeric::odeint;\n',...
    '\n',...
    'const double sigma = 10.0;\n',...
    'const double R = 28.0;\n',...
    'const double b = 8.0 / 3.0;\n',...
    '\n',...
    'typedef boost::array< double , 3 > state_type;\n',...
    '\n',...
    'void lorenz( const state_type &x , state_type &dxdt , double t )\n',...
    '{\n',...
    '    dxdt[0] = sigma * ( x[1] - x[0] );\n',...
    '    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];\n',...
    '    dxdt[2] = -b * x[2] + x[0] * x[1];\n',...
    '}\n',...
    '\n',...
    'int main(int argc, char **argv)\n',...
    '{\n',...
    '    state_type x = { 10.0 , 1.0 , 1.0 }; // initial conditions\n',...
    '    integrate( lorenz , x , 0.0 , 25.0 , 0.1  );\n',...
    '}\n',...
    '\n']);
fclose(fout);

isOk = true;

try
    eval(['mex -silent -c try.cpp -I"',boostDir,'/include" -L"',boostDir,'/lib"'])
catch
    isOk = false;
end

delete('try.*')

if ~isOk
    error('Unable to find odeint library');
else
    disp('Odeint library found');
end


disp(' ');
disp('Compiling mex files...')



% First compile libraries
cd c_file

mkdir bin;

infoFile = dir('src');

totalString = [];

for i=1:numel(infoFile)
    tmpFile = infoFile(i).name;
    if numel(tmpFile) > 4
        if strcmp(tmpFile(end-3:end),'.cpp')
            eval(['mex -silent  -c src/',tmpFile,' -outdir bin']);
            totalString = [totalString,' bin/',tmpFile(1:end-3),'o'];
        end
    end
    
end

if isunix
    status = system(['ar rcs libCEPAGE.a ',totalString]);
elseif ispc
    status = system(['ar rcs libCEPAGE.lib ',totalString]);
else
    error('Unsopported operative system');
end
% Then compiel C file
% boostDir = getCEPAGEPar();
% boostDir = boostDir.boostDir;
mex -silent  -c eulero.cpp -L. -lCEPAGE 
disp('eulero.cpp compiled')
mex -silent -c euleroEvents.cpp -L. -lCEPAGE 
disp('euleroEvents.cpp compiled')
eval(['mex -c odeint.cpp -silent -L. -lCEPAGE -I"',boostDir,'/include" -L"',boostDir,'/lib"']);
disp('odeint.cpp compiled')
eval(['mex -c odeintEvents.cpp -silent -L. -lCEPAGE -I"',boostDir,'/include" -L"',boostDir,'/lib"']);
disp('odeintEvents.cpp compiled')
warning on

cd ..

disp(' ')
disp('Done.')
disp(' ')

disp(' ')
disp('Adding folders to path...')
disp(' ')

addpath([pwd,'/classes'])
disp([pwd,'/classes'])
addpath([pwd,'/c_file'])
disp([pwd,'/c_file'])
addpath(genpath([pwd,'/function']))
disp([pwd,'/function'])

disp(' ')
disp('Done.')
disp(' ')


disp(' ')
disp('Installation complete.');




