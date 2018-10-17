function flow = vsf( jsonFile0, jsonFile1, resultDir, customParam )
%VSF Run VSF on Sim data

global VSFDir
binFile = './StereoFlow';

if nargin < 4
    customParam = struct();
end

iniFile = writeVsfIni(jsonFile0, jsonFile1, resultDir, customParam);

currDir = pwd();
cd(VSFDir);
command = sprintf('%s %s',binFile,iniFile);
disp(command);
system(command,'-echo');
cd(currDir);

flow = convertVsf(jsonFile0, resultDir);

end

