function cocolibDisparity( lffile, customParam )
%COCOLIBDISPARITY Call Cocolib to estimate disparity from light fields

global CocolibDir
binFile = fullfile(CocolibDir, 'lightfields/lightfields');

[pathstr, fname, ~] = fileparts(lffile);
epiFile = fullfile(pathstr,['disparity_epi_' fname]);
filterFile = fullfile(pathstr,['disparity_epi_filter_' fname]);

if nargin < 2
    customParam = struct();
end

writeCocolibConfig(lffile, customParam);

command = sprintf('%s -config %s',binFile,epiFile);
disp(command);
system(command,'-echo');
command = sprintf('%s -config %s',binFile,filterFile);
disp(command);
system(command,'-echo');

delete(fullfile(pathstr,fname,'disparity_epi','*.h5'));

end

