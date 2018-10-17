function [ list ] = listCamera( part )
%LISTCAMERA List all the available camera parts

global LFOptDir
files = dir(fullfile(LFOptDir,part,'*.json'));
list = cell(size(files));
for i = 1:numel(files)
    [~,fName,~] = fileparts(files(i).name);
    list{i} = fName;
end

end

