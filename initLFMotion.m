function initLFMotion
%INITLFMOTION Initialize the path settings for the LFMotion library
fprintf('Setting up the paths and global variables for LFMotion...');
p = mfilename('fullpath');
[root, ~, ~] = fileparts(p);
addpath(root);
subdirs = {'LFSceneFlow', 'utils'};
for i = 1:length(subdirs)
    addpath(genpath(fullfile(root,subdirs{i})));
end

% This assumes the light field dataset is placed at ./data.
% Modify this if you want to put them elsewhere.
global LFTopDir
LFTopDir = fullfile(root, 'data');

% If you want to generate depth maps by yourself,
% replace this with your cocolib root path
global CocolibDir
CocolibDir = '/home/sizhuo/Documents/cocolib-release-6';

global DEBUG
DEBUG = false;

fprintf('Done.\n');

end

