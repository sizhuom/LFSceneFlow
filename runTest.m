%% Run this script to test the performance of the CLG method
% choose the dataset to test among: {'card', 'threecards', 'hand', 
% 'card-pile', 'flash', 'desktop', 'plant', 'mug', 'pillow', 'shakehand',
% 'wavehand'}
dataset = 'wavehand';

% Subaperture image are already extracted and included in the dataset.
% If you want to work on you own dataset:
% Please use LightField_GeoCalibration_ver2 to calibrate your Lytro Illum
% camera and extract the subaperture images. Put them in a similar folder
% structure as the included dataset.

% Depth maps are already included in the dataset, but if you want to run
% Cocolib by yourself, uncomment the following line:
% preprocess(['data/' dataset]);

% Run CLG
runCLG(['data/' dataset '/param.json']);