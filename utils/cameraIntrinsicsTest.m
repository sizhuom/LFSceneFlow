%% Unrectified
% Input filenames
LFTopDirectory = '/media/sizhuo/Seagate Backup Plus Drive/LFData';
% directory = fullfile(LFTopDirectory, 'Cameras/B5150805530/Cal26_100_33');
directory = fullfile(LFTopDirectory, 'Cameras/B5150805530/Cal51_inf_100');
filepath = fullfile(directory, 'CalInfo.json');

% Load the calibration info
CalInfo = LFReadMetadata(filepath);
EstCamIntrinsicsH = CalInfo.EstCamIntrinsicsH;

% Assemble the (i, j, k, l) coordinates
% The 4D structure of the LF array is organized as [j, i, l, k].
sz = [15 15 434 625];
[Y, X, V, U] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3), 1:sz(4));
NPixels = length(X(:));
XYUV = [ X(:)'; Y(:)'; U(:)'; V(:)'; ones(1, NPixels); ];

% Transform into (x, y, u, v) coordinates
% (s, t, u, v) in [Dansereau et al.]
XYUV = EstCamIntrinsicsH * XYUV;
X = reshape(XYUV(1, :), sz);
Y = reshape(XYUV(2, :), sz);
U = reshape(XYUV(3, :), sz);
V = reshape(XYUV(4, :), sz);

%% Rectified
% Input filenames
LFTopDirectory = '/media/sizhuo/Seagate Backup Plus Drive/LFData';
directory = fullfile(LFTopDirectory, 'Images/Motion/trans-160826-rectified/');
filenames = arrayfun(@(id) sprintf('IMG_%04d__Decoded.mat', id), [401], 'UniformOutput', false);

% Print the instrinsics
Hs = cell(numel(filenames),1);
for i = 1:numel(filenames)
    load(fullfile(directory,filenames{i}));
    Hs{i} = RectOptions.RectCamIntrinsicsH;
    disp(Hs{i});
end

%% Synthetic
H = genIntrinsics();