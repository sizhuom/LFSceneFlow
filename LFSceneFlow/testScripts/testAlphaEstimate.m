global LFTopDir;
workDir = fullfile(LFTopDir, 'Images/Sim/0505-objmotion');
subDir = '.';
% logFileName = 'result-0510-test';

lfFilePrefix = 'frame-480';
motionFile = 'motion.mat';

% logFile = fopen(fullfile(workDir, subDir, [logFileName '.txt']), 'a');
lf0 = rgb2gray(im2double(imread(fullfile(workDir, subDir, [lfFilePrefix sprintf('_%04d.png', 0)]))));
lf0param = LFReadMetadata(fullfile(workDir, subDir, [lfFilePrefix sprintf('_%04d.json', 0)]));
lf0 = raw2LF(lf0, lf0param.camParam.resol);
lf1 = rgb2gray(im2double(imread(fullfile(workDir, subDir, [lfFilePrefix sprintf('_%04d.png', 1)]))));
lf1 = raw2LF(lf1, lf0param.camParam.resol);

sigma = 2;

% pre-filtering with Gaussian
% TODO: combine with the deriv filters
if sigma > 0
    sigma_rad = ceil(sigma*1.5);
    lf0 = prefilter(lf0, sigma, sigma_rad);
    lf1 = prefilter(lf1, sigma, sigma_rad);
else
    sigma_rad = 0;
end

LFs = cat(5, lf0, lf1);

fprintf('Compute the derivatives...');
tic;
% deriv filters
% three-point stencil
h1 = [-1 0 1] / 2;
% d1 = fspecial('gaussian',[1 3], 0.5);
% h1 = conv(d1, h1);

% five-point stencil
h2 = [1 -8 0 8 -1]/12;
% d2 = fspecial('gaussian',[1 5], 1);
% h2 = conv(d2, h2);
[Lt, Li, Lj, Lk, Ll] = partialDeriv(LFs, h1, h2);
NPixels = numel(Lt);
Lt = reshape(Lt(:,:,:,:), [NPixels, 1]);
Li = reshape(Li(:,:,:,:), [NPixels, 1]);
Lj = reshape(Lj(:,:,:,:), [NPixels, 1]);
Lk = reshape(Lk(:,:,:,:), [NPixels, 1]);
Ll = reshape(Ll(:,:,:,:), [NPixels, 1]);
toc;

%% estimate alpha
alpha = (Li ./ Lk + Lj ./ Ll) / 2;
alpha_f = LF2Raw(reshape(alpha,size(lf0)));