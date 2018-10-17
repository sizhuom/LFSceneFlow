function [ S, V ] = svdTest( LF0, HM, param )
%SVDTEST SVD test on the motion matrix LHM
fprintf('Initializing...\n');
tic;

% Process the parameters
if (nargin < 4)
    param = struct();
end
if (~isfield(param, 'thresh'))
    thresh = [0 Inf];
else
    if (numel(param.thresh) == 1)
        thresh = [param.thresh Inf];
    else
        thresh = param.thresh;
    end
end
if (~isfield(param, 'extn'))
    extn = NaN(6, 1);
else
    extn = reshape(param.extn, [6 1]);
end
if (~isfield(param, 'patchu'))
    patchu = 1:size(LF0,4); 
else
    patchu = param.patchu;
end
if (~isfield(param, 'patchv'))
    patchv = 1:size(LF0,3);
else
    patchv = param.patchv;
end
if (~isfield(param, 'sigma'))
    sigma = 8;
else
    sigma = param.sigma;
end

% Convert the first light field to grayscale
LF0 = LF0(:,:,patchv,patchu,:);

toc;

% Solve for the motion for every subsequent frame
fprintf('Processing Frame ...\n');
tic;
% LF1 = LF1(:,:,patchv,patchu,:);
% pre-filtering with Gaussian
% TODO: combine with the deriv filters
if sigma > 0
    sigma_rad = ceil(sigma*1.5);
    LF0 = prefilter(LF0, sigma, sigma_rad);
%     LF1 = prefilter(LF1, sigma, sigma_rad);
else
    sigma_rad = 0;
end

% LFs = cat(5, LF0, LF1);
LFs = LF0;
toc;

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
NPixels = numel(Li);
% Lt = reshape(Lt(:,:,:,:), [NPixels, 1]);
Li = reshape(Li(:,:,:,:), [NPixels, 1]);
Lj = reshape(Lj(:,:,:,:), [NPixels, 1]);
Lk = reshape(Lk(:,:,:,:), [NPixels, 1]);
Ll = reshape(Ll(:,:,:,:), [NPixels, 1]);
% Compute the ray mask: rays that are used in the computation
rayMask = true(size(LF0));
b1 = floor(length(h1) / 2);
b2 = max(sigma_rad, floor(length(h2) / 2));
rayMask(1:b1,:,:,:) = false; rayMask(end-b1+1:end,:,:,:) = false;
rayMask(:,1:b1,:,:) = false; rayMask(:,end-b1+1:end,:,:) = false;
rayMask(:,:,1:b2,:) = false; rayMask(:,:,end-b2+1:end,:) = false;
rayMask(:,:,:,1:b2) = false; rayMask(:,:,:,end-b2+1:end) = false;
if isfield(param, 'window')
    h = param.window(1)/2;
    w = param.window(2)/2;
    uc = (1+size(LF0,4)) / 2;
    vc = (1+size(LF0,3)) / 2;
    ul = round(uc-w+0.5);
    uh = round(uc+w-0.5);
    vl = round(vc-h+0.5);
    vh = round(vc+h-0.5);
    rayMask(:,:,vl:vh,ul:uh) = false;
end
rayMask = rayMask(:);
toc;

fprintf('Compute the coefficient matrix LH^-1M...');
tic;
LM = zeros(NPixels, 6);
for j = 1 : NPixels
    LM(j, :) = [Li(j) Lj(j) Lk(j) Ll(j)] * HM(:, :, j);
end
toc;

fprintf('Perform SVD...');
tic;
A = LM(rayMask,:);
[~,S,V] = svd(A, 'econ');
toc;

end

