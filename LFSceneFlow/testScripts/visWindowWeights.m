% visualize window weights used in sceneFlowLK, for debugging
%% load file
global LFTopDir;
workDir = fullfile(LFTopDir, 'Images/Sim', '0601-sceneflow');
lfFilePrefix = 'frame-480';

lf0 = rgb2gray(im2double(imread(fullfile(workDir, [lfFilePrefix sprintf('_%04d.png', 0)]))));
lf0param = LFReadMetadata(fullfile(workDir, [lfFilePrefix sprintf('_%04d.json', 0)]));
lf0 = raw2LF(lf0, lf0param.camParam.resol);
lf1 = rgb2gray(im2double(imread(fullfile(workDir, [lfFilePrefix sprintf('_%04d.png', 1)]))));
lf1 = raw2LF(lf1, lf0param.camParam.resol);

if isfield(lf0param.camParam, 'H')
    Hnew = lf0param.camParam.H;
else
    Hnew = genIntrinsics2(lf0param.camParam.resol, lf0param.camParam.apert,...
        lf0param.camParam.fov, lf0param.camParam.fLen);
end
if exist('H', 'var') && ~isequal(H, Hnew)
    clear HM
    tic;
    HM = LFMotionMatrix(Hnew, size(lf0));
    toc;
    fprintf('Finished precomputing HM.\n');
elseif ~exist('HM', 'var')
    tic;
    HM = LFMotionMatrix(Hnew, size(lf0));
    toc;
    fprintf('Finished precomputing HM.\n');
end
H = Hnew;

alphaMap = calcGTAlpha(workDir, 'cc-frame-480_0000.png', 0);

[gtFlowx, gtFlowy, gtFlowz] = calcGTFlow(workDir, 'cc-frame-480_0000.png', 1);

%% init
% parameters
sigma = 2;

fprintf('Initializing...\n');

% Convert the first light field to grayscale
lf0 = lf0(:,:,:,:,:);

% Solve for the motion for every subsequent frame
fprintf('Processing Frame ...\n');
tic;
lf1 = lf1(:,:,:,:,:);
% pre-filtering with Gaussian
% NOTE: also filter the alphaMap?
if sigma > 0
    sigma_rad = ceil(sigma*1.5);
    lf0 = prefilter(lf0, sigma, sigma_rad);
    lf1 = prefilter(lf1, sigma, sigma_rad);
%     alphaMap = prefilter(alphaMap, sigma, sigma_rad);
else
    sigma_rad = 0;
end

LFs = cat(5, lf0, lf1);
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
NPixels = numel(Lt);
Lt = reshape(Lt(:,:,:,:), [NPixels, 1]);
Li = reshape(Li(:,:,:,:), [NPixels, 1]);
Lj = reshape(Lj(:,:,:,:), [NPixels, 1]);
Lk = reshape(Lk(:,:,:,:), [NPixels, 1]);
Ll = reshape(Ll(:,:,:,:), [NPixels, 1]);
toc;

fprintf('Compute the coefficient matrix LH^-1M...');
tic;
LM = zeros(NPixels, 3);
for j = 1 : NPixels
    LM(j, :) = [Li(j) Lj(j) Lk(j) Ll(j)] * HM(:, 1:3, j);
end
toc;

%% main loop
% parameters
halfWindow = [7 7 4 4];
useGaussWeight = true;
useShearedGauss = true;
occSigma = 0.1;
bilSigma = 0;

i = round((1+size(lf0,1))/2);
j = round((1+size(lf0,2))/2);
k = round((1+size(lf0,3))/2);
l = round((1+size(lf0,4))/2);

f1 = figure('Name','Central Subaperture');
sai = squeeze(lf0(i,j,:,:));
imshow(sai);
f2 = figure('Name','Horizontal EPI');
f3 = figure('Name','Weight Map');
f4 = figure('Name','Residual Map');
f5 = figure('Name','Alpha Map at t0');

while true
    % Compute the window weights
    if useShearedGauss
        windowSigma = (halfWindow+0.5)./3;
        % sheared window
        [I,J,K,L] = ndgrid(-halfWindow(1):+halfWindow(1),...
            -halfWindow(2):+halfWindow(2),...
            -halfWindow(3):+halfWindow(3),...
            -halfWindow(4):+halfWindow(4));
        W = exp(-(I/windowSigma(1)).^2/2-(J/windowSigma(2)).^2/2-(K/windowSigma(3)).^2/2-(L/windowSigma(4)).^2/2);
        K = K - alphaMap(i,j,k,l) * I;
        L = L - alphaMap(i,j,k,l) * J;
        I = I+i;
        J = J+j;
        K = min(max(round(K)+k,1),size(lf0,3));
        L = min(max(round(L)+l,1),size(lf0,4));
        ind = sub2ind(size(lf0), I,J,K,L);
        ind = ind(:);
        % occlusion term
        W = W(:);
        if occSigma > 0
            W = W .* exp(-((alphaMap(ind)-alphaMap(i,j,k,l))/occSigma).^2/2);
        end
        % bilateral term
        if bilSigma > 0
            W = W .* exp(-((lf0(ind)-lf0(i,j,k,l))/bilSigma).^2/2);
        end
        % solve the linear system
        A = LM(ind, :);
        b = -Lt(ind);
        ATW = bsxfun(@times,A',W');
        ATA = ATW*A;
        ATb = ATW*b;
        if (rank(ATA) == 3)
            flow = inv(ATA)*ATb;
        else
            flow = NaN(3,1);
        end
    else
        if useGaussWeight
            [I,J,K,L] = ndgrid(-halfWindow(1):+halfWindow(1),...
                -halfWindow(2):+halfWindow(2),...
                -halfWindow(3):+halfWindow(3),...
                -halfWindow(4):+halfWindow(4));
            windowSigma = (halfWindow+0.5)./3;
            W0 = exp(-(I/windowSigma(1)).^2/2-(J/windowSigma(2)).^2/2-(K/windowSigma(3)).^2/2-(L/windowSigma(4)).^2/2);
        else
            W0 = ones(size(2*halfWindow+1));
        end
        % Gaussian window
        [I,J,K,L] = ndgrid(i-halfWindow(1):i+halfWindow(1),...
            j-halfWindow(2):j+halfWindow(2),...
            k-halfWindow(3):k+halfWindow(3),...
            l-halfWindow(4):l+halfWindow(4));
        K = min(max(K,1),size(lf0,3));
        L = min(max(L,1),size(lf0,4));
        ind = sub2ind(size(lf0), I, J, K, L);
        ind = ind(:);
        W = W0(:);
        % occlusion term
        if occSigma > 0
            W = W .* exp(-((alphaMap(ind)-alphaMap(i,j,k,l))/occSigma).^2/2);
        end
        % bilateral term
        if bilSigma > 0
            W = W .* exp(-((lf0(ind)-lf0(i,j,k,l))/bilSigma).^2/2);
        end
        % solve the linear system
        A = LM(ind, :);
        b = -Lt(ind);
        ATW = bsxfun(@times,A',W');
        ATA = ATW*A;
        ATb = ATW*b;
        if (rank(ATA) == 3)
            flow = inv(ATA)*ATb;
        else
            flow = NaN(3,1);
        end
    end
    
    % Draw the window
    ul = [L(1+halfWindow(1),1,1+halfWindow(3),1),J(1+halfWindow(1),1,1+halfWindow(3),1)];
    ur = [L(1+halfWindow(1),1,1+halfWindow(3),1+2*halfWindow(4)),J(1+halfWindow(1),1,1+halfWindow(3),1+2*halfWindow(4))];
    dl = [L(1+halfWindow(1),1+2*halfWindow(2),1+halfWindow(3),1),J(1+halfWindow(1),1+2*halfWindow(2),1+halfWindow(3),1)];
    dr = [L(1+halfWindow(1),1+2*halfWindow(2),1+halfWindow(3),1+2*halfWindow(4)),J(1+halfWindow(1),1+2*halfWindow(2),1+halfWindow(3),1+2*halfWindow(4))];
    
    figure(f2); 
    epi = squeeze(lf0(i,:,k,:));
    imshow(epi);
    figure(f2); % need this line, or it will draw on f3 (why?)
    hold on;
    plot([ul(1) ur(1) dr(1) dl(1) ul(1)],...
        [ul(2) ur(2) dr(2) dl(2) ul(2)], '-r');
    hold off;
    
    % Draw the weights
    weightMap = zeros(size(lf0));
    weightMap(ind) = W;
%     weightEpi = squeeze(weightMap(i,:,k,:));
    weightEpi = squeeze(sum(sum(weightMap,1),3));
    figure(f3);
    imagesc(weightEpi);
    pbaspect([size(lf0,4) size(lf0,2) 1]);
    colorbar;
    hold on;
    plot([ul(1) ur(1) dr(1) dl(1) ul(1)],...
        [ul(2) ur(2) dr(2) dl(2) ul(2)], '-r');
    hold off;
    
    % Draw the residuals
    residualMap = zeros(size(lf0));
    gtFlow = [gtFlowx(i,j,k,l); gtFlowy(i,j,k,l); gtFlowz(i,j,k,l)];
    residualMap(ind) = W.*(A*gtFlow-b).^2;
    residualEpi = squeeze(residualMap(i,:,k,:));
    figure(f4);
    imagesc(residualEpi);
    pbaspect([size(lf0,4) size(lf0,2) 1]);
    colorbar;
    hold on;
    plot([ul(1) ur(1) dr(1) dl(1) ul(1)],...
        [ul(2) ur(2) dr(2) dl(2) ul(2)], '-r');
    hold off;
    
    % Draw the alpha map
    alphaEpi = squeeze(alphaMap(i,:,k,:));
    figure(f5);
    imagesc(alphaEpi);
    pbaspect([size(lf0,4) size(lf0,2) 1]);
    colorbar;
    hold on;
    plot([ul(1) ur(1) dr(1) dl(1) ul(1)],...
        [ul(2) ur(2) dr(2) dl(2) ul(2)], '-r');
    hold off;
    
    % Display info
    fprintf('Coordinates: (%d, %d, %d, %d)\n', j, i, l, k);
    fprintf('Alpha: %g\n', alphaMap(i,j,k,l));
    fprintf('Weigts: %g (min), %g (max)\n', min(weightEpi(:)), max(weightEpi(:)));
    fprintf('Flow: (%g, %g, %g)\n', flow(1), flow(2), flow(3));
    
    % Pick next point
    figure(f1);
    [l, k] = ginput(1);
    l = round(l);
    if l < 1
        l = 1;
    elseif l > 320
        l = 320;
    end
    k = round(k);
    if k < 1
        k = 1;
    elseif k > 240
        k = 240;
    end
end