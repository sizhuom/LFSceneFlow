global LFTopDir
workDir = fullfile(LFTopDir, 'Images/Sim/0706-sceneflow');
lfFilePrefix = 'frame-480';
motionFile = 'motion.mat';

logFile = fopen(fullfile(workDir, [logFileName '.txt']), 'a');
lf0 = rgb2gray(im2double(imread(fullfile(workDir, [lfFilePrefix sprintf('_%04d.png', 0)]))));
lf0param = LFReadMetadata(fullfile(workDir, [lfFilePrefix sprintf('_%04d.json', 0)]));
lf0 = raw2LF(lf0, lf0param.camParam.resol);
lf1 = rgb2gray(im2double(imread(fullfile(workDir, [lfFilePrefix sprintf('_%04d.png', 1)]))));
lf1 = raw2LF(lf1, lf0param.camParam.resol);

if isfield(lf0param.camParam, 'H')
    H = lf0param.camParam.H;
else
    H = genIntrinsics2(lf0param.camParam.resol, lf0param.camParam.apert,...
        lf0param.camParam.fov, lf0param.camParam.fLen);
end

intr = struct('H', H, 'sz', size(lf0), 'S', 1000, 'D', 1000);

load(fullfile(workDir,motionFile), 'dT');

[gtFlowx, gtFlowy, gtFlowz] = calcGTFlow(workDir,...
    sprintf('cc-%s_%04d.png',lfFilePrefix,0),1);
gtFlowx = gtFlowx * intr.S;
gtFlowy = gtFlowy * intr.S;
gtFlowz = gtFlowz * intr.S;
P =cat(5, gtFlowx, gtFlowy, gtFlowz);

[Y, X, V, U] = ndgrid(1:intr.sz(1), 1:intr.sz(2), 1:intr.sz(3), 1:intr.sz(4));
Q = premultHM(P, intr, X, Y, U, V);

X2 = X + Q(:,:,:,:,1);
Y2 = Y + Q(:,:,:,:,2);
U2 = U + Q(:,:,:,:,3);
V2 = V + Q(:,:,:,:,4);

tic;
lf1interpn = interpn(lf1, Y2, X2, V2, U2);
toc;

tic;
bs = bsarray(lf1, 'degree', 2);
toc;
tic;
lf1interp4 = interp4(bs, X2, Y2, U2, V2);
toc;

errn = abs(lf1interpn-lf0);
errn = errn(~isnan(errn));
errn = mean(errn(:))
err4 = abs(lf1interp4-lf0);
err4 = err4(~isnan(err4));
err4 = mean(err4(:))