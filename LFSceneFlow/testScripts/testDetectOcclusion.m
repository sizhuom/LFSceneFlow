workDir = fullfile(LFTopDir, 'Images/Sim', '0829-hand');
subDir = '.';
lfFilePrefix = 'frame-480';
motionFile = 'motion.mat';

lf0 = rgb2gray(im2double(imread(fullfile(workDir, subDir, [lfFilePrefix sprintf('_%04d.png', 0)]))));
lf0param = LFReadMetadata(fullfile(workDir, subDir, [lfFilePrefix sprintf('_%04d.json', 0)]));
lf0 = raw2LF(lf0, lf0param.camParam.resol);
if isfield(lf0param.camParam, 'H')
    H = lf0param.camParam.H;
else
    H = genIntrinsics2(lf0param.camParam.resol, lf0param.camParam.apert,...
        lf0param.camParam.fov, lf0param.camParam.fLen);
end

intr = struct('H', H, 'sz', size(lf0), 'S', 1000, 'D', 1000);

lf1 = rgb2gray(im2double(imread(fullfile(workDir, subDir, [lfFilePrefix sprintf('_%04d.png', 1)]))));
lf1 = raw2LF(lf1, lf0param.camParam.resol);

images = cat(5, lf0, lf1);

load(fullfile(workDir,subDir,sprintf('gtflow_%04d.mat',1)));
gtFlow = cat(5, gtFlowx, gtFlowy, gtFlowz) * intr.S;

tic;
occ = detect_occlusion(gtFlow, images, intr);
toc;