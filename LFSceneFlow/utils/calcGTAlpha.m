function [ gtAlpha ] = calcGTAlpha( workDir, ccFilename, frameNum )
%CALCGTALPHA Calculate the ground trugh alpha map
% alpha = Lx/Lu, -1/alpha is the slope in the EPI (x-u plot)

load(fullfile(workDir, 'motion.mat'), 'dT', 'T0');
cdFile = fopen(fullfile(workDir,'colorData.txt'),'r');
colors = fscanf(cdFile,'%f,',[3,Inf]);

imcc = imread(fullfile(workDir, ccFilename));
[~,name,~] = fileparts(ccFilename);
param = LFReadMetadata(fullfile(workDir, [name '.json']));
sz = param.camParam.resol;
fLen = param.camParam.fLen;
focus = param.camParam.focus;
F = 1/(1/fLen-1/focus);
[H, ~] = genIntrinsics2(param.camParam.resol,param.camParam.apert,...
    param.camParam.fov, param.camParam.fLen);

lfcc = raw2LF(imcc, sz);

if isa(imcc, 'uint8')
    colors = colors * 255;
elseif isa(imcc, 'uint16')
    colors = colors * 65535;
end

gtAlpha = ones(sz) * (H(3,1)/H(1,1));
for i = 1:size(dT,1)
    mask = lfcc(:,:,:,:,1)==colors(1,i)&lfcc(:,:,:,:,2)==colors(2,i)...
        &lfcc(:,:,:,:,3)==colors(3,i);
    if frameNum == 0
        Z = T0(i,3);
    else
        Z = T0(i,3)+dT(i,3,frameNum);
    end
    gtAlpha(mask) = (H(1,1)+Z*H(3,1))/(H(1,3)+Z*H(3,3));
end

end

