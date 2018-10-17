function [ gtFlowX, gtFlowY, gtFlowZ ] = calcGTFlow( workDir, ccFilename, frameNum )
%CALCGTFLOW Calculate ground truth flow from simulated results

load(fullfile(workDir, 'motion.mat'), 'dT');
cdFile = fopen(fullfile(workDir,'colorData.txt'),'r');
colors = fscanf(cdFile,'%f,',[3,Inf]);

imcc = imread(fullfile(workDir, ccFilename));
[~,name,~] = fileparts(ccFilename);
param = LFReadMetadata(fullfile(workDir, [name '.json']));
sz = param.camParam.resol;
lfcc = raw2LF(imcc, sz);

if isa(imcc, 'uint8')
    colors = colors * 255;
elseif isa(imcc, 'uint16')
    colors = colors * 65535;
end

gtFlowX = zeros(sz);
gtFlowY = zeros(sz);
gtFlowZ = zeros(sz);
for i = 1:size(dT,1)
    mask = lfcc(:,:,:,:,1)==colors(1,i)&lfcc(:,:,:,:,2)==colors(2,i)...
        &lfcc(:,:,:,:,3)==colors(3,i);
    gtFlowX(mask) = dT(i,1,frameNum);
    gtFlowY(mask) = dT(i,2,frameNum);
    gtFlowZ(mask) = dT(i,3,frameNum);
end

end

