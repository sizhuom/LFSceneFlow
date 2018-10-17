function [ occMap ] = calcGTOcc( workDir, ccFilename0, ccFilename1, frameNum )
%CALCGTOCC Calculate ground truth occlusion map for simulated data
% occMap(i) = true if the pixel i is NOT occluded

load(fullfile(workDir, 'motion.mat'), 'dT', 'T0');
cdFile = fopen(fullfile(workDir,'colorData.txt'),'r');
colors = fscanf(cdFile,'%f,',[3,Inf]) * 255;

[ gtFlowX, gtFlowY, gtFlowZ ] = calcGTFlow( workDir, ccFilename0, frameNum );
sz = size(gtFlowX);
[Y,X,V,U] = ndgrid(1:size(gtFlowX,1),1:size(gtFlowX,2),1:size(gtFlowX,3),...
    1:size(gtFlowX,4));
dx = gtFlowX - U .* gtFlowZ;
dy = gtFlowY - V .* gtFlowZ;
dx(dx>0) = ceil(dx(dx>0)); % round towards larger absolute value
dx(dx<0) = floor(dx(dx<0));
dy(dy>0) = ceil(dy(dy>0));
dy(dy<0) = floor(dy(dy<0));
newX = X + dx;
newY = Y + dy;

im0 = imread(fullfile(workDir,ccFilename0));
cc0 = raw2LF(im0, sz);
cc0r = cc0(:,:,:,:,1);
cc0g = cc0(:,:,:,:,2);
cc0b = cc0(:,:,:,:,3);
im1 = imread(fullfile(workDir,ccFilename1));
cc1 = raw2LF(im1, sz);
cc1r = cc1(:,:,:,:,1);
cc1g = cc1(:,:,:,:,2);
cc1b = cc1(:,:,:,:,3);

newX(newX < 1) = 1;
newX(newX > sz(2)) = sz(2);
newY(newY < 1) = 1;
newY(newY > sz(1)) = sz(1);
linInd = sub2ind(sz, newY, newX, V, U);
occMap = cc0r == cc1r(linInd) & cc0g == cc1g(linInd) ...
    & cc0b == cc1b(linInd);

end

