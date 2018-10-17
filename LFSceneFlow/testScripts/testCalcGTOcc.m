%% compute occ map
global LFTopDir
workDir = fullfile(LFTopDir, 'Images/Sim', '0601-sceneflow');
ccFilename0 = 'cc-frame-480_0000.png';
ccFilename1 = 'cc-frame-480_0001.png';
frameNum = 1;
gtOcc = calcGTOcc(workDir, ccFilename0, ccFilename1, frameNum);

%% visualize
gtOccCentral = centralSub(gtOcc);
imagesc(gtOccCentral); colorbar;