global LFTopDir;
workDir = fullfile(LFTopDir, 'Images/Real/1026-lion-60/result-1103-CC-sor-py-halfdown-a9-seg-gnc/result-1103-CC-sor-py-halfdown-a9-seg-gnc-iters')
% matFile = 'result-0605-SG-Occ-1-7-4-occ0.1.mat';
% load(fullfile(workDir, matFile));

listing = dir(fullfile(workDir, 'result-*.mat'));
for i=1:length(listing)
    clear flowx flowy flowz
    
    [~,fname,~] = fileparts(listing(i).name);
    load(fullfile(workDir, [fname '.mat']));
    if ~exist('flowx', 'var')
        if exist('opq', 'var');
            flowx = centralSub(opq(:,:,:,:,1));
            flowy = centralSub(opq(:,:,:,:,2));
            flowz = centralSub(opq(:,:,:,:,3));
        else
            continue
        end
    end
    f = plotFlowHSV(flowx,flowy,flowz);
    
    saveas(f, fullfile(workDir, [fname '-hsv.png']));
    fprintf('Saved as %s.\n', [fname '-hsv.png']);
    close(f);
end

