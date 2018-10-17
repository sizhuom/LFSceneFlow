global LFTopDir;
workDir = fullfile(LFTopDir, 'Images/Real/0225-pillow/result-0225-CLGp1-py-sor-halfdown-a9-owSize0/');
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
    f = plotFlowXYZ(flowx,flowy,flowz,6);
    
    outName = [fname '-xyz-fixedRange.png'];
    saveas(f, fullfile(workDir, outName));
    fprintf('Saved as %s.\n', outName);
    close(f);
end

