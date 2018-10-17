global LFTopDir;
workDir = fullfile(LFTopDir, 'Images/Real/1009-candy-00')
% matFile = 'result-0605-SG-Occ-1-7-4-occ0.1.mat';
% load(fullfile(workDir, matFile));

listing = dir(fullfile(workDir, 'result-*.mat'));
for i=1:length(listing)
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
    f = plotFlowErrorRel(flowx,flowy,flowz,gtFlowx,gtFlowy,gtFlowz);
    
    savefig(fullfile(workDir, [fname '-rel.fig']));
    saveas(f, fullfile(workDir, [fname '-rel.png']));
    fprintf('Saved as %s.\n', [fname '-rel.fig']);
    close(f);
    
%     f = plotFlowAE(flowx,flowy,flowz,gtFlowx,gtFlowy,gtFlowz,0.05);
%     
%     savefig(fullfile(workDir, [fname '-ae-clamped.fig']));
%     saveas(f, fullfile(workDir, [fname '-ae-clamped.png']));
%     fprintf('Saved as %s.\n', [fname '-ae-clamped.fig']);
%     pause(1);
%     close(f);
    clear flowx flowy flowz
end

