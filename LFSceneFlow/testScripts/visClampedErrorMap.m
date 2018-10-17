global LFTopDir;
workDir = fullfile(LFTopDir, 'Images/Sim', '0611-sceneflow');
maxval = 10;
% matFile = 'result-0605-SG-Occ-1-7-4-occ0.1.mat';
% load(fullfile(workDir, matFile));

listing = dir(fullfile(workDir, 'result-*.mat'));
for i=1:length(listing)
    [~,fname,~] = fileparts(listing(i).name);
    load(fullfile(workDir, [fname '.mat']));
    if ~exist('flowx', 'var')
        continue
    end
    f = plotFlowError(flowx,flowy,flowz,gtFlowx,gtFlowy,gtFlowz,maxval);
    
    savefig(fullfile(workDir, [fname '-clamped.fig']));
    saveas(f, fullfile(workDir, [fname '-clamped.png']));
    fprintf('Saved as %s.\n', [fname '-clamped.fig']);
    pause(1);
    close(f);
    
%     f = plotFlowAE(flowx,flowy,flowz,gtFlowx,gtFlowy,gtFlowz,0.05);
%     
%     savefig(fullfile(workDir, [fname '-ae-clamped.fig']));
%     saveas(f, fullfile(workDir, [fname '-ae-clamped.png']));
%     fprintf('Saved as %s.\n', [fname '-ae-clamped.fig']);
%     pause(1);
%     close(f);
end

