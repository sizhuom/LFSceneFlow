function [flowx,flowy,flowz] = runCLGInner( workDir,logFileName,frame0,frame1,param,initFlow,gtFlow)
%RUNCLGINNER Run one iteration of CLG. This function is called twice for 
%each complete runCLG call.

%% Set up log files
if param.saveResult
    resultDir = fullfile(workDir, logFileName);
    if ~exist(resultDir, 'dir')
        mkdir(resultDir);
    end
    diary(fullfile(resultDir,[logFileName,'_diary.txt']));
    logFile = fopen(fullfile(resultDir, [logFileName '.txt']), 'a');
end

%% Compute scene flow
% Load file and pre-process
fprintf('Reading frame %d: %s\n', 0, frame0);
[lf0, H] = gcReadSubImages(fullfile(workDir, frame0));
lf0 = LFRGB2Gray(im2double(lf0(:,:,:,:,:)));

patchCoord = [
    3, 11;
    3, 11;
    1, size(lf0,3);
    1, size(lf0,4);
    ];
[lf0, H] = makePatch(lf0, patchCoord, H);
ratio = [1 1 0.5 0.5];
[lf0, H] = LFResize(lf0, ratio, H);

intr = struct('H', H, 'sz', size(lf0), 'S', 1, 'D', 1, 'scaleZ', 1);

if isfield(param, 'debug') && param.debug
    debug_dir = fullfile(resultDir,'iters');
    if ~exist(debug_dir, 'dir')
        mkdir(debug_dir);
    end
    param.debug_dir = debug_dir;
end
cellParam = makeCellParam(param);

fprintf('Reading frame %d: %s\n', 1, frame1);
lf1 = gcReadSubImages(fullfile(workDir, frame1));
lf1 = LFRGB2Gray(im2double(lf1(:,:,:,:,:)));
[lf1, ~] = makePatch(lf1, patchCoord);
lf1 = LFResize(lf1, ratio);

fprintf('Computing flow...\n');
method = 'clg2d';

extra = struct();
if ~isempty(initFlow)
    extra.xyflow = initFlow(:,:,1:2);
    %     extra.initFlow = initFlow;
end
alphaMap = double(readCocolibDisparity(fullfile(workDir,frame0,'disparity_epi_filter/disparity_estimate_filtered_all_views.h5')));
alphaMap = makePatch(alphaMap,patchCoord);
alphaMap = LFResize(alphaMap, ratio, [], 'nearest');
alphaMap = alphaMap * ratio(3);
extra.alphaMap = alphaMap;

% Compute scene flow by CLG method
timerVal = tic;
flow = estimate_flow_interface(lf0, lf1, intr, method, cellParam, extra);
% flow = rand(size(lf0,3), size(lf0,4), 3);
timeElapsed = toc(timerVal);
flowc = flow;

flowx = flowc(:,:,1);
flowy = flowc(:,:,2);
flowz = flowc(:,:,3) / intr.scaleZ;
mask = ~isnan(flowx);

% Save results
if param.saveResult
    % Plot the flow
    fprintf(logFile, 'time: %g seconds\n', timeElapsed);
    
    dataPrefix = sprintf('%s',logFileName);
    f2 = plotFlowHSV(flowx,flowy,flowz);
    saveas(f2, fullfile(resultDir, [dataPrefix '-hsv.png']));% Compare with ground truth
    close(f2);
    f3 = plotFlowXYZ(flowx,flowy,flowz);
    saveas(f3, fullfile(resultDir, [dataPrefix '-xyz.png']));% Compare with ground truth
    close(f3);
    
    % Save the flow as mat file, write logs
    if isempty(gtFlow)
        save(fullfile(resultDir, [dataPrefix '.mat']),...
            'flowx', 'flowy', 'flowz', 'flow', 'param');
    else
        if size(gtFlow,1) == 1
            gtFlowx = gtFlow(:,:,1) * ones(size(flowx));
            gtFlowy = gtFlow(:,:,2) * ones(size(flowx));
            gtFlowz = gtFlow(:,:,3) * ones(size(flowx));
        else
            gtFlowx = imresize(gtFlow(:,:,1),size(flowx),'nearest');
            gtFlowy = imresize(gtFlow(:,:,2),size(flowx),'nearest');
            gtFlowz = imresize(gtFlow(:,:,3),size(flowx),'nearest');
        end
        save(fullfile(resultDir, [dataPrefix '.mat']),...
            'gtFlowx', 'gtFlowy', 'gtFlowz', 'flowx', 'flowy', 'flowz', 'param');
        
        % Plot error
        f = plotFlowErrorRel(flowx,flowy,flowz,gtFlowx,gtFlowy,gtFlowz);
%         savefig(fullfile(resultDir, [dataPrefix '-rel1.fig']));
        saveas(f, fullfile(resultDir, [dataPrefix '-rel1.png']));
        close(f);
        
        % Compute MSE
        flowx(~mask) = 0;
        flowy(~mask) = 0;
        flowz(~mask) = 0;
        errx = flowx-gtFlowx;
        erry = flowy-gtFlowy;
        errz = flowz-gtFlowz;
        msex = sum(sum(errx.^2))/numel(flowx);
        msey = sum(sum(erry.^2))/numel(flowy);
        msez = sum(sum(errz.^2))/numel(flowz);
        fprintf(logFile, 'MSE(x): %.16g\n', msex);
        fprintf(logFile, 'MSE(y): %.16g\n', msey);
        fprintf(logFile, 'MSE(z): %.16g\n', msez);
        ee = sqrt(errx.^2+erry.^2+errz.^2);
        aee = mean(ee(:));
        fprintf(logFile, 'AEE: %.16g\n', aee);
        ae = acos((flowx.*gtFlowx+flowy.*gtFlowy+flowz.*gtFlowz+1)./...
            sqrt((flowx.^2+flowy.^2+flowz.^2+1).*(gtFlowx.^2+gtFlowy.^2+gtFlowz.^2+1)));
        aae = mean(ae(:));
        fprintf(logFile, 'AAE: %.16g\n', aae);
    end
    
    fclose(logFile);
    diary off;
end

end

function cellParam = makeCellParam(param)
fn = fieldnames(param);
cellParam = cell(numel(fn)*2, 1);
for i = 1:numel(fn)
    cellParam{2*i-1} = fn{i};
    cellParam{2*i} = param.(fn{i});
end
end