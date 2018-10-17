function runCLGSim( workDir,logFileName,frameNo0,frameNo1,c,param,initFlow,gtFlow,alphaMap )
%RUNCLGSIM Run CLG on simulated data

if nargin < 9
    alphaMap = [];
end

resultDir = fullfile(workDir, logFileName);
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

diary(fullfile(resultDir,[logFileName,'_diary.txt']));

%%
logFile = fopen(fullfile(resultDir, [logFileName '.txt']), 'a');

listing = dir(fullfile(workDir, 'frame*.json'));
fprintf('Reading frame %d: %s\n', 0, listing(frameNo0).name);
[lf0, lfParam] = LFRead(fullfile(workDir, listing(frameNo0).name));
lf0 = LFRGB2Gray(im2double(lf0(:,:,:,:,:)));
H = diag([1000 1000 1 1 1]) * lfParam.camParam.H;

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

if getParamField(param, 'debug')
    debug_dir = fullfile(resultDir,'iters');
    if ~exist(debug_dir, 'dir')
        mkdir(debug_dir);
    end
    param{end+1} = 'debug_dir';
    param{end+1} = debug_dir;
end

fprintf('Reading frame %d: %s\n', 1, listing(frameNo1).name);
lf1 = LFRead(fullfile(workDir, listing(frameNo1).name));
lf1 = LFRGB2Gray(im2double(lf1(:,:,:,:,:)));
[lf1, ~] = makePatch(lf1, patchCoord);
lf1 = LFResize(lf1, ratio);

fprintf('Computing flow...\n');
method = 'clg2d';

[~,lf0name,~] = fileparts(listing(frameNo0).name);
extra = struct();
if ~isempty(initFlow)
    extra.xyflow = initFlow(:,:,1:2);
    %     extra.initFlow = initFlow;
end
if isempty(alphaMap)
    alphaMap = double(readCocolibDisparity(fullfile(workDir,lf0name,'disparity_epi_filter/disparity_estimate_filtered_all_views.h5')));
end
alphaMap = makePatch(alphaMap,patchCoord);
alphaMap = LFResize(alphaMap, ratio, [], 'nearest');
alphaMap = alphaMap * ratio(3);
extra.alphaMap = alphaMap;

timerVal = tic;
flow = estimate_flow_interface(lf0, lf1, intr, method, param, extra);
timeElapsed = toc(timerVal);
flowc = flow;
fprintf(logFile, 'frame: %d->%d\n', frameNo0, frameNo1);
fprintf(logFile, '%s-%02d: %s\n',...
    logFileName,c,cell2str(param));
fprintf(logFile, '%s\n', struct2str(intr));
fprintf(logFile, 'time: %g seconds\n', timeElapsed);

flowx = flowc(:,:,1);
flowy = flowc(:,:,2);
flowz = flowc(:,:,3) / intr.scaleZ;
mask = ~isnan(flowx);

% Plot results
dataPrefix = sprintf('%s-%04d',logFileName, c);
f2 = plotFlowHSV(flowx,flowy,flowz);
saveas(f2, fullfile(resultDir, [dataPrefix '-hsv.png']));% Compare with ground truth
close(f2);
f3 = plotFlowXYZ(flowx,flowy,flowz);
saveas(f3, fullfile(resultDir, [dataPrefix '-xyz.png']));% Compare with ground truth
close(f3);

% Save results
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
    savefig(fullfile(resultDir, [dataPrefix '-rel1.fig']));
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

function value = getParamField(param, fieldname)
value = [];
for i = 1:2:numel(param)
    if strcmp(param{i},fieldname)
        value = param{i+1};
        break
    end
end
end