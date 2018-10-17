function runCLG( jsonFile )
%RUNCLG Run CLG method from here

% First iteration
workDir = fileparts(jsonFile);
param = LFReadMetadata(jsonFile);
if isfield(param, 'p2')
    logFileNameP1 = 'result-CLGp1';
else
    logFileNameP1 = 'result-CLG';
end
if isfield(param, 'gtFlow')
    load(fullfile(workDir,param.gtFlow), 'gtFlowx', 'gtFlowy', 'gtFlowz');
    gtFlow = cat(3, gtFlowx, gtFlowy, gtFlowz);
else
    gtFlow = [];
end
[flowx, flowy, flowz] = runCLGInner(workDir, logFileNameP1, param.frame0,...
    param.frame1, param.p1, [], gtFlow);

% Second iteration: use XY-flow from first iteration
if isfield(param, 'p2')
    logFileNameP2 = 'result-CLG';
    initFlow = cat(3, flowx, flowy, flowz);
    runCLGInner(workDir, logFileNameP2, param.frame0, param.frame1, param.p2,...
        initFlow, gtFlow);
end

end

