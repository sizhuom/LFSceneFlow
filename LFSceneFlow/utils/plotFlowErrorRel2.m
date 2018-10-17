function f = plotFlowErrorRel2( flowx, flowy, flowz, gtFlowx, gtFlowy, gtFlowz )
%PLOTFLOWERRORREL2 Plot the relative L2 error map of the scene flow
% (relative to the maximum maginitude of the flow)

% Compute error
xMap = flowx-gtFlowx;
yMap = flowy-gtFlowy;
zMap = flowz-gtFlowz;
errMap = sqrt(xMap.^2+yMap.^2+zMap.^2);

% Clamp
magMap = sqrt(gtFlowx.^2+gtFlowy.^2+gtFlowz.^2);
maxval = max(magMap(:));
errMap = errMap / maxval;
errMap(errMap==Inf | errMap==-Inf) = NaN;
errMapClean = errMap(~isnan(errMap));
meanErr = mean(errMapClean(:));
medianErr = median(errMapClean(:));
errMap(errMap > 1) = 1;

% Replace NaN with white, as in
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/140607
nanMap = isnan(errMap);
errMap(nanMap) = 1.1;
f = figure(); imagesc(errMap,[0 1.1]); colorbar;
colordata = colormap;
colordata(end,:) = [1 1 1];
colormap(colordata);
title(sprintf('mean=%g,median=%g',meanErr,medianErr));

end

