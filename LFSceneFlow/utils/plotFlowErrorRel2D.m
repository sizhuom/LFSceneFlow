function f = plotFlowErrorRel2D( flowx, flowy, gtFlowx, gtFlowy )
%PLOTFLOWERRORREL2D Plot the relative L2 error map of the scene flow

% Compute error
xMap = flowx-gtFlowx;
yMap = flowy-gtFlowy;
errMap = sqrt(xMap.^2+yMap.^2);

% Clamp
magMap = sqrt(gtFlowx.^2+gtFlowy.^2);
errMap = errMap ./ magMap;
errMap(errMap==Inf | errMap==-Inf) = NaN;
errMapClean = errMap(~isnan(errMap));
meanErr = mean(errMapClean(:));
medianErr = median(errMapClean(:));
errMap(errMap > 1) = 1;

% Replace NaN with white, as in
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/140607
nanMap = isnan(errMap);
errMap(nanMap) = 1.1;
axisColor = 'black';
f = figure(); imagesc(errMap,[0 1.02]); colorbar('XColor',axisColor,'YColor',axisColor);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'TickLength',[0 0]);
colordata = colormap;
colordata(end,:) = [1 1 1];
colormap(colordata);
pbaspect([size(flowx,2) size(flowy,1) 1]);
set(gca,'XColor',axisColor);
set(gca,'YColor',axisColor);
title('Relative Error','Color',axisColor);
text(0.5,-0.1,...
    sprintf('mean=%g,median=%g',meanErr,medianErr),...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'center', 'Units', 'normalized');

end

