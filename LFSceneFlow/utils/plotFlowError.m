function f = plotFlowError( flowx, flowy, flowz, gtFlowx, gtFlowy, gtFlowz, maxval )
%PLOTFLOWERROR Plot the L2 error map of the scene flow, NaN as white,
%clamped at maxval if given

% Compute error
xMap = flowx-gtFlowx;
yMap = flowy-gtFlowy;
zMap = flowz-gtFlowz;
errMap = sqrt(xMap.^2+yMap.^2+zMap.^2);

% Clamp
if nargin == 7
    errMap(errMap > maxval) = maxval;
else
    maxval = max(errMap(:));
    if maxval == 0
        maxval = 0.1;
    end
end

% Replace NaN with white, as in
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/140607
nanMap = isnan(errMap);
if (sum(nanMap(:)) > 0)
    errMap(nanMap) = maxval + maxval/10;
    f = figure(); imagesc(errMap,[0 maxval*1.1]); colorbar;
    colordata = colormap;
    colordata(end,:) = [1 1 1];
    colormap(colordata);
else
    f = figure(); imagesc(errMap,[0 maxval*1.1]); colorbar;
end

end

