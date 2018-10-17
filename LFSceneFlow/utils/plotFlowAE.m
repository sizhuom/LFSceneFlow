function [ f, ae ] = plotFlowAE( flowx, flowy, flowz, gtFlowx, gtFlowy, gtFlowz, maxval )
%PLOTFLOWAAE Plot Average Angular Error (AAE) for the scene flow

ae = acos((flowx.*gtFlowx+flowy.*gtFlowy+flowz.*gtFlowz+1)./...
    sqrt((flowx.^2+flowy.^2+flowz.^2+1).*(gtFlowx.^2+gtFlowy.^2+gtFlowz.^2+1)));

% Clamp
if nargin == 7
    ae(ae > maxval) = maxval;
else
    maxval = max(ae(:));
    if maxval == 0
        maxval = 0.1;
    end
end

% Replace NaN with white, as in
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/140607
nanMap = isnan(ae);
if (sum(nanMap(:)) > 0)
    ae(nanMap) = maxval + maxval/10;
    f = figure(); imagesc(ae,[0 maxval*1.1]); colorbar;
    colordata = colormap;
    colordata(end,:) = [1 1 1];
    colormap(colordata);
else
    f = figure(); imagesc(ae,[0 maxval*1.1]); colorbar;
end

end

