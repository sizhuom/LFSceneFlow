function f = plotFlowXYZ( flowx, flowy, flowz, mag )
%PLOTFLOWXYZ Plot the x,y,z components of scene flow separately

f = figure; 
if nargin < 4
    mag = max([max(abs(flowx(:))),max(abs(flowy(:))),max(abs(flowz(:)))]);
end

subplot_tight(1,3,1); imagesc(flowx,[-mag mag]); colorbar; colormap(hot(256));
pbaspect([size(flowx,2) size(flowy,1) 1]);
title('X Motion');
set(gca,'XTick',[]);
set(gca,'XTickLabel',[]);
set(gca,'YTick',[]);
set(gca,'YTickLabel',[]);

subplot_tight(1,3,2); imagesc(flowy,[-mag mag]); colorbar; colormap(hot(256));
pbaspect([size(flowx,2) size(flowy,1) 1]);
title('Y Motion');
set(gca,'XTick',[]);
set(gca,'XTickLabel',[]);
set(gca,'YTick',[]);
set(gca,'YTickLabel',[]);

subplot_tight(1,3,3); imagesc(flowz,[-mag mag]); colorbar; colormap(hot(256));
pbaspect([size(flowx,2) size(flowy,1) 1]);
title('Z Motion');
set(gca,'XTick',[]);
set(gca,'XTickLabel',[]);
set(gca,'YTick',[]);
set(gca,'YTickLabel',[]);

set(gcf, 'Position', [100, 100, 2000, 400])
end

