function plotFlowRGB( flowx, flowy, flowz )
%PLOTFLOWRGB Plot the flow in RGB

map = zeros([size(flowx) 3]);
map(:,:,1) = (flowx - min(flowx(:))) / (max(flowx(:)) - min(flowx(:)));
map(:,:,2) = (flowy - min(flowy(:))) / (max(flowy(:)) - min(flowy(:)));
map(:,:,3) = (flowz - min(flowz(:))) / (max(flowz(:)) - min(flowz(:)));

figure; imshow(map);

end

