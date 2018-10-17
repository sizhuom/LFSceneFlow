function [ P ] = compute_xyweight_pyramid3( xyflow, sigma, r, nL, ratio )
%COMPUTE_XYWEIGHT_PYRAMID3 Compute a pyramid of weight map based on xy flow

P   = cell(nL,1);
if ~isempty(xyflow)
    P{1}= computeWeightMap(xyflow, sigma, r);
    sz = [size(xyflow,1) size(xyflow,2) size(xyflow,3) size(xyflow,4)];
    
    for m = 2:nL
        tmp = resample_flow(xyflow, floor(sz*ratio));
        tmpr = round(r * (ratio^(m-1)));
        tmpr(tmpr<1) = 1;
        P{m} = computeWeightMap(tmp, sigma, tmpr);
    end
else
    for m = 1:nL
        P{m} = [];
    end
end


end

function weight = computeWeightMap(xyflow, sigma, r)

[~,Xx,Xy,Xu,Xv] = partialDeriv(xyflow(:,:,:,:,1));
[~,Yx,Yy,Yu,Yv] = partialDeriv(xyflow(:,:,:,:,2));

gradSq = Xx.^2+Xy.^2+Xu.^2+Xv.^2+Yx.^2+Yy.^2+Yu.^2+Yv.^2;
% gradSq = gradSq / max(gradSq(:));
% weight = exp(-gradSq / (sigma^2));
weight = 1./(1 + gradSq / (sigma^2));
weight = imerode(weight, ones(r(1),r(1),r(2),r(2)));

end