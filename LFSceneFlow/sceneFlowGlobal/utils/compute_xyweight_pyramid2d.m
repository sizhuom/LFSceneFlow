function [ P ] = compute_xyweight_pyramid2d( xyflow, sigma, r, nL, ratio, functionType )
%COMPUTE_XYWEIGHT_PYRAMID2D Compute a pyramid of weight map based on xy flow

if nargin < 6
    functionType = 2;
end

P   = cell(nL,1);
if ~isempty(xyflow)
    if functionType == 1
        P{1} = computeWeightMap(xyflow, sigma, r);
    else
        P{1} = computeWeightMap2(xyflow, sigma, r);
    end
    sz = [size(xyflow,1) size(xyflow,2)];
    
    for m = 2:nL
        tmp = resample_flow2d(xyflow, floor(sz*ratio),'bicubic');
        tmpr = round(r * (ratio^(m-1)));
        tmpr(tmpr<1) = 1;
        if functionType == 1
            P{m} = computeWeightMap(tmp, sigma, tmpr);
        else
            P{m} = computeWeightMap2(tmp, sigma, tmpr);
        end
    end
else
    for m = 1:nL
        P{m} = [];
    end
end


end

function weight = computeWeightMap(xyflow, sigma, r)
h = [1 -8 0 8 -1]/12;
Xu = imfilter(xyflow(:,:,1), h,  'corr', 'symmetric', 'same');
Xv = imfilter(xyflow(:,:,1), h', 'corr', 'symmetric', 'same');
Yu = imfilter(xyflow(:,:,2), h,  'corr', 'symmetric', 'same');
Yv = imfilter(xyflow(:,:,2), h', 'corr', 'symmetric', 'same');

gradSq = Xu.^2+Xv.^2+Yu.^2+Yv.^2;
gradSq = gradSq / max(gradSq(:));
weight = exp(-gradSq / (sigma^2));
weight = imerode(weight, ones(r));

end

function weight = computeWeightMap2(xyflow, sigma, r)
h = [1 -8 0 8 -1]/12;
Xu = imfilter(xyflow(:,:,1), h,  'corr', 'symmetric', 'same');
Xv = imfilter(xyflow(:,:,1), h', 'corr', 'symmetric', 'same');
Yu = imfilter(xyflow(:,:,2), h,  'corr', 'symmetric', 'same');
Yv = imfilter(xyflow(:,:,2), h', 'corr', 'symmetric', 'same');

gradSq = Xu.^2+Xv.^2+Yu.^2+Yv.^2;
weight = 1./(1 + gradSq / (sigma^2));
weight = imerode(weight, ones(r));

end