function [ P ] = compute_aweight_pyramid2d( pyramid_alphaCen, sigma, r, nL, ratio, functionType )
%COMPUTE_XYWEIGHT_PYRAMID2D Compute a pyramid of weight map based on xy flow

if nargin < 6
    functionType = 2;
end

P   = cell(nL,1);
if functionType == 1
    P{1}= computeWeightMap(1./pyramid_alphaCen{1}, sigma, r);
else
    P{1} = computeWeightMap2(1./pyramid_alphaCen{1},sigma,r);
end

for m = 2:nL
    tmpr = round(r * (ratio^(m-1)));
    tmpr(tmpr<1) = 1;
    if functionType == 1
        P{m} = computeWeightMap(1./pyramid_alphaCen{m}, sigma, tmpr);
    else
        P{m} = computeWeightMap2(1./pyramid_alphaCen{m}, sigma, tmpr);
    end
end


end

function weight = computeWeightMap(alphaMap, sigma, r)
h = [1 -8 0 8 -1]/12;
Xu = imfilter(alphaMap, h,  'corr', 'symmetric', 'same');
Xv = imfilter(alphaMap, h', 'corr', 'symmetric', 'same');

gradSq = Xu.^2+Xv.^2;
weight = exp(-gradSq / (sigma^2));
weight = imerode(weight, ones(r));

end

function weight = computeWeightMap2(alphaMap, sigma, r)
h = [1 -8 0 8 -1]/12;
Xu = imfilter(alphaMap, h,  'corr', 'symmetric', 'same');
Xv = imfilter(alphaMap, h', 'corr', 'symmetric', 'same');

gradSq = Xu.^2+Xv.^2;
weight = 1 ./ (1 + gradSq / (sigma^2));
weight = imerode(weight, ones(r));

end