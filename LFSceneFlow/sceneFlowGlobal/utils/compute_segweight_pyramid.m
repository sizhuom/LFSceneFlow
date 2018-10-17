function [ P ] = compute_segweight_pyramid( segMap, sigma, r, nL, ratio )
%COMPUTE_SEGWEIGHT_PYRAMID Compute a pyramid of weight map based on
%segmentation

P   = cell(nL,1);
if ~isempty(segMap)
    P{1}= computeWeightMap(segMap, sigma, r);
    sz = [size(segMap,1) size(segMap,2) size(segMap,3) size(segMap,4)];
    
    for m = 2:nL
        tmp = resample_flow(segMap, floor(sz*ratio), 'nearest');
        tmpr = round(r * (ratio^(m-1)));
        tmpr(tmpr<3) = 3;
        P{m} = computeWeightMap(tmp, sigma, tmpr);
    end
else
    for m = 1:nL
        P{m} = [];
    end
end


end

function weight = computeWeightMap(segMap, sigma, r)

minMap = imerode(segMap, ones(r(1),r(1),r(2),r(2)));
maxMap = imdilate(segMap, ones(r(1),r(1),r(2),r(2)));
diffMap = double(minMap ~= maxMap);

weight = exp(-diffMap / (sigma^2));

end