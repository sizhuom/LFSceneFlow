function [ P ] = compute_mask_pyramid( mask, nL, ratio )
%COMPUTE_MASK_PYRAMID Compute a pyramid of logical masks

P   = cell(nL,1);
P{1}= mask;

for m = 2:nL
    P{m} = resizeMask(P{m-1}, ratio);
end

