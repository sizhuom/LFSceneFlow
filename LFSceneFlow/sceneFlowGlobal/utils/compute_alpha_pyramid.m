function [ P ] = compute_alpha_pyramid( alpha, L, ratio )
%COMPUTE_ALPHA_PYRAMID Compute a pyramid of alpha maps

P = cell(L, 1);
if isempty(alpha)
    for i = 1:numel(P)
        P{i} = [];
    end
else
    P{1} = alpha;
    for i = 2:numel(P)
        P{i} = LFResize(P{i-1}, ratio, [], 'nearest');
    end
end

end

