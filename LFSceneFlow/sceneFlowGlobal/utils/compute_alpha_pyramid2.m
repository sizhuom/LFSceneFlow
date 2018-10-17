function [ P ] = compute_alpha_pyramid2( alpha, pyramid_images )
%COMPUTE_ALPHA_PYRAMID2 Compute a pyramid of alpha maps

P = cell(numel(pyramid_images), 1);
alphaCen = centralSub(alpha);
for i = 1:numel(P)
    P{i} = imresize(alphaCen, [size(pyramid_images{i},3),size(pyramid_images{i},4)]);
end

end

