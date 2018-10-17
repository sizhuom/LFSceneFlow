function [ P ] = compute_alphaCen_pyramid( alphaCen, pyramid_images )
%COMPUTE_ALPHACEN_PYRAMID2 Compute a pyramid of alpha maps

P = cell(numel(pyramid_images), 1);
for i = 1:numel(P)
    P{i} = imresize(alphaCen, [size(pyramid_images{i},3),size(pyramid_images{i},4)], 'nearest');
end

end

