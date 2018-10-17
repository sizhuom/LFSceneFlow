function [ It, Ix, Iy, Iu, Iv ] = partialDeriv_old( images, deriv_filter_xy, deriv_filter_uv, b )
%PARTIALDERIV_OLD Partial derivative

% h1: deriv_filter for x,y
if nargin >= 2 && ~isempty(deriv_filter_xy)
    h1 = deriv_filter_xy;
else
    h1 = [-1 0 1] / 2;
end
% h2: deriv_filter for u,v
if nargin >= 3 && ~isempty(deriv_filter_uv)
    h2 = deriv_filter_uv;
else
    h2 = [1 -8 0 8 -1]/12; % used in Wedel etal "improved TV L1"
end;
% blending factor for temporal average
if nargin < 4
    b = 0.5;
end

assert(size(images, 6) == 1); % only accepts Grayscale images
img1 = images(:,:,:,:,1);

Ix = imfilter(img1, h1,  'corr', 'symmetric', 'same');
Iy = imfilter(img1, h1', 'corr', 'symmetric', 'same');
Iu = imfilter(img1, reshape(h2,[1 1 1 length(h2)]),...
    'corr', 'symmetric', 'same');
Iv = imfilter(img1, reshape(h2,[1 1 length(h2) 1]),...
    'corr', 'symmetric', 'same');
if size(images, 5) > 1
    img2 = images(:,:,:,:,2);
    It = img2 - img1;
    
    % temporal average
    Ix2 = imfilter(img2, h1,  'corr', 'symmetric', 'same');
    Iy2 = imfilter(img2, h1', 'corr', 'symmetric', 'same');
    Iu2 = imfilter(img2, reshape(h2,[1 1 1 length(h2)]),...
        'corr', 'symmetric', 'same');
    Iv2 = imfilter(img2, reshape(h2,[1 1 length(h2) 1]),...
        'corr', 'symmetric', 'same');
    
    Ix = b*Ix + (1-b)*Ix2;
    Iy = b*Iy + (1-b)*Iy2;
    Iu = b*Iu + (1-b)*Iu2;
    Iv = b*Iv + (1-b)*Iv2;
else
    It = [];
end

end

