function [ lf ] = prefilter( lf, sigma, rad )
%PREFILTER Gaussian prefilter the LF 

if nargin < 3
    rad = ceil(sigma * 1.5);
end
if sigma == 0
    return;
end
l = 2*rad + 1;
x = -rad:rad;
k = 1/sqrt(2*pi)/sigma*exp(-x.^2/2/sigma^2);
k = k/sum(k);
h = reshape(k, [1 1 l 1]);
lf = imfilter(lf, h, 'corr', 'symmetric', 'same');
h = reshape(k, [1 1 1 l]);
lf = imfilter(lf, h, 'corr', 'symmetric', 'same');

% l = 2*rad + 1;
% h = fspecial('gaussian', l, sigma);
% h = reshape(h, [1 1 l l]);
% lf = imfilter(lf, h, 'corr', 'symmetric', 'same');
    
end

