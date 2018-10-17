function [ lf ] = gaussFilt4D( lf, sigma, rad )
%GAUSSFILT4D 4D Gaussian filtering on a light field

if nargin < 3
    rad = ceil(sigma * 1.5);
end

for i = 1:4
    if rad(i) == 0
        continue;
    end
    l = 2*rad(i) + 1;
    x = -rad(i):rad(i);
    k = 1/sqrt(2*pi)/sigma(i)*exp(-x.^2/2/sigma(i)^2);
    k = k/sum(k);
    s = [1 1 1 1];
    s(i) = l;
    h = reshape(k, s);
    lf = imfilter(lf, h, 'corr', 'symmetric', 'same');
end

end

