function [ lf1, H1 ] = LFResize( lf0, ratio, H0, method )
%LFRESIZE Naive resize function for light fields

if nargin < 4
    method = 'linear';
end

if (numel(ratio) == 1)
    ratio = ratio * ones(1,4);
end

% Prefiltering if downsampling
sigma = sqrt(1./(2*ratio));
sigma(ratio>=1) = 0;
lf0 = gaussFilt4D(lf0, sigma);

% Compute the intrinsic matrix
sz0 = [size(lf0,1) size(lf0,2) size(lf0,3) size(lf0,4)];
sz1 = floor(sz0 .* ratio);
a = 1 ./ ratio;
b = (sz0+1)/2 - (sz1+1)/2.*a;
a = [a(2) a(1) a(4) a(3)];
b = [b(2) b(1) b(4) b(3)];
Hinc = [
    diag(a) b(:);
    zeros(1,4) 1;
    ];
if nargout > 1
    H1 = H0 * Hinc;
end

% Resample the LF
[y,x,v,u] = ndgrid(1:sz1(1),1:sz1(2),1:sz1(3),1:sz1(4));
xyuv = [x(:)'; y(:)'; u(:)'; v(:)'; ones(1,numel(x))];
xyuv = Hinc * xyuv;
x = reshape(xyuv(1,:), sz1);
y = reshape(xyuv(2,:), sz1);
u = reshape(xyuv(3,:), sz1);
v = reshape(xyuv(4,:), sz1);
x(x<1) = 1; x(x>sz0(2)) = sz0(2);
y(y<1) = 1; y(y>sz0(1)) = sz0(1);
u(u<1) = 1; u(u>sz0(4)) = sz0(4);
v(v<1) = 1; v(v>sz0(3)) = sz0(3);
lf1 = zeros([sz1 size(lf0,5)]);
for c = 1:size(lf0,5)
    lf1(:,:,:,:,c) = interpn(lf0(:,:,:,:,c),y,x,v,u,method);
end

end

