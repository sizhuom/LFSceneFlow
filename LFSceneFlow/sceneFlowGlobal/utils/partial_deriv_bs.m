function [It, IX, IY, IZ] = partial_deriv_bs(images, bs_images, intr, flow_prev, ...
    interpolation_method, b)
%PARTIAL_DERIV_BS Summary of this function goes here
%   Detailed explanation goes here

if nargin < 6
    b = 0.6;    % recommended in Wedal etal "improved TV-L1" 2008
end;

% SY   = size(images, 1);
% SX   = size(images, 2);
% SV   = size(images, 3);
% SU   = size(images, 4);
sz = get(bs_images(1),'datasize');
SY = sz(1);
SX = sz(2);
SV = sz(3);
SU = sz(4);

flow_proj = premultHM(flow_prev, intr);

[y,x,v,u]   = ndgrid(1:SY,1:SX,1:SV,1:SU);
x2          = x + flow_proj(:,:,:,:,1);
y2          = y + flow_proj(:,:,:,:,2);
u2          = u + flow_proj(:,:,:,:,3);
v2          = v + flow_proj(:,:,:,:,4);

[It, Ix, Iy, Iu, Iv] = partialDerivBS(images, bs_images, flow_proj, interpolation_method,...
    b);

% Back-project to 3D spatial derivatives
IXYZ = postmultHM(cat(5,Ix,Iy,Iu,Iv),intr,x2,y2,u2,v2);
IX = IXYZ(:,:,:,:,1);
IY = IXYZ(:,:,:,:,2);
IZ = IXYZ(:,:,:,:,3);
        
end

