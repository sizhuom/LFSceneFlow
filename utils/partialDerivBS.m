function [ It, Ix, Iy, Iu, Iv ] = partialDerivBS( images, bs_images, flow_prev, ...
    interpolation_method, b )
%PARTIALDERIVBS Compute partial derivatives using B-splines
%TODO: implement interp4 for B-splines for better accuracy

% blending factor for temporal average
if nargin < 5
    b = 0.6;
end

% SY   = size(images, 1);
% SX   = size(images, 2);
% SV   = size(images, 3);
% SU   = size(images, 4);
sz = get(bs_images(1),'datasize');
SY = sz(1);
SX = sz(2);
SV = sz(3);
SU = sz(4);

[y,x,v,u]   = ndgrid(1:SY,1:SX,1:SV,1:SU);
x2          = x + flow_prev(:,:,:,:,1);
y2          = y + flow_prev(:,:,:,:,2);
u2          = u + flow_prev(:,:,:,:,3);
v2          = v + flow_prev(:,:,:,:,4);

% Record out of boundary pixels
B = (x2>SX) | (x2<1) | (y2>SY) | (y2<1) ...
    | (u2>SU) | (u2<1) | (v2>SV) | (v2<1);

if strcmp(interpolation_method, 'bi-linear') || strcmp(interpolation_method, 'cubic')
    
    if strcmp(interpolation_method, 'bi-linear')
        method = 'linear';
    else
        method = 'cubic';
    end;
    
    % Matlab built-in
    % Gray-level
%     warpIm  = interpn(images(:,:,:,:,2),y2,x2,v2,u2,method);
    warpIm = interp4(bs_images(2),x2,y2,u2,v2);
    It      = warpIm - images(:,:,:,:,1);
    
    if nargout == 5
        
        % First compute derivative then warp
        % TODO: size(images, 6) > 1?
        bsdy = partial(bs_images(2), 1);
        bsdx = partial(bs_images(2), 2);
        bsdv = partial(bs_images(2), 3);
        bsdu = partial(bs_images(2), 4);
%         I2x = indirectFilter(bsdx);
%         I2y = indirectFilter(bsdy);
%         I2u = indirectFilter(bsdu);
%         I2v = indirectFilter(bsdv);
%         
%         Ix  = interpn(I2x,y2,x2,v2,u2,method);
%         Iy  = interpn(I2y,y2,x2,v2,u2,method);
%         Iu  = interpn(I2u,y2,x2,v2,u2,method);
%         Iv  = interpn(I2v,y2,x2,v2,u2,method);
        Ix = interp4(bsdx, x2, y2, u2, v2);
        Iy = interp4(bsdy, x2, y2, u2, v2);
        Iu = interp4(bsdu, x2, y2, u2, v2);
        Iv = interp4(bsdv, x2, y2, u2, v2);
        
        % Temporal average
        bsdy = partial(bs_images(1), 1);
        bsdx = partial(bs_images(1), 2);
        bsdv = partial(bs_images(1), 3);
        bsdu = partial(bs_images(1), 4);
        I1x = indirectFilter(bsdx);
        I1y = indirectFilter(bsdy);
        I1u = indirectFilter(bsdu);
        I1v = indirectFilter(bsdv);
        
        Ix  = b*Ix+(1-b)*I1x;
        Iy  = b*Iy+(1-b)*I1y;
        Iu  = b*Iu+(1-b)*I1u;
        Iv  = b*Iv+(1-b)*I1v;
        
    end;
    
    % Disable those out-of-boundary pixels in warping
    It(B)   = 0;
    if nargout == 5
        Ix(B) = 0;
        Iy(B) = 0;
        Iu(B) = 0;
        Iv(B) = 0;
    end;
    
else
    error('partial_deriv: unknown interpolation method!');
end;
end

