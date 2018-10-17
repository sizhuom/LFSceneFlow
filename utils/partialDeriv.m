function [ It, Ix, Iy, Iu, Iv ] = partialDeriv( images, flow_prev, ...
    interpolation_method, deriv_filter_xy, deriv_filter_uv, ...
    interp_filter_xy, interp_filter_uv, b, mask )
%PARTIALDERIV Partial derivative

if nargin < 2 || isempty(flow_prev)
    flow_prev = zeros([size(images,1) size(images,2) size(images,3)...
        size(images,4) 4]);
end
if nargin < 3
    interpolation_method = 'bi-linear';
end;

% h1: deriv_filter for x,y
if nargin >= 4
    h1 = deriv_filter_xy;
else
    h1 = [-1 0 1] / 2;
end
% h2: deriv_filter for u,v
if nargin >= 5
    h2 = deriv_filter_uv;
else
    h2 = [1 -8 0 8 -1]/12; % used in Wedel etal "improved TV L1"
end;
% k1: interp_filter for x,y
if nargin >= 6
    k1 = interp_filter_xy;
else
    k1 = [1];
end;
% k2: interp_filter for u,v
if nargin >= 7
    k2 = interp_filter_uv;
else
    k2 = [1];
end

if nargin < 8
    b = 0.6;    % recommended in Wedal etal "improved TV-L1" 2008
end;

if nargin < 9
    mask = [];
end

SY   = size(images, 1);
SX   = size(images, 2);
SV   = size(images, 3);
SU   = size(images, 4);

[y,x,v,u]   = ndgrid(1:SY,1:SX,1:SV,1:SU);
x2          = x + flow_prev(:,:,:,:,1);
y2          = y + flow_prev(:,:,:,:,2);
u2          = u + flow_prev(:,:,:,:,3);
v2          = v + flow_prev(:,:,:,:,4);

% Record out of boundary pixels
B = (x2>SX) | (x2<1) | (y2>SY) | (y2<1) ...
    | (u2>SU) | (u2<1) | (v2>SV) | (v2<1);

if size(images, 6) ~= 1
    B = repmat(B, [1 1 1 1 size(images, 5)]);
end;

if size(images, 5) == 1
    images = repmat(images, [1 1 1 1 2]);
end
if size(images, 6) == 1
    img1 = images(:,:,:,:,1);
    img2 = images(:,:,:,:,2);
else
    img1 = images(:,:,:,:,:,1);
    img2 = images(:,:,:,:,:,2);
end;

if strcmp(interpolation_method, 'bi-cubic')
    % TODO: adapt this?
    error('bi-cubic interpolation is not available!');
    %
    %     %h = [-0.5, 0, 0.5]; % consistent with bi-cubic interpolation
    %
    %     % Bicubic interpolation
    %     if nargout == 1
    %
    %         if size(images, 4) == 1
    %             % gray-level
    %             warpIm = interp2_bicubic(images(:,:,2),x2,y2, h);
    %         else
    %             % color
    %             warpIm  = zeros(size(images(:,:,:,1)));
    %             for j = 1:size(images,3)
    %                 warpIm(:,:,j) = interp2_bicubic(images(:,:,j, 2),x2,y2, h);
    %             end;
    %         end;
    %
    %     elseif nargout == 3
    %
    %         if size(images, 4) == 1
    %             % gray-level
    %             [warpIm Ix Iy] = interp2_bicubic(images(:,:,2),x2,y2, h);
    %         else
    %             % color
    %             warpIm  = zeros(size(images(:,:,:,1)));
    %             Ix      = warpIm;
    %             Iy      = warpIm;
    %             for j = 1:size(images,3)
    %                 [warpIm(:,:,j) Ix(:,:,j) Iy(:,:,j)] = interp2_bicubic(images(:,:,j,2),x2,y2, h);
    %             end;
    %         end;
    %
    %     else
    %         error('partial_deriv: number of output wrong!');
    %     end;
    %
    %     indx        = isnan(warpIm);
    %     if size(images, 4) == 1
    %         It          = warpIm - images(:,:,1);
    %     else
    %         It          = warpIm - images(:,:,:,1);
    %     end;
    %
    %     % Disable those out-of-boundary pixels in warping
    %     It(indx)    = 0;
    %     if nargout == 3
    %
    %         % Temporal average
    %         I1x = imfilter(img1, h,  'corr', 'symmetric', 'same');  %
    %         I1y = imfilter(img1, h', 'corr', 'symmetric', 'same');
    %
    %         % b = 0.6;    % recommended in Wedal etal "improved TV-L1" 2008
    %         % b = 0.5;
    % %         b = 1;
    %
    %
    %         Ix  = b*Ix+(1-b)*I1x;
    %         Iy  = b*Iy+(1-b)*I1y;
    %
    %         Ix(indx) = 0;
    %         Iy(indx) = 0;
    %     end;
    %
elseif strcmp(interpolation_method, 'bi-linear') || strcmp(interpolation_method, 'cubic')
    
    if strcmp(interpolation_method, 'bi-linear')
        method = 'linear';
    else
        method = 'cubic';
    end;
    
    % Matlab built-in
    if size(images, 6) == 1
        % Gray-level
        warpIm  = interpn(images(:,:,:,:,2),y2,x2,v2,u2,method);
        It      = warpIm - images(:,:,:,:,1);
        
    else
        % Color
        warpIm  = zeros(size(images(:,:,:,:,:,1)));
        for j = 1:size(images,5)
            warpIm(:,:,:,:,j) = interpn(images(:,:,:,:,j,2),y2,x2,v2,u2,method, NaN);
        end;
        It      = warpIm - images(:,:,:,:,:,1);
        
    end;
    
    % spline-based bicubic interpolation code with incorrect derivatives (below)
    %         if size(images, 4) == 1
    %             % gray-level
    %             warpIm = interp2_bicubic(images(:,:,2),x2,y2, h);
    %             It      = warpIm - images(:,:,1);
    %         else
    %             % color
    %             warpIm  = zeros(size(images(:,:,:,1)));
    %             for j = 1:size(images,3)
    %                 warpIm(:,:,j) = interp2_bicubic(images(:,:,j, 2),x2,y2, h);
    %             end;
    %             It      = warpIm - images(:,:,:,1);
    %         end;
    
    if nargout == 5
        
        % First compute derivative then warp
        I2x = imfilter(img2, h1,  'corr', 'symmetric', 'same');
        I2y = imfilter(img2, h1', 'corr', 'symmetric', 'same');
        I2u = imfilter(img2, reshape(h2,[1 1 1 length(h2)]),...
            'corr', 'symmetric', 'same');
        I2v = imfilter(img2, reshape(h2,[1 1 length(h2) 1]),...
            'corr', 'symmetric', 'same');
        
        % interp filter
        if numel(k1) > 1
            I2x = imfilter(I2x, k1', 'corr', 'symmetric', 'same');
            I2x = imfilter(I2x, reshape(k2,[1 1 length(k2) 1]),...
                'corr', 'symmetric', 'same');
            I2x = imfilter(I2x, reshape(k2,[1 1 1 length(k2)]),...
                'corr', 'symmetric', 'same');
            I2y = imfilter(I2y, k1, 'corr', 'symmetric', 'same');
            I2y = imfilter(I2y, reshape(k2,[1 1 length(k2) 1]),...
                'corr', 'symmetric', 'same');
            I2y = imfilter(I2y, reshape(k2,[1 1 1 length(k2)]),...
                'corr', 'symmetric', 'same');
            I2u = imfilter(I2u, k1, 'corr', 'symmetric', 'same');
            I2u = imfilter(I2u, k1', 'corr', 'symmetric', 'same');
            I2u = imfilter(I2u, reshape(k2,[1 1 length(k2) 1]),...
                'corr', 'symmetric', 'same');
            I2v = imfilter(I2v, k1, 'corr', 'symmetric', 'same');
            I2v = imfilter(I2v, k1', 'corr', 'symmetric', 'same');
            I2v = imfilter(I2v, reshape(k2,[1 1 1 length(k2)]),...
                'corr', 'symmetric', 'same');
        end
        
        if size(images, 6) == 1
            % Gray-level
            Ix  = interpn(I2x,y2,x2,v2,u2,method);
            Iy  = interpn(I2y,y2,x2,v2,u2,method);
            Iu  = interpn(I2u,y2,x2,v2,u2,method);
            Iv  = interpn(I2v,y2,x2,v2,u2,method);
            
            %             % spline-based bicubic interpolation code with incorrect derivatives
            %             Ix  = interp2_bicubic(I2x,x2,y2, h);
            %             Iy  = interp2_bicubic(I2y,x2,y2, h);
            
        else
            % Color
            Ix  = zeros(size(images(:,:,:,:,:,1)));
            Iy  = zeros(size(images(:,:,:,:,:,1)));
            Iu  = zeros(size(images(:,:,:,:,:,1)));
            Iv  = zeros(size(images(:,:,:,:,:,1)));
            
            for j = 1:size(images,5)
                Ix(:,:,:,:,j)  = interpn(I2x(:,:,:,:,j),y2,x2,v2,u2,method);
                Iy(:,:,:,:,j)  = interpn(I2y(:,:,:,:,j),y2,x2,v2,u2,method);
                Iu(:,:,:,:,j)  = interpn(I2u(:,:,:,:,j),y2,x2,v2,u2,method);
                Iv(:,:,:,:,j)  = interpn(I2v(:,:,:,:,j),y2,x2,v2,u2,method);
            end;
            
        end;
        
        % warp then derivative
        %         if size(images, 4) == 1
        %             tmp           = images(:,:,1);
        %         else
        %             tmp           = images(:,:,:,1);
        %         end;
        %         warpIm(B) = tmp(B);
        %         Ix        = imfilter(warpIm, h,  'corr', 'symmetric', 'same');  %
        %         Iy        = imfilter(warpIm, h', 'corr', 'symmetric', 'same');
        % end of warp then derivative
        
        
        % Temporal average
        I1x = imfilter(img1, h1,  'corr', 'symmetric', 'same');  %
        I1y = imfilter(img1, h1', 'corr', 'symmetric', 'same');
        I1u = imfilter(img1, reshape(h2,[1 1 1 length(h2)]),...
            'corr', 'symmetric', 'same');
        I1v = imfilter(img1, reshape(h2,[1 1 length(h2) 1]),...
            'corr', 'symmetric', 'same');
        
        % interp filter
        if numel(k1) > 1
            I1x = imfilter(I1x, k1', 'corr', 'symmetric', 'same');
            I1x = imfilter(I1x, reshape(k2,[1 1 length(k2) 1]),...
                'corr', 'symmetric', 'same');
            I1x = imfilter(I1x, reshape(k2,[1 1 1 length(k2)]),...
                'corr', 'symmetric', 'same');
            I1y = imfilter(I1y, k1, 'corr', 'symmetric', 'same');
            I1y = imfilter(I1y, reshape(k2,[1 1 length(k2) 1]),...
                'corr', 'symmetric', 'same');
            I1y = imfilter(I1y, reshape(k2,[1 1 1 length(k2)]),...
                'corr', 'symmetric', 'same');
            I1u = imfilter(I1u, k1, 'corr', 'symmetric', 'same');
            I1u = imfilter(I1u, k1', 'corr', 'symmetric', 'same');
            I1u = imfilter(I1u, reshape(k2,[1 1 length(k2) 1]),...
                'corr', 'symmetric', 'same');
            I1v = imfilter(I1v, k1, 'corr', 'symmetric', 'same');
            I1v = imfilter(I1v, k1', 'corr', 'symmetric', 'same');
            I1v = imfilter(I1v, reshape(k2,[1 1 1 length(k2)]),...
                'corr', 'symmetric', 'same');
        end
        
        Ix  = b*Ix+(1-b)*I1x;
        Iy  = b*Iy+(1-b)*I1y;
        Iu  = b*Iu+(1-b)*I1u;
        Iv  = b*Iv+(1-b)*I1v;
        
        
    end;
    
    % Disable those out-of-boundary pixels in warping
    It(B)   = 0;
    It(~mask) = 0;
    if nargout == 5
        Ix(B) = 0;
        Iy(B) = 0;
        Iu(B) = 0;
        Iv(B) = 0;
        Ix(~mask) = 0;
        Iy(~mask) = 0;
        Iu(~mask) = 0;
        Iv(~mask) = 0;
    end;
    
else
    error('partial_deriv: unknown interpolation method!');
end;

end

