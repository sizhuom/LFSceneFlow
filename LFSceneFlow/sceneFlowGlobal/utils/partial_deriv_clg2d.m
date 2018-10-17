function [It, IX, IY, IZ] = partial_deriv_clg2d(images, intr, flow_prev, ...
    owSize, alphaCen, interpolation_method, deriv_filter_xy, deriv_filter_uv, ...
    interp_filter_xy, interp_filter_uv, b)

%PARTIAL_DERIV   Spatio-temporal derivatives
%   P = PARTIAL_DERIV(IMAGES, INIT) computes the spatio-temporal derivatives
%
%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2007-11-30 $
%   $Revision $
%
% Copyright 2007-2010, Brown University, Providence, RI. USA
%
%                          All Rights Reserved
%
% All commercial use of this software, whether direct or indirect, is
% strictly prohibited including, without limitation, incorporation into in
% a commercial product, use in a commercial service, or production of other
% artifacts for commercial purposes.
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for research purposes is hereby granted without fee,
% provided that the above copyright notice appears in all copies and that
% both that copyright notice and this permission notice appear in
% supporting documentation, and that the name of the author and Brown
% University not be used in advertising or publicity pertaining to
% distribution of the software without specific, written prior permission.
%
% For commercial uses contact the Technology Venture Office of Brown University
%
% THE AUTHOR AND BROWN UNIVERSITY DISCLAIM ALL WARRANTIES WITH REGARD TO
% THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
% FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR OR
% BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
% DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
% PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
% ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
% THIS SOFTWARE.

if nargin < 3
    flow_prev = zeros([intr.sz 3]);
end

if nargin < 6
    interpolation_method = 'bi-linear';
end;

% h1: deriv_filter for x,y
if nargin >= 7 && ~isempty(deriv_filter_xy)
    h1 = deriv_filter_xy;
else
    h1 = [-1 0 1] / 2;
end
% h2: deriv_filter for u,v
if nargin >= 8 && ~isempty(deriv_filter_uv)
    h2 = deriv_filter_uv;
else
    h2 = [1 -8 0 8 -1]/12; % used in Wedel etal "improved TV L1"
end;
% k1: interp_filter for x,y
if nargin >= 9 && ~isempty(interp_filter_xy)
    k1 = interp_filter_xy;
else
    k1 = [1];
end;
% k2: interp_filter for u,v
if nargin >= 10 && ~isempty(interp_filter_uv)
    k2 = interp_filter_uv;
else
    k2 = [1];
end

if nargin < 11
    b = 0.6;    % recommended in Wedal etal "improved TV-L1" 2008
end;

SY   = size(images, 1);
SX   = size(images, 2);
SV   = size(images, 3);
SU   = size(images, 4);
SOW = prod(2*owSize+1);

CY = (1+SY)/2;
CX = (1+SX)/2;

[winY,winX,winV,winU] = ndgrid(CY-owSize(1):CY+owSize(1),...
    CX-owSize(2):CX+owSize(2),...
    -owSize(3):owSize(3),-owSize(4):owSize(4));
x1 = repmat(reshape(winX,1,1,[]),SV,SU,1);
y1 = repmat(reshape(winY,1,1,[]),SV,SU,1);
u1 = repmat(reshape(winU,1,1,[]),SV,SU,1);
v1 = repmat(reshape(winV,1,1,[]),SV,SU,1);

[cenV,cenU,~] = ndgrid(1:SV,1:SU,1:SOW);
u1 = u1+cenU-bsxfun(@times,alphaCen,x1-CX);
v1 = v1+cenV-bsxfun(@times,alphaCen,y1-CY);
clear winY winX winV winU cenV cenU

flow_prev_rep = repmat(reshape(flow_prev,SV,SU,1,3),1,1,SOW,1);
flow_proj = premultHM(flow_prev_rep,intr,x1,y1,u1,v1);
x2          = x1 + flow_proj(:,:,:,1);
y2          = y1 + flow_proj(:,:,:,2);
u2          = u1 + flow_proj(:,:,:,3);
v2          = v1 + flow_proj(:,:,:,4);

% Record out of boundary pixels
B = (x1>SX) | (x1<1) | (y1>SY) | (y1<1) ...
    | (u1>SU) | (u1<1) | (v1>SV) | (v1<1) ...
    | (x2>SX) | (x2<1) | (y2>SY) | (y2<1) ...
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
        im1s = interpn(images(:,:,:,:,1),y1,x1,v1,u1,method);
        It      = warpIm - im1s;
        
    else
        % Color
        error('partial_deriv_clg2d: support for color images not implemented.');
%         warpIm  = zeros(size(images(:,:,:,:,:,1)));
%         for j = 1:size(images,5)
%             warpIm(:,:,:,:,j) = interpn(images(:,:,:,:,j,2),y2,x2,v2,u2,method, NaN);
%         end;
%         It      = warpIm - images(:,:,:,:,:,1);
        
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
    
    if nargout > 1
        
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
            error('partial_deriv_clg2d: support for color images not implemented.');
%             Ix  = zeros(size(images(:,:,:,:,:,1)));
%             Iy  = zeros(size(images(:,:,:,:,:,1)));
%             Iu  = zeros(size(images(:,:,:,:,:,1)));
%             Iv  = zeros(size(images(:,:,:,:,:,1)));
%             
%             for j = 1:size(images,5)
%                 Ix(:,:,:,:,j)  = interpn(I2x(:,:,:,:,j),y2,x2,v2,u2,method);
%                 Iy(:,:,:,:,j)  = interpn(I2y(:,:,:,:,j),y2,x2,v2,u2,method);
%                 Iu(:,:,:,:,j)  = interpn(I2u(:,:,:,:,j),y2,x2,v2,u2,method);
%                 Iv(:,:,:,:,j)  = interpn(I2v(:,:,:,:,j),y2,x2,v2,u2,method);
%             end;
            
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
        
        I1x  = interpn(I1x,y1,x1,v1,u1,method);
        I1y  = interpn(I1y,y1,x1,v1,u1,method);
        I1u  = interpn(I1u,y1,x1,v1,u1,method);
        I1v  = interpn(I1v,y1,x1,v1,u1,method);
            
        Ix  = b*Ix+(1-b)*I1x;
        Iy  = b*Iy+(1-b)*I1y;
        Iu  = b*Iu+(1-b)*I1u;
        Iv  = b*Iv+(1-b)*I1v;
        
        
    end;
    
    % Disable those out-of-boundary pixels in warping
    It(B)   = 0;
    if nargout > 1
        Ix(B) = 0;
        Iy(B) = 0;
        Iu(B) = 0;
        Iv(B) = 0;
    end;
    
else
    error('partial_deriv: unknown interpolation method!');
end;

% Back-project to 3D spatial derivatives
if nargout > 1
    if size(images, 6) == 1 % grayscale
        IXYZ = postmultHM(cat(4,Ix,Iy,Iu,Iv),intr,x2,y2,u2,v2);
        IX = IXYZ(:,:,:,1);
        IY = IXYZ(:,:,:,2);
        IZ = IXYZ(:,:,:,3);
    else % color
        error('partial_deriv_clg2d: support for color images not implemented.');
        %     SC = size(images, 5);
        %     IX = zeros(SY,SX,SV,SU,SC);
        %     IY = zeros(SY,SX,SV,SU,SC);
        %     IZ = zeros(SY,SX,SV,SU,SC);
        %     for c = 1:SC
        %         IXYZ = postmultHM(cat(5,Ix(:,:,:,:,c),Iy(:,:,:,:,c),...
        %             Iu(:,:,:,:,c),Iv(:,:,:,:,c)),intr,x2,y2,u2,v2);
        %         IX(:,:,:,:,c) = IXYZ(:,:,:,:,1);
        %         IY(:,:,:,:,c) = IXYZ(:,:,:,:,2);
        %         IZ(:,:,:,:,c) = IXYZ(:,:,:,:,3);
        %     end
    end
end

