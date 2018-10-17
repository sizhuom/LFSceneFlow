function [ It, IX, IY, IZ ] = dbgPartialDeriv( images, intr, flow_prev )
%DBGPARTIALDERIV Just for debugging, 0623-sceneflow

global idealParams LFTopDir
% 
%% Compute the ideal derivatives
% SY   = intr.sz(1);
% SX   = intr.sz(2);
% SV   = intr.sz(3);
% SU   = intr.sz(4);
% 
% flow_proj = premultHM(flow_prev, intr);
% 
% [y,x,v,u]   = ndgrid(1:SY,1:SX,1:SV,1:SU);
% x2          = x + flow_proj(:,:,:,:,1);
% y2          = y + flow_proj(:,:,:,:,2);
% u2          = u + flow_proj(:,:,:,:,3);
% v2          = v + flow_proj(:,:,:,:,4);
% 
% % Record out of boundary pixels
% B = (x2>SX) | (x2<1) | (y2>SY) | (y2<1) ...
%     | (u2>SU) | (u2<1) | (v2>SV) | (v2<1);
% 
% if size(images, 6) ~= 1
%     B = repmat(B, [1 1 1 1 size(images, 5)]);
% end;
% 
% if size(images, 6) == 1
%     img1 = images(:,:,:,:,1);
%     img2 = images(:,:,:,:,2);
% else
%     img1 = images(:,:,:,:,:,1);
%     img2 = images(:,:,:,:,:,2);
% end;
% 
% % Matlab built-in
% method = 'linear';
% if size(images, 6) == 1
%     % Gray-level
%     warpIm  = interpn(images(:,:,:,:,2),y2,x2,v2,u2,method);
%     It      = warpIm - images(:,:,:,:,1);
%     
% else
%     % Color
%     warpIm  = zeros(size(images(:,:,:,:,:,1)));
%     for j = 1:size(images,5)
%         warpIm(:,:,:,:,j) = interpn(images(:,:,:,:,j,2),y2,x2,v2,u2,method, NaN);
%     end;
%     It      = warpIm - images(:,:,:,:,:,1);
%     
% end;
% 
% % spline-based bicubic interpolation code with incorrect derivatives (below)
% %         if size(images, 4) == 1
% %             % gray-level
% %             warpIm = interp2_bicubic(images(:,:,2),x2,y2, h);
% %             It      = warpIm - images(:,:,1);
% %         else
% %             % color
% %             warpIm  = zeros(size(images(:,:,:,1)));
% %             for j = 1:size(images,3)
% %                 warpIm(:,:,j) = interp2_bicubic(images(:,:,j, 2),x2,y2, h);
% %             end;
% %             It      = warpIm - images(:,:,:,1);
% %         end;
% 
% if nargout == 4
%     
%     % First compute derivative then warp
%     NPixels = numel(img1);
%     XYUV = [x2(:)';y2(:)';u2(:)';v2(:)';];
%     Hl = intr.H(1:4, 1:4) * diag([intr.S intr.S intr.D intr.D]); % linear part of H
%     Ho = intr.H(1:4, 5) .* [intr.S intr.S intr.D intr.D]'; % offset part of H
%     XYUV = bsxfun(@plus, Hl * XYUV, Ho);
%     
%     Ix = zeros(NPixels,1);
%     Iy = zeros(NPixels,1);
%     Iu = zeros(NPixels,1);
%     Iv = zeros(NPixels,1);
%     for i = 1:NPixels
%         grad = LFIdeal(XYUV(:,i),1,idealParams);
%         Ix(i) = grad(1);
%         Iy(i) = grad(2);
%         Iu(i) = grad(3);
%         Iv(i) = grad(4);
%     end    
%     
%     % Back-project to 3D spatial derivatives
%     P = cat(5,Ix,Iy,Iu,Iv);
%     Px = reshape(P, [size(P,1)*size(P,2)*size(P,3)*size(P,4) 4]);
%     gX = Px(:, 1);
%     gY = Px(:, 2);
%     gZ = -gX .* XYUV(3,:)'/intr.D - gY .* XYUV(4,:)'/intr.D;
%     Q = [gX gY gZ];
%     Q = reshape(Q, [size(P,1) size(P,2) size(P,3) size(P,4) 3]);
%     IX = Q(:,:,:,:,1);
%     IY = Q(:,:,:,:,2);
%     IZ = Q(:,:,:,:,3);
%     
% end;
% 
% % Disable those out-of-boundary pixels in warping
% It(B)   = 0;
% if nargout == 4
%     IX(B) = 0;
%     IY(B) = 0;
%     IZ(B) = 0;
% end;

%% Load the ideal derivatives
load(fullfile(LFTopDir, 'Images/Sim/0623-sceneflow/frame-480_0000-id.mat'));
It2 = It * 255;
IX = IX * 255;
IY = IY * 255;
IZ = IZ * 255;

%% Use numerical It
It = images(:,:,:,:,2) - images(:,:,:,:,1);

end

