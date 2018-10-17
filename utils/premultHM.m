function [ Q ] = premultHM( P, intr, X, Y, U, V )
%PREMULTHM Pre-multiply by HM

if nargin < 3
    [Y, X, V, U] = ndgrid(1:intr.sz(1), 1:intr.sz(2), 1:intr.sz(3), 1:intr.sz(4)); 
end

nPixels = numel(X);
XYUV = [ X(:)'; Y(:)'; U(:)'; V(:)'];

Hl = intr.H(1:4, 1:4) * diag([intr.S intr.S intr.D intr.D]); % linear part of H
Ho = intr.H(1:4, 5) .* [intr.S intr.S intr.D intr.D]'; % offset part of H
XYUV = bsxfun(@plus, Hl * XYUV, Ho);

% assert(size(P,5) == 3); % only translation is supported
dX = P(1:nPixels) - XYUV(3,:)/intr.D .* P(2*nPixels+1:3*nPixels) / intr.scaleZ;
dY = P(nPixels+1:2*nPixels) - XYUV(4,:)/intr.D .* P(2*nPixels+1:3*nPixels) / intr.scaleZ;
dXYUV = [dX; dY; zeros(2, nPixels)];
Q = Hl \ dXYUV;
szP = size(P);
Q = reshape(Q', [szP(1:end-1) 4]);
    
end

