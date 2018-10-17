function [ Q ] = projFlowXY( P, intr )
%PROJFLOWXY Project 3D scene flow onto 2D xy space

[Y, X, V, U] = ndgrid(round((1+intr.sz(1))/2),round((1+intr.sz(2))/2),...
    1:intr.sz(3), 1:intr.sz(4));
nPixels = numel(X);
XYUV = [ X(:)'; Y(:)'; U(:)'; V(:)'];

Hl = intr.H(1:4, 1:4) * diag([intr.S intr.S intr.D intr.D]); % linear part of H
Ho = intr.H(1:4, 5) .* [intr.S intr.S intr.D intr.D]'; % offset part of H
XYUV = bsxfun(@plus, Hl * XYUV, Ho);

dX = P(1:nPixels) - XYUV(3,:)/intr.D .* P(2*nPixels+1:3*nPixels) / intr.scaleZ;
dY = P(nPixels+1:2*nPixels) - XYUV(4,:)/intr.D .* P(2*nPixels+1:3*nPixels) / intr.scaleZ;
dXY = [dX; dY];
Q = reshape(dXY', [size(P,1) size(P,2) 2]);

end

