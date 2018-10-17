function [ Q ] = postmultHM( P, intr, X, Y, U, V, full6D )
%POSTMULTHM Post-multipy by HM

if nargin < 3
    [Y, X, V, U] = ndgrid(1:intr.sz(1), 1:intr.sz(2), 1:intr.sz(3), 1:intr.sz(4)); 
end
if nargin < 7
    full6D = false;
end

XYUV = [ X(:)'; Y(:)'; U(:)'; V(:)'];

Hl = intr.H(1:4, 1:4) * diag([intr.S intr.S intr.D intr.D]); % linear part of H
Ho = intr.H(1:4, 5) .* [intr.S intr.S intr.D intr.D]'; % offset part of H
XYUV = bsxfun(@plus, Hl * XYUV, Ho);

assert(~full6D); % only translation is supported
Px = reshape(P, [], 4);
Hlinv = inv(Hl);
Px = Px * Hlinv(:,1:2);
gX = Px(:, 1);
gY = Px(:, 2);
gZ = (-gX .* XYUV(3,:)'/intr.D - gY .* XYUV(4,:)'/intr.D) / intr.scaleZ;
Q = [gX gY gZ];
szP = size(P);
Q = reshape(Q, [szP(1:end-1) 3]);

end

