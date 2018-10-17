function show4D( q, sz )
%SHOW4D Visualize a 4D quantity by rescale the values in [0,1]

if nargin > 1
    q = reshape(q, sz);
else
    sz = size(q);
end

% scaling
minQ = min(q(:));
rangeQ = max(q(:)) - minQ;
q = (q - minQ) / rangeQ;

% central view
figure;
imshow(squeeze(q(ceil(sz(1) / 2), ceil(sz(2) / 2), :, :)));

% circular vis
LFDispVidCirc(q);

end

