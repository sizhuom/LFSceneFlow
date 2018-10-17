function showDerivative( d, sz )
%SHOWDERIVATIVE helper function which displays the derivative image
im = reshape(d, sz(1:4));
if max(im(:)) > 1
    im = (im + 256) / 512;
else
    im = (im + 1) / 2;
end
imshow(squeeze(im(ceil(sz(1) / 2), ceil(sz(2) / 2), :, :)));


end

