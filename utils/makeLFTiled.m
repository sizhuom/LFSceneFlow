function tiled = makeLFTiled( LF )
%SHOWLF Display subaperture images of a light field

method = 'bilinear';
if (size(LF, 5) == 4)
    LF = LF(:,:,:,:,1:3);
end

[sy,sx,sv,su,sc] = size(LF);
scale = 100 / su;
su2 = ceil(su * scale);
sv2 = ceil(sv * scale);
tiled = zeros(sy*sv2, sx*su2, sc);

for x = 1:sx
    for y = 1:sy
        im = imresize(squeeze(LF(y, x, :, :, :)), scale, method);
        tiled((y-1)*sv2+1:y*sv2,(x-1)*su2+1:x*su2,:) = im;
    end
end

mx = max(tiled(:));
if mx > 256
    tiled = double(tiled) / 65536;
elseif mx > 1
    tiled = double(tiled) / 256;
end

% for x = 1:sx
%     tiled(:,x*su2,:) = 1;
% end
% for y = 1:sy
%     tiled(y*sv2,:,:) = 1;
% end

end

