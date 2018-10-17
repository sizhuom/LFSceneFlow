function showLF( LF )
%SHOWLF Display subaperture images of a light field

relrad = 0.6; % only show the center region within a relative radius
nx = 1; % show a grid of (2*nx+1)x(2*ny+1)
ny = 1;

if isa(LF, 'double') && max(LF(:)) > 1
    LF = LF / 256;
end

sy = size(LF, 1);
sx = size(LF, 2);
centery = floor((1+sy) / 2);
centerx = floor((1+sx) / 2);
rady = (centery - 1) * relrad;
radx = (centerx - 1) * relrad;
stepy = floor(rady / ny);
stepx = floor(radx / nx);
idy = centery-stepy*ny:stepy:centery+stepy*ny;
idx = centerx-stepx*nx:stepx:centerx+stepx*nx;

for i = 1:2*ny+1
    for j = 1:2*nx+1
        subplot(2*ny+1, 2*nx+1, (i-1)*(2*nx+1)+j);
        if (size(LF, 5) > 1)
            imshow(squeeze(LF(idy(i),idx(j),:,:,1:3)));
        else
            imshow(squeeze(LF(idy(i),idx(j),:,:)));
    end
end


end