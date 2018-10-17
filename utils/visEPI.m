function visEPI( lf, mode )
%VISEPI Interactive visualization of EPI

if nargin < 2
    mode = 0;
end

f1 = figure;
imshow(centralSub(lf));

u = round((size(lf,4)+1)/2);
v = round((size(lf,3)+1)/2);
f2 = figure;
if mode == 0
    epi = squeeze(lf(round((size(lf,1)+1)/2),:,v,:,:));
    imshow(epi);
else
    epi = squeeze(lf(:,round((size(lf,2)+1)/2),:,u,:));
    imshow(epi);
end
    
while true
    figure(f1);
    [u,v] = ginput(1);
    u = roundClamp(u, size(lf, 4));
    v = roundClamp(v, size(lf, 3));
    
    figure(f2);
    if mode == 0
        epi = squeeze(lf(round((size(lf,1)+1)/2),:,v,:,:));
        imshow(epi);
    else
        epi = squeeze(lf(:,round((size(lf,2)+1)/2),:,u,:));
        imshow(epi);
    end
end

end

function r = roundClamp(x, sz)
r = round(x);
if r < 1
    r = 1;
elseif r > sz
    r = sz;
end
end
