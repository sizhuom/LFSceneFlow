function [ resized, Hnew ] = LFResize2D( LF, scale, method, H )
%LFRESIZE2D Resize the u-v dimension of a 4D light field

if (nargin < 3)
    method = 'bicubic';
end

if (size(scale, 2) == 1)
    sy = size(LF, 1);
    sx = size(LF, 2);
    sv = ceil(size(LF, 3) * scale); % round? ceil?
    su = ceil(size(LF, 4) * scale);
    ch = size(LF, 5);
    resized = zeros([sy sx sv su ch]);
    
    for j = 1:sy
        for i = 1:sx
            im = imresize(squeeze(LF(j, i, :, :, :)), scale, method);
            resized(j, i, :, :, :) = im;
        end
    end
    
    if nargout > 1
        dH = [1 0 0 0 0;
            0 1 0 0 0;
            0 0 1/scale 0 1/2-1/(2*scale);
            0 0 0 1/scale 1/2-1/(2*scale);
            0 0 0 0 1];
        Hnew = H * dH;
    end
else
    sy = size(LF, 1);
    sx = size(LF, 2);
    sv = scale(1);
    su = scale(2);
    ch = size(LF, 5);
    resized = zeros([sy sx sv su ch]);
    
    for j = 1:sy
        for i = 1:sx
            im = imresize(squeeze(LF(j, i, :, :, :)), scale, method);
            resized(j, i, :, :, :) = im;
        end
    end
    
    if nargout > 1
        rv = scale(1)/size(LF,3);
        ru = scale(2)/size(LF,4);
        dH = [1 0 0 0 0;
            0 1 0 0 0;
            0 0 1/ru 0 1/2-1/(2*ru);
            0 0 0 1/rv 1/2-1/(2*rv);
            0 0 0 0 1];
        Hnew = H * dH;
    end
end


end

