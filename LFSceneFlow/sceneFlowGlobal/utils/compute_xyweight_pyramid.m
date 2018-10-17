function [ P ] = compute_xyweight_pyramid( xyflow, sigma, r, smooth_sigma, nL, ratio )
%COMPUTE_XYWEIGHT_PYRAMID Compute a pyramid of weight map based on xy flow

P   = cell(nL,1);
if ~isempty(xyflow)
    P{1}= computeWeightMap(xyflow, sigma, r);
    tmp = P{1};
    
    for m = 2:nL
        % Gaussian filtering
        tmp = gaussFilt4D(tmp, smooth_sigma*ones(4));
        sz1 = [size(tmp,1) size(tmp,2) size(tmp,3) size(tmp,4)];
        sz2 = floor(sz1 * ratio);
        tmp2 = zeros(sz2);
        
        % Compute the intrinsic matrix for the new level
        a = 1 / ratio;
        b = (sz1+1)/2 - (sz2+1)/2*a;
        b = [b(2) b(1) b(4) b(3)];
        Hinc = [
            eye(4)*a b(:);
            zeros(1,4) 1;
            ];
        
        % Subsample the LF
        [y,x,v,u] = ndgrid(1:sz2(1),1:sz2(2),1:sz2(3),1:sz2(4));
        xyuv = [x(:)'; y(:)'; u(:)'; v(:)'; ones(1,numel(x))];
        xyuv = Hinc * xyuv;
        x = reshape(xyuv(1,:), sz2);
        y = reshape(xyuv(2,:), sz2);
        u = reshape(xyuv(3,:), sz2);
        v = reshape(xyuv(4,:), sz2);
        tmp2(:,:,:,:) = interpn(tmp(:,:,:,:),y,x,v,u,'linear');
        
        P{m} = tmp2;
        tmp = tmp2;
    end
else
    for m = 1:nL
        P{m} = [];
    end
end


end

function weight = computeWeightMap(xyflow, sigma, r)

[~,Xx,Xy,Xu,Xv] = partialDeriv(xyflow(:,:,:,:,1));
[~,Yx,Yy,Yu,Yv] = partialDeriv(xyflow(:,:,:,:,2));

gradSq = Xx.^2+Xy.^2+Xu.^2+Xv.^2+Yx.^2+Yy.^2+Yu.^2+Yv.^2;
gradSq = gradSq / max(gradSq(:));
weight = exp(-gradSq / (sigma^2));
weight = imerode(weight, ones(r,r,r,r));

end