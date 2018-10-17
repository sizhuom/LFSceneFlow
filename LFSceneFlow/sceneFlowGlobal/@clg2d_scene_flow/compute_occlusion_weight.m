function [ w ] = compute_occlusion_weight( this )
%COMPUTE_OCCLUSION_WEIGHT Compute weight due to occlusion detected on EPI

owSize = this.owSize;
alphaMap = this.alphaMap;
SY   = size(this.images, 1);
SX   = size(this.images, 2);
SV   = size(this.images, 3);
SU   = size(this.images, 4);
SOW = prod(2*owSize+1);

CY = (1+SY)/2;
CX = (1+SX)/2;

[winY,winX,winV,winU] = ndgrid(CY-owSize(1):CY+owSize(1),...
    CX-owSize(2):CX+owSize(2),...
    -owSize(3):owSize(3),-owSize(4):owSize(4));
x1 = repmat(reshape(winX,1,1,[]),SV,SU,1);
y1 = repmat(reshape(winY,1,1,[]),SV,SU,1);
u1 = repmat(reshape(winU,1,1,[]),SV,SU,1);
v1 = repmat(reshape(winV,1,1,[]),SV,SU,1);

[cenV,cenU,~] = ndgrid(1:SV,1:SU,1:SOW);
alphaCen = centralSub(alphaMap);
u1 = u1+cenU-bsxfun(@times,alphaCen,x1-CX);
v1 = v1+cenV-bsxfun(@times,alphaCen,y1-CY);
clear winY winX winV winU cenV cenU

invAlpha = interpn(1./alphaMap,y1,x1,v1,u1,'nearest');
cenIdx = sub2ind(2*owSize+1,1+owSize(1),1+owSize(2),1+owSize(3),1+owSize(4));
w = exp(-bsxfun(@minus,invAlpha,invAlpha(:,:,cenIdx)).^2/(2*this.owOccSigma.^2));
w(isnan(w)) = 0;

end

