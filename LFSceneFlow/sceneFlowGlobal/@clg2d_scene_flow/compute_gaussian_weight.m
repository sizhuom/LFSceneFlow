function [ w ] = compute_gaussian_weight( this )
%COMPUTE_GAUSSIAN_WEIGHT Compute Gaussian weight for the oriented window

[Y,X,V,U] = ndgrid(-this.owSize(1):this.owSize(1),...
    -this.owSize(2):this.owSize(2),-this.owSize(3):this.owSize(3),...
    -this.owSize(4):this.owSize(4));

w = exp(-Y.^2/(2*this.owSigma(1)^2)-X.^2/(2*this.owSigma(2)^2)...
    -V.^2/(2*this.owSigma(3)^2)-U.^2/(2*this.owSigma(4)^2));
w = w(:);
w = w / sum(w);

end

