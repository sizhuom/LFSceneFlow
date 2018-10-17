function [ mfae ] = evalMFAngularError( flowx, flowy, flowz, gtFlowx, gtFlowy, gtFlowz, borderSize )
%EVALMFANGULARERROR Evaluate median of foreground angular error

if nargin < 7
    borderSize = 0;
end
ae = acos((flowx.*gtFlowx+flowy.*gtFlowy+flowz.*gtFlowz+1)./...
    sqrt((flowx.^2+flowy.^2+flowz.^2+1).*(gtFlowx.^2+gtFlowy.^2+gtFlowz.^2+1)));
ae = ae / pi * 180;

mask1 = gtFlowx~=0 | gtFlowy~=0 | gtFlowz~=0;
mask2 = zeros(size(flowx));
mask2(1+borderSize:end-borderSize,1+borderSize:end-borderSize) = 1;
mask = mask1 & mask2;
ae = ae(mask);
mfae = median(ae(:));

end

