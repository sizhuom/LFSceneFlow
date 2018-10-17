function [ aae ] = evalAae( flowx, flowy, flowz, gtFlowx, gtFlowy, gtFlowz, borderSize )
%EVALAAE Evaluate average angular error

if nargin < 7
    borderSize = 0;
end
ae = acos((flowx.*gtFlowx+flowy.*gtFlowy+flowz.*gtFlowz+1)./...
    sqrt((flowx.^2+flowy.^2+flowz.^2+1).*(gtFlowx.^2+gtFlowy.^2+gtFlowz.^2+1)));
ae = ae / pi * 180;
ae = ae(1+borderSize:end-borderSize,1+borderSize:end-borderSize);
aae = mean(ae(:));

end

