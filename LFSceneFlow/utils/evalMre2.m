function [ mre ] = evalMre2( flowx, flowy, flowz, gtFlowx, gtFlowy, gtFlowz, borderSize )
%EVALMRE2 Evaluate mean relative error (relative to maximum magnitude)

if nargin < 7
    borderSize = 0;
end
xMap = flowx-gtFlowx;
yMap = flowy-gtFlowy;
zMap = flowz-gtFlowz;
errMap = sqrt(xMap.^2+yMap.^2+zMap.^2);

magMap = sqrt(gtFlowx.^2+gtFlowy.^2+gtFlowz.^2);
errMap = errMap ./ max(magMap(:));
errMap = errMap(1+borderSize:end-borderSize,1+borderSize:end-borderSize);
mre = mean(errMap(:));


end

