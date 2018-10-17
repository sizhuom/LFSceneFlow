function [ mre ] = evalMre( flowx, flowy, flowz, gtFlowx, gtFlowy, gtFlowz, borderSize )
%EVALMRE Evaluate mean relative error

if nargin < 7
    borderSize = 0;
end
xMap = flowx-gtFlowx;
yMap = flowy-gtFlowy;
zMap = flowz-gtFlowz;
errMap = sqrt(xMap.^2+yMap.^2+zMap.^2);

magMap = sqrt(gtFlowx.^2+gtFlowy.^2+gtFlowz.^2);
errMap = errMap ./ magMap;
errMap = errMap(1+borderSize:end-borderSize,1+borderSize:end-borderSize);
errMap = errMap(~isnan(errMap));
mre = mean(errMap(:));


end

