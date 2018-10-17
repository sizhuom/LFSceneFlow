function [ mae ] = evalMae( flowx, flowy, flowz, gtFlowx, gtFlowy, gtFlowz, borderSize )
%EVALMAE Evaluate mean absolute error for X, Y, Z flows separately

if nargin < 7
    borderSize = 0;
end
xMap = abs(flowx-gtFlowx);
yMap = abs(flowy-gtFlowy);
zMap = abs(flowz-gtFlowz);

mae = zeros(1,3);
xMap = xMap(1+borderSize:end-borderSize,1+borderSize:end-borderSize);
mae(1) = mean(xMap(:));
yMap = yMap(1+borderSize:end-borderSize,1+borderSize:end-borderSize);
mae(2) = mean(yMap(:));
zMap = zMap(1+borderSize:end-borderSize,1+borderSize:end-borderSize);
mae(3) = mean(zMap(:));


end

