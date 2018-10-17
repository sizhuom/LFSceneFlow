function [ mae ] = evalMFAbsoluteError( flowx, flowy, flowz, gtFlowx, gtFlowy, gtFlowz, borderSize )
%EVALMFABSOLUTEERROR Evaluate mean foreground absolute error for X, Y, Z flows separately

if nargin < 7
    borderSize = 0;
end
xMap = abs(flowx-gtFlowx);
yMap = abs(flowy-gtFlowy);
zMap = abs(flowz-gtFlowz);

mask1 = gtFlowx~=0 | gtFlowy~=0 | gtFlowz~=0;
mask2 = zeros(size(flowx));
mask2(1+borderSize:end-borderSize,1+borderSize:end-borderSize) = 1;
mask = mask1 & mask2;

mae = zeros(1,3);
xMap = xMap(mask);
mae(1) = median(xMap(:));
yMap = yMap(mask);
mae(2) = median(yMap(:));
zMap = zMap(mask);
mae(3) = median(zMap(:));


end

