function [ rmse ] = evalRmseXYZ( flowx, flowy, flowz, gtFlowx, gtFlowy, gtFlowz, borderSize )
%EVALRMSEXYZ Evaluate root mean squared error for X,Y,Z separately

if nargin < 7
    borderSize = 0;
end
rmse = zeros(1,3);

se = (flowx-gtFlowx).^2;
se = se(1+borderSize:end-borderSize,1+borderSize:end-borderSize);
rmse(1) = sqrt(mean(se(:)));

se = (flowy-gtFlowy).^2;
se = se(1+borderSize:end-borderSize,1+borderSize:end-borderSize);
rmse(2) = sqrt(mean(se(:)));

se = (flowz-gtFlowz).^2;
se = se(1+borderSize:end-borderSize,1+borderSize:end-borderSize);
rmse(3) = sqrt(mean(se(:)));

end

