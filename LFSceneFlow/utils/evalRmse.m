function [ rmse ] = evalRmse( flowx, flowy, flowz, gtFlowx, gtFlowy, gtFlowz, borderSize )
%EVALRMSE Evaluate root mean squared error

if nargin < 7
    borderSize = 0;
end
se = (flowx-gtFlowx).^2 + (flowy-gtFlowy).^2 + (flowz-gtFlowz).^2;
se = se(1+borderSize:end-borderSize,1+borderSize:end-borderSize);
rmse = sqrt(mean(se(:)));

end

