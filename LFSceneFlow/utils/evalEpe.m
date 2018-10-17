function [ epe ] = evalEpe( flowx, flowy, flowz, gtFlowx, gtFlowy, gtFlowz, borderSize )
%EVALEPE Evaluate endpoint error

if nargin < 7
    borderSize = 0;
end
se = (flowx-gtFlowx).^2 + (flowy-gtFlowy).^2 + (flowz-gtFlowz).^2;
se = se(1+borderSize:end-borderSize,1+borderSize:end-borderSize);
epe = mean(sqrt(se(:)));

end

