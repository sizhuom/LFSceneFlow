function [ alpha ] = estimateAlpha( Li, Lj, Lk, Ll, sz, param )
%ESTIMATEALPHA Estimate the alpha (effective depth) map using EPI gradients
% alpha = Lx/Lu, -1/alpha is the slope in the EPI (x-u plot)

thresh = 5;
sigma = [0.5 0.5 2 2];

if (nargin > 5)
    fields = fieldnames(param);
    for i = 1:numel(fields)
        eval([fields{i} '=param.' fields{i} ';']);
    end
end

% Estimate by computing the EPI gradients
alpha = (Li ./ Lk + Lj ./ Ll) / 2;
alpha = reshape(alpha,sz);

% Post-filtering
gaussFilt4D(alpha, sigma);

% Remove nan and large values
alpha(isnan(alpha)) = 0;
alpha(alpha<-1 | alpha>thresh) = 0;

end

