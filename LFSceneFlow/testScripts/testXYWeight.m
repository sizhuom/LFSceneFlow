load('/media/sizhuo/Seagate/LFData/Images/Real/1017-card-60-bg/result-1017-CCXY-sor-py-halfdown-a9/result-1017-CCXY-sor-py-halfdown-a9-01.mat');
xyflow = flow;
% load('/media/sizhuo/Seagate/LFData/Images/Real/1002-card-60/mask0002.mat');
% xyflow = cat(5,double(mask),zeros(size(mask)));
sigma = 0.3;
[~,Xx,Xy,Xu,Xv] = partialDeriv(xyflow(:,:,:,:,1));
[~,Yx,Yy,Yu,Yv] = partialDeriv(xyflow(:,:,:,:,2));

gradSq = Xx.^2+Xy.^2+Xu.^2+Xv.^2+Yx.^2+Yy.^2+Yu.^2+Yv.^2;
gradSq = gradSq / max(gradSq(:));
weight = exp(-gradSq / (sigma^2));
se = strel(ones(5,5,5,5));
wr = imerode(weight, se);

% gradL1 = abs(Xx)+abs(Xy)+abs(Xu)+abs(Xv)+abs(Yx)+abs(Yy)+abs(Yu)+abs(Yv);
% gradL1 = gradL1 / max(gradL1(:));
% weight = exp(-gradL1.^2 / (sigma^2));
% se = strel(ones(5,5,5,5));
% wr = imerode(weight, se);

% flowNorm = sqrt(xyflow(:,:,:,:,1).^2 + xyflow(:,:,:,:,2).^2);
% thresh = 0.2;
% flowNorm(flowNorm<thresh) = 0;
% flowNorm(flowNorm>=thresh) = 1;
% [~,Fx,Fy,Fu,Fv] = partialDeriv(flowNorm);
% flowGradSq = Fx.^2+Fy.^2+Fu.^2+Fv.^2;
% flowGradSq = flowGradSq / max(flowGradSq(:));
% weight = exp(-flowGradSq / (sigma^2));
% se = strel(ones(5,5,5,5));
% wr = imerode(weight, se);