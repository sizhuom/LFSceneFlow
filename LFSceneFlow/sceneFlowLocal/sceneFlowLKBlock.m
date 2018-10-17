function [ flow ] = sceneFlowLKBlock( LF0, LF1, intr, param )
%SCENEFLOWLKBLOCK Compute the scene flow between two light fields

% init
% Process the parameters
sigma = [1 1 0 0];
maxIters = 3;
tempBlend = 0.6;
if (nargin > 3)
    fields = fieldnames(param);
    for i = 1:numel(fields)
        eval([fields{i} '=param.' fields{i} ';']);
    end
end

% pre-filtering with Gaussian
% TODO: combine with the deriv filters
sigma_rad = ceil(sigma*1.5);
LF0 = gaussFilt4D(LF0, sigma, sigma_rad);
LF1 = gaussFilt4D(LF1, sigma, sigma_rad);

[~, Li0, Lj0, Lk0, Ll0] = partialDeriv(LF0);
[~, Li1, Lj1, Lk1, Ll1] = partialDeriv(LF1);

[I,J,K,L] = ndgrid(1:size(LF0,1),1:size(LF0,2),1:size(LF0,3),1:size(LF0,4));
ind = sub2ind(size(LF0), I, J, K, L);
ind = ind(:);
flowv = zeros(3,1);
for iter = 1:maxIters
    flowPrev = repmat(reshape(flowv,[1 1 1 1 3]),...
        size(I));
    flowProj = premultHM(flowPrev,intr,J,I,L,K);
    I2 = I + flowProj(:,:,:,:,2);
    J2 = J + flowProj(:,:,:,:,1);
    K2 = K + flowProj(:,:,:,:,4);
    L2 = L + flowProj(:,:,:,:,3);
    
    warpIm  = interpn(LF1,I2,J2,K2,L2,'cubic');
    Lt      = warpIm(:) - LF0(ind);
    Lt(isnan(Lt)) = 0;
    
    Li = interpn(Li1,I2,J2,K2,L2,'cubic');
    Lj = interpn(Lj1,I2,J2,K2,L2,'cubic');
    Lk = interpn(Lk1,I2,J2,K2,L2,'cubic');
    Ll = interpn(Ll1,I2,J2,K2,L2,'cubic');
    Li = tempBlend*Li + (1-tempBlend)*reshape(Li0(ind),size(Li));
    Lj = tempBlend*Lj + (1-tempBlend)*reshape(Lj0(ind),size(Lj));
    Lk = tempBlend*Lk + (1-tempBlend)*reshape(Lk0(ind),size(Lk));
    Ll = tempBlend*Ll + (1-tempBlend)*reshape(Ll0(ind),size(Ll));
    Li(isnan(Li)) = 0;
    Lj(isnan(Lj)) = 0;
    Lk(isnan(Lk)) = 0;
    Ll(isnan(Ll)) = 0;
    
    P = cat(5, Li, Lj, Lk, Ll);
    LM = postmultHM(P, intr, J2, I2, L2, K2);
    LM = reshape(LM, [], 3);
    
    % solve the linear system
    A = LM;
    b = -Lt;
    flowv = flowv + A\b;
end
flow = flowv;

end

