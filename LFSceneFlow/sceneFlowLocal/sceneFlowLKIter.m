function [ flow ] = sceneFlowLKIter( LF0, LF1, intr, param )
%SCENEFLOWLKITER Compute the scene flow between two light fields

% init
% Process the parameters
sigma = [0 0 0 0];
halfWindow = [1 1 1 1];
useGaussWeight = true;
useShearedGauss = false;
occSigma = 0;
bilSigma = 0;
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

% pad image to simplify computation
LF0p = padarray(LF0, [0,0,halfWindow(3),halfWindow(4)], 'replicate');
LF1p = padarray(LF1, [0,0,halfWindow(3),halfWindow(4)], 'replicate');
Hinc = [
    1 0 0 0 0;
    0 1 0 0 0;
    0 0 1 0 -halfWindow(4);
    0 0 0 1 -halfWindow(3);
    0 0 0 0 1 ];
intrp = intr;
intrp.sz = size(LF0p);
intrp.H = intrp.H * Hinc;

flow = zeros(size(LF0,3), size(LF0,4), 3);

[~, Li0, Lj0, Lk0, Ll0] = partialDeriv(LF0p);
[~, Li1, Lj1, Lk1, Ll1] = partialDeriv(LF1p);

i = round((1+size(LF0,1))/2);
j = round((1+size(LF0,2))/2);
charCount = 0;
% estimate alpha, if no external alpha map is given
if useShearedGauss || occSigma > 0
    if ~exist('alphaMap', 'var')
        alphaMap = estimateAlpha(Li0,Lj0,Lk0,Ll0,size(LF0));
        alphaMap = padarray(alphaMap, [0,0,halfWindow(3),halfWindow(4)], 'replicate');
    end
end
if useShearedGauss
    windowSigma = (halfWindow+0.5)./3;
    % precompute sheared window weights: ignore sub-pixel offsets
    [I,J,K,L] = ndgrid(-halfWindow(1):+halfWindow(1),...
        -halfWindow(2):+halfWindow(2),...
        -halfWindow(3):+halfWindow(3),...
        -halfWindow(4):+halfWindow(4));
    windowSigma = (halfWindow+0.5)./3;
    W0 = exp(-(I/windowSigma(1)).^2/2-(J/windowSigma(2)).^2/2-(K/windowSigma(3)).^2/2-(L/windowSigma(4)).^2/2);
    for k = 1+halfWindow(3):size(LF0,3)+halfWindow(3)
        fprintf(repmat('\b',1,charCount));
        charCount = fprintf('Computing: %d/%d...\n',...
            k-halfWindow(3),size(LF0,3));
        for l = 1+halfWindow(4):size(LF0,4)+halfWindow(4)
            % sheared window
            [I,J,K,L] = ndgrid(-halfWindow(1):+halfWindow(1),...
                -halfWindow(2):+halfWindow(2),...
                -halfWindow(3):+halfWindow(3),...
                -halfWindow(4):+halfWindow(4));
            K = K - alphaMap(i,j,k,l) * I;
            L = L - alphaMap(i,j,k,l) * J;
            I = I+i;
            J = J+j;
            K = min(max(round(K)+k,1),size(LF0p,3));
            L = min(max(round(L)+l,1),size(LF0p,4));
            ind = sub2ind(size(LF0p), I,J,K,L);
            ind = ind(:);
            % occlusion term
            W = W0(:);
            if occSigma > 0
                W = W .* exp(-((alphaMap(ind)-alphaMap(i,j,k,l))/occSigma).^2/2);
            end
            % bilateral term
            if bilSigma > 0
                W = W .* exp(-((LF0p(ind)-LF0p(i,j,k,l))/bilSigma).^2/2);
            end
            % solve the linear system
            A = LM(ind, :);
            b = -Lt(ind);
            ATW = bsxfun(@times,A',W');
            ATA = ATW*A;
            ATb = ATW*b;
            if (rank(ATA) == 3)
                flow(k-halfWindow(3),l-halfWindow(4),:) = inv(ATA)*ATb;
            else
                flow(k-halfWindow(3),l-halfWindow(4),:) = NaN;
            end
        end
    end
else
    if useGaussWeight
        [I,J,K,L] = ndgrid(-halfWindow(1):+halfWindow(1),...
            -halfWindow(2):+halfWindow(2),...
            -halfWindow(3):+halfWindow(3),...
            -halfWindow(4):+halfWindow(4));
        windowSigma = (halfWindow+0.5)./3;
        W0 = exp(-(I/windowSigma(1)).^2/2-(J/windowSigma(2)).^2/2-(K/windowSigma(3)).^2/2-(L/windowSigma(4)).^2/2);
    else
        W0 = ones(size(2*halfWindow+1));
    end
    for k = 1+halfWindow(3):size(LF0,3)+halfWindow(3)
        fprintf(repmat('\b',1,charCount));
        charCount = fprintf('Computing: %d/%d...\n',...
            k-halfWindow(3),size(LF0,3));
        for l = 1+halfWindow(4):size(LF0,4)+halfWindow(4)
            % Gaussian window
            [I,J,K,L] = ndgrid(i-halfWindow(1):i+halfWindow(1),...
                j-halfWindow(2):j+halfWindow(2),...
                k-halfWindow(3):k+halfWindow(3),...
                l-halfWindow(4):l+halfWindow(4));
%             K = min(max(K,1),size(LF0,3)); % ?
%             L = min(max(L,1),size(LF0,4));
            ind = sub2ind(size(LF0p), I, J, K, L);
            ind = ind(:);
            flowv = zeros(3,1);
            for iter = 1:maxIters
                if isnan(flowv)
                    break;
                end
                flowPrev = repmat(reshape(flowv,[1 1 1 1 3]),...
                    size(I));
                flowProj = premultHM(flowPrev,intrp,J,I,L,K);
                I2 = I + flowProj(:,:,:,:,2);
                J2 = J + flowProj(:,:,:,:,1);
                K2 = K + flowProj(:,:,:,:,4);
                L2 = L + flowProj(:,:,:,:,3);
                
                warpIm  = interpn(LF1p,I2,J2,K2,L2,'cubic');
                Lt      = warpIm(:) - LF0p(ind);
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
                LM = postmultHM(P, intrp, J2, I2, L2, K2);
                LM = reshape(LM, [], 3);
                
                W = W0(:);
                % occlusion term
                if occSigma > 0
                    alphas = interpn(alphaMap,I2,J2,K2,L2,'cubic');
                    W = W .* exp(-((alphas(:)-alphaMap(i,j,k,l))/occSigma).^2/2);
                    W(isnan(W)) = 0;
                end
                % bilateral term
                if bilSigma > 0
                    values = interpn(LF0p,I2,J2,K2,L2,'cubic');
                    W = W .* exp(-((values(:)-LF0p(i,j,k,l))/bilSigma).^2/2);
                    W(isnan(W)) = 0;
                end
                % solve the linear system
                A = LM;
                b = -Lt;
                ATW = bsxfun(@times,A',W');
                ATA = ATW*A;
                ATb = ATW*b;
                if (rank(ATA) == 3)
                    flowv = flowv + inv(ATA)*ATb;
                else
                    flowv = NaN;
                    break;
                end
            end
            flow(k-halfWindow(3),l-halfWindow(4),:) = flowv;
        end
    end
end

end

