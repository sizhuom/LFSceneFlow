function [ flow ] = sceneFlowLK( LF0, LF1, intr, param )
%SCENEFLOW Compute the scene flow between two light fields

% init
% Process the parameters
sigma = [0 0 0 0];
halfWindow = [1 1 1 1];
useGaussWeight = true;
useShearedGauss = false;
occSigma = 0;
bilSigma = 0;
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

LFs = cat(5, LF0, LF1);
flow = zeros(size(LF0,3), size(LF0,4), 3);

[Lt, Li, Lj, Lk, Ll] = partialDeriv(LFs);
NPixels = numel(Lt);
Lt = reshape(Lt(:,:,:,:), [NPixels, 1]);
Li = reshape(Li(:,:,:,:), [NPixels, 1]);
Lj = reshape(Lj(:,:,:,:), [NPixels, 1]);
Lk = reshape(Lk(:,:,:,:), [NPixels, 1]);
Ll = reshape(Ll(:,:,:,:), [NPixels, 1]);

P = cat(5, Li, Lj, Lk, Ll);
LM = postmultHM(P, intr);
LM = reshape(LM, [NPixels 3]);

i = round((1+size(LF0,1))/2);
j = round((1+size(LF0,2))/2);
charCount = 0;
% estimate alpha, if no external alpha map is given
if useShearedGauss || occSigma > 0
    if ~exist('alphaMap', 'var')
        alphaMap = estimateAlpha(Li,Lj,Lk,Ll,size(LF0));
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
    for k = 1:size(LF0,3)
        fprintf(repmat('\b',1,charCount));
        charCount = fprintf('Computing: %d/%d...\n',...
            k,size(LF0,3));
        for l = 1:size(LF0,4)
            % sheared window
            [I,J,K,L] = ndgrid(-halfWindow(1):+halfWindow(1),...
                -halfWindow(2):+halfWindow(2),...
                -halfWindow(3):+halfWindow(3),...
                -halfWindow(4):+halfWindow(4));
            K = K - alphaMap(i,j,k,l) * I;
            L = L - alphaMap(i,j,k,l) * J;
            I = I+i;
            J = J+j;
            K = min(max(round(K)+k,1),size(LF0,3));
            L = min(max(round(L)+l,1),size(LF0,4));
            ind = sub2ind(size(LF0), I,J,K,L);
            ind = ind(:);
            % occlusion term
            W = W0(:);
            if occSigma > 0
                W = W .* exp(-((alphaMap(ind)-alphaMap(i,j,k,l))/occSigma).^2/2);
            end
            % bilateral term
            if bilSigma > 0
                W = W .* exp(-((LF0(ind)-LF0(i,j,k,l))/bilSigma).^2/2);
            end
            % solve the linear system
            A = LM(ind, :);
            b = -Lt(ind);
            ATW = bsxfun(@times,A',W');
            ATA = ATW*A;
            ATb = ATW*b;
            if (rank(ATA) == 3)
                flow(k,l,:) = inv(ATA)*ATb;
            else
                flow(k,l,:) = NaN;
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
    for k = 1:size(LF0,3)
        fprintf(repmat('\b',1,charCount));
        charCount = fprintf('Computing: %d/%d...\n',...
            k,size(LF0,3));
        for l = 1:size(LF0,4)
            % Gaussian window
            [I,J,K,L] = ndgrid(i-halfWindow(1):i+halfWindow(1),...
                j-halfWindow(2):j+halfWindow(2),...
                k-halfWindow(3):k+halfWindow(3),...
                l-halfWindow(4):l+halfWindow(4));
            K = min(max(K,1),size(LF0,3));
            L = min(max(L,1),size(LF0,4));
            ind = sub2ind(size(LF0), I, J, K, L);
            ind = ind(:);
            W = W0(:);
            % occlusion term
            if occSigma > 0
                W = W .* exp(-((alphaMap(ind)-alphaMap(i,j,k,l))/occSigma).^2/2);
            end
            % bilateral term
            if bilSigma > 0
                W = W .* exp(-((LF0(ind)-LF0(i,j,k,l))/bilSigma).^2/2);
            end
            % solve the linear system
            A = LM(ind, :);
            b = -Lt(ind);
            ATW = bsxfun(@times,A',W');
            ATA = ATW*A;
            ATb = ATW*b;
            if (rank(ATA) == 3)
                flow(k,l,:) = inv(ATA)*ATb;
            else
                flow(k,l,:) = NaN;
            end
        end
    end
end

end

