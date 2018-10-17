function u = compute_flow_base(this, u)

%COMPUTE_FLOW_BASE   Base function for computing flow field
%   UV = COMPUTE_FLOW_BASE(THIS, INIT) computes the flow field UV with
%   algorithm THIS and the initialization INIT.
%
%   This is a member function of the class 'hs_optical_flow'.
%
%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2007-11-30 $
%   $Revision: $
%
% Copyright 2007-2010, Brown University, Providence, RI. USA
%
%                          All Rights Reserved
%
% All commercial use of this software, whether direct or indirect, is
% strictly prohibited including, without limitation, incorporation into in
% a commercial product, use in a commercial service, or production of other
% artifacts for commercial purposes.
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for research purposes is hereby granted without fee,
% provided that the above copyright notice appears in all copies and that
% both that copyright notice and this permission notice appear in
% supporting documentation, and that the name of the author and Brown
% University not be used in advertising or publicity pertaining to
% distribution of the software without specific, written prior permission.
%
% For commercial uses contact the Technology Venture Office of Brown University
%
% THE AUTHOR AND BROWN UNIVERSITY DISCLAIM ALL WARRANTIES WITH REGARD TO
% THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
% FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR OR
% BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
% DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
% PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
% ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
% THIS SOFTWARE.

GRAD_IS_ZERO = 1e-10;

sz = [size(u,1) size(u,2) size(u,3) size(u,4)];
l_t = this.lambda * this.theta;
lz_t = this.lambdaZ * this.theta;

w = this.xyweightMap;

% initialize p
p1 = zeros([sz 4]);
p2 = zeros([sz 4]);
p3 = zeros([sz 4]);

% Iterate flow computation
for i = 1:this.max_warping_iters
    u0 = u;
    
    if this.display
        fprintf('energy of solution \t%3.3e\n', -evaluate_log_posterior(this, u));
    end;
    
    % Compute warped image derivatives\
    [It, IX, IY, IZ] = partial_deriv(this.images, this.intr, u, this.interpolation_method);
    
    % store |grad(I1)|^2
    grad = IX.^2 + IY.^2 + IZ.^2;
    
    % compute the constant part of the rho function
    rho_c = It - IX .* u(:,:,:,:,1) - IY .* u(:,:,:,:,2) - IZ .* u(:,:,:,:,3);
    
    n = 0;
    error = Inf;
    while (error > this.epsilon^2 && n < this.max_iters)
        n = n + 1;
        % estimate the values of the variable (v1, v2)
        % (thresholding operator TH)
        rho = rho_c + IX .* u(:,:,:,:,1) + IY .* u(:,:,:,:,2) + IZ .* u(:,:,:,:,3);
        
        d1 = zeros(sz);
        d2 = zeros(sz);
        d3 = zeros(sz);
        
        mask1 = rho < -this.theta * grad;
        mask2 = rho > this.theta * grad;
        mask3 = ~(mask1 | mask2) & (grad >= GRAD_IS_ZERO);
        
        d1(mask1) = this.theta * IX(mask1);
        d2(mask1) = this.theta * IY(mask1);
        d3(mask1) = this.theta * IZ(mask1);
        d1(mask2) = -this.theta * IX(mask2);
        d2(mask2) = -this.theta * IY(mask2);
        d3(mask2) = -this.theta * IZ(mask2);
        fi = -rho ./ grad;
        d1(mask3) = fi(mask3) .* IX(mask3);
        d2(mask3) = fi(mask3) .* IY(mask3);
        d3(mask3) = fi(mask3) .* IZ(mask3);
        
        v(:,:,:,:,1) = u(:,:,:,:,1) + d1;
        v(:,:,:,:,2) = u(:,:,:,:,2) + d2;
        v(:,:,:,:,3) = u(:,:,:,:,3) + d3;
        
        % compute the divergence of the dual variable (p1, p2, p3)
        if isempty(w)
            div_p1 = divergence(p1);
            div_p2 = divergence(p2);
            div_p3 = divergence(p3);
        else
            div_p1 = divergence(bsxfun(@times,p1,w));
            div_p2 = divergence(bsxfun(@times,p2,w));
            div_p3 = divergence(bsxfun(@times,p3,w));
        end
        
        % estimate the values of the optical flow (u1, u2)
        uk = u;
        u = v + cat(5, l_t*div_p1, l_t*div_p2, lz_t*div_p3);
        error = u - uk;
        error = sum(error(:).^2);
        
        % compute the gradient of the optical flow (Du1, Du2)
        [u1x, u1y, u1u, u1v] = forward_gradient(u(:,:,:,:,1));
        [u2x, u2y, u2u, u2v] = forward_gradient(u(:,:,:,:,2));
        [u3x, u3y, u3u, u3v] = forward_gradient(u(:,:,:,:,3));
        if ~isempty(w)
            u1x = w .* u1x;
            u1y = w .* u1y;
            u1u = w .* u1u;
            u1v = w .* u1v;
            u2x = w .* u2x;
            u2y = w .* u2y;
            u2u = w .* u2u;
            u2v = w .* u2v;
            u3x = w .* u3x;
            u3y = w .* u3y;
            u3u = w .* u3u;
            u3v = w .* u3v;
        end
        
        % estimate the values of the dual variable
        taut = this.tau / l_t;
        tautz = this.tau / lz_t;
        g1 = sqrt(u1x.^2 + u1y.^2 + u1u.^2 + u1v.^2);
        g2 = sqrt(u2x.^2 + u2y.^2 + u2u.^2 + u2v.^2);
        g3 = sqrt(u3x.^2 + u3y.^2 + u3u.^2 + u3v.^2);
        ng1 = 1 + taut * g1;
        ng2 = 1 + taut * g2;
        ng3 = 1 + tautz * g3;
        p1 = bsxfun(@rdivide, p1 + taut * cat(5,u1x,u1y,u1u,u1v), ng1);
        p2 = bsxfun(@rdivide, p2 + taut * cat(5,u2x,u2y,u2u,u2v), ng2);
        p3 = bsxfun(@rdivide, p3 + tautz * cat(5,u3x,u3y,u3u,u3v), ng3);
    end
    
    x = u - u0;
    % Print status information
    if this.display
        disp(['--Iteration: ', num2str(i), '    (', ...
            num2str(norm(x(:))), ')'])
    end;
    
    % debug: save result after each iteration
    if this.debug
        save(findNextTmp(this.debug_dir,'result-*.mat'), 'opq');
        fprintf('iteration data saved.\n');
    end
    
    % Terminate iteration early if flow doesn't change substantially
    if (length(this.lambda) == 1 && norm(x(:)) < 1E-3)
        break
    end
    
    % If limiting the incremental flow to [-1, 1] is requested, do so
    if (this.limit_update)
        x(x > 1)  = 1;
        x(x < -1) = -1;
    end
    
    u = u0 + x;
    
    % Perform median filtering to remove outliers
    if ~isempty(this.median_filter_size)
        for m = 1:this.mf_iter; % extensive MF filtering
            LFMedfilt2(opq, this.median_filter_size);
        end;
    end;
    
    % debug: save result after each iteration
%     opq = u;
%     save(findNextTmp('~/Documents','result-*.mat'), 'opq');
    
    % Terminate early if the flow_operator doesn't require multiple
    % interations
    %     if (~iterative)
    %       break;
    %     end
    
end
