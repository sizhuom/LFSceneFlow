function L = evaluate_log_posterior(this, opq)
%EVALUATE_LOG_POSTERIOR computes the log-posterior (negative energy) of the
%   flow fields UV
%   Actually only proportional to the log posterior since the variance of neither the
%   spatial nor the data terms is considered
%
% This is a member function of the class 'classic_nl_optical_flow'.
%
% Authors: Deqing Sun, Department of Computer Science, Brown University
% Contact: dqsun@cs.brown.edu
% $Date: $
% $Revision: $
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

% Spatial term
S = this.spatial_filters;
p = 0;

for i = 1:length(S)
    
%     o_ = conv2(opq(:,:,1), S{i}, 'valid');
%     p_ = conv2(opq(:,:,2), S{i}, 'valid');
%     q_ = conv2(opq(:,:,3), S{i}, 'valid');
    o_ = conv2(opq(:,:,1), S{i}, 'same');
    p_ = conv2(opq(:,:,2), S{i}, 'same');
    q_ = conv2(opq(:,:,3), S{i}, 'same');
    
    if isempty(this.weightMap)
        if isa(this.rho_spatial_u{i}, 'robust_function')
            
            p = p - this.lambda*sum(evaluate(this.rho_spatial_u{i}, o_(:)))...
                - this.lambda*sum(evaluate(this.rho_spatial_v{i}, p_(:)))...
                - this.lambdaZ*sum(evaluate(this.rho_spatial_w{i}, q_(:)));
            
        elseif isa(this.rho_spatial_u{i}, 'gsm_density')
            
            p   = p + this.lambda*sum(evaluate_log(this.rho_spatial_u{i}, o_(:)'))...
                + this.lambda*sum(evaluate_log(this.rho_spatial_v{i}, p_(:)'))...
                + this.lambdaZ*sum(evaluate_log(this.rho_spatial_w{i}, q_(:)'));
            
        else
            error('evaluate_log_posterior: unknown rho function!');
        end;
    else
        w1 = this.weightMap(:);
        if isa(this.rho_spatial_u{i}, 'robust_function')
            
            p = p - this.lambda*sum(w1.*evaluate(this.rho_spatial_u{i}, o_(:)))...
                - this.lambda*sum(w1.*evaluate(this.rho_spatial_v{i}, p_(:)))...
                - this.lambdaZ*sum(w1.*evaluate(this.rho_spatial_w{i}, q_(:)));
            
        elseif isa(this.rho_spatial_u{i}, 'gsm_density')
            
            p   = p + this.lambda*sum(w1.*evaluate_log(this.rho_spatial_u{i}, o_(:)'))...
                + this.lambda*sum(w1.*evaluate_log(this.rho_spatial_v{i}, p_(:)'))...
                + this.lambdaZ*sum(w1.*evaluate_log(this.rho_spatial_w{i}, q_(:)'));
            
        else
            error('evaluate_log_posterior: unknown rho function!');
        end;
    end;
end;

% likelihood
It = partial_deriv_clg2d(this.images, this.intr, opq, ...
    this.owSize, this.alphaCen, this.interpolation_method,...
    this.deriv_filter_xy, this.deriv_filter_uv, ...
    this.interp_filter_xy, this.interp_filter_uv,...
    0.6);

if isa(this.rho_data, 'robust_function')
    
    l   = -evaluate(this.rho_data, It(:));
    
elseif isa(this.rho_data, 'gsm_density')
    
    l   = evaluate_log(this.rho_data, It(:)');
    
else
    error('evaluate_log_posterior: unknown rho function!');
end;

% w2 = ones(1,1,size(It,3)) / size(It,3);
w2 = reshape(compute_gaussian_weight(this),1,1,[]);
if this.owOccSigma > 0
    wocc = compute_occlusion_weight(this);
    w2 = bsxfun(@times,w2,wocc);
    w2 = bsxfun(@rdivide,w2,sum(w2,3));
end

l = reshape(l,size(It));
l = bsxfun(@times,l,w2);
l = sum(l(:));

% zero-motion bias term
zm = sum(reshape(opq.^2,[],1)) * this.mu;

L = p + l + zm;

if this.display
    fprintf('spatial\t%3.2e\tdata\t%3.2e\tzero-motion\t%3.2e\n', p, l, zm);
end;
