function [A, b, params, iterative] = flow_operator(this, opq)
%FLOW_OPERATOR   Linear flow operator (equation) for flow estimation
%   [A, b] = FLOW_OPERATOR(THIS, UV, INIT)
%   returns a linear flow operator (equation) of the form A * x = b.  The
%   flow equation is linearized around UV with the initialization INIT
%   (e.g. from a previous pyramid level).
%
%   [A, b, PARAMS, ITER] = FLOW_OPERATOR(...) returns optional parameters
%   PARAMS that are to be passed into a linear equation solver and a flag
%   ITER that indicates whether solving for the flow requires multiple
%   iterations of linearizing.
%
%   This is a member function of the class 'hs_optical_flow'.
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

% Compute spatial and temporal partial derivatives
if this.use_bs
    [It,IX,IY,IZ] = partial_deriv_bs(this.images, this.bs_images, this.intr, opq, ...
        this.interpolation_method);
else
    [It,IX,IY,IZ] = partial_deriv(this.images, this.intr, opq, ...
        this.interpolation_method, this.deriv_filter_xy, this.deriv_filter_uv,...
        this.interp_filter_xy, this.interp_filter_uv);
end
% [It,IX,IY,IZ] = dbgPartialDeriv(this.images, this.intr, opq); % debugging

sz        = [size(IX, 1) size(IX, 2) size(IX, 3) size(IX, 4)];
npixels   = prod(sz);

% Spatial/prior term
if isempty(this.xyweightMap)
    L     = ndlaplacian(4);     % Laplacian operator
    % notice the difference between 'replicate' and 'valid';
    % should not be a problem?
    F     = make_imfilter_mat(L, sz, 'replicate', 'same');
%     F     = make_convn_mat(L, sz, 'valid', 'sameswap');
    
    % Replicate operator for u and v
    M     = [this.lambda*F, sparse(npixels, 2*npixels);
        sparse(npixels, npixels), this.lambda*F, sparse(npixels, npixels);
        sparse(npixels, 2*npixels), this.lambdaZ*F];
else
    S = this.spatial_filters;
    
    FO = sparse(npixels, npixels);
    FP = sparse(npixels, npixels);
    FQ = sparse(npixels, npixels);
    for i = 1:length(S)
        
        FMi = make_convn_mat(S{i}, sz, 'valid', 'sameswap');
        Fi  = FMi';
        W = spdiags(this.xyweightMap(:),0,npixels,npixels);
        FO        = FO+ Fi*W*FMi;
        FP        = FP+ Fi*W*FMi;
        FQ        = FQ+ Fi*W*FMi;
        
    end;
    clear FMi Fi o_ p_ q_
    
    M = [-this.lambda*FO, sparse(npixels, 2*npixels);
        sparse(npixels, npixels), -this.lambda*FP, sparse(npixels, npixels);
        sparse(npixels, 2*npixels), -this.lambdaZ*FQ];
    clear FO FP FQ FR;
end

% For color processing
IX2 = mean(IX.^2, 5);
IY2 = mean(IY.^2, 5);
IZ2 = mean(IZ.^2, 5);
IXY = mean(IX.*IY, 5);
IXZ = mean(IX.*IZ, 5);
IYZ = mean(IY.*IZ, 5);
ItX = mean(It.*IX, 5);
ItY = mean(It.*IY, 5);
ItZ = mean(It.*IZ, 5);

doo   = spdiags(IX2(:), 0, npixels, npixels);
dpp   = spdiags(IY2(:), 0, npixels, npixels);
dqq   = spdiags(IZ2(:), 0, npixels, npixels);
dop   = spdiags(IXY(:), 0, npixels, npixels);
doq   = spdiags(IXZ(:), 0, npixels, npixels);
dpq   = spdiags(IYZ(:), 0, npixels, npixels);

% Compute the operator
A     = [doo dop doq; dop dpp dpq; doq dpq dqq]/this.sigmaD2 - M/this.sigmaS2;
b     = M*opq(:)/this.sigmaS2 - [ItX(:); ItY(:); ItZ(:)]/this.sigmaD2;

% No auxiliary parameters
params    = [];
iterative = true;

% If the non-linear weights are non-uniform, we have to iterate
%   if (max(pp_s(:)) - min(pp_s(:)) < 1E-6 && ...
%       max(pp_d(:)) - min(pp_d(:)) < 1E-6)
%     iterative = false;
%   else
%     iterative = true;
%   end