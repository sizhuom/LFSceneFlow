function [A, b, params, iterative] = flow_operator(this, opq, dopq, It, IX, IY, IZ)

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

sz        = [size(IX,1) size(IX,2) size(IX,3) size(IX,4)];
npixels   = prod(sz);

% spatial term
S = this.spatial_filters;

FO = sparse(npixels, npixels);
FP = sparse(npixels, npixels);
FQ = sparse(npixels, npixels);
for i = 1:length(S)
    
    FMi = make_convn_mat(S{i}, sz, 'valid', 'sameswap');
    Fi  = FMi';
    
    % Use flow increment to update the nonlinearity
    o_        = FMi*reshape(opq(:, :, :, :, 1)+dopq(:, :, :, :, 1), [npixels 1]);
    p_        = FMi*reshape(opq(:, :, :, :, 2)+dopq(:, :, :, :, 2), [npixels 1]);
    q_        = FMi*reshape(opq(:, :, :, :, 3)+dopq(:, :, :, :, 3), [npixels 1]);
    
    if isa(this.rho_spatial_u{i}, 'robust_function')
        pp_so     = deriv_over_x(this.rho_spatial_u{i}, o_);
        pp_sp     = deriv_over_x(this.rho_spatial_v{i}, p_);
        pp_sq     = deriv_over_x(this.rho_spatial_w{i}, q_);
    elseif isa(this.rho_spatial_u{i}, 'gsm_density')
        pp_so     = -evaluate_log_grad_over_x(this.rho_spatial_u{i}, o_')';
        pp_sp     = -evaluate_log_grad_over_x(this.rho_spatial_v{i}, p_')';
        pp_sq     = -evaluate_log_grad_over_x(this.rho_spatial_w{i}, q_')';
    else
        error('evaluate_log_posterior: unknown rho function!');
    end;
    
    if isempty(this.weightMap)
        FO        = FO+ Fi*spdiags(pp_so, 0, npixels, npixels)*FMi;
        FP        = FP+ Fi*spdiags(pp_sp, 0, npixels, npixels)*FMi;
        FQ        = FQ+ Fi*spdiags(pp_sq, 0, npixels, npixels)*FMi;
    else
        W = spdiags(this.weightMap(:),0,npixels,npixels);
        FO        = FO+ Fi*W*spdiags(pp_so, 0, npixels, npixels)*FMi;
        FP        = FP+ Fi*W*spdiags(pp_sp, 0, npixels, npixels)*FMi;
        FQ        = FQ+ Fi*W*spdiags(pp_sq, 0, npixels, npixels)*FMi;
    end
    
end;
clear FMi Fi o_ p_ q_

M = [-this.lambda*FO, sparse(npixels, 2*npixels);
    sparse(npixels, npixels), -this.lambda*FP, sparse(npixels, npixels);
    sparse(npixels, 2*npixels), -this.lambdaZ*FQ];
clear FO FP FQ FR;

% Data term
IX2 = IX.^2;
IY2 = IY.^2;
IZ2 = IZ.^2;
IXY = IX.*IY;
IXZ = IX.*IZ;
IYZ = IY.*IZ;
ItX = It.*IX;
ItY = It.*IY;
ItZ = It.*IZ;

% Filter the structure tensors
if this.useCLG
    IX2 = gaussFilt4D(IX2,this.sigmaTS);
    IY2 = gaussFilt4D(IY2,this.sigmaTS);
    IZ2 = gaussFilt4D(IZ2,this.sigmaTS);
    IXY = gaussFilt4D(IXY,this.sigmaTS);
    IXZ = gaussFilt4D(IXZ,this.sigmaTS);
    IYZ = gaussFilt4D(IYZ,this.sigmaTS);
    ItX = gaussFilt4D(ItX,this.sigmaTS);
    ItY = gaussFilt4D(ItY,this.sigmaTS);
    ItZ = gaussFilt4D(ItZ,this.sigmaTS);
end

% Perform linearization - note the change in It
It = It + IX.*repmat(dopq(:,:,:,:,1), [1 1 1 1 size(It,5)]) ...
    + IY.*repmat(dopq(:,:,:,:,2), [1 1 1 1 size(It,5)]) ...
    + IZ.*repmat(dopq(:,:,:,:,3), [1 1 1 1 size(It,5)]);


if isa(this.rho_data, 'robust_function')
    pp_d  = deriv_over_x(this.rho_data, It(:));
elseif isa(this.rho_data, 'gsm_density')
    pp_d = -evaluate_log_grad_over_x(this.rho_data, It(:)')';
else
    error('flow_operator: unknown rho function!');
end;

tmp = pp_d.*IX2(:);
doo = spdiags(tmp, 0, npixels, npixels);
tmp = pp_d.*IY2(:);
dpp = spdiags(tmp, 0, npixels, npixels);
tmp = pp_d.*IZ2(:);
dqq = spdiags(tmp, 0, npixels, npixels);
tmp = pp_d.*IXY(:);
ddop = spdiags(tmp, 0, npixels, npixels);
tmp = pp_d.*IXZ(:);
ddoq = spdiags(tmp, 0, npixels, npixels);
tmp = pp_d.*IYZ(:);
ddpq = spdiags(tmp, 0, npixels, npixels);
clear IX2 IY2 IZ2 IXY IXZ IYZ tmp;

A = [doo ddop ddoq;
    ddop dpp ddpq;
    ddoq ddpq dqq];
clear doo ddop ddoq dpp ddpq dqq;
A = A - M;

% right hand side
b =  M * opq(:) - ...
    [pp_d.*ItX(:); pp_d.*ItY(:); pp_d.*ItZ(:)];

A = A / max(pp_d) / 100; % scale down both A and b, for some reason it helps
b = b / max(pp_d) / 100;

% No auxiliary parameters
params    = [];

% save('cc.mat');
% If the non-linear weights are non-uniform, do more linearization
if (max(pp_so(:)) - min(pp_so(:)) < 1E-6 && ...
        max(pp_sp(:)) - min(pp_sp(:)) < 1E-6 && ...
        max(pp_sq(:)) - min(pp_sq(:)) < 1E-6 && ...
        max(pp_d(:)) - min(pp_d(:)) < 1E-6)
    iterative = false;
else
    iterative = true;
end