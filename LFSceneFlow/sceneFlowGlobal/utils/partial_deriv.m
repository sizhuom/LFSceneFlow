function [It, IX, IY, IZ] = partial_deriv(images, intr, flow_prev, ...
    interpolation_method, deriv_filter_xy, deriv_filter_uv, ...
    interp_filter_xy, interp_filter_uv, b, mask)

%PARTIAL_DERIV   Spatio-temporal derivatives
%   P = PARTIAL_DERIV(IMAGES, INIT) computes the spatio-temporal derivatives
%
%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2007-11-30 $
%   $Revision $
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

if nargin < 3
    flow_prev = zeros([intr.sz 3]);
end

if nargin < 4
    interpolation_method = 'bi-linear';
end;

% h1: deriv_filter for x,y
if nargin >= 5 && ~isempty(deriv_filter_xy)
    h1 = deriv_filter_xy;
else
    h1 = [-1 0 1] / 2;
end
% h2: deriv_filter for u,v
if nargin >= 6 && ~isempty(deriv_filter_uv)
    h2 = deriv_filter_uv;
else
    h2 = [1 -8 0 8 -1]/12; % used in Wedel etal "improved TV L1"
end;
% k1: interp_filter for x,y
if nargin >= 7 && ~isempty(interp_filter_xy)
    k1 = interp_filter_xy;
else
    k1 = [1];
end;
% k2: interp_filter for u,v
if nargin >= 8 && ~isempty(interp_filter_uv)
    k2 = interp_filter_uv;
else
    k2 = [1];
end

if nargin < 9
    b = 0.6;    % recommended in Wedal etal "improved TV-L1" 2008
end;

if nargin < 10
    mask = [];
end

SY   = size(images, 1);
SX   = size(images, 2);
SV   = size(images, 3);
SU   = size(images, 4);

flow_proj = premultHM(flow_prev, intr);

[y,x,v,u]   = ndgrid(1:SY,1:SX,1:SV,1:SU);
x2          = x + flow_proj(:,:,:,:,1);
y2          = y + flow_proj(:,:,:,:,2);
u2          = u + flow_proj(:,:,:,:,3);
v2          = v + flow_proj(:,:,:,:,4);

[It,Ix,Iy,Iu,Iv] = partialDeriv(images,flow_proj,interpolation_method,...
    h1,h2,k1,k2,b,mask);

% Back-project to 3D spatial derivatives
if size(images, 6) == 1 % grayscale
    IXYZ = postmultHM(cat(5,Ix,Iy,Iu,Iv),intr,x2,y2,u2,v2);
    IX = IXYZ(:,:,:,:,1);
    IY = IXYZ(:,:,:,:,2);
    IZ = IXYZ(:,:,:,:,3);
else % color
    SC = size(images, 5);
    IX = zeros(SY,SX,SV,SU,SC);
    IY = zeros(SY,SX,SV,SU,SC);
    IZ = zeros(SY,SX,SV,SU,SC);
    for c = 1:SC
        IXYZ = postmultHM(cat(5,Ix(:,:,:,:,c),Iy(:,:,:,:,c),...
            Iu(:,:,:,:,c),Iv(:,:,:,:,c)),intr,x2,y2,u2,v2);
        IX(:,:,:,:,c) = IXYZ(:,:,:,:,1);
        IY(:,:,:,:,c) = IXYZ(:,:,:,:,2);
        IZ(:,:,:,:,c) = IXYZ(:,:,:,:,3);
    end
end

