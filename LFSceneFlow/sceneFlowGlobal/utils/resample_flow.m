function out = resample_flow(opq, sz, method)
% function out = resample_flow(uv, factor, method)
%RESAMPLE_FLOW   Resample flow field
%   OUT = RESAMPLE_FLOW(IN, FACTOR[, METHOD]) resamples (resizes) the flow
%   field IN using a factor of FACTOR.  The optional argument METHOD
%   specifies the interpolation method ('bilinear' (default) or
%   'bicubic').
%
%   This is a private member function of the class 'clg_2d_optical_flow'.
%
%   Author:  Stefan Roth, Department of Computer Science, TU Darmstadt
%   Contact: sroth@cs.tu-darmstadt.de
%   $Date$
%   $Revision$

% Copyright 2004-2007, Brown University, Providence, RI. USA
% Copyright 2007-2010 TU Darmstadt, Darmstadt, Germany.
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

% Make bilinear the default method
if (nargin < 3)
    method = 'bilinear';
    %     method = 'bicubic';
end

sz1 = size(opq); sz1 = sz1(1:4);
ratio = sz ./ sz1;
a = 1 ./ ratio;
b = (sz1+1)/2 - (sz+1)/2.*a;
a = [a(2) a(1) a(4) a(3)];
b = [b(2) b(1) b(4) b(3)];
Hinc = [
    diag(a) b(:);
    zeros(1,4) 1;
    ];

[y,x,v,u] = ndgrid(1:sz(1),1:sz(2),1:sz(3),1:sz(4));
xyuv = [x(:)'; y(:)'; u(:)'; v(:)'; ones(1,numel(x))];
xyuv = Hinc * xyuv;
x = reshape(xyuv(1,:), sz);
y = reshape(xyuv(2,:), sz);
u = reshape(xyuv(3,:), sz);
v = reshape(xyuv(4,:), sz);
x(x<1) = 1; x(x>sz1(2)) = sz1(2);
y(y<1) = 1; y(y>sz1(1)) = sz1(1);
u(u<1) = 1; u(u>sz1(4)) = sz1(4);
v(v<1) = 1; v(v>sz1(3)) = sz1(3);
out = zeros([sz size(opq,5)]);
for c = 1:size(opq,5)
    out(:,:,:,:,c) = interpn(opq(:,:,:,:,c),y,x,v,u,method);
end




