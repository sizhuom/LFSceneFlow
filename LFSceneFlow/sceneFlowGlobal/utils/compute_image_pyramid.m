function [P, intrP] = compute_image_pyramid(img, intr, smooth_sigma, nL, ratio)
%%  COMPUTE_IMAGE_PYRAMID computes nL level image pyramid of the input image IMG using filter F

%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2007-10-10 $
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
P   = cell(nL,1);
intrP = cell(nL,1);
tmp = img;
P{1}= tmp;
intrP{1} = intr;

for m = 2:nL
    % Gaussian filtering
    tmp = gaussFilt4D(tmp, smooth_sigma*ones(4));
    sz1 = [size(tmp,1) size(tmp,2) size(tmp,3) size(tmp,4)];
    sz2 = floor(sz1 * ratio);
    tmp2 = zeros([sz2 size(img,5)]);
    
    % Compute the intrinsic matrix for the new level
    a = 1 / ratio;
    b = (sz1+1)/2 - (sz2+1)/2*a;
    b = [b(2) b(1) b(4) b(3)];
    Hinc = [
        eye(4)*a b(:);
        zeros(1,4) 1;
        ];
    intrP{m} = intrP{m-1};
    intrP{m}.sz = sz2;
    intrP{m}.H = intrP{m-1}.H * Hinc;
    
    % Subsample the LF
    [y,x,v,u] = ndgrid(1:sz2(1),1:sz2(2),1:sz2(3),1:sz2(4));
    xyuv = [x(:)'; y(:)'; u(:)'; v(:)'; ones(1,numel(x))];
    xyuv = Hinc * xyuv;
    x = reshape(xyuv(1,:), sz2);
    y = reshape(xyuv(2,:), sz2);
    u = reshape(xyuv(3,:), sz2);
    v = reshape(xyuv(4,:), sz2);
    for c = 1:size(img,5)
        tmp2(:,:,:,:,c) = interpn(tmp(:,:,:,:,c),y,x,v,u,'linear');
    end
        
    P{m} = tmp2;
    tmp = tmp2;
end;







