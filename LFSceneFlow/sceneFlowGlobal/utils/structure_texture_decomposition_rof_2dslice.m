function [texture, structure] = structure_texture_decomposition_rof_2dslice(lf, theta, nIters, alp)
%
% Decompose the input IMAGE into structure and texture parts using the
% Rudin-Osher-Fatemi method. The final output is a linear combination 
% of the decomposed texture and the structure parts.
%
% According to Wedel etal "An Improved Algorithm for TV-L1 optical flow"
%       equations (8)-(10)
%
% Test code
% [im1, im2, tu, tv] = read_image_flow_tune_para(4,0);
% [t1 s1] = structure_texture_decomposition_rof(im1); 
% [t2 s2] = structure_texture_decomposition_rof(im2); 
% indx = ~isnan(tu) | ~isnan(tv);
% uv = cat(3, tu, tv);
% uv(isnan(uv)) = 0;
% figure; imshow(-abs(partial_deriv(cat(3,im1, im2), uv)).*indx, []); title('original');
% figure; imshow(-abs(partial_deriv(cat(3,s1, s2), uv)).*indx, []);  title('structure');
% figure; imshow(-abs(partial_deriv(cat(3,t1, t2), uv)).*indx, []);  title('texture');
% tmp = partial_deriv(cat(3,t1, t2, uv);
% figure; imshow(t1, []); figure; imshow(s1, []);

%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2009 $
%
% Copyright 2009-2010, Brown University, Providence, RI. USA
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

if nargin == 1
    theta   = 1/8; 
    nIters  = 100;
    alp     = 0.95;  % alp = 0.75 results in 4:1
end;

SY = size(lf, 1);
SX = size(lf, 2);
texture = zeros(size(lf));
if nargout > 1
    structure = zeros(size(lf));
end

for y = 1:SY
    for x = 1:SX
        if nargout < 2
            [texture(y,x,:,:,:), structure(y,x,:,:,:)] = structure_texture_decomposition_rof(squeeze(lf(y,x,:,:,:)),theta,nIters,alp);
        else
            texture(y,x,:,:,:) = structure_texture_decomposition_rof(squeeze(lf(y,x,:,:,:)),theta,nIters,alp);
        end
    end
end