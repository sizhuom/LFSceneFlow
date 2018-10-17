function occ = detect_occlusion2(op, images, intr, sigma_d, sigma_i)

% Detect occlusion regions using flow divergence and brightness constancy
% error

% the output taks continous value 
% close to 0: occluded
% close to 1: nonoccluded

% according to Appendix A.2 Sand & Teller "Particle Video" IJCV 2008

%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2009$
%   $Revision $
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


if nargin == 3
    sigma_d = 0.3; % divergence
    sigma_i = 20;  % intensity
end;

op_pad = cat(5,op,zeros(size(op,1),size(op,2),size(op,3),size(op,4)));
flow_proj = premultHM(op_pad, intr);
[~,Xx,~,~,~] = partialDeriv(flow_proj(:,:,:,:,1));
[~,~,Yy,~,~] = partialDeriv(flow_proj(:,:,:,:,2));
[~,~,~,Uu,~] = partialDeriv(flow_proj(:,:,:,:,3));
[~,~,~,~,Vv] = partialDeriv(flow_proj(:,:,:,:,4));
div = Xx + Yy + Uu + Vv;

div(div>0) =0;

It = partialDeriv(images, flow_proj);
occ = exp(-div.^2/2/sigma_d^2).*exp(-It.^2/2/ sigma_i^2);

