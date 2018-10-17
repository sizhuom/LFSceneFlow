function this = parse_input_parameter(this, param)
%PARSE_INPUT_PARAMETER parses and set parameters in the struct PARAM
%
%   This is a member function of the class 'hs_optical_flow'. 

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
    
if mod(length(param), 2) ~=0
    error('Parse_input_parameter: Input parameters must be given in pairs (name and value)!');
end;

i = 1;

while (i <= length(param))
    
%     if ischar(param{i+1})
%         param{i+1} = str2num(param{i+1});
%     end;
    
    switch lower(param{i})
        
        case 'lambda'
            this.lambda         = param{i+1};
            
        case 'pyramid_levels'
            this.pyramid_levels  = param{i+1};
            
        case 'pyramid_spacing'
            this.pyramid_spacing = param{i+1};           
            
        case 'max_warping_iters'
            this.max_warping_iters       = param{i+1};
            
        case 'median_filter_size' 
            this.median_filter_size = [param{i+1} param{i+1}];
            
        case 'image_filter'
            if strcmpi(param{i+1}, 'gaussian')
                this.image_filter = fspecial('gaussian', [15 15], 0.4);            
            elseif strcmpi(param{i+1}, 'LoG')
                this.image_filter = - fspecial('gaussian', [15 15], 1);
                this.image_filter(8,8) = 1 + this.image_filter(8,8);                
            elseif strcmpi(param{i+1}, 'DoG')    
                this.image_filter =  fspecial('gaussian', [15 15], 0.4) - fspecial('gaussian', [15 15], 1);
            else
                this.image_filter = [];
            end;    
            
        case 'sigmap'
            this.sigmaP = param{i+1};
            
        case 'use_bs'
            this.use_bs = param{i+1};
            
        case 'lambdaz'
            this.lambdaZ         = param{i+1};
            
        case 'pcg_iters'
            this.pcg_iters = param{i+1};
            
        case 'deriv_filter_xy'
            this.deriv_filter_xy = param{i+1};
            
        case 'deriv_filter_uv'
            this.deriv_filter_uv = param{i+1};
            
        case 'interp_filter_xy'
            this.interp_filter_xy = param{i+1};
            
        case 'interp_filter_uv'
            this.interp_filter_uv = param{i+1};
            
        case 'solver'
            this.solver = param{i+1};
            
        case 'sor_max_iters'
            this.sor_max_iters = param{i+1};
            
        case 'sor_tol'
            this.sor_tol = param{i+1};
            
        case 'debug'
            this.debug = param{i+1};
            
        case 'debug_dir'
            this.debug_dir = param{i+1};
            
        case 'limit_update'
            this.limit_update = param{i+1};
            
        case 'xyweightsigma'
            this.xyweightSigma = param{i+1};
            
        case 'xyr'
            this.xyr = param{i+1};
    end;
    
    i = i+2;
end;
