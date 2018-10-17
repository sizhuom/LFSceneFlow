function this = parse_input_parameter(this, param)
%PARSE_INPUT_PARAMETER parses and set parameters in the struct PARAM
%
%   This is a member function of the class 'classic_nl_optical_flow'. 
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


if length(param) ==1
    param = param{1};
end;

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
            this.lambda_q       = param{i+1};           
            
        case 'pyramid_levels'
            this.pyramid_levels = param{i+1};
            this.auto_level     = false;
            
        case 'pyramid_spacing'
            this.pyramid_spacing        = param{i+1};
            
        case 'gnc_pyramid_levels'
            this.gnc_pyramid_levels     = param{i+1};
            
        case 'gnc_pyramid_spacing'    
            this.gnc_pyramid_spacing    = param{i+1};
            
        case 'sigma_d'
            for j = 1:length(this.spatial_filters);
                this.rho_spatial_u{j}.param   = param{i+1};
                this.rho_spatial_v{j}.param   = param{i+1};
            end;
            
        case 'sigma_s'
            this.rho_data.param         = param{i+1};
            
        case 'sigmap'
            this.sigmaP = param{i+1};
            
        case 'max_iters'
            this.max_iters = param{i+1};
            
        case 'gnc_iters'
            this.gnc_iters = param{i+1};
            
        case 'sigmats'
            this.sigmaTS = param{i+1};
            
        case 'lambdaz'
            this.lambdaZ         = param{i+1};
            this.lambdaZ_q       = param{i+1};    
            
        case 'median_filter_size'
            this.median_filter_size = [param{i+1} param{i+1}];
            
        case 'pcg_iters'
            this.pcg_iters = param{i+1};
            
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
            
        case 'segsigma'
            this.xyweightSigma = param{i+1};
            
        case 'segr'
            this.xyr = param{i+1};
            
        case 'deriv_filter_xy'
            this.deriv_filter_xy = param{i+1};
            
        case 'deriv_filter_uv'
            this.deriv_filter_uv = param{i+1};
            
        case 'interp_filter_xy'
            this.interp_filter_xy = param{i+1};
            
        case 'interp_filter_uv'
            this.interp_filter_uv = param{i+1};
            
    end;
    
    i = i+2;
end;