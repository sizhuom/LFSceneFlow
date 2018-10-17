function this = hs_xy_scene_flow(varargin)
%
%HS_OPTICAL_FLOW   Optical flow computation with Horn & Schunck method
%       B.K.P. Horn and B.G. Schunck. Determining optical flow. 
%       Artificial Intelligence, 16:185-203, Aug. 1981.
%       http://people.csail.mit.edu/bkph/papers/Optical_Flow_OPT.pdf
%       
%   HS_OPTICAL_FLOW([IMGS]) constructs a HS optical flow object
%   with the optional image sequence IMGS ([n x m x 2] array). 
%   HS_OPTICAL_FLOW(O) constructs HS optical flow object by copying O.
%  
%   This is a member function of the class 'hs_optical_flow'. 

%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2007-10-30 $
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

narginchk(0,1);
  
  switch (length(varargin))
    case 0
        
      this.images               = [];             
      this.intr                 = struct();   % intrinsics for light fields
      this.lambda               = 80;   %5e-2;
      this.lambdaZ         = 80;    % weight for z motion
      this.lambda_q             = 80;   %not used, for consistency with other program
      this.gnc_iters            = 1;    %not used, for consistency with other program
      this.pyramid_levels       = 2; %
      
      this.pyramid_spacing      = 2; %         % default to the prolongation/restriction filter
      
      
      this.max_warping_iters    = 10;           % # of warping/linearization per pyramid level
      this.median_filter_size   = [];
      this.texture              = false;        % apply to the image pyramid;     
      
      this.solver               = 'sor';  % pcg, sor      
      this.sor_max_iters   = 2000;       % 100 seems sufficient
      this.sor_tol = 1e-3; 
      this.pcg_iters       = 2000;
      
      this.interpolation_method = 'bi-linear';      % 'bi-cubic', 'cubic', 'bi-linear'
      this.deriv_filter_xy      = [-1 0 1]/2;
      this.deriv_filter_uv      = [1 -8 0 8 -1]/12;
      this.interp_filter_xy     = [1 2 1]/4;
      this.interp_filter_uv     = [1 2 1]/4;
      this.display              = true;
      this.limit_update         = false;
      this.sigmaD2              = 100;                 % data term
      this.sigmaS2              = 100;                 % spatial term
      
      this.sigmaP               = [0 0 0 0];                 % sigma of Gaussian for performing "texture" decomposition
      this.weightP              = 1;                 % images -weightPxGaussian*images

      this.mf_iter              = 1;
      
      method = 'quadratic'; 
      this.spatial_filters = {[1 -1], [1; -1],...
          reshape([1 -1], [1 1 1 2]), reshape([1 -1], [1 1 2 1])};
      for i = 1:length(this.spatial_filters);
          this.rho_spatial_u{i}   = robust_function(method, 1); % 0.1
          this.rho_spatial_v{i}   = robust_function(method, 1);
      end;
      this.rho_data        = robust_function(method, 1); % 6.3
      
      this.color_images     = [];
      
      this.use_bs = false;  % use b-spline for interpolation
      this.bs_degree = 3;
      this.bs_lambda = 0;
      this.bs_images = [];
      
      this.mask = []; % only compute where mask=true
      
      this.xyflow = []; % used for segmentation
      this.xyweightSigma = 0.3;
      this.xyweightMap = [];
      
      this.debug = false;
      this.debug_dir = '~/Documents';
      
      this = class(this, 'hs_xy_scene_flow');         

      
    case 1
      if isa(varargin{1}, 'hs_xy_scene_flow')
        this = other;
        
      else    
          this = hs_xy_scene_flow;
      end
      
    otherwise


      error('Incompatible arguments!');
      
  end
