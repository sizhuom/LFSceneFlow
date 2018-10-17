function opq = compute_flow_base(this, opq)
%%
%COMPUTE_FLOW_BASE   Base function for computing flow field
%   UV = COMPUTE_FLOW_BASE(THIS, INIT) computes the flow field UV with
%   algorithm THIS and the initialization INIT.
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

% TODO: adapt this part to enable GNC
% Construct quadratic formulation
qua_this          = this;
qua_this.lambda   = this.lambda_q;
qua_this.lambdaZ  = this.lambdaZ_q;

% Spatial
if isa(this.rho_spatial_u{1}, 'robust_function')
    for i = 1:length(this.rho_spatial_u)
        a = this.rho_spatial_u{i}.param;
        qua_this.rho_spatial_u{i}   = robust_function('quadratic', a(1));
        a = this.rho_spatial_v{i}.param;
        qua_this.rho_spatial_v{i}   = robust_function('quadratic', a(1));
        a = this.rho_spatial_w{i}.param;
        qua_this.rho_spatial_w{i}   = robust_function('quadratic', a(1));
    end;
elseif isa(this.rho_spatial_u{1}, 'gsm_density')
    for i = 1:length(this.rho_spatial_u)
        qua_this.rho_spatial_u{i}   = robust_function('quadratic', sqrt(1/this.rho_spatial_u{i}.precision));
        qua_this.rho_spatial_v{i}   = robust_function('quadratic', sqrt(1/this.rho_spatial_v{i}.precision));
        qua_this.rho_spatial_w{i}   = robust_function('quadratic', sqrt(1/this.rho_spatial_w{i}.precision));
    end;
else
    error('evaluate_log_posterior: unknown rho function!');
end;

% Data
if isa(qua_this.rho_data, 'robust_function')
    a = this.rho_data.param;
    qua_this.rho_data        = robust_function('quadratic', a(1));
elseif isa(qua_this.rho_data, 'gsm_density')
    qua_this.rho_data        = robust_function('quadratic', sqrt(1/this.rho_data.precision));
else
    error('evaluate_log_posterior: unknown rho function!');
end;

% Iterate flow computation
for i = 1:this.max_iters
    
    dopq = zeros(size(opq));
    
    % Compute spatial and temporal partial derivatives
    [It, IX, IY, IZ] = partial_deriv(this.images, this.intr, opq, this.interpolation_method,...
        this.deriv_filter_xy, this.deriv_filter_uv, ...
        this.interp_filter_xy, this.interp_filter_uv,...
        0.6, this.mask);
    
    for j = 1:this.max_linear
%         evaluate_energy_map(this,opq);
        fprintf('energy of solution \t%3.3e\n', -evaluate_log_posterior(this, opq));
        
        % Every linearization step updates the nonlinearity using the
        % previous flow increments
        
        % Compute linear flow operator
        if this.alpha == 1
            [A, b, parm, iterative] = ...
                flow_operator(qua_this, opq, dopq, It, IX, IY, IZ);
            
        elseif this.alpha > 0
            [A, b] = ...
                flow_operator(qua_this, opq, dopq, It, IX, IY, IZ);
            [A1, b1, parm, iterative] = ...
                flow_operator(this, opq, dopq, It, IX, IY, IZ);
            A = this.alpha * A + (1-this.alpha) * A1;
            b = this.alpha * b + (1-this.alpha) * b1;
            
        elseif this.alpha == 0
            [A, b, parm, iterative] = ...
                flow_operator(this, opq, dopq, It, IX, IY, IZ);
            
        else
            error('flow_operator@classic_nl_optical_flow: wrong gnc parameter alpha %3.2e', this.alpha);
        end;
        
        % Invoke the selected linear equation solver
        switch (lower(this.solver))
            case 'backslash'
                x = reshape(A \ b, size(opq));
            case 'sor'
                % Use complied mex file (may need to compile utils/mex/sor.pp)
%                 [x, flag, res, n] = sor(A', b, 1.9, this.sor_max_iters, 1E-2, opq(:));
                % Use matlab code from netlib
                [x, err, iter, flag]  = sor(A, zeros(size(b)), b, 1.9, this.sor_max_iters, this.sor_tol);
                relres = norm(A*x-b) / norm(b);
                x = reshape(x, size(opq));
                fprintf('flag: %d\terror: %g\tres: %g\titer: %d\n', flag, err, relres, iter);
                
            case 'bicgstab'
                [x,flag] = reshape(bicgstab(A, b, 1E-3, 200, [], [], opq(:)), size(opq)); %, parm
            case 'pcg'
                alpha = max(sum(abs(A),2)./diag(A))-2;
                L = ichol(A, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
                [x, flag, relres, iter] = pcg(A,b, 1e-6, this.pcg_iters, L, L');
%                 [x, flag, relres, iter] = pcg(A,b, 1e-6, 2000);
                x        = reshape(x, size(opq));
                fprintf('flag: %d\tres: %g\titer: %d\n', flag, relres, iter);
            otherwise
                error('Invalid solver!')
        end
        
        % If limiting the incremental flow to [-1, 1] is requested, do so
        if (this.limit_update)
            limit = this.intr.H(1,1)*this.intr.S;
            limitZ = this.intr.H(1,1)*this.intr.S*this.intr.scaleZ;
            
            xyFlow = x(:,:,:,:,1:2);
            xyFlow(xyFlow > limit)  = limit;
            xyFlow(xyFlow < -limit) = -limit;
            
            zFlow = x(:,:,:,:,3);
            zFlow(zFlow > limitZ) = limitZ;
            zFlow(zFlow < -limitZ) = -limitZ;
            
            x = cat(5, xyFlow, zFlow);
        end
        
        % Print status information (change betwen flow increment = flow
        %   increment at first linearization step)
        if this.display
            disp(['--Iteration: ', num2str(i), '   ', num2str(j), '   (norm of flow increment   ', ...
                num2str(norm(x(:)-dopq(:))), ')'])
        end;
        
        % Terminate iteration early if flow doesn't change substantially
        %         if norm(x(:)-duv(:)) < 1E-3
        %             break;
        %         end
        
        dopq = x;
        
        opq0 = opq;
        opq  = opq+dopq;
        
        timerVal = tic;
        if ~isempty(this.median_filter_size)
            %Compute weighted median solved by Li & Osher formula
            occ = detect_occlusion(opq, this.images, this.intr);
            for mi = 1:size(opq,1)
                for mj = 1:size(opq,2)
                    opq(mi,mj,:,:,:) = denoise_color_weighted_medfilt3_2d(...
                        squeeze(opq(mi,mj,:,:,:)), ...
                        squeeze(this.color_images(mi,mj,:,:,:)), ...
                        squeeze(occ(mi,mj,:,:,:)),...
                        this.area_hsz, this.median_filter_size, this.sigma_i, this.fullVersion);
                end
            end
        end;
%         % Perform median filtering to remove outliers
%         if ~isempty(this.median_filter_size)
%             LFMedfilt2(opq, this.median_filter_size);
%         end;
        timeElapsed = toc(timerVal);
        fprintf('median filter: %f seconds.\n', timeElapsed);
        
        dopq = opq - opq0;
        opq  = opq0;
        
    end;
    
    % Update flow fileds
    opq = opq + dopq;
    
    % debug: save result after each iteration
    if this.debug
        save(findNextTmp(this.debug_dir,'result-*.mat'), 'opq');
        fprintf('iteration data saved.\n');
    end
    
end
