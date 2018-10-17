function opq = compute_flow_base(this, opq)

%COMPUTE_FLOW_BASE   Base function for computing flow field
%   UV = COMPUTE_FLOW_BASE(THIS, INIT) computes the flow field UV with
%   algorithm THIS and the initialization INIT.
%
%   This is a member function of the class 'hs_optical_flow'.
%
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

% Iterate flow computation
for i = 1:this.max_warping_iters
fprintf('energy of solution \t%3.3e\n', -evaluate_log_posterior(this, opq));
    
% [It,IX,IY,IZ] = partial_deriv(this.images, this.intr, opq, this.interpolation_method, this.deriv_filter);
% fprintf('convex energy of solution (before) \t%3.3e\n',-evaluate_log_posterior_convex(this, opq, zeros(size(opq)), It, IX, IY, IZ));

    % Compute linear flow operator
    [A, b, parm, iterative] = ...
        flow_operator(this, opq);
    
    % Invoke the selected linear equation solver
    switch (lower(this.solver))
        case 'backslash'
            x = reshape(A \ b, size(opq));
        case 'sor'
%             [x, flag, res, n] = sor(A', b, 1.9, this.sor_max_iters, 1E-2);
%             x = reshape(x, size(opq));
%             fprintf('%d %d %d  ', flag, res, n);
            % Use matlab code from netlib
            [x, err, iter, flag]  = sor(A, zeros(size(b)), b, 1.9, this.sor_max_iters, this.sor_tol);
            relres = norm(A*x-b) / norm(b);
            x = reshape(x, size(opq));
            fprintf('flag: %d\terror: %g\tres: %g\titer: %d\n', flag, err, relres, iter);
        case 'bicgstab'
            x = reshape(bicgstab(A, b, 1E-3, 200, [], [], opq(:), parm), size(opq));
        case 'pcg'
            alpha = max(sum(abs(A),2)./diag(A))-2;
            L = ichol(A, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
            [x, flag, relres, iter] = pcg(A,b, [], this.pcg_iters, L, L');
%             [x, flag, relres, iter] = pcg(A,b, 1e-6, 2000);
            x        = reshape(x, size(opq));
            fprintf('flag: %d\tres: %g\titer: %d\n', flag, relres, iter);
        otherwise
            error('Invalid solver!')
    end
    
% fprintf('convex energy of solution (after) \t%3.3e\n',-evaluate_log_posterior_convex(this, opq, x, It, IX, IY, IZ));

    % Print status information
    if this.display
        disp(['--Iteration: ', num2str(i), '    (', ...
            num2str(norm(x(:))), ')'])
    end;
    
    % Terminate iteration early if flow doesn't change substantially
    if (length(this.lambda) == 1 && norm(x(:)) < 1E-3)
        break
    end
    
    % If limiting the incremental flow to [-1, 1] is requested, do so
    limit = this.intr.H(1,1)*this.intr.S;
    if (this.limit_update)
        x(x > limit)  = limit;
        x(x < -limit) = -limit;
    end
    
    opq = opq + x;
    
    % Perform median filtering to remove outliers
%     timerVal = tic;
    if ~isempty(this.median_filter_size)
        for m = 1:this.mf_iter; % extensive MF filtering
            LFMedfilt2(opq, this.median_filter_size);
        end;
    end;
%     timeElapsed = toc(timerVal);
%     fprintf('median filter: %f seconds.\n', timeElapsed);
    
    % debug: save result after each iteration
    if this.debug
        save(findNextTmp(this.debug_dir,'result-*.mat'), 'opq');
        fprintf('iteration data saved.\n');
    end
    
    % Terminate early if the flow_operator doesn't require multiple
    % interations
    if (~iterative)
        break;
    end
    
end
