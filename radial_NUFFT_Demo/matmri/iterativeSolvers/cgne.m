function [x, resSqAll, mseAll, xnormAll, xdiffAll, stopThresh] = cgne(Ain,bin,x0,NitMax,opt)
    % Use conjugate gradient method to solve Ax = b in a least squares sense. 
    %    i.e., "conjugate gradient on the normal equations", cjne
    %
    % x = cgne(A,b,x0,maxIterations,options)
    %
    %   A must have a transpose that can be evaluated using A'*b, or as a function with A(b,'transp')
    %   A operates on x via A*x or A(x,'notransp')
    %
    % (c) Corey Baron 2020
    %
    
    % Set options
    if nargin<4 || isempty(NitMax)
        % Maximum number of iterations allowed
        NitMax = 1000;
    end
    if nargin<5 || ~isfield(opt,'expResN')
        % Set to non-zero to explicitely compute residual every expResN iterations. 
        %   also restarts conj grad iterations by setting beta to 0
        % Can reduce accumulation of roundoff errors
        opt.expResN = 100; 
    end
    if nargin<5 || ~isfield(opt,'nseThreshFact')
        % Set factor to increase estimated noise by for stopping iterations
        % based on the discrepancy principle. 
        % For stopping criterion, see a great review here: 
        % https://doi-org.proxy1.lib.uwo.ca/10.3846/1392-6292.2007.12.61-70
        opt.nseThreshFact = 1.0; 
    end
    if nargin<5 || ~isfield(opt,'nseExtraIt')
        % Stopping on noise threshold seems to often stop a little too
        % early, so we add some extra iterations.
        opt.nseExtraIt = 4; 
    end
    if nargin<5 || ~isfield(opt,'stopOnResInc')
        % Stop on nth iteration where ||res|| increases, where n is the
        % value that opt.stopOnResInc is set to. This works well for single
        % shot spiral MRI. 
        % Similar to methods used for CT, https://doi.org/10.1016/S0096-3003(98)10007-3
        opt.stopOnResInc = 0; 
    end
    if nargin<5 || ~isfield(opt,'resIncThresh')
        % fractional threshold for stopOnResInc. Should be close to 1
        opt.resIncThresh = 1.0; 
    end
    if nargin<5 || ~isfield(opt,'stopOnXdifInc')
        % Stop on nth iteration where ||x_{k+1}-x_k|| increases, where n is the
        % value that opt.stopOnXdifInc is set to. 
        % Similar to methods used for CT, https://doi.org/10.1016/S0096-3003(98)10007-3
        opt.stopOnXdifInc = 0; 
    end
    if nargin<5 || ~isfield(opt,'noiseVar')
        % Variance of the noise expected in the input bin
        opt.noiseVar = []; 
    end
    if nargin<5 || ~isfield(opt,'plotting')
        % Show plots of progress
        opt.plotting = 0; 
    end
    if nargin<5 || ~isfield(opt,'plottingReshape')
        % Reshape x for plotting images. Only 2D matrix allowed.
        opt.plottingReshape = []; 
    end
    if nargin<5 || ~isfield(opt,'gtruth')
        % Useful for simulations. Allows computation of mse per iteration
        opt.gtruth = []; 
    end

    if (nargout > 3) 
        findxnorm = 1;
    else
        findxnorm = 0;
    end
    
    % Account for different ways of supplying A
    if isa(Ain,'function_handle')
        A = Ain;
    else
        A = @(x,transp) Asub(x,transp,Ain);
    end
    
    % Set default starting guess
    if nargin<3 || isempty(x0)
        x0 = A(bin,'transp');
    end

    % Prep for stopping based on measured noise
    stopThresh = [];
    if (~isempty(opt.noiseVar) && opt.noiseVar>0)
        % Create a synthetic noise vector
        if ~isreal(bin)
            % The scaling by 0.5 presumes the noise was estimated from
            % complex data
            nseVecIn = sqrt(opt.noiseVar/2)*(randn(size(bin)) + 1i*randn(size(bin)));
        else
            nseVecIn = sqrt(opt.noiseVar)*randn(size(bin));
        end
        % Pass the noise through A' to find the expected noise in b
        nse = A(nseVecIn,'transp');
        % Define the stopping criteria
        stopThresh = opt.nseThreshFact * nse(:)' * nse(:);
    end
    
    % Initialize first iteration
    b = A(bin,'transp');
    xold = x0;
    res = b - A(A(xold, 'notransp'), 'transp');
    p = res;
    res_sq_old = res(:)' * res(:);
    nit = 0;
    if opt.plotting
        nf = figure;
        Ncol = ceil(sqrt(NitMax));
        Nrow = ceil(NitMax/Ncol);
        if ~isempty(opt.gtruth)
            nf2 = figure;
        end
    end
    
    % Run iterations
    resSqAll = zeros(NitMax+1,1);
    resSqAll(1) = res_sq_old;
    xnormAll = zeros(NitMax+1,1);
    if findxnorm
        xnormAll(1) = x0(:)'*x0(:);
    end
    xdiffAll = zeros(NitMax+1,1);
    xdiffAll(1) = xnormAll(1);
    mseAll = zeros(NitMax+1,1);
    if ~isempty(opt.gtruth)
        tmp = opt.gtruth(:)-x0(:);
        mseAll(1) = tmp'*tmp;
    end
    finished = 0;
    resIncsTotal = 0;
    xdifIncsTotal = 0;
    lastxdifInc = 0; % consider sequential climbing xdiff as a single increase
    while ~finished
        Ap = A(A(p, 'notransp'), 'transp');
        alph = res_sq_old / (p(:)' * Ap(:));
        x = xold + alph*p;
        if mod(nit,opt.expResN) == 0
            % Avoid accumulation of rounding errors when doing many iterations
            res = b - A(A(x, 'notransp'), 'transp');
            betaFact = 0;
        else
            res = res - alph*Ap;
            betaFact = 1;
        end
        res_sq_new = res(:)' * res(:);
        resSqAll(nit+2) = res_sq_new; 
        if findxnorm
            xnormAll(nit+2) = x(:)'*x(:);
        end
        xdiff = x-xold;
        xdiffAll(nit+2) = xdiff(:)'*xdiff(:);
        if ~isempty(opt.gtruth)
            tmp = opt.gtruth(:)-x(:);
            mseAll(nit+2) = tmp'*tmp;
        end
        % Plotting
        if opt.plotting
            figure(nf);
            subplot(Nrow,Ncol,nit+1);
            tmp = abs(xdiff);
            if ~isempty(opt.plottingReshape)
                tmp = reshape(tmp, opt.plottingReshape);
                imagesc(tmp);
            else
                plot(tmp);
            end
            hold('all')
            if ~isempty(opt.gtruth) && ~isempty(opt.plottingReshape)
                figure(nf2)
                subplot(Nrow,Ncol,nit+1);
                imagesc(reshape(abs(tmp), opt.plottingReshape));
                colorbar;
                hold('all')
            end
        end
        % Check if converged 
        if (res_sq_new/resSqAll(1) < 1e-10)
            % Here the residual is tiny
            finished = 1;
        end
        if (~isempty(stopThresh) && (res_sq_new < stopThresh)) 
            % Here the residual dropped below the residual expected from
            % noise.
            NitMax = nit + opt.nseExtraIt;
            stopThresh = 0;
        end
        if opt.stopOnResInc>0 && (resSqAll(nit+2)/resSqAll(nit+1) > opt.resIncThresh)
            resIncsTotal = resIncsTotal + 1;
            if resIncsTotal >= opt.stopOnResInc
                finished = 1;
            end
        end
        if opt.stopOnXdifInc>0 && (xdiffAll(nit+2) > xdiffAll(nit+1)) 
            if ~lastxdifInc
                xdifIncsTotal = xdifIncsTotal + 1;
                lastxdifInc = 1;
                if xdifIncsTotal >= opt.stopOnXdifInc
                    finished = 1;
                end
            end
        else
            lastxdifInc = 0;
        end
        % Finish up current iteration
        p = res + betaFact * (res_sq_new / res_sq_old) * p;
        res_sq_old = res_sq_new;
        xold = x;
        nit = nit+1;
        if nit >= NitMax
            finished = 1;
        end
    end
    resSqAll = resSqAll(1:(nit+1));
    xnormAll = xnormAll(1:(nit+1));
    mseAll = mseAll(1:(nit+1));
    
end

function y = Asub(x,transp,Ain)
    if strcmp(transp,'notransp')
        y = Ain*x;
    else
        y = Ain'*x;
    end
end