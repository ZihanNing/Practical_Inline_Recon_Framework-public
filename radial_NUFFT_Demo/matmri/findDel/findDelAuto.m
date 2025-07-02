function [delSk, delSk_perIt] = findDelAuto(delSk0,data_in,tdwell,phs_spha,phs_conc,datatime,phs_sph_grid,b0map,Crcvr,maxNit_cgne,delJumpFact,numCoarseSearch)
    % Find optimal delay for traj data with respect to MRI raw data
    % Presumes a 2D slice 
    % Units for all times should be consistent (tdwell [time], datatime [time], b0map [rad/time])
    %
    % Inputs:
    %   delSk0:     starting guess for delay
    %   data_in:    MRI raw data
    %   tdwell:     dwell time for trajectory samples. Time for trajectory
    %               is presumed to start at 0.
    %   phs_spha:   trajectory using a higher order encoding model. See
    %               sampHighOrder class in matmri repository
    %               (https://doi.org/10.5281/zenodo.4495477)
    %   phs_conc:   phase from concomitant gradient fields. See sampHighOrder.
    %   datatime:   times for each data point in data_in
    %   phs_sph_grd: x, y, and z positions for each voxel. See sampHighOrder.
    %   b0map:      b0map. See sampHighOrder.
    %   Crcvr:      receiver sensitivity profiles. First three dims should
    %               have same size as b0map, with one extra 4th dim to
    %               count the number of receivers.
    %   maxNit_cgne: maximum number of iterations to use for each
    %                update of the image
    %   numCourseSearch: number of samples to perform a course search over
    %               e.g., for a value of 2, the algorithm finds which of
    %               the delays n*(datatime(2)-datatime(1)) gives the lowest
    %               residual when finding the image, where n spans from
    %               -numCourseSearch to +numCourseSearch. The iterative
    %               optimization then starts from that delay.
    %
    % (c) Corey Baron 2021

    % Options
    if nargin<12
        numCoarseSearch = 2;
    end
    if nargin<11
        delJumpFact = 3;
    end
    
    fprintf('findDelAuto:\n')
    
    % Basic prep computations
    numsamps = numCoarseSearch; % number of samples to try in initial course search
    opt_cgne.stopOnResInc = 1;
    R = rcvrOp(Crcvr,0); % Create receiver operator
    interpThresh = 0;    % Disable interpolation for iterations, since it could affect the results.
    dtsamp = datatime(2)-datatime(1);

    % Find numerical derivative of phs_spha and phs_conc
    kern = [1/8 1/4 0 -1/4 -1/8]';
    dspha_dt = conv2(phs_spha(:,:), kern, 'same')/tdwell; % units rad/{1,m,m^2,m^3}/time (spatial units depend on spatial order of term)
    dconc_dt = conv2(phs_conc(:,:), kern, 'same')/tdwell;

    % Peform course search
    resvec_min = inf;
    if numCoarseSearch > 0
        fprintf('    Course search (N=%d)', 2*numsamps+1);
        for courseTest = -numsamps:numsamps
            % Interpolate to match datatimes
            delSk = delSk0+courseTest*dtsamp;
            [phs_spha_a,phs_conc_a] = interpTrajTime(phs_spha,phs_conc,tdwell,delSk,datatime);

            % Update image
            S = sampHighOrder(b0map, datatime(:), phs_spha_a', phs_conc_a', phs_sph_grid, [], [], [], interpThresh);
            opFunc = @(x,transp) mrSampFunc(x,transp,S,R);
            % Use the same number of iterations for all recons
            opt_cgne_c.stopOnResInc = 0;
            [~, resvec] = cgne(opFunc,data_in(:),[],maxNit_cgne,opt_cgne_c); 
            if resvec(end)<resvec_min
                delSk_min = delSk;
                resvec_min = resvec(end);
            end
            fprintf('.');
        end
        fprintf('diff in delay from initial guess after course search: %.1f\n', delSk_min-delSk0)
        delSk = delSk_min;
    else
        delSk = delSk0;
    end
    
    % Get ready to start iterations
    ntry = 0;
    delChange = inf;
    delThresh = 0.005; % us
    fprintf('    diff in del from initial guess per iteration: ');
    delSk_perIt = zeros(1000,1);
    delSk_perIt(1) = delSk;
    while abs(delChange) > delThresh
        % Interpolate to match datatimes
        [phs_spha_a,phs_conc_a] = interpTrajTime(phs_spha,phs_conc,tdwell,delSk,datatime);
        [dspha_dt_a,dconc_dt_a] = interpTrajTime(dspha_dt,dconc_dt,tdwell,delSk,datatime);

        % Update image
        S = sampHighOrder(b0map, datatime(:), phs_spha_a', phs_conc_a', phs_sph_grid, [], [], [], interpThresh);
        opFunc = @(x,transp) mrSampFunc(x,transp,S,R);
        if ntry>0
            % Use the same number of iterations for all recons
            opt_cgne.stopOnResInc = 0;
            maxNit_cgne = length(resvec)-1;
        end
        [x, resvec] = cgne(opFunc,data_in,[],maxNit_cgne,opt_cgne);
        
        % Update delay
        bhat = data_in - opFunc(x,'notransp');
        clear opFunc
        S = S.setPhiDiv(dspha_dt_a',dconc_dt_a');
        opFunc = @(x,transp) mrSampFunc(x,transp,S,R);  
        bhathat = opFunc(x,'notransp');
        clear opFunc
        delChange = bhathat(:)\bhat(:);
        delChange = delJumpFact*real(delChange);
        delSk = delSk + delChange;
        fprintf('  %.2f', delSk - delSk0)
        if (ntry>0) && (sign(delChangePrev) ~= sign(delChange)) && (delJumpFact>1)
            delJumpFact = delJumpFact/2;
        end
        delChangePrev = delChange;
        ntry = ntry+1;
        delSk_perIt(ntry+1) = delSk;
        if ntry > 50
            warning('findDel did not converge after %d iterations\n', ntry)
            break;
        end
    end
    delSk_perIt = delSk_perIt(1:ntry+1);
    fprintf('\n    Done\n');
end
