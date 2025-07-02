function bfunc = basisFuncConcGrad(x,y,z,n,subInds)
    % Outputs basis function for index n of concomitant gradient terms
    % n = index of phs_conc in sampHighOrder.m
    if (nargin > 4) && ~isempty(subInds)
        x = x(subInds(1):subInds(2));
        y = y(subInds(1):subInds(2));
        z = z(subInds(1):subInds(2));
    end
    szx = size(x);
    bfuncOut = zeros([numel(x),numel(n)],'like',x);
    for n_ind = 1:numel(n)
        n_a = n(n_ind);
        switch n_a
        case 1
            bfunc = z.*z;
        case 2
            bfunc = x.*x + y.*y;
        case 3                
            bfunc = x.*z;
        case 4                
            bfunc = y.*z;
        otherwise
            error('Basis function undefined.')
        end
        bfuncOut(:,n_ind) = bfunc(:);
    end
    bfuncOut = squeeze(reshape(bfuncOut,[szx,size(n)]));
end