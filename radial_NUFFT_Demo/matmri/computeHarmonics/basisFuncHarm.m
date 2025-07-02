function bfuncOut = basisFuncHarm(x,y,z,n,subInds)
    % Outputs basis function for index n of spherical harmonics
    % n = index of phs_spha in sampHighOrder.m
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
            bfunc = ones(size(x), 'like', x);
        case 2
            bfunc = x;
        case 3                
            bfunc = y;
        case 4                
            bfunc = z;
        case 5                
            bfunc = (x.*y);
        case 6               
            bfunc = (z.*y);
        case 7                
            bfunc = ((3.*z.^2 - (x.^2 + y.^2 + z.^2)));
        case 8                
            bfunc = (x.*z);
        case 9                
            bfunc = (x.^2 - y.^2);
        case 10               
            bfunc = (3.*y.*x.^2 - y.^3);
        case 11                
            bfunc = (x.*z.*y);
        case 12               
            bfunc = ((5.*z.^2 - (x.^2 + y.^2 + z.^2)).*y);
        case 13               
            bfunc = (5*z.^3 - 3.*z.*(x.^2 + y.^2 + z.^2));
        case 14                
            bfunc = ((5.*z.^2 - (x.^2 + y.^2 + z.^2)) .* x);
        case 15                
            bfunc = (x.^2.*z - y.^2.*z);
        case 16                
            bfunc = (x.^3 - 3.*x.*y.^2);
        otherwise
            error('Basis function undefined.')
        end
        bfuncOut(:,n_ind) = bfunc(:);
    end
    bfuncOut = squeeze(reshape(bfuncOut,[szx,size(n)]));
end