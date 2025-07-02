function y = mrSampFunc(x,transp,samplingOp,R,opt)
% Perform MRI sampling.
%
% Usage: y = mrSampFunc(x,transp,samplingOp,R,opt)
%     y = mrSampFunc(x,'notransp',...) performs MRI sampling for
%     object-domain input x. y = mrSampFunc(x,'transp',...) performs
%     adjoint of MRI sampling for data samples x.
%
% Inputs:
%   x: when transp = 'transp', x is object-domain image data
%      when transp = 'notransp', x is sampled data (i.e., k-space samples)
%   transp: see description for x samplingOp: Applies MR sampling. Either a
%   mask for Cartesian sampling, or an object of the nufftOp class for
%   non-Cartesian sampling. Note that mask input for sampling opt must be
%   real-valued. For Cartesian, fft's are performed over ndims(samplingOp)
%   R: Applies receiver sensitivities. An object of the rcvrOp class. opt:
%   options struct with fields (okay to define only subset)
%       opt.imNFull (default = []): object domain data size. If not empty,
%         a vector input of x is assumed for transp = 'notransp', which is
%         then reshaped to opt.imNFull. For transp = 'transp', y is
%         vectorized as a final step of the function. 
%       opt.daNFull (default = []): sampling domain data size. If not 
%         empty, a vector input of x is assumed for transp = 'transp',
%         which is then reshaped to opt.daNFull. For transp = 'notransp', y
%         is vectorized as a final step of the function. 
%       opt.tikReg (default = []):  performs tikhonov regularization with 
%         weighting factor opt.tikReg, assumging mrSampFunc() is supplied
%         to a solver like lsqr().
%       opt.sumDims (default = []): performs sum along the listed 
%         dimensions during forward operation, after application of R.
%         Useful for SMS recons, or ESPIRiT-based R with multiple sets
%         of maps.
%  Forward operation converts image domain to data samples. Adjoint
%  converts data samples to image domain
%
% (c) Corey Baron 2020

% Parse inputs
if nargin<1
    tests;
    return;
end
if nargin<4
    R = [];
end
if nargin<5 || ~isfield(opt,'imNFull')
    opt.imNFull = [];
elseif ~isempty(opt.imNFull)
    opt.imNFull = [opt.imNFull,1];
end
if nargin<5 || ~isfield(opt,'daNFull')
    opt.daNFull = [];
elseif ~isempty(opt.daNFull)
    opt.daNFull = [opt.daNFull,1];
end
if nargin<5 || ~isfield(opt,'tikReg')
    opt.tikReg = [];
end
if nargin<5 || ~isfield(opt,'sumDims')
    opt.sumDims = [];
end
if nargin<5 || ~isfield(opt,'hardReal')
    opt.hardReal = 0;
end

% Check incompatibilities
if ~isempty(opt.tikReg) && (isempty(opt.imNFull) || isempty(opt.daNFull))
    error('Tikhonov regularization requires opt.imNFull and opt.daNFull')
end

switch transp
    case 'notransp'
        if opt.hardReal
            x = real(x);
        end
        if ~isempty(opt.tikReg)
            if numel(opt.tikReg) == 2
                x_tik = opt.tikReg(1)*real(x) + 1i*opt.tikReg(2)*imag(x);
            else
                x_tik = opt.tikReg*x;
            end
        end
        if ~isempty(opt.imNFull)
            x = reshape(x,opt.imNFull);
        end
        if ~isempty(R)
            x = R*x;
        end
        if ~isempty(opt.sumDims)
            % Note that the adjoint of this is not observed in transp case because of implicit repmat in matlab
            for n=1:length(opt.sumDims)
                x = sum(x,opt.sumDims(n));
            end
        end
        if isa(samplingOp,'nufftOp') || isa(samplingOp,'sampHighOrder')
            x = samplingOp*x;
        else
            nDim = ndims(samplingOp);
            if nDim==2 && size(samplingOp,2)==1
                nDim = 1;
            end
            x = samplingOp.*ifftnc(x,nDim);
        end
        if ~isempty(opt.daNFull)
            x = x(:);
        end
        if ~isempty(opt.tikReg)
            x = [x; x_tik];
        end
    case 'transp'
        if ~isempty(opt.tikReg)
            x_tik = x(prod(opt.daNFull)+1:end);
            x = x(1:prod(opt.daNFull));
        end
        if ~isempty(opt.daNFull)
            x = reshape(x,opt.daNFull);
        end
        if isa(samplingOp,'nufftOp') || isa(samplingOp,'sampHighOrder')
            x = samplingOp'*x;
        else
            nDim = ndims(samplingOp);
            if nDim==2 && size(samplingOp,2)==1
                nDim = 1;
            end
            x = samplingOp.*x;
            x = fftnc(x,nDim);
        end
        if ~isempty(R)
            x = R'*x;
        end
        if ~isempty(opt.imNFull)
            x = x(:);
        end
        if ~isempty(opt.tikReg)
            if numel(opt.tikReg) == 2
                x_tik = opt.tikReg(1)*real(x_tik) + 1i*opt.tikReg(2)*imag(x_tik);
            else
                x_tik = opt.tikReg*x_tik;
            end
            x = x + x_tik;
        end
        if opt.hardReal
            x = real(x);
        end
end   
y = x;
end

function tests
    N = 32;
    NR = 2;
    R = rcvrOp(rand([N N NR]),false);
    x = randn([N N]) + 1i*randn([N N]);
    for sumDims = true % [false,true]
        for vectorize_im = [false,true]
            for vectorize_da = [false,true]
                if vectorize_im && vectorize_da
                    tikTests = [false,true];
                else
                    tikTests = false;
                end
                for tikReg = tikTests
                    if sumDims
                        samplingOp = double(rand([N,1])>0.5);
                        y = randn([N 1 NR]) + 1i*randn([N 1 NR]);
                        opt.sumDims = 2;
                    else
                        samplingOp = double(rand([N,N])>0.5);
                        y = randn([N N NR]) + 1i*randn([N N NR]);
                        opt.sumDims = [];
                    end
                    if vectorize_im
                        opt.imNFull = [N N];
                        x_a = x(:);
                    else
                        opt.imNFull = [];
                        x_a = x;
                    end
                    if vectorize_da
                        opt.daNFull = size(y);
                        y_a = y(:);
                    else
                        opt.daNFull = [];
                        y_a = y;
                    end
                    if tikReg
                        opt.tikReg = 0.01;
                        y_a = [y_a; x_a];
                    else
                        opt.tikReg = [];
                    end
                    Fx = mrSampFunc(x_a,'notransp',samplingOp,R,opt);
                    Fy = mrSampFunc(y_a,'transp',samplingOp,R,opt);
                    d1 = dot(x_a(:),Fy(:));
                    d2 = dot(Fx(:),y_a(:));
                    assert(abs(d1-d2)/min(abs(d1),abs(d2)) < 1e-8, 'Adjoint test failed.')
                end
            end
        end
    end
    fprintf('mrSampFunc unit test success!\n')
end