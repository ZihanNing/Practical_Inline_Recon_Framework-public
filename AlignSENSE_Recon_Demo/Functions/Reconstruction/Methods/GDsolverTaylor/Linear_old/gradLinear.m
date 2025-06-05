
function [gradFi, gradRe] = gradLinear(y,x,E,EH,W,sep)

%GRADLINEAR   Computes the gradient of a least squares 
%problem with different regularisations w.r.t. the linaer model of pose-dependent B0 variations. 
%   [EN,W]=GRADLINEAR(Y,X,E,EH,{W},{SEP})
%   * Y is the measured data
%   * X is the reconstructed data
%   * E is the encoding structure
%   * EH is the decoding structure (used for weighted least squares)
%   * {W} is the weight to use for certain regularisers
%   * {SEP} whether to separate the data fidelity and the regularisation gradient terms
%
if nargin < 5 || isempty(W); W = ones(size(x), 'like',real(x)); applyWeight = 0; else applyWeight = 1; end %applyWeight used for efficiency lambda_susc
if nargin < 6 || isempty(sep); sep=(nargout>1);end
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'BlSzStates');BlSzStates = E.Dl.Optim.BlSzStates; else BlSzStates=E.bS(1) ;end% Block size for coefficient gradient calculation

%%%INITIALISE
NStates  = E.NMs; 
N = size(E.Dl.D);
    
%%%DATA FIDELITY
gradFi = zeros(N, 'like',E.Dl.D);
res=y-encode(x,E);
for a=1:BlSzStates:NStates %Chunk states otwerwise is has size [Nx Nstates] and this runs out of memmory depending on NX and NStates
    E.oS(1) = a; E.bS(1) = BlSzStates; E.dS(1) = min(a+BlSzStates-1,NStates);
    vA=a:min(a+E.bS(1)-1,NStates);

    notsum = 1; % if 0, decode step will sum over motion states
    [~,EHres] =  decode(res,EH,E,notsum); % motion states in the 5th dimension
        %rotPar = E.Dl.fun(dynInd(E.Tr, {vA}, 5)); %Rotation parameters used in the Taylor model
    rotPar = convertRotation(dynInd(E.Tr, {vA,E.Dl.d}, 5:6) , 'rad','deg');
    EHres = bsxfun(@times, EHres , rotPar ); %might have last segment that is issue here
    xToMultiply = exp(1i*angle(x));%Play with this term to investigate effect of signal voids on D update
    if a==1 && ~isequal(x,xToMultiply);fprintf('gradientLinear:: taking magnitude out of x.\n');end
    temp = bsxfun(@times , EHres, conj(xToMultiply) );
    gradFi = gradFi + multDimSum(real( +1i * 2*pi*E.Dl.TE * temp) , 5 ) / numel(y) ; % numel(y) since in objective function relative error 
end
    
%%%REGULARISATION - only L2 gradients as L1 incorporated in (F)ISTA
gradRe = [];
dim=4; %dimension where to store the different regularisation terms

if isfield(E.Dl,'Reg') && isfield(E.Dl.Reg,'lambda_smooth') && ~isequal(E.Dl.Reg.lambda_smooth,0)%Smoothness
    gradRe = cat(dim,gradRe, E.Dl.Reg.lambda_smooth * ( FiniteDiff( W .*  FiniteDiff(E.Dl.D,1),1 ,0) + ...
                                                        FiniteDiff( W .*  FiniteDiff(E.Dl.D,2),2 ,0) + ...
                                                        FiniteDiff( W .*  FiniteDiff(E.Dl.D,3),3 ,0) ) / numel(E.Dl.D));
end

if isfield(E.Dl,'Reg') && isfield(E.Dl.Reg,'lambda_ridge') && ~isequal(E.Dl.Reg.lambda_ridge,0)%Ridge 
    gradRe = cat(dim,gradRe, E.Dl.Reg.lambda_ridge * E.Dl.D /numel(E.Dl.D));
end

if isfield(E.Dl,'Reg') && isfield(E.Dl.Reg,'lambda_susc') && ~isequal(E.Dl.Reg.lambda_susc,0)%Implicit susceptibility model
    K = linearTerms(N,pi,6);% create linear terms

    part2 = dynInd(K,2,6) .* fftn(fftshift(dynInd(E.Dl.D,1,6))) + dynInd(K,1,6) .* fftn(fftshift(dynInd(E.Dl.D,2,6)));
    if applyWeight; part2 = fftn(fftshift(abs(W).*ifftshift(ifftn(part2))));end %avoid two additional fft's if isequal(W,ones(size(W))==1
    
    gradTemp = zeros(N, 'like', real(E.Dl.D)); %make sure gradient real
    for i = 1:N(6)
        temp = ifftshift( ifftn( conj(dynInd(K,(N(6)+1)-i,6)).* part2)) ;
        temp = real(temp);
        gradTemp = dynInd(gradTemp, i, 6, temp);
    end
    gradRe = cat(dim,gradRe, E.Dl.Reg.lambda_susc * gradTemp / numel(E.Dl.D));
end

%%%SUM
if sep==0 && ~isempty(gradRe); gradFi = gradFi + multDimSum(gradRe,dim); end

end