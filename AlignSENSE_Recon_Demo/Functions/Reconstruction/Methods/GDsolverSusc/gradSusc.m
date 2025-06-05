
function [gradFi, gradRe] = gradSusc(y,x,E,EH,W,sep)

%GRADSUSC  Computes the gradient of a least squares 
%problem with different regularisations w.r.t. the susceptibility model of pose-dependent B0 variations. 
%   [EN,W]=GRADSUSC(Y,X,E,EH,{W},{SEP})
%   * Y is the measured data
%   * X is the reconstructed data
%   * E is the encoding structure
%   * EH is the decoding structure (used for weighted least squares)
%   * {W} is the weight to use for certain regularisers
%   * {SEP} whether to separate the data fidelity loss and the regularisatoin gradient terms
%

if nargin < 5 || isempty(W); W = ones(size(x), 'like',real(x)); applyWeight = 0; else applyWeight = 1; end %applyWeight used for efficiency lambda_susc
if nargin < 6 || isempty(sep); sep=(nargout>1);end
if isfield(E.Ds,'Optim') && isfield(E.Ds.Optim, 'BlSzStates');BlSzStates = E.Ds.Optim.BlSzStates; else BlSzStates=E.bS(1) ;end% Block size for coefficient gradient calculation

%%%INITIALISE
NStates  = E.NMs; 
N = size(E.Ds.chi);
    
%%%DATA FIDELITY
gradFi = zeros(N, 'like',E.Ds.chi);
res=y-encode(x,E);
for a=1:BlSzStates:NStates %Chunk states otwerwise is has size [Nx Nstates] and this runs out of memmory depending on NX and NStates
    E.oS(1) = a; E.bS(1) = BlSzStates; E.dS(1) = min(a+BlSzStates-1,NStates);
    vA=a:min(a+E.bS(1)-1,NStates);

    notsum = 1; % if 0, decode step will sum over motion states
    [~,EHres] = decode(res,EH,E,notsum); % motion states in the 5th dimension
    
    xToMultiply =  x;%exp(1i*angle(x));%Play with this term to investigate effect of signal voids on D update
    if a==1 && ~isequal(x,xToMultiply);fprintf('gradSusc:: taking magnitude out of x.\n');end
    
    EHres = bsxfun(@times, EHres, conj(xToMultiply) );
    
    %%% Padding to match 
    EHres = padArrayND(EHres,E.Ds.padDim,1,0,'both'); 

    %%% Fourier Transform
    for i=1:3; EHres = fftGPU(EHres,i);end
    
    %%% Apply Kernels
    if ~isempty(E.Ds.kernelStruct.Kroll) && ~isempty(E.Ds.kernelStruct.Kpitch)
        rotPar = E.Ds.f(dynInd(E.Ds.Tr, {vA}, 5)); %Rotation parameters used in the Taylor model
        kernelToApply = bsxfun(@times, cat( 6, E.Ds.kernelStruct.Kroll,E.Ds.kernelStruct.Kpitch) , rotPar ); %might have last segment that is issue here
        kernelToApply = multDimSum( kernelToApply, 6);
    else
        kernelToApply=0;
    end
    
    if ~isempty(E.Ds.kernelStruct.Kequilibrium)
        kernelToApply =  bsxfun(@plus, kernelToApply , E.Ds.kernelStruct.Kequilibrium );%Add in k-space
    end
    
    %Filter
    kernelToApply =  bsxfun(@times, fftshift(E.Ds.kernelStruct.H),kernelToApply );
    for i=1:3; kernelToApply = ifftshiftGPU(kernelToApply,i);end%EHres not shifted so shift back
    
    EHres = bsxfun(@times, EHres, kernelToApply ) ;
        
    %%% Inverse Fourier Transform
    for i=1:3; EHres = ifftGPU(EHres,i);end

    %%% Inverse Padding 
    %EHres = padArrayND(EHres,E.Ds.padDim,0,[],'both'); 
    
    gradFi = gradFi + multDimSum(real( +1i * 2*pi*E.Ds.TE * EHres) , 5 );% / numel(y) ; % numel(y) since in objective function relative error 
end
   
%%%REGULARISATION - only L2 gradients as L1 incorporated in (F)ISTA
gradRe = [];
dim=4; %dimension where to store the different regularisation terms

if isfield(E.Ds,'Reg') && isfield(E.Ds.Reg,'lambda_smooth') && ~isequal(E.Ds.Reg.lambda_smooth,0)%Smoothness
    gradRe = cat(dim,gradRe, E.Ds.Reg.lambda_smooth * ( FiniteDiff( W .*  FiniteDiff(E.Ds.chi,1),1 ,0) + ...
                                                        FiniteDiff( W .*  FiniteDiff(E.Ds.chi,2),2 ,0) + ...
                                                        FiniteDiff( W .*  FiniteDiff(E.Ds.chi,3),3 ,0) ) / numel(E.Ds.chi)) ;
end

if isfield(E.Ds,'Reg') && isfield(E.Ds.Reg,'lambda_ridge') && ~isequal(E.Ds.Reg.lambda_ridge,0)%Ridge 
    gradRe = cat(dim,gradRe, E.Ds.Reg.lambda_ridge * E.Ds.chi /numel(E.Ds.chi));
end

if isfield(E.Ds,'Reg') && isfield(E.Ds.Reg,'lambda_tissueConfidence') && ~isequal(E.Ds.Reg.lambda_tissueConfidence,0)%Ridge 
    gradRe = cat(dim,gradRe, E.Ds.Reg.lambda_tissueConfidence * W .* (E.Ds.chi - E.Ds.Reg.tissueSusc) /numel(E.Ds.chi));
end

%%%SUM
if sep==0 && ~isempty(gradRe); gradFi = gradFi + multDimSum(gradRe,dim); end

end