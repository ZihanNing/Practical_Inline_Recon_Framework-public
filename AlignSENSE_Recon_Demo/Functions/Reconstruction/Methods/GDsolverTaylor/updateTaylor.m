
function E = updateTaylor (E, x, N, res, parXB, Labels)

%UPDATETAYLOR  initialises and updates the encoding structure E with the
%necessary information to estimate the Taylor model at different resolutions.
%   [EN,W]=UPDATETAYLOR(E,X,N,RES,PARXB,LABELS)
%   * E is the encoding structure
%   * X is the reconstructed data. Currently only used to determine the class to initialise Taylor maps with.
%   * N is the size of the reconstructed data
%   * RES is the resolution of the reconstruction level
%   * PARXB is the structure containing the parameters for B0 estimation
%   * LABELS is the structure containing the labels of the acquisition
%

if length(res)>1; res = res(1);end

if ~isfield(E, 'Dl')%Dl is the structure within E that takes care of the Taylor model
    fprintf('updateTaylor:: Creating Taylor coeffients for B0 modelling\n')
    
    % Initialise
    if ~isfield( parXB, 'orderTaylor') || isempty(parXB.orderTaylor); parXB.orderTaylor = 1;end
    yawIdx = 4; rollIdx=5; pitchIdx=6;%Reverse order as in sincRigidTransformation convention/parameters in RAS
    [E.Dl.f, E.Dl.fDeriv] = higherOrderFunHandle(parXB.orderTaylor,  yawIdx, rollIdx , pitchIdx);
    nMaps = numel( E.Dl.f( ones(1,1,1,1,1,6))); %See how many of the 6 motion parameters are store in f(T)
    E.Dl.D =  zeros([N , 1, 1, nMaps], 'like',real(x) );%The Taylor maps are stored in the 6th dimension
     
    % Parameters
    if parXB.orderTaylor == 1; E.Dl.d = 5:6; elseif parXB.orderTaylor>1; E.Dl.d = 4:6;end %Indices that extract the relevant motion parameters for the Taylor model.
    E.Dl.TE = Labels.TE * 1e-3;%in seconds
    
    % Regularisation
    if isfield(parXB,'Reg'); E.Dl.Reg = parXB.Reg; else E.Dl.Reg =[];end
    
    % Optimisation
    if isfield(parXB,'Optim'); E.Dl.Optim = parXB.Optim; else E.Dl.Optim =[];end
    if ~isfield( E.Dl.Optim, 'alphaList'); E.Dl.Optim.alphaList=10.^(-1* [-2:0.7:6]);end
    if strcmp(Labels.FatShiftDir,'F');E.Dl.Optim.Mbg=1:floor((1-parXB.redFOV)*N(3));elseif strcmp(Labels.FatShiftDir,'H');E.Dl.Optim.Mbg=floor(1+parXB.redFOV*N(3)):N(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.Dl.Optim.Mbg=1:N(3);end
    
    % Constrain 
    if isfield(parXB,'C'); E.Dl.C = parXB.C; else E.Dl.C =[];end
    if ~isempty(parXB.C.filterTaylor)
        gpuIn=isa(x,'gpuArray');
        if isfield ( parXB.C.filterTaylor, 'sp'); E.Dl.C.filterTaylor.sp = parXB.C.filterTaylor.sp; else E.Dl.C.filterTaylor.sp = 10;end%Store on E.Dl.C so that it can be used at other resolution levels
        if isfield ( parXB.C.filterTaylor, 'gibbsRinging'); E.Dl.C.filterTaylor.gibbsRinging = parXB.C.filterTaylor.gibbsRinging; else E.Dl.C.filterTaylor.gibbsRinging = 0;end
        E.Dl.C.H = buildFilter(2*N(1:3),'tukeyIso',(res/E.Dl.C.filterTaylor.sp)*ones(1,3),gpuIn,E.Dl.C.filterTaylor.gibbsRinging,1);%Filter in Cosine domain 
    end
    if isfield(parXB.C,'deSH'); E.Dl.C.deSH = parXB.C.deSH; else E.Dl.C.deSH = 0; end
    if isfield(parXB.C,'unwrapTaylor'); E.Dl.C.unwrapTaylor = parXB.C.unwrapTaylor; else E.Dl.C.unwrapTaylor = 0; end
    if isfield(parXB.C,'denoiseTaylor'); E.Dl.C.denoiseTaylor = parXB.C.denoiseTaylor; else E.Dl.C.denoiseTaylor = 0; end
             
else
    fprintf('updateTaylor:: Updating Taylor coeffients for B0 modelling\n')

    %%% Resample
    E.Dl.D=resampling(E.Dl.D,N,0,2*ones(1,3));
        
    %%% Constrain 
    if ~isempty(parXB.C.filterTaylor)
        gpuIn=isa(x,'gpuArray');
        E.Dl.C.H = buildFilter(2*N(1:3),'tukeyIso',(res/E.Dl.C.filterTaylor.sp)*ones(1,3),gpuIn,E.Dl.C.filterTaylor.gibbsRinging,1);%Filter in Cosine domain    
    end
    %%% Change order if needed
    nMaps = size(E.Dl.D,6);% How many in model now
    
    [E.Dl.f, E.Dl.fDeriv] = higherOrderFunHandle(parXB.orderTaylor);
    nMapsModel = numel( E.Dl.f( ones(1,1,1,1,1,6)));% How many there should be in (orderTaylor might have changed)
    if strcmp(Labels.FatShiftDir,'F');E.Dl.Optim.Mbg=1:floor((1-parXB.redFOV)*N(3));elseif strcmp(Labels.FatShiftDir,'H');E.Dl.Optim.Mbg=floor(1+parXB.redFOV*N(3)):N(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.Dl.Optim.Mbg=1:N(3);end
    
    if nMaps ~= nMapsModel
        fprintf('updateTaylor:: Changing the Taylor expansion applied to motion-induced B0 model.\n')
        Dnew = zeros([N , 1, 1, nMapsModel], 'like', real(E.Dl.D) );
        Dnew = dynInd( Dnew, 1:nMaps, 6, E.Dl.D);%Replace first terms with the ones coming from the previous order
        E.Dl.D = Dnew;
    end
    
end

