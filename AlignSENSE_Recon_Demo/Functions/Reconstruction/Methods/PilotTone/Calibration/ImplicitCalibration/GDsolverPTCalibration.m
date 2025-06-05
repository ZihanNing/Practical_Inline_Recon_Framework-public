
function [E, PT, En] = GDsolverPTCalibration(y,E,EH,A,C,R,x,nIt, tol, deb, PT)

%GDSOLVERPTCALIBRATION   Performs a GD-based pseudoinverse reconstruction 
%   [E, En]=CGSOLVER(Y,E,EH,{X},{NIT},{TOL},{DEB})
%   * Y is the array of measures
%   * E is an encoding structure
%   * EH is a decoding structure
%   * {X} is an initial condition for the array to be reconstructed
%   * {NIT} is the maximum number of iterations
%   * {TOL} is the tolerance
%   * {DEB} indicates whether to print information about convergence

if nargin<5 || isempty(nIt);nIt=300;end
if nargin<6 || isempty(tol);tol=1e-5;end
if nargin<7 || isempty(deb);deb=0;end
if nargin<8 || isempty(PT);PT=[];end

if isempty(PT); warning('GDsolverPTCalibration:: No PT data provided: returned');return; end

%%% INITIALISE
%General parameters
gpu = isa(x,'gpuArray');
NX = size(x);
NStates  = E.NMs; 
NSegments  = E.NSe; %inlcudes the outer segment (shutter)

%(F)ISTA parameters
if isfield(PT,'Optim') && isfield(PT.Optim, 'lineSearch');   lineSearch=PT.Optim.lineSearch;else lineSearch=1; end%line search flag
if isfield(PT,'Optim') && isfield(PT.Optim, 'alphaList');   alphaList=PT.Optim.alphaList;else alphaList=10.^(-1* [-5:0.7:15]); end%line search step size
if isfield(PT,'Reg') && isfield(PT.Reg, 'lambda_sparse');   lambda_sparse=PT.Reg.lambda_sparse;else lambda_sparse=0; end
%if isfield(PT,'Optim') && isfield(PT.Optim, 'Mbg');     Mbg = PT.Optim.Mbg; else Mbg=1:NX(3);end

%Energy
if isfield(EH,'Mbe');Mbe = EH.Mbe;end %In contrast to motion estimation, use whole readout dimension - NOT USED FOR NOW
EH.Mbe = zeros([1 1 NX(3)],'like',real(x));EH.Mbe(E.nF)=1;%Because also used in LMSolver and gradPTCalibration
if isfield(EH,'We')&& ~isempty(EH.We)
    We = EH.We; %set aside
    %EH.Mbe= bsxfun(@times, We, EH.Mbe ); %Make sure these weights are also used in computeEnergy.m -> NOT HERE BECAUSE NOT WEIGHTED GRADIENT YET
end
EnPrevious = energyPTCalibration(y, x, E, EH);
En = EnPrevious; % keep track of energy

pTimeResGrOrig = PT.pTimeResGr;
AbOrig = PT.Ab;

%% ===========  TEMP SIMULATION
rng(22);
NChaTest=1;
paramSelect=1:6; paramsDiscard = dynInd(1:6,~ismember(1:6,paramSelect),2);
tol=0;
nIt = 30;

PT.NChaPT = NChaTest;
PT.pTimeResGr = dynInd(pTimeResGrOrig,1:NChaTest,1);

PT.Ab = dynInd(AbOrig,1:NChaTest,2);

PT.pTimeResGr = abs(PT.pTimeResGr);
PTGT=PT;
amp = 200;
PTGT.Ab = amp*rand(size(PT.Ab))/NChaTest;
PTGT.Ab = dynInd(PTGT.Ab, 4:6, 1, pi/1000/180*dynInd(PTGT.Ab,4:6,1));
PTGT.Ab = dynInd(PTGT.Ab, paramsDiscard, 1, 0);

%PTGT.Ab = ones(size(PT.Ab));
TGT=pilotTonePrediction (PTGT.Ab, PT.pTimeResGr, 'PT2T','backward', PT.NChaPT, size(E.Tr,6));
E.Tr=TGT; 
y = encode(x,E);
if isequal(paramSelect,1:3)
    PT.Ab = PTGT.Ab + amp/10*[ ones([1 3]) 0*ones([1 3])]'.*rand([1 6])';
elseif isequal(paramSelect,4:6)
    PT.Ab = PTGT.Ab + amp/10*[ 0*ones([1 3]) pi/180*ones([1 3])]'.*rand([1 6])';
else
    PT.Ab = PTGT.Ab + amp/10*[ 80*ones([1 3]) 0*pi/180*ones([1 3])]'.*rand([1 6])';
end
E.Tr = pilotTonePrediction (PT.Ab, PT.pTimeResGr, 'PT2T','backward', PT.NChaPT, size(E.Tr,6));
EnPrevious = energyPTCalibration(y, x, E, EH);
En=[];
visMotion(TGT,[],[],0,[],[],[],[],[],65);
visMotion(TGT-E.Tr,[],[],0,[],[],[],[],[],66);
linePlotA( y, E, EH, x, PT.Ab, PTGT.Ab, [], [], 100, 11, PT.pTimeResGr );

%%% ITERATE
for n=1:max(nIt)    
    Ab = PT.Ab; % Set original value aside so no need to duplicate E
    
    %%%GRADIENT
    grad = gradPTCalibration(y,E,EH,x,C,[],PT);%TODO: include weighting
    grad = dynInd( grad, paramsDiscard, 1, 0);
    %linePlotA(y, E, EH, x, PT.Ab, PT.Ab - 1*grad, [], [], 100, 11, PT.pTimeResGr );
% 
%    for alpha = alphaList
%         %%% Possible update 
%         update = - alpha*grad; %Only selecting one Taylor term (all others to zero)
% 
%         PT.Ab = Ab + update;
%         E.Tr = pilotTonePrediction(PT.Ab, PT.pTimeResGr, 'PT2T','backward', PT.NChaPT, size(E.Tr,6));
%         %E.Tr = restrictTransform(E.Tr);
%         
%         %%% Check possible energy reduction 
%         EnNew = energyPTCalibration(y,x,E,EH,[],0);         
%         if ~lineSearch || (EnNew < EnPrevious)
%             EnPrevious = EnNew;
%             break;
%         end 
% 
%         if alpha == alphaList(end) 
%             if deb >=2; fprintf('   GD iteration %d: no step size found\n   --> set to zero and GDsolver aborted\n', n);end
%             PT.Ab = Ab; % Restore with original value
%             E.Tr = pilotTonePrediction (PT.Ab , PT.pTimeResGr, 'PT2T','backward', PT.NChaPT, size(E.Tr,6));
%             return %don't rescale step size and leave GDsolver
%         end
%     end
    Etemp = LMsolver(y,E,x,C);
    E.Tr = Etemp.Tr;
    
    PT.Ab-PTGT.Ab
    multDimMea(abs(PT.Ab-PTGT.Ab))
    visMotion(TGT-E.Tr,[],[],0,[],[],[],[],[],67);
    
    %%%CHECK CONVERGENCE
%     En = cat(2, En, EnNew); 
%     if EnNew<tol;fprintf('   Convergence %.2f reached after %d iterations.\n',tol, n);break;end % reached convergence   
%     if n==max(nIt) && deb==2;fprintf('   GD solver terminated without reaching convergence\n');end
% 
%     if deb
%         h=figure(89); 
%         set(h,'color','w');
%         plot(1:length(En), En);
%         xlabel('Number of iterations [#]');
%         ylabel('Residuals [a.u.]');
%         title('Convergenve plot for GDsolverPTCalibration');
%     end
end

%%
%%%SET BACK PARAMETERS
if exist('Mbe', 'var');EH.Mbe=Mbe;end
if exist('We', 'var');EH.We=We;end

end