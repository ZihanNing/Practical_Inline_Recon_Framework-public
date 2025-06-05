
function rec = B0Algorithm(rec)

%B0ALGORITHM   Sets the default parameters for the B0 model used in the reconstruction
%   REC=B0ALGORITHM(REC)
%   * REC is a reconstruction structure without the algorithmic parameters. 
%   At this stage it may contain the naming information (.Names), the 
%   status of the reconstruction (.Fail), the .lab information (.Par), the 
%   fixed plan information (.Plan) and the dynamic plan information (.Dyn) 
%   ** REC is a reconstruction structure with filled algorithmic 
%   information (.Alg)
%

%General parameters
if ~isfield(rec.Alg,'parXB');rec.Alg.parXB=[];end
if ~isfield(rec.Alg.parXB,'useTaylor');rec.Alg.parXB.useTaylor=0;end%Whether to use the Taylor model
if ~isfield(rec.Alg.parXB,'orderTaylor');rec.Alg.parXB.orderTaylor=1;end%Order of Taylor model to use. Implementation goes up to 2nd order.
if ~isfield(rec.Alg.parXB,'weightedEst');rec.Alg.parXB.weightedEst=0;end%To use shot weigths from outlier detection (see DISORDER) for B0 estimation as well.
if ~isfield(rec.Alg.parXB,'redFOV');rec.Alg.parXB.redFOV=rec.Alg.parXT.redFOV;end%To exclude part of the FOV in the readout dimension when calculating the energy for the Taylor model update.
if ~isfield(rec.Alg.parXB,'nGD_Taylor');rec.Alg.parXB.nGD_Taylor=2;end %Number of iteration of the GD optimisation within the alternating optimisation.
if ~isfield(rec.Alg.parXB,'taylorDelay');rec.Alg.parXB.taylorDelay=6;end %To de-activate the GD optimisation for the first 6 outer iterations (as motion estimates might not be accurate yet).
if ~isfield(rec.Alg.parXB,'intraVoxDeph');rec.Alg.parXB.intraVoxDeph=0;end %To take into account intra-voxel dephasing in the forward model
if ~isfield(rec.Alg.parXB,'useSH');rec.Alg.parXB.useSH=0;end 
if ~isfield(rec.Alg.parXB,'unwrapTaylor');rec.Alg.parXB.unwrapTaylor=0;end 

%Optimisation parameters for the Gradient Descent based optimisation
if ~isfield(rec.Alg.parXB,'Optim');rec.Alg.parXB.Optim=[];end
if ~isfield(rec.Alg.parXB.Optim,'FISTA');rec.Alg.parXB.Optim.FISTA=0;end%Whether to use the FISTA algorithm for L1 regularisation
if ~isfield(rec.Alg.parXB.Optim,'RMSProp');rec.Alg.parXB.Optim.RMSProp=0;end%Whether to use the RMSProp algorithm version of the GD optimisation
if ~isfield(rec.Alg.parXB.Optim,'beta');rec.Alg.parXB.Optim.beta=0.9;end%Beta hyperparameter for the RMSProp algorithm
if ~isfield(rec.Alg.parXB.Optim,'alphaList');rec.Alg.parXB.Optim.alphaList=10.^(-1* [-3:0.7:6]);end%Line search step size

%Regularisation
if ~isfield(rec.Alg.parXB,'Reg');rec.Alg.parXB.Reg=[];end
if ~isfield(rec.Alg.parXB.Reg,'lambda_sparse');rec.Alg.parXB.Reg.lambda_sparse=0;end%Weight for the L1 regularisation of the Taylor maps
if ~isfield(rec.Alg.parXB.Reg,'lambda_smooth');rec.Alg.parXB.Reg.lambda_smooth=0;end%Weight for the L2 regularisation of the finite differences of the Taylor maps
if ~isfield(rec.Alg.parXB.Reg,'lambda_ridge');rec.Alg.parXB.Reg.lambda_ridge=0;end%Weight for the ridge regularisation of the Taylor maps

%Constrain
if ~isfield(rec.Alg.parXB,'C');rec.Alg.parXB.C=[];end
if ~isfield(rec.Alg.parXB.C,'filterTaylor');rec.Alg.parXB.C.filterTaylor=[];end %Filter structure to constrain Taylor maps. Structure should contain spacing filterTaylor.sp and Gibbs ringing factor filterTaylor.gibbsRinging
if ~isfield(rec.Alg.parXB.C,'DupLim');rec.Alg.parXB.C.DupLim=inf;end%Maximum absolute update of the Taylor maps that will be thresholded. In units of Hz/degree.
if ~isfield(rec.Alg.parXB.C,'DLim');rec.Alg.parXB.C.DLim=inf;end%Maximum value of the Taylor maps that will be thresholded. In units of Hz/degree.





