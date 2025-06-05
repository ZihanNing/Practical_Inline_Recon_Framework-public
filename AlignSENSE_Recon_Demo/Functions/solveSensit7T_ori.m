
function rec=solveSensit7T_ori(rec, recB)
    
    addpath(genpath('Functions_SENSE')); % ZN: to enable original espirit, use SENSE version support function folder
    
    rec.x=rec.y;
    % ZN: way of estB
    rec.Alg.bodyCoilEstimationMethod = 'BCC';
    % espirit related parameters
    rec.Alg.ThreeDEspirit=1;%To use 3D ESPIRIT
    %if rec.Alg.ThreeDEspirit;rec.Alg.parE.NC=2*ones(1,3);else rec.Alg.parE.NC=2*ones(1,2);end%Resolution (mm) of calibration area to compute compression
    if rec.Alg.ThreeDEspirit;rec.Alg.parE.NC=1.5*ones(1,3);else rec.Alg.parE.NC=1.5*ones(1,2);end%Resolution (mm) of calibration area to compute compression
    if rec.Alg.ThreeDEspirit;rec.Alg.parE.eigSc=[0.98 0.3];else rec.Alg.parE.eigSc=[0.9 0.3];end%Possible cut-off for the eigenmaps in soft SENSE, used the first for mask extraction by thresholding the eigenmaps%Default was [0.85 0.3]
    if rec.Alg.ThreeDEspirit;rec.Alg.parE.Ksph=200;else rec.Alg.parE.Ksph=50;end%Number of points for spherical calibration area, unless 0, it overrides other K's%Default was 200
    %if rec.Alg.ThreeDEspirit;rec.Alg.parE.factScreePoint=0.125;else rec.Alg.parE.factScreePoint=0.25;end%Factor over threshold computed with scree point
    if rec.Alg.ThreeDEspirit;rec.Alg.parE.factScreePoint=0.5;else rec.Alg.parE.factScreePoint=1;end%Factor over threshold computed with scree point
    %if rec.Alg.ThreeDEspirit;rec.Alg.parE.factScreePoint=2;else rec.Alg.parE.factScreePoint=4;end%Factor over threshold computed with scree point
    rec.Alg.parE.mirr=[0 0 0];%Whether to mirror along a given dimension%Default was [8 8 8]
    rec.Alg.parE.virCo=0;%Flag to use the virtual coil to normalize the maps%Default was 1, use 5 to get the Siemens contrast 
    %rec.Alg.parE.virCo=0;%Flag to use the virtual coil to normalize the maps%Default was 1, use 5 to get the Siemens contrast 
    rec.Alg.parE.factorBody=1;%Factor for virtual body coil in ESPIRIT (it was 1)
    rec.Alg.sensitivityEstimationMethod = 'ESPIRIT';

%     if rec.Alg.averagePolaritySensitivities==2;rec.S=mean(rec.S,11);end
    if ~isfield(rec,'B');estB=1;else estB=0;end
    rec.B=[];rec.M=[];rec.W=[];
    size(rec.y) % ZN: if you use seperate ACS, then no need to extract first echo
%     rec.S=dynInd(rec.S,1,8);%We extract first echo 
    
    for s=1:size(rec.y,11)
        S=dynInd(rec.y,s,11);
        if estB
            if strcmp(rec.Alg.bodyCoilEstimationMethod,'RSoS');B=sqrt(normm(S,[],4));
            elseif strcmp(rec.Alg.bodyCoilEstimationMethod,'CC');B=compressCoils(S,1);
            elseif strcmp(rec.Alg.bodyCoilEstimationMethod,'BCC');B=blockCompressCoils(S,1);
            elseif strcmp(rec.Alg.bodyCoilEstimationMethod,'CCPh');B=sqrt(normm(S,[],4)).*signz(compressCoils(S,1));
            else error('Virtual body coil estimation method %s not contemplated',rec.Alg.bodyCoilEstimationMethod);
            end
        end
        if strcmp(rec.Alg.sensitivityEstimationMethod,'Filter')
            S=gather(S);B=gather(B);AcqDelta=rec.Enc.S.AcqDelta;maskThS=rec.Alg.maskThS;
            %save('/mnt/nfs/home/lcordero/Matlab/Data/dat.mat','S','B','AcqDelta','maskThS');
            %exit
            [S,M]=solveSFilt(B,S,rec.Enc.S.AcqDelta,rec.Alg.maskThS);
            
        elseif strcmp(rec.Alg.sensitivityEstimationMethod,'Standard')
            rph.x=B;rph.y=S;rph.Enc.AcqVoxelSize=rec.Enc.S.AcqDelta;rph.Alg.parS=rec.Alg.parS;
            rph=solveS(rph);
            S=rph.S;M=rph.M;rph=[];
        elseif strcmp(rec.Alg.sensitivityEstimationMethod,'ESPIRIT')
            rph.x=B;rph.y=S;rph.Enc.AcqVoxelSize=rec.Enc.AcqVoxelSize;rph.Alg.parE=rec.Alg.parE;
            %rph.y=rph.y.*signz(B).*conj(signz(dynInd(rph.y,1,4)));%Phase correction
            %visSegment(angle(dynInd(rph.y,1,4)))
            rph=solveESPIRIT_lucilio(rph);
            S=rph.S;M=single(rph.W>rec.Alg.parE.eigSc(1));W=rph.Wfull;rph=[];
            rec.W=cat(11,rec.W,W);
        else
            error('Sensitivity estimation method %s not contemplated',rec.Alg.sensitivityEstimationMethod);
        end
        rec.S=dynInd(rec.y,s,11,S);
        rec.B=cat(11,rec.B,B);
        rec.M=cat(11,rec.M,M);
    end
%     if rec.Alg.averagePolaritySensitivities==1;rec.S=mean(rec.S,11);rec.B=mean(rec.B,11);end
    rec.M=max(rec.M,[],11);
    if useGPU;rec.S=gpuArray(rec.S);rec.x=gpuArray(rec.x);end
    rec.x=sum(conj(rec.S).*rec.x,4)./normm(rec.S,[],4);  
    rec.x=resPop(rec.x,11,[],5);rec.S=resPop(rec.S,11,[],5);rec.B=resPop(rec.B,11,[],5);
    if ~isempty(rec.W);rec.W=resPop(rec.W,11,[],5);end
%     if size(rec.xS,5)>1
%         rec.PS=sqrt(dynInd(rec.xS,2,5).*conj(dynInd(rec.xS,1,5)));rec.PS=gather(rec.PS);
%         rec.Geom.PS.MS=rec.Geom.S.MS;rec.Geom.PS.MT=rec.Geom.S.MT;
%     end
    rec.S=gather(rec.S);rec.x=gather(rec.x);rec.M = gather(rec.M);

end