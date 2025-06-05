function [rec] = assignSensitivities_SPACE(rec, recS, solve_espirit)

NS=size(recS.S);NS=NS(1:3);
NY=rec.Enc.FOVSize; 

%Filter
if solve_espirit == 1
    gibbsRinging=1; % used to be 1
    if gibbsRinging~=0 && gibbsRinging<=1;HY=buildFilter(2*NS,'tukeyIso',0.5,0,gibbsRinging,1);elseif gibbsRinging==-1;HY=buildFilter(2*NS,'CubicBSpline',[],0,1);else HY=[];end
    if ~isempty(HY)
        rec.S=filtering(recS.S,HY,1);
        rec.W=abs(filtering(recS.W,HY,1));
        HY=buildFilter(NS,'tukeyIso',1,0,gibbsRinging);
        rec.M=abs(filtering(abs(recS.x),HY)/20);           
    else            
        rec.S=recS.S;rec.W=recS.W;rec.M=recS.x;
    end
    
    %Map to right FOV
    MTS = recS.Par.Mine.APhiRec;MTy = rec.Par.Mine.APhiRec;
    rec.S = mapVolume (recS.S, ones(NY(1:3)) , MTS, MTy,[],[],'spline');
    %rec.W = mapVolume (recS.W, ones(NY(1:3)) , MTS, MTy,[],[],'spline');
    recS.M = gather(recS.M);
%     rec.M = mapVolume (recS.M, ones(NY(1:3)) , MTS, MTy,[],[],'spline');recS.M=single(recS.M>0.5); % ZN: previous gadgetron configuration for MPRAGE or SWI
    rec.M = mapVolume (rec.M, ones(NY(1:3)), MTS, MTy,[],[],'spline'); % ZN: current configuration for SPACE
    recS=[];

else
    gibbsRinging=1; % used to be 1
    if gibbsRinging~=0 && gibbsRinging<=1;HY=buildFilter(2*NS,'tukeyIso',0.5,0,gibbsRinging,1);elseif gibbsRinging==-1;HY=buildFilter(2*NS,'CubicBSpline',[],0,1);else HY=[];end
    if ~isempty(HY)
        rec.S=filtering(recS.S,HY,1);
%         rec.W=abs(filtering(recS.W,HY,1));
        HY=buildFilter(NS,'tukeyIso',1,0,gibbsRinging);
        rec.M=abs(filtering(abs(recS.x),HY)/20);           
    else            
        rec.S=recS.S;rec.W=recS.W;rec.M=recS.x;
    end
    
    %Map to right FOV
    MTS = recS.Par.Mine.APhiRec;MTy = rec.Par.Mine.APhiRec;
    rec.S = mapVolume (rec.S, ones(NY(1:3)) , MTS, MTy,[],[],'spline');
    %rec.W = mapVolume (recS.W, ones(NY(1:3)) , MTS, MTy,[],[],'spline');
    rec.M = mapVolume (rec.M, ones(NY(1:3)) , MTS, MTy,[],[],'spline');
    rec.ref_y = mapVolume(RSOS(recS.y),ones(NY(1:3)),MTS,MTy,[],[],'spline'); % by ZN, for masking
    recS=[];
end







