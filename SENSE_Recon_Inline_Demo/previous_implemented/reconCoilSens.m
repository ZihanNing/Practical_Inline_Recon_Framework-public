function rec=reconCoilSens(rec)
    rec.xS=rec.S;
    if rec.Alg.averagePolaritySensitivities==2;rec.S=mean(rec.S,11);end
    if ~isfield(rec,'B');estB=1;else estB=0;end
    rec.B=[];rec.M=[];rec.W=[];
    size(rec.S)
    rec.S=dynInd(rec.S,1,8);%We extract first echo 
    
    for s=1:size(rec.S,11)
        S=dynInd(rec.S,s,11);
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
            rph.x=B;rph.y=S;rph.Enc.AcqVoxelSize=rec.Enc.S.AcqDelta;rph.Alg.parE=rec.Alg.parE;
            %rph.y=rph.y.*signz(B).*conj(signz(dynInd(rph.y,1,4)));%Phase correction
            %visSegment(angle(dynInd(rph.y,1,4)))
            rph=solveESPIRIT(rph);
            S=rph.S;M=single(rph.W>rec.Alg.parE.eigSc(1));W=rph.Wfull;rph=[];
            rec.W=cat(11,rec.W,W);
        else
            error('Sensitivity estimation method %s not contemplated',rec.Alg.sensitivityEstimationMethod);
        end
        rec.S=dynInd(rec.S,s,11,S);
        rec.B=cat(11,rec.B,B);
        rec.M=cat(11,rec.M,M);
    end
    if rec.Alg.averagePolaritySensitivities==1;rec.S=mean(rec.S,11);rec.B=mean(rec.B,11);end
    rec.M=max(rec.M,[],11);
    if useGPU;rec.S=gpuArray(rec.S);rec.xS=gpuArray(rec.xS);end
    rec.xS=sum(conj(rec.S).*rec.xS,4)./normm(rec.S,[],4);  
    rec.xS=resPop(rec.xS,11,[],5);rec.S=resPop(rec.S,11,[],5);rec.B=resPop(rec.B,11,[],5);
    if ~isempty(rec.W);rec.W=resPop(rec.W,11,[],5);end
    if size(rec.xS,5)>1
        rec.PS=sqrt(dynInd(rec.xS,2,5).*conj(dynInd(rec.xS,1,5)));rec.PS=gather(rec.PS);
        rec.Geom.PS.MS=rec.Geom.S.MS;rec.Geom.PS.MT=rec.Geom.S.MT;
    end
    rec.S=gather(rec.S);rec.xS=gather(rec.xS);
end