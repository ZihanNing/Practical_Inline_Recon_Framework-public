function rec=PEregrid(rec,di)

%PEREGRID   Forwards and backwards regridding of multiple phase encode information
%   REC=PEREGRID(REC)
%   * REC is a reconstruction structure. At this stage it may contain the naming information (rec.Names), the status of the reconstruction (.Fail), the
%   .lab information (rec.Par), the fixed plan information (rec.Plan), the dynamic plan information (rec.Dyn), the data information
%   (rec.(rec.Plan.Types)), the information for correction (rec.Corr.(rec.Plan.Types)), the informaton for further sorting 
%   (rec.Assign(rec.Plan.Types)) and the encoding information (rec.Enc)
%   * REC is a reconstruction structure with transformed (or not, depending on TYP) coordinates
%   * DI is the regridding direction (it defaults to forward)
%

if nargin<2 || isempty(di);di=1;end
if ~isfield(rec,'x');return;end

gpu=isa(rec.x,'gpuArray');

typ2Rec=rec.Dyn.Typ2Rec;
if ~isempty(rec.Par.Mine.Signs)    
    if isfield(rec,'Mprov')
        if gpu;rec.Mprov=gpuArray(rec.Mprov);end
        rec.M=cat(4,rec.M,rec.Mprov);
    end
    NPE=size(rec.Par.Mine.Signs,3);
    indPEV=cell(1,NPE);
    voInds=rec.Plan.Par2Rec{11}+1;%Volumetric indexes
    if ~isempty(rec.Par.Mine.diInfo) && any(voInds>size(rec.Par.Mine.diInfo,1))
        fprintf('More volumes in the data than in the direction file for %s\n',rec.Names.Name);
        fprintf('DATA NOT WRITTEN. LAST MESSAGE SHOULD PROVIDE THE REASON\n');
        rec.Fail=1;return;
    end
    if ~isempty(rec.Par.Mine.diInfo) && size(rec.Par.Mine.diInfo,2)>=6;voInds=voInds(rec.Par.Mine.diInfo(voInds,6)==1);end
    if isfield(rec,'x');voInds=voInds(1:size(rec.x,4));end
    if isfield(rec.Par.Mine,'InitVol');voInds=voInds+rec.Par.Mine.InitVol;end
    if di
        %RESAMPLING AND MASKING
        if any(abs(rec.Par.Mine.Signs(1,2,:)))==1
            M=size(rec.x);M=M(1:2);
            O=rec.Alg.OverDec(1:2);O(O>1)=1;O=abs(O);
            R=round(max(O)*M./O);
            N=max(R)*ones(1,2);
            resampleAndMask; 
        end

        for n=typ2Rec';datTyp=rec.Plan.Types{n};rec.(datTyp)=gather(rec.(datTyp));end
                
        %PERMUTE AND FLIP
        for n=typ2Rec';datTyp=rec.Plan.Types{n};
            if gpu;rec.(datTyp)=gpuArray(rec.(datTyp));end
            c=1;
            for p=1:NPE
                %if ~ismember(n,[7:9 19:21]);indPEV=find(rec.Par.Mine.diInfo(voInds,5)==rec.Par.Mine.pedsUn(p));di=4;elseif ismember(n,[9 19:21]);indPEV=p;di=4;else indPEV=p;di=5;end
                if ~ismember(n,[7:9 19:21 27]);indPEV=find(rec.Par.Mine.diInfo(voInds,5)==p);di=4;
                elseif ismember(n,19:21);indPEV=p;di=4;
                elseif ismember(n,9);indPEV=find(rec.Par.Mine.pedsUn==p);di=4;
                elseif size(rec.(datTyp),5)==NPE;indPEV=p;di=5;
                elseif ismember(p,rec.Par.Mine.pedsUn);indPEV=c;di=5;c=c+1;
                else indPEV=[];di=5;
                end 
                if ~isempty(indPEV)
                    perm=(1:ndims(rec.(datTyp)))';        
                    perm(1:2)=rec.Par.Mine.Signs(:,:,p)*perm(1:2);
                    for s=1:size(rec.(datTyp),6)
                        rec.(datTyp)=dynInd(rec.(datTyp),{indPEV,s},[di 6],permute(dynInd(rec.(datTyp),{indPEV,s},[di 6]),abs(perm')));
                        for m=1:2
                            if perm(m)<0
                                rec.(datTyp)=dynInd(rec.(datTyp),{indPEV,s},[di 6],circshift(ifftGPU(ifftshift(flip(fftshift(fftGPU(dynInd(rec.(datTyp),{indPEV,s},[di 6]),m),m),m),m),m),-1,m));%It may not be necessary to perform so many fftshifts
                            end
                        end
                    end
                end
            end
            magnitudeAndThreshold;
            if n~=8;rec.(datTyp)=gather(rec.(datTyp));end
        end        
    else
        %FLIP AND PERMUTE      
        for n=typ2Rec';datTyp=rec.Plan.Types{n};
            c=1;
            for p=1:NPE
                if ~ismember(n,[7:9 19:21 27]);indPEV=find(rec.Par.Mine.diInfo(voInds,5)==p);di=4;
                elseif ismember(n,19:21);indPEV=p;di=4;
                elseif ismember(n,9);indPEV=find(rec.Par.Mine.pedsUn==p);di=4;
                elseif size(rec.(datTyp),5)==NPE;indPEV=p;di=5;
                elseif ismember(p,rec.Par.Mine.pedsUn);indPEV=c;di=5;c=c+1;
                else indPEV=[];di=5;
                end                 
                if ~isempty(indPEV)                   
                    perm=(1:ndims(rec.(datTyp)))';
                    perm(1:2)=rec.Par.Mine.Signs(:,:,p)*perm(1:2); 
                    for s=1:size(rec.(datTyp),6)
                        for m=1:2
                            if perm(m)<0;rec.(datTyp)=dynInd(rec.(datTyp),{indPEV,s},[di 6],ifftGPU(ifftshift(flip(fftshift(fftGPU(circshift(dynInd(rec.(datTyp),{indPEV,s},[di 6]),1,m),m),m),m),m),m));end%It may not be necessary to perform so many fftshifts
                        end
                        rec.(datTyp)=dynInd(rec.(datTyp),{indPEV,s},[di 6],permute(dynInd(rec.(datTyp),{indPEV,s},[di 6]),abs(perm')));
                    end
                end
            end
            magnitudeAndThreshold;
        end        
        %RESAMPLING AND MASKING
        if any(abs(rec.Par.Mine.Signs(1,2,:)))==1
            R=rec.Enc.FOVSize(1:2);
            O=rec.Alg.OverDec(1:2);O(O>1)=1;O=abs(O);
            N=round(max(O)*R./O);
            resampleAndMask;
        end        
    end
    if isfield(rec,'Mprov');rec=rmfield(rec,'Mprov');rec.M=dynInd(rec.M,2,4);end
end

function resampleAndMask     
     for n=typ2Rec';datTyp=rec.Plan.Types{n}; 
         if n==8;rec.(datTyp)=single(abs(rec.(datTyp))>1e-12);end            
         if max(M)>1.5*min(M);rec.(datTyp)=gather(rec.(datTyp));gpu=0;end%Potential memory problems
         if di
            rec.(datTyp)=resampling(rec.(datTyp),R,2);%Oversampling
            rec.(datTyp)=resampling(rec.(datTyp),N);%Resampling          
         else
            rec.(datTyp)=resampling(rec.(datTyp),N);%Resampling
            rec.(datTyp)=resampling(rec.(datTyp),R,2);%Oversampling
         end
         magnitudeAndThreshold;
     end     
     for n=typ2Rec';datTyp=rec.Plan.Types{n};
         if ~ismember(n,[7:8 12 19:21 27])
             c=1;
             for p=1:NPE
                 if ~ismember(n,[7:9 19:21 27]);indPEV=find(rec.Par.Mine.diInfo(voInds,5)==p);di=4;
                 elseif ismember(n,19:21);indPEV=[p 2];di=4:5;
                 elseif ismember(n,9);indPEV=find(rec.Par.Mine.pedsUn==p);di=4;
                 elseif size(rec.(datTyp),5)==NPE;indPEV=p;di=5;
                 elseif ismember(p,rec.Par.Mine.pedsUn);indPEV=c;di=5;c=c+1;
                 else indPEV=[];di=5;
                 end               
                 if ~isempty(indPEV)
                    if size(rec.M,5)==4;rec.(datTyp)=dynInd(rec.(datTyp),indPEV,di,bsxfun(@times,dynInd(rec.(datTyp),indPEV,di),dynInd(rec.M,[1 p],4:5)));
                    elseif size(rec.M,5)==length(rec.Par.Mine.pedsUn);rec.(datTyp)=dynInd(rec.(datTyp),indPEV,di,bsxfun(@times,dynInd(rec.(datTyp),indPEV,di),dynInd(rec.M,[1 indPEV],4:5)));
                    else error('Unknown PE configuration at PEregrid. Number of reformatted masks: %d / Number of PEs: %d',size(rec.M,5),length(rec.Par.Mine.pedsUn));
                    end
                 end
             end
         end
     end
end

function magnitudeAndThreshold
     if ismember(n,[8:10 19:21 27]);rec.(datTyp)=abs(rec.(datTyp));end
     if n==8;rec.(datTyp)=single(rec.(datTyp)>0.5);end
end

end