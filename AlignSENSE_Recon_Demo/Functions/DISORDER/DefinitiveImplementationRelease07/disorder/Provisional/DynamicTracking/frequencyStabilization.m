function rec=frequencyStabilization(rec,appl)

% FREQUENCYSTABILIZATION performs temporal tracking of volumes
%   REC=FREQUENCYSTABILIZATION(REC)
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (rec.Names), the status of the reconstruction 
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), and the data 
%   information (rec.(rec.Plan.Types))
%   * REC is a reconstruction structure with information in tracked coordinates
%

if nargin<2 || isempty(appl);appl=0;end

if rec.Fail || rec.Par.Mine.Modal~=9;return;end%Only run if modality is fMRI

NE=length(rec.Par.Labels.TE);
if rec.Par.Labels.TE(end)<rec.Par.Labels.TE(1);NE=1;end%Sometimes it appears there is a second TE while there isn't
if NE==1;TE=rec.Par.Labels.TE(1);elseif NE==2;TE=diff(rec.Par.Labels.TE);else error('Distortion reversal only working for single or double echo data');end     

gpu=rec.Dyn.GPU;
if appl==0
    if gpu;gpuF=2;else gpuF=0;end
    if gpu;rec.x=gpuArray(rec.x);end

    AdHocArray=rec.Par.Mine.AdHocArray;
    assert(AdHocArray(1)~=102,'SIR reconstruction is not supported any longer'); 

    x=abs(rec.x);
    rec.x=gather(rec.x);
    N=size(x);N(end+1:4)=1;
    if NE>1 && numDims(rec.x)<=4;x=reshape(x,[N(1:3) NE N(4)/NE]);end
    if NE==1 && numDims(rec.x)>=5;x=reshape(x,[N(1:3) prod(N(4:end))]);end

    ND=numDims(x);ND=max(ND,4);
    if ND==5;x=permute(x,[1 2 3 5 4]);end

    subPyr=4;
    mirr=[0 0 0 1];mirr(end+1:ND)=0;
    T=computeFrequencyStabilization(x,subPyr,mirr);
    T{2}(:);%OVERRIDDING THE FUNCTION!
    [rec.w,ESFOV]=applyFrequencyStabilization(rec.x,T,rec,gpu);
    rec.w=gather(rec.w);
    rec.Dyn.Typ2Rec(rec.Dyn.Typ2Rec==12)=[];
    if isfield(rec,'E')
        rec.Dyn.Typ2Wri(28)=1;
        rec.E=gather(abs(applyFrequencyStabilization(rec.E,T,rec,gpu)));%KEY INSTRUCTION
    end
    %TE
    rec.F=convertB0Field(gather(T{2}),TE,ESFOV,'pix','Hz',N(2));
    rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,24);
    rec.Dyn.Typ2Wri(23:24)=1;
else        
    %ES
    NPE=length(rec.Enc.kGrid{2});
    if isfield(rec.Par.Mine,'ES')
        ES=rec.Par.Mine.ES;
    else    
        WFS=rec.Par.Labels.WFS;
        ES=1000*(WFS/(3*42.576*3.35*(NPE+1)));%Using NPE+1 is really strange but has provided best match (not complete with ES read from GoalC parameters...   
        %echo spacing in msec = 1000 * (water-fat shift (per pixel)/(water-fat shift (in Hz) * echo train length))
        %    echo train length (etl) = EPI factor + 1
        %    water-fat-shift (Hz) = fieldstrength (T) * water-fat difference (ppm) * resonance frequency (MHz/T)
        %    water-fat difference (ppm) = 3.35 [2]
        %    resonance frequency (MHz/T) = 42.576  (for 1H; see Bernstein pg. 960)
    end
    %fprintf('Echo spacing (ms): %.2f\n',ES);
    N=size(rec.E);
    ESFOV=ES*NPE/N(2);

    T{2}=rec.F;
    T{2}=convertB0Field(T{2},TE,ESFOV,'Hz','pix',N(2));
    rec.E=gather(applyFrequencyStabilization(rec.E,T,rec,gpu));
end
