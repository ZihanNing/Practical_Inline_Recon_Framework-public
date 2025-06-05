function rec=reverseDistortion(rec,TE,norm,magn,appl)

% REVERSEDISTORTION reverses the distortion from B0 field
%   REC=REVERSEDISTORTION(REC,{TE},{NORM},{MAGN})
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (rec.Names), the status of the reconstruction 
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), and the data 
%   information (rec.(rec.Plan.Types))
%   * {TE} is to input a TE for field expressed in radians, otherwise it is
%   obtained from the reconstruction structure metadata
%   * {NORM} is to normalize the field estimates
%   * {MAGN} is to take the magnitude of the output data
%   ** REC is a reconstruction structure with information in undistorted 
%   coordinates (rec.u) 
%

if rec.Fail;return;end

if nargin<3 || isempty(norm);norm=0;end
if nargin<4 || isempty(magn);magn=1;end
if nargin<5 || isempty(appl);appl=0;end

fieldInHz=(nargin>=2) && ~isempty(TE);

NE=length(rec.Par.Labels.TE);
if rec.Par.Labels.TE(end)<rec.Par.Labels.TE(1);NE=1;end%Sometimes it appears there is a second TE while there isn't 
if NE==1 && rec.Par.Mine.Modal==10 && ~norm;return;end

gpu=rec.Dyn.GPU;

if appl
    rec.Dyn.Typ2Rec=[11;28];
    typInt=28;
    if gpu;rec.E=gpuArray(rec.E);end
else
    typInt=13;
    if gpu;rec.u=gpuArray(rec.u);end
end
if gpu;rec.B=gpuArray(rec.B);end

if (rec.Dyn.Typ2Wri(11)==1 && (rec.Alg.parU.useUndist && ~norm)) || (rec.Alg.DistoSensi && norm) || appl    
    %PERMUTE
    typ2Rec=rec.Dyn.Typ2Rec;
    if ~norm
        for n=typ2Rec';datTyp=rec.Plan.Types{n};
            rec.(datTyp)=permute(rec.(datTyp),[1 2 3 5 4]);
        end
    end

    overs=2;%Oversampling in resolution for reversing the gradient
    %TE
    if ~fieldInHz       
        if NE==1;TE=rec.Par.Labels.TE(1);else TE=rec.Par.Labels.TE(2)-rec.Par.Labels.TE(1);end%Currently for more than one echoe we assume the difference is the distance between the first and the second
    end    

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
        %From the PPE code
        %define MGG_GAMMA_1H (42577.46778)
        %define MGG_PPM_WATER_FAT_SHIFT (3.4)
        %Then: 
        %ES=1000*(WFS/(3*42.57746778*3.4*(NPE+1)));
        %It may also perhaps be necessary to fix the field to real one, as
        %it may be different from nominal 3T   
        %In summary (see below)
        %General Philips 3T: 3*42.576*3.35=427.8888
        %Achieva scanners: 434.215
        %From PPE: 434.2902
    end
    fprintf('Echo spacing (ms): %.2f\n',ES);
    N=size(rec.(rec.Plan.Types{typInt}));
    NRes=N;NRes(2)=overs*N(2);NRes(end+1:4)=1;
    ESFOV=ES*NPE/N(2);
    NB=size(rec.B);NB(end+1:4)=1;

    %FIELD FILTER
    %if fieldInHz;H=buildFilter(N(1:2),'tukeyIso',0.1*ones(1:2),gpu,rec.Alg.parU.GibbsRingi);else H=buildFilter(N(1:2),'tukeyIso',ones(1:2),gpu,rec.Alg.parU.GibbsRingi);end%This probably producing best results in fMRI    
    if fieldInHz;H=buildFilter(N(1:2),'tukeyIso',ones(1:2),gpu,rec.Alg.parU.GibbsRingi);else H=buildFilter(N(1:2),'tukeyIso',ones(1:2),gpu,rec.Alg.parU.GibbsRingi);end
    if ~appl;rec.B=filtering(rec.B,H);end;H=[];
    %CONVERT B0 FROM RAD TO 1/FOV
    if ~fieldInHz;orD='rad';else orD='Hz';end
    rec.B=convertB0Field(rec.B,TE,ESFOV,orD,'res');

    %APODIZE
    if fieldInHz
        %HA=buildFilter(N(2),'tukey',1,gpu,1);%This was shown probably to produce best results in fMRI
        HA=buildFilter(N(2),'tukey',1,gpu,0.1);
        HA=fftshift(permute(HA,[2 1]));
        if ~appl;rec.B=bsxfun(@times,rec.B,HA);end;HA=[];
    end

    %INVERT THE FIELD
    if norm;nor=rec.(rec.Plan.Types{typInt});nor(:)=1;end
    NT=size(rec.(rec.Plan.Types{typInt}),5);
    if NT>1000
        blSzN=100;
        for s=1:blSzN:NT;vS=s:min(s+blSzN-1,NT);             
            if size(rec.B,5)==NT;Bs=dynInd(rec.B,vS,5);else Bs=rec.B;end
            NRes(5)=length(vS);N(5)=length(vS);
            rec.(rec.Plan.Types{typInt})=dynInd(rec.(rec.Plan.Types{typInt}),vS,5,resampling(reverseB0Field(resampling(dynInd(rec.(rec.Plan.Types{typInt}),vS,5),NRes),real(resampling(repmat(Bs,[ones(1,3) NRes(4)/NB(4)]),NRes(1:4)))),N));
            if isfield(rec,'E') && ~appl;rec.E=dynInd(rec.E,vS,5,abs(resampling(reverseB0Field(resampling(dynInd(rec.E,vS,5),NRes),real(resampling(repmat(Bs,[ones(1,3) NRes(4)/NB(4)]),NRes(1:4)))),N)));end%KEY INSTRUCTION
            if norm;nor=dynInd(nor,vS,5,resampling(reverseB0Field(resampling(dynInd(nor,vS,5),NRes),real(resampling(repmat(Bs,[ones(1,3) NRes(4)/NB(4)]),NRes(1:4)))),N));end   
        end
        if norm;rec.(rec.Plan.Types{typInt})=rec.(rec.Plan.Types{typInt})./nor;end
    else
        rec.(rec.Plan.Types{typInt})=resampling(reverseB0Field(resampling(rec.(rec.Plan.Types{typInt}),NRes),real(resampling(repmat(rec.B,[ones(1,3) NRes(4)/NB(4)]),NRes(1:4)))),N);
        rec.(rec.Plan.Types{typInt})=gather(rec.(rec.Plan.Types{typInt}));
        if isfield(rec,'E') && ~appl
            if gpu;rec.E=gpuArray(rec.E);end
            rec.E=gather(abs(resampling(reverseB0Field(resampling(rec.E,NRes),real(resampling(repmat(rec.B,[ones(1,3) NRes(4)/NB(4)]),NRes(1:4)))),N)));%KEY INSTRUCTION
        end
        if norm;nor=resampling(reverseB0Field(resampling(nor,NRes),real(resampling(repmat(rec.B,[ones(1,3) NRes(4)/NB(4)]),NRes(1:4)))),N);rec.(rec.Plan.Types{typInt})=rec.(rec.Plan.Types{typInt})./nor;end   
    end
    if magn;rec.(rec.Plan.Types{typInt})=abs(rec.(rec.Plan.Types{typInt}));end      
    rec.B=convertB0Field(rec.B,TE,ESFOV,'res','Hz');
    rec.B=real(rec.B);

    %PERMUTE BACK
    typ2Rec=rec.Dyn.Typ2Rec;
    if ~norm
        for n=typ2Rec';datTyp=rec.Plan.Types{n};
            rec.(datTyp)=permute(rec.(datTyp),[1 2 3 5 4]);
        end
    end
    rec.B=gather(rec.B);
end

% 
% % calc_echo_spacing_philips_mod.m
% %
% % Formula from http://www.spinozacentre.nl/wiki/index.php/NeuroWiki:Current_developments#B0_correction
% % modified after suggestions from Pieter Buur and Nikolaus Weiskopf; also thanks to Lawrie MacKay.
% % Please check results carefully, if possible with physicist of the Philips scanner where data were acquired from.
% % 
% % For Philips Achieva, the best equation for calculating "effective" echo spacing 
% % (corrected for SENSE factor) in milliseconds would seem to be:
% %  
% % echo spacing in msec = (1000 * wfs)/(434.215 * (etl)) (etl = EPI factor
% % + 1)
% % 
% % We have:
% % - matrix size in the phase encoding direction 
% % - EPI factor (EPI factor + 1 = echo train length (ETL))
% % - water-fat shift per pixel
% % - field strength in Tesla
% 
% % % We calculate:
% % - total bandwidth
% % - water-fat shift in Hz
% % - bandwidth per pixel in Hz
% % - echo spacing in sec
% % - echo spacing in msec
% 
% % We assume:
% % - that the SENSE factor is incorporated in Philips' EPI factor
% % - that the water-fat difference in parts-per-million (ppm) is: 3.35 (from Haacke)
% % 
% % Use for own risk. Please don't use for clinical applications. 
% % Hester Breman, 23-08-2013; latest update: 19-11-14
%  
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ variables, please adapt for each dataset (can be found in *.PAR file) ~~~~~~~~
% 
% matrixsize_phase_enc_dir    = 80; % matrix size in y direction (if that is phase encoding direction)
% epifactor 					= 43; % EPI factor
% water_fat_shift_pixel 		= 12.0; % water-fat shift per pixel
% fieldstrength_tesla 		= 3.0;  % magnetic field strength (T)
% 
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ not change below line ~~~~~~~~
% 
% water_fat_diff_ppm          = 3.35;  
% resonance_freq_mhz_tesla    = 42.576; % gyromagnetic ratio for proton (1H)
% echo_train_length           = epifactor + 1 
% water_fat_shift_hz          = fieldstrength_tesla * water_fat_diff_ppm * resonance_freq_mhz_tesla % water_fat_shift_hz 3T = 434.215 Hz
% BW_hz_pixel                 = water_fat_shift_hz / water_fat_shift_pixel
% totBW                       = BW_hz_pixel * echo_train_length
% echo_spacing_sec            = 1/totBW 
% echo_spacing_msec           = echo_spacing_sec * 1000
% 
% 
