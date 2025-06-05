function [x,ESFOV,TE]=applyFrequencyStabilization(x,T,rec,gpu,NEqLines)

% APPLYFREQUENCYSTABILIZATION applies the computed frequency stabilization
% to the data
%   [X,ESFOV]=APPLYFREQUENCYSTABILIZATION(X,T,REC,{GPU}) 
%   * X is an array that needs to be stabilized, 4th dimension are dynamics 
%   and fifth dimension are components
%   * T are some stabilization shifts
%   * REC is a reconstruction structure with information about
%   stabilization
%   {GPU} determines whether to use gpu computations 
%   {NEQLINES} is the number of equivalent lines for if this is applied to
%   subsampled data
%   ** X is a stabilized array
%   ** ESFOV is the effective echo spacing for the reconstruction FOV
%

if nargin<4 || isempty(gpu);gpu=isa(x,'gpuArray');end

NE=length(rec.Par.Labels.TE);
if rec.Par.Labels.TE(end)<rec.Par.Labels.TE(1);NE=1;end%Sometimes it appears there is a second TE while there isn't

%TE
if NE==1;TE=rec.Par.Labels.TE(1);elseif NE==2;TE=diff(rec.Par.Labels.TE);else error('Distortion reversal only working for single or double echo data');end          
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
N=size(x);
if nargin<5 || isempty(NEqLines);NEqLines=N(2);end
ESFOV=ES*NPE/NEqLines;

B=convertB0Field(T{2},TE,ESFOV,'pix','rad',N(2));

if gpu;x=gpuArray(x);end
x=bsxfun(@times,exp(-1i*B),x);
%if numel(rec.x)>4e8;rec.x=gather(rec.x);T{2}=gather(T{2});end    
x=shifting(x,T);
