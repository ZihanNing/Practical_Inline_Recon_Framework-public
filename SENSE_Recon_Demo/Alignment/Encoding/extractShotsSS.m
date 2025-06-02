function [ktraj,NSweeps,NSamples]=extractShotsSS(ktraj)

%EXTRACTSHOTSSS generates the shot sampling structure for a steady-state 
%sequence
%   [KTRAJ,SUPSHOTS,NSWEEPS,NSAMPLES]=EXTRACTSHOTS(KTRAJ,SUPSHOTSFACT,RESPYR)
%   * KTRAJ is the spectral trajectory (NEchosxNShotsx2 matrix)
%   ** KTRAJ is the reshaped spectral trajectory (NEchosxNShotsx2 matrix)
%   ** NSWEEPS is the number of motion states at this resolution level
%   ** NSAMPLES is the number of samples per motion state at this resolution
%   level
%   ** SUPSHOTS is the grouping factor of motion states at this resolution
%   level (not used currently)
%

%DETECT THE NUMBER OF SWEEPS
supShots=supShotsFact;
ktrajs=abs(fct(ktraj));
N=size(ktrajs);
[~,iAs]=max(dynInd(ktrajs,1:floor(N(1)/2),1),[],1);
iM=min(iAs,[],3);  
NProf=size(ktraj,1);
NSamples=NProf/(iM-1);
NSweeps=round(NProf/NSamples);
fprintf('Number of detected steady state sweeps: %d\n',NSweeps);
NSweeps=max(round(NSweeps/supShots),1);
NSamples=NProf/NSweeps;
NSweeps=round(NProf/NSamples);
NRepeats=size(ktraj,4);

%GROUP THE SAMPLES
cont=1;
ktrajaux=inf*single(ones([ceil(NSamples) NSweeps 2 NRepeats]));
for n=1:NSweeps
    nelShot=round(n*NSamples)-round((n-1)*NSamples);
    ktrajaux(1:nelShot,n,:)=ktraj(cont:cont+nelShot-1,1,:);
    cont=cont+nelShot;
end
ktraj=ktrajaux;
