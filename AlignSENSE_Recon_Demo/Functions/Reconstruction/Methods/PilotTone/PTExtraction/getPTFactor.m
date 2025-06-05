

function Factor = getPTFactor(fPilotTone, fCentre , baseResol , BWPixel, offsetFOV, resolRO)

% fPilotTone in MHz
% fCentre in MHz
% baseResol in #
% BWPixel in Hz/pixel (see console)

if nargin < 5 || isempty(offsetFOV); offsetFOV = 0;end %FOV offset from isocentre in mm
if nargin < 6 || isempty(resolRO); resolRO = inf; end %resolution in mm in RO direction

BWMm = BWPixel/resolRO;%Bandwidth in Hz/milimeter  

Factor = (1e6*fPilotTone  - ...
         1e6*fCentre - ...
         BWMm * offsetFOV)  / ... 
         (baseResol * BWPixel );
     
setLowerPart = Factor<0;
if setLowerPart; fprintf('Pilot Tone signal was set to lower part of the RO FOV.\n');end

end
