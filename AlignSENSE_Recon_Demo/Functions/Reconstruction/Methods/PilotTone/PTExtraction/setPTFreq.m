

function fPilotTone = setPTFreq (fCentre , baseResol , BWPixel, Factor, offsetFOV, resolRO, slice, PE1)

% fCentre in Hz
% baseResol in #
% BWPixel in Hz/pixel (see console)

if nargin < 4 || isempty(Factor); Factor = .95;end %How far into oversampled FOV (0 at edge FOV and 1 at edge 2*FOV)
if nargin < 5 || isempty(offsetFOV); offsetFOV = 0;end %FOV offset from isocentre in mm
if nargin < 6 || isempty(resolRO); resolRO = inf; end %resolution in mm in RO direction

BWMm = BWPixel/resolRO;%Bandwidth in Hz/mm

setLowerPart = Factor<0;
if setLowerPart; fprintf('Pilot Tone signal set to lower part of the RO FOV.\n');end

fPilotTone = fCentre + ...%In Hz
             BWMm * offsetFOV + ...
             baseResol * BWPixel * Factor ;

fPilotTone = fPilotTone;%In Hz

end

