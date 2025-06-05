

function [D] = unwrapLinear ( D ,x, explicit, TE, theta) 

%UNWRAPLINEAR   Tries to invert the sign from estimated linear maps of the pose-dependent fields w.r.t motion parameters first unwrapping the
%   individual phase maps and then fitting the linear model again.
%   [X,T]=UNWRAPLINEAR(D,TE,{THETA})
%   * D is the linaer term (in Hz/degree)
%   * TE the echo time (in seconds)
%   * {THETA} the angles (in degree) to use for 
%   ** D The 'unwrapped' linear term
%
N  = size(D);

if nargin<2 || isempty(x); x = ones(N(1:3), 'like', D); end%in seconds
if nargin<3 || isempty(explicit); explicit = 0; end%in seconds
if nargin<4 || isempty(TE); TE = +5e-3; end%in seconds
if nargin<5 || isempty(theta); theta = convertRotation([-15:4:15 ]','deg','rad');end%in radians

%%% Check if on pathname
Folder='/home/ybr19/Software/DISORDER/DefinitiveImplementationRelease07/disorder/Provisional/DynamicDistortion';%This contains old version of unwrapper and will throw error!
pathCell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
  onPath = any(strcmpi(Folder, pathCell));
else
  onPath = any(strcmp(Folder, pathCell));
end
if onPath; rmpath(genpath(Folder));end

%%% Magnitude for the phase unwrapping
if isequal( x, ones(N));warning('unwrapLinear.m:: No magnitude image provided in D. Might result in suboptimal unwrapping performance.');end

if explicit == 0 %Go via motion to phase and unwrap in time
    fprintf('Unwrapping linear term implicitely.\n');

    theta = reshape(theta, [numel(theta), 1]);
    theta = permute(theta, [2:4 1]);

    %Phase = 2*pi* TE * bsxfun(@times, theta, D);
    PhaseWr = exp( +1i * 2*pi* TE * bsxfun(@times, theta, D));

    NS = length(theta);

    for i = 1:NS
        z = CNCGUnwrapping(x.*dynInd(PhaseWr,i,4) , ones(1,numDims(dynInd(PhaseWr,i,4)) ),'MagnitudeGradient4','LSIt'); %! Ideally has a magnitude as well!
        z = z/2/pi/TE; % in Hz

        if i==1; Z = z(:).'; else Z = cat(1, Z, z(:).' );end %Stack like this to make LS easy later
    end
    z=[]; PhaseWr =[];

    % plotND([],angle(PhaseWr),[-pi pi],[],0,{0:2,1},0.3, [], [], 301);
    % plotND([],PhaseUnwr,[],[],0,{0:2,1},0.3, [], [], 401);

    % LS interpolation of the linear D term ( Z = T D)
    D = theta(:)\Z; Z=[];
    D = reshape(D, N); %in Hz/rad
    
elseif explicit ==1
    fprintf('Unwrapping linear term explicitely.\n');
    %Scale to phase domain
    rangeIn = [-5 5]; rangeOut = [-pi pi];
    D = rescaleND(D, rangeOut, rangeIn);
    
    x = ones(N(1:3), 'like', D);
    D = CNCGUnwrapping(bsxfun(@times, x,exp(+1i .* D)), ones(1,numDims(D)),'MagnitudeGradient4','LSIt');
    %Scale back
    D = rescaleND(D, rangeIn, rangeOut);
end

