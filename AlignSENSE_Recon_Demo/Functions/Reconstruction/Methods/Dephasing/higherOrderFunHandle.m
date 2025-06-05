
function [f, fDeriv] = higherOrderFunHandle(n, yawIdx, rollIdx , pitchIdx)

%HIGHERORDERFUNCTIONHANDLE  returns a function handle to extract the relevant motion
%parameters from the motion parameters for the B0 Taylor model.
%   [F,FDERIV]=HIGHERORDERFUNCTIONHANDLE(N,{YAWIDX},{ROLLIDX},{PITCHIDX})
%   * N is the order of the Taylor model
%   * {YAWIDX} is the index where the yaw rotation is stored in the 6 motion parameters
%   * {ROLLIDX} is the index where the roll rotation is stored in the 6 motion parameters
%   * {PITCHIDX} is the index where the pitch rotation is stored in the 6 motion parameters
%   ** F is the function handle to extract the relevant motion parameters
%   ** FDERIV is the function handle to extract the relevant motion parameters whem computing the derivatives - not used currently
%

if nargin <2 || isempty(yawIdx); yawIdx = 4;end
if nargin <3 || isempty(rollIdx); rollIdx = 5;end
if nargin <4 || isempty(pitchIdx); pitchIdx  = 6;end

if n==0
    f = @(T) cat(6,0);
    fDeriv = @(T) cat(7, cat(6,0),cat(6,0),cat(6,0));
    
elseif n ==1 %1st order Taylor expansion (Linear model)
    
    % Rotation terms that multiply with the Taylor maps
    f = @(T) cat(6, dynInd(T, rollIdx, 6), ...
                    dynInd(T, pitchIdx, 6));
    
    % Derivatives used in LMsolver to add contribution to Jacobian
    fDeriv = @(T) cat(7, cat(6,0,0), ... 
                         cat(6,1,0), ...
                         cat(6,0,1) ) ;
    
elseif n==2 %2nd order Taylor expansion
    % Rotation terms that multiply with the Taylor maps
    f = @(T) cat(6, dynInd(T, rollIdx, 6) - dynInd(T, pitchIdx, 6) .* dynInd(T, yawIdx, 6), ...
                    dynInd(T, pitchIdx, 6) + dynInd(T, rollIdx, 6) .* dynInd(T, yawIdx, 6), ...
                    dynInd(T, pitchIdx, 6) .* dynInd(T, rollIdx, 6), ...
                    dynInd(T, rollIdx, 6) .^2, ...
                    dynInd(T, pitchIdx, 6) .^2 );
    % Derivatives used in LMsolver to add contribution to Jacobian            
    fDeriv = @(T) cat(7, cat(6,-dynInd(T, pitchIdx, 6), dynInd(T, rollIdx, 6),0*dynInd(T, rollIdx, 6),0*dynInd(T, rollIdx, 6),0*dynInd(T, rollIdx, 6)), ... % Yaw derivative
                         cat(6, ones(size(dynInd(T, rollIdx, 6))), dynInd(T, yawIdx, 6),dynInd(T, pitchIdx, 6), 2* dynInd(T, rollIdx, 6),0*dynInd(T, rollIdx, 6)), ...% Roll derivative
                         cat(6,-dynInd(T, yawIdx, 6),ones(size(dynInd(T, rollIdx, 6))),dynInd(T, rollIdx, 6),0*dynInd(T, rollIdx, 6), 2* dynInd(T, pitchIdx, 6)) ) ; % Pitch derivative
else
    error('higherOrderFunHandle:: Order %d not implemented', n);
end

