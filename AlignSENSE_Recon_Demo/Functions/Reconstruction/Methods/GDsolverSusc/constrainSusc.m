
function [E, chiCon ] = constrainSusc(E,chi) %D as separate input so you can also use constrain on e.g. gradient - in that case don't update E.Db.cr
    
    if nargin<2; chiCon = E.Ds.chi; passToC=1; else chiCon=chi; passToC=0;end 
    
    % Filtering
    if isfield(E.Ds, 'C') && isfield(E.Ds.C, 'H') && ~isempty(E.Ds.C.H)
        fprintf('Filtering susceptibility distribution.\n');
        for i =1:size(chiCon,6); chiCon = real( dynInd(chiCon,i,6, filtering(dynInd(chiCon,i,6), E.Ds.C.H, 1)));end %Filter in cosine-domain
    end
   
    
    
end