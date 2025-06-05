
function [E, Dcon] = constrainTaylor(E,D)

%CONSTRAINTAYLOR  Constrains the Taylor maps of the B0 model. 
%   [E,DCON]=CONSTRAINTAYLOR(E,{D})
%   * E is the encoding structure, which can contain the Taylor maps in E.Dl.D that will be used for constraining.
%   * {D} are Taylor maps to be contrained. If provided, the maps in E.Dl will not be used
%   ** E the encoding structure, with optionally the original maps in E.Dl.D
%   ** DCON contrained Taylor maps
%    

if nargin<2; Dcon = E.Dl.D; passToC=1; else Dcon=D; passToC=0;end 

% Filtering
if isfield(E.Dl, 'C') && isfield(E.Dl.C, 'H') && ~isempty(E.Dl.C.H)
    fprintf('Filtering Taylor terms.\n');
    for i =1:size(Dcon,6); Dcon = real( dynInd(Dcon,i,6, filtering(dynInd(Dcon,i,6), E.Dl.C.H, 1)));end %Filter in cosine-domain
end
   
% % Take away lower order basis functions
% if isfield(E.Dl,'C') && isfield(E.Dl.C,'deSH') && E.Dl.C.deSH==1 && isfield(E,'Db') 
%     fprintf('Taking out SH terms.\n');
%     for i=1:size(Dcon,6) %Different Taylor terms
%         cSH = E.Db.B \ reshape( dynInd( Dcon, i,6), [size(E.Db.B,1) 1]); % Least squares interpolate with basis functions
%         cSH = permute(cSH,[2 3 4 5 6 1]);%Permute coefficients in 6th dimension
%         [~,p] = dephaseBasis( E.Db.B, cSH, dynInd(size(Dcon),1:3,2),E.Db.TE); % Remove lower order term from linear field
%         Dcon = dynInd(Dcon,i,6, bsxfun(@minus, dynInd(Dcon, i,6), p/2/pi/E.Db.TE ));
% 
%         if passToC 
%             E.Db.cr = bsxfun(@plus, E.Db.cr, bsxfun(@times, cSH, dynInd( E.Dl.f(E.Tr) , i, 6)) ); % add coefficients 
%         end
%     end
% end

end