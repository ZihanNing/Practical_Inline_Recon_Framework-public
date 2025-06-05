
function E = TransformCoef(E, Tf)
%used in LMsolver to correct for SH coefficients when substracting a meanT

NStates = size(E.Db.cr,5);
cr_new = zeros(size(E.Db.cr),'like', E.Db.cr);
NX = [ size(E.Fof{1},1) , size(E.Fof{2},1) , size(E.Fof{3},1)];

tic
for s = 1:NStates
    
    c_temp = dynInd(E.Db.cr,s,5);
    [~,f] = dephaseBasis(E.Db.B, c_temp, NX);
    f = real(sincRigidTransform(f,Tf,1,E.Fof,E.Fob));
    
    c_temp_new = E.Db.B\reshape(f,[numel(f),1]);%this is the transformed one, not called fnew to save memory
    c_temp_new = permute(c_temp_new, [2 3 4 5 6 1]);
    
    cr_new = dynInd(cr_new, s, 5, c_temp_new);

end
toc

%%% ASSIGN NEW COEFFICIENTS
E.Db.cr = cr_new;