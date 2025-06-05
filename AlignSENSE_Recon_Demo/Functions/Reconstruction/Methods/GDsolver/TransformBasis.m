

function Db = TransformBasis(Db, T, Tf, Fof, Fob)
%used in LMsolver to correct for SH coefficients when substracting a meanT
assert(isfield(Db,'B'),'No basis functions present in encodinf structur')

NCoefs = size(Db.B,2);
NX = Db.NX;
NXflatten = [prod(NX),1];

for idx = 1:NCoefs
   
    %%% TRANSFORM BASIS FUNCTION
    b = dynInd( Db.B, idx,2);%extract ones basis
    b = reshape(b,NX);%reshape to image size
    b = real(sincRigidTransform(b,Tf,1,Fof,Fob));%transform
    b = reshape(b, NXflatten);%flatten
    
    %%%STORE
    Db.B = dynInd( Db.B, idx,2, b);    
end

%%% FLAG
Db.transformed = 1;
if isempty(Db.transformHist)
    Db.transformHist = T ;
else
    Tall = cat(5,Db.transformHist,T);
    di = [1 1]; %both forward transforms
    Db.transformHist = compositeTransform(Tall, di);%T applied after transformHist
end

end

