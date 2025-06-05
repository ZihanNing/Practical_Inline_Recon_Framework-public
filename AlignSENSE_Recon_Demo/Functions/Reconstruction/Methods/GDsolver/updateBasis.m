

function [E] = updateBasis(E, parXB, NX)
if parXB.useSH==2; E.Db = E.Dbs;end
SHorder=parXB.SHorder;

transformBasis=0;
if isfield(E,'Db') && isfield(E.Db,'transformed') %we have created basis function in past and have transformed it in LMsolver
    if E.Db.transformed==1
        transformBasis=1;
        T = E.Db.transformHist;
    end
end

%%%CREATE NEW BASIS
[E.Db.B, E.Db.Bidx] = SH_basis(NX, SHorder);
E.Db.NX = NX;

%%%TRANSFORM BASIS - so rotation basis functions taken into account
if transformBasis
    Tf = precomputeFactorsSincRigidTransform(E.kG,E.rkG,T,1,0,1,1);%Transform factors 
    E.Db = TransformBasis(E.Db,T,Tf, E.Fof, E.Fob); 
end

%%%RESET TRANFORM PARAMETERS
E.Db.transformed = 0;
E.Db.transformHist = [];

if parXB.useSH==2; E.Dbs = E.Db;E = rmfield(E,'Db'); end
