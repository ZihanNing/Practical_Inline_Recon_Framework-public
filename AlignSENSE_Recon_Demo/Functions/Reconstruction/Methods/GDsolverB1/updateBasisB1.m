

function [E] = updateBasisB1(E, parXB, NX)

transformBasis=0;
if isfield(E,'B1m') && isfield(E.B1m,'transformed') %we have created basis function in past and have transformed it in LMsolver
    if E.B1m.transformed==1
        transformBasis=1;
        T = E.B1m.transformHist;
    end
end

%%%CREATE NEW BASIS
%[E.B1m.B, E.B1m.Bidx] = PolynomialBasis(NX, parXB.B1order );
[E.B1m.B, E.B1m.Bidx] = SH_basis(NX, parXB.B1order );
%E.B1m.B=extractDCTAtoms(NX,parXB.B1order);
E.B1m.NX = NX;

%%%TRANSFORM BASIS - so rotation basis functions taken into account
if transformBasis
    Tf = precomputeFactorsSincRigidTransform(E.kG,E.rkG,T,1,0,1,1);%Transform factors 
    E.B1m = TransformBasis(E.B1m,T,Tf, E.Fof, E.Fob); 
end

%%%RESET TRANFORM PARAMETERS
E.B1m.transformed = 0;
E.B1m.transformHist = [];
