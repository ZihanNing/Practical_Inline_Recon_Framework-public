function x=solidHarmonics(x,C,S,R0,GRef)

% SOLIDHARMONICS solid harmonics from cartesian coordinates
%   X=SOLIDHARMONICS(x,M,N,GPU)
%   X are Cartesian coordinates
%   ALPHA are the solid harmonic parameters
%   BETA are the solid harmonic parameters
%   RO is the reference radious
%   GREF is the reference gradient
%   X are the solid harmonics
%

R = sqrt(sum(x.^2,4));
R=R+eps;
Theta=dynInd(x,3,4)./R;
Theta(Theta>1)=1;
Theta(Theta<-1)=-1;
Theta=acos(Theta);
Phi=atan2(dynInd(x,2,4)./R,dynInd(x,1,4)./R);
x=[];

N=size(R);
nmax = size(C,1)-1;

B=zeros(N,'like',R);
for n = 0:nmax
    if any(C(n+1,:)~=0) || any(S(n+1,:)~=0)
        P=real(legendre(n,cos(Theta)));   
        if n~=0;P=permute(P,[2 3 4 1]);end
   
        F=(R/R0).^n;
        for m=0:min(n,size(C,2)-1)
            if C(n+1,m+1)~=0 || S(n+1,m+1)~=0
                F2=C(n+1,m+1)*cos(m*Phi)+S(n+1,m+1)*sin(m*Phi);      
                B=B+F.*P(:,:,:,m+1).*F2;
            end
        end
    end
end
B(isnan(B))=0;
x=B*GRef; % Normalize to meters
