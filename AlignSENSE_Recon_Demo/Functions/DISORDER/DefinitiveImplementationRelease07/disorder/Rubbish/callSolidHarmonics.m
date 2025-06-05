function x=callSolidHarmonics(x,C,S,R0)

% CALLSOLIDHARMONICS solid harmonics call
%   X=CALLSOLIDHARMONICS(X,ALPHA,BETA)
%   X are Cartesian coordinates
%   ALPHA are the solid harmonic parameters for each direction
%   BETA are the solid harmonic parameters for each direction
%   X are the solid harmonics
%

%R0=250/1000;%(m for 3T, see file create_displacement_table_for_philips
%R0=R0;

Gref(1)= R0/C{1}(2,2);
Gref(2)= R0/S{2}(2,2);
Gref(3)= R0/C{3}(2,1);

fprintf('\n\nCreating displacement table ...\n');


fprintf('start calculations ...\n');

N=size(x);
B=zeros(N,'like',x);
for m=1:3
    B(:,:,:,m)=solidHarmonics(x,C{m},S{m},R0,Gref(m));
end
B(:,:,:,1:2)=-B(:,:,:,1:2);
x=B;
%x=-B;
