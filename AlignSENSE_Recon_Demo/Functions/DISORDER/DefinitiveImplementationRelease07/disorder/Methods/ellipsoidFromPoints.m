function [par,v,chi2]=ellipsoidFromPoints(X)

%ELLIPSOIDFROMPOINTS fits an ellipsoid to a set of xyz data points. The code is
%based on a modification of the following implementation:
%https://uk.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit:
%Copyright (c) 2015, Yury Petrov
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the 
%   distribution
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
%IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
%THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
%PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
%CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
%EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
%PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
%PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
%LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
%NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%   [CENTER,RADII,EVECS,V,CHI2]=ELLIPSOIDFROMPOINTS(X)
%   * PAR are the 9 ellipsoid parameters. The parameters are specified as 
%   c1,c2,c3,r1,r2,r3,th1,th2,th3 (center,radious,rotation, respectively 
%   in pixels and degrees). They correspond to the parameters used to
%   generate the ellipsoid in the BUILDGEOMETRY method
%   * V are the 10 parameters describing the ellipsoid /conic
%   algebraically: Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 
%   2Iz + J = 0
%   * CHI2 is the residual sum of squared errors (chi^2), this chi2 is in 
%   the coordinate frame in which the ellipsoid is a unit sphere.
%

assert(size(X,2)==3,'Input data must have three columns and it has %d',size(X,2));
    
x=X(:,1);y=X(:,2);z=X(:,3);


% fit ellipsoid in the form Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx +
% 2Hy + 2Iz + J = 0 and A + B + C = 3 constraint removing one extra
% parameter
D=[x.*x+y.*y-2*z.*z x.*x+z.*z-2*y.*y 2*x.*y 2*x.*z 2*y.*z 2*x 2*y 2*z 1+0*x];

% solve the normal system of equations
d2=x.*x+y.*y+z.*z; % the RHS of the llsq problem (y's)
u=(D'*D)\(D'*d2);% solution to the normal equations

% find the residual sum of errors
% chi2=sum((1-(D*u)./d2).^2); % this chi2 is in the coordinate frame in which the ellipsoid is a unit sphere.

% find the ellipsoid parameters
% convert back to the conventional algebraic form
v=[u(1)+u(2)-1;u(1)-2*u(2)-1;u(2)-2*u(1)-1;u(3:9)];
%v=v';

% form the algebraic form of the ellipsoid
A=[v(1) v(4) v(5) v(7);
   v(4) v(2) v(6) v(8);
   v(5) v(6) v(3) v(9);
   v(7) v(8) v(9) v(10)];
% find the center of the ellipsoid
center=-A(1:3,1:3)\v(7:9);
% form the corresponding translation matrix
T=eye(4,'like',X);
T(4,1:3)=center';
% translate to the center
R=T*A*T';
% solve the eigenproblem
[evecs,evals]=eig(R(1:3,1:3)/-R(4,4));

radii=sqrt(1./diag(abs(evals)));
sgns=sign(diag(evals));
radii=radii.*sgns;
if det(evecs)<-0.99;evecs=-evecs;end
[c,r,U]=parUnaFun({center,radii,evecs},@gather);

a=SpinCalc('DCMtoEA123',U',1e-3,0);
par=[c' r' a];

% calculate difference of the fitted points from the actual data normalized by the conic radii
d=bsxfun(@minus,x,center');%shift data to origin
d=d*evecs; % rotate to cardinal axes of the conic;
d=bsxfun(@rdivide,d,radii');%normalize to the conic radii
chi2=sum(abs(1-sum(bsxfun(@times,d.^2,sgns'),2)));

% if big enough, normalize to the more conventional form with constant term = -1
if abs(v(end))>1e-6;v=-v/v(end);else v=-sign(v(end))*v;end

[v,chi2]=parUnaFun({v,chi2},@gather);

