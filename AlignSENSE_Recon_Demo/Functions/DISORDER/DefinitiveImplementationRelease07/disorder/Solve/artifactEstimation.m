function [E,EH]=artifactEstimation(E,EH,y,sigma,c,x)

%ARTIFACTESTIMATION   Estimates the artifact matrix in [1] JM Peeters and M
%Fuderer, "SENSE with improved tolerance to inaccuracies in coil
%sensitivity maps," Magn Reson Med, 69:1665-1669, 2013
%   [E,EH]=ARTIFACTESTIMATION(E,EH,{Y},{SIGMA},{C},{X})
%   * E is an encoding structure
%   * EH is a decoding structure
%   * {Y} is the measured data
%   * {SIGMA} serves to weight the level of artifacts, it defaults to 0.5
%   * {C} is the chi2-factor
%   * {X} is the reconstructed data
%   ** E is an updated encoding structure
%   ** EH is an updated decoding structure
%

if nargin<4 || isempty(sigma);sigma=0.5;end
if nargin<5;c=[];end
if nargin<6;x=[];end

if ~isfield(E,'Sf');return;end%There must be some sensitivities

if isempty(c) && ~isempty(x) && ~isempty(y);c=solveC(x,y,E,EH);end
S=E.Sf;
if ~isempty(c);E.Sf=sqrt(max(c-median(c(:)),0)).*abs(E.Sf);else E.Sf=abs(E.Sf).^2;end%Use of chi2 maps to estimate errors

NS=size(E.Sf);NS(end+1:3)=1;NS=NS(1:3);
x=ones(NS,'like',E.Sf);
E.Zf=[];

%THE PROBLEM IS HERE!!
%EH=rmfield(EH,'Bb');%SENSE 3+shift
%E=rmfield(E,'Bf');%SENSE 3+shift

x=encode(x,E);

Sf=E.Sf;
if isfield(E,'Sf');E=rmfield(E,'Sf');end
if isfield(EH,'Mb');EH=rmfield(EH,'Mb');end


x=abs(decode(x,EH,E));
rho2=normm(y)./normm(S);

E.Zf=sigma^2*rho2*max(x-Sf,0);
E.Zf=1./(1+E.Zf);
E.Sf=S;

