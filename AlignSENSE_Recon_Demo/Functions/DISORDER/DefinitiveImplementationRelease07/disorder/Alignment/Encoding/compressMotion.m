function [Tou,F,NE,cou,couIdx]=compressMotion(Tin,F,res,parXT,cin,disableGrouping)

%COMPRESSMOTION   Bins the motion states using a Haar wavelet 
%   decomposition for quicker reconstructions
%   [TOU,F,NE]=COMPRESSMOTION(TIN,F,RES,PARXT)
%   * TIN are the input transforms for the motion states
%   * F are the Fourier encoding matrices with the different motion
%   states arranged as cell column elements
%   * RES is the normalized, isotropic equivalent data resolution
%   * PARXT are the parameters for aligned reconstruction
%   * CIN are the input coefficients for the motion states
%   * DISABLEGROUPING disables the grouping
%   ** TOU are the transforms for the compressed motion states
%   ** F are the Fourier encoding matrices with the different compressed
%   motion states arranged as cell column elements
%   ** NE are the starting indexes of the Fourier encoding matrices for the
%   different compressed motion states
%

if nargin < 5 || isempty(cin); cin = zeros(size(Tin)); end % Added by YB - coefficients of basis
if nargin < 6 || isempty(disableGrouping); disableGrouping = 0; end 
compress_thres = 10; % ZN

%INITIALIZATION
NT=size(Tin);
Nc = size(cin);Nc(end+1:6)=1;
NSt=length(F{1});%Number of states
Nre=NSt-NT(5);%Number of states not associated with motion but with boundary conditions of the spectrum (samples outside the elliptic shutter). Corresponds to the number of repeats

Tou=Tin;
cou = cin;% Basis coefficients
couIdx = ones(Nc(1:5) ,'like',cin);
%COMPRESSION
perm=[6 5 1 2 3 4];
if NT(5)>1
    %TRANSFORM COMPRESSION
    Test=permute(Tin,perm);
    NL=floor(log2(prod(NT(1:5))));%Number of levels of wavelet decomposition %YB: log2 so that last level is on motion state resolution
    C=cell(1,NT(6));L=cell(1,NT(6));
    for m=1:NT(6)%Number of parameters (6)     %YB: Doing the Haar decomposition for every motion parameter indicidually 
        if m<4;th=parXT.traLimX(1)*res;else th=convertRotation(parXT.traLimX(2)*res,'deg','rad');end%Binning factor depends on resolution
        [C{m},L{m}]=wavedec(Test(m,:),NL,'haar');
        C{m}(abs(C{m})<th)=0;
        Test(m,:) = waverec(C{m},L{m},'haar');
    end
    Test=ipermute(Test,perm);    

    %SAMPLE GROUPING
    n=1;
    NB=0; %YB: Number of bins
    indGroup=single(ones([1 NT(5)]));
    while n<=NT(5)
        Tcur=dynInd(Test,n,5);
        m=1;
        while m+n<=NT(5)
            Tnex=dynInd(Test,n+m,5);
            if any(abs(Tnex(:)-Tcur(:))>(-1)^disableGrouping * compress_thres * 1e-6);break;end %YB: when setting this negative, you will not group anything
            m=m+1;
        end
        indGroup(n:n+m-1)=m;
        n=m+n;
        NB=NB+1;
    end
    
    %MOTION GROUPING    
    Tou=single(zeros([NT(1:4) NB NT(6)]));%YB: Size of T is also going to change, not only grouping indices!
    cou = zeros([Nc(1:4) NB Nc(6)],'like',cin);
    %Inner samples
    n=1;nb=1;     
    while n<=NT(5)
        nG=n:n+indGroup(n)-1;
        if ~disableGrouping
            Tou=dynInd(Tou,nb,5,mean(dynInd(Test,nG,5),5));
        else
            Tou=dynInd(Tou,nb,5,mean(dynInd(Tin,nG,5),5));%YB: added since in simulations, I want the T not to be altered at all, and Test could change a little bit (negligible in-vivo)
        end
        cou = dynInd(cou,nb,5,mean(dynInd(cin,nG,5),5));
        couIdx = dynInd(couIdx,nG,5,nb); % nb*ones(size(nG))
        nGr=nG(nG>nb);        
        for f=1:2;F{f}{nb}=cat(1,F{f}{nG});F{f}(nGr)={[]};end
        n=n+indGroup(n);
        nb=nb+1;
    end    
    %Outer samples and clean non used motion states / harmonics
    nbG=nb:nb+Nre-1;
    nG=n:n+Nre-1;nGr=nG(nG>nb+Nre-1);
    %------YB: dealing with empty shots for motion as well ----START
    idxEmptyShots = cellfun(@isempty,F{1});
    idxEmptyShotsTest = cellfun(@isempty,F{2});assert(isequal(idxEmptyShotsTest,idxEmptyShots),'error --> look at difference F{1} and F{2}');
    idxShots = ~idxEmptyShots;
    idxShots = idxShots(1:NT(5));
    Tou = dynInd(Tou,idxShots,5);
    cou = dynInd(cou,idxShots,5);
    couIdx = dynInd(couIdx,idxShots,5);
    %------YB: dealing with empty shots for motion as well ----START
    for f=1:2
        F{f}(nbG)=F{f}(nG);F{f}(nGr)={[]};
        F{f}=F{f}(~cellfun(@isempty,F{f}));
    end
end

%INDEX GROUPING
NE=single(zeros(1,length(F{1})+1));
for n=1:length(F{1});NE(n+1)=size(F{1}{n},1);end
NE=cumsum(NE);
