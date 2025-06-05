function [CParent,CconvParent,C,Cconv,Parent,Childs,S]=subdivideClusters(C,Cconv,X,xx,dim,Nmin,typFeat,S)

%SUBDIVIDECLUSTERS  subdivides an existing set of row/column clusters 
%   [CPARENT,CCONVPARENT,C,CCONV,PARENT,CHILDS]=SUBDIVIDECLUSTERS(C,CCONV,X,XX,GLOBALFEATURE,DIM)
%   * C is the current set of clusters in the data
%   * CCONV indicates whether convergence has been achieved for a given
%   cluster
%   * X is the denoised matrix
%   * XX are the spatial coordinates
%   * DIM is the dimension for clustering
%   * NMIN is the minimum number of points per cluster
%   * TYPFEAT is the type of feature to compute. One of the following:
%   'NormEnerg', 'Alpha'
%   * S are the singular values
%   * {CPARENT} are the parent clusters
%   * {CCONVPARENT} indicates whether convergence has been achieved for the
%   parent clusters
%   ** {C} is the current set of clusters in the data
%   ** {CCONV} indicates whether convergence has been achieved for a given
%   cluster
%   ** {PARENT} is the index of the parent for current clusters
%   ** {CHILDS} are the indexes of the childs for a given parent
%   ** {S} are the singular values
%

if isempty(gcp);parpool;end
opts = statset('UseParallel',1);

MN=size(X);
NC=length(C);
if ~iscell(S);F=computeFeature(X,typFeat,S);
else SParent=S;S=cell(1,2*NC);
end
CParent=C;
CconvParent=Cconv;
C=cell(1,2*NC);Cconv=nan(1,2*NC);Parent=nan(1,2*NC);Childs=cell(1,NC);
c=1;
for n=1:NC
    if CconvParent(n)%Already converged
        C{c}=CParent{n};Cconv(c)=1;Parent(c)=n;
        if iscell(S);S{c}=SParent{n};end
        Childs{n}=c;
        c=c+1;
    else
        if iscell(S);Fs=computeFeature(dynInd(X,CParent{n},dim),typFeat,SParent{n});else Fs=dynInd(F,CParent{n},dim);end
        if ~isempty(xx);Fs=cat(setdiff(1:2,dim),Fs,dynInd(xx,CParent{n},dim));end%We add the spatial distances
        we=ones(size(Fs,setdiff(1:2,dim))-1,1)/(30*MN(setdiff(1:2,dim)));
        if dim==1;on=on';end
        wF=cat(setdiff(1:2,dim),1,we);
        Fs=bsxfun(@times,wF,Fs);
        if dim==2;Fs=Fs';end
        Clev=kmeans(Fs,2,'Options',opts,'Replicates',8);
        if sum(Clev==1)<Nmin || sum(Clev==2)<Nmin;
            fprintf('Number of clusters: %d / %d. Limit: %d\n',sum(Clev==1),sum(Clev==2),Nmin);
            C{c}=CParent{n};Cconv(c)=1;Parent(c)=n;
            if iscell(S);S{c}=SParent{n};end
            Childs{n}=c;CconvParent(n)=1;
            c=c+1;
        else
            Childs{n}=[c c+1];
            C{c}=CParent{n}(Clev==1);Cconv(c)=0;Parent(c)=n;
            if iscell(S);S{c}=SParent{n};end
            c=c+1;
            C{c}=CParent{n}(Clev==2);Cconv(c)=0;Parent(c)=n;
            if iscell(S);S{c}=SParent{n};end
            c=c+1;
        end
    end
end
notAssign=find(isnan(Cconv));
Cconv(notAssign)=[];
C(notAssign)=[];
Parent(notAssign)=[];
if iscell(S);S(notAssign)=[];end

function F=computeFeature(F,typFeat,S)
    if strcmp(typFeat,'NormEnerg')
        F=bsxfun(@times,F,1./(sqrt(sum(abs(F).^2,dim))+eps));%Normalize
        F=sqrt(sum(abs(F).^2,setdiff(1:2,dim))/MN(setdiff(1:2,dim)));
    elseif strcmp(typFeat,'Alpha')
        F=sqrt(sum(abs(F).^2,setdiff(1:2,dim)));
    elseif strcmp(typFeat,'AlphaW')%Alpha from William Leeb paper
        if dim==2;S=S(:);end
        F=sqrt(sum(abs(bsxfun(@times,S,F)).^2,setdiff(1:2,dim)));
    else
        error('Not contemplated feature: %s\n',typFeat);
    end
end

end
