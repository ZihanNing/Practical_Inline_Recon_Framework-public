function [A,echoLim]=extractShots(ktraj,kRange,sub,subTree,maxSubTree,sup)

%EXTRACTSHOTS generates the shot sampling structure for a given acquisition
%   [A,ECHOLIM]=EXTRACTSHOTS(KTRAJ,KRANGE,{SUB},{SUBTREE},{MAXSUBTREE},{SUP})
%   * KTRAJ is the spectral trajectory (NEchosxNShotsx2 matrix)
%   * KRANGE is the spectral coverage (2x2 matrix: Rows: PE dimensions. Columns: min/max range)
%   * {SUB} indicates whether to divide the spectrum in two for each shot (either 1 for full spectrum or 2 for halved spectrum, defaults to 1)
%   * {SUBTREE} specifies the tree level of within shot subdivisions for each half of the spectrum for this iteration (defaults to 1)
%   * {MAXSUBTREE} defines the maximum number of subdivisions of each half of the spectrum that are going to be used by the algorithm (defaults to 1)
%   * {SUP} specifies grouping of shots for motion estimation (defaults to 1) 
%   ** A is a binary matrix with the shot structure in the fifth dimension
%   ** ECHOLIM are the echo limits in within shot subdivisions
%

if nargin<3 || isempty(sub);sub=1;end
if nargin<4 || isempty(subTree);subTree=1;end
if nargin<5 || isempty(maxSubTree);maxSubTree=1;end
if nargin<6 || isempty(sup);sup=1;end

%ERROR CONTROL
assert(~(subTree>1 && sub==1),'%s within shot subdivisions have been specified without the spectrum being halved',subTree);
assert(~(maxSubTree>1 && sub==1),'%s maximum within shot subdivisions have been specified without the spectrum being halved',maxSubTree);

%DEFINITION OF THE SCHEME
NEchos=size(ktraj,1);NShots=size(ktraj,2);NDimens=size(ktraj,3);
if mod(NShots,sup)~=0
    fprintf('Shot grouping is only implemented for divisors of the number of shots\n');
else
    NShots=NShots/sup;
    NEchos=NEchos*sup;
    ktraj=reshape(ktraj,[NEchos NShots NDimens]);
end
kShift=floor((diff(kRange,1,2)+1)/2);
kSize=diff(kRange,1,2)+1;
if NDimens==1;kSize(2,1)=1;end
if subTree==1;maxSubTree=1;end

A=single(zeros([kSize' subTree sub NShots]));
NmaxSubTree=sub*maxSubTree;
NEchosEff=NEchos/NmaxSubTree;%Approximate number of echoes per motion state for the targeted temporal subdivision
ww=1:NmaxSubTree;
echoLim=round([NEchosEff*(ww-1)+1;NEchosEff*ww]);%Echo limits for the targeted within shot subdivision

%GROUPING OF SAMPLES
for ns=1:NShots
    for ws=1:sub
        if ws==1%Right-hand side, decreasing from low to high frequencies
            contwt=subTree;
            for wt=maxSubTree:-1:1
                for ne=echoLim(1,(ws-1)*maxSubTree+wt):echoLim(2,(ws-1)*maxSubTree+wt)  
                    if NDimens==2
                        if ktraj(ne,ns,1)~=inf;A(ktraj(ne,ns,1)+kShift(1)+1,ktraj(ne,ns,2)+kShift(2)+1,contwt,ws,ns)=1;end
                    else
                        if ktraj(ne,ns,1)~=inf;A(ktraj(ne,ns,1)+kShift(1)+1,1,contwt,ws,ns)=1;end
                    end
                end
                if maxSubTree-wt+1<subTree;contwt=contwt-1;end
            end
        else%Left-hand side, increasing from low to high frequencies
            contwt=1;
            for wt=1:maxSubTree
                for ne=echoLim(1,(ws-1)*maxSubTree+wt):echoLim(2,(ws-1)*maxSubTree+wt)         
                    if NDimens==2
                        if ktraj(ne,ns,1)~=inf;A(ktraj(ne,ns,1)+kShift(1)+1,ktraj(ne,ns,2)+kShift(2)+1,contwt,ws,ns)=1;end
                    else
                        if ktraj(ne,ns,1)~=inf;A(ktraj(ne,ns,1)+kShift(1)+1,1,contwt,ws,ns)=1;end
                    end
                end
                if wt<subTree;contwt=contwt+1;end
            end
        end            
    end
end
for n=1:2;A=ifftshift(A,n);end
