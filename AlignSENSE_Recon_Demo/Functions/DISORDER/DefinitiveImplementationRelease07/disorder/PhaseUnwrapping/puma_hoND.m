function [unwph,k,iter,erglist] = puma_hoND(psi,p,varargin)
%puma_ho   ND extension of graph optimization based phase unwrapping algorithm.
%   [unwph,iter,erglist] = puma_ho(psi,p,'potential',potential,'cliques',cliques, 'qualitymaps', qualitymaps,
%   'schedule',schedule)
%   Unwrapps the observed modulo-2pi phase, casting the problem as an energy minimization via graph mincut
%   calculation. The algorithm is described in "Phase Unwrapping via Graph Cuts" submitted to IEEE IP, October, 2005.
%   Herein we make a generalization, by allowing cliques of any order (not
%   just 1st order) and N-D images / clique weights
%
%   Authors: Jose Bioucas-Dias, Gonzalo Valadao 
%             Lucilio Cordero-Grande (ND extension / clique weights)
%
%   Last change: Goncalo Valadao (19/9/2012 23h46m)
%   Last change: Lucilio Cordero-Grande (~01/10/2017 00h00m)
%
% ======================== REQUIRED INPUT PARAMETERS ======================
% Parameter             Values
% name                  and description
% =========================================================================
%
% psi                   (double) The wrapped phase image.
% p                     (double) It defines the clique potential exponent (>0).
%
% ======================== OPTIONAL INPUT PARAMETERS ======================
% Parameter             Values
% name                  and description
% =========================================================================
%
% potential             (1x1 struct array) This struct array has 2 fields:
% potential.quantized   ('yes','no') Default: 'no'.
% potential.threshold   (double) it defines a region over which the
%                        potential grows quadratically. By default is pi.
%
% c               (nxN double matrix) Each row defines the
%                       "displacement vector" corresponding to each clique.
%                       The columns correspond to the dimension of the
%                       image along which the clique operates. By default
%                       is [1 0;0 1] (first order cliques, 2D image).
%
% q           (size(psi) x n (nocliques) double array). The quality matrices
%                       may take values between 0 and 1 (value 1: discontinuity presence;
%                       value 0: discontinuity absence).
%                       There is one quality matrix per clique type. By default there is
%                       discontinuity absence. A quality value corresponding to a certain
%                       clique must be signalled in the pixel corresponding to the end of
%                       the clique displacement vector (for each pair of pixels).
%
% schedule              (double vector) This vector contains a schedule of jump sizes.
%
% verbose               ('yes', 'no') -> display the unwrapped  phase along  
%                        the iterations.
%                        Default = 'yes'
%
%
% ========================= OUTPUT PARAMETERS =============================
% Parameter             Values
% name                  and description
% =========================================================================
% unwph                 (double array) This is the unwrapped phase image.
% iter                  (double) This is the number of iterations the algorithm runs through
% erglist               (double vector) This is a sequence of the energies of the unwrapped
%                       phases along the algorithm run.
%
% =========================== EXAMPLES ====================================
%       Note: the optional arguments must be provided in pairs string+value;
%       those pairs may be entered in any order.
%
%       potential.quantized = 'no'; potential.threshold = 0.5;
%       [unwph,iter,erglist] = puma_ho(psi,2,'potential',potential)
%
%       potential.quantized = 'yes'; potential.threshold = 2;
%       c = [1 0; 0 1; -1 1];
%       [unwph,iter,erglist] = puma_ho(psi,1,'c',c,'potential',potential)
%
%       potential.quantized = 'no';
%       potential.threshold = 0.1; c = [1 1];
%       q = ones(size(psi,1),size(psi,2))
%       [unwph,iter,erglist] = puma_ho(psi,p,'potential',potential,'c',c,'q',q)

% ========================== REFERENCES ===================================
%   For reference see:
%   J. Bioucas-Dias and G. Valad�o, "Phase Unwrapping via Graph Cuts"
%   IEEE Transactions Image Processing, 2007 (to appear).
%   The algorithm here coded corresponds to a generalization for any
%   c set (not only vertical and horizontal).
%
%   J. Bioucas-Dias and J. Leit�o, "The ZpiM Algorithm for Interferometric Image Reconstruction
%   in SAR/SAS", IEEE Transactions Image Processing, vol. 20, no. Y, 2001.
%
%   The technique here employed is also based on the article:
%   V. kolmogorov and R. Zabih, "What Energy Functions can be Minimized via Graph Cuts?",
%   European Conference on Computer Vision, May 2002.
% =========================================================================
%
%
%
% Modification: 
%     
%    1 - Fix a bug in the way discontinties were deal with 
%        (Gon�alo Valad�o, Sep.,2012)
%
%    2 - Change in the potential default parameters:
%        potential.quantized = 'no';
%        potential.threshold = pi;
%        (J Bioucas-Dias, Sep.,2012)
%
%    3 - Introdution of the verbose input parameter:
%        verbose = 'yes' -> display (iter, enerrg_actual, jump_size)
%                           and display the unwrapped  phase along
%                           the iterations.
%        Default = 'yes'
%        
%    
%   
%
% ------------------------------------------------------------------
% Author: Gon�alo Valad�o & Jose Bioucas-Dias, 2007
%
%

%
%% -------------------------------------------------------------------------
%
% Copyright (July, 2007):        Gon�alo Valad�o (gvaladao@lx.it.pt) 
%                                Jos� Bioucas-Dias (bioucas@lx.it.pt)
%
% PUMA_HO is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ---------------------------------------------------------------------

gpu=isa(psi,'gpuArray');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Default values
pot.quant='no';
pot.thres=pi;
ND=numDims(psi);
c=eye(ND,ND);
w=ones(ND,1);
NC=size(c);
N=size(psi);
q=single(zeros(N));k=q;
v=ones(1,ND+1);v(end)=NC(1);
q=repmat(q,v);
qual=0;
init=0;
weig=0;
schedule=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Error out if there are not at least the two required input arguments
assert(nargin-length(varargin)==2,'Wrong number of required parameters');
assert(rem(length(varargin),2)==0,'Optional parameters should always go by pairs');

% Read the optional parameters
if ~isempty(varargin)
    for l=1:2:(length(varargin)-1)        
        switch varargin{l}% change the value of parameter
            case 'pot';pot = varargin{l+1};% potential definition               
            case 'c';c = varargin{l+1};% cliques to consider
            case 'w';w = varargin{l+1};weig=1;% weights of cliques
            case 'q';q=varargin{l+1};qual=1;% the quality maps
            case 'k0';k=varargin{l+1};init=1;% the quality maps
            case 'schedule';schedule=varargin{l+1};% jump size schedule
            otherwise;error(['Unrecognized parameter: ''' varargin{n} '''']);% Something wrong with the parameter string
        end;
    end
end;
NC=size(c);
assert(ND==NC(2),'Dimensions of the input image (%d) have to match the number of columns of cliques (%d)',ND,NC);
v=ones(1,ND+1);v(end)=NC(1);
assert(~(qual==1 && size(q,ND+1)~=NC(1)),'q must be a ND+1 matrix whos ND+1 size (%d) is equal to no. cliques (%d). Each hyperplane on qualitymaps corresponds to a clique',size(q,ND+1),NC(1));
assert(~(weig==1 && size(w,1)~=NC(1)),'Weight size (%d) should match clique size (%d)',size(w,1),NC(1));
assert(~(init==1 && any(size(k)~=size(psi))),'k0 size (%s) must match psi size (%s)',sprintf(' %d',size(k)),sprintf(' %d',size(psi)));

k_aux=k;
iter=0;
erg_list=[];

if qual==0;q=single(zeros(N));q=repmat(q,v);end

if weig==0;w=ones(NC(1),1);end
perm=1:ND+1;perm(1)=ND+1;perm(ND+1)=1;
w = bsxfun(@times,permute(w,perm),1-q);

maxCl = max(abs(c(:)));
ROI=single(ones(N));
if gpu;ROI=gpuArray(ROI);end
ROI=padarray(ROI,(maxCl+1)*ones(1,ND));

ext_w=padarray(w,(maxCl+1)*[ones(1,ND) 0]);
permA=1:ND+1;permA=circshift(permA,[0 1]);
ext_w=permute(ext_w,permA);
NN=size(ext_w);NN(end+1:3)=1;

ext_psi=padarray(psi,(maxCl+1)*ones(1,ND));

vr=cell(1,ND);for m=1:ND;vr{m}=maxCl+2:N(m)+1+maxCl;end

for jump_size=schedule
    impr=1;
    erg_pre=energy_hoND(k,ext_psi,ROI,p,c,ext_w,pot);
    %fprintf('Energy: %.3f\n',erg_previous);
    while impr
        iter = iter + 1;
        erg_list=horzcat(erg_list,erg_pre);
        remain=[];
        ext_k=padarray(k,(maxCl+1)*ones(1,ND));
                
        sourcefinal=zeros(N);sinkfinal=zeros(N);                        
        for t=1:NC(1) 
            ext_wCur=reshape(dynInd(ext_w,t,1),NN(2:end));
            ROIPos=circshift(ROI,-c(t,:)).*ROI;ROINeg=circshift(ROI,c(t,:)).*ROI;            
      
            dk= ext_k-circshift(ext_k,c(t,:));          
            dpsi = circshift(ext_psi,c(t,:))-ext_psi;
                     
            a=(2*pi*dk-dpsi).*ROINeg;         
            A=clique_energy_ho(abs(a),p,pot).*ROINeg.*ext_wCur;
            D=A;
            C=clique_energy_ho(abs(2*pi*jump_size+a),p,pot).*ROINeg.*ext_wCur;
            B=clique_energy_ho(abs(-2*pi*jump_size+a),p,pot).*ROINeg.*ext_wCur;            
            
            source = circshift((C-A).*((C-A)>0),-c(t,:)).*ROIPos+((D-C).*((D-C)>0)).*ROINeg;
            sink = circshift((A-C).*((A-C)>0),-c(t,:)).*ROIPos+((C-D).*((C-D)>0)).*ROINeg;

            aux1=B+C-A-D;
            [source,sink,aux1,ROIPos,ROINeg]=parUnaFun({source,sink,aux1,ROIPos,ROINeg},@dynInd,vr,1:ND);
            [start,end]=parUnaFun({ROIPos~=0,ROINeg~=0},@find);                     
            zArr=zeros([length(endd) 1]);            
            aux2=horzcat(start,endd,double(aux1(endd).*(aux1(endd)>0)),zArr);%To use the double in aux1 solved a tricky bug as otherwise indexing was wrong because of conversion to float precision!!!!
            aux2(aux2(:,3)==0,:)=[];
            remain=vertcat(remain,aux2);
            
            sourcefinal=sourcefinal+source;
            sinkfinal=sinkfinal+sink;
            
        end
        sourcesink=horzcat((1:prod(N))',double(sourcefinal(:)),double(sinkfinal(:)));%To use the double in sourcefinal/sinkfinal solved a tricky bug as otherwise indexing was wrong because of conversion to float precision!!!!      
        
        typ='Dou';%Necessary for big graphs (perhaps also 'Int' can be used, to be checked).
        %typ='Flo';
        %typ='Int';
        if strcmp(typ,'Int')        
            NL=65536;        
            [remain,sourcesink]=quantizeEnergy(remain,sourcesink,NL);
        end
        
        % KAPPA RELABELING
        [~,cutside] = mincut(gather(sourcesink),gather(remain),typ);%the flow is not retrieved
        k_aux(cutside(:,1)+1)=k(cutside(:,1)+1)+(1-cutside(:,2))*jump_size;       
     
        erg_act=energy_hoND(k_aux,ext_psi,ROI,p,c,ext_w,pot);
        %fprintf('Energy: %.3f\n',erg_actual);
        if erg_act<erg_pre
            erg_pre=erg_act;
            k=k_aux;
        else
            impr=0;
            unwph=2*pi*k+psi;
        end
    end
end

function [remain,sourcesink]=quantizeEnergy(remain,sourcesink,NL)
    maxPot=0;
    for n=3:4
        maxAux=max(remain(:,n));
        maxPot=max(maxPot,maxAux);
    end
    for n=2:3
        maxAux=max(sourcesink(:,n));
        maxPot=max(maxPot,maxAux);
    end

    remain(:,3:4)=round((NL-1)*(remain(:,3:4)/maxPot));
    sourcesink(:,2:3)=round((NL-1)*(sourcesink(:,2:3)/maxPot));       
end

end

