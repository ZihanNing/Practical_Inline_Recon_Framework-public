function svrWriteData(svr,suff,x,spa,res,suffSVR,infor,ma,gr)

%SVRWRITEDATA   Writes svr data
%   SVR=SVRWRITEDATA(SVR,SUFF,X)
%   SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   SUFF is the suffix of the data to be written
%   X is the data to be writen
%   SPA indicates whether to write information in material coordinates
%   RES manipulates the original resolution. Only has effect if SPA equals
%   0
%   SUFFSVR is the suffix identifying the SVR run
%   INFOR is a flag to write surrogate information to a corresponding subfolder
%   MA is a flag to indicate that we are interpolating a mask
%   GR serves to use Gibbs ringing for resampling
%

if nargin<4 || isempty(spa);spa=0;end
if nargin<5 || isempty(res);res=1;end
if nargin<6 || isempty(suffSVR);suffSVR='';end
%if nargin<6 || isempty(suffSVR);suffSVR='';end
if nargin<7 || isempty(infor);infor=0;end
if nargin<8 || isempty(ma);ma=0;end
if nargin<9 || isempty(gr);gr=0;end

gpu=isa(x,'gpuArray');
if ~infor;pathModa=svr.ParSVR.pathModal;else pathModa=svr.ParSVR.pathSurro;end
if ~exist(pathModa,'dir');mkdir(pathModa);end
suff=strcat(suffSVR,suff);
if spa==0    
    MS=svr.MSS;
    MT=svr.MTT;
    N=size(x);N(end+1:3)=1;N=N(1:3);
    Nr=round(N./(res./MS));
    if ma;x=morphFourier(x,4*ones(1,3),svr.MSS,ones(1,3),1);end
    if gr~=0
        x=filtering(x,buildFilter(2*N,'tukeyIso',ones(1,3),gpu,gr,1),1);   
        x=real(resampling(x,Nr,0,ones(1,3)));
    else%Simply conventient to avoid memory problems
        x=real(resampling(x,Nr));
    end
    if ma;x=mapToZeroOne(x,0.5);end
    if res~=1;x=max(x,0);end
    for n=1:3;MT(n,n)=MT(n,n)*res/svr.MSS(n);end
    MS(:)=res;
    fileName=fullfile(pathModa,svr.rec{1}.Names.Name);    
    writeNII(fileName,{suff},{x},{MS},{MT}); 
elseif spa==1
    if ~iscell(svr.MT)
        svr.MT=num2cell(svr.MT,1:2);
        svr.MS=num2cell(svr.MS,1:2);
    end    
    for v=1:svr.NV
        fileName=fullfile(pathModa,svr.rec{v}.Names.Name);        
        writeNII(fileName,{suff},{x{v}},{svr.MS{v}},{svr.MT{v}}); 
    end
else
    fileName=fullfile(pathModa,svr.rec{1}.Names.Name);
    writeNII(fileName,{suff},{x},{svr.MSB1},{svr.MTB1});
end
