function svr=svrGenerateSoftMask(svr,mode)

%SVRGENERATESOFTMASK   Generates a soft from a hard mask
%   SVR=SVRGENERATESOFTMASK(SVR)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%

if nargin<2 || isempty(mode);mode=1;end

if mode==1
    dist=4;%THIS IS A HARDCODED PARAMETER
    svr.My=svr.El;
    for v=1:svr.NV;svr.My{v}=morphFourier(svr.El{v},dist*ones(1,3),svr.MSS,ones(1,3),1);end
elseif mode==2
    dist=4;%THIS IS A HARDCODED PARAMETER
    svr.M=svr.Mx;
    svr.M=morphFourier(svr.M,dist*ones(1,3),svr.MSS,ones(1,3),1);
elseif mode==3  
    dist=-4;%THIS IS A HARDCODED PARAMETER
    svr.M=svr.Mx;
    svr.M=morphFourier(svr.M,dist*ones(1,3),svr.MSS,zeros(1,3));
    svr.M=single(svr.M>0.5);
    %svr.M=morphFourier(morphFourier(svr.M,-dist*ones(1,3),svr.MSS,zeros(1,3)),dist*ones(1,3),svr.MSS,zeros(1,3),1);%First erosion then distance function for soft masking including only voxels inside
elseif mode==4
    dist=4;%THIS IS A HARDCODED PARAMETER
    svr.M=svr.Elx;
    svr.M=morphFourier(svr.M,dist*ones(1,3),svr.MSS,ones(1,3),1);
end