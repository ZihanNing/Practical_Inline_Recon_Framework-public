function svr=svrResiduals(svr)

%SVRRESIDUALS   Computes the residuals of the reconstruction
%   SVR=SVRSLICEWEIGHTS(SVR)
%   SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)

svr.xx=svr.x;
svr=svrEncode(svr,0);
for v=1:svr.NV;svr.E{v}=svr.yy{v};end
svr.xx=svr.M;
svr=svrEncode(svr,0);
for v=1:svr.NV;svr.My{v}=real(svr.yy{v});end
    

