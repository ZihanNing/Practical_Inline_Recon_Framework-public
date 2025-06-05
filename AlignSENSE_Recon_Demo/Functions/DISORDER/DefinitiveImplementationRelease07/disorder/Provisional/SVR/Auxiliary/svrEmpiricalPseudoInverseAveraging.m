function svr=svrEmpiricalPseudoInverseAveraging(svr,filt,arr)

%SVREMPIRICALPSEUDOINVERSEAVERAGING   Performs empirical pseudo-inverse of
%the data and the mask
%   SVR=SVREMPIRICALPSEUDOINVERSEAVERAGING(SVR)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%

if nargin<2 || isempty(filt);filt=1;end
if nargin<3 || isempty(arr);arr=0;end

svr.yy=svr.y;
svr=svrDecode(svr,1,filt,arr);
xx=svr.xx;
svrWriteData(svr,sprintf('Step%02d-StacksSpatialCoordinates',svr.ParSVR.Step),xx,[],[],[],1);svr.ParSVR.Step=svr.ParSVR.Step+1;
svr.yy=svr.My;
svr=svrDecode(svr,1);
MMxx=mapToZeroOne(svr.xx,0.05);
svrWriteData(svr,sprintf('Step%02d-StacksSpatialCoordinatesMask',svr.ParSVR.Step),MMxx,[],[],[],1);svr.ParSVR.Step=svr.ParSVR.Step+1;
svr.x=mean(xx,4);
if isfield(svr,'EnViews');svr.x=svr.x*svr.NV/sum(svr.EnViews);end
svrWriteData(svr,sprintf('Step%02d-SpatialCoordinates',svr.ParSVR.Step),svr.x,[],[],[],1);svr.ParSVR.Step=svr.ParSVR.Step+1;
svr.Mx=mean(MMxx,4);
if isfield(svr,'EnViews');svr.Mx=svr.Mx*svr.NV/sum(svr.EnViews);end
svrWriteData(svr,sprintf('Step%02d-SpatialCoordinatesMask',svr.ParSVR.Step),svr.Mx,[],[],[],1);svr.ParSVR.Step=svr.ParSVR.Step+1;

if svr.ParSVR.Visualize
    svrVisualization(num2cell(xx,1:3),num2cell(MMxx,1:3),sprintf('Stacks empirical pseudoinverse step %d',svr.ParSVR.Step-4),round(mean(svr.MSS(:))));
    svrVisualization(svr.x,svr.Mx,sprintf('Empirical pseudoinverse step %d',svr.ParSVR.Step-4),round(mean(svr.MSS(:))));
end
