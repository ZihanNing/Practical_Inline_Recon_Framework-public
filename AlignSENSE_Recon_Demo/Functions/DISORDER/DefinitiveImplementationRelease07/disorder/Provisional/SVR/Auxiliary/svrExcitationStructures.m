function svr=svrExcitationStructures(svr)

%SVREXCITATIONSTRUCTURES   Sets up the excitation structures for SVR: 
%Slice orders, number of packages, number of slices per package and 
%indexes of slices in package
%   SVR=SVREXCITATIONSTRUCTURES(SVR)
%   * SVR is a svr structure containing different views (svr.y), 
%   spacings (svr.MS), orientations (svr.MT) and reconstruction parameters 
%   (svr.rec)
%   ** SVR is a svr structure containing different views (svr.y), spacings 
%   (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%

svr.NV=length(svr.y);
svr.MPack=cell(svr.NV);
for v=1:svr.NV    
    typ=svr.rec{v}.Plan.Typ==1;
    coi=svr.rec{v}.Plan.Assign{4}==0;
    kyv=svr.rec{v}.Plan.Assign{2}==0;
    ave=svr.rec{v}.Plan.Assign{12}==0;
    dyn=svr.rec{v}.Plan.Assign{5}==svr.dynamic{v}-1;
    
    svr.slOrd{v}=svr.rec{v}.Plan.Assign{8}(typ & coi & kyv & ave & dyn);    
    if svr.rec{v}.Par.Mine.AdHocArray(1)==101;svr.MBFactor(v)=svr.rec{v}.Par.Mine.AdHocArray(4);else svr.MBFactor(v)=1;end
    svr.NSlices{v}=size(svr.y{v},3);    
    svr.P(v)=numel(find(diff(svr.slOrd{v})<0))+1;
    svr.MPack{v}=real(zeros([1 1 svr.NSlices{v} 1 svr.P(v)],'like',svr.y{v}));%Mask of packs       
    
    svr.slPerPack{v}=find(diff(svr.slOrd{v})<0);
    svr.slPerPack{v}=cat(1,svr.slPerPack{v}(1),diff(svr.slPerPack{v}),length(svr.slOrd{v})-svr.slPerPack{v}(end));
    cont=1;   
    for p=1:svr.P(v)        
        svr.MSlices{v}{p}=real(zeros([1 1 svr.NSlices{v} 1 1 svr.slPerPack{v}(p)],'like',svr.y{v}));%3rd Dim slices-6th Dim excitations
        svr.MSlicesP{v}{p}=real(zeros([1 1 svr.NSlices{v} 1 1 svr.slPerPack{v}(p)],'like',svr.y{v}));%3rd Dim slices-6th Dim excitations
        svr.MSlicesPlus{v}{p}=real(zeros([1 1 svr.NSlices{v} 1 1 svr.slPerPack{v}(p)],'like',svr.y{v}));%3rd Dim slices-6th Dim excitations
        svr.MSlicesMinu{v}{p}=real(zeros([1 1 svr.NSlices{v} 1 1 svr.slPerPack{v}(p)],'like',svr.y{v}));%3rd Dim slices-6th Dim excitations
        svr.slInPack{v}{p}=svr.slOrd{v}(cont:cont+svr.slPerPack{v}(p)-1);
        for m=1:svr.MBFactor(v)
            svr.MPack{v}(1,1,svr.slInPack{v}{p}+1+(m-1)*svr.NSlices{v}/svr.MBFactor(v),1,p)=1;%3rd Dim slices-5th Dim packages
            for s=1:svr.slPerPack{v}(p)
                svr.MSlices{v}{p}(1,1,svr.slInPack{v}{p}(s)+1+(m-1)*svr.NSlices{v}/svr.MBFactor(v),1,1,s)=1;
                svr.MSlicesP{v}{p}(1,1,svr.slInPack{v}{p}(s)+1+(m-1)*svr.NSlices{v}/svr.MBFactor(v),1,1,s)=1;
                if s~=svr.slPerPack{v}(p);svr.MSlicesPlus{v}{p}(1,1,svr.slInPack{v}{p}(s+1)+1+(m-1)*svr.NSlices{v}/svr.MBFactor(v),1,1,s)=0.5;end%%IF NON-MULTIBANDED WE LINK!
                if s~=1;svr.MSlicesMinu{v}{p}(1,1,svr.slInPack{v}{p}(s-1)+1+(m-1)*svr.NSlices{v}/svr.MBFactor(v),1,1,s)=0.5;end%%IF NON-MULTIBANDED WE LINK!        
            end                     
        end
        %if svr.MBFactor(v)==1;svr.MSlicesP{v}{p}=svr.MSlicesP{v}{p}+svr.MSlicesPlus{v}{p}+svr.MSlicesMinu{v}{p};end
        cont=cont+svr.slPerPack{v}(p);
    end
end