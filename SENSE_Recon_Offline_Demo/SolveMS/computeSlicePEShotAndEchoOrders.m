function [SliceOrder,PEOrder,ShotOrder,EchoOrder,SlicePositions,SliceOffsets,iSlicePositionsMin,SliceOrderReduced,ShotOrderReduced]=computeSlicePEShotAndEchoOrders(Par,NX)

%COMPUTESLICEPESHOTANDECHOORDERS   Computes Slice, PE, Shot and Echo
%ordering in time for Siemens data, as well as SlicePositions for
%reordering data in the slice dimension, and Slice and Shot orders grouped
%by shots
%   [SLICEORDER,PEORDER,SHOTORDER,ECHOORDER,SLICEPOSITIONS,SLICEOFFSETS,SLICEORDERREDUCED,SHOTORDERREDUCED]=COMPUTESLICEPESHOTANDECHOORDERS(PAR)
%   * PAR is a Twix.Image structure
%   * NX is the size of the data
%   * SLICEORDER is the slice order in time
%   * PEORDER is the PE order in time (with PEs from 1 to N in fftshifted
%   space)
%   * SHOTORDER is the shot order in time
%   * ECHOORDER is the echo order in time
%   * SLICEPOSITIONS serves to reorder slices from the way they are stored 
%   in Twix
%   * SLICEOFFSETS slice offsets with respect to the isocenter in radians
%   * SLICEORDERREDUCED is the slice order when grouping k-space samples
%   shotwise
%   * SHOTORDERREDUCED is the shot order when grouping k-space samples
%   shotwise
%

ND=16;NX(end+1:ND)=1;

Lin=Par.Lin;Sli=Par.Sli;Seg=Par.Seg;slicePos=Par.slicePos;Rep=Par.Rep;
NP=length(Lin);fprintf('Number of profiles: %d\n',NP);
NE=Par.NSeg;fprintf('Number of echoes: %d\n',NE);
NKS=length(unique(Lin));fprintf('Number of sampled k-space lines: %d\n',NKS);
NKT=Par.NLin;fprintf('Number of targeted k-space lines: %d\n',NKT);
NSl=length(unique(Sli));fprintf('Number of slices: %d\n',NSl);
NRep=length(unique(Rep));fprintf('Number of repeats: %d\n',NRep);

minLin=min(Lin);minSli=min(Sli);minRep=min(Rep);minSeg=min(Seg);%Minimum line, minimum slice, minimum repeat and minimum echo
if NX(11)==Par.NSeg
    iminLin=find(Lin==minLin & Rep==minRep & Seg==minSeg);iminSli=find(Sli==minSli & Rep==minRep & Seg==minSeg);
else 
    iminLin=find(Lin==minLin & Rep==minRep);iminSli=find(Sli==minSli & Rep==minRep);%Positions in time where minimum line and minimum slice have been sampled
end
%In case of repeated indexes:
[~,ia]=unique(Sli(iminLin));iminLin=iminLin(ia);
SliList=Sli(iminLin);%List of slices
SliceLocations=slicePos(:,iminLin);%Offset of list of slices (spatial coordinates)
SliceNormal=SliceLocations(1:3,SliList==2)-SliceLocations(1:3,SliList==1);%Slice normal
SliceNormalMagnitude=norm(SliceNormal);
SliceNormal=SliceNormal/SliceNormalMagnitude;%Slice normal unit vector
if size(SliceNormal,2)~=1;warning('Slice normal not clearly identified');end
SlicePositions=sum(SliceLocations(1:3,:).*SliceNormal(1:3,1),1);%Slice locations along the normal
NSl=size(SliceLocations,2);
SliceFOV=(max(SlicePositions)-min(SlicePositions))*NSl/(NSl-1);
SliceOffsets=2*pi*SlicePositions/SliceFOV;%Offset in radians
[SlicePositionsMin,iSlicePositionsMin]=min(SlicePositions);%Minimum and index of minimum slice location along the normal
iSlicePositionsMin=iminLin(iSlicePositionsMin);%Index of minimum slice location in list of acquired profiles
SlicePositions=SlicePositions-SlicePositionsMin;%Positive slice locations (i.e., with origin in left-most slice along the normal)
sliPmin2=mink(SlicePositions,2);%Two left-most slices along the normal
SlicePositions=round(SlicePositions/sliPmin2(2))+1;%Normalized slice locations starting at 1, i.e., indexes of slice locations
SliceOrder=SlicePositions(Sli);%Actual mapping to order of slices in space
PEOrder=Lin;%Phase encoding order
EchoOrder=Seg;%Echo order

if NX(11)==Par.NSeg
    PEsl=Lin(Sli==1 & Rep==1 & Seg==1);%Phase encodes for a given slice and repeat
    SEsl=Seg(Sli==1 & Rep==1 & Seg==1);%Segments for a given slice and repeat
else
    PEsl=Lin(Sli==1 & Rep==1);%Phase encodes for a given slice and repeat
    SEsl=Seg(Sli==1 & Rep==1);%Segments for a given slice and repeat
end
NEperShot=diff([0 find(diff(SEsl)<0) NKS]);%Number of echoes per shot
NSh=length(NEperShot);fprintf('Number of shots: %d\n',NSh);%Number of shots
if any(NEperShot<0);
    warning('Number of echoes per shot not clearly identified');
    ShotOrder=[];
    SliceOrderReduced=[];
    ShotOrderReduced=[];
else
    ShotsKS=repelem(1:NSh,NEperShot);%Shots for sampled k-space lines (careful, in the temporal order they are sampled!)
    ShotsKT=zeros(1,NKT);ShotsKT(PEsl)=ShotsKS;%Shots for targeted k-space lines
    ShotOrder=ShotsKT(Lin);%Shot order

    for e=1:NE;if sum(Seg==e)==NSh*NSl;break;end;end%Identifying an echo that has been sampled for all shots and slices
    SliceOrderReduced=SliceOrder(Seg==e);%Slice order when grouping k-space samples belonging to same shot
    ShotOrderReduced=ShotOrder(Seg==e);%Shot order when grouping k-space samples belonging to same shot
end
