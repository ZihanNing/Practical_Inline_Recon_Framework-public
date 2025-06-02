function [MS,MT]=generateNIIInformation(Par)

%GENERATENIIINFORMATION   Generates the geometric information necessary to
%write NII files
%   [MS,MT]=GENERATENIIINFORMATION(PAR)
%   * PAR are the parameters of a reconstruction structure
%   ** MS is the spacing
%   ** MT is the orientation
%

MT=dynInd(Par.Mine.APhiRec,Par.Mine.Nat,3);

MS=sqrt(sum(MT(1:3,1:3).^2,1));%MS=rec.Enc.AcqVoxelSize;

MTT=[-1 0 0 0;0 -1 0 0;0 0 1 0;0 0 0 1];
MT=MTT*MT;
if ismember(Par.Mine.Modal,9:10);MS(4)=Par.Labels.RepetitionTime(1)/1000;else MS(4)=1;end
