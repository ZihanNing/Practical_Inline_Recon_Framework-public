
function [kIndex, idx, hitMat, timeMat] = PESamplesFromRecPTTest(rec)

%ACQUISITION STRUCTURE
NY = size(rec.y);
if ~isfield(rec.Enc,'kRange');rec.Enc.kRange={[-NY(1) NY(1)-1],[-NY(2)/2 NY(2)/2-1],[-NY(3)/2 NY(3)/2-1]};end %YB: RO-PE-SL. Factor 2 comes from the oversampling in the readout direction (1st dimension here)

perm=1:12; perm(1:3)=[2 3 1]; %If you follow solveXTB logic
NProfs = length(rec.Assign.z{2});
kTraj=zeros([NProfs 2],'single');
for n=1:2;kTraj(:,n)=rec.Assign.z{perm(n)}(:);end

kRange=single(zeros(2,2));
for n=1:2;kRange(n,:)=rec.Enc.kRange{perm(n)};end
kShift=(floor((diff(kRange,1,2)+1)/2)+1)';
kIndex=bsxfun(@plus,kTraj,kShift);
%%%take fft (rec.y) and then fftfift - now ready to use sampling

rec.y = permute(rec.y, perm);
NY = size(rec.y);
idx = sub2ind(NY(1:2),kIndex(:,1),kIndex(:,2))';
  
timeMat = zeros([prod(NY(1:2)) 1]);
timeMat(idx) = 1:length(idx);
timeMat = resSub(timeMat, 1, NY(1:2)); 
hitMat = single(timeMat>0);
