function H=sliceProfile(N,SlTh,SlOv,parX)

%SLICEPROFILE   Constructs the slice profile filter using a top-hat
%approximation to the continuous profile
%   H=SLICEPROFILE(N,SLTH,SLOV,PARX) constructs a slice profile filter based
%   on slice thickness and slice overlap
%   * N is the number of slices
%   * SLTH is the slice thickness
%   * SLOV is the slice overlap
%   * PARX are the parameters of the reconstruction
%   ** H is the slice filter in Fourier space
%

H=zeros([1 N],'single');
SlSp=SlTh-SlOv;
H(1)=1;
if SlSp>SlTh
    fprintf('Slice thickness: %6.3f\n',SlTh);
    fprintf('Slice gap: %6.3f\n',SlGp);
    error('The slice overlap pattern has not been implemented in the reconstruction'); 
end
if parX.threeD==1
    H(1)=SlSp;
    SlRes=SlOv;
    inda=2;
    indb=0;
    while SlRes>0
        if SlRes-2*SlSp>0
            H(inda)=SlSp;
            H(end-indb)=SlSp;
            inda=inda+1;
            indb=indb+1;
        else
            H(inda)=SlRes/2;
            H(end-indb)=SlRes/2;
        end
        SlRes=SlRes-2*SlSp;
    end
    H=H/SlTh;
end
H=fftGPU(H,2);
H=permute(H,[1 3 2]);   
