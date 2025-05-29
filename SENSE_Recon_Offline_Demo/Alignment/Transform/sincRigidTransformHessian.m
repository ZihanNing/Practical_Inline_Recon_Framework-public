function H=sincRigidTransformHessian(xB,GB,GC,et,etg,eth,mH,F,FH)

%SINCRIGIDTRANSFORMHESSIAN   Computes the Hessian of sinc-interpolated 
%3D rigid transforms
%   G=SINCRIGIDTRANSFORMHESSIAN(XB,GB,GC,ET,ETG,ETH,MH,{F},{FH}) 
%   obtains the Hessian of the transform of the images
%   * XB contains the images before the first, second and third rotations 
%   and before the translation
%   * GB contains gradient of the first, second and third rotations before 
%   applying the translation
%   * GC contains the gradient of the first, first and second rotations 
%   before applying the second, third and third rotations respectively
%   * ET are the transform factors
%   * ETG are the transform gradient factors
%   * ETH are the transform Hessian factors
%   * MH indexes the Hessian term to be computed
%   * {F} contains discrete Fourier transform matrices
%   * {FH} contains inverse discrete Fourier transform matrices
%   ** G is the Hessian of the transformed image
%

if nargin<8 || isempty(F);F={[],[],[]};end
if nargin<9 || isempty(FH);FH={[],[],[]};end

%Translation parameters
for m=1:6
    if m==mH
        x{1}=bsxfun(@times,xB{4},eth{1}{m});
        if iscell(et{1})  
            for n=1:3;x{1}=bsxfun(@times,x{1},et{1}{n});end
        end
        for n=1:3;x{1}=ifftGPU(x{1},n,FH{n});end  
        H=x{1};
    end
end

%Translation-rotation parameters
%First rotation
for n=1:3
    for m=1:3
        if (6+3*(n-1)+m)==mH
            x{1}=bsxfun(@times,GB{n},etg{1}{m});
            if iscell(et{1})  
                for o=1:3;x{1}=bsxfun(@times,x{1},et{1}{o});end
            end
            for o=1:3;x{1}=ifftGPU(x{1},o,FH{o});end   
            H=x{1};
        end
    end 
end

%Rotation cross terms
%First-second
if mH==16
    x{1}=bsxfun(@times,GC{1},et{2}{2});
    x{2}=bsxfun(@times,GC{1},etg{2}{2});
    for m=1:2
        x{m}=ifftGPU(x{m},3,FH{3});
        x{m}=fftGPU(x{m},1,F{1});
    end
    x{2}=bsxfun(@times,etg{3}{2},x{1})+bsxfun(@times,et{3}{2},x{2});
    x{1}=bsxfun(@times,et{3}{2},x{1});
    for m=1:2
        x{m}=ifftGPU(x{m},1,FH{1});
        x{m}=fftGPU(x{m},3,F{3});
    end
    x{1}=bsxfun(@times,etg{2}{2},x{1})+bsxfun(@times,et{2}{2},x{2});
    x{1}=ifftGPU(x{1},3,FH{3});
    
    x{1}=fftGPU(x{1},2,F{2});
    x{1}=bsxfun(@times,x{1},et{2}{3});
    x{1}=ifftGPU(x{1},2,FH{2});
    x{1}=fftGPU(x{1},3,F{3});
    x{1}=bsxfun(@times,x{1},et{3}{3});
    x{1}=ifftGPU(x{1},3,FH{3});
    x{1}=fftGPU(x{1},2,F{2});
    x{1}=bsxfun(@times,x{1},et{2}{3});

    for m=1:2:3;x{1}=fftGPU(x{1},m,F{m});end
    if ~iscell(et{1})
        x{1}=bsxfun(@times,x{1},et{1});
    else
        for m=1:3;x{1}=bsxfun(@times,x{1},et{1}{m});end
    end
    for m=3:-1:1;x{1}=ifftGPU(x{1},m,FH{m});end
    H=x{1};
end


%First-third
if mH==17
    x{2}=bsxfun(@times,GC{2},etg{2}{3});
    x{1}=bsxfun(@times,GC{2},et{2}{3});
    for m=1:2
        x{m}=ifftGPU(x{m},2,FH{2});
        x{m}=fftGPU(x{m},3,F{3});
    end
    x{2}=bsxfun(@times,etg{3}{3},x{1})+bsxfun(@times,et{3}{3},x{2});
    x{1}=bsxfun(@times,et{3}{3},x{1});
    for m=1:2
        x{m}=ifftGPU(x{m},3,FH{3});
        x{m}=fftGPU(x{m},2,F{2});
    end
    x{1}=bsxfun(@times,etg{2}{3},x{1})+bsxfun(@times,et{2}{3},x{2});

    for m=1:2:3;x{1}=fftGPU(x{1},m,F{m});end
    if ~iscell(et{1})
        x{1}=bsxfun(@times,x{1},et{1});
    else
        for m=1:3;x{1}=bsxfun(@times,x{1},et{1}{m});end
    end    
    for m=3:-1:1;x{1}=ifftGPU(x{1},m,FH{m});end
    H=x{1};
end

%Second-third
if mH==18
    x{2}=bsxfun(@times,GC{3},etg{2}{3});
    x{1}=bsxfun(@times,GC{3},et{2}{3});
    for m=1:2
        x{m}=ifftGPU(x{m},2,FH{2});
        x{m}=fftGPU(x{m},3,F{3});
    end
    x{2}=bsxfun(@times,etg{3}{3},x{1})+bsxfun(@times,et{3}{3},x{2});
    x{1}=bsxfun(@times,et{3}{3},x{1});
    for m=1:2
        x{m}=ifftGPU(x{m},3,FH{3});
        x{m}=fftGPU(x{m},2,F{2});
    end
    x{1}=bsxfun(@times,etg{2}{3},x{1})+bsxfun(@times,et{2}{3},x{2});

    for m=1:2:3;x{1}=fftGPU(x{1},m,F{m});end
    if ~iscell(et{1})
        x{1}=bsxfun(@times,x{1},et{1});
    else
        for m=1:3;x{1}=bsxfun(@times,x{1},et{1}{m});end
    end
    for m=3:-1:1;x{1}=ifftGPU(x{1},m,FH{m});end
    H=x{1};
end


%Rotation second order
%First rotation
if mH==19
    x{1}=bsxfun(@times,xB{1},et{2}{1});
    x{2}=bsxfun(@times,xB{1},eth{2}{1});
    x{3}=bsxfun(@times,xB{1},etg{2}{1});
    x{4}=bsxfun(@times,xB{1},etg{2}{1});
    for m=1:4
        x{m}=ifftGPU(x{m},1,FH{1});
        x{m}=fftGPU(x{m},2,F{2});
    end
    x{2}=bsxfun(@times,eth{3}{1},x{1})+bsxfun(@times,et{3}{1},x{2});
    x{5}=bsxfun(@times,etg{3}{1},x{1})+bsxfun(@times,et{3}{1},x{3});
    x{6}=bsxfun(@times,etg{3}{1},x{1})+bsxfun(@times,et{3}{1},x{4});    
    x{1}=bsxfun(@times,et{3}{1},x{1});
    x{3}=bsxfun(@times,etg{3}{1},x{3});
    x{4}=bsxfun(@times,etg{3}{1},x{4});
    for m=1:6
        x{m}=ifftGPU(x{m},2,FH{2});
        x{m}=fftGPU(x{m},1,F{1});
    end
    x{1}=bsxfun(@times,eth{2}{1},x{1})+bsxfun(@times,et{2}{1},x{2});
    x{2}=bsxfun(@times,etg{2}{1},x{5})+bsxfun(@times,et{2}{1},x{3});
    x{3}=bsxfun(@times,etg{2}{1},x{6})+bsxfun(@times,et{2}{1},x{4});
    x{1}=x{1}+x{2}+x{3};
    x{1}=ifftGPU(x{1},1,FH{1});
        
    if any(et{5}{2}(:)==1);x{1}=bsxfun(@times,x{1},1-et{5}{2})+bsxfun(@times,flipping(x{1},et{4}{2}),et{5}{2});end
    x{1}=fftGPU(x{1},3,F{3});
    x{1}=bsxfun(@times,x{1},et{2}{2});
    x{1}=ifftGPU(x{1},3,FH{3});
    x{1}=fftGPU(x{1},1,F{1});
    x{1}=bsxfun(@times,x{1},et{3}{2});
    x{1}=ifftGPU(x{1},1,FH{1});
    x{1}=fftGPU(x{1},3,F{3});
    x{1}=bsxfun(@times,x{1},et{2}{2});
    x{1}=ifftGPU(x{1},3,FH{3});
    
    if any(et{5}{3}(:)==1);x{1}=bsxfun(@times,x{1},1-et{5}{3})+bsxfun(@times,flipping(x{1},et{4}{3}),et{5}{3});end
    x{1}=fftGPU(x{1},2,F{2});
    x{1}=bsxfun(@times,x{1},et{2}{3});
    x{1}=ifftGPU(x{1},2,FH{2});
    x{1}=fftGPU(x{1},3,F{3});
    x{1}=bsxfun(@times,x{1},et{3}{3});
    x{1}=ifftGPU(x{1},3,FH{3});
    x{1}=fftGPU(x{1},2,F{2});
    x{1}=bsxfun(@times,x{1},et{2}{3});

    for m=1:2:3;x{1}=fftGPU(x{1},m,F{m});end
    if ~iscell(et{1})
        x{1}=bsxfun(@times,x{1},et{1});
    else
        for m=1:3;x{1}=bsxfun(@times,x{1},et{1}{m});end
    end
    for m=3:-1:1;x{1}=ifftGPU(x{1},m,FH{m});end
    H=x{1};
end


%Second rotation
if mH==20
    x{1}=bsxfun(@times,xB{2},et{2}{2});
    x{2}=bsxfun(@times,xB{2},eth{2}{2});
    x{3}=bsxfun(@times,xB{2},etg{2}{2});
    x{4}=bsxfun(@times,xB{2},etg{2}{2});
    for m=1:4
        x{m}=ifftGPU(x{m},3,FH{3});
        x{m}=fftGPU(x{m},1,F{1});
    end
    x{2}=bsxfun(@times,eth{3}{2},x{1})+bsxfun(@times,et{3}{2},x{2});
    x{5}=bsxfun(@times,etg{3}{2},x{1})+bsxfun(@times,et{3}{2},x{3});
    x{6}=bsxfun(@times,etg{3}{2},x{1})+bsxfun(@times,et{3}{2},x{4});    
    x{1}=bsxfun(@times,et{3}{2},x{1});
    x{3}=bsxfun(@times,etg{3}{2},x{3});
    x{4}=bsxfun(@times,etg{3}{2},x{4});
    for m=1:6
        x{m}=ifftGPU(x{m},1,FH{1});
        x{m}=fftGPU(x{m},3,F{3});
    end
    x{1}=bsxfun(@times,eth{2}{2},x{1})+bsxfun(@times,et{2}{2},x{2});
    x{2}=bsxfun(@times,etg{2}{2},x{5})+bsxfun(@times,et{2}{2},x{3});
    x{3}=bsxfun(@times,etg{2}{2},x{6})+bsxfun(@times,et{2}{2},x{4});
    x{1}=x{1}+x{2}+x{3};
    x{1}=ifftGPU(x{1},3,FH{3});
        
    if any(et{5}{3}(:)==1);x{1}=bsxfun(@times,x{1},1-et{5}{3})+bsxfun(@times,flipping(x{1},et{4}{3}),et{5}{3});end
    x{1}=fftGPU(x{1},2,F{2});
    x{1}=bsxfun(@times,x{1},et{2}{3});
    x{1}=ifftGPU(x{1},2,FH{2});
    x{1}=fftGPU(x{1},3,F{3});
    x{1}=bsxfun(@times,x{1},et{3}{3});
    x{1}=ifftGPU(x{1},3,FH{3});
    x{1}=fftGPU(x{1},2,F{2});
    x{1}=bsxfun(@times,x{1},et{2}{3});

    for m=1:2:3;x{1}=fftGPU(x{1},m,F{m});end
    if ~iscell(et{1})
        x{1}=bsxfun(@times,x{1},et{1});
    else
        for m=1:3;x{1}=bsxfun(@times,x{1},et{1}{m});end
    end
    for m=3:-1:1;x{1}=ifftGPU(x{1},m,FH{m});end
    H=x{1};
end

%Third rotation
if mH==21
    x{1}=bsxfun(@times,xB{3},et{2}{3});
    x{2}=bsxfun(@times,xB{3},eth{2}{3});
    x{3}=bsxfun(@times,xB{3},etg{2}{3});
    x{4}=bsxfun(@times,xB{3},etg{2}{3});
    for m=1:4
        x{m}=ifftGPU(x{m},2,FH{2});
        x{m}=fftGPU(x{m},3,F{3});
    end
    x{2}=bsxfun(@times,eth{3}{3},x{1})+bsxfun(@times,et{3}{3},x{2});
    x{5}=bsxfun(@times,etg{3}{3},x{1})+bsxfun(@times,et{3}{3},x{3});
    x{6}=bsxfun(@times,etg{3}{3},x{1})+bsxfun(@times,et{3}{3},x{4});    
    x{1}=bsxfun(@times,et{3}{3},x{1});
    x{3}=bsxfun(@times,etg{3}{3},x{3});
    x{4}=bsxfun(@times,etg{3}{3},x{4});
    for m=1:6
        x{m}=ifftGPU(x{m},3,FH{3});
        x{m}=fftGPU(x{m},2,F{2});
    end
    x{1}=bsxfun(@times,eth{2}{3},x{1})+bsxfun(@times,et{2}{3},x{2});
    x{2}=bsxfun(@times,etg{2}{3},x{5})+bsxfun(@times,et{2}{3},x{3});
    x{3}=bsxfun(@times,etg{2}{3},x{6})+bsxfun(@times,et{2}{3},x{4});
    x{1}=x{1}+x{2}+x{3};

    for m=1:2:3;x{1}=fftGPU(x{1},m,F{m});end
    if ~iscell(et{1})
        x{1}=bsxfun(@times,x{1},et{1});
    else
        for m=1:3;x{1}=bsxfun(@times,x{1},et{1}{m});end
    end
    for m=3:-1:1;x{1}=ifftGPU(x{1},m,FH{m});end
    H=x{1};
end
