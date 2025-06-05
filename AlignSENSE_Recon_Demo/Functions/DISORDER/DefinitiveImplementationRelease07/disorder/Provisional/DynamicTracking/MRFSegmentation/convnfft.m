function x=convnfft(x,y)

%CONVNFFT performs multiple N-D convolutions using the fft. It replicates 
%the behaviour of convn(a,b,'same') but for the singleton dimensions of x,
%for which the convolution with the elements in y is performed separately
%   X=CONVNFFT(X,Y)
%   * X is an N-D array
%   * Y is a kernel array
%   * X is the convolved array
%

N=size(x);
ND=numDims(x);
M=size(y);M(end+1:ND)=1;
if ND==0;x=bsxfun(@times,x,y);return;end

N=N(1:ND);M=M(1:ND);
O=N+M-1;
y=padarray(y,ceil((O-M)/2),0,'pre');
y=padarray(y,floor((O-M)/2),0,'post');
x=padarray(x,O-N,0,'post');
for n=1:ND 
    if mod(M(n),2)==0;y=fftshift(y,n);else y=ifftshift(y,n);end
    x=fftGPU(x,n);y=fftGPU(y,n);
end
x=bsxfun(@times,x,y);
nV=cell(1,ND);
for n=1:ND
    x=ifftGPU(x,n);
    nV{n}=1:N(n);
end
x=dynInd(x,nV,1:ND);
   