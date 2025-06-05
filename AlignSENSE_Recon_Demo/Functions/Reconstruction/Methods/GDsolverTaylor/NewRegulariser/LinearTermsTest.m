
N=101;
x = 1:N;
x = x.^2+3;
X = fft(x);
disp = 44;
lin = (1:N(1)) - ( floor(N(1)/2)+1 ); 
lin = ifftshift(lin)

Xdisp = X .* exp(2*pi*1i*disp*lin/N);
xdisp = ifft(Xdisp);

x
xdisp
figure;plot(x)
figure;plot(xdisp)

