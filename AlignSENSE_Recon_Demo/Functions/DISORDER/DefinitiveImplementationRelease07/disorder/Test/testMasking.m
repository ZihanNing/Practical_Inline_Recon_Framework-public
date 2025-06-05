gamma=0.25;
N=10000;
M=gamma*N;
x=randn(M,N)/sqrt(N);
y=eig(x*x');
figure;
histogram(y,200,'Normalization','pdf')
hold on
z=0.05:0.05:4;w=marcenkoPastur(z,gamma);
plot(z,w);
xb=x;
xb(1:2:end,:)=0;
y=eig(xb*xb');
figure;
histogram(y,200,'Normalization','pdf')
hold on
w=marcenkoPastur(z,gamma/2)/2;
plot(z,w)
w=marcenkoPastur(z,2/gamma)/gamma;
plot(z,w)
y=eig(x'*x);
figure;
histogram(y,200,'Normalization','pdf')
hold on
w=marcenkoPastur(z,1/gamma);
plot(z,w)