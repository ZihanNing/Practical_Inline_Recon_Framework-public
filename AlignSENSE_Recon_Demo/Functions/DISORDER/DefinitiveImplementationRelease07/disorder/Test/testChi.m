

x=0:0.001:32;
ko=16;
figure
for p=[0.7 1 1.4 2]
    k=ko;
    %k=ko^(2/p);

    logy=(1-k/2)*log(2)-log(p)-gammaln(k/2)+(k/p-1)*log(x)-x.^(2/p)/2;
    %y=exp(logy);
    y=((2/p)/(2^(k/2)*gamma(k/2)))*x.^(k/p-1).*exp(-x.^(2/p)/2);
    
    (sum(y(1:end-1).*diff(x))+sum(y(2:end).*diff(x)))/2
    
    hold on
    plot(x,y)
end
legend('p=0.7','p=1','p=1.4','p=2')