
load inpRobLoss

    figure
    plot(slR)
    figure
    plot(slW)
    
[b,c]=shapeParameterEstimation(slR,slW);
fprintf('Scale parameter: %.4f\n',c);
fprintf('Shape parameter: %.4f\n',b);

%figure
%plot(x,'*')
%figure
%plot(W,'*')
