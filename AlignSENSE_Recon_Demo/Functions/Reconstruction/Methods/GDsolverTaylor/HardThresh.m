function x = HardThresh(y,t)
% SoftThresh -- Apply Hard Threshold 
%  Usage 
%    x = HardThresh(y,t)
%  Inputs 
%    y     Noisy Data 
%    t     Threshold
%  Outputs 
%    x     

    x = y;
	x(abs(x)<t) = 0;
    
end