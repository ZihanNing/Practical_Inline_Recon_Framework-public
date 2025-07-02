function wout = softthresh(win,lambda)
% Perform soft thresholding, which is the proximal operator of the l1 norm

wout = sign(win) .* max(abs(win)-lambda, 0);

end

