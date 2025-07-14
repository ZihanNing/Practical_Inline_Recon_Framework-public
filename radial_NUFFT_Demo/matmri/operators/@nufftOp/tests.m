function tests(obj)
% (c) Corey Baron 2020

N0 = 32;
if gpuDeviceCount>0
    useGPU = 1;
else
    useGPU = 0;
end
for nD = 1:3 
    % Check adjoint for basic cases
    imN = repmat(N0, 1, nD);
    imN = imN + (0:4:4*(nD-1)); % Make each dim have different size
    kloc = (-imN(1)/2:imN(1)/2-1)'/imN(1); % Use cartesian kloc to enable comparison to fft
    switch nD
        case 2
            kloc2 = (-imN(2)/2:imN(2)/2-1)'/imN(2);
            [k1,k2] = ndgrid(kloc,kloc2);
            kloc = [k1(:),k2(:)];
            clear k1 k2
        case 3
            kloc2 = (-imN(2)/2:imN(2)/2-1)'/imN(2);
            kloc3 = (-imN(3)/2:imN(3)/2-1)'/imN(3);
            [k1,k2,k3] = ndgrid(kloc,kloc2,kloc3);
            kloc = [k1(:),k2(:),k3(:)];
            clear k1 k2 k3
    end
    imk = [size(kloc,1), 1];
    dcf = ones(imk);
    S = nufftOp(imN, kloc, dcf, useGPU, 1.5, 6);
    x = randn([imN,1]) + 1i*randn([imN,1]) + randn;
    y = randn([imk,1]) + 1i*randn([imk,1]) + randn;
    Sx = S*x;
    Sy = S'*y;
    d1 = dot(x(:),Sy(:));
    d2 = dot(Sx(:),y(:));
    assert(abs(d1-d2)/min(abs(d1),abs(d2)) < 1e-8, 'Adjoint test failed. nD = %d', nD)

    % Test repetitions over loops
    x = repmat(x, [ones(1,nD),2,2]);
    y = repmat(y, [1,2,2]);
    for loopDim=nD:nD+1
        S.loopDim = loopDim;
        Sx = S*x;
        Sy = S'*y;
        d1 = dot(x(:),Sy(:));
        d2 = dot(Sx(:),y(:));
        assert(abs(d1-d2)/min(abs(d1),abs(d2)) < 1e-8, 'Adjoint test on repetitions failed. nD = %d', nD)
        d1 = mean(reshape(abs(diff(Sx,[],1+1)),1,[])) + ...;
            mean(reshape(abs(diff(Sx,[],1+2)),1,[]));
        d2 = mean(reshape(abs(diff(Sy,[],nD+1)),1,[])) + ...;
            mean(reshape(abs(diff(Sy,[],nD+2)),1,[]));
        assert(mean(d1)+mean(d2) < 1e-8, 'Adjoint test failed. nD = %d', nD)
    end
    clear x y Sx Sy
    
    % Test circular shift property by shifting kloc outside nominal
    % [-0.5,0.5) range
    y = randn([imk,1]) + 1i*randn([imk,1]) + randn;
    S = nufftOp(imN, kloc, dcf, useGPU,1.5,6);
    Sy1 = S'*y;
    S = nufftOp(imN, kloc+0.25, dcf, useGPU,1.5,6);
    Sy2 = S'*y;
    assert(abs(mean(abs(Sy1(:))-abs(Sy2(:)))) < 1e-5, 'Circ shift test abs failed. nD = %d', nD)
    clear Sy1 Sy2
    
    % Point spread function tests
    x = zeros([imN,1]);
    if nD==1
        x(end/2+1) = 1;
    elseif nD==2
        x(end/2+1,end/2+1) = 1;
    elseif nD==3
        x(end/2+1,end/2+1,end/2+1) = 1;
    end
    S = nufftOp(imN, kloc, dcf, useGPU, 1.5, 6);
    Sx1 = S*x;
    testval = mean(abs(angle(Sx1(:))));
    assert(testval < 1e-6, 'Point spread function phase test failed. nD = %d', nD)
    testval = std(abs(Sx1(:)))/mean(abs(Sx1(:))); 
    assert(testval < 1e-4, 'Point spread function abs test failed. nD = %d', nD)
    S = nufftOp(imN, kloc, dcf, useGPU, 1.375, 6);
    Sx2 = S*x;
    S = nufftOp(imN, kloc, dcf, useGPU, 1.375, 4);
    Sx3 = S*x;
    testval = [mean(abs(Sx1(:))), mean(abs(Sx2(:))), mean(abs(Sx3(:)))];
    assert(max(abs(diff(testval))) < 1e-3, 'Scaling test failed. nD = %d', nD)  
    clear Sx1 Sx2 Sx3
    
    % Comparison to fft
    im = phantom(min(imN)); 
    imN_a = [imN, 1, 1];
    for nD_a=1:3
        im = padcrop(im,imN_a(nD_a),nD_a);
    end
    S = nufftOp(imN, kloc, dcf, useGPU, 1.5, 6);
    k1 = reshape(S*im, [imN,1]);
    im1 = fftnc(k1);
    testval = im1(:)-im(:);
    testval = sqrt(sum(testval.*conj(testval)))/numel(testval);
    assert(testval < 1e-4, 'Comparison to fft failed. nD = %d', nD)
end

% Test Toeplitz
for nD = 1:3
    imN = repmat(64,[1 nD]);
    kloc = [];
    for n=1:length(imN)
        kloc = cat(2,kloc, rand(0.25*prod(imN),1)-0.5);
    end
    dcf = rand(size(kloc,1),1);
    useGPU = 0;
    S = nufftOp(imN, kloc, dcf, useGPU,1.5,6);
    S = S.prepToep;
    % Gen data
    x = randn([imN,2]) + 1i*randn([imN,2]);
    S0 = nufftOp(imN, kloc, dcf, useGPU, 2, 8); % "truth" case
    y0 = S0'*(S0*x);
    y1 = S*x;
    % Test
    testval = y0-y1;
    testval = sqrt(mean(testval(:).*conj(testval(:))));
    assert(testval < 1e-3, 'Toeplitz failed. nD = %d', nD)
end

% Test dcf computation
N0 = 64;
dk_base = 1/N0;
kloc_diff = linspace(dk_base/2,dk_base,2*N0);
kloc = cumsum(kloc_diff);
kloc = kloc/max(kloc)-0.5;
S = nufftOp(N0,kloc');
S = S.findDcf;
dcf = S.dcf(1+10:end-1-10);
dcf_theory = diff(kloc(1+10:end-10)');
normVal = sum(dcf(round(2*N0/4):round(2*N0*3/4)))/sum(dcf_theory(round(2*N0/4):round(2*N0*3/4)));
costval = norm(dcf/normVal - dcf_theory)/norm(dcf_theory);
assert(costval < 0.025, 'dcf computation failed')

fprintf('nufftOp unit test success!\n')

end
