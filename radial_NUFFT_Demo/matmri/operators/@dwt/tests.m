function tests(obj)
% (c) Corey Baron 2022
% TODO: anything that involves the "extra" field in waveletObj is not currently tested

levToTest = {...
    [2],...
    [2 2],...
    [2 2 2],...
    [2 0 1],...
    [0 2]};
families = {'db1', 'db2'};

rng('default');
im = repmat(phantom(32), [1 1 32]);
im = im .* reshape(randn(size(im,3),1), [1 1 size(im,3)]);
    
% Test that dwt_adjoint(dwt(im)) = im for a variety of cases. Should be
% true for both decimated and undecimated transforms.
for useGPU = [0,1]
    for isDec = [0,1]
        for nFam = 1:length(families)
            for nLev = 1:length(levToTest)
                W = dwt(levToTest{nLev},size(im),isDec,families{nFam},useGPU);
                wim = W*im;
                im_a = W'*wim;
                assert(norm(im(:)-im_a(:))/numel(im) < 1e-8, ...
                    sprintf('Adjoint test failed for useGPU = %d, isDec = %d, nFam = %d, nLev = %d\n', useGPU, isDec, nFam, nLev))
            end
        end
    end
end

% Test operations between multiple wavelet tranform objects. Compare the
% wavelet formatted version with the element wise version.
levToTest = {[2 2 2]};
for useGPU = [0,1]
    for isDec = [0,1]
        for nFam = 1:length(families)
            for nLev = 1:length(levToTest)
                for opType = 1:3
                    
                    W = dwt(levToTest{nLev},size(im),isDec,families{nFam},useGPU);
                    wim = W*im;
                    switch opType
                        case 1
                            % A single scalar
                            wim_b = randn;
                            wim_c = wim_b;
                        case 2
                            % Full wavelet obj
                            wim_b = W * randn(size(im));
                            wim_c = wim_b;
                        case 3 
                            % A wavelet obj with same structure as wim, but each level contains only a single scalar
                            wim_b = waveletObj(1, wim);
                            wim_c = wim;
                            wim_b.low = randn;
                            wim_c.low(:) = wim_b.low;
                            for i = 1:length(wim_b.high)
                                for m = 1:length(wim_b.high{i})
                                    wim_b.high{i}{m} = randn;
                                    wim_c.high{i}{m}(:) = wim_b.high{i}{m};
                                end
                            end
                            if ~isempty(wim_b.extra)
                                wim_b.extra = randn;
                                wim_c.extra(:) = wim_b.extra;
                            end
                    end
                    wim2   = wim - wim_b; % This tests uminus, minus, and plus
                    wim2_b = wim(:) - wim_c(:);
                    assert(norm(wim2(:)-wim2_b(:))/numel(wim2(:)) < 1e-8, ...
                        sprintf('Subtraction test failed for useGPU = %d, isDec = %d, nFam = %d, nLev = %d, opType = %d\n', useGPU, isDec, nFam, nLev, opType))
                    wim2   = wim .* wim_b;
                    wim2_b = wim(:) .* wim_c(:);
                    assert(norm(wim2(:)-wim2_b(:))/numel(wim2(:)) < 1e-8, ...
                        sprintf('Multiplication test failed for useGPU = %d, isDec = %d, nFam = %d, nLev = %d, opType = %d\n', useGPU, isDec, nFam, nLev, opType))
                    wim2   = max(wim, wim_b);
                    wim2_b = max(wim(:), wim_c(:));
                    assert(norm(wim2(:)-wim2_b(:))/numel(wim2(:)) < 1e-8, ...
                        sprintf('Max test failed for useGPU = %d, isDec = %d, nFam = %d, nLev = %d, opType = %d\n', useGPU, isDec, nFam, nLev, opType))
                    wim2   = min(wim, wim_b);
                    wim2_b = min(wim(:), wim_c(:));
                    assert(norm(wim2(:)-wim2_b(:))/numel(wim2(:)) < 1e-8, ...
                        sprintf('Min test failed for useGPU = %d, isDec = %d, nFam = %d, nLev = %d, opType = %d\n', useGPU, isDec, nFam, nLev, opType))
                    wim2   = abs(wim);
                    wim2_b = abs(wim(:));
                    assert(norm(wim2(:)-wim2_b(:))/numel(wim2(:)) < 1e-8, ...
                        sprintf('abs test failed for useGPU = %d, isDec = %d, nFam = %d, nLev = %d, opType = %d\n', useGPU, isDec, nFam, nLev, opType))
                    wim2   = sign(wim);
                    wim2_b = sign(wim(:));
                    assert(norm(wim2(:)-wim2_b(:))/numel(wim2(:)) < 1e-8, ...
                        sprintf('sign test failed for useGPU = %d, isDec = %d, nFam = %d, nLev = %d, opType = %d\n', useGPU, isDec, nFam, nLev, opType))
                    
                end
            end
        end
    end
end


fprintf('dwt unit test success!\n')
return;










