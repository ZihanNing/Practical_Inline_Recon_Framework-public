function tests(obj)
% (c) Corey Baron 2020

N = 32;
NR = 4;
for nD = 1:3
    szIm = N*ones(1,nD);
    input = randn([szIm, NR]) + 1i*randn([szIm, NR]);
    x = randn([szIm,1]) + 1i*randn([szIm,1]);
    y = randn([szIm, NR]) + 1i*randn([szIm, NR]);
    for isCalDat = [false,true]
        if isCalDat
            allApN = [N/2,N];
        else
            allApN = N;
        end
        for apN = allApN
            apN_a = apN*ones(1,nD);
            for doLoop = [false,true]
                R = rcvrOp(input,isCalDat,apN_a);
                R.doLoop = doLoop;
                Rx = R*x;
                Ry = R'*y;
                d1 = dot(x(:),Ry(:));
                d2 = dot(Rx(:),y(:));
                assert(abs(d1-d2)/min(abs(d1),abs(d2)) < 1e-8, 'Adjoint test failed. nD = %d', nD)
                if ~doLoop
                    Rx_a = Rx;
                else
                    assert( sum(abs(Rx_a(:)-Rx(:))) < 1e-8, 'doLoop test failed, nD = %d', nD)
                end
            end
        end
    end
end

fprintf('rcvrOp unit test success!\n')

end
