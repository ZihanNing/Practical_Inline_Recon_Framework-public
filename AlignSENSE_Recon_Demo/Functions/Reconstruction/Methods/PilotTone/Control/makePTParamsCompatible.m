
function rec = makePTParamsCompatible( rec, resPyr, resIso, estT, estB)


L = length(resPyr);

if length( rec.Alg.parXT.PT.usePT) < L
    rec.Alg.parXT.PT.usePT = cat( 2, ...
                                  rec.Alg.parXT.PT.usePT, ...
                                  rec.Alg.parXT.PT.usePT(end)*ones([1 L-length( rec.Alg.parXT.PT.usePT)])  );
end

if length( rec.Alg.parXT.PT.subDiv) < L
    rec.Alg.parXT.PT.subDiv = cat( 2, ...
                                  rec.Alg.parXT.PT.subDiv, ...
                                  rec.Alg.parXT.PT.subDiv(end)*ones([1 L-length( rec.Alg.parXT.PT.subDiv)])  );
end

if length( rec.Alg.parXT.PT.nItActivate) < L
    rec.Alg.parXT.PT.nItActivate = cat( 2, ...
                                  rec.Alg.parXT.PT.nItActivate, ...
                                  rec.Alg.parXT.PT.nItActivate(end)*ones([1 L-length( rec.Alg.parXT.PT.nItActivate)])  );
end

