function bhat=estbeta(m1,m2)

% ESTBETA Estimate the beta parameter of a generalized Gaussian
%   BETA=ESTBETA(M1,M2) estimates de beta parameter of a generalized
%   Gaussian (based on moments). The code is based on Minh N. Do, Dec. 
%   1999 version
%   M1 is the mean absolute value
%   M2 is the variance
%   BHAT is the estimated beta


% beta = F^(-1)(m1^2 / m2);
x=0.05*(1:100)';
f=[2.4663e-05 .004612 .026349 .062937 .10606 .15013 .19234 .23155 .26742 .3 .32951 .35624 .38047 .40247 .42251 .44079 .45753 .47287 .48699 .5 .51202 .52316 .53349 .5431 .55206 .56042 .56823 .57556 .58243 .58888 .59495 .60068 .60607 .61118 .616 .62057 .6249 .62901 .63291 .63662 .64015 .64352 .64672 .64978 .65271 .6555 .65817 .66073 .66317 .66552 .66777 .66993 .672 .674 .67591 .67775 .67953 .68123 .68288 .68446 .68599 .68747 .68889 .69026 .69159 .69287 .69412 .69532 .69648 .6976 .69869 .69974 .70076 .70175 .70271 .70365 .70455 .70543 .70628 .70711 .70791 .70869 .70945 .71019 .71091 .71161 .71229 .71295 .71359 .71422 .71483 .71543 .71601 .71657 .71712 .71766 .71819 .7187 .7192 .71968];

val=(m1.^2)./m2;
val(val<=f(1))=f(1)+1e-5;
val(val>=f(end))=f(end)-1e-4;
bhat=reshape(interp1(f,x,val(:)),size(val));