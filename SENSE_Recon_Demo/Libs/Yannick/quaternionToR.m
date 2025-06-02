function rotMatrix =  quaternionToR(Q)
%     """
%     Covert a quaternion into a full three-dimensional rotation matrix.
%  
%     Input
%     :param Q: A 4 element array representing the quaternion (q0,q1,q2,q3) 
%  
%     Output
%     :return: A 3x3 element matrix representing the full 3D rotation matrix. 
%              This rotation matrix converts a point in the local reference 
%              frame to a point in the global reference frame.
%     """
    %%%# Extract the values from Q
    if length(Q)==4
        a = Q(1);
        b = Q(2);
        c = Q(3);
        d = Q(4);
    elseif length(Q)==3
        a = sqrt( 1 -  ( Q(1).^2 + Q(2).^2 + Q(3).^2));
        b = Q(1);
        c = Q(2);
        d = Q(3);
    else
        error('quaternionToR:: Input not right.')
    end
     
    if any(~isreal([a b c d]))
        warning('quaternionToR:: complex element found. Real part taken.')
        a = real(a); b=real(b); c=real(c); d=real(d);
    end
     
    %%%# First row of the rotation matrix
    r00 = a^2+b^2-c^2-d^2;
    r01 = 2 * (b * c - a * d);
    r02 = 2 * (b * d + a * c);
     
    %%%# Second row of the rotation matrix
    r10 = 2 * (b * c + a * d);
    r11 = a^2+c^2-b^2-d^2;
    r12 = 2 * (c * d - a * b);
     
    %%%# Third row of the rotation matrix
    r20 = 2 * (b * d - a * c);
    r21 = 2 * (c * d + a * b);
    r22 = a^2+d^2-b^2-c^2;
     
    %%%# 3x3 rotation matrix
    rotMatrix =    [      [r00, r01, r02];
                           [r10, r11, r12];
                           [r20, r21, r22]  ];
%Should be the same as below                            
%rot_matrixTest = SpinCalc('QtoDCM',[Q(2:4) Q(1)],.1,1)' %Transpose because
%of DCM convention