
function [euler] = T_to_Euler (T , euler_order)

    %%% Put the two T parameters in right order in a matrix multiplication

    di = 0; % Backwards
    R_1 = T_to_R ( dynInd( T, 1,4),  di); 

    di = 1; % Forward
    R_2 = T_to_R ( dynInd( T, 2,4) , di);

    R = R_2*R_1;

    %%% Convert to Euler angles in the specified order
    tol = 0.1;
    ichk = 1;

    if isequal( euler_order , 'XYZ'); euler_order = '123';
    elseif isequal( euler_order , 'XZY'); euler_order = '132';
    elseif isequal( euler_order , 'YXZ'); euler_order = '213';
    elseif isequal( euler_order , 'YZX'); euler_order = '231';
    elseif isequal( euler_order , 'ZXY'); euler_order = '312';
    elseif isequal( euler_order , 'ZYX'); euler_order = '321';
    end

    euler = SpinCalc(strcat('DCMtoEA', euler_order), inv(R(1:3, 1:3)), tol, ichk); % DCMtoEA converts from orthonormal rotation matrix to euler angles
    % you have to take the inverse of R - see SpinCalc.m :  DCM - 3x3xN multidimensional matrix which pre-multiplies a coordinate
%              frame column vector to calculate its coordinates in the desired 
%              new frame.  YB: This is the inverse of our definition!
    
    %%% Make sure rotations are in the range [-180 180]
    idx = euler > 180;
    euler(idx) = euler(idx) - 360;

    idx = euler < -180;
    euler(idx) = euler(idx) + 360;
    
    