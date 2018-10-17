function [ Hnew ] = resizeIntrinsics( H, scale )
%RESIZEINTRINSICS Modify the intrinsic matrix due to a light field scaling

a = 1 / scale;
b = (1 + a) / 2 - a;
Hadd = [
    1 0 0 0 0;
    0 1 0 0 0;
    0 0 a 0 b;
    0 0 0 a b;
    0 0 0 0 1;
    ];

Hnew = H * Hadd;

end

