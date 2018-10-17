function [ H, params ] = genIntrinsics2( resol, apert, fov, fLen )
%GENINTRINSICS2 A better way to generate the camera intrinsic matrix

% generate the physical parameters
assert(resol(1)==resol(2));
mMratio = 0.63; % size ratio between sensor plane and main lens, assume fixed
params.np = resol(1);
params.nu = resol(4);
params.nv = resol(3);
params.sM = apert;
params.sp = params.sM * mMratio / params.nu / params.np;
params.dM = params.sM * mMratio / 2 * cot(fov/180*pi/2);
params.dm = params.dM / params.sM * (params.sp*params.np);
params.sm = params.sp * params.np * params.dM / (params.dM+params.dm);
params.fM = fLen;
params.D = 1;

% disp(params);

H = genIntrinsics(params);

end

