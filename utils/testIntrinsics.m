inParams.np = 14;
inParams.sp = 9e-6;
inParams.nu = 100;
inParams.nv = 100;
inParams.sm = 125e-6;
inParams.sM = 80e-3 / 4;
inParams.dm = 500e-6;
inParams.dM = 80e-3;
inParams.fM = 80e-3;
inParams.D = 1;

H = genIntrinsics(inParams);
disp(H);