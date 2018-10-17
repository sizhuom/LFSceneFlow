function [ flow ] = convertVsf( jsonFile0, resultDir )
%CONVERTVSF Convert D0,Dt,U,V to 3D scene flow

D0 = -csvread(fullfile(resultDir, 'D0_SF_final.csv'));
D1 = -csvread(fullfile(resultDir, 'Dt_SF_final.csv'));
dI = csvread(fullfile(resultDir, 'U_SF_final.csv'));
dJ = csvread(fullfile(resultDir, 'V_SF_final.csv'));

lfparam = LFReadMetadata(jsonFile0);
if strcmp(lfparam.type, 'stereo2')
    K0 = inv(lfparam.camParam.K0);
    K1 = inv(lfparam.camParam.K1);
    Do = K0(1,3) - K1(1,3);
else
    Do = 0;
end
K = inv(lfparam.camParam.K);
Kinv = lfparam.camParam.K;
SU = lfparam.camParam.resol(2);
SV = lfparam.camParam.resol(1);

f = K(1,1);
b = lfparam.camParam.baseline;
Z0 = f * b ./ (D0+Do);
Z1 = f * b ./ (D1+Do);
[I0, J0] = meshgrid(1:SU, 1:SV);
I1 = I0 + dI;
J1 = J0 + dJ;

U0 = I0 * Kinv(1,1) + Kinv(1,3);
V0 = J0 * Kinv(2,2) + Kinv(2,3);
X0 = U0 .* Z0;
Y0 = V0 .* Z0;

Z1w = interp2(Z1,I1,J1,'bicubic',NaN);
U1 = I1 * Kinv(1,1) + Kinv(1,3);
V1 = J1 * Kinv(2,2) + Kinv(2,3);
X1 = U1 .* Z1w;
Y1 = V1 .* Z1w;

flowx = X1 - X0;
flowy = Y1 - Y0;
flowz = Z1w - Z0;
flow = cat(3, flowx, flowy, flowz);

end

