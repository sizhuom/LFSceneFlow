function [ H ] = genIntrinsicsP( params )
%GENINTRINSICSP Generate the camera intrinsic matrix for the perspective
%part

% Parameters
np = 14;
sp = 9e-6;
nu = 100;
nv = 100;
sm = 125e-6;
sM = 80e-3 / 4; % unused in the formula
dm = 500e-6;
dM = sM*dm/(sp*np);
fM = 80e-3;
D = 1;
SA = [0 0]; % subaperture location
if (nargin > 0)
    fields = fieldnames(params);
    for i = 1:numel(fields)
        eval([fields{i} '=params.' fields{i} ';']);
    end
end

% H1: location of each pixel
H1 = [sp 0 0;
      0 sp 0;
      0 0 1;];

% center H1
p_cen = [(1+np*nu)/2,(1+np*nv)/2,1]';
offset = [0;0;1] - H1 * p_cen;
H1(1:2,3) = offset(1:2);

% H2: ray before refraction
H2 = [0 0 SA(1);
      0 0 SA(2);
      -1/(dm+dM) 0 SA(1)/(dm+dM);
      0 -1/(dm+dM) SA(2)/(dm+dM);
      0 0 1];
  
% H3: ray after refraction
H3 = [1 0 0 0 0;
      0 1 0 0 0;
      -1/fM 0 1 0 0;
      0 -1/fM 0 1 0;
      0 0 0 0 1;];
  
% H4: light field representation (virtual plane at D)
H4 = [1 0 0 0 0;
      0 1 0 0 0;
      0 0 D 0 0;
      0 0 0 D 0;
      0 0 0 0 1];
  
% Combine them together
% The axes need to be reversed
% Because the pixels are mirrored after read from the sensors
H = H4*H3*H2*H1;
H = diag([-1,-1,-1,-1,1]) * H;

end

