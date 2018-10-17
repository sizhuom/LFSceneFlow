function [ H ] = genIntrinsics( params )
%GENINTRINSICS Generate the camera intrinsic matrix

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
vp = false;
if (nargin > 0)
    fields = fieldnames(params);
    for i = 1:numel(fields)
        eval([fields{i} '=params.' fields{i} ';']);
    end
end

if ~vp
    % H1: rel to abs
    H1 = [1 0 np 0 0;
        0 1 0 np 0;
        0 0 1 0 0;
        0 0 0 1 0;
        0 0 0 0 1;];
    
    % H2: abs to light field ray
    H2 = [sp 0 0 0 0;
        0 sp 0 0 0;
        0 0 sm 0 0;
        0 0 0 sm 0;
        0 0 0 0 1;];
    
    % add the translational term for H21
    % such that the central pixel maps to (0,0,0,0,1);
    H21 = H2 * H1;
    p_cen = [(1+np)/2,(1+np)/2,(1+nu)/2,(1+nv)/2,1]';
    offset = [0;0;0;0;1] - H21 * p_cen;
    H21(1:4,5) = offset(1:4);
else
    % Another way to formalize H21
    % Assume postprocessing to interpolate "virtual" pixels
    H21 = [ sp 0 sm 0 0;
        0 sp 0 sm 0;
        0 0 sm 0 0;
        0 0 0 sm 0;
        0 0 0 0 1];
    p_cen = [(1+np)/2,(1+np)/2,(1+nu)/2,(1+nv)/2,1]';
    offset = [0;0;0;0;1] - H21 * p_cen;
    H21(1:4,5) = offset(1:4);
end

% H3: express the ray as position (x,y,0) and direction (u,v,1)
idm = 1 / dm;
H3 = [1 0 0 0 0;
      0 1 0 0 0;
      -idm 0 idm 0 0;
      0 -idm 0 idm 0;
      0 0 0 0 1];
  
% H4: propagate to the main lens
H4 = [1 0 dm+dM 0 0;
      0 1 0 dm+dM 0;
      0 0 1 0 0;
      0 0 0 1 0;
      0 0 0 0 1];
  
% H5: apply the main lens
ifM = 1 / fM;
H5 = [1 0 0 0 0;
      0 1 0 0 0;
      -ifM 0 1 0 0;
      0 -ifM 0 1 0;
      0 0 0 0 1];
  
% H6: convert back to light field ray representation
% (two-plane relative)
H6 = [1 0 0 0 0;
      0 1 0 0 0;
      0 0 D 0 0;
      0 0 0 D 0;
      0 0 0 0 1];
  
% Combine them together
% The axes need to be reversed
% Because the pixels are mirrored after read from the sensors
H = H6*H5*H4*H3*H21;
H = diag([-1,-1,-1,-1,1]) * H;

end

