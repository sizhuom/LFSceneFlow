global LFTopDir;
projName = '0611-sceneflow';
workDir = fullfile(LFTopDir, 'Images/Sim', projName);
subDir = '.';
lfFilePrefix = 'frame-480';
motionFile = 'motion.mat';

lf0param = LFReadMetadata(fullfile(workDir, subDir, [lfFilePrefix sprintf('_%04d.json', 0)]));

SY = lf0param.camParam.resol(1); 
SX = lf0param.camParam.resol(2);
SV = lf0param.camParam.resol(3);
SU = lf0param.camParam.resol(4);
if isfield(lf0param.camParam, 'H')
    H = lf0param.camParam.H;
else
    H = genIntrinsics2(lf0param.camParam.resol, lf0param.camParam.apert,...
        lf0param.camParam.fov, lf0param.camParam.fLen);
end
HM = LFMotionMatrix(H, lf0param.camParam.resol, 1000, 1000);

flow_prev = rand(SY,SX,SV,SU,3);

%% Method 1 
tic;
flow_proj = zeros(SY*SX*SV*SU,4);
nPixels = SY*SX*SV*SU;
for ind = 1:nPixels
    flow_proj(ind,:) = HM(:,1:3,ind) * flow_prev(ind:nPixels:end)';
end
flow_proj = reshape(flow_proj, [SY SX SV SU 4]);
toc;

%% Method 2
tic;
intr = struct('H', H, 'sz', [SY SX SV SU], 'S', 1000, 'D', 1000);
flow_proj2 = premultHM(flow_prev, intr);
toc;

%% Method 1
tic;
Ix = rand(SY,SX,SV,SU);
Iy = rand(SY,SX,SV,SU);
Iu = rand(SY,SX,SV,SU);
Iv = rand(SY,SX,SV,SU);
IX = zeros(SY,SX,SV,SU);
IY = zeros(SY,SX,SV,SU);
IZ = zeros(SY,SX,SV,SU);
for ind = 1:size(HM,3)
    f = [Ix(ind) Iy(ind) Iu(ind) Iv(ind)] * HM(:,1:3,ind);
    IX(ind) = f(1);
    IY(ind) = f(2);
    IZ(ind) = f(3);
end
toc;

%% Method 2
tic;
IXYZ = postmultHM(cat(5,Ix,Iy,Iu,Iv), intr);
IX2 = IXYZ(:,:,:,:,1);
IY2 = IXYZ(:,:,:,:,2);
IZ2 = IXYZ(:,:,:,:,3);
errX = IX2 - IX; disp(min(errX(:))); disp(max(errX(:)));
errY = IY2 - IY; disp(min(errY(:))); disp(max(errY(:)));
errZ = IZ2 - IZ; disp(min(errZ(:))); disp(max(errZ(:)));
toc;