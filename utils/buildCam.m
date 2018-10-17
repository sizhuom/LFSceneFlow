function [ camParam ] = buildCam( sensor, objective, lenslet, focus, fLen )
%BUILDCAM Compute the intrinsics based on real camera parameters

param.np = floor(lenslet.pitch / sensor.pixSz);
param.sp = sensor.pixSz;
param.nu = floor(min(sensor.resol(1)*sensor.pixSz,lenslet.width)/lenslet.pitch);
param.nv = floor(min(sensor.resol(2)*sensor.pixSz,lenslet.height)/lenslet.pitch);
param.sm = lenslet.pitch;
param.dm = lenslet.fLen;

if (nargin > 4)
    if numel(objective.fLen) == 1
        if fLen ~= objective.fLen
            error('Error: Specified focal length %fm not available.', fLen);
        else
            param.fM = fLen;
        end
    else
        if fLen >= objective.fLen(1) && fLen <= objective.fLen(2)
            param.fM = fLen;
        else
            error('Error: Specified focal length %fm not available.', fLen);
        end
    end
else
    param.fM = objective.fLen(end);
end
param.dM = 1/(1/param.fM-1/focus);
if (param.dM == Inf || param.dM <= 0)
    error('Error: Invalid objective lens distance: %fm.', param.dM);
end

fNum = param.dm/param.sm;
param.sM = param.dM / fNum;
if numel(objective.maxApert) == 1
    if fNum < objective.maxApert
        error('Error: Required aperture f/%f is not available.', fNum);
    end
else
    if fNum < objective.maxApert(1)
        error('Error: Required aperture f/%f is not available.', fNum);
    end
    if fNum < objective.maxApert(2)
        warning('Required aperture f/%f may not be available. Please check.', fNum);
    end
end
if fNum > objective.minApert
    error('Error: Required aperture f/%f is not available.', fNum);
end
   
param.vp = true;
% param.D = param.fM;
H = genIntrinsics(param);

camParam.H = H;
camParam.resol = [param.np param.np param.nv param.nu];
camParam.apert = param.sM;
camParam.fov = atan(param.sm*param.nu/2/param.dM)*2/pi*180;
camParam.fLen = param.fM;


end

