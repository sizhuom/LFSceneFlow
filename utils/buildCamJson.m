function [ H ] = buildCamJson( sFile, oFile, lFile, focus, fLen )
%BUILDCAMJSON Build cam from a .json profile
global LFOptDir
sensor = LFReadMetadata(fullfile(LFOptDir,'sensor',[sFile '.json']));
objective = LFReadMetadata(fullfile(LFOptDir,'objective',[oFile '.json']));
lenslet = LFReadMetadata(fullfile(LFOptDir,'lenslet',[lFile '.json']));
if nargin == 5
    H = buildCam(sensor, objective, lenslet, focus, fLen);
else
    H = buildCam(sensor, objective, lenslet, focus);
end

end

