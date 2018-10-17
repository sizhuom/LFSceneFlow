function extractRawImage(workDir, devignetting)
%EXTRACTRAWIMAGE Extract png images from .RAW, .LFR Lytro Illum images

global LFTopDir
WhiteImageDatabasePath = fullfile(LFTopDir,'Cameras/WhiteImageDatabase.mat');
if nargin < 2
    devignetting = true;
end

files = dir(fullfile(workDir,'*.RAW'));
for i = 1:numel(files)
    raw = LFReadRaw(fullfile(workDir,files(i).name), '10bit');
    im = demosaic(raw, 'grbg');
    imd = double(im) / 1023;
    [~,fname,~] = fileparts(files(i).name);
    imwrite(imd, fullfile(workDir,[fname '.png']));
end

files2 = dir(fullfile(workDir,'*.LFR'));
for i =1:numel(files2)
    lfp = LFReadLFP(fullfile(workDir,files2(i).name));
    raw = lfp.RawImg;
    LFMetadata = lfp.Metadata;
    BlackLevel = LFMetadata.image.pixelFormat.black.gr;
    WhiteLevel = LFMetadata.image.pixelFormat.white.gr;
    im = (double(raw) - BlackLevel) ./ (WhiteLevel - BlackLevel);
    if devignetting
        DesiredCam = struct('CamSerial', lfp.Serials.camera.serialNumber, ...
            'ZoomStep', LFMetadata.devices.lens.zoomStep, ...
            'FocusStep', LFMetadata.devices.lens.focusStep );
        WhiteImageInfo = LFSelectFromDatabase( DesiredCam, WhiteImageDatabasePath );
        [WhiteImagePath,WhiteImageFname,~] = fileparts(WhiteImageInfo.Fname);
        whiteImage = LFReadRaw(fullfile(LFTopDir,'Cameras',WhiteImagePath,[WhiteImageFname '.RAW']), '10bit');
        whiteImage = (double(whiteImage) - BlackLevel) ./ (WhiteLevel - BlackLevel);
        im = im ./ whiteImage;
    end
    im = im2uint16(im);
    im = demosaic(im, 'grbg');
    [~,fname,~] = fileparts(files2(i).name);
    imwrite(im, fullfile(workDir,[fname '.png']));
end