function preprocess( workDir )
%PREPROCESS Preprocessing the dataset: generate .json metafile and create
%alpha maps using cocolib

listing = dir(fullfile(workDir, 'IMG_*.LFR'));
for i = 1:numel(listing)
    lffile = fullfile(workDir,listing(i).name);
    [~,lfname,~] = fileparts(listing(i).name);
    [lf,H] = gcReadSubImages(lffile);
    camParam = struct('resol',[size(lf,1) size(lf,2) size(lf,3) size(lf,4)],...
        'H', H);
    param = struct('type','illum','camParam',camParam);
    LFWriteMetadata(fullfile(workDir,[lfname '.json']),param);
end
for i = 1:numel(listing)
    lffile = fullfile(workDir,listing(i).name);
    [~,lfname,~] = fileparts(listing(i).name);
    customParam = struct();
    customParam.lf_path = fullfile(workDir,'SubAperture/');
    customParam.lf_name = [lfname,'_Sub_%d_%d.png'];
    cocolibDisparity(lffile,customParam);
end

end

