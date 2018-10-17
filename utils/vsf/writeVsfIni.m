function iniFile = writeVsfIni( jsonFile0, jsonFile1, resultDir, customParam )
%WRITEVSFINI Write vsf .ini file

[pathstr, fname0, ~] = fileparts(jsonFile0);
[~, fname1, ~] = fileparts(jsonFile1);
lfParam = LFReadMetadata(jsonFile0);
sz = lfParam.camParam.resol;

% Global
secG.Alpha = 250;
secG.AlphaOF = 250;
secG.AlphaST = 250;
secG.AlphaSF = 250;
secG.Lambda = 0.05;
secG.Mu = 0.05;
secG.Gamma = 1.5;

secG.PyramidLevelOF = 18;
secG.PyramidLevelST = 12;
secG.PyramidLevelSF = 12;

secG.PyramidLevelFinal = 1;
secG.PyramidEta = 0.9;
secG.LogFile = fullfile(resultDir, 'vsf.log');
secG.DoStereo = 1;
secG.InitStereoBP = 1;
secG.DoSceneFlow = 1;
secG.Il0_FileName = fullfile(pathstr, [fname0 '-f0.png']);
secG.Ilt_FileName = fullfile(pathstr, [fname1 '-f0.png']);
secG.Ir0_FileName = fullfile(pathstr, [fname0 '-f1.png']);
secG.Irt_FileName = fullfile(pathstr, [fname1 '-f1.png']);
secG.ImagePyramidFileSuffix = 'png';
% secG.Il0PyramidFilePrefix = fullfile(resultDir, 'Debug', 'Il0_');
% secG.IltPyramidFilePrefix = fullfile(resultDir, 'Debug', 'Ilt_');
% secG.Ir0PyramidFilePrefix = fullfile(resultDir, 'Debug', 'Ir0_');
% secG.IrtPyramidFilePrefix = fullfile(resultDir, 'Debug', 'Irt_');

% Optic Flow
% secOF.OF_I0_FileName = fullfile(resultDir, 'Debug', 'I0.png');
% secOF.OF_It_FileName = fullfile(resultDir, 'Debug', 'It.png');
% secOF.OF_It_Warped_FileName = fullfile(resultDir, 'Debug', 'It_warped.png');
% secOF.OF_U_increment_FileName = fullfile(resultDir, 'Debug', 'delta_U_OF.png');
% secOF.OF_V_increment_FileName = fullfile(resultDir, 'Debug', 'delta_V_OF.png');
% secOF.OF_U_FileName = fullfile(resultDir, 'Debug', 'U_OF.png');
% secOF.OF_V_FileName = fullfile(resultDir, 'Debug', 'V_OF.png');
% secOF.OF_disactivated_FileName = fullfile(resultDir, 'Debug', 'disactivated_f.png');
secOF.OF_U_Final_FileName = fullfile(resultDir, 'U_OF_final.png');
secOF.OF_V_Final_FileName = fullfile(resultDir, 'V_OF_final.png');
secOF.OF_Ur_Final_FileName = fullfile(resultDir, 'Ur_OF_final.png');
secOF.OF_Vr_Final_FileName = fullfile(resultDir, 'Vr_OF_final.png');

% Stereo
secST.ST_Il0_FileName = fullfile(resultDir, 'Il0.pgm');
secST.ST_Ir0_FileName = fullfile(resultDir, 'Ir0.pgm');
secST.ST_D0_Init_Middlebury_FileName = fullfile(resultDir, 'D0_ST_init_Middlebury.pgm');
secST.ST_D0_FileName = fullfile(resultDir, 'D0_ST_final.png');

% Scene Flow
% secSF.SF_Dt_Init_FileName = fullfile(resultDir, 'Debug', 'Dt_SF_Init.png');
% secSF.SF_U_increment_FileName = fullfile(resultDir, 'Debug', 'delta_U_SF.png');
% secSF.SF_V_increment_FileName = fullfile(resultDir, 'Debug', 'delta_V_SF.png');
% secSF.SF_D0_increment_FileName = fullfile(resultDir, 'Debug', 'delta_D0_SF.png');
% secSF.SF_Dt_increment_FileName = fullfile(resultDir, 'Debug', 'delta_Dt_SF.png');
% secSF.SF_U_FileName = fullfile(resultDir, 'Debug', 'U_SF.png');
% secSF.SF_V_FileName = fullfile(resultDir, 'Debug', 'V_SF.png');
% secSF.SF_D0_FileName = fullfile(resultDir, 'Debug', 'D0_SF.png');
% secSF.SF_Dt_FileName = fullfile(resultDir, 'Debug', 'Dt_SF.png');
% secSF.SF_disactivated_stereoinit_FileName = fullfile(resultDir, 'Debug', 'disactivated_sti.png');
% secSF.SF_disactivated_stereo_FileName = fullfile(resultDir, 'Debug', 'disactivated_st.png');
% secSF.SF_disactivated_flowleft_FileName = fullfile(resultDir, 'Debug', 'disactivated_fl.png');
% secSF.SF_disactivated_flowright_FileName = fullfile(resultDir, 'Debug', 'disactivated_fr.png');
% secSF.SF_Il0_FileName = fullfile(resultDir, 'Debug', 'Il0.png');
% secSF.SF_Ilt_FileName = fullfile(resultDir, 'Debug', 'Ilt.png');
% secSF.SF_Ir0_FileName = fullfile(resultDir, 'Debug', 'Ir0.png');
% secSF.SF_Irt_FileName = fullfile(resultDir, 'Debug', 'Irt.png');
% secSF.SF_Ilt_Warped_FileName = fullfile(resultDir, 'Debug', 'Ilt_warped.png');
% secSF.SF_Ir0_Warped_FileName = fullfile(resultDir, 'Debug', 'Ir0_warped.png');
% secSF.SF_Irt_Warped_FileName = fullfile(resultDir, 'Debug', 'Irt_warped.png');
secSF.SF_U_Final_FileName = fullfile(resultDir, 'U_SF_final.csv');
secSF.SF_V_Final_FileName = fullfile(resultDir, 'V_SF_final.csv');
secSF.SF_D0_Final_FileName = fullfile(resultDir, 'D0_SF_final.csv');
secSF.SF_Dt_Final_FileName = fullfile(resultDir, 'Dt_SF_final.csv');
secSF.SF_D0_Final_Middlebury_FileName = fullfile(resultDir, 'D0_SF_final_Middlebury.pgm');
secSF.SF_Steps = 1;

% Apply custom parameters
if nargin > 3
    if isfield(customParam, 'secG')
        secG = replaceField(secG, customParam.secG);
    end
    if isfield(customParam, 'secOF')
        secOF = replaceField(secOF, customParam.secOF);
    end
    if isfield(customParam, 'secST')
        secST = replaceField(secST, customParam.secST);
    end
    if isfield(customParam, 'secSF')
        secSF = replaceField(secSF, customParam.secSF);
    end
end

% Create .ini file
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end
if ~exist(fullfile(resultDir, 'Debug'), 'dir')
    mkdir(resultDir, 'Debug');
end
iniFile = fullfile(resultDir,'vsf.ini');
fid = fopen(iniFile, 'w');
fprintf(fid, '[Global]\n');
writeSection(fid, secG);
fprintf(fid, '[OpticFlow]\n');
writeSection(fid, secOF);
fprintf(fid, '[Stereo]\n');
writeSection(fid, secST);
fprintf(fid, '[SceneFlow]\n');
writeSection(fid, secSF);
fclose(fid);

end

function param = replaceField(param, customParam)
fields = fieldnames(customParam);
for i = 1:numel(fields)
    if isfield(param, fields{i})
        param.(fields{i}) = customParam.(fields{i});
    end
end
end

function writeSection(fid, param)
fields = fieldnames(param);
for i = 1:numel(fields)
    if ischar(param.(fields{i}))
        fprintf(fid,'%s = %s\n',fields{i},param.(fields{i}));
    else
        fprintf(fid,'%s = %s\n',fields{i},num2str(param.(fields{i})));
    end
end
fprintf(fid, '\n');
end

