function writeCocolibConfig( lffile, customParam )
%WRITECOCOLIBCONFIG Write cocolib config file

if nargin < 2
    customParam = struct();
end

[pathstr, fname, ~] = fileparts(lffile);
lfParam = LFReadMetadata(fullfile(pathstr,[fname,'.json']));
sz = lfParam.camParam.resol;

% config for epi_disparity
epiParam = struct();
epiParam.lf_format = 'Grid';
epiParam.lf_mode = 'rgb';
epiParam.lf_reverse_s = 1;
epiParam.lf_reverse_rowcol = 1;
epiParam.lf_reverse_t = 1;
epiParam.stereo_err_max	= 10.0;

epiParam.lf_path = fullfile(pathstr,[fname '/']);
epiParam.lf_name = 'Sub_%02i_%02i.png';
epiParam.S = sz(2);
epiParam.smin = 1;
epiParam.smax = sz(2);
epiParam.T = sz(1);
epiParam.tmin = 1;
epiParam.tmax = sz(1);

epiParam.W = sz(4);
epiParam.H = sz(3);

epiParam.algorithm = 'epi_disparity';

epiParam.dataset_dx = 'disparity_h';
epiParam.dataset_dy = 'disparity_v';
epiParam.dataset_fx = 'coherence_h';
epiParam.dataset_fy = 'coherence_v';

% number of disparity labels and disparity range
epiParam.dmin = 0;
epiParam.dmax = 2.0;
epiParam.L = 64;

% structure tensor config
epiParam.sigma = 1.0;
epiParam.tau = 0.5;

epiParam.outdir = fullfile(pathstr,fname,'disparity_epi');

if nargin > 1
    fields = fieldnames(customParam);
    for i = 1:numel(fields)
        if isfield(epiParam, fields{i})
            epiParam.(fields{i}) = customParam.(fields{i});
        end
    end
end

writeConfigFile(epiParam,fullfile(pathstr,['disparity_epi_' fname]));

% config for epi_disparity_filter
filterParam = struct();
filterParam.lf_format = 'Grid';
filterParam.lf_mode = 'rgb';
filterParam.lf_reverse_s = 1;
filterParam.lf_reverse_rowcol = 1;
filterParam.lf_reverse_t = 1;
filterParam.stereo_err_max	= 10.0;

filterParam.lf_path = fullfile(pathstr,[fname '/']);
filterParam.lf_name = 'Sub_%02i_%02i.png';
filterParam.S = sz(2);
filterParam.smin = 1;
filterParam.smax = sz(2);
filterParam.T = sz(1);
filterParam.tmin = 1;
filterParam.tmax = sz(1);

filterParam.W = sz(4);
filterParam.H = sz(3);

filterParam.algorithm = 'disparity_epi_filter';

filterParam.all_views = 1;

% number of disparity labels and disparity range
filterParam.dmin = 0;
filterParam.dmax = 2.0;
filterParam.L = 64;

% smoothness factor
filterParam.lambda = 0.75;

% disparity map filter
filterParam.disp_filter = 'tv';
filterParam.disp_filter_iter = 200;

filterParam.outdir = fullfile(pathstr,fname,'disparity_epi_filter');

filterParam.lf_name_channels = fullfile(pathstr,fname,'disparity_epi','disparity_estimate_channels.h5');
filterParam.dataset_dx = 'disparity_h';
filterParam.dataset_dy = 'disparity_v';
filterParam.dataset_fx = 'coherence_h';
filterParam.dataset_fy = 'coherence_v';

if nargin > 1
    fields = fieldnames(customParam);
    for i = 1:numel(fields)
        if isfield(filterParam, fields{i})
            filterParam.(fields{i}) = customParam.(fields{i});
        end
    end
end

writeConfigFile(filterParam,fullfile(pathstr,['disparity_epi_filter_' fname]));

end

function writeConfigFile(param,filename)
fid = fopen(filename,'w');
fields = fieldnames(param);
for i = 1:numel(fields)
    if ischar(param.(fields{i}))
        fprintf(fid,'%s %s\n',fields{i},param.(fields{i}));
    else
        fprintf(fid,'%s %s\n',fields{i},num2str(param.(fields{i})));
    end
end
fclose(fid);
end

