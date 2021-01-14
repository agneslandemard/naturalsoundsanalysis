function data = LoadDenoisedData(load_params)

global analysis_path data_path

% Define default parameters
if ~isfield(load_params,'n_ics')
    load_params.n_ics = 8;
end
if ~isfield(load_params,'name')
    load_params.name = 'myVersion';
end

proc_path = [analysis_path load_params.experiment '_' load_params.name ];

% Load hemispheres names
data = load([proc_path '/D_norm.mat'],'hemis');

% Load and recompose data based on runs
D = load([proc_path '/D_denoised_K' num2str(load_params.n_ics) '.mat']);
data.SrecoFullTC = D.D_denoised;
data.SrecoTC = snm(D.D_denoised,4);

% Define other parameters
data.experiment = load_params.experiment;
data.n_hemis = length(data.hemis);
data.data_suffix = '_Anat&Param';
X = load([data_path  data.hemis{1} data.data_suffix '.mat'],'param');
data.param = X.param; clear X
data.data_path = proc_path;
data.si = D.si; % hemisphere identifier

% Define response window and time-averaged response
if ~isfield(load_params,'RespWindow')
    load_params.RespWindow = data.param.proc.MeanInterval;
end
data.Sreco = snm(data.SrecoTC(data.param.exp.PreStimSilence+load_params.RespWindow,:,:),1);
data.SrecoFull = snm(data.SrecoFullTC(data.param.exp.PreStimSilence+load_params.RespWindow,:,:,:),1);
data.RespWindow = load_params.RespWindow;
end
