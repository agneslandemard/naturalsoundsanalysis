%% Compute NSE maps and NSE as function of distance to PAC (Fig 2E&F)

% Load data
P.experiment = 'natural';
P.n_ics = NbICs;
P.name = version_name;
D = LoadDenoisedData(P);

% 2E : maps
% define color scale
clear Color
Color.ColorAxis = [0 1];
Color.cm = cmap_from_name('lightblue-to-yellow1');

% load distances to pac
load([additional_path '/Coordinates/distances_to_pac_' P.experiment '.mat'])
load([additional_path '/Coordinates/annuli-both-lowfreq.mat']);

figure;
n_conds = D.param.snd.Nmodels;
local_nse = nan(n_conds,length(distances),D.n_hemis);

hemi = 4;
for md = 2 : n_conds
    
    % NSE between selected model and natural 
    nsm = NSE_noise_corrected(D.SrecoFull(D.param.snd.idxSound(:,1),:,1),D.SrecoFull(D.param.snd.idxSound(:,1),:,2),...
        D.SrecoFull(D.param.snd.idxSound(:,md),:,1),D.SrecoFull(D.param.snd.idxSound(:,md),:,2),1);
    
    for j = 1 : D.n_hemis
        % pick out one subject
        xi = D.si == j;
        nse_single_subj = nsm(xi);
      
       for d = 1 : length(distances)
            local_nse(md,d,j) = nanmedian(nse_single_subj(distance_to_pac(xi)==distances(d)));
        end
        
    end

    xi = D.si == hemi;
    nse_single_subj = nsm(xi);
    X = load([data_path  D.hemis{hemi} D.data_suffix '.mat'],'param');
    Q = X.param;
   
    subplot(2,n_conds-1,md-1)
    PlotTopView(nse_single_subj,Q,Color);
    title(regexprep(D.param.snd.models{md},'synth_',''))

    subplot(2,n_conds-1,md-1+(n_conds-1))
    hold all
    plot(distances(2:end)./10,squeeze(local_nse(md,1:end-1,:)),'color',0.6.*[1 1 1])
    plot(distances(2:end)./10,snm(local_nse(md,1:end-1,:),3),'k','LineWidth',2)
    
    if md == 2
        xlabel('Distance from center of PAC (mm)')
        ylabel('NSE (noise-corrected)')
    end
    
    ylim([0 1])
    yticks(0:0.2:1)
    xticks(distances(1:2:end-1)./10)
    xlim([0 distances(end-1)./10])
    if md == 5
        plot(distances(2:end-1)./10,annuli_stat_allsubjects,'color',[232 191 88]./255)
        plot(distances(2:end-1)./10,snm(annuli_stat_allsubjects,1),'LineWidth',2,'color',[173 133 43]./255)
    end
end

%% Plot top view of distance bins

clear Color
Color.ColorAxis = distances;
Color.cm = cbrewer('div','Spectral',length(distances),'pchip');
Color.cm = Color.cm(end:-1:1,:);
figure;
PlotTopView(distance_to_pac(xi),Q,Color);
cbh = colorbar;
cbh.Ticks = linspace(0,1,length(distances));
cbh.TickLabels = distances;
