%%%  Plot figure 4 panels D, E and G (ferret data, experiment II)

% Load data from experiment II
P.experiment = 'vocalization';
P.n_ics = 8;
P.name = 'myVersion';

D = LoadDenoisedData(P);

A_Idx = find(contains(D.hemis,'A'));
T_Idx = find(contains(D.hemis,'T'));

D.SrecoA = D.Sreco(:,ismember(D.si,A_Idx));
D.SrecoT = D.Sreco(:,ismember(D.si,T_Idx));

D.SrecoFullA = D.SrecoFull(:,ismember(D.si,A_Idx),:);
D.SrecoFullT = D.SrecoFull(:,ismember(D.si,T_Idx),:);

load([additional_path '/Coordinates/distances_to_pac_' P.experiment '.mat'])

%% plot NSE Map (panel D)

hemis_to_plot = 1:D.n_hemis;

n_conds = D.param.snd.Nmodels; % will compare natural sounds to this model (spectrotemporal)

% define color scale
clear Color
Color.ColorAxis = [0 1];
Color.cm = cmap_from_name('lightblue-to-yellow1');

figure('Position',[136 354 1305 415]);
n_cols = 4;
nsm = NSE_noise_corrected(D.SrecoFull(D.param.snd.idxSound(:,1),:,1),D.SrecoFull(D.param.snd.idxSound(:,1),:,2),...
    D.SrecoFull(D.param.snd.idxSound(:,n_conds),:,1),D.SrecoFull(D.param.snd.idxSound(:,n_conds),:,2),1);

for h = 1: length(hemis_to_plot)
    hemi = hemis_to_plot(h);
    
    xi = D.si == hemi;
    nse_single_subj = nsm(xi);
    X = load([data_path D.hemis{hemi} D.data_suffix '.mat'],'param');
    Q = X.param;
   
    subplot(length(hemis_to_plot),n_cols,1+(h-1)*n_cols)
    PlotTopView(nse_single_subj,Q,Color);
    title([D.hemis{hemi} ', mean= ' num2str(mean(nse_single_subj),2) ])
end

 %% Plot difference maps (panel E)

SoundType = {'Ferrets','Speech','Music'};
m = 2;
clear Color
Color.ColorAxis = [-m m];
Color.cm = cmap_from_name('cbrewer-blue-red');

denom_by_pix = nanstd(D.Sreco(D.param.snd.idxSound(:,[1 n_conds]),:),[],1);
diff_all = D.Sreco(D.param.snd.idxSound(:,1),:,1) - D.Sreco(D.param.snd.idxSound(:,n_conds),:,1);
    
for sd = 1:length(SoundType)

    SdsIdx = SelectSounds(SoundType{sd},D.param);
    
    for h = 1: length(hemis_to_plot)
        hemi = hemis_to_plot(h);
        
        % pick out one subject
        xi = D.si == hemi;
        diff_single_subj = nanmean(diff_all(SdsIdx,xi),1)./denom_by_pix(xi);
        
        X = load([data_path D.hemis{hemi} D.data_suffix '.mat'],'param');
        Q = X.param;
        
        subplot(length(hemis_to_plot),n_cols,1+sd+(h-1)*n_cols)
        PlotTopView(diff_single_subj,Q,Color);
        title([SoundType{sd} ' ' D.hemis{hemi}])
    end
end

%% Plot NSE as function of distance to PAC, by category (panel G)

plot_type = 'bounded';
SoundType = {'Ferrets','Speech','Music'};

% Ferret A
dist = load([additional_path '/Coordinates/distances_to_pac_' P.experiment '.mat']);
dist.distance_to_pac = dist.distance_to_pac(ismember(D.si,A_Idx));
nse_by_dis = NSE_disFromPAC(D.SrecoFullA,D.si(ismember(D.si,A_Idx)),D.param,dist,SoundType);

colors = [0.949    0.604	0.722  ;...
        0.1882    0.7490    0.6196 ;...
        0.0549    0.3020    0.5843 ];
    
figure('Position',[489 462 500 336]);
axf1 = subplot(1,2,1);
hold all
for sd = 1:length(SoundType)
   switch plot_type
        case 'normal'
            plot((dist.distances+2.5)./10,snm(nse_by_dis(sd,:,:,1),3),'color',colors(sd,:),'LineWIdth',2);
            plot((dist.distances+2.5)./10,squeeze(nse_by_dis(sd,:,:,1)),'color',colors(sd,:),'LineWIdth',0.3);
        case 'bounded'
            boundedline((dist.distances+2.5)./10,snm(nse_by_dis(sd,:,:,1),3),snm(nse_by_dis(sd,:,:,2),3),'cmap',colors(sd,:),'alpha');
   end
end
xlabel('Distance to PAC')
title('Ferret A')
yticks(0:0.2:0.4)
ylim([0 0.4])

yticks(0:0.5:1)
ylim([0 1])
xticks(dist.distances(1:2:end)./10)

% Ferret T
dist = load([additional_path '/Coordinates/distances_to_pac_' P.experiment '.mat']);
dist.distance_to_pac = dist.distance_to_pac(ismember(D.si,T_Idx));
nse_by_dis = NSE_disFromPAC(D.SrecoFullT,D.si(ismember(D.si,T_Idx))-1,D.param,dist,SoundType);

colors = [0.949    0.604	0.722  ;...
        0.1882    0.7490    0.6196 ;...
        0.0549    0.3020    0.5843 ];

axf2 = subplot(1,2,2);
hold all
for sd = 1:length(SoundType)
   switch plot_type
        case 'normal'
            plot((dist.distances+2.5)./10,snm(nse_by_dis(sd,:,:,1),3),'color',colors(sd,:),'LineWIdth',2);
            plot((dist.distances+2.5)./10,squeeze(nse_by_dis(sd,:,:,1)),'color',colors(sd,:),'LineWIdth',0.3);
        case 'bounded'
            boundedline((dist.distances+2.5)./10,snm(nse_by_dis(sd,:,:,1),3),snm(nse_by_dis(sd,:,:,2),3),'cmap',colors(sd,:),'alpha');
   end
end
xlabel('Distance to PAC')
title('Ferret T')
yticks(0:0.2:0.4)
ylim([0 0.4])

yticks(0:0.5:1)
ylim([0 1])
xticks(dist.distances(1:2:end)./10)
 