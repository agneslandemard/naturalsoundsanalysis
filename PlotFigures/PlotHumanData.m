%% Compare ferrets and humans, for exp I
% - Difference nat vs synth maps, by category of sounds, for ferret in exp
% I, and humans (Fig 4E)
% - NSE by category of sounds, as function of distance to center of PAC,
% for ferrets in exp I and humans (Fig 4F)
% - Boxplot summary and quantification (Fig S8D)

%% Load ferret data

P.experiment = 'natural';
P.n_ics = NbICs;
P.name = version_name;
D = LoadDenoisedData(P);

% Separate ferrets
A_Idx = find(contains(D.hemis,'A'));
T_Idx = find(contains(D.hemis,'T'));

D.SrecoFullA = D.SrecoFull(:,ismember(D.si,A_Idx),:);
D.SrecoFullT = D.SrecoFull(:,ismember(D.si,T_Idx),:);

%% Load human data
clear cat

load([additional_path 'Human/human_data.mat'],'human_R','stim_names','fullData','surf_test_retest_nse');
xi = [1 5]; % models to consider
R = nanmean(human_R(:,[1 5],:,:),4);
R = reshape(R,size(R,1)*size(R,2),size(R,3));

% get back component weights from data
LH_allSubjectsRep1 = permute(snm(fullData.G.grid_data{2}(:,:,:,[2 1],:,1:2:end),6),[1 2 5 3 4]);
RH_allSubjectsRep1 = permute(snm(fullData.G.grid_data{1}(:,:,:,[2 1],:,1:2:end),6),[1 2 5 3 4]);
[xl, yl, n_subj, n_sds, n_mds] = size(LH_allSubjectsRep1);
[xr, yr, ~, ~, ~] = size(RH_allSubjectsRep1);
Dor_allSubjectsRep1 = cat(1,reshape(LH_allSubjectsRep1,xl*yl,n_subj,n_sds*n_mds),reshape(RH_allSubjectsRep1,xr*yr,n_subj,n_sds*n_mds));
Dor_allSubjectsRep1 = reshape(Dor_allSubjectsRep1,(xl*yl+xr*yr)*n_subj,n_sds*n_mds)';
%( LH + RH )*n_subj
W_allSubjectsRep1 = pinv(R)*Dor_allSubjectsRep1;

LH_allSubjectsRep2 = permute(snm(fullData.G.grid_data{2}(:,:,:,[2 1],:,2:2:end),6),[1 2 5 3 4]);
RH_allSubjectsRep2 = permute(snm(fullData.G.grid_data{1}(:,:,:,[2 1],:,2:2:end),6),[1 2 5 3 4]);
Dor_allSubjectsRep2 = cat(1,reshape(LH_allSubjectsRep2,xl*yl,n_subj,n_sds*n_mds),reshape(RH_allSubjectsRep2,xr*yr,n_subj,n_sds*n_mds));
Dor_allSubjectsRep2 = reshape(Dor_allSubjectsRep2,(xl*yl+xr*yr)*n_subj,n_sds*n_mds)';
%( LH + RH )*n_subj
W_allSubjectsRep2 = pinv(R)*Dor_allSubjectsRep2;

% Project on those weights
Dreco_tmp_allSubjects1 = reshape(R*W_allSubjectsRep1,n_sds, n_mds ,(xl*yl+xr*yr)*n_subj);
Dreco_tmp_allSubjects2 = reshape(R*W_allSubjectsRep2,n_sds, n_mds ,(xl*yl+xr*yr)*n_subj);

D.SrecoFull_human_allSubjects = nan(length(D.param.snd.idxSound(:)),(xl*yl+xr*yr)*n_subj,2);

% Match order of sounds
for md = 1:length(xi)
    for sd = 1:length(D.param.snd.SoundList)
        matchSD = find(strcmp(stim_names,D.param.snd.SoundList(sd)));
        if ~isempty(matchSD)
            D.SrecoFull_human_allSubjects(D.param.snd.idxSound(sd,xi(md)),:,1) = Dreco_tmp_allSubjects1(matchSD,md,:);
            D.SrecoFull_human_allSubjects(D.param.snd.idxSound(sd,xi(md)),:,2) = Dreco_tmp_allSubjects2(matchSD,md,:);
                     
        else
            disp('Sound not found')
        end
    end
end

% remove not reliable voxels
TestRetestH = NSE(D.SrecoFull_human_allSubjects(:,:,1),D.SrecoFull_human_allSubjects(:,:,2),1);
D.SrecoFull_human_allSubjects(:,TestRetestH>0.4,:) = nan;

%% Plot human maps of difference nat vs synth
% fig 4E, right

% get properties for plotting
toPlot = fullData.G;

SoundType = {'Speech','Music','Others','All'};
hemi = 2; %1=RH, 2=LH
plotWhat = 'diff_resp';
% can also view 'nse_across_sounds' which matches the plots of NSE as funciton to
% distance to PAC computed later on

% Compute difference nat vs synth 
[~, diff_all_sds] = NSE_disFromPAC(D.SrecoFull_human_allSubjects,...
    [],D.param,[],SoundType,plotWhat);

% colormap properties
clear Color
m = 2;
Color.ColorAxis = m*[-1 1];
Color.cm = cmap_from_name('cbrewer-blue-red');

for sd = 1:length(SoundType)
    
    SdsIdx = SelectSounds(SoundType{sd},D.param);
    
    % for selected sounds
    diff_all= snm(diff_all_sds(SdsIdx,:),1);
    
    % average across subjects
    diff_all = nanmean(reshape(diff_all,length(diff_all)/n_subj,n_subj),2);
    
    % prepare plotting
    D_LH = reshape(diff_all(1:xl*yl),xl,yl);
    D_RH = reshape(diff_all(xl*yl+1:end),xr,yr);
    
    toPlot.grid_data{1} = D_RH;
    toPlot.grid_data{2} = D_LH;
    
    surftoPlot = grid2surface(toPlot);
    surftoPlot(surf_test_retest_nse(:,1)>0.4,1) = NaN;
    surftoPlot(surf_test_retest_nse(:,2)>0.4,2) = NaN;
    
    plot_fsaverage_1D_overlay(surftoPlot(:,hemi),'lh',Color.cm ,Color.ColorAxis);
    set(gcf,'Position',[100+(sd-1)*300 487 310 303])
    set(gcf, 'Renderer', 'opengl');
    
end


%% Plot ferret maps of difference nat vs synth, for exp I
% fig 4E, row 1

SoundType = {'Ferrets','Speech','Music','Others'};

hemistoplot = 4; %1:D.n_hemis;
figure('Position',[136 172 868 597]);
n_cols = length(SoundType);

plotWhat = 'diff_resp';
% can also view 'nse_across_sounds' which matches the plots of NSE as funciton to
% distance to PAC computed later on

% Compute difference nat vs synth 
[~, diff_all_sds] = NSE_disFromPAC(D.SrecoFull,...
    [],D.param,[],SoundType,plotWhat);

% colormap properties (same colormap as for humans)
clear Color
m = 2;
Color.ColorAxis = m*[-1 1];
Color.cm = cmap_from_name('cbrewer-blue-red');

for sd = 1:length(SoundType)

   SdsIdx = SelectSounds(SoundType{sd},D.param);
    
    % for selected sounds
    diff_all = snm(diff_all_sds(SdsIdx,:),1);
 
    for h = 1 : length(hemistoplot)
        hemi = hemistoplot(h);
        xi = D.si == hemi;
        diff_single_subj = diff_all(xi);
        
        X = load([data_path D.hemis{hemi} D.data_suffix '.mat'],'param');
        Q = X.param;
        
        subplot(length(hemistoplot),n_cols,sd+(h-1)*n_cols)
        PlotTopView(diff_single_subj,Q,Color);
        if h == 1
            title([SoundType{sd}])
        end

    end
end


%% Plot distance from pac for both species
% Fig 4F

plot_type = 'bounded';

% FERRET 1
SoundType = {'Ferrets','Speech','Music','Others'};
dist = load([additional_path 'Coordinates/distances_to_pac_' P.experiment '.mat']);
dist.distance_to_pac = dist.distance_to_pac(ismember(D.si,A_Idx));
nse_by_dis = NSE_disFromPAC(D.SrecoFullA,D.si(ismember(D.si,A_Idx)),D.param,dist,SoundType);

colors = [0.949    0.604	0.722  ;...
        0.1882    0.7490    0.6196 ;...
        0.0549    0.3020    0.5843 ;...
        snm(D.param.plt.SoundColors(SelectSounds('Others',D.param),:),1)];
    
figure('Position',[489 462 500 336]);
axf1 = subplot(1,3,1);
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
yticks(0:0.5:1)
ylim([0 1.4])
xticks(dist.distances(1:2:end)./10)

% FERRET 2
SoundType = {'Ferrets','Speech','Music','Others'};
dist = load([additional_path 'Coordinates/distances_to_pac_' P.experiment '.mat']);

dist.distance_to_pac = dist.distance_to_pac(ismember(D.si,T_Idx));
nse_by_dis = NSE_disFromPAC(D.SrecoFullT,D.si(ismember(D.si,T_Idx))-2,D.param,dist,SoundType);

colors = [0.949    0.604	0.722  ;...
        0.1882    0.7490    0.6196 ;...
        0.0549    0.3020    0.5843 ;...
        snm(D.param.plt.SoundColors(SelectSounds('Others',D.param),:),1)];
    
axf2 = subplot(1,3,2);
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
yticks(0:0.5:1)
ylim([0 1.4])
xticks(dist.distances(1:2:end)./10)

% HUMANS
clear dist
SoundType = {'Speech','Music','Others'};
dist = load([additional_path 'Coordinates/distances_to_pac_human.mat']);
dist.distance_to_pac = dist.distance_to_pac(:);
dist.distances = dist.distances(1:end-2);

SubjIdx = mat2vec(repmat(1:5,[xl*yl+xr*yr 1]));
nse_by_dis = NSE_disFromPAC(D.SrecoFull_human_allSubjects,SubjIdx,D.param,dist,SoundType);

colors = [0.1882    0.7490    0.6196 ;...
    0.0549    0.3020    0.5843 ;...
    snm(D.param.plt.SoundColors(SelectSounds('Others',D.param),:),1)];

axh = subplot(1,3,3);
hold all
for sd = 1:length(SoundType)
    switch plot_type
        case 'normal'
            
            plot(dist.distances+2.5,snm(nse_by_dis(sd,:,:,1),3),'color',colors(sd,:),'LineWIdth',2);
            plot(dist.distances+2.5,squeeze(nse_by_dis(sd,:,:,1)),'color',colors(sd,:),'LineWIdth',0.3);
        case 'bounded'
            boundedline(dist.distances+2.5,snm(nse_by_dis(sd,:,:,1),3),snm(nse_by_dis(sd,:,:,2),3),'cmap',colors(sd,:),'alpha');
    end
end
xlabel('Distance to PAC')
title('Human')
yticks(0:0.5:1)
ylim([0 1.4])
xticks(dist.distances(1:2:end))
linkaxes([axf1 axf2 axh],'y');

%% Boxplot summary and quantification (Fig S8D)

method = 'ncNSE_bycat';
SoundType = {'Ferrets','Speech','Music','Others'};
sbj.Subjects = {'SrecoFullA','SrecoFullT','SrecoFull_human_allSubjects'};
sbj.Subject_list = {'Ferret A','Ferret T','Human'};
 
[corr_pvals_indiv] = PlotNSEBoxplots(D,sbj,SoundType,method);

