%%% Predict human data based on ferret components
% parameters of analysis
clear I;
I.target = 'human';
% if I.target  = 'human' -> fig S6
% if I.target = 'ferret' -> fig S7

%% Load and format components

clear LP;
LP.experiment = 'natural';
LP.name = version_name;
D = LoadComponents(LP);

% load human components
load([additional_path 'Human/HumanComponentsForPrediction.mat'],'human_R','human_R1','human_R2');

% load ferret components
if exist([analysis_path LP.experiment '_' LP.name '/FerretComponentsForPrediction.mat'],'file')
    load([analysis_path LP.experiment '_' LP.name '/FerretComponentsForPrediction.mat'],'ferret_R','ferret_R1','ferret_R2');
else
    % need to run this script first to format components correctly
    open FormatFerretComponentsForPrediction
end

X = load([data_path D.hemis{1} D.data_suffix],'param');
P = X.param;

% keep only sounds that are common to both species
ns = find(isnan(nanmean(human_R,[2 3])));
SoundstoKeep = setdiff(1:P.snd.Nsound,ns);
ModelstoKeep = [1 5]; % models to consider

switch I.target
    case 'human'
        target_R = human_R(SoundstoKeep,ModelstoKeep,:);
        target_R1 = human_R1(SoundstoKeep,ModelstoKeep,:);
        target_R2 = human_R2(SoundstoKeep,ModelstoKeep,:);
        
        predictor_R = ferret_R(SoundstoKeep,ModelstoKeep,:);
        predictor_R1 = ferret_R1(SoundstoKeep,ModelstoKeep,:);
        predictor_R2 = ferret_R2(SoundstoKeep,ModelstoKeep,:);
        
        CompOrder = 1:size(target_R,3);
    case 'ferret'
        target_R = ferret_R(SoundstoKeep,ModelstoKeep,:);  
        target_R1 = ferret_R1(SoundstoKeep,ModelstoKeep,:);
        target_R2 = ferret_R2(SoundstoKeep,ModelstoKeep,:);
   
        predictor_R = human_R(SoundstoKeep,ModelstoKeep,:);
        predictor_R1 = human_R1(SoundstoKeep,ModelstoKeep,:);
        predictor_R2 = human_R2(SoundstoKeep,ModelstoKeep,:);
        
        CompOrder = 1:size(target_R,3);
        % can be modified
    otherwise
        error('Target must be human or ferret not %s', I.target);
end

%% Quantify variance in target for synthetic and natural-synthetic parts of components response

switch I.target
    case 'human'
        load([additional_path 'Human/human_data.mat'],'human_R','stim_names','fullData');
        human_R = snm(human_R(:,ModelstoKeep,:,:),4);
        LH = snm(fullData.G.grid_data{1}(:,:,:,[2 1],:,:),[5 6]);
        RH = snm(fullData.G.grid_data{2}(:,:,:,[2 1],:,:),[5 6]);
        s = size(LH);
        D = cat(1,reshape(LH,s(1)*s(2),s(3)*s(4)),reshape(LH,s(1)*s(2),s(3)*s(4)))';
        W = pinv(reshape(human_R,s(3)*s(4),size(human_R,3)))*D;
        
        RMS_w = sqrt(nanmean(W.^2,2));

    case 'ferret'
        RMS_w = sqrt(nanmean(D.W_ica_train.^2,2));
        
end

variance = nan(size(target_R,3),2);

target_R1_n = nan(size(target_R1));
target_R2_n = nan(size(target_R2));

for i = 1 : size(target_R,3)
    target_R1_n(:,:,i) = target_R1(:,:,i).*RMS_w(i);
    target_R2_n(:,:,i) = target_R2(:,:,i).*RMS_w(i);
    
    % natural
    variance(i,1) = 0.25*(var(target_R1_n(:,end,i)+target_R2_n(:,end,i))-var(target_R1_n(:,end,i)-target_R2_n(:,end,i)));
    % difference natural -synthetic
    variance(i,2) = 0.25*(var((target_R1_n(:,1,i)-target_R1_n(:,end,i))+(target_R2_n(:,1,i)-target_R2_n(:,end,i)))-var((target_R1_n(:,1,i)-target_R1_n(:,end,i))-(target_R2_n(:,1,i)-target_R2_n(:,end,i))));
 
end


%% Prediction of target components

[n_stims_per_cond, n_conds, n_target_components, ~] = size(target_R);
[~, ~, n_predictor_components,~] = size(predictor_R);
assert(size(target_R,1) == size(predictor_R,1))
assert(size(target_R,2) == size(predictor_R,2))

Y = reshape(target_R, [n_stims_per_cond*n_conds, n_target_components]);
X = reshape(predictor_R, [n_stims_per_cond*n_conds, n_predictor_components]);

% folds used for cross-validation
ResetRandStream2(1);
n_folds = 9;
folds = subdivide(n_stims_per_cond, n_folds)';
folds = repmat(folds, 1, n_conds);
folds = folds(randperm(n_stims_per_cond),:);
folds_unwrap = reshape(folds, [n_stims_per_cond*n_conds, 1]);

% prediction
Yh = regress_predictions_from_3way_crossval(X, Y, ...
    'test_folds', folds_unwrap, 'train_folds', folds_unwrap, ...
    'demean_feats',true , 'std_feats', true, 'method', 'ridge','K', 2.^(-100:100));

% split out conditions
target_prediction = reshape(Yh, [n_stims_per_cond, n_conds, n_target_components]);

%% Plot measured vs predicted component responses, for synthetic and difference natural - synthetic

figh = figure;
n_rows = 2;
n_cols = n_target_components;

set(figh, 'Position', [0, 0, 200*n_cols, 200*n_rows]);
for ic = 1:n_target_components
    
    i = CompOrder(ic);
    
    % Plot measured vs predicted response to synthetic
    subplot(n_rows, n_cols, ic);
    hold on;
    X = cat(2,[target_prediction(:,end,i); target_R(:,end,i)],[target_prediction(:,1,i)-target_prediction(:,end,i);target_R(:,1,i)-target_R(:,end,i)]);
    maxrange = max(diff([min(X); max(X)]));
    bounds = min([target_prediction(:,end,i); target_R(:,end,i)])+[0 maxrange];
    bounds = bounds + [-1 1]*diff(bounds)*0.1;
    plot(bounds, bounds, 'r--','LineWidth', 2);
    scatter(target_prediction(:,end,i), target_R(:,end,i),[],P.plt.SoundColors(SoundstoKeep,:),'filled')
    xlabel('Predicted');
    ylabel('Measured');
    title([upper(I.target(1)) num2str(i)]);
    xlim(bounds);
    ylim(bounds);
    axis equal tight
    if bounds(1)<0
        vline(0,'k--')
        hline(0,'k--')
    end
    
    
    % Plot measured vs predicted difference between natural and synthetic
    subplot(n_rows, n_cols, ic+n_cols);
    hold on;
    bounds=min([target_prediction(:,1,i)-target_prediction(:,end,i); target_R(:,1,i)-target_R(:,end,i)])+[0 maxrange];
    bounds = bounds + [-1 1]*diff(bounds)*0.1;
    plot(bounds, bounds, 'r--','LineWidth', 2);
    scatter(target_prediction(:,1,i)-target_prediction(:,end,i), target_R(:,1,i)-target_R(:,end,i),[],P.plt.SoundColors(SoundstoKeep,:),'filled')
    xlabel('Predicted');
    ylabel('Measured');
    xlim(bounds);
    ylim(bounds);
    axis equal tight
    
    if bounds(1)<0
        vline(0,'k--')
        hline(0,'k--')
    end
end


%% Prediction using independent splits, for quantification
% Fig S6D

splits = {target_R1,target_R2;
    predictor_R1,predictor_R2};

Nbs = 5; % number of bootstraps

[n_stims_per_cond, n_conds, n_target_components, ~] = size(target_R1);
        
ResetRandStream2(1);
nc_nse_all = nan(n_target_components,2,2,Nbs);
nc_corr_all = nan(n_target_components,2,2,Nbs);

for spl = 1:2
    
    target_R1_tmp=splits{1,1};
    predictor_R_tmp=splits{2,spl};
    
    [~, ~, n_predictor_components,~] = size(predictor_R_tmp);
    assert(size(target_R1_tmp,1)==size(predictor_R_tmp,1))
    assert(size(target_R1_tmp,2)==size(predictor_R_tmp,2))
    
    Y = reshape(target_R1_tmp, [n_stims_per_cond*n_conds, n_target_components]);
    X = reshape(predictor_R_tmp, [n_stims_per_cond*n_conds, n_predictor_components]);
    
    for bs=1:Nbs
        % folds used for cross-validation
        
        n_folds = 9;
        folds = subdivide(n_stims_per_cond, n_folds)';
        folds = repmat(folds, 1, n_conds);
        folds = folds(randperm(n_stims_per_cond),:);
        folds_unwrap = reshape(folds, [n_stims_per_cond*n_conds, 1]);
        
        % prediction
        Yh = regress_predictions_from_3way_crossval(X, Y, ...
            'test_folds', folds_unwrap, 'train_folds', folds_unwrap, ...
            'demean_feats', true, 'std_feats', true, 'method','ridge','K', 2.^(-100:100));
        
        % split out conditions
        target_prediction1 = reshape(Yh, [n_stims_per_cond, n_conds, n_target_components]);
        
        % split 2
        target_R2_tmp=splits{1,2};
        predictor_R_tmp=splits{2,setdiff(1:2,spl)};
        
        [n_stims_per_cond, n_conds, n_target_components, ~] = size(target_R2);
        [~, ~, n_predictor_components,~] = size(predictor_R_tmp);
        assert(size(target_R2_tmp,1)==size(predictor_R_tmp,1))
        assert(size(target_R2_tmp,2)==size(predictor_R_tmp,2))
        
        Y = reshape(target_R2_tmp, [n_stims_per_cond*n_conds, n_target_components]);
        X = reshape(predictor_R_tmp, [n_stims_per_cond*n_conds, n_predictor_components]);
        
        % prediction
        Yh = regress_predictions_from_3way_crossval(X, Y, ...
            'test_folds', folds_unwrap, 'train_folds', folds_unwrap, ...
            'demean_feats', true, 'std_feats', true, 'method', 'ridge','K', 2.^(-100:100));
        
        % split out conditions
        target_prediction2 = reshape(Yh, [n_stims_per_cond, n_conds, n_target_components]);
       
        % Noise corrected NSE
        nse_tmp = nan(n_target_components,2);
          % only mm=shared
        nse_tmp(:,1) = NSE_noise_corrected(squeeze(target_prediction1(:,end,:)),squeeze(target_prediction2(:,end,:)),...
            squeeze(target_R1(:,end,:)),squeeze(target_R2(:,end,:)),1);
        
          % orig - mm
        nse_tmp(:,2) = NSE_noise_corrected(squeeze(target_prediction1(:,1,:)-target_prediction1(:,end,:)),squeeze(target_prediction2(:,1,:)-target_prediction2(:,end,:)),...
            squeeze(target_R1(:,1,:)-target_R1(:,end,:)),squeeze(target_R2(:,1,:)-target_R2(:,end,:)),1);
        
        nc_nse_all(:,:,spl,bs)=(1-nse_tmp).^2;
        
        
        % Noise corrected correlation
        corr_tmp=nan(n_target_components,2);
        for ic = 1:n_target_components
            corr_tmp(ic,1) = noise_corrected_correlation(squeeze(target_prediction1(:,end,ic)),squeeze(target_prediction2(:,end,ic)),...
                squeeze(target_R1(:,end,ic)), squeeze(target_R2(:,end,ic)));
            
            corr_tmp(ic,2) = noise_corrected_correlation(squeeze(target_prediction1(:,1,ic)-target_prediction1(:,end,ic)), squeeze(target_prediction2(:,1,ic)-target_prediction2(:,end,ic)),...
                squeeze(target_R1(:,1,ic)-target_R1(:,end,ic)),squeeze(target_R2(:,1,ic)-target_R2(:,end,ic)));
        end
        nc_corr_all(:,:,spl,bs) = corr_tmp.^2;

    end
end

%% Plot summary of variance explained for synthetic and natural-synthetic
% parts of components' responses.

metric = 'nse';
switch metric
    case 'nse'
        metric_id = nc_nse_all;
    case 'corr'
        metric_id = nc_corr_all;
end
 
met = snm(metric_id,[3 4]);
nnvariance = variance./sum(variance(:));
xpos = @(i,j) (i-1)*3+j;
 
error = cat(3,abs(prctile(snm(metric_id,3),2.5,3)-snm(metric_id,[3 4])),prctile(snm(metric_id,3),97.5,3)-snm(metric_id,[3 4]));
figure;
subplot(1,3,1:2)
for i = 1 : n_target_components
  	ic = CompOrder(i);
    hold all
    
    % Variance of synthetic
    b = bar(xpos(i,1),nnvariance(ic,1));
    b.EdgeColor = 'k';
    b.FaceColor = [1 1 1];
    
    % Predicted variance
    c = bar(xpos(i,2),met(ic,1).*nnvariance(ic,1));
    c.FaceColor = [1 1 1]*0.8;
    errorbar(xpos(i,2),met(ic,1).*nnvariance(ic,1),error(ic,1).*nnvariance(ic,1),'k')
    
    % Variance of difference natural - synthetic 
    b = bar(xpos(i,3*(n_target_components+1)),nnvariance(ic,2));
    b.EdgeColor = 'k';
    b.FaceColor = [1 1 1];
    
    % Predicted variance
    d = bar(xpos(i,3*(n_target_components+1)+1),met(ic,2).*nnvariance(ic,2));
    d.FaceColor = [1 1 1]*0.3;
    errorbar(xpos(i,3*(n_target_components+1)+1),met(ic,2).*nnvariance(ic,2),error(ic,2).*nnvariance(ic,2),'k')
    
    
end
sgtitle(['Predict ' I.target ', metric : ' metric])

xticks([xpos(1:n_target_components,1) xpos(1:n_target_components,3*(n_target_components+1))])
xticklabels(repmat(1:n_target_components,1,2))
xlabel([I.target ' components'])
ylabel('Variance (normalized)')


% same but averaged across components
nnvariance = variance./sum(mean(variance,1),2);

xpost = @(i,j) (i-1)*3+j;
error = cat(2,abs(prctile(snm(metric_id,[1 3]),2.5,2)-snm(metric_id,[ 1 3 4])'),prctile(snm(metric_id,[1 3]),97.5,2)-snm(metric_id,[1 3 4])');
subplot(1,3,3)
hold all
b = bar(1,snm(nnvariance(:,1),1));
b.EdgeColor = 'k';
b.FaceColor = [1 1 1];

c = bar(2,snm(met(:,1).*nnvariance(:,1),1));
c.FaceColor = [1 1 1]*0.8;
errorbar(2,snm(met(:,1).*nnvariance(:,1),1),error(1,1).*snm(nnvariance(:,1),1),error(1,2).*snm(nnvariance(:,1),1),'k')%,'horizontal')

b = bar(4,snm(nnvariance(:,2),1));
b.EdgeColor = 'k';
b.FaceColor = [1 1 1];

c = bar(5,snm(met(:,2).*nnvariance(:,2),1));
c.FaceColor = [1 1 1]*0.8;
errorbar(5,snm(met(:,2).*nnvariance(:,2),1),error(2,1).*snm(nnvariance(:,2),1),error(2,2).*snm(nnvariance(:,2),1),'k')%,'horizontal')

ylim([0 1])
xticks([])
ylabel('Variance (normalized)')
xlabel('Average across components')
legend([b,c,d],'Variance of each part of IC','mm variance explained','orig-mm variance explained','Location','NorthOutside'); 
   