%% Plot Figure 3, panels 3 B-E
% Three example components for experiment I

%% Load data 

clear LP;
LP.experiment = 'natural';
LP.n_ics = NbICs;
LP.name = version_name;

D =  LoadComponents(LP);

X = load([data_path D.hemis{1} D.data_suffix '.mat'],'param');
P = X.param;
sound_colors = permute(reshape(P.plt.SoundColors(P.snd.idxSound(:,1:D.n_conds),:), [D.n_stims_per_cond, D.n_conds, 3]), [3, 1, 2]);
symbols = {'o','x', 'square', 'diamond', '+'};

boundsW = quantile(abs(D.W_ica_train'), 0.99)' * [-1 1];

% load cortical model for sounds
CM = load([additional_path 'CorticalModel/CM_natural.mat'],'mods','freqs','freqs_list','specmods_list','tempmods_list','stims','conditions');

%% Define to plot

ICstoPlot = [7 5 3];
hemi = 4;

%% Plot figure 3 B-E

figure;

for i = 1:length(ICstoPlot)
    
    % Select component
    ic = ICstoPlot(i);
    R = D.R_ica_all(:,:,ic);
    
    % Component weights for one hemisphere
    subplot(length(ICstoPlot),5,(i-1)*5+1) 
    xi = D.si == hemi;
    W_single_subj = D.W_ica_train(:,xi);
    X = load([data_path D.hemis{hemi} D.data_suffix '.mat'],'param');
    P = X.param;
    Color.ColorAxis = boundsW(ic,:);
    Color.cm = cmap_from_name('cbrewer-blue-red');
    PlotTopView(W_single_subj(ic,:)',P,Color)
    
    % Responses to natural vs synthetic sounds
    subplot(length(ICstoPlot),5,(i-1)*5+2)
    bounds = [min(R(:)), max(R(:))];
    bounds = bounds + 0.1*[-1 1]*diff(bounds);
    hold on
    for k = 1:D.n_stims_per_cond
        scatter(R(k,1), R(k,D.n_conds),10,sound_colors(:, k, D.n_conds)','filled');
        xlim(bounds); ylim(bounds);
    end
    xlabel('Natural Sounds');
    ylabel('spectrotemporal');
    title(['NSE = ' num2str(NSE(R(:,1), R(:,D.n_conds)),2)])
    plot(bounds, bounds, 'r:', 'LineWidth', 0.3)
    axis equal tight
     
    
    % Test-retest 
    subplot(length(ICstoPlot),5,(i-1)*5+3)
    hold all
    for md = 1:D.n_conds-1:D.n_conds
        cl = sound_colors(:,:,md);
        scatter(mat2vec(D.R_ica_test_even(:,md,ic)),mat2vec(D.R_ica_test_odd(:,md,ic)),10,cl(:,:)',symbols{md},'LineWidth',0.2);
    end
    mi = min([mat2vec(D.R_ica_test_even(:,:,ic));mat2vec(D.R_ica_test_odd(:,:,ic))]);
    ma = max([mat2vec(D.R_ica_test_even(:,:,ic)); mat2vec(D.R_ica_test_odd(:,:,ic))]);
    line([mi ma],[mi ma],'LineStyle',':','Color','k','LineWidth',0.3)
    title(['NSE = ' num2str(NSE(D.R_ica_test_even(:,:,ic),mat2vec(D.R_ica_test_odd(:,:,ic))),2)])
    axis equal tight
    xlabel('test')
    ylabel('retest')
    
    % Correlation with frequency content
    subplot(length(ICstoPlot),5,(i-1)*5+4)
    n_freqbins = length(CM.freqs_list);
    r = corr(CM.freqs,R(:),'rows','pairwise');
    min_freq = 4; % remove lowest frequencies
    plot(r(min_freq:end), 'k-', 'LineWidth', 0.75);
    hold on;
    plot(zeros(n_freqbins-min_freq+1,1), 'r:', 'LineWidth', 0.3);
    bounds = [min(r(min_freq:end)) max(r(min_freq:end))];
    bounds = bounds.*(1/diff(bounds));
    ylim(bounds);
    xlim([1, n_freqbins-min_freq+1]);
    set(gca, 'XTick', 1:2:n_freqbins-min_freq+1, 'XTickLabel', round(CM.freqs_list(min_freq:2:n_freqbins)));
    xlabel('Frequency'); ylabel('Correlation');
    set(gca,'DataAspectRatio',[4 1 1])
    
    % Correlation with modulation content
    subplot(length(ICstoPlot),5,(i-1)*5+5)
    n_spec = length(CM.specmods_list);
    n_temp = length(CM.tempmods_list);
    r = corr(CM.mods(:,:),R(:),'rows','pairwise');
    r = reshape(r,n_spec,n_temp);
    corr_scale = cmap_from_name('cbrewer-blue-red');
    imagesc(flipud(r), 0.7*[-1, 1]);
    axis equal tight
    colormap(corr_scale);
    flip_yvals = fliplr(CM.specmods_list);
    set(gca, ...
        'YTick', 2:2:n_spec, ...
        'XTick', 2:2:n_temp,...
        'YTickLabel', flip_yvals(2:2:n_spec), ...
    'XTickLabel', CM.tempmods_list(2:2:n_temp));
    xlabel('Temp (Hz)'); ylabel('Spec (cyc/oct)');
 
    
end

