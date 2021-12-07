%% Load data

% % specify which organoids to analyse
ORGANOIDS = {"L1", 0, 0, 0}; % only test data is included
TREATMENTS = {"C", "D50"};

% make empty cell arrays for results
spike_mat = cell(length(TREATMENTS),length(ORGANOIDS));
spike_times = cell(length(TREATMENTS),length(ORGANOIDS));
burst_edges = cell(length(TREATMENTS),length(ORGANOIDS));
pop_act = cell(length(TREATMENTS),length(ORGANOIDS));
pop_act_times = cell(length(TREATMENTS),length(ORGANOIDS));

% for each array
for organoid = 1:length(ORGANOIDS)
    
    if isstring(ORGANOIDS{organoid})
        
        for treatment = 1:length(TREATMENTS)
            
            % % load data
            
            % load spike data
            t_spk_mat_cont = load(sprintf("./data/%s_t_spk_mat_sorted.mat", TREATMENTS{treatment}));
            spike_mat{treatment,organoid} = t_spk_mat_cont.t_spk_mat;
            spike_times{treatment,organoid} = t_spk_mat_cont.spike_times;
            
            % load burst edges
            edge_cont = load(sprintf("./data/%s_BurstEdgesFull.mat", TREATMENTS{treatment}));
            burst_edges{treatment,organoid} = edge_cont.edges;
            
            % load population activity data and store
            pop_act_cont = load(sprintf("./data/processed/%s_population_activity.mat", TREATMENTS{treatment}));
            pop_act{treatment,organoid} = pop_act_cont.pop_act;
            pop_act_times{treatment,organoid} = pop_act_cont.pop_act_times;
            
        end % treatment
        
    end % if
    
end % organoid

%% Plot figure 4

% define plot parameters
ARRAY_OI = 1; % specify array of interest for plot (panel A, B, D) 

START = -200; % specify plot range for panel A
END = 500; % specify plot range for panel A
SPIKE_MAT_RANGE = [-100, 350]; % specify region indicated with blue dotted lines in A used to calculate B
STEP = 10; % specify step size for panel D

% make lists of start and end values for panel D
start_vals = linspace(START,0,(-START/STEP)+1);
end_vals = linspace(0,END,(END/STEP)+1);

% make save_path
save_path = "./figures/Fig_4";

% plot results (include save_path variable to save results)
plot_fig4(spike_mat, burst_edges, pop_act, pop_act_times, ARRAY_OI, ...
    start_vals, end_vals, SPIKE_MAT_RANGE)%, save_path)

%% Functions

function plot_fig4(spike_mat, burst_edges, pop_act, pop_act_times, array_oi, ...
    start_vals, end_vals, spike_mat_range, save_path)

global organoid_names

% initiate figure
cm = figure(1);
clf

% adjust size of figure
set(gcf,'PaperPositionMode','auto')
set(cm, 'Position', [0 20 1800 600])
set(cm, 'Renderer', 'painters')

% make empty result array
burst_pop_vecs_plot = cell(1,size(spike_mat,1));
dissim_mat = cell(1,size(spike_mat,1));
frac_diff_mat = cell(1,size(spike_mat,2));
av_window_score = zeros(1,size(spike_mat,2));
av_dissim_score = zeros(size(spike_mat));
SEM_dissim_score = zeros(size(spike_mat));

% % pre process data

% for each array
for array = 1:size(spike_mat, 2)
    
    if ~isempty(spike_mat{1,array})
        
        % make empty cell array
        av_pop_vec_diff = cell(1,size(spike_mat, 1));
        
        % for each treatment
        for treatment = 1:size(spike_mat, 1)
            
            % make empty result array
            av_pop_vec_diff{treatment} = NaN(length(start_vals), length(end_vals));
            
            % for each start val
            for s_val = 1:length(start_vals)
                
                % for each end val
                for e_val = 1:length(end_vals)
                    
                    % determine analysis edges
                    analysis_edges = obtain_edges(burst_edges{treatment,array}, ...
                        [start_vals(s_val), end_vals(e_val)], pop_act{treatment,array}, ...
                        pop_act_times{treatment,array});
                    
                    % make empty result array for burst pop vectors
                    burst_pop_vecs = zeros(size(analysis_edges, 1), ...
                        end_vals(e_val)-start_vals(s_val)+1);
                    
                    % sum spiking activity over all channels
                    sum_spike = sum(spike_mat{treatment, array}, 2);
                    smooth_spike = smoothdata(sum_spike, "movmean", 5);
                    
                    % smooth summed spiking activity
                    pop_vec = smooth_spike.*1000;
                    
                    % for each burst
                    for burst = 1:size(analysis_edges, 1)
                        
                        burst_pop_vec = pop_vec(analysis_edges(burst,1):...
                            analysis_edges(burst,2));
                        
                        % cut population vection for specified time
                        burst_pop_vecs(burst, :) = burst_pop_vec;
                        
                    end % burst
                    
                    % store burst pop vecs to be plotted
                    if array == array_oi && s_val == 1 && e_val == length(end_vals)
                        burst_pop_vecs_plot{treatment} = burst_pop_vecs;
                    end % if
                    
                    % % compute burst to burst similarity
                    
                    % make empty result matrix
                    pop_vec_diff = NaN(size(burst_pop_vecs,1),size(burst_pop_vecs,1));
                    
                    % make array with comp_bursts
                    comp_bursts = 1:size(burst_pop_vecs, 1);
                    
                    % for each ref burst
                    for ref_burst = 1:size(burst_pop_vecs, 1)
                        
                        % remove ref_burst from comp bursts
                        comp_bursts(comp_bursts == ref_burst) = [];
                        
                        % for each comp_burst
                        for comp_burst = comp_bursts
                            
                            % subtract comp burst pop vec from ref burst
                            burst_diff = burst_pop_vecs(ref_burst,:) - burst_pop_vecs(comp_burst,:);
                            
                            % compute difference and store
                            pop_vec_diff(comp_burst, ref_burst) = mean(abs(burst_diff));
                            
                        end % comp_burst
                        
                    end % ref_burst
                    
                    av_pop_vec_diff{treatment}(s_val, e_val) = mean(pop_vec_diff, "all", "omitnan");
                    
                    % if results for this time window should be plotted
                    if start_vals(s_val) == spike_mat_range(1) && ...
                            end_vals(e_val) == spike_mat_range(2)
                        
                        % store average and SEM of dissim score
                        av_dissim_score(treatment, array) = ...
                            mean(pop_vec_diff, "all", "omitnan");
                        SEM_dissim_score(treatment, array) = ...
                            std(pop_vec_diff, [], "all", "omitnan")/...
                            sqrt(sum(~isnan(pop_vec_diff),"all"));
                        
                        % if this is array of interest
                        if array == array_oi
                            dissim_mat{treatment} = pop_vec_diff;
                        end % if
                        
                    end % if
                    
                end % e_val
            end % s_val
        end % treatment
        
        % compute fractional difference
        frac_diff_mat{array} = (av_pop_vec_diff{1}-av_pop_vec_diff{2}) ./ ...
            (av_pop_vec_diff{1}+av_pop_vec_diff{2});
        
        % take average of window score
        av_window_score(array) = mean(frac_diff_mat{array}, "all", "omitnan");
        
        % remove NaN values for plotting
        frac_diff_mat{array}(isnan(frac_diff_mat{array})) = 0;
        
    end % if
    
end % array

% % plot results of population vectors for array_oi

% make subplot
subplot(2,4,1)

% plot population vectors per burst
plot_burst_pop_vecs(burst_pop_vecs_plot{1}, [start_vals(1), end_vals(end)], ...
    spike_mat_range)
title("Control", "FontSize", 16)

% make subplot
subplot(2,4,2)

% plot population vectors per burst
plot_burst_pop_vecs(burst_pop_vecs_plot{2}, [start_vals(1), end_vals(end)], ...
    spike_mat_range)
title("Diazepam", "FontSize", 16)

% % plot results of dissimilarity matrix for array_oi

% determine maximum dissimilarity score
max_dissim = max([max(dissim_mat{1}, [], "all"), max(dissim_mat{2}, [], "all")]);

% make subplot
subplot(2,8,[5,6])

% plot dissimilarity matrix
cb = plot_dissim_mat(dissim_mat{1}, max_dissim);
cb.Visible = "off";
title("Control", "FontSize", 16)

% make subplot
subplot(2,8,[7,8])

% plot dissimilarity matrix
cb = plot_dissim_mat(dissim_mat{2}, max_dissim);
title("Diazepam", "FontSize", 16)

% % plot results for all arrays

% make subplot
subplot(2,4,5)

plot_all_array_bar(av_dissim_score', SEM_dissim_score')

% % plot results for all time windows

% make subplot
sph = subplot(2,16,[21:27]);

plot_window_range_heatmap(frac_diff_mat{array_oi}, start_vals, end_vals, ...
    spike_mat_range, sph)

% % plot time window results over all arrays

% make subplot
subplot(2,4,8)

plot_window_array_bar(av_window_score)


% if save_path is included in function call
if exist("save_path", "var")

    "saving"
    
    % save figure
    print(save_path, "-dpng", "-r950")
    
end % if

    % % nested functions
    function plot_window_array_bar(av_window_score)
        
        bar(av_window_score, "b")
        
        xticklabels(organoid_names)
        ylabel("(C-D)/(C+D)")
        
        ax = gca;
        ax.FontSize = 14;
  
    
    end % fun plot_window_array_bar

    % % %
    
    function plot_window_range_heatmap(frac_diff_mat, start_vals, end_vals, ...
            spike_mat_range, sph)
        
        hold on
        
        % plot results
        imagesc(end_vals, start_vals, frac_diff_mat)
        
        % mark example on map
        scatter(spike_mat_range(2), spike_mat_range(1), 36, "k", "o", "filled")
        
        axis ij
        
        colormap(sph, flipud(redblue))
        caxis([-(max(abs(frac_diff_mat), [], "all")*10)/10, (max(abs(frac_diff_mat), [], "all")*10)/10])
        cb = colorbar;
        
        ylabel(cb, "(C-D)/(C+D)", "FontSize", 14)
        ylabel("Start time (ms)")
        xlabel("End time (ms)")
        
        xlim([end_vals(1), end_vals(end)])
        ylim([start_vals(1), start_vals(end)])
        
        ax = gca;
        ax.FontSize = 14;
        
    end % fun plot_window_range_heatmap

    % % %
    
    function plot_all_array_bar(av_dissim_score, err_bars)
        
        hb = bar(av_dissim_score);
        
        hold on
        
        % Find the number of groups and the number of bars in each group
        [ngroups, nbars] = size(av_dissim_score);
        
        % Calculate the width for each bar group
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        
        % Set the position of each error bar in the centre of the main bar
        for i = 1:nbars
            
            % Calculate center of each bar
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x, av_dissim_score(:,i), err_bars(:,i), 'k', 'linestyle', 'none');
            
        end
        
        hold off
        
        set(hb(1),'FaceColor',[0,0,1])
        set(hb(2),'FaceColor',[1,0,0])
        
        legend([hb(1), hb(2)], {"Control", "Diazepam"}, "Location", "SouthOutside")
        
        xticklabels(organoid_names)
        ylabel(["Pop vec diff (ms^{-1})", "(mean +- SEM)"])
        
        ax = gca;
        ax.FontSize = 14;
        
    end % fun plot_all_array_bar

    % % %
    
    function cb = plot_dissim_mat(dissim_mat, max_dissim)
        
        imagesc(dissim_mat)
        caxis([0,ceil(max_dissim/100)*100])
        
        colormap hot
        cb = colorbar;
        ylabel(cb, "Pop vec diff (ms^{-1})", "FontSize", 14)
        
        xlabel("Burst")
        ylabel("Burst")
        
        ax = gca;
        ax.FontSize = 14;
        
    end % fun plot_dissim_mat

    % % %
    
    function plot_burst_pop_vecs(burst_pop_vecs, plot_range, example_markers)
        
        hold on
                
        % for each burst
        for burst = 1:size(burst_pop_vecs,1)
            
            % plot pop vec
            plot([plot_range(1):plot_range(2)], burst_pop_vecs(burst,:), "Color", [0,0,0, 0.05])
            
        end % burst
        
        xline(example_markers(1), "b--", "LineWidth", 2);
        xline(example_markers(2), "b--", "LineWidth", 2);
        
        ylabel("Pop. firing rate (Hz)")
        
        % compute sem over time
        pop_vec_sem = std(burst_pop_vecs)/sqrt(size(burst_pop_vecs,1));
                
        yyaxis right
        plot([plot_range(1):plot_range(2)], pop_vec_sem, "r--", "LineWidth", 2)
        ax = gca;
        ax.YColor = 'r';
        ylabel("SEM over bursts")
        ylim([0,150])
        
        xlabel("Time (ms)")
        xlim([plot_range(1), plot_range(2)])

        ax.FontSize = 14;
        
    end % fun plot_burst_pop_vecs

end % fun pop_vec_clus

%% other functions

function edges = obtain_edges(burst_edges, window, pop_act, pop_act_times)
 
% make empty result array
edges = zeros(size(burst_edges));

% for each burst
for burst = 1:size(burst_edges,1)
    
    % find min and max index value for MUA data
    [~,MUA_min_cut] = min(abs(pop_act_times-burst_edges(burst,1)));
    [~,MUA_max_cut] = min(abs(pop_act_times-burst_edges(burst,2)));
    
    % find index for MUA peak
    [~,MUA_peak_index] = max(pop_act(MUA_min_cut:MUA_max_cut));
    
    % obtain timepoint of MUA peak
    peak_time = pop_act_times(MUA_peak_index+MUA_min_cut)*1000;
    
    % specify cut region
    min_cut = peak_time + window(1);
    max_cut = peak_time + window(2);
    
    edges(burst, :) = [min_cut, max_cut];
    
end % burst

end % fun obtain_edges

% % %

function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
%   Adam Auton, 9th October 2009
if nargin < 1, m = size(get(gcf,'colormap'),1); end
if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end
c = [r g b]; 
end % fun redblue


