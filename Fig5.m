%% Load data and and set parameters

% % This script requires the following toolbox:                         % %
% % http://strategic.mit.edu/downloads.php?page=matlab_networks         % %


TREATMENT = "C"; % choose C = control or D50 = diazepam
DIRECTIONAL = true; % specify whether directional adjacency matrix should be generated

% column channel receives signals from row channels %
% row channel sends signals to column channels %

dT = 0.02; % specify delta T in s
MIN_SPIKES = 5; % specify minimum number of spikes per unit

MIN_COEF_THRESH = 0.35; % specify STTC coefficient threshold
MEAN_THRESH = 0.002; % specify threshold for mean spike time latency (in s)
PVAL_THRESH = 0.1; % specify threshold for assigning latency distributions multimodality
FWHM_THRESH = 0.015; % specify threshold for FWHM of latency distribution (in s)

NBOOT = 500; % specify number of boot strap samples for Hartigans' dip test 

% % load data
t_spk_mat_cont = load(sprintf("./data/%s_t_spk_mat_sorted.mat", TREATMENT));
spike_mat = t_spk_mat_cont.t_spk_mat;
coordinates = t_spk_mat_cont.xy_raw;
spike_times = t_spk_mat_cont.spike_times;

% % determine recording time range and channels of interest
RANGE_OI = [0, size(spike_mat, 1)/1000]; % set start and end of recording in s
CHANNELS_OI = [1:size(spike_mat,2)]; % specify channels to be analysed in vector    

%% Processes the data

% make empty result arrays (dim1 = comp_channels, dim2 = ref_channels)
tile_coefficients_whole = NaN(size(spike_mat,2), size(spike_mat,2));
mean_latencies = NaN(size(spike_mat,2), size(spike_mat,2));
pvals_per_edge = NaN(size(spike_mat,2), size(spike_mat,2));
fwhm_per_edge = NaN(size(spike_mat,2), size(spike_mat,2));

% make empty result array for number of spikes per channel
channel_numSpikes = zeros(1,size(spike_mat,2));

% % run tileCoef analysis between individual channels for whole time range

% make vector with all comparison channels
comp_channels = CHANNELS_OI;

% for each reference channel in channels of interest
for ref_channel = CHANNELS_OI
    
    % remove reference channel from comparison channels array
    comp_channels(comp_channels == ref_channel) = [];
    
    % obtain spikes within range of interest for reference channel
    ref_spike_mask = spike_times{ref_channel} >= RANGE_OI(1) & ...
        spike_times{ref_channel} < RANGE_OI(2);
    ref_spikes = spike_times{ref_channel}(ref_spike_mask);
    
    % obtain number of spikes within burst range for reference channel
    N_ref_spikes = length(ref_spikes);
    
    % for each comparison channel
    for comp_channel = comp_channels
        
        % obtain spikes within range of interest for reference channel
        comp_spike_mask = spike_times{comp_channel} >= RANGE_OI(1) & ...
            spike_times{comp_channel} < RANGE_OI(2);
        comp_spikes = spike_times{comp_channel}(comp_spike_mask);
        
        % obtain number of spikes within burst range for reference channel
        N_comp_spikes = length(comp_spikes);
        
        if N_ref_spikes >= MIN_SPIKES && N_comp_spikes >= MIN_SPIKES
            
            % calculate the spike timing tiling coefficient between the
            % reference and comparison channel
            sttc_val = sttc(N_ref_spikes, N_comp_spikes, dT, ...
                RANGE_OI, ref_spikes, comp_spikes);
            
            % if sttc is larger than threshold
            if sttc_val > MIN_COEF_THRESH
                
                % obtain latencies between spike trains
                latencies = extract_latencies(ref_spikes, comp_spikes);
                
                % remove latencies larger than dT
                latencies(abs(latencies) > dT) = [];
                
                if abs(mean(latencies)) > MEAN_THRESH
                    
                    % sort time differences
                    pdf = sort(latencies);
                    
                    % perform Hartigans test for multimodality
                    [dip, p_value, xlow, xup] = HartigansDipSignifTest(pdf, NBOOT);
                    
                    if p_value > PVAL_THRESH
                        
                        % fit normal distribution to time differences
                        pd = fitdist(latencies',"Normal");
                        
                        % extract mean and standard deviation from fit
                        std = pd.sigma;
                        
                        % compute FWHM
                        fwhm = 2*sqrt(2*log(2))*std;
                        
                        if fwhm < FWHM_THRESH
                            
                            % store mean latency
                            mean_latencies(comp_channel, ref_channel) = mean(latencies);
                            
                            % store p_value
                            pvals_per_edge(comp_channel, ref_channel) = p_value;
                            
                            % store fwhm
                            fwhm_per_edge(comp_channel, ref_channel) = fwhm;
                            
                            if DIRECTIONAL == true
                                
                                if mean_latencies(comp_channel, ref_channel) > MEAN_THRESH
                                    
                                    % store the sttc value
                                    tile_coefficients_whole(ref_channel, comp_channel) = sttc_val;
                                    
                                elseif mean_latencies(comp_channel, ref_channel) < -MEAN_THRESH
                                    
                                    % store the sttc value
                                    tile_coefficients_whole(comp_channel, ref_channel) = sttc_val;
                                    
                                end % if
                                
                            else
                                
                                if abs(mean_latencies(comp_channel, ref_channel)) > MEAN_THRESH
                                    
                                    % store the sttc value
                                    tile_coefficients_whole(comp_channel, ref_channel) = sttc_val;
                                    
                                end % if
                            end % if
                        end % if
                    end % if
                end % if
            end % if
        end % if
        
    end % comp_channel
end % ref_channel

%% Save results

if DIRECTIONAL == true
    
    % set file name
    var_save_path = sprintf("./data/processed/%s_tile_coef_mat_whole_dir_dT_%.0f_minSpikes_%.0f", ...
        TREATMENT, dT*1000, MIN_SPIKES);
  
else
    
    % set file name
    var_save_path = sprintf("./data/processed/%s_tile_coef_mat_whole_dT_%.0f_minSpikes_%.0f", ...
        TREATMENT, dT*1000, MIN_SPIKES);

end % if

% save results
save(var_save_path, "tile_coefficients_whole", "mean_latencies", ...
    "pvals_per_edge", "fwhm_per_edge", "dT", "MIN_SPIKES", "MIN_COEF_THRESH", ...
    "MEAN_THRESH", "PVAL_THRESH", "FWHM_THRESH")

%% Plot figure 5 A-C

% set plotting parameters
REF_OI = 26; % specify reference unit number
COMP_OI = 32; % specify comparison unit number

BIN_WIDTH = 0.003; % specify bin width for plotting (in s)
PLOT_WIDTH = 0.03; % specify number of seconds before and after 0 to show in plot

SPIKE_TRAIN_RANGE = [0,120]; % specify time range for which to plot spike train (s)

ref_spikes = spike_times{REF_OI};
comp_spikes = spike_times{COMP_OI};

MIN_dT = 0.0005; % specify lowest delta t value to use for analysis (in s)
MAX_dT = 0.150; % specify highest delta t value to use for analysis (in s)
STEP_dT = 0.0005; % specify stepsize in delta t range (in s) 
 
% set delta t range
dT_range = linspace(MIN_dT, MAX_dT, 1+(MAX_dT-MIN_dT)/STEP_dT);
sttc_range = zeros(1,length(dT_range));

% compute STTC for delta t range
for i = 1:length(dT_range)
    sttc_range(i) = sttc(length(ref_spikes), length(comp_spikes), dT_range(i), ...
                [0, size(spike_mat, 1)/1000], ref_spikes, comp_spikes);
end % i

% compute latencies for unit pair
latencies = extract_latencies(ref_spikes, comp_spikes);

% make save_path
save_path = "./figures/Fig_5ABC";

% plot results (include save_path variable to save results)
plot_fig5ABC(sttc_range, dT_range, latencies, BIN_WIDTH, ...
    PLOT_WIDTH, ref_spikes, comp_spikes, SPIKE_TRAIN_RANGE, save_path)
 

%% Plot figure 5 D-F

NO_PAIRS = 90; % specify number of pairs with highest tile coefficient to plot
STTC_THRESH = 0.35; % specify edge strenght threshold for calculating node degree
FRAC_THRESH = 0.8; % specify in/out degree fraction for assigning receiver/sender
SUBTRACT = 0.35; % specify line thickness that should be subtracted from STTC value
LC = true; % specify whether only largest component channels should be plotted

% specify coordinates of example nodes to be plotted
NODE_MARKERS = [3342.5 , 122.5 ; ... % sender
                3027.5 , 122.5 ; ... % broker
                3657.5 , 507.5];     % receiver
                
% specify window sizes for connectivity map
XMIN = -100;
XMAX = 4000;YMIN = -100;
YMAX = 2200;

% load spike matrices
C_cont = load(sprintf("./data/processed/C_tile_coef_mat_whole_dir_dT_%.0f_minSpikes_%.0f.mat", ...
        dT*1000, MIN_SPIKES));
sttc_whole1 = C_cont.tile_coefficients_whole;
D50_cont = load(sprintf("./data/processed/D50_tile_coef_mat_whole_dir_dT_%.0f_minSpikes_%.0f.mat", ...
        dT*1000, MIN_SPIKES));
sttc_whole2 = D50_cont.tile_coefficients_whole;

% make save_path
save_path = "./figures/Fig_5DEF";

% plot results (include save_path variable to save results)
plot_fig5DEF(sttc_whole1, sttc_whole2, coordinates, STTC_THRESH, FRAC_THRESH, ...
    NO_PAIRS, SUBTRACT, LC, NODE_MARKERS, [XMIN, XMAX, YMIN, YMAX], save_path)

%% Functions

function plot_fig5DEF(sttc1, sttc2, coordinates, sttc_thresh, frac_thresh, ...
    no_pairs, subtract, LC, node_markers, sizes, save_path)

% initiate figure
cm = figure(2);
clf

% adjust size of figure
set(gcf,'PaperPositionMode','auto')
set(cm, 'Position', [0 20 1500 1000])
set(cm, 'Renderer', 'painters')

% % process data

% binarize matrices
bin_sttc1 = sttc1;
bin_sttc1(isnan(bin_sttc1)) = 0;
bin_sttc1(bin_sttc1<sttc_thresh) = 0;
bin_sttc1(bin_sttc1>=sttc_thresh) = 1;

bin_sttc2 = sttc2;
bin_sttc2(isnan(bin_sttc2)) = 0;
bin_sttc2(bin_sttc2<sttc_thresh) = 0;
bin_sttc2(bin_sttc2>=sttc_thresh) = 1;

% select LC units
LC_units1 = 1:size(sttc1,1);
LC_units2 = 1:size(sttc2,1);

if LC == true
    
    % extract the giant components channels of two sided binary matrix
    [~,LC_units1] = giant_component(symmetrize(bin_sttc1));
    [~,LC_units2] = giant_component(symmetrize(bin_sttc2));
    
    % obtain gc adjacency matrices
    bin_sttc1 = bin_sttc1(LC_units1,LC_units1);
    bin_sttc2 = bin_sttc2(LC_units2,LC_units2);
    sttc1 = sttc1(LC_units1,LC_units1);
    sttc2 = sttc2(LC_units2,LC_units2);
    
end % if

% obtain degree per channel
[tot_deg1, in_deg1, out_deg1] = degrees(bin_sttc1);
[tot_deg2, in_deg2, out_deg2] = degrees(bin_sttc2);

% calculate in out ratio and store
ratios1 = (in_deg1-out_deg1)./tot_deg1;
ratios2 = (in_deg2-out_deg2)./tot_deg2;


% % plot results
% initiate subplot for map of first treatment
subplot_tight(2, 2, 1, [0.05, 0.01])

% plot connectivity map
plot_connectivity_map(sttc1, tot_deg1, ratios1, ...
    coordinates(LC_units1,:), ...
    frac_thresh, no_pairs, subtract, node_markers, sizes);

% initiate subplot for map of second treatment
subplot_tight(2, 2, 2, [0.05, 0.01])

% plot connectivity map
plot_connectivity_map(sttc2, tot_deg2, ratios2, ...
    coordinates(LC_units2,:), ...
    frac_thresh, no_pairs, subtract, [], sizes);

% initiate subplot for sender node
subplot_tight(4, 6, 13, [0.02, 0.001])

% plot zoomed in with connections
plot_zoomed_dir_node(sttc1, coordinates(LC_units1,:), ...
    node_markers(1,:), sttc_thresh, [1,0,0], false)

% initiate subplot for receiver node
subplot_tight(4, 6, 14, [0.02, 0.001])

% plot zoomed in with connections
plot_zoomed_dir_node(sttc1, coordinates(LC_units1,:), ...
    node_markers(3,:), sttc_thresh, [0,0,1], true)

% initiate subplot for broker node
subplot_tight(4, 6, 15, [0.02, 0.001])

% plot zoomed in with connections
plot_zoomed_dir_node(sttc1, coordinates(LC_units1,:), ...
    node_markers(2,:), sttc_thresh, [0.5,0.5,0.5], false)


% if save_path is included in function call
if exist("save_path", "var")

    "saving"
    
    % save figure
    print(save_path, "-dpng", "-r600")
    
end % if


    % nested functions
    function plot_zoomed_dir_node(adj, coordinates, node_coordinates, thresh, color, leg)
        
        % column channel receives signals from row channels %
        % row channel sends signals to column channels %
        
        % find channel index of node to plot
        channel_index = find(coordinates(:,1) == node_coordinates(1) & coordinates(:,2) == node_coordinates(2));
        channel_index = channel_index(1);
        
        hold on
        
        % plot channel location
        plot(coordinates(channel_index,1), coordinates(channel_index,2), "o", ...
            "MarkerEdgeColor", color, "MarkerFaceColor", color, ...
            "MarkerSize", 20)
        
        % for each incoming signal
        for in_sig = 1:size(adj,1)
            
            % if signal is larger than thresh
            if adj(in_sig, channel_index) > thresh
                
                % plot line 
                plot([coordinates(channel_index,1), coordinates(in_sig,1)],...
                    [coordinates(channel_index,2), coordinates(in_sig,2)],...
                    "-", "Color", [0,0,1,0.5], "LineWidth", 2) 
                
            end % if
        end % in_sig
        
        % for each outgoing signal
        for out_sig = 1:size(adj,1)
            
            % if signal is larger than thresh
            if adj(channel_index, out_sig) > thresh
                
                % plot line 
                plot([coordinates(channel_index,1), coordinates(out_sig,1)],...
                    [coordinates(channel_index,2), coordinates(out_sig,2)],...
                    "-", "Color", [1,0,0,0.5], "LineWidth", 2) 
                
            end % if
        end % in_sig
        
        % plot lines for legend
        leg_line1 = plot([-1000,-1000],[-950,-950], "-", "Color", [0,0,1], "LineWidth", 3);
        leg_line2 = plot([-1000,-1000],[-950,-950], "-", "Color", [1,0,0], "LineWidth", 3);
        
        % add legend
        hleg = legend([leg_line2, leg_line1], {"Outgoing", "Incoming"}, "Location", ...
            "SouthOutside", "Orientation", "vertical", "FontSize", 14);
            
        if leg == false
            hleg.Visible = "off";
        end
        
        hold off
        
        % adjust axis
        axis equal
        xlim([coordinates(channel_index,1)-100, coordinates(channel_index,1)+100])
        ylim([coordinates(channel_index,2)-100, coordinates(channel_index,2)+100])
        xticks([])
        yticks([])
        set(gca, 'YDir','reverse')
        
        ax = gca;
        ax.FontSize = 14;
        
    end % fun plot_zoomed_dir_node
    
    % % %
    
    function plot_connectivity_map(adj, tot_deg, ratios, coordinates, frac_thresh, ...
            no_pairs, subtract, node_markers, sizes)
   
        % specify multiplier for plot sizes
        WIDTH_MULTIPLY = 12;
        BUB_MULTIPLY = 6;
        LINE_COLOR = [0.25,0.25,0.25];
        BUB_COLOR = [0.5,0.5,0.5];
        TRANSPARENCY = 0.2;
        
        % make empty result arrays for sender and receiver nodes
        senders = [];
        receivers = [];
        brokers = [];
        
        % % plot results
        hold on 
        
        % for each channel
        for channel = 1:length(ratios)
            
            % if out_degree is larger than in_degree
            if ratios(channel) > -frac_thresh && ratios(channel) < frac_thresh
                
                % plot channel location with size of total degree
                plot(coordinates(channel,1), coordinates(channel,2), "o", ...
                    "MarkerEdgeColor", [0.5,0.5,0.5], "MarkerFaceColor", [0.5,0.5,0.5], "MarkerSize", 5)
                
            end % if
            
        end % channel
        
        % for each channel
        for channel = 1:length(ratios)
            
            % if out_degree is larger than in_degree
            if ratios(channel) < -frac_thresh 
                
                % store channel
                senders = [senders; channel];
                
                % plot channel location with size of total degree
                plot(coordinates(channel,1), coordinates(channel,2), "o", ...
                    "MarkerEdgeColor", [1,0,0], "MarkerFaceColor", [1,0,0], ...
                    "MarkerSize", sqrt(tot_deg(channel)*BUB_MULTIPLY))
      
                
            elseif ratios(channel) > frac_thresh
                
                % store channel
                receivers = [receivers; channel];
                
                % plot channel location with size of total degree
                plot(coordinates(channel,1), coordinates(channel,2), "o", ...
                    "MarkerEdgeColor", [0,0,1], "MarkerFaceColor", [0,0,1], ...
                    "MarkerSize", sqrt(tot_deg(channel)*BUB_MULTIPLY))
                
            end % if
        end % channel        
        
        if ~isempty(node_markers)
                        
            % plot locations of nodes of interest
            plot(node_markers(1,1), node_markers(1,2), "o", ...
                "LineWidth", 2, "MarkerEdgeColor", "k", "MarkerFaceColor", [1,0,0], ...
                "MarkerSize", sqrt(20*BUB_MULTIPLY))
            
            plot(node_markers(2,1), node_markers(2,2), "o", ...
                "LineWidth", 2, "MarkerEdgeColor", "k", "MarkerFaceColor", [0.5,0.5,0.5], ...
                "MarkerSize", sqrt(20*BUB_MULTIPLY))
            
            plot(node_markers(3,1), node_markers(3,2), "o", ...
                "LineWidth", 2, "MarkerEdgeColor", "k", "MarkerFaceColor", [0,0,1], ...
                "MarkerSize", sqrt(20*BUB_MULTIPLY))
            
        end
        
        % column channel receives signals from row channels %
        % row channel sends signals to column channels %
        
        % % select top x signals from senders and receivers
        % make adj matrix of all outgoing signals from senders
        send_matrix = adj(senders,:);
        send_coordinates = coordinates(senders,:);
        
        % sort all coefficients from high to low
        sorted_send_coefs = sort(send_matrix(~isnan(send_matrix)), "descend");
        
        % obtain threshold based on no_pairs to plot
        send_threshold = sorted_send_coefs(no_pairs);
        
        % make adj matrix of all outgoing signals from senders
        receive_matrix = adj(:,receivers);
        receive_coordinates = coordinates(receivers,:);
        
        % sort all coefficients from high to low
        sorted_receive_coefs = sort(receive_matrix(~isnan(receive_matrix)), "descend");
        
        % obtain threshold based on no_pairs to plot
        receive_threshold = sorted_receive_coefs(no_pairs);
        
        
        % for each reference channel in send matrix
        for ref_channel =  1:size(send_matrix,2)
            
            % for each comparison channel
            for comp_channel = 1:size(send_matrix,1)
                
                % if tile coefficient lies above threshold
                if send_matrix(comp_channel, ref_channel) > abs(send_threshold)
                        
                    % plot line between ref_channel and comp_channel with thickness of
                    % tile_coef
                    plot([coordinates(ref_channel,1),send_coordinates(comp_channel,1)],...
                        [coordinates(ref_channel,2),send_coordinates(comp_channel,2)],...
                        "-", "Color", [LINE_COLOR, TRANSPARENCY], "LineWidth", ...
                        WIDTH_MULTIPLY*(send_matrix(comp_channel, ref_channel)-subtract)) % *send_matrix(comp_channel, ref_channel)
                    
                end % if
                
            end % comp_channel
        end % ref_channel
        
        % for each reference channel in receive matrix
        for ref_channel =  1:size(receive_matrix,2)
            
            % for each comparison channel
            for comp_channel = 1:size(receive_matrix,1)
                
                % if tile coefficient lies above threshold
                if receive_matrix(comp_channel, ref_channel) > abs(receive_threshold)
                        
                    % plot line between ref_channel and comp_channel with thickness of
                    % tile_coef
                    plot([receive_coordinates(ref_channel,1),coordinates(comp_channel,1)],...
                        [receive_coordinates(ref_channel,2),coordinates(comp_channel,2)],...
                        "-", "Color", [LINE_COLOR, TRANSPARENCY], "LineWidth", ...
                        WIDTH_MULTIPLY*(receive_matrix(comp_channel, ref_channel)-subtract))
                    
                end % if
                
            end % comp_channel
        end % ref_channel    
            
        % plot lines to be used for legend
        leg_line1 = plot([-1000,-1000],[-950,-950], "-", "Color", [LINE_COLOR, 2*TRANSPARENCY*0.6], ...
            "LineWidth", WIDTH_MULTIPLY*(0.6-subtract));
        leg_line2 = plot([-1000,-1000],[-950,-950], "-", "Color", [LINE_COLOR 2*TRANSPARENCY*0.8], ...
            "LineWidth", WIDTH_MULTIPLY*(0.8-subtract));
        leg_line3 = plot([-1000,-1000],[-950,-950], "-", "Color", [LINE_COLOR 2*TRANSPARENCY*1], ...
            "LineWidth", WIDTH_MULTIPLY*(1-subtract));
        
        % plot circles to be used for second legend
        leg_bub1 = plot(-1000, -1000, "o", "MarkerEdgeColor", "k", ...
            "MarkerFaceColor", [1,1,1], "MarkerSize", sqrt(5*BUB_MULTIPLY));
        leg_bub2 = plot(-1000, -1000, "o", "MarkerEdgeColor", "k", ...
            "MarkerFaceColor", [1,1,1], "MarkerSize", sqrt(20*BUB_MULTIPLY));
        leg_bub3 = plot(-1000, -1000, "o", "MarkerEdgeColor", "k", ...
            "MarkerFaceColor", [1,1,1], "MarkerSize", sqrt(50*BUB_MULTIPLY));
        
        % add scale bar
        plot([0,500],[2100,2100], "k", "LineWidth", 5)
        
        % add legend
        Hleg = legend([leg_line1, leg_line2, leg_line3, leg_bub1, ...
            leg_bub2, leg_bub3], {"0.6", "0.8", "1", "5", ...
            "20", "50"}, "Location", "SouthOutside", ...
            "Orientation", "horizontal");
        
        % add legend titles
        title(Hleg,"Tile coefficient              Node degree")
        
        % adjust axis
        axis equal
        xlim([sizes(1), sizes(2)])
        ylim([sizes(3), sizes(4)])
        xticks([])
        yticks([])
        set(gca, 'YDir','reverse')
        
        ax = gca;
        ax.FontSize = 14;
            
        hold off

    end % fun plot_degree_map

end % fun plot_fig5DEF

% % %

function plot_fig5ABC(sttc_range, dT_range, latencies, bin_width, plot_width, ...
    ref_spikes, comp_spikes, spike_train_range, save_path)

% initiate figure
cm = figure(1);
clf

% adjust size of figure
set(gcf,'PaperPositionMode','auto')
set(cm, 'Position', [0 20 1000 300])
set(cm, 'Renderer', 'painters')


% initiate first subplot for spike train
subplot(3,3,4)

% plot spike trains
plot_spike_trains(ref_spikes, comp_spikes, spike_train_range)


% initiate subplot for spike time latency histogram
subplot(3,3,[2,5,8])

% plot spike time latency histogram
plot_time_diff_histogram(latencies, bin_width, plot_width)


% initiate subplot for sttc range
subplot(3,3,[3,6,9])

% plot sttc range results
make_single_sttc_range_plot(sttc_range, dT_range)

% if save_path is included in function call
if exist("save_path", "var")

    "saving"
    
    % save figure
    print(save_path, "-dpng")
end % if


    % nested functions
    function plot_spike_trains(ref_spikes, comp_spikes, range)

        % plot results
        hold on
        plot([ref_spikes; ref_spikes], ...
            [ones(size(ref_spikes, 2), 1), 2*ones(size(ref_spikes, 2), 1)]', ...
            "color", [0.5, 0.5, 0.5]);
        plot([comp_spikes; comp_spikes], ...
            [2.5*ones(size(comp_spikes, 2), 1), 3.5*ones(size(comp_spikes, 2), 1)]', ...
            "color", [0.5, 0.5, 0.5]);
        
        ax = gca;
        ax.FontSize = 14;
        
        xlim([range(1), range(2)])
        
        % improve axis
        yticks([])
        xticks([])
        ylim([0, 3.6])
            
    end % fun plot_spike_trains

    % % %
    
    function plot_time_diff_histogram(time_differences, bin_width, high_cut)
        
        low_cut = -high_cut;
        
        % determine bin edges
        bin_edges = linspace(low_cut, high_cut, ((high_cut-low_cut)/bin_width)+1);
        
        % plot histogram with specified bin edges
        histogram(time_differences, bin_edges, 'FaceColor', [0 0.4470 0.7410])
        
        ax = gca;
        ax.FontSize = 14;
        
        % adjust axes
        xlim([low_cut,high_cut])
        ylim([0,ceil(1.1*max(histcounts(time_differences, bin_edges)))])
        
        if high_cut == 0.02
            
            xticks(linspace(low_cut, high_cut, 9))
            xticklabels(linspace(low_cut*1000, high_cut*1000, 9))
            xlabel("Time difference (ms)", "fontsize", 16)
            
        elseif high_cut < 0.2
            
            xticks(linspace(low_cut, high_cut, 7))
            xticklabels(linspace(low_cut*1000, high_cut*1000, 7))
            xlabel("Time difference (ms)", "fontsize", 16)
            
        else
            
            xlabel("Time difference (s)", "fontsize", 16)
            
        end % if
        
        ylabel("Number of spike pairs", "fontsize", 16)
        
    end % fun plot_time_diff_histogram
    
    % % %
    
    function make_single_sttc_range_plot(sttc_range, dT_range)
        
        hold on
        
        % plot sttc range
        plot(dT_range, sttc_range, "Color", "k", "LineWidth", 2);

        % plot vertical line to indicate time cutoff
        xline(0.02, "k--", "LineWidth", 2, 'HandleVisibility', 'off');
            
        % adjust axis
        ax = gca;
        ax.FontSize = 14;
        
        xlabel("Δt (ms)", "fontsize", 16)
        ylabel("STTC", "fontsize", 16)
        xticks(linspace(0,max(dT_range), 1+max(dT_range)*1000/30))
        xticklabels(linspace(0,max(dT_range)*1000, 1+max(dT_range)*1000/30))
        ylim([0,1])

    end % fun make_sttc_range_plot

end % fun plot_fig5ABC

% % %

function latencies = extract_latencies(ref_spike_times, comp_spike_times)

% make empty result array
latencies = [];

% for each spike time in the reference channel
for ref_time = 1:length(ref_spike_times)
    
    % find which spike in the comparison channel lies closest
    closest_comp = comp_spike_times(find(abs(ref_spike_times(ref_time) - ...
        comp_spike_times) == min(abs(ref_spike_times(ref_time) - comp_spike_times))));
    
    % obtain time difference
    time_diff = closest_comp(1)-ref_spike_times(ref_time);
    
    % if the time difference is positive
    if time_diff > 0
        
        % determine if next ref spike is closer than closest comp spike
        if ref_time == length(ref_spike_times) || ...
            ref_spike_times(ref_time+1) - ref_spike_times(ref_time) > time_diff
            
            % store result
            latencies = [latencies, time_diff];
            
        end % if
        
    elseif time_diff < 0
        
        % determine if next ref spike is closer than closest comp spike
        if ref_time == 1 || ref_spike_times(ref_time-1) - ref_spike_times(ref_time) < time_diff
            
            % store result
            latencies = [latencies, time_diff];
            
        end % if
    end % if

end % ref_time

end % fun extract_latencies

% % %

function tileCoef = sttc(N1v, N2v, dtv, Time, spike_times_1, spike_times_2)

%sttc calculates the spike timing tiling coefficient 
%   Originally written in C by Cutts and Eglen 
% https://github.com/CCutts/Detecting_pairwise_correlations_in_spike_trains/blob/master/spike_time_tiling_coefficient.c
%   See their 2014 paper on this 
% https://www.ncbi.nlm.nih.gov/pubmed/25339742
%   Tranlsated to matlab by Tim Sit (sitpakhang@gmail.com)
% https://github.com/Timothysit/organoids/blob/master/correlation_analysis/sttc.m
%   See their paper on applying tileCoef to MEA recordings of neural
%   organoids
% https://www.nature.com/articles/s41593-019-0350-2

% INPUTS 
    % N1v           | The number of spikes in electrode 1 (double) 
    % N2v           | The number of spikes in electrode 2 (double)
    % dtv           | The delay (in seconds) (double) 
    % Time          | 2 x 1 vector containing the start time and end time 
    % of the recording (seconds), so that Time(2) - Time(1) = length of 
    % recording
    % spike_times_1 | The spike times in electrode 1 (in seconds)  (vector)
    % spike_times_2 | the spikes times in electrode 2 (in seconds) (vector)

% OUTPUT

    % tileCoef | The tiling coefficient
 

N1 = N1v; % need to think about these gloabl var names
N2 = N2v; 
dt = dtv; 

if N1 == 0 || N2 == 0 
   %  index = R_NaN; % I think this just means NaN values
   index = NaN; 
   
else 
    T = Time(2) - Time(1); 
    TA = run_T(N1, dt, Time(1), Time(2), spike_times_1); 
    TA = TA / T; 
    TB = run_T(N2, dt, Time(1), Time(2), spike_times_2); 
    TB = TB / T; 
    PA = run_P(N1, N2, dt, spike_times_1, spike_times_2); 
    PA = PA / N1; 
    PB = run_P(N2, N1, dt, spike_times_2, spike_times_1); 
    PB = PB / N2; 
    index = 0.5 * (PA - TB) / (1 - TB * PA) + 0.5 * (PB - TA) / (1 - TA * PB);
end 

tileCoef = index; 


    function Nab = run_P(N1, N2, dt, spike_times_1, spike_times_2)
        Nab = 0; 
        j = 1; % change to 1 for 1 indexing
        % also note the switch from 0 to 1 indexing
        for i = 1:N1
            while j <= N2 % changed to <= for 1 indexing
                if abs(spike_times_1(i) - spike_times_2(j)) <= dt
                    Nab = Nab + 1; 
                    break 
                elseif spike_times_2(j) > spike_times_1(i) 
                    break
                else 
                    j = j + 1;
                end 
            end 
        end 
    end


    function time_A = run_T(N1v, dtv, startv, endv, spike_times_1)
        dt = dtv; 
        start = startv; 
        endvv = endv; % end is not a valid variable name in MATLAB 
        tempN = N1v; % changed N1 into tempN as nested function variables are declared
        % globally. This is problematic when N1v and N2v are different
        % values
        
        % maximum
        time_A = 2 * tempN * dt; 
        
        % if just one spike in train 
        if tempN == 1
            
            if spike_times_1(1) - start < dt 
                time_A = time_A - start + spike_times_1(1) - dt; 
            elseif spike_times_1(1) + dt > endvv 
                time_A = time_A - spike_times_1(1) - dt + endvv; 
            end
        
            % if more than one spike in train 
        else 
            i = 1; % added by TS
            while i < tempN % switched from N1 - 1, to take account of 1 indexing
            
                diff = spike_times_1(i+1) - spike_times_1(i); 
                
                if diff < 2 * dt 
                    % subtract overlap 
                    time_A = time_A -2 * dt + diff; 
                end 
                 
                i = i + 1; 
            end 
            
            % check if spikes are within dt of the start and/or end, if so
            % just need to subtract overlap of first and/or last spike as
            % all within-train overlaps have been accounted for 
            
            if spike_times_1(1) - start < dt 
                time_A = time_A - start + spike_times_1(1) - dt; 
            end 
            
            if endvv - spike_times_1(tempN) < dt % switched from N1 - 1 to N1 to for 1 indexing
                time_A = time_A - spike_times_1(tempN) - dt + endvv; 
            end  
        end       
    end 
end % fun tileCoef

% % %

function		[dip, p_value, xlow,xup]=HartigansDipSignifTest(xpdf,nboot)

%  function		[dip,p_value,xlow,xup]=HartigansDipSignifTest(xpdf,nboot)
%
% calculates Hartigan's DIP statistic and its significance for the empirical p.d.f  XPDF (vector of sample values)
% This routine calls the matlab routine 'HartigansDipTest' that actually calculates the DIP
% NBOOT is the user-supplied sample size of boot-strap
% Code by F. Mechler (27 August 2002)

% calculate the DIP statistic from the empirical pdf
[dip,xlow,xup, ifault, gcm, lcm, mn, mj]=HartigansDipTest(xpdf);
N=length(xpdf);

% calculate a bootstrap sample of size NBOOT of the dip statistic for a uniform pdf of sample size N (the same as empirical pdf)
boot_dip=[];
for i=1:nboot
   unifpdfboot=sort(unifrnd(0,1,1,N));
   [unif_dip]=HartigansDipTest(unifpdfboot);
   boot_dip=[boot_dip; unif_dip];
end;
boot_dip=sort(boot_dip);
p_value=sum(dip<boot_dip)/nboot;

% % Plot Boot-strap sample and the DIP statistic of the empirical pdf
% figure(1); clf;
% [hy,hx]=hist(boot_dip); 
% bar(hx,hy,'k'); hold on;
% plot([dip dip],[0 max(hy)*1.1],'r:');

end % fun HartigansDipSignifTest

% % % 

function	[dip,xl,xu, ifault, gcm, lcm, mn, mj] = HartigansDipTest(xpdf)

% function	[dip,xl,xu, ifault, gcm, lcm, mn, mj]=HartigansDipTest(xpdf)
%
% This is a direct translation by F. Mechler (August 27 2002)
% into MATLAB from the original FORTRAN code of Hartigan's Subroutine DIPTST algorithm 
% Ref: Algorithm AS 217 APPL. STATIST. (1985) Vol. 34. No.3 pg 322-325
%
% Appended by F. Mechler (September 2 2002) to deal with a perfectly unimodal input
% This check the original Hartigan algorithm omitted, which leads to an infinite cycle
%
% HartigansDipTest, like DIPTST, does the dip calculation for an ordered vector XPDF using
% the greatest convex minorant (gcm) and the least concave majorant (lcm),
% skipping through the data using the change points of these distributions.
% It returns the 'DIP' statistic, and 7 more optional results, which include
% the modal interval (XL,XU), ann error flag IFAULT (>0 flags an error)
% as well as the minorant and majorant fits GCM, LCM, and the corresponding support indices MN, and MJ

% sort X in increasing order in column vector
x=sort(xpdf(:));
N=length(x);
mn=zeros(size(x));
mj=zeros(size(x));
lcm=zeros(size(x));
gcm=zeros(size(x));
ifault=0;

% Check that N is positive
if (N<=0) 
   ifault=1;
   fprintf(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);
   return;
end;

% Check if N is one
if (N==1)
   xl=x(1);
   xu=x(N);
   dip=0.0;
   ifault=2;
   fprintf(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);
   return;
end;

if (N>1)
   % Check that X is sorted
   if (x ~= sort(x))
      ifault=3;
      fprintf(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);
      return;
   end;
   % Check for all values of X identical OR for case 1<N<4
   if ~((x(N)>x(1)) & (4<=N))
      xl=x(1);
      xu=x(N);
      dip=0.0;
      ifault=4;
      fprintf(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);
      return;
   end;
end;

% Check if X is perfectly unimodal
% Hartigan's original DIPTST algorithm did not check for this condition
% and DIPTST runs into infinite cycle for a unimodal input
% The condition that the input is unimodal is equivalent to having 
% at most 1 sign change in the second derivative of the input p.d.f.
xsign=-sign(diff(diff(x)));
% This condition check below works even 
% if the unimodal p.d.f. has its mode in the very first or last point of the input 
% because then the boolean argument is Empty Matrix, and ANY returns 1 for an Empty Matrix
posi=find(xsign>0);
negi=find(xsign<0);
if isempty(posi) | isempty(negi) | all(posi<min(negi))
   % A unimodal function is its own best unimodal approximation, with a zero corresponding dip
   xl=x(1);
   xu=x(N);
   dip=0.0;
   ifault=5;
	%fprintf(1,'\n  The input is a perfectly UNIMODAL input function\n');
   return;
end;

% LOW  contains the index of the current estimate of the lower end of the modal interval
% HIGH contains the index of the current estimate of the upper end of the modal interval
fn=N;
low=1;
high=N;
dip=1/fn;
xl=x(low);
xu=x(high);

% establish the indices over which combination is necessary for the convex minorant fit
mn(1)=1;
for j=2:N
   mn(j)=j-1;
   % here is the beginning of a while loop
   mnj=mn(j);
   mnmnj=mn(mnj);
   a=mnj-mnmnj;
   b=j-mnj;
   while ~( (mnj==1) | ((x(j)-x(mnj))*a < (x(mnj)-x(mnmnj))*b))
      mn(j)=mnmnj;
      mnj=mn(j);
      mnmnj=mn(mnj);
      a=mnj-mnmnj;
      b=j-mnj;
   end;   % here is the end of the while loop
end; % end  for j=2:N

% establish the indices over which combination is necessary for the concave majorant fit
mj(N)=N;
na=N-1;
for jk=1:na
   k=N-jk;
   mj(k)=k+1;
   % here is the beginning of a while loop
   mjk=mj(k);
   mjmjk=mj(mjk);
   a=mjk-mjmjk;
   b=k-mjk;
   while ~( (mjk==N) | ((x(k)-x(mjk))*a < (x(mjk)-x(mjmjk))*b))
      mj(k)=mjmjk;
      mjk=mj(k);
      mjmjk=mj(mjk);
      a=mjk-mjmjk;
      b=k-mjk;
   end;   % here is the end of the while loop
end; % end  for jk=1:na

itarate_flag = 1;

% start the cycling of great RECYCLE
while itarate_flag 

% collect the change points for the GCM from HIGH to LOW
% CODE BREAK POINT 40
ic=1;
gcm(1)=high;
igcm1=gcm(ic);
ic=ic+1;
gcm(ic)=mn(igcm1);
while(gcm(ic) > low)
   igcm1=gcm(ic);
   ic=ic+1;
   gcm(ic)=mn(igcm1);
end;
icx=ic;

% collect the change points for the LCM from LOW to HIGH
ic=1;
lcm(1)=low;
lcm1=lcm(ic);
ic=ic+1;
lcm(ic)=mj(lcm1);
while(lcm(ic) < high)
   lcm1=lcm(ic);
   ic=ic+1;
   lcm(ic)=mj(lcm1);
end;
icv=ic;

% ICX, IX, IG are counters for the convex minorant
% ICV, IV, IH are counters for the concave majorant
ig=icx;
ih=icv;

% find the largest distance greater than 'DIP' between the GCM and the LCM from low to high
ix=icx-1;
iv=2;
d=0.0;

% Either GOTO CODE BREAK POINT 65 OR ELSE GOTO CODE BREAK POINT 50;
if ~(icx~=2 | icv~=2)
   d=1.0/fn;
else
   iterate_BP50=1;
   while iterate_BP50
		% CODE BREAK POINT 50
		igcmx=gcm(ix);
      lcmiv=lcm(iv);
      if ~(igcmx > lcmiv)
         % if the next point of either the GCM or LCM is from the LCM then calculate distance here
         % OTHERWISE, GOTO BREAK POINT 55
         lcmiv1=lcm(iv-1);
         a=lcmiv-lcmiv1;
         b=igcmx-lcmiv1-1;
         dx=(x(igcmx)-x(lcmiv1))*a/(fn*(x(lcmiv)-x(lcmiv1)))-b/fn;
         ix=ix-1;
         if(dx < d) 
            goto60 = 1; 
         else
            d=dx;
            ig=ix+1;
            ih=iv;
            goto60 = 1;
         end;
      else
         % if the next point of either the GCM or LCM is from the GCM then calculate distance here
         % CODE BREAK POINT 55
         lcmiv=lcm(iv);
         igcm=gcm(ix);
         igcm1=gcm(ix+1);
         a=lcmiv-igcm1+1;
         b=igcm-igcm1;
         dx=a/fn-((x(lcmiv)-x(igcm1))*b)/(fn*(x(igcm)-x(igcm1)));
         iv=iv+1;
         if ~(dx < d) 
            d=dx;
            ig=ix+1;
            ih=iv-1;
         end;
         goto60 = 1;
      end;
      
      if goto60
         % CODE BREAK POINT 60
         if (ix < 1) ix=1; end;
         if (iv > icv) iv=icv; end;
         iterate_BP50 = (gcm(ix) ~= lcm(iv)); 
      end;
   end; % End of WHILE iterate_BP50
end; % End of ELSE (IF ~(icx~=2 | icv~=2)) i.e., either GOTO CODE BREAK POINT 65 OR ELSE GOTO CODE BREAK POINT 50

% CODE BREAK POINT 65
itarate_flag = ~(d < dip);
if itarate_flag
% if itarate_flag is true, then continue calculations and the great iteration cycle
% if itarate_flag is NOT true, then stop calculations here, and break out of great iteration cycle to BREAK POINT 100
   
% calculate the DIPs for the corrent LOW and HIGH

% the DIP for the convex minorant
dl=0.0;
% if not true, go to CODE BREAK POINT 80
if (ig ~= icx)
   icxa=icx-1;
   for j=ig:icxa
      temp=1.0/fn;
   	jb=gcm(j+1);
      je=gcm(j);
      % if not true either, go to CODE BREAK POINT 74
      if ~(je-jb <= 1)
         if~(x(je)==x(jb))
            a=(je-jb);
            const=a/(fn*(x(je)-x(jb)));
            for jr=jb:je
               b=jr-jb+1;
               t=b/fn-(x(jr)-x(jb))*const;
               if (t>temp) temp=t; end;
            end;
         end;
      end;
      % CODE BREAK POINT 74
      if (dl < temp) dl=temp; end;
   end;
end;

% the DIP for the concave majorant
% CODE BREAK POINT 80
du=0.0;
% if not true, go to CODE BREAK POINT 90
if ~(ih==icv)
   icva=icv-1;
   for k=ih:icva
      temp=1.0/fn;
      kb=lcm(k);
      ke=lcm(k+1);
      % if not true either, go to CODE BREAK POINT 86
      if ~(ke-kb <= 1)
         if ~(x(ke)==x(kb))
            a=ke-kb;
            const=a/(fn*(x(ke)-x(kb)));
            for kr=kb:ke
               b=kr-kb-1;
               t=(x(kr)-x(kb))*const-b/fn;
               if (t>temp) temp=t; end;
            end;
         end;
      end;
      % CODE BREAK POINT 86
      if (du < temp) du=temp; end;
   end;
end;

% determine the current maximum
% CODE BREAK POINT 90
dipnew=dl;
if (du > dl) dipnew=du; end;
if (dip < dipnew) dip=dipnew; end;
low=gcm(ig);
high=lcm(ih);      

end; % end of IF(itarate_flag) CODE from BREAK POINT 65

% return to CODE BREAK POINT 40 or break out of great RECYCLE;
end; % end of WHILE of great RECYCLE

% CODE BREAK POINT 100
dip=0.5*dip;
xl=x(low);
xu=x(high);

end % fun HartigansDipTest

% % %

function vargout=subplot_tight(m, n, p, margins, varargin)
% subplot_tight
% A subplot function substitude with margins user tunabble parameter.
%
% Syntax
%  h=subplot_tight(m, n, p);
%  h=subplot_tight(m, n, p, margins);
%  h=subplot_tight(m, n, p, margins, subplotArgs...);
%
% Description
% Our goal is to grant the user the ability to define the margins between neighbouring
%  subplots. Unfotrtunately Matlab subplot function lacks this functionality, and the
%  margins between subplots can reach 40% of figure area, which is pretty lavish. While at
%  the begining the function was implememnted as wrapper function for Matlab function
%  subplot, it was modified due to axes del;etion resulting from what Matlab subplot
%  detected as overlapping. Therefore, the current implmenetation makes no use of Matlab
%  subplot function, using axes instead. This can be problematic, as axis and subplot
%  parameters are quie different. Set isWrapper to "True" to return to wrapper mode, which
%  fully supports subplot format.
%
% Input arguments (defaults exist):
%   margins- two elements vector [vertical,horizontal] defining the margins between
%        neighbouring axes. Default value is 0.04
%
% Output arguments
%   same as subplot- none, or axes handle according to function call.
%
% Issues & Comments
%  - Note that if additional elements are used in order to be passed to subplot, margins
%     parameter must be defined. For default margins value use empty element- [].
%  - 
%
% Example
% close all;
% img=imread('peppers.png');
% figSubplotH=figure('Name', 'subplot');
% figSubplotTightH=figure('Name', 'subplot_tight');
% nElems=17;
% subplotRows=ceil(sqrt(nElems)-1);
% subplotRows=max(1, subplotRows);
% subplotCols=ceil(nElems/subplotRows);
% for iElem=1:nElems
%    figure(figSubplotH);
%    subplot(subplotRows, subplotCols, iElem);
%    imshow(img);
%    figure(figSubplotTightH);
%    subplot_tight(subplotRows, subplotCols, iElem, [0.0001]);
%    imshow(img);
% end
%
% See also
%  - subplot
%
% Revision history
% First version: Nikolay S. 2011-03-29.
% Last update:   Nikolay S. 2012-05-24.
%
% *List of Changes:*
% 2012-05-24
%  Non wrapping mode (based on axes command) added, to deal with an issue of disappearing
%     subplots occuring with massive axes.
% Default params
isWrapper=false;
if (nargin<4) || isempty(margins)
    margins=[0.04,0.04]; % default margins value- 4% of figure
end
if length(margins)==1
    margins(2)=margins;
end
%note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
[subplot_col,subplot_row]=ind2sub([n,m],p);  
height=(1-(m+1)*margins(1))/m; % single subplot height
width=(1-(n+1)*margins(2))/n;  % single subplot width
% note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
subplot_cols=1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot 
subplot_rows=1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot   
merged_height=subplot_rows*( height+margins(1) )- margins(1);   % merged subplot height
merged_width= subplot_cols*( width +margins(2) )- margins(2);   % merged subplot width
merged_bottom=(m-max(subplot_row))*(height+margins(1)) +margins(1); % merged subplot bottom position
merged_left=min(subplot_col)*(width+margins(2))-width;              % merged subplot left position
pos=[merged_left, merged_bottom, merged_width, merged_height];
if isWrapper
   h=subplot(m, n, p, varargin{:}, 'Units', 'Normalized', 'Position', pos);
else
   h=axes('Position', pos, varargin{:});
end
if nargout==1
   vargout=h;
end

end % fun subplot_tight
