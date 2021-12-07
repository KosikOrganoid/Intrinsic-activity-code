%% load data

% % specify which organoids to analyse
ORGANOIDS = {"L1", 0, 0}; % only test data is included
TREATMENTS = {"C", "D50"};

% make empty cell arrays for results
sttc_whole = cell(length(TREATMENTS),length(ORGANOIDS));
coordinates = cell(1,length(ORGANOIDS));

% for each array
for organoid = 1:length(ORGANOIDS)
    
    if isstring(ORGANOIDS{organoid})
        
        for treatment = 1:length(TREATMENTS)
            
            % load STTC data (from Fig5.m script)
            sttc_cont = load(sprintf("./data/processed/%s_tile_coef_mat_whole_dir_dT_20_minSpikes_5.mat", ...
                TREATMENTS{treatment}));
            sttc_whole{treatment, organoid} = sttc_cont.tile_coefficients_whole;
            
        end % treatment
        
        % load coordinates
        t_spk_mat_cont = load("./data/C_t_spk_mat_sorted.mat");
        coordinates{organoid} = t_spk_mat_cont.xy_raw;
        
    end % if
    
end % organoid

%% Process data

STTC_THRESH = 0.35;

edges_mat = cell(size(sttc_whole));
all_edges = cell(size(sttc_whole));

% for each organoid
for organoid = 1:size(sttc_whole,2)
    
    if ~isempty(sttc_whole{1,organoid})
        
        organoid
        
        % for each treatment
        for treatment = 1:size(sttc_whole,1)
            
            % binarize matrix
            bin_sttc_whole = sttc_whole{treatment, organoid};
            bin_sttc_whole(isnan(bin_sttc_whole)) = 0;
            bin_sttc_whole(bin_sttc_whole<STTC_THRESH) = 0;
            bin_sttc_whole(bin_sttc_whole>=STTC_THRESH) = 1;
            
            % extract the giant components channels of two sided binary matrix
            [~,LC_units] = giant_component(symmetrize(bin_sttc_whole));
            
            % determine non-LC channels
            non_LC = ~ismember(LC_units, 1:size(bin_sttc_whole,1));

            % copy sttc matrix and store results
            edges_mat{treatment, organoid} = sttc_whole{treatment, organoid};
            
            % remove non LC units
            edges_mat{treatment, organoid}(non_LC,:) = 0;
            edges_mat{treatment, organoid}(:,non_LC) = 0;
            
            % remove below threshold values 
            edges_mat{treatment, organoid}(edges_mat{treatment, organoid} ...
                < STTC_THRESH) = 0;
            
            % obtain list of all edge strengths and store
            all_edges{treatment, organoid} = edges_mat{treatment, organoid} ...
                (edges_mat{treatment, organoid} > STTC_THRESH);
            
        end % treatment
        
    end % if
    
end % array

%% Plot figure 6

% define plot parameters
ARRAY_OI = 1; % specify array of interest for plot 
BIN_SIZE = 0.05; % specify size to bin STTC values
TOP_EDGES = 90; % specify number of pairs with highest tile coefficient to plot

% define symbols to be used for plotting
symbols = ["o", "d"];
colors = ["b", "r"];

% define gamma distribution
modelfun = @(b,x) b(3) .* x(:,1) .^ (b(1)-1) .* exp(-x(:,1) ./ (b(2))) ./ (b(2).^(b(1)).*gamma(b(1)));

% specify bin edges
bin_edges = linspace(0.35,1,1+(1-0.35)/BIN_SIZE);

% make save path
save_path = "./figures/Fig_6";

% plot figure
plot_fig6(edges_mat, all_edges, coordinates, bin_edges, symbols, colors, ...
    modelfun, ARRAY_OI, TOP_EDGES, save_path)

%% functions

function plot_fig6(edges_mat, all_edges, coordinates, bin_edges, ...
    symbols, colors, modelfun, array_oi, top_edges, save_path)

% initiate figure
cm = figure(1);
clf

% adjust size of figure
set(gcf,'PaperPositionMode','auto')
set(cm, 'Position', [0 20 1200 1000])
set(cm, 'Renderer', 'painters')

% % plot edge strength distribution for array_oi
subplot(3,3,1)

[fits, xfits, lh] = plot_edge_dist(all_edges, array_oi, bin_edges, modelfun, symbols, colors);
legend(lh, "Control", "Diazepam", "Location", "NorthEast")


% % plot fractional change in edge strength distribution for all arrays
subplot(3,3,2)

plot_edge_dist_frac_change(all_edges, fits, xfits, bin_edges, {[1,0,0] ; [0,0,1] ; [0.4660 0.6740 0.1880]}) 


% % process data based on which edges are present for each recording
[edge_strengths_both, edge_strengths_c, edge_strengths_d, ...
    coor_both, coor_c, coor_d] = process_data(edges_mat, coordinates);


% % make barplot with number of edges in each group
subplot(3,3,3)

plot_group_count_bar(edge_strengths_both, edge_strengths_c, edge_strengths_d)

legend("Both C and D", "Only C", "Only D", "Location", "EastOutside")
    

% % plot distributions of individual group
subplot(3,3,4)

[~, ~, lh] = plot_edge_dist({edge_strengths_both{array_oi}(:,1); edge_strengths_both{array_oi}(:,2)}, 1, bin_edges, ...
    modelfun, symbols, colors);
legend(lh, "Control", "Diazepam", "Location", "NorthEast")
title(sprintf("%.0f edges", size(edge_strengths_both{array_oi},1)))

subplot(3,3,5)

[~, ~, lh] = plot_edge_dist(edge_strengths_c(array_oi), 1, bin_edges, modelfun, ...
    symbols(1), colors(1));
title(sprintf("%.0f edges", length(edge_strengths_c{array_oi})))

subplot(3,3,6)

[~, ~, lh] = plot_edge_dist(edge_strengths_d(array_oi), 1, bin_edges, modelfun, ...
    symbols(2), colors(2));
title(sprintf("%.0f edges", length(edge_strengths_d{array_oi})))


% % plot maps of individual group
subplot(3,3,7)

plot_map(coordinates{array_oi}, coor_both{array_oi}, edge_strengths_both{array_oi}, top_edges)

subplot(3,3,8)

plot_map(coordinates{array_oi}, coor_c{array_oi}, edge_strengths_c{array_oi}, top_edges)

subplot(3,3,9)

plot_map(coordinates{array_oi}, coor_d{array_oi}, edge_strengths_d{array_oi}, top_edges)


% if save_path is included in function call
if exist("save_path", "var")

    "saving"
    
    % save figure
    print(save_path, "-dpng", "-r950")
    
end % if


    % % nested functions
    function plot_map(all_coordinates, group_coordinates, edge_strengths, no_pairs)
        
        WIDTH_MULTIPLY = 5;
        SUBTRACT = 0.3499;
        
        % sort all strengths from high to low
        sorted_strenghts = sort(edge_strengths, "descend");
        
        % obtain threshold based on no_pairs to plot
        threshold = sorted_strenghts(abs(no_pairs));
        
        hold on
            
        % for each channel
        for channel = 1:size(all_coordinates,1)
            
            % plot channel location
            plot(all_coordinates(channel,1), all_coordinates(channel,2), "o", ...
                "MarkerEdgeColor", [0.5 0.5 0.5], "MarkerFaceColor", [0.5 0.5 0.5], ...
                "MarkerSize", 2)
            
        end % channel
        
        % for each edge
        for edge = 1:length(edge_strengths)
            
            % if edge lies below threshold
            if edge_strengths(edge) < threshold
                
                % plot line between ref_channel and comp_channel with thickness of tile_coef
                plot([group_coordinates(edge,1),group_coordinates(edge,3)],...
                    [group_coordinates(edge,2),group_coordinates(edge,4)],...
                    "-", "Color", [0.5 0.5 0.5 0.5], "LineWidth", ...
                    WIDTH_MULTIPLY*(0.4-SUBTRACT))
                
            end % if
            
        end % edge
        
        % for each edge
        for edge = 1:length(edge_strengths)
            
            % if edge lies above threshold
            if edge_strengths(edge) >= threshold
                
                % plot line between ref_channel and comp_channel with thickness of
                % tile_coef
                plot([group_coordinates(edge,1),group_coordinates(edge,3)],...
                    [group_coordinates(edge,2),group_coordinates(edge,4)],...
                    "-", "Color", [1 0 0 0.5], "LineWidth", ...
                    WIDTH_MULTIPLY*((edge_strengths(edge)-SUBTRACT)))
                
            end % if
            
        end % edge
        
        % plot lines to be used for legend
        leg_line1 = plot([-1000,-1000],[-950,-950], "-", "Color", [1 0 0 0.5], ...
            "LineWidth", WIDTH_MULTIPLY*(0.6-SUBTRACT));
        leg_line2 = plot([-1000,-1000],[-950,-950], "-", "Color", [1 0 0 0.5], ...
            "LineWidth", WIDTH_MULTIPLY*(0.8-SUBTRACT));
        leg_line3 = plot([-1000,-1000],[-950,-950], "-", "Color", [1 0 0 0.5], ...
            "LineWidth", WIDTH_MULTIPLY*(1-SUBTRACT));
             
        % add scale bar
        plot([0,500],[2100,2100], "k", "LineWidth", 5)
        
        % add legend
        Hleg = legend([leg_line1, leg_line2, leg_line3], {"0.6", "0.8", "1"}, ...
            "Location", "SouthOutside", "Orientation", "horizontal");
        
        % adjust axis
        axis equal
        set(gca, 'YDir','reverse')
        xlim([-100, 4000])
        ylim([-100, 2200])
        xticks([])
        yticks([])
        ax = gca;
        ax.FontSize = 14;
          
        hold off
        
    end % fun plot_connectivity_map

    % % %
    
    function [fits, xfits, lh] = plot_edge_dist(all_edges, array_oi, bin_edges, modelfun, symbols, colors)
        
        % make empty array for legend handles
        lh = NaN(1,size(all_edges,1));
        
        % make empty cell array for fitted models
        fits = cell(size(all_edges));
        
        % define bin centers
        bin_centers = bin_edges(2:end)-0.5*diff(bin_edges);
        
        hold on
        
        % for each array
        for array = 1:size(all_edges, 2)
            
            if ~isempty(all_edges{1,array})
                
                % for each treatment
                for treatment = 1:size(all_edges,1)
                    
                    % determine initial parameter estimates
                    est_initial = [((std(all_edges{treatment,array}))^2)/ ...
                        mean(all_edges{treatment,array}), (std(all_edges{treatment, ...
                        array})/(((std(all_edges{treatment,array}))^2)/mean( ...
                        all_edges{treatment,array})))^2, 1];
                    
                    % obtain hist counts
                    counts = histcounts(all_edges{treatment,array}, bin_edges, 'Normalization','probability');
                    
                    % combine x and y values in table
                    tbl = table(bin_centers', counts');
                    
                    try_bool = true;
                    
                    % while model is not properly fitted
                    while try_bool == true
                        try
                            % fit model
                            mdl = fitnlm(tbl, modelfun, est_initial);
                            try_bool = false;
                        catch
                            est_initial = est_initial*2;
                        end
                        
                    end % while try_bool
                    
                    % Extract the coefficient values from the the model object.
                    coefficients = mdl.Coefficients{:, 'Estimate'};
                    
                    % get x and y values for fitted model
                    xFitted = linspace(min(bin_centers), max(bin_centers), 301)';
                    yFitted = modelfun(coefficients,xFitted);
                    
                    if array == array_oi
                        
                        % plot results
                        lh(treatment) = scatter(bin_centers, counts, 50, symbols(treatment), ...
                            "MarkerEdgeColor",  colors(treatment), "LineWidth", 2); 
                        
                        % plot fitted model
                        plot(xFitted, yFitted, "Color", colors(treatment), "LineWidth", 1, "LineStyle", "--"); 
                        
                    end % if
                    
                    % store results of fitted model
                    fits{treatment, array} = yFitted;
                    xfits = xFitted;
                    
                end % treatment
                
            end % if
            
        end % array
        
        % adjust figure layout
        xlim([bin_edges(1),bin_edges(end)])
        
        xlabel("STTC")
        ylabel("Fraction of edges")

        ax = gca;
        ax.FontSize = 14;
        
    end % fun plot_edge_dist

    % % %
    
    function plot_edge_dist_frac_change(all_edges, fits, xfits, bin_edges, line_color)
        
        hold on
        
        % for each treatment
        for array = 1:size(all_edges,2)
            
            if ~isempty(fits{1,array})
                
                % plot results
                plot(xfits, (fits{2,array}-fits{1,array})./(fits{1,array}+fits{2,array}), ...
                    "LineWidth", 2, "Color", line_color{array});
                
            end % if
            
        end % treatment
        
        % add line at 0
        yline(0, "k--");
        
        % adjust figure layout
        ylim([-1,1])
        xlim([bin_edges(1),bin_edges(end)])
        
        xlabel("STTC")
        ylabel(["Fractional change", "(D-C)/(D+C)"])

        ax = gca;
        ax.FontSize = 14;
        
    end % fun plot_edge_dist_frac_change

    % % %

    function plot_group_count_bar(edge_strengths_both, edge_strengths_c, edge_strengths_d)
        
        % make empty result array
        bar_results = zeros(3, length(edge_strengths_both));
        
        % for each array
        for array = 1:length(edge_strengths_both)
            
            % compute number of edges per group
            num_both = size(edge_strengths_both{array}, 1);
            num_c = size(edge_strengths_c{array}, 1);
            num_d = size(edge_strengths_d{array}, 1);
        
            % store fraction of edges per group in bar_results array
            bar_results(1, array) = num_both / sum([num_both, num_c, num_d]);
            bar_results(2, array) = num_c / sum([num_both, num_c, num_d]);
            bar_results(3, array) = num_d / sum([num_both, num_c, num_d]);

        end % treatment
        
        % plot barplot
        hb = bar(bar_results', "stacked");
        
        set(hb,{'FaceColor'},{[0.4940 0.1840 0.5560] ; [0,0,1] ; [1,0,0]});

        xticklabels(["L1", "L2", "L3"])
        ylabel("Fraction of edges")

        ax = gca;
        ax.FontSize = 14;
        
    end % fun plot_group_count_bar

    % % %
    
    function [edge_strengths_both, edge_strengths_c, edge_strengths_d, ...
            coor_both, coor_c, coor_d] = process_data(edges_mat, coordinates)
        
        % make empty result cell arrays
        edge_strengths_both = cell(1,size(edges_mat,2));
        edge_strengths_c = cell(1,size(edges_mat,2));
        edge_strengths_d = cell(1,size(edges_mat,2));
        
        coor_both = cell(1,size(edges_mat,2));
        coor_c = cell(1,size(edges_mat,2));
        coor_d = cell(1,size(edges_mat,2));
        
        % for each array
        for array = 1:size(edges_mat,2)
            
            % initiate counters
            counter_both = 1;
            counter_c = 1;
            counter_d = 1;
            
            strength_array_both = [];
            strength_array_c = [];
            strength_array_d = [];
            coor_array_both = [];
            coor_array_c = [];
            coor_array_d = [];
            
            % for each channel row
            for row_chan = 1:size(edges_mat{1,array}, 1)
                
                % for each channel column
                for col_chan = 1:size(edges_mat{1,array}, 1)
                    
                    % if there is an edge for control and diazepam
                    if edges_mat{1,array}(row_chan, col_chan) >= 0.35 && edges_mat{2,array}(row_chan, col_chan) >= 0.35
                        
                        % store edge value
                        strength_array_both(counter_both, 1) = edges_mat{1,array}(row_chan, col_chan);
                        strength_array_both(counter_both, 2) = edges_mat{2,array}(row_chan, col_chan);
                        
                        coor_array_both(counter_both, [1,2]) = coordinates{array}(row_chan, :);
                        coor_array_both(counter_both, [3,4]) = coordinates{array}(col_chan, :);
                        
                        % add one to counter_c
                        counter_both = counter_both + 1;
                        
                        % if there is only an edge for control
                    elseif edges_mat{1,array}(row_chan, col_chan) >= 0.35
                        
                        % store edge value
                        strength_array_c(counter_c, 1) = edges_mat{1,array}(row_chan, col_chan);
                        
                        coor_array_c(counter_c, [1,2]) = coordinates{array}(row_chan, :);
                        coor_array_c(counter_c, [3,4]) = coordinates{array}(col_chan, :);
                        
                        % add one to counter_c
                        counter_c = counter_c + 1;
                        
                        % if there is only an edge for diazepam
                    elseif edges_mat{2,array}(row_chan, col_chan) >= 0.35
                        
                        % store edge value
                        strength_array_d(counter_d, 1) = edges_mat{2,array}(row_chan, col_chan);
                        
                        coor_array_d(counter_d, [1,2]) = coordinates{array}(row_chan, :);
                        coor_array_d(counter_d, [3,4]) = coordinates{array}(col_chan, :);
                        
                        % add one to counter_c
                        counter_d = counter_d + 1;
                        
                    end % if
                    
                end % col_chan
            end % row_chan
            
            % store results in cell for array
            edge_strengths_both{array} = strength_array_both;
            edge_strengths_c{array} = strength_array_c;
            edge_strengths_d{array} = strength_array_d;
            
            coor_both{array} = coor_array_both;
            coor_c{array} = coor_array_c;
            coor_d{array} = coor_array_d;
            
        end % array
        
    end % fun process_data

end % fun plot_fig6
