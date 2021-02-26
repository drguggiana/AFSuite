%% Quantify cluster relationships
clearvars
close all
Paths
%% Load the meta structure

average_str = load(metapost_path);
average_str = average_str.average_str;
%% Define some constants

% target_groups = fields(average_str);
target_groups = {'pre'};
num_subfolders = length(target_groups);

% load the labels
labels = load(constants_path,'labels');
labels = labels.labels;

af_labels = labels.af;
celltype_labels = labels.celltype;

% define the figure path
figure_path = fullfile(fig_path,'clusterAverages');

% get the total number of regions
num_regions = length(af_labels);
%% Count the occurrences


% get the number of fish
num_fish = size(average_str.pre.meta,1);

% allocate memory for all the fish
cluster_matrix = cell(num_fish,3);
% for all the fish
for fish = 1:num_fish

    % load the celltype info
    cell_info = average_str.pre.meta(fish).cell_type;
    % load the region info
    region_info = average_str.pre.meta(fish).AF_info;
    % for all celltypes
    for celltype = 1:3
        % get the number of clusters
        clu_num = average_str.pre.meta(fish).clu_num(celltype);
        % allocate memory for the cluster numbers
        cluster_by_area = zeros(clu_num,num_regions);
        
        % get the region info for this celltype
        cell_region = region_info(cell_info==celltype-1,1);
        % load the idx of the celltype
        idx = average_str.pre.meta(fish).idx_cell{celltype};
        
        % for all the regions
        for region = 1:num_regions
            % for all the clusters
            for clu = 1:clu_num
                % count the occurrences of this region
                cluster_by_area(clu,region) = ...
                    sum((idx==clu)&(cell_region==str2double(af_labels(region).number)));
            end
        end
        % save the matrix
        cluster_matrix{fish,celltype} = cluster_by_area;
        
    end
end
%% Average across fish and plot
close all
% for all the celltypes
for celltype = 1:3
    % average
    average_matrix = nanmean(cat(3,cluster_matrix{:,celltype}),3);
    % plot
    figure
    imagesc(log10((average_matrix)))
    title(celltype_labels{celltype})
    set(gca,'TickLength',[0 0],'XTick',1:num_regions,'XTickLabel',{af_labels.name},...
        'XTickLabelRotation',45,'TickLabelInterpreter','None')
    ylabel('Clusters')
    colorbar
end
autoArrangeFigures

%% Correlate the patterns

% allocate memory for the results
corr_matrix = cell(num_fish,3);
% for all the celltypes
for celltype = 1:3
    
    % for all the fish
    for fish = 1:num_fish
        % calculate the correlation matrix
        corr_matrix{fish,celltype} = corr(cluster_matrix{fish,celltype});
    end
end
%% Plot the averages

close all
% for all the celltypes
for celltype = 1:3
    % average
    average_matrix = nanmean(cat(3,corr_matrix{:,celltype}),3);
    % plot
    figure
    imagesc(average_matrix)
    title(celltype_labels{celltype})
    set(gca,'TickLength',[0 0],'XTick',1:num_regions,'XTickLabel',{af_labels.name},...
        'XTickLabelRotation',45,'TickLabelInterpreter','None')
    set(gca,'YTick',1:num_regions,'YTickLabel',{af_labels.name})
    axis square
    axis tight
    colorbar
end
autoArrangeFigures