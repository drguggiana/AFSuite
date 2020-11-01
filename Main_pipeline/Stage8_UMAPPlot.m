%% clean up
clearvars
close all
Paths
%% Load the files and define paths

% define the groups to cluster
target_groups = {'pre','postcontrol'};

% assemble the overall path
main_path = fullfile(analysis_path,'Meta_files');

% % load the paths
% [group_folders,num_subfolders] = group_loader(target_groups,main_path);

% load the constants
constants = load(constants_path,'constants');
constants = constants.constants;

% load the labels
labels = load(constants_path,'labels');
labels = labels.labels;

af_labels = labels.af;
celltype_labels = labels.celltype;

% define the figure saving path
figure_path = fullfile(fig_path,'UMAP');
%% load the average structure (SLOW)
average_str = load(meta_path);
average_str = average_str.average_str;
%% Plot the regions on the map
close all

% define the colormap
cmap = jet(length(af_labels));
% turn the af labels into a vector
af_numbers = str2double({af_labels.number});
% allocate a vector for the labels
af_vector = zeros(max(af_numbers),1);
% for all the labels
for labels = 1:length(af_numbers)
    af_vector(af_numbers(labels)) = labels;
end
% for all the groups
for group = 1:length(target_groups)
    
    % get the region info
    region_all = cat(1,average_str.(target_groups{group}).meta.AF_info);
    celltype_all = cat(1,average_str.(target_groups{group}).meta.cell_type);
    % for all the celltypes
    for celltype = 1:3
        figure
        % get the reduced data
        reduced_data = average_str.(target_groups{group}).reduced_data{celltype};
        
        % get the regions
        region_celltype = region_all(celltype_all==(celltype-1),1);
        
        % remove the nans
        nan_vector = ~isnan(region_celltype);
        reduced_data = reduced_data(nan_vector,:);
        region_celltype = region_celltype(nan_vector);
        
        % subsample
        % define the subsampling factor
        subsample_factor = 20;
        reduced_data = reduced_data(1:subsample_factor:end,:);
        region_celltype = region_celltype(1:subsample_factor:end);
        
        % plot
        gscatter(reduced_data(:,1), reduced_data(:,2),af_vector(region_celltype),cmap,'.',10)
%         scatter(reduced_data(:,1), reduced_data(:,2),5,region_celltype,'.')
%         colormap(colorcube)
        axis square
        % get the current regions
        current_regions = unique(af_vector(region_celltype));
        legend({af_labels(current_regions).name},'Interpreter','None','Location','bestoutside')
       
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = figure_path;
        fig_set(1).fig_name = strjoin({'Regions',target_groups{group},celltype_labels{celltype},'.png'},'_');
        fig_set(1).fig_size = 5;
        fig_set(1).Title = strjoin({target_groups{group},celltype_labels{celltype}},'_');
%         fig_set(1).cmap = hsv;
        
        h = style_figure(gcf,fig_set);
    end
end

autoArrangeFigures
%% Plot the fish
close all
% for all the groups
for group = 1:length(target_groups)
    
    % get the number of fish
    fish_num = size(average_str.(target_groups{group}).meta,1);
    % allocate memory for the fish vector
    fish_all = cell(fish_num,1);
    % define the colormap
    cmap = jet(fish_num);
    % for all the fish
    for fish = 1:fish_num
        fish_all{fish} = zeros(size(average_str.(target_groups{group}).meta(fish).cell_type))+fish;
    end
    % concatenate to assemble
    fish_all = vertcat(fish_all{:});
    % get the region info
%     fish_all = cat(1,average_str.(target_groups{group}).meta.AF_info);
    celltype_all = cat(1,average_str.(target_groups{group}).meta.cell_type);
    % for all the celltypes
    for celltype = 1:3
        figure
        
        % get the reduced data
        reduced_data = average_str.(target_groups{group}).reduced_data{celltype};
        
        % get the regions
        fish_celltype = fish_all(celltype_all==(celltype-1),1);
        % subsample
        % define the subsampling factor
        subsample_factor = 20;
        reduced_data = reduced_data(1:subsample_factor:end,:);
        fish_celltype = fish_celltype(1:subsample_factor:end,:);

        % plot
        gscatter(reduced_data(:,1), reduced_data(:,2),fish_celltype,cmap,'.',10)
%         scatter(reduced_data(:,1), reduced_data(:,2),10,fish_celltype,'.')
        colormap(hsv)
        axis square
        
        legend('Location','bestoutside')
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = figure_path;
        fig_set(1).fig_name = strjoin({'Fish',target_groups{group},celltype_labels{celltype},'.png'},'_');
        fig_set(1).fig_size = 5;
        fig_set(1).Title = strjoin({target_groups{group},celltype_labels{celltype}},'_');
%         fig_set(1).cmap = hsv;
        
        h = style_figure(gcf,fig_set);
    end
end

autoArrangeFigures
%% Plot the clusters
close all

% for all the groups
for group = 1:length(target_groups)
    
    % for all the celltypes
    for celltype = 1:3
        figure
        % get the reduced data
        reduced_data = average_str.(target_groups{group}).reduced_data{celltype};
        
        % get the clusters
        % get the number of fish
        fish_num = size(average_str.(target_groups{group}).meta,1);
        % allocate memory for the cluster idx
        idx_celltype = cell(fish_num,1);
        % for all the fish
        for fish = 1:fish_num
            idx_celltype{fish} = average_str.(target_groups{group}).meta(fish).idx_cell{celltype};
        end
        % concatenate the indexes
        idx_celltype = vertcat(idx_celltype{:});
        
        % define the colormap
        cmap = jet(length(unique(idx_celltype)));
        % subsample
        % define the subsampling factor
        subsample_factor = 20;
        reduced_data = reduced_data(1:subsample_factor:end,:);
        idx_celltype = idx_celltype(1:subsample_factor:end,:);
        
        % plot
        gscatter(reduced_data(:,1), reduced_data(:,2),idx_celltype,cmap,'.',10)
%         scatter(reduced_data(:,1), reduced_data(:,2),30,idx_celltype,'.')
        colormap(hsv)
        axis square
        
        legend('Location','bestoutside')
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = figure_path;
        fig_set(1).fig_name = strjoin({'Cluster',target_groups{group},celltype_labels{celltype},'.png'},'_');
        fig_set(1).fig_size = 5;
        fig_set(1).Title = strjoin({target_groups{group},celltype_labels{celltype}},'_');
%         fig_set(1).cmap = hsv;
        
        h = style_figure(gcf,fig_set);
    end
end

autoArrangeFigures
%% Plot the on/off index

close all

% define the colormap
cmap = parula(2);
% for all the groups
for group = 1:length(target_groups)
    
    % for all the celltypes
    for celltype = 1:3
        figure
        % get the reduced data
        reduced_data = average_str.(target_groups{group}).reduced_data{celltype};
        
        % get the onoff
        % get the number of fish
        fish_num = size(average_str.(target_groups{group}).meta,1);
        % allocate memory for the cluster idx
        onoff_celltype = cell(fish_num,1);
        % for all the fish
        for fish = 1:fish_num
            onoff_celltype{fish} = average_str.(target_groups{group}).meta(fish).onoff{celltype};
        end
        % concatenate the indexes
        onoff_celltype = vertcat(onoff_celltype{:});
        
        % subsample
        % define the subsampling factor
        subsample_factor = 20;
        reduced_data = reduced_data(1:subsample_factor:end,:);
        onoff_celltype = onoff_celltype(1:subsample_factor:end,:);
        
        % plot
        gscatter(reduced_data(:,1), reduced_data(:,2),onoff_celltype,cmap,'.',10)
%         scatter(reduced_data(:,1), reduced_data(:,2),30,onoff_celltype,'.')
%         colormap(viridis)
        axis square
        legend({'Off','On'},'Location','bestoutside')
        
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = figure_path;
        fig_set(1).fig_name = strjoin({'OnOff',target_groups{group},celltype_labels{celltype},'.png'},'_');
        fig_set(1).fig_size = 5;
        fig_set(1).Title = strjoin({target_groups{group},celltype_labels{celltype}},'_');
        fig_set(1).cmap = parula;
        
        h = style_figure(gcf,fig_set);
    end
end

autoArrangeFigures
%% Plot the pref_dir index

close all
% define the colormap
cmap = hsv(8);
% for all the groups
for group = 1:length(target_groups)
    
    % for all the celltypes
    for celltype = 1:3
        figure
        % get the reduced data
        reduced_data = average_str.(target_groups{group}).reduced_data{celltype};
        
        % get the onoff
        % get the number of fish
        fish_num = size(average_str.(target_groups{group}).meta,1);
        % allocate memory for the cluster idx
        prefdir_celltype = cell(fish_num,1);
        % for all the fish
        for fish = 1:fish_num
            prefdir_celltype{fish} = average_str.(target_groups{group}).meta(fish).prefdir{celltype};
        end
        % concatenate the indexes
        prefdir_celltype = vertcat(prefdir_celltype{:});
        
        % subsample
        % define the subsampling factor
        subsample_factor = 20;
        reduced_data = reduced_data(1:subsample_factor:end,:);
        prefdir_celltype = prefdir_celltype(1:subsample_factor:end,:);
        % plot
        gscatter(reduced_data(:,1), reduced_data(:,2),prefdir_celltype,cmap,'.',10)
%         scatter(reduced_data(:,1), reduced_data(:,2),30,onoff_celltype,'.')
%         colormap(viridis)
        axis square
%         legend({'Off','On'},'Location','bestoutside')
        
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = figure_path;
        fig_set(1).fig_name = strjoin({'Prefdir',target_groups{group},celltype_labels{celltype},'.png'},'_');
        fig_set(1).fig_size = 5;
        fig_set(1).Title = strjoin({target_groups{group},celltype_labels{celltype}},'_');
        fig_set(1).cmap = parula;
        
        h = style_figure(gcf,fig_set);
    end
end

autoArrangeFigures
%% Plot the signal strength

close all

% for all the groups
for group = 1:length(target_groups)
    
    % for all the celltypes
    for celltype = 1:3
        figure
        % get the reduced data
        reduced_data = average_str.(target_groups{group}).reduced_data{celltype};
        
        % get the onoff
        % get the number of fish
        fish_num = size(average_str.(target_groups{group}).meta,1);
        % allocate memory for the cluster idx
        signal_celltype = cell(fish_num,1);
        % for all the fish
        for fish = 1:fish_num
            signal_celltype{fish} = average_str.(target_groups{group}).meta(fish).signal{celltype};
        end
        % concatenate the indexes
        signal_celltype = vertcat(signal_celltype{:});
        
        % subsample
        % define the subsampling factor
        subsample_factor = 20;
        reduced_data = reduced_data(1:subsample_factor:end,:);
        signal_celltype = signal_celltype(1:subsample_factor:end,:);
        % plot
%         gscatter(reduced_data(:,1), reduced_data(:,2),signal_celltype,cmap,'.',10)
        scatter(reduced_data(:,1), reduced_data(:,2),30,log(signal_celltype),'.')
%         colormap(viridis)
        axis square
%         legend({'Off','On'},'Location','bestoutside')
        
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = figure_path;
        fig_set(1).fig_name = strjoin({'Signal',target_groups{group},celltype_labels{celltype},'.png'},'_');
        fig_set(1).fig_size = 5;
        fig_set(1).Title = strjoin({target_groups{group},celltype_labels{celltype}},'_');
        fig_set(1).cmap = parula;
        
        h = style_figure(gcf,fig_set);
    end
end

autoArrangeFigures