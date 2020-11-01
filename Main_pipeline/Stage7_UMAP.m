%% Cluster each AF
% This one clusters all AFs, all fish and all conditions together

%this script takes the folder you select (as many AFs as you want) and
%clusters each one of them (including all the fish you put in the folder
%and that also match the keyword, either pre or post, you enter in the
%input dialog). The clustered files are saved for next stage.
%% clean up
clearvars
close all
Paths
%% Load the files and define paths

% initialize rng
rng(1)

% define the groups to cluster
target_groups = {'pre','postcontrol'};

% assemble the overall path
main_path = fullfile(analysis_path,'Meta_files');

% load the paths
[group_folders,num_subfolders] = group_loader(target_groups,main_path);

% load the constants
constants = load(constants_path,'constants');
constants = constants.constants;

% load the average structure
average_str = load(meta_path);
average_str = average_str.average_str;

% for all the subfolders
for group = 1:num_subfolders
    % load the traces
    [main_cell,conc_trace_all,fish_ori_all,num_fish,name_cell] = main_loader(group_folders(group),1);

    %% Extract the singular values for activity of the directional stimuli
    % define the type of extraction
    stimTypeNum = 4;
    % get the SVD values
    [~,tconc_trace, stim_length, stim_num, dsi_clas_final] = fish_svd(conc_trace_all,...
        stimTypeNum,1,cat(1,average_str.(target_groups{group}).meta.cell_type));
    
    % clean up memory
    clear conc_trace_all    
    %% Load the snr info
%     % allocate memory for the info
%     snr_info = cell(num_fish,1);
%     % for all the fish
%     for fish = 1:num_fish
%         snr_info{fish} = main_cell{fish}.snr_mat(main_cell{fish}.snr_vec,:);
%     end
%     % concatenate the values
%     snr_info = vertcat(snr_info{:});
    %% Run the sPCA and cluster the selected data
    
    % restructure tconc_trace for the parallel loop
    feature_cell = cell(3,1);
    index_cell = cell(3,1);
%     snr_cell = cell(3,1);
    for celltype = 1:3
        feature_cell{celltype} = tconc_trace(dsi_clas_final==(celltype-1),1:sum(stim_length(1:end-1)));
        index_cell{celltype} = tconc_trace(dsi_clas_final==(celltype-1),sum(stim_length(1:end-1))+1:end);
%         snr_cell{celltype} = snr_info(dsi_clas_final==(celltype-1),:);
    end

    % allocate a cell to store the umap data
    umap_cell = cell(3,1);
    reduced_cell = cell(3,1);
    % for all 3 types of cells (PARFOR GOES HERE)
    for celltype = 1:3
  
        % define the bounds of each stimulus
        bounds_top = [1 cumsum(stim_length(1:end-2))+1];
        bounds_bottom = cumsum(stim_length(1:end-1));
        bounds = [bounds_top;bounds_bottom];
        % define the number of components per stimulus
        K = ones(1,stim_num).*5;
        % define the time bins per sparse component
        t_bins = ones(stim_num,1).*10;
        % define which stimuli to evaluate with spca vs regular pca
        pca_vec = ones(stim_num,1).*1;
        
        % SPCA and cluster
        [~,~,~,pca_output,~,feature_matrix] =...
            sPCA_GMM_cluster(feature_cell{celltype},bounds,K,t_bins,pca_vec,[],[],2);
        
        % column normalize the indexes
        current_index = index_cell{celltype};
        current_index = (current_index - mean(current_index,2))./std(current_index,0,2);
        % concatenate the index data to the PCA results and run umap
        [reduced_cell{celltype},umap_cell{celltype}] = run_umap(horzcat(feature_matrix,current_index),...
            'min_dist',0.7,'n_neighbors',20);
        

    end
    %% Save analysis output
    
    % save the umap model and reduced coordinates in the structure
    average_str.(target_groups{group}).umap_model = umap_cell;
    average_str.(target_groups{group}).reduced_data = reduced_cell;
    
end
%% Save the structure file again
% define the save path
save_path = fullfile(analysis_path,'Meta_files','all_meta.mat');
save(save_path,'average_str','-v7.3')