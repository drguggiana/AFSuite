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

% define the groups to cluster
target_groups = {'precontrol'};

% assemble the overall path
main_path = fullfile(analysis_path,'Meta_files');

% load the paths
[group_folders,num_subfolders] = group_loader(target_groups,main_path);

% for all the subfolders
for group = 1:num_subfolders
    % load the traces
    [main_cell,conc_trace_all,fish_ori_all,num_fish] = main_loader(group_folders(group),1);

    %% Extract the singular values for activity of the directional stimuli
    % define the type of extraction
    stimTypeNum = 4;
    % get the SVD values
    [~,tconc_trace, stim_length, stim_num, dsi_clas_final] = fish_svd(conc_trace_all,stimTypeNum);
    
    % clean up memory
    clear conc_trace_all    
    %% Run the sPCA on selected data
    
    % restructure tconc_trace for the parallel loop
    feature_cell = cell(3,1);
    index_cell = cell(3,1);
    for celltype = 1:3
        feature_cell{celltype} = tconc_trace(dsi_clas_final==(celltype-1),1:sum(stim_length(1:end-1)));
        index_cell{celltype} = tconc_trace(dsi_clas_final==(celltype-1),sum(stim_length(1:end-1))+1:end);
    end
    % allocate memory for the clustering results
    gmm_str = struct([]);
    % create the pool of workers
    if group == 1
        worker_pool = parpool;
    end
    
    % for all 3 types of cells
    parfor celltype = 1:3
        % define (maybe load later) the cluster numbers to try
        clu_vec = [5 10 20 30 50 70 100];
  
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
        % define the number of replicates to run
        replicates = 10;
        % if the cluster number is not appropriate, remove it
        temp_clu_vec = clu_vec(clu_vec<size(index_cell{celltype},1));
        
        % SPCA and cluster
        [gmm_str(celltype).idx_clu,gmm_str(celltype).GMM,gmm_str(celltype).clu_num,~,gmm_str(celltype).bic] =...
            sPCA_GMM_cluster(feature_cell{celltype},bounds,K,t_bins,pca_vec,temp_clu_vec,replicates,4,...
            index_cell{celltype});

    end
    % if it's the last run, delete the pool
    if group == num_subfolders
        delete(worker_pool)
    end
    %% Save analysis output
    
    % update the main_str for each fish
    % for all the fish
    for fish = 1:num_fish
        % get the corresponding structure
        main_str = main_cell{fish};
        % get the type, indexes and cluster number corresponding to that fish
        main_str.cell_type = dsi_clas_final(fish_ori_all==fish);
        main_str.clu_num = vertcat(gmm_str.clu_num);
        
        % allocate memory to store the indexes
        temp_cell = cell(3,1);
        % for all the celltypes
        for celltype = 1:3
            % get the fish id vector for the type
            type_ori = fish_ori_all(dsi_clas_final==(celltype-1));
            temp_cell{celltype} = gmm_str(celltype).idx_clu(type_ori==fish);
        end
        % store in the main structure
        main_str.idx_cell = temp_cell;
        
        % save the structure
        save(name_cell{fish},'main_str');
    end
    % also save the clustering structure
    % assemble the name
    clustering_name = strjoin({'Clustering',target_groups{group}},'_');
    % save the file
    save(fullfile(clustering_path,clustering_name),'gmm_str')
    
end