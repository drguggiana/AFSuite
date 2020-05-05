%% clean up
clearvars
close all
Paths
%% Load the traces
% define the groups to cluster
target_groups = {'precontrol'};

% assemble the overall path
main_path = fullfile(analysis_path,'Meta_files');

% get the folders in there
group_folders = dir(main_path);
group_folders = group_folders(3:end);
% get the number of subfolders
num_subfolders = length(group_folders);
keep_vector = zeros(num_subfolders,1);
% for all the groups
for group = 1:num_subfolders
    if any(strcmp(target_groups,group_folders(group).name))
        keep_vector(group) = 1;
    else
        keep_vector(group) = 0;
    end
end
% select the subfolders
group_folders = group_folders(keep_vector==1);
% for all the subfolders
for group = 1:num_subfolders
    % get the files in the subfolder
    tar_path_temp = dir(fullfile(group_folders(group).folder,group_folders(group).name,'*.mat'));
    % get the number of experiments
    num_data = length(tar_path_temp);
    % allocate memory for the full paths
    name_cell = cell(num_data,1);
    % turn into a cell
    for folders = 1:num_data
        name_cell{folders} = fullfile(tar_path_temp(folders).folder,tar_path_temp(folders).name);
    end
    % remove the empties
    name_cell = name_cell(~cellfun(@isempty,name_cell));
    
    % for now, filter the list of the names for the 4dpf fish
    num_data = length(name_cell);
    
    % allocate memory for the traces
    conc_trace_all = cell(num_data,1);
    % allocate memory for the origins
    fish_ori_all = cell(num_data,1);
    % allocate memory to store the main cells
    main_cell = cell(num_data,1);
    % for all the fish in the group
    for fish = 1:num_data
        fprintf(strjoin({'Current group:',num2str(group),'current fish:',num2str(fish),...
            'out of',num2str(num_data),'\r\n'},'_'))
        % get the data for this fish
        main_str = load(name_cell{fish},'main_str');
        main_str = main_str.main_str;
        % store the structure
        main_cell{fish} = main_str;
        % get the dimensions of the trace file
        seed_num = size(main_str.fish_ori,1);
        stim_num = main_str.stim_num;
        time_num = main_str.time_num;
        
        % load the traces
        temp_path = main_str.path_trial_ave;
        fID = fopen(temp_path);
        temp_data = fread(fID,'double');
        fclose(fID);
        % save in the cell
        conc_trace_all{fish} = reshape(temp_data,seed_num,time_num,stim_num);
        % save the fish ori
        fish_ori_all{fish} = ones(seed_num,1).*fish;
    end
    % concatenate the fish and the ori
    conc_trace_all = cat(1,conc_trace_all{:});
    fish_ori_all = cat(1,fish_ori_all{:});
end