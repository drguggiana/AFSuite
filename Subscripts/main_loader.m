function [main_cell,conc_trace_all,fish_ori_all,num_fish] = main_loader(group_path,flag)
% Load the main (and alternatively the traces) from the given folder

% get the files in the subfolder
% tar_path_temp = dir(fullfile(group_folders(group).folder,group_folders(group).name,'*.mat'));
tar_path_temp = dir(fullfile(group_path.folder,group_path.name,'*.mat'));
% get the number of experiments
num_fish = length(tar_path_temp);
% allocate memory for the full paths
name_cell = cell(num_fish,1);
% turn into a cell
for folders = 1:num_fish
    name_cell{folders} = fullfile(tar_path_temp(folders).folder,tar_path_temp(folders).name);
end
% remove the empties
name_cell = name_cell(~cellfun(@isempty,name_cell));

% for now, filter the list of the names for the 4dpf fish
num_fish = length(name_cell);

% allocate memory for the traces
conc_trace_all = cell(num_fish,1);
% allocate memory for the origins
fish_ori_all = cell(num_fish,1);
% allocate memory to store the main cells
main_cell = cell(num_fish,1);
% for all the fish in the group
for fish = 1:num_fish
    fprintf(strjoin({'Current fish:',num2str(fish),...
        'out of',num2str(num_fish),'\r\n'},'_'))
    % get the data for this fish
    main_str = load(name_cell{fish},'main_str');
    main_str = main_str.main_str;
    % store the structure
    main_cell{fish} = main_str;
    % get the dimensions of the trace file
    seed_num = size(main_str.fish_ori,1);
    stim_num = main_str.stim_num;
    time_num = main_str.time_num;
    
    % if the flag is 1, also load the traces
    if flag == 1
        % load the traces
        temp_path = main_str.path_trial_ave;
        fID = fopen(temp_path);
        temp_data = fread(fID,'double');
        fclose(fID);
        % save in the cell
        conc_trace_all{fish} = reshape(temp_data,seed_num,time_num,stim_num);
    elseif flag == 2
        % load the traces
        temp_path = main_str.path_single_reps;
        fID = fopen(temp_path);
        temp_data = fread(fID,'double');
        fclose(fID);
        % save in the cell
        conc_trace_all{fish} = reshape(temp_data,seed_num,time_num,stim_num,main_str.rep_num);
    else
        conc_trace_all{fish} = [];
    end
    % save the fish ori
    fish_ori_all{fish} = ones(seed_num,1).*fish;
end
% concatenate the fish and the ori
conc_trace_all = cat(1,conc_trace_all{:});
fish_ori_all = cat(1,fish_ori_all{:});