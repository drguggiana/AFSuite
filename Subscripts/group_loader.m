function [group_folders,num_subfolders] = group_loader(target_groups,main_path)
% Load the paths to the individual fish on the groups given

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
% update the number of subfolders
num_subfolders = length(group_folders);