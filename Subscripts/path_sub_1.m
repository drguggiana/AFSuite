function [path_cell] = path_sub_1(tar_path)


%get the files and folders in this path
curr_list = dir(tar_path);
%get rid of the dots
curr_list = curr_list(3:end);
%get the number of directories in the list
num_subdir = sum(vertcat(curr_list(:).isdir));
%if there are directories in the list
if num_subdir > 0
    
    %isolate the folders
    subdir_list = curr_list(vertcat(curr_list(:).isdir)==1);
    %get the parent path
    parent_path = subdir_list(1).folder;
    %allocate memory to store their subdirs
    subdir_cell = cell(num_subdir,1);
    %for all of them
    for subdir = 1:num_subdir
        %assemble the path to the subfolder
        sub_path = strcat(parent_path,'\',subdir_list(subdir).name);
        %run this function recursively
        subdir_cell{subdir} = path_sub_1(sub_path);
    end
    %concatenate the sublists to generate the output
    path_cell = vertcat(subdir_cell{:});
else
    %output the last part of the input folder since there are no folders deeper
    path_cell = {tar_path};
end
