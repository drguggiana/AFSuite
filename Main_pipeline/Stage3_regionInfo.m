%% AF folder

%this script will grab the files from the snr script that you select (so
%from stage 2) and then split the traces across different AF folders. This
%is specified in the list of AFs below. It also draws info from the AF
%files obtained with the GUI, so at this stage you need that done before
%running this script.
%% clean up
clearvars
close all
Paths
%% Load the files and define paths

% assemble the overall path
main_path = fullfile(analysis_path,'Meta_files');

% get the folders in there
group_folders = dir(main_path);
group_folders = group_folders(3:end);
% get the number of subfolders
num_subfolders = length(group_folders);

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
        if contains(tar_path_temp(folders).name,{'4dpf','Fish'})
            continue
        else
            name_cell{folders} = fullfile(tar_path_temp(folders).folder,tar_path_temp(folders).name);
        end
    end
    % remove the empties
    name_cell = name_cell(~cellfun(@isempty,name_cell));
    
    % for now, filter the list of the names for the 4dpf fish
    num_data = length(name_cell); 
    %% Define the AFs in each dataset

    %load the AF list from the pipeline constants
    % constants_path = 'E:\Behavioral data\Matlab\AF_proc\Clemens_suite\20170827_Software_pipeline\subfunctions and scripts\pipeline_constants.mat';
    af_list = load(constants_path,'af_list');
    af_list = af_list.af_list;
    %% Load the AF files

    %the structure of the files is that per each fish there are several files
    %corresponding to the different imaging sessions. Each one of these files
    %contains an nxm cell where n is z and m is AF. Within this cell there are
    %1x2 cells that contain a px2 vector with the coordinates of the
    %corresponding ROI and a hxw logical matrix containing the ROI in mask form

    %allocate memory to store the matrices
    af_cell = cell(num_data,1);

    %the fish names
    fish_name = cell(num_data,1);
    %and to store the list cells
    af_sublist = cell(num_data,1);

    %for all the fish
    for fish = 1:num_data


        %load the file name of the fish to find out which folders to select
        [~,f_name,~] = fileparts(name_cell{fish});
        %split the name to get the folder names (fish and pre/post)
        f_parts = strsplit(f_name,'_');
        % edit the name to match the folders
        if strcmp(f_parts{1},'postcontrol') == 1
            f_parts{1} = 'control';
        elseif strcmp(f_parts{1},'precontrol') == 1
            f_parts{1} = 'pre';
        end
            
        %assemble the path to load files from
        f_path = fullfile(af_path,f_parts{2},f_parts{1});
        %get the files in the path
        af_files = dir(fullfile(f_path,'*.mat'));
        af_files = {af_files.name};
        %get the number of files
        file_num = size(af_files,2);
        %allocate memory for the file contents
        af_subcell = cell(file_num,1);
        %for all the files
        for files = 1:file_num
            %load the file
            % Clemens changed those 2 lines of code, so it will only load mat
            % files and not any tif
            af_subcell{files} = load(fullfile(f_path,af_files{files}),'a_store');

            af_subcell{files} = af_subcell{files}.a_store;
        end

        %also search for the list cell corresponding to this fish
        fish_name{fish} = strcat(f_parts{2},'_',f_parts{1});
        %extract just the names from the list
        %get the number of items in the list
        i_number = size(af_list,1);
        %allocate memory for just the names
        i_names = cell(i_number,1);
        %for all the fish
        for list_i = 1:i_number
            i_names{list_i} = af_list{list_i}{1};
        end
        list_c = strcmp(fish_name{fish},i_names);

        %store the af info in the cell
        af_sublist{fish} = af_list{list_c}{2};
        %concatenate and add to the main cell
        af_cell{fish} = af_subcell;
    end
    %% Combine the AF and seed coordinate info to yield a matrix with AF and z

    %allocate memory for the matrices per fish
    af_fish = cell(num_data,1);

    %for all of the fish
    for fish = 1:num_data
        
        fprintf(strjoin({'Current group:',num2str(group),'current fish:',num2str(fish),...
            'out of',num2str(num_data),'\r\n'},'_'))
        % get the data for this fish
        main_str = load(name_cell{fish},'main_str');
        main_str = main_str.main_str;
        % get the seed info for this fish
        seed_data = main_str.seed_data;
        %assemble a volume with the labeled regions for each fish

        %get the number of files
        file_num = size(af_cell{fish},1);
        %allocate memory to store the temporary matrix
        af_filescell = cell(file_num,1);

        %for each file
        for files = 1:file_num
            %determine the number of z sections
            z_num = size(af_cell{fish}{files},1);
            %also determine the number of AFs
            af_num = size(af_cell{fish}{files},2);
            %for all the z sections
            for z = 1:z_num
                %if the cell is empty, skip to the next iterations
                if isempty(af_cell{fish}{files}{z})
                    continue
                else
                    %determine the dimensions of the image
                    im_height = size(af_cell{fish}{files}{z}{2},1);
                    im_width = size(af_cell{fish}{files}{z}{2},2);
                    break
                end
            end
            %allocate memory for the volume
            af_vol = zeros(im_height,im_width,z_num);
            %for every z section
            for z = 1:z_num
                %before performing the assignment of the traces, correct for AF
                %overlap
                %allocate memory for a temp frame
                temp_frame = zeros(im_height,im_width);
                %for all the AFs
                for afc = 1:af_num
                    %if there are no ROIs there
                    if isempty(af_cell{fish}{files}{z,afc})
                        continue
                    end
                    %accumulate the frames, so I can find the places with
                    %overlap
                    temp_frame = temp_frame + af_cell{fish}{files}{z,afc}{2};
                end

                %find the places with overlap
                overlap_vec = find(temp_frame>1);
                %if there is overlap
                if ~isempty(overlap_vec)
                    %for all the AFs
                    for afc = 1:af_num
                         %if there are no ROIs there
                        if isempty(af_cell{fish}{files}{z,afc})
                            continue
                        end
                        %eliminate those voxels from every af
                        af_cell{fish}{files}{z,afc}{2}(overlap_vec) = 0;
                    end
                end
                %now actually assign the voxels to the respective AFs
                %for all the AFs
                for afc = 1:af_num
                    %if there are no ROIs there
                    if isempty(af_cell{fish}{files}{z,afc})
                        continue
                    end
                    %load the volume
                    af_vol(:,:,z) = af_vol(:,:,z) + af_sublist{fish}{files}(afc).*af_cell{fish}{files}{z,afc}{2};
                end
            end
    
            % get the seed_data for this file
            seed_file = seed_data(vertcat(seed_data.experiment)==files);
            % turn the coordinates into a seed_num x3 matrix
            seed_3d = round(horzcat(vertcat(seed_file.centroid),vertcat(seed_file.z)));
            % get the seed number
            seed_num = size(seed_3d,1);

            %allocate memory for the output matrix
            af_filemat = zeros(seed_num,2);
            %load the z
            af_filemat(:,2) = seed_3d(:,3);

            %and index the corresponding volume to get the AF
            %for all the seeds
            for seeds = 1:seed_num
                %load the points
                p = seed_3d(seeds,:);
                %and load the AF for the point
                af_filemat(seeds,1) = af_vol(p(2),p(1),p(3));
                %if the AF is zero, scan the surroundings of the points
                if af_filemat(seeds,1) == 0
                    %allocate a vector to contain the surroundings
                    val_vec = zeros(9,1);
                    vec_c = 1;
                    for x = -1:1
                        for y = -1:1
                            val_vec(vec_c) = af_vol(p(2)+x,p(1)+y,p(3));
                            vec_c = vec_c + 1;
                        end
                    end
                    %define the target AF as the most common for the seed
                    af_filemat(seeds,1) = mode(val_vec(val_vec~=0));
                end

            end
            %store the matrix
            af_filescell{files} = af_filemat;
        end

        %concatenate all the files and store in the final output cell
        af_fish{fish} = cat(1,af_filescell{:});
        %% Add the AF data to the meta data and save
        
        % write the field in the structure
        main_str.AF_info = af_fish{fish};
        % save the file where it was
        save(name_cell{fish},'main_str');
    end
end