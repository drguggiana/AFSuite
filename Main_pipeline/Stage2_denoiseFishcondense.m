%% Remove low SNR traces and combine files for a single fish

%this script will take the output from the previous one (you gotta select
%the pre or post folder of the desired fish) and it'll combine all of the
%files for that fish, also removing the low snr traces
%% clean up
clearvars
close all
Paths
%% Load the files and define paths

%get the folder where the image files are
% tar_path_all_pre = uipickfiles('FilterSpec','E:\Behavioral data\Matlab\AF_proc\Clemens_suite\pipeline_test\Stage1\*.mat');
tar_path_all_pre = uipickfiles('FilterSpec',fullfile(analysis_path,'Stage1'));

% define the save path
save_path = fullfile(analysis_path,'Stage2');
% %get the file paths from the selected folder
% tar_path_all = path_sub_1(tar_path_all{1});
%allocate memory for the path list
tar_path_all = cell(size(tar_path_all_pre));

%for each one of the folders selected, get the subfolders
for subdir = 1:length(tar_path_all)
    tar_path_all{subdir} = path_sub_1(tar_path_all_pre{subdir});
    %get a vector with the dirs containing REP
    rep_dir = ~cellfun(@isempty,strfind(tar_path_all{subdir},'REP'));
    %if REP files are present, leave only those
    if sum(rep_dir) > 0
        tar_path_all{subdir} = tar_path_all{subdir}(rep_dir);
    end
end
%concatenate the full list of files
tar_path_all = vertcat(tar_path_all{:});

%get the number of experiments selected
num_exp = length(tar_path_all);

%define the list of labels to sort the files
label_list = {'_traces.mat'};
for exps = 1:num_exp
    % %get the number of each type of data (the round is to avoid the nonscalar
    % %warning for the for loop)
    % num_data = round(num_exp./length(label_list));
    %get a list of file names in the directory
    file_list = dir(strcat(tar_path_all{exps},'\*.mat'));
    file_list = {file_list.name};
    %get the number of trace files in the directory
    num_data = size(dir(strcat(tar_path_all{exps},'\*',label_list{1})),1);
    %allocate memory for the different types of files
    name_cell = cell(num_data,length(label_list));

    %for the types of files
    for f_type = 1:length(label_list)

        %get the coordinates of the file names
        name_map = file_list(~cellfun(@isempty,strfind(file_list,label_list{f_type})));
        %for all the names
        for names = 1:size(name_map,2)
            %store them in the corresponding layer of the name cell
            name_cell(names,f_type) = fullfile(tar_path_all(exps),name_map{names});
        end
    end
    %% Define/load constants

    %get the number of time points
    time_num = load(name_cell{1},'time_num');
    time_num = time_num.time_num;
    %get the number of stimuli
    stim_num2 = load(name_cell{1},'stim_num');
    stim_num2 = stim_num2.stim_num;
    %% Get the snr threshold

    %load the thresholding parameters from the constants file
%     constants_path = 'E:\Behavioral data\Matlab\AF_proc\Clemens_suite\20170827_Software_pipeline\subfunctions and scripts\pipeline_constants.mat';
    
    percentile_constant = load(constants_path,'percentile_constant');
    percentile_constant = percentile_constant.percentile_constant;
    stimulus_constant = load(constants_path,'stimulus_constant');
    stimulus_constant = stimulus_constant.stimulus_constant;
    
    
    %allocate memory to store the snr for the entire fish
    snr_all = cell(num_data,1);
    %for all the files
    for files = 1:num_data
        %load and concatenate the snr for this file
        snr_temp = load(name_cell{files,1},'snr_cell');
        snr_all{files} = cat(1,snr_temp.snr_cell{:});
    end
    %concatenate the snr for all of the files
    snr_mat = cat(1,snr_all{:});

    %allocate memory for storing the thresholds
    thres_vec = zeros(stim_num2,1);

    %for all of the stimuli
    for stim = 1:stim_num2
        %calculate the 25th percentile and store as the threshold
        thres_vec(stim) = prctile(snr_mat(:,stim),percentile_constant);
    end

    % figure
    % plot(thres_vec)
    %% Load the traces 
    
    %define a list of the stimuli to keep
    % keep_stim = [1 3 13 14 10 12];
    keep_stim = [];
    %load and concatenate all the traces

    %allocate memory for the actual traces
    conc_trace = cell(num_data,1);
    % allocate memory for the single reps
    reps_trace = cell(num_data,1);
    % allocate memory for the seed coordinates and area
    seed_data = cell(num_data,1);
    % save the number of seeds per fish also
    seeds_per_fish = cell(num_data,1);
    %for all the traces
    for fish = 1:num_data
        %load the conc_trace variable into the trace_all cell
        conc_trace{fish} = load(name_cell{fish,1},'conc_trace');
        %extract the matrix of values from the struct
        conc_trace{fish} = conc_trace{fish}.conc_trace;
        %load the conc_trace variable into the trace_all cell
        reps_trace{fish} = load(name_cell{fish,1},'reps_trace');
        %extract the matrix of values from the struct
        reps_trace{fish} = reps_trace{fish}.reps_trace;
        % save the seeds per fish
        seeds_per_fish{fish} = size(conc_trace{fish},1);
        
        % load the coordinates
        seed_concat = load(name_cell{fish,1},'seed_concat');
        seed_concat = seed_concat.seed_concat;
        
        z_seed = load(name_cell{fish,1},'z_seed');
        z_seed = z_seed.z_seed;
        
        if ~isempty(keep_stim)
            %load a temp matrix for calculations
            temp_mat = trace_all{fish};
            %allocate memory for the stimuli to keep
            keep_mat = zeros(size(temp_mat,1),time_num*length(keep_stim));

            %initialize a counter for indexing the matrix
            keepc = 1;
            %for all the stimuli to keep
            for kstim = keep_stim
                %store the desired stim in the matrix
                keep_mat(:,keepc:keepc+time_num-1) = ...
                    temp_mat(:,(kstim-1)*time_num+1:kstim*time_num);
                %update the counter
                keepc = keepc + time_num;
            end
            %load the modified matrix into the original storage
            conc_trace{fish} = keep_mat;

            %if it's the first fish
            if fish == 1
                %redefine the number of stimuli
                stim_num2 = length(keep_stim);
                %also modify the color matrix
                col_out = col_out(keep_stim,:,:);
            end
        end
    end

    %concatenate the responses
    conc_trace = vertcat(conc_trace{:});
    reps_trace = vertcat(reps_trace{:});
    %and create a vector with the fish of origin for each seed and the original
    %seed number
    fish_ori = zeros(size(conc_trace,1),2);
    %initialize a counter for the index within the seeds
    fish_count = 1;
    %for all the experiments
    for fish = 1:num_data
        %get the number of seeds in the fish
%         seed_num = size(trace_all{fish},1);
        seed_num = seeds_per_fish{fish};

        %insert it in the final vector
        fish_ori(fish_count:fish_count+seed_num-1,1) = fish;
        %also insert the number of each seed
        fish_ori(fish_count:fish_count+seed_num-1,2) = 1:seed_num;
        %update the counter
        fish_count = seed_num + fish_count;
    end
    %% Threshold the traces based on the snr calculation above

    %the idea is to exclude any trace that has a value under the threshold,
    %under any stimulus

    %turn the NaNs (probably from traces with zeros across the board)into 0
    snr_mat(isnan(snr_mat)) = 0;
    %get the binary version of the traces
    snr_bin = bsxfun(@gt,snr_mat,thres_vec');

    thres_res = zeros(stim_num2,1);
    for t = 1:stim_num2
        %to do this, AND the columns of the snr matrix
        snr_vec = logical(sum(snr_bin,2)<t);
        %get the number of traces excluded
        thres_res(t) = sum(snr_vec);
    end

    %based on the results above, I'll set the threshold to be at least 10
    %stimuli with significant signal, which excludes around 11% of the traces
    snr_vec = logical(sum(snr_bin,2)>stimulus_constant);

    %adapt the conc_trace and fish_ori matrices 
    conc_trace = conc_trace(snr_vec,:);
    fish_ori = fish_ori(snr_vec,:);
    % also the single reps
    reps_trace = reps_trace(snr_vec,:,:,:);    
    %% Save analysis output
    
    % assemble the path for the structure
    % parse the Stage 1 path
    name_parts = strsplit(tar_path_all{exps},'\');
    fish_name = name_parts{8};
    condition = name_parts{9};
    file_name = strjoin({condition,fish_name,'meta'},'_');
    
    
    %define the save path
%     save_path = 'E:\Behavioral data\Matlab\AF_proc\Clemens_suite\pipeline_test\Stage2\';

%     save_var = 1;
%     if save_var == 1
% 
%     %     %get the root of the save name
%     %     [ori_name,~] = uiputfile(strcat(save_path,'*.*'));
%         %aseemble the path
%         full_name = strsplit(tar_path_all{exps},'\');
% 
%         %save the clustering output
%         save_clu = strcat(full_name{8},'_',full_name{9},'_wholefish.mat');
%         save(fullfile(save_path,save_clu),'conc_trace','stim_num2','snr_vec','-v7.3')
%     end
    %% Assemble the structure with the metadata
    %% Save the structure
end