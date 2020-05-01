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
tar_path_all_pre = uipickfiles('FilterSpec',fullfile(analysis_path,'Stage1'));

% define the save path
save_path = fullfile(analysis_path,'Meta_files');

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
    
    fprintf(strjoin({'Current experiment:',num2str(exps),'out of',num2str(num_exp),'\r\n'},'_'))
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
    stim_num = load(name_cell{1},'stim_num');
    stim_num = stim_num.stim_num;
    %% Get the snr threshold

    %load the thresholding parameters from the constants file    
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
    thres_vec = zeros(stim_num,1);

    %for all of the stimuli
    for stim = 1:stim_num
        %calculate the 25th percentile and store as the threshold
        thres_vec(stim) = prctile(snr_mat(:,stim),percentile_constant);
    end
    %% Load the traces 
    
    %define a list of the stimuli to keep
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
        
        % concatenate the coordinate information
        % for all the seeds
        for seeds = 1:length(z_seed)
            seed_concat(seeds).z = z_seed(seeds);
            seed_concat(seeds).experiment = fish;
        end
        % store in the cell
        seed_data{fish} = seed_concat;
        
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
                stim_num = length(keep_stim);
                %also modify the color matrix
                col_out = col_out(keep_stim,:,:);
            end
        end
    end

    %concatenate the responses
    conc_trace = vertcat(conc_trace{:});
    reps_trace = vertcat(reps_trace{:});
    seed_data = vertcat(seed_data{:});
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

    thres_res = zeros(stim_num,1);
    for t = 1:stim_num
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
    % and the coordinates
    seed_data = seed_data(snr_vec);
    %% Save analysis output
    
    % assemble the path for the traces
    % parse the Stage 1 path
    name_parts = strsplit(tar_path_all{exps},'\');
    fish_name = name_parts{8};
    condition = name_parts{9};
    
    file_name_trial_ave = strjoin({condition,fish_name,'trial_ave'},'_');
    path_trial_ave = fullfile(traces_path,file_name_trial_ave);
    file_name_single_reps = strjoin({condition,fish_name,'single_reps'},'_');
    path_single_reps = fullfile(traces_path,file_name_single_reps);
    
    % save the files
    file_ID = fopen(path_single_reps,'w');
    fwrite(file_ID,conc_trace,'double');
    fclose(file_ID);
    
    file_ID = fopen(path_trial_ave,'w');
    fwrite(file_ID,reps_trace,'double');
    fclose(file_ID);
    
    %% Load the average stack
    ave_stack = load(name_cell{fish,1},'ave_stack');
    ave_stack = ave_stack.ave_stack;
    %% Assemble the structure with the metadata
    main_str = struct([]);
    
    % fish name
    main_str(1).fish_name = fish_name;
    % condition
    main_str(1).condition = condition;
    % path to the traces file
    main_str(1).path_trial_ave = path_trial_ave;
    % path to the single reps
    main_str(1).path_single_reps = path_single_reps;
    % fish ori
    main_str(1).fish_ori = fish_ori;
    % coordinate info
    main_str(1).seed_data = seed_data;
    % snr
    main_str(1).snr_vec = snr_vec;
    % average stack
    main_str(1).ave_stack = ave_stack;
    % time num
    main_str(1).time_num = time_num;
    % stim num
    main_str(1).stim_num = stim_num;
    % number of reps
    main_str(1).rep_num = size(reps_trace,4);
    % list of files included in the experiment (matching the seed_data)
    main_str(1).files = name_cell;
    
    %% Save the structure
    file_name_meta = strjoin({condition,fish_name,'meta'},'_');
    
    save(fullfile(save_path,condition,file_name_meta),'main_str')
    
end