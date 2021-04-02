%% Predict stimulus on a whole fish basis
clearvars
close all
Paths
%% Define some constants

% define the number of reps per classifier
num_classreps = 5;
% define the flags
subsample_flag = 0;
label_flag = 0;
roi_flag = 1;
roi_percentage = 0.95;

% assemble the parameter structure
params = struct([]);
params(1).num_classreps = num_classreps;
params(1).subsample_flag = subsample_flag;
params(1).label_flag = label_flag;
params(1).roi_flag = roi_flag;
params(1).roi_percentage = roi_percentage;

% load the labels
labels = load(constants_path,'labels');
labels = labels.labels;

af_labels = labels.af;

% get the total number of regions
num_regions = length(af_labels);
%% Load the data

% assemble the overall path
main_path = fullfile(analysis_path,'Meta_files');

% load the paths
[group_folders,num_subfolders] = group_loader('pre',main_path);

% load the traces
[main_cell,conc_trace_all,fish_ori_all,num_fish,name_cell] = main_loader(group_folders(1),2);
%% Use the neural activity to determine which stimulus was being presented
tic
% allocate memory to save the results
confusion_cell = cell(num_fish,3,num_regions+1);
randomtime_cell = cell(num_fish,3,num_regions+1,num_classreps);
randomroi_cell = cell(num_fish,3,num_regions+1,num_classreps);

% load the celltype and region info
cell_info = cell(num_fish,1);
region_info = cell(num_fish,1);
for fish = 1:num_fish
    cell_info{fish} = main_cell{fish}.cell_type;
    % load the region info
    region_info{fish} = main_cell{fish}.AF_info(:,1);
end

cell_info = vertcat(cell_info{:});
region_info = vertcat(region_info{:});
% turn nans into 0
region_info(isnan(region_info)) = 0;
% for all the regions
for region = 1:num_regions+1
    
    % for all the fish
    for celltype = 1:3
        
        % remove the target AF from the data
        if region <= num_regions
            % get the selection vector
            selection_vector = (region_info~=str2double(af_labels(region).number))&region_info>0;
            af_fish_ori = fish_ori_all(selection_vector);
            af_cell_info = cell_info(selection_vector);
            af_conc_trace = conc_trace_all(selection_vector,:,:,:);
        else
            af_fish_ori = fish_ori_all;
            af_cell_info = cell_info;
            af_conc_trace = conc_trace_all;
        end
        % get the number of traces to use per fish/celltype
        fish_numbers = af_fish_ori(af_cell_info==celltype-1);
        counts = count_occurrences(fish_numbers,1:num_fish);
        min_perfish = min(counts);
        
        % for all celltypes
        for fish = 1:num_fish
            disp(strcat('Fish:',num2str(fish),'type:',num2str(celltype),'region:',num2str(region)))
            
            % get the traces 
            conc_trace = af_conc_trace(af_fish_ori==fish&af_cell_info==celltype-1,:,:,:);
            % select the target time span
            conc_trace = conc_trace(:,21:68,:,:);
            % allocate memory for the reps
            rep_results = cell(num_classreps,2);
            % for all the classreps
            for classr = 1:num_classreps
                
                % subsample based on the stimulus
                [subsampled_traces,random_time] = subsample_data(conc_trace,subsample_flag);
                % subsample based on the rois
                [subsampled_traces,random_rois] = subsample_roi(subsampled_traces,roi_flag,min_perfish,roi_percentage);
                % generate the label vector
                label_vector = label_maker(subsampled_traces,label_flag);
                % reshape the data and labels to 2D
                subsampled_traces = reshape(subsampled_traces,size(subsampled_traces,1),[]);
                label_vector = label_vector(:);
                % standardize the data
                subsampled_traces = normalize(subsampled_traces');
                % set up the template
%                 t = templateSVM('KernelFunction','linear','Standardize',1);
                t = templateLinear('Lambda',1e-4,'Regularization','lasso');
                % train the classifier
                model = fitcecoc(subsampled_traces,label_vector,'Prior','Uniform','Learners',t,...
                    'KFold',5,'Coding','onevsall');
                
                % get conf matrix
                label_predicted = kfoldPredict(model);
                conf_matrix = confusionmat(label_vector,label_predicted);
                % save for averaging
                rep_results{classr,1} = label_predicted;
                rep_results{classr,2} = conf_matrix;
                
                % save the randomization vectors
                randomtime_cell{fish,celltype,region,classr} = random_time;
                randomroi_cell{fish,celltype,region,classr} = random_rois;
            end
            % store
            confusion_cell{fish,celltype,region} = rep_results;

        end
    end
end

% store the randomization in the structure
params(1).random_rois = randomroi_cell;
params(1).random_time = randomtime_cell;

toc
%% Save the results

% define the save path
save_name = strjoin({'ClassificationPerFish','classreps',num2str(num_classreps),...
    'tsubsample',num2str(subsample_flag),'rsubsample',num2str(roi_flag),'label',num2str(label_flag)},'_');
save_path = fullfile(classification_path,strcat(save_name,'.mat'));

% save the results
save(save_path,'confusion_cell','params')