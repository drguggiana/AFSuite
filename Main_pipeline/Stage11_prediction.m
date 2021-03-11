%% Predict stimulus
clearvars
close all
Paths
%% Define some constants

% define the number of reps per classifier
num_classreps = 5;
% define the flags
subsample_flag = 1;
label_flag = 0;
roi_flag = 1;

% assemble the parameter structure
params = struct([]);
params(1).num_classreps = num_classreps;
params(1).subsample_flag = subsample_flag;
params(1).label_flag = label_flag;
params(1).roi_flag = roi_flag;

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
confusion_cell = cell(num_fish,3,num_regions);
randomtime_cell = cell(num_fish,3,num_regions,num_classreps);
randomroi_cell = cell(num_fish,3,num_regions,num_classreps);
% for all the fish
for fish = 1:num_fish
    
    % get the traces corresponding to this fish
    conc_trace = conc_trace_all(fish_ori_all==fish,:,:,:);
    
    % load the celltype info
    cell_info = main_cell{fish}.cell_type;
    % load the region info
    region_info = main_cell{fish}.AF_info(:,1);
    
    % for all celltypes
    for celltype = 1:3
        
        % get the region info for this celltype
        cell_region = region_info(cell_info==celltype-1,1);
        % turn NaN into 0
        cell_region(isnan(cell_region)) = 0;
        % calculate the minimum number of ROIs across regions
        af_numbers = cellfun(@str2double,{af_labels(:).number});
        counts = count_occurrences(cell_region,[0 af_numbers]);
        min_perregion = min(counts(counts>0));
        
        % for all the regions
        for region = 1:num_regions
            
            disp(strcat('Fish:',num2str(fish),',type:',num2str(celltype),',region:',num2str(region)))
            % get the traces for the region
            current_traces = conc_trace(cell_region==str2double(af_labels(region).number),:,:,:);
            % if there are no traces in this region, skip and record a nan
            if size(current_traces,1) == 0
                confusion_cell{fish,celltype,region} = NaN;
                continue
            end
            % select only the stimulus period
            current_traces = current_traces(:,21:70,:,:);
            
            % allocate memory for the rep results
            rep_results = cell(num_classreps,2);
            
            % for all the classreps
            for classr = 1:num_classreps
                
                % subsample based on the stimulus
                [subsampled_traces,random_time] = subsample_data(current_traces,subsample_flag);
                % subsample based on the rois
                [subsampled_traces,random_rois] = subsample_roi(subsampled_traces,roi_flag,min_perregion);
                % generate the label vector
                label_vector = label_maker(subsampled_traces,label_flag);
                % reshape the data and labels to 2D
                subsampled_traces = reshape(subsampled_traces,size(subsampled_traces,1),[]);
                label_vector = label_vector(:);
                % train the classifier
                model = fitcecoc(subsampled_traces',label_vector,'Prior','Uniform','Learners','svm',...
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
            % save
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
save_name = strjoin({'Classification','classreps',num2str(num_classreps),...
    'tsubsample',num2str(subsample_flag),'rsubsample',num2str(roi_flag),'label',num2str(label_flag)},'_');
save_path = fullfile(classification_path,strcat(save_name,'.mat'));

% save the results
save(save_path,'confusion_cell','params')