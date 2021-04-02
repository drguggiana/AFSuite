function [subsampled_data,random_vector] = subsample_data(data_in,flag)
% function to subsample the data according to stimulus selection

% select the subsampling to use
switch flag
    case 0 % consider all stimuli, so only 95% subsampling
        roi_num = size(data_in,1);
        time_num = size(data_in,2);
        stim_num = size(data_in,3);
        rep_num = size(data_in,4);
        % get the target number of time bins
        target_time = ceil(time_num.*0.95);
        
        % allocate for the downsampled matrices
        subsampled_data = zeros(roi_num,target_time,stim_num,rep_num);
        
        % allocate memory to save the randomization vector
        random_vector = cell(stim_num,rep_num);
        
        % for all the categories
        for stim = 1:stim_num

            % for all the reps
            for reps = 1:rep_num
                % get the data
                cat_data = reshape(data_in(:,:,stim,reps),roi_num,[]);
                % get the random vector
                random_vector{stim,reps} = randperm(size(cat_data,2),target_time);
                % subsample
                subsampled_data(:,:,stim,reps) = cat_data(:,random_vector{stim,reps});
            end
        end
    case 1 % subsample to match the stim categories (95% of min)
        % get the info from the data_in
        roi_num = size(data_in,1);
        time_num = size(data_in,2);
        rep_num = size(data_in,4);
        % define the indexes of the categories
        cat_index = {[1 10],2:9,11:18,19:26};
        
        % get the number of categories
        cat_number = length(cat_index);
        % assuming all stimuli last the same, get the number of time points
        % from the category with the least indexes
        stim_percat = cellfun(@length,cat_index);
        
        % get the target number of time bins
        target_time = ceil(min(stim_percat).*time_num.*0.95);
        
        % allocate for the downsampled matrices
        subsampled_data = zeros(roi_num,target_time,cat_number,rep_num);
        
        % allocate memory to save the randomization vector
        random_vector = cell(cat_number,rep_num);
        
        % for all the categories
        for cats = 1:cat_number

            % for all the reps
            for reps = 1:rep_num
                % get the data
                cat_data = reshape(data_in(:,:,cat_index{cats},reps),roi_num,[]);
                % get the random vector
                random_vector{cats,reps} = randperm(size(cat_data,2),target_time);
                % subsample
                subsampled_data(:,:,cats,reps) = cat_data(:,random_vector{cats,reps});
            end
        end      
    case 2 % split the categories into 4 time bins
        % define the time bin number
        bin_number = 4;
        % get the info from the data_in
        roi_num = size(data_in,1);
        time_num = size(data_in,2);
        rep_num = size(data_in,4);
        % define the indexes of the categories
        cat_index = {[1 10],2:9,11:18,19:26};
        % assuming all stimuli last the same, get the number of time points
        % from the category with the least indexes
        stim_percat = cellfun(@length,cat_index);
        
        % get the target number of time bins
        target_time = ceil(min(stim_percat).*time_num/bin_number.*0.95);
        % get the indexes for the time
        time_index = [1:time_num/bin_number:time_num;...
            time_num/bin_number:time_num/bin_number:time_num];
        
        % get the number of categories
        cat_number = length(cat_index)*size(time_index,2);
        
        % allocate for the downsampled matrices
        subsampled_data = zeros(roi_num,target_time,cat_number,rep_num);
        
        % allocate memory to save the randomization vector
        random_vector = cell(length(cat_index),size(time_index,2),rep_num);
        
        % initialize a counter for the "stimuli" (times and categories)
        count = 1;
        

        % for all the categories
        for cats = 1:length(cat_index)
            % for all the time bins
            for times = 1:size(time_index,2)
            
                % for all the reps
                for reps = 1:rep_num
                    % get the data
                    cat_data = reshape(data_in(:,time_index(1,times):time_index(2,times),...
                        cat_index{cats},reps),roi_num,[]);
                    % get the random vector
                    random_vector{cats,times,reps} = randperm(size(cat_data,2),target_time);
                    % subsample
                    subsampled_data(:,:,count,reps) = cat_data(:,random_vector{cats,times,reps});
                end
                % update the counter
                count = count + 1;
            end
        end
end