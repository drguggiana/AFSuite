function [subsampled_data,random_vector] = subsample_data(data_in,flag)
% function to subsample the data according to stimulus selection

% select the subsampling to use
switch flag
    case 0 % consider all stimuli, so no subsampling
        subsampled_data = data_in;
        random_vector = [];
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
end