function label_vector = label_maker(data_in,flag)
% Function to make labels for the training of neural classifiers


% select the corresponding case
switch flag
    case 0 % label all stimuli separately
        
        % get the info from the data_in
        roi_num = size(data_in,1);
        time_num = size(data_in,2);
        stim_num = size(data_in,3);
        rep_num = size(data_in,4);
        
        % allocate the label matrix
        label_vector = zeros(time_num,stim_num,rep_num);
        % for all the stimuli
        for stim = 1:stim_num
            % generate the label
            labels = zeros(time_num,rep_num)+stim;
            % assign to the output
            label_vector(:,stim,:) = labels;
        end
        
%     case 1 % label the categories of stimuli
        
end


