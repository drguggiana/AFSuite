function [subsampled_rois,random_vector] = subsample_roi(data_in,flag,varargin)
% subsample ROIs

% select the subsampling method
switch flag
    case 0 % don't subsample
        subsampled_rois = data_in;
        random_vector = [];
    case 1 % subsample to 95% of the minimum number
        % get the number from varargin
        min_number = ceil(varargin{1}*0.95);
        
        % get the randomization vector
        random_vector = randperm(size(data_in,1),min_number);
        % subsample the data
        subsampled_rois = data_in(random_vector,:,:,:);
end