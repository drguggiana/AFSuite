function [subsampled_rois,random_vector] = subsample_roi(data_in,flag,varargin)
% subsample ROIs

% select the subsampling method
switch flag
    case 0 % don't subsample
        subsampled_rois = data_in;
        random_vector = [];
    case 1 % subsample to 95% of the minimum number
        % get the number from varargin
        min_number = ceil(varargin{1}*varargin{2});
        
        % get the randomization vector
        random_vector = randperm(size(data_in,1),min_number);
        % subsample the data
        subsampled_rois = data_in(random_vector,:,:,:);
    case 2 % subsample all AFs present to the same number of ROIs
        % get the number from varargin
%         max_number = ceil(varargin{3}*varargin{2});
        
        % get the current region
        current_region = varargin{1};
        total_rois = varargin{3};
        
%         % get the number of regions
%         num_regions = size(total_rois,1);
        % get the number to remove
        target_remove = max(total_rois)-total_rois(current_region);

        % get the number of traces remaining
        traces_remaining = round((size(data_in,1) - target_remove).*varargin{2});
%         % get a list of the regions
%         region_list = unique(region_info);
%         % get the number of regions
%         num_regions = length(region_list);
        % create the random vector
        random_vector = randperm(size(data_in,1),traces_remaining);
        % take the subselection
        subsampled_rois = data_in(random_vector,:,:,:);
%         % allocate memory to store the output traces
%         subsampled_rois = zeros(traces_remaining,num_regions,size(data_in,2));
%         % allocate memory for the random vectors
%         random_vector = zeros(traces_remaining,num_regions);
%         % for all the regions
%         for region = 1:num_regions
%             % get the traces in this region
%             region_traces = data_in(region_info==region_list(region),:);
%             % get the random vector
%             random_temp = randperm(size(region_traces,1),min_number);
%             % store it
%             random_vector(:,region) = random_temp;
%             % get the random sample and store
%             subsampled_rois(:,region,:) = region_traces(random_temp,:);
%         end
        % remove the region dimension
%         subsampled_rois = permute(reshape(permute(subsampled_rois,[3 4 5 1 2]),[],max_number*num_regions),[4 1 2 3]);
%         random_vector = reshape(random_vector,[],1);
end