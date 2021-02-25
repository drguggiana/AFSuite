% Model the interaction between regions

clearvars
close all
Paths
%% Load the meta structure

average_str = load(metapost_path);
average_str = average_str.average_str;
%% Define some constants

% target_groups = fields(average_str);
target_groups = {'pre'};
num_subfolders = length(target_groups);

% load the labels
labels = load(constants_path,'labels');
labels = labels.labels;

af_labels = labels.af;
celltype_labels = labels.celltype;

% define the figure path
figure_path = fullfile(fig_path,'clusterAverages');

% get the total number of regions
num_regions = length(af_labels);
%% Run the modelling 
% TODO: use the clusters per region
close all
% % define the regions to use
% tectum_regions = {'L-TcN','R-TcN','L-TcP','R-TcP','L-Cb','R-Cb','L-Hb','R-Hb','L-Pt','R-Pt'};
% tectum_numbers = 1:10;
% af_regions = {'AF4','AF5','AF8','AF9','AF10'};
% af_numbers = [1 2 5 6 7];

% num_datasets = size(data,2);
num_datasets = 3;

% allocate memory to store the model data
model_cell = cell(num_datasets,1);
% allocate memory to store the region information
region_cell = cell(num_datasets,2);
% for all the data sets
for datas = 1:3 %(ori, dir and ns)
%     % get the region info
%     region_info = data(datas).anatomy_info(:,1);
%     % define the regions to be considered (depending on the stimulus protocol)
%     if contains(data(datas).name,{'syn','Syn'})
%         region_list = af_regions;
%         region_numbers = af_numbers;
%     else
%         region_list = tectum_regions;
%         region_numbers = tectum_numbers;
% %     end
%     % get the region numbers 
%     num_data = length(region_numbers);

    %get all the pairwise combinations of the regions
    region_comb = [nchoosek(1:num_regions,2);fliplr(nchoosek(1:num_regions,2))];

    %get the number of combs
    num_comb = size(region_comb,1);
    
    % allocate memory for the model data within this dataset
    current_models = cell(num_comb,1);
% 
%     % define the period of interest (0 pre, 1 stim, 2 post, 3 pre-post)
%     period = 1;
%     % get the target period labeled with ones
%     rest_all = period_of_interest(period,data(datas).stim_num,1);
    %for all the combinations
    for combs = 1:num_comb
        %concatenate the clusters from both animals involved
%         tar1 = data(datas).region_clusters(region_comb(combs,1)).clu_ave;
%         tar2 = data(datas).region_clusters(region_comb(combs,2)).clu_ave;
        name1 = af_labels(region_comb(combs,1)).name;
        tar1 = average_str.pre.region_ave(datas).(name1);
        tar1 = reshape(tar1,size(tar1,1),[]);
        tar1 = tar1(sum(isnan(tar1),2)~=size(tar1,2),:);
        
        name2 = af_labels(region_comb(combs,2)).name;
        tar2 = average_str.pre.region_ave(datas).(name2);
        tar2 = reshape(tar2,size(tar2,1),[]);
        tar2 = tar2(sum(isnan(tar2),2)~=size(tar2,2),:);

        
        % if either of them is empty, put a nan in the cell and skip
        if isempty(tar1) || isempty(tar2)
            current_models{combs} = NaN;
            continue
        end
        
        %for all the averages in 1, calculate models from the raw traces in 2
        %allocate memory to store the model results
        model_para = cell(size(tar1,1),1);
        %for all the averages
        for clu = 1:size(tar1,1)
            fprintf(strcat('Current comb:',num2str(combs),'Current clu: ',num2str(clu),'\r\n'))
            
            model_para{clu} = fitrlinear(tar2',tar1(clu,:)','CrossVal','on');

        end
        % fill up the dataset cell
        current_models{combs} = model_para;

    end
    % fill up the overall cell
    model_cell{datas} = current_models;
end
%% Plot the results

close all

% for all the data sets
for datas = 1:3
    figure
%     % define the regions to be considered (depending on the stimulus protocol)
%     if contains(data(datas).name,{'syn','Syn'})
%         region_list = af_regions;
%         region_numbers = af_numbers;
%     else
%         region_list = tectum_regions;
%         region_numbers = tectum_numbers;
%     end
%     % get the region numbers 
%     num_data = length(region_numbers);
    %get all the pairwise combinations of the regions
    region_comb = [nchoosek(1:num_regions,2);fliplr(nchoosek(1:num_regions,2))];

    %get the number of combs
    num_comb = size(region_comb,1);
    % allocate memory for the combination matrix
    combination_matrix = zeros(num_regions);
    % for all the combinations
    for combs = 1:num_comb
        % get the coordinates of the combination
        target = region_comb(combs,1);
        source = region_comb(combs,2);
        % get the number of clusters in the x_coord
        name1 = af_labels(region_comb(combs,1)).name;
        tar1 = average_str.pre.region_ave(datas).(name1);
        tar1 = reshape(tar1,size(tar1,1),[]);
        tar1 = tar1(sum(isnan(tar1),2)~=size(tar1,2),:);
        clu_num = size(tar1,1);
        % if the model of interest is a nan, put a nan and skip the
        % iteration
        if ~iscell(model_cell{datas}{combs}) && isnan(model_cell{datas}{combs})
            combination_matrix(target,source) = nan;
            continue
        end
        % calculate the average of the cluster losses for this combination
        % for all the clusters
        for clu = 1:clu_num
            combination_matrix(source,target) = ...
                combination_matrix(source,target) + ...
                kfoldLoss(model_cell{datas}{combs}{clu})/clu_num;
        end
    end
    % plot the matrix
    imagesc(log(combination_matrix))
%     if contains(data(datas).name,{'syn','Syn'})
%         set(gca,'XTick',1:num_data,'XTickLabel',region_list,'XTickLabelRotation',45)
%     else
%         set(gca,'XTick',1:2:num_data,'XTickLabel',{'TcN','TcP','Pt','Hb','Cb'},'XTickLabelRotation',45)
%     end
%     set(gca,'YTick',1:num_data,'YTickLabel',region_list)
%     set(gca,'CLim',[0.65 1])
    title(celltype_labels{datas})
    set(gca,'XTick',1:num_regions,'XTickLabels',{af_labels.name},...
        'YTick',1:num_regions,'YTickLabels',{af_labels.name},...
        'TickLabelInterpreter','None','TickLength',[0 0],...
        'XTickLabelRotation',45)
    axis square
    c = colorbar;
%     % create the settings
%     fig_set = struct([]);
%     
%     fig_set(1).fig_path = fig_path;
%     fig_set(1).fig_name = strjoin({'modelMatrix',data(datas).name,'.eps'},'_');
%     fig_set(1).fig_size = 5;
%     fig_set(1).colorbar = 1;
%     fig_set(1).colorbar_label = 'Goodness of Fit';
%     fig_set(1).box = 'on';
%     fig_set(1).LineWidth = 0.05;
%     
%     h = style_figure(gcf,fig_set);
    
end
autoArrangeFigures