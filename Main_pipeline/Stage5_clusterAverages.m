%% clean up
clearvars
close all
Paths
%% Load the traces

% define the figure path
figure_path = fullfile(fig_path,'clusterAverages');

% define the save path
save_path = fullfile(analysis_path,'Meta_files','all_meta.mat');

% load the labels
labels = load(constants_path,'labels');
labels = labels.labels;

af_labels = labels.af;
celltype_labels = labels.celltype;

% define the groups to plot
target_groups = {'pre','postcontrol'};

% assemble the overall path
main_path = fullfile(analysis_path,'Meta_files');

% load the paths
[group_folders,num_subfolders] = group_loader(target_groups,main_path);

% allocate memory to store the compiled results
% group_struct = cell(num_subfolders,4);
average_str = struct([]);
% for all the subfolders
for group = 1:num_subfolders
    % load the traces
    [main_cell,~,~,num_fish] = main_loader(group_folders(group),0);
    
    % allocate a temp structure
    temp_str = ([]);
    % allocate memory for the averages
    cluster_averages = cell(num_fish,3);
    % allocate memory for the region averages
    region_cell = cell(num_fish,3);
    
    % for all the fish
    for fish = 1:num_fish
        fprintf(strjoin({'Current fish:',num2str(fish),'out of',num2str(num_fish),'\r\n'},'_'))
        % get the indexes
        idx_cell = main_cell{fish}.idx_cell;
        % get the number of clusters
        clu_num = main_cell{fish}.clu_num;
        % get the cell type
        cell_type_info = main_cell{fish}.cell_type;
        % get the dimensions of the trace file
        seed_num = size(main_cell{fish}.fish_ori,1);
        stim_num = main_cell{fish}.stim_num;
        time_num = main_cell{fish}.time_num;
        % load the traces
        fID = fopen(main_cell{fish}.path_trial_ave);
        conc_trace = reshape(fread(fID,'double'),seed_num,time_num,stim_num);
        fclose(fID);
        
        % load the anatomy info
        anatomy_info = main_cell{fish}.AF_info;
        % get a list of the regions
        region_list = unique(anatomy_info(:,1));
        % remove the NaN entries
        region_list = region_list(~isnan(region_list));
        % get the number of regions
        region_number = length(region_list);
        
        % get the cluster averages
        % for all the cell types
        for celltype = 1:3
            % get the cluster number
            temp_clu_num = clu_num(celltype);
            % allocate memory for the average matrix
            cluster_average_matrix = zeros(temp_clu_num,time_num,stim_num);
            region_average_matrix = struct([]);

            % get the traces for this cell type
            type_traces = conc_trace(cell_type_info==(celltype-1),:,:);
            type_regions = anatomy_info(cell_type_info==(celltype-1),1);
            % for all the regions
            for region = 1:region_number
                region_temp_matrix = zeros(temp_clu_num,time_num,stim_num);
                % for all the clusters
                for clu = 1:temp_clu_num
                    % if it's the first region, fill in the overall average
                    % also
                    if region == 1
                        % get the traces for this cluster average and store
                        cluster_average_matrix(clu,:,:) = mean(type_traces(idx_cell{celltype}==clu,:,:),1);  
                    end
                    % fill in the region matrix
                    region_temp_matrix(clu,:,:) = nanmean(type_traces(idx_cell{celltype}==clu&...
                            type_regions==region_list(region),:,:),1);
  
    %                 % for all the regions
    %                 for region = 1:region_number
    % %                     region_average_matrix(clu,:,:,region) = mean(type_traces(idx_cell{celltype}==clu&...
    % %                         type_regions==region_list(region),:,:),1);
    %                     region_average_matrix(                    = mean(type_traces(idx_cell{celltype}==clu&...
    %                         type_regions==region_list(region),:,:),1);
    %                 end
                end
                % find the region name
                region_name = strcmp({af_labels.number},num2str(region_list(region)));
                % append to the structure
                region_average_matrix(1).(af_labels(region_name).name) = region_temp_matrix;
                
            end
            % store in the main cell
            cluster_averages{fish,celltype} = cluster_average_matrix;
            region_cell{fish,celltype} = region_average_matrix;
        end

    end
    %% Average across animals
    
    % average across fish
    fish_average = cell(3,1);
%     region_average = cell(3,1);
    region_average = struct([]);
    % for all the cell types
    for celltype = 1:3
        % average across all fish
        fish_average{celltype} = mean(cat(4,cluster_averages{:,celltype}),4);
        % for all the regions
        for region = 1:length(af_labels)
            % get the current region name
            current_region = af_labels(region).name;
            % allocate memory to store the temp region data
            temp_region = cell(num_fish,1);
            % for all the fish
            for fish = 1:num_fish
                % check if the field is in this fish
                if sum(contains(fields(region_cell{fish,celltype}),current_region))==0
                    continue
                end
                % get the region for this fish
                temp_region{fish} = region_cell{fish,celltype}.(current_region);
            end
            % average across fish and store
            region_average(celltype).(current_region) = nanmean(cat(4,temp_region{:}),4);
        end
    end
    %% Plot the averages
    
    close all
    for celltype = 1:3
        figure
        % get the traces
        traces = reshape(fish_average{celltype},[],time_num*stim_num);
        % sort and plot
        imagesc(normr_1(sort_traces(traces),0))
        title(strjoin({target_groups{group},celltype_labels{celltype}},'_'),'Interpreter','None')
        set(gca,'TickLength',[0 0])
        % define the path and save
        file_path = strjoin({'clusterAverages',target_groups{group},celltype_labels{celltype},'.png'},'_');
        print(fullfile(figure_path,file_path), '-dpng','-r600')
        
    end
    autoArrangeFigures
    %% Plot the cluster distribution per area
    close all
    
    % allocate memory for the anatomy breakdown
    anatomy_cell = cell(num_fish,3);
    % get the total number of regions
    num_regions = length(af_labels);
    % get the number of clusters (should be the same for the whole group)
    clu_num = main_cell{1}.clu_num;
    
    % for all the fish
    for fish = 1:num_fish
        % load the anatomy info
        anatomy_info = main_cell{fish}.AF_info;
        % load the indexes
        index_info = main_cell{fish}.idx_cell;
        % load the cell type info
        celltype_info = main_cell{fish}.cell_type;
        % for all the cell types
        for celltype = 1:3
            % get the anatomy for this cell type
            curr_anatomy = anatomy_info(celltype_info==celltype-1,1);
            % allocate memory to store the results
            anatomy_matrix = zeros(num_regions,clu_num(celltype));
            % for all the clusters
            for clu = 1:clu_num(celltype)
                % for all the regions
                for region = 1:num_regions
                    % add the number of traces
                    anatomy_matrix(region,clu) = anatomy_matrix(region,clu) +...
                        sum(curr_anatomy==str2double(af_labels(region).number)&index_info{celltype}==clu);
                    
                end
            end
            % store the resulting matrix
            anatomy_cell{fish,celltype} = anatomy_matrix;
        end
    end
    %% Average the anatomy and plot
    close all
    % for all the cell types
    for celltype = 1:3
        average_anatomy = mean(cat(3,anatomy_cell{:,celltype}),3);
        figure
        imagesc(normr_1(average_anatomy,1))
        set(gca,'YTick',1:num_regions,'YTickLabels',{af_labels.name})
        title(strjoin({target_groups{group},celltype_labels{celltype}},'_'),'Interpreter','None')
        set(gca,'TickLength',[0 0])
        % define the path and save
        file_path = strjoin({'clusterPerRegion',target_groups{group},celltype_labels{celltype},'.png'},'_');
        print(fullfile(figure_path,file_path), '-dpng','-r600')
    end
    
    autoArrangeFigures
    
    % store the averages
    temp_str(1).cluster_ave_all = fish_average;
    temp_str(1).region_counts = anatomy_cell;
    temp_str(1).meta = vertcat(main_cell{:});
    temp_str(1).region_ave = region_average;
    temp_str(1).region_ave_perfish = region_cell;
    
    average_str(1).(target_groups{group}) = temp_str;
end
%% Save the structure

save(save_path,'average_str','-v7.3')