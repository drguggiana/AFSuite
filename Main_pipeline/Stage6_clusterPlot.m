%% Ablation comparisons

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
%% Plot the cluster averages

% define the target regions
target_regions = {'AF5','AF8','AF9d','AF9v'};
num_targets = length(target_regions);

close all
% for all celltypes
for celltype = 1:3
    figure
    counter = 1;
    % get the data
    main1 = average_str.(target_groups{1}).region_ave(celltype);
    % for all the regions
    for region = 1:length(af_labels)
        
        
        % get the current region name
        current_region = af_labels(region).name;
        if sum(contains(target_regions,current_region))==0
            continue
        end
        % extract the region of interest
        average1 = main1.(current_region);

        
        % reshape them to 2D
        average1 = sort_traces(reshape(average1,size(average1,1),[]));
        

        %             % nan the non-significant values
        %             correlation(pval<0.05) = NaN;
        
        subplot(round(sqrt(num_targets)),ceil(sqrt(num_targets)),counter)
        counter = counter + 1;
        imagesc(normr_1(average1,0))
        set(gca,'CLim',[0 1],'TickLength',[0 0],'YTick',[],'XTick',[])
        xlabel(current_region)
        axis tight

%         % plot the correlation
%         subplot(4,4,region)
%         imagesc(correlation)
%         axis equal
%         axis tight
%         set(gca,'CLim',[0 1],'TickLength',[0 0],'YTick',[],'XTick',[])
%         title(strjoin({target_groups{combo_vector(combo,1)},...
%             target_groups{combo_vector(combo,2)},current_region},' '),...
%             'Interpreter','None')
    end
    sgtitle(celltype_labels{celltype})
    
end
%% Quantify correspondence between clusters
close all
% define the type of extraction
stimTypeNum = 4;
% get the number of fish in the pre group
pre_number = size(average_str.pre.meta,1);


    
% for all the celltypes
for celltype = 1:3
    per_fish = zeros(pre_number,num_regions,8,2);
    % for all the fish in the pre group
    for fish = 1:pre_number
        
        fprintf(strjoin({'Current fish:',num2str(fish),'out of',num2str(pre_number),'\r\n'},' '))
        % get the fish name
        fish_name = average_str.pre.meta(fish).fish_name;
%         corr_fig = figure;
%         histo_fig = figure;
%         % set up a counter for the matches
%         match_counter = struct([]);
        
        
        % get the fish averages per region
        pre_averages = average_str.pre.region_ave_perfish{fish,celltype};
        
        % get the matching averages
        
        % get the counts too
        pre_counts = average_str.pre.region_counts;
        
        % compare the averages
        % for all the regions
        for region = 1:length(af_labels)
            % get the current region name
            current_region = af_labels(region).name;
            
            % check if the region is in the fish
            if contains(current_region,fields(pre_averages)) == 0
                continue
            end
            % extract the region of interest
            average1 = pre_averages.(current_region);
            
            average1(isnan(average1)) = 0;
%             if isempty(average1) || isempty(average2)
%                 continue
%             end
            
%             % reshape them to 2D
%             average1 = sort_traces(reshape(average1,size(average1,1),[]));
%             average2 = sort_traces(reshape(average2,size(average2,1),[]));
% 
%             % get the correlation matrix between them and plot
%             [correlation,pval] = corr(average1',average2','rows','pairwise');

            % run the svd on the averages
            [~,pre_svd] = fish_svd(average1,stimTypeNum);

            % separate the cells based on pref dir/ori and on/off
            [~,pref_dir] = max(pre_svd(:,51:58),[],2);
            on_off = (pre_svd(:,225)<pre_svd(:,226))+1;
            % collapse into a matrix with only on/off and dir
            temp_pre = zeros(8,2);
            
            % for all the clusters
            for clu = 1:size(pref_dir,1)
                temp_pre(pref_dir(clu),on_off(clu)) = ...
                    temp_pre(pref_dir(clu),on_off(clu)) + pre_counts{fish,celltype}(region,clu);
            end
              
            
            % store the overall result per animal
            per_fish(fish,region,:,:) = temp_pre;
            
        end
%         % save the correlation figure
%         figure(corr_fig)

        
        
%         % define the path and save
%         set(gcf,'Color','w')
%         file_path = fullfile(figure_path,'PerFish',strjoin({'correlationPerRegionPerFish','pre',...
%             match_group,celltype_labels{celltype},fish_name,'.png'},'_'));
%         export_fig(file_path,'-r600')
        

    end
    % average across fish and plot
    figure
    sgtitle(celltype_labels{celltype})
    % for all the regions
    for region = 1:length(af_labels)
        subplot(round(sqrt(length(af_labels))),ceil(sqrt(length(af_labels))),region)
        
        imagesc(squeeze(mean(per_fish(:,region,:,:),1)))
        xlabel(af_labels(region).name,'Interpreter','None')
        set(gca,'TickLength',[0 0],'XTick',[],'YTick',[])
        colormap(viridis)
    end
    

    
    
end

%% Compare the cluster averages

close all
% get the combinations and their number
combo_vector = nchoosek(1:num_subfolders,2);
num_comb = size(combo_vector,1);

% for all celltypes
for celltype = 1:3
    figure
    % for all the combinations
    for combo = 1:num_comb
        % get the data
        average1 = average_str.(target_groups{combo_vector(combo,1)}).cluster_ave_all{celltype};
        average2 = average_str.(target_groups{combo_vector(combo,2)}).cluster_ave_all{celltype};
        
        % reshape them to 2D
        average1 = sort_traces(reshape(average1,size(average1,1),[]));
        average2 = sort_traces(reshape(average2,size(average2,1),[]));
        % 
        
        % get the correlation matrix between them and plot
        [correlation,pval] = corr(average1',average2');
        
%         % nan the non-significant values
%         correlation(pval<0.05) = NaN;

        % calculate the mean
        mean_corr = nanmean(correlation(:));
        % plot the correlation
        subplot(1,3,combo)
        imagesc(correlation)
        title(strjoin({num2str(mean_corr)},' '))
        
        sgtitle(celltype_labels{celltype})
        ylabel(target_groups{combo_vector(combo,1)})
        xlabel(target_groups{combo_vector(combo,2)})
        axis equal
        axis tight
    end
    
end
autoArrangeFigures
%% Compare the region profiles
close all
% get the combinations and their number
combo_vector = nchoosek(1:num_subfolders,2);
num_comb = size(combo_vector,1);

% for all celltypes
for celltype = 1:3
    % for all the combinations
    for combo = 1:num_comb
        figure
        % get the data
        main1 = average_str.(target_groups{combo_vector(combo,1)}).region_ave(celltype);
        main2 = average_str.(target_groups{combo_vector(combo,2)}).region_ave(celltype);
        % for all the regions
        for region = 1:length(af_labels)
            % get the current region name
            current_region = af_labels(region).name;
            % extract the region of interest
            average1 = main1.(current_region);
            average2 = main2.(current_region);
            if isempty(average1) || isempty(average2)
                continue
            end

            % reshape them to 2D
            average1 = sort_traces(reshape(average1,size(average1,1),[]));
            average2 = sort_traces(reshape(average2,size(average2,1),[]));

            % get the correlation matrix between them and plot
            [correlation,pval] = corr(average1',average2');

%             % nan the non-significant values
%             correlation(pval<0.05) = NaN;

            % plot the correlation
            subplot(4,4,region)
            imagesc(correlation)
            axis equal
            axis tight
            set(gca,'CLim',[0 1],'TickLength',[0 0],'YTick',[],'XTick',[])
            title(strjoin({target_groups{combo_vector(combo,1)},...
                target_groups{combo_vector(combo,2)},current_region},' '),...
                'Interpreter','None')
        end
        sgtitle(celltype_labels{celltype})
        % define the path and save
        file_path = strjoin({'correlationPerRegion',target_groups{combo_vector(combo,1)},...
                target_groups{combo_vector(combo,2)},celltype_labels{celltype},'.png'},'_');
        print(fullfile(figure_path,file_path), '-dpng','-r600')

        
    end
    
end
autoArrangeFigures
%% Compare cluster average on a fish by fish basis

close all

% get the number of fish in the pre group
pre_number = size(average_str.pre.meta,1);

% allocate memory for the average correlations
average_correlations = struct([]);
average_correlations(1).post = zeros(size(average_str.post.meta,1),num_regions,3);
average_correlations(1).postcontrol = zeros(size(average_str.postcontrol.meta,1),num_regions,3);

% for all the fish in the pre group
for fish = 1:pre_number
    % get the fish name
    fish_name = average_str.pre.meta(fish).fish_name;
    % see if there's a match in post or postcontrol
    match_idx = find(strcmp(fish_name,{average_str.post.meta.fish_name}));
    match_group = 'post';
    % if it's not in post, check postcontrol
    if isempty(match_idx)
        match_idx = find(strcmp(fish_name,{average_str.postcontrol.meta.fish_name}));
        match_group = 'postcontrol';
    end
    % if still no match, skip
    if isempty(match_idx)
        continue
    end
    
    % for all the celltypes
    for celltype = 2%:3
        corr_fig = figure;
        histo_fig = figure;
        % set up a counter for the matches
        match_counter = struct([]);
        
        
        % get the fish averages per region
        pre_averages = average_str.pre.region_ave_perfish{fish,celltype};
        % get the matching averages
        match_averages = average_str.(match_group).region_ave_perfish{match_idx,celltype};
        
        % compare the averages
        % for all the regions
        for region = 1:length(af_labels)
            % get the current region name
            current_region = af_labels(region).name;
            
            % check if the region is in the fish
            if contains(current_region,fields(pre_averages)) == 0 ||...
                contains(current_region,fields(match_averages)) == 0
                continue
            end
            % extract the region of interest
            average1 = pre_averages.(current_region);
            average2 = match_averages.(current_region);
            if isempty(average1) || isempty(average2)
                continue
            end

            % reshape them to 2D
            average1 = sort_traces(reshape(average1,size(average1,1),[]));
            average2 = sort_traces(reshape(average2,size(average2,1),[]));

            % get the correlation matrix between them and plot
            [correlation,pval] = corr(average1',average2','rows','pairwise');
            
            % save the average correlation
            average_correlations.(match_group)(fish,region,celltype) = nanmean(nanmean(tril(correlation,-1))); 

%             % nan the non-significant values
%             correlation(pval<0.05) = NaN;
            % select the correlation figure
            figure(corr_fig)
            % plot the correlation
            subplot(4,4,region)
            imagesc(correlation)
            axis equal
            axis tight
            set(gca,'CLim',[-1 1],'TickLength',[0 0],'YTick',[],'XTick',[])
            title(strjoin({'pre',match_group,current_region},' '),...
                'Interpreter','None')
            
            % select the histogram figure
            figure(histo_fig)
            % plot overlapping histograms
            subplot(4,4,region)
            histogram(average1(:))
            
            hold on
            histogram(average2(:))
            set(gca,'TickLength',[0 0])
            title(strjoin({'pre',match_group,current_region},' '),...
                'Interpreter','None')
            
            
        end
        % save the correlation figure
        figure(corr_fig)
        sgtitle(celltype_labels{celltype})
        colormap(viridis)
        % define the path and save
        set(gcf,'Color','w')
        file_path = fullfile(figure_path,'PerFish',strjoin({'correlationPerRegionPerFish','pre',...
            match_group,celltype_labels{celltype},fish_name,'.png'},'_'));
        export_fig(file_path,'-r600')
        

    end
    
    
end

% plot the average correlation
figure
subplot(1,2,1)
imagesc(average_correlations.post(:,:,2))
set(gca,'XTick',1:num_regions,'XTickLabels',{af_labels.name},'XtickLabelRotation',45,'TickLength',[0 0])
axis equal
axis tight
subplot(1,2,2)
imagesc(average_correlations.postcontrol(:,:,2))
set(gca,'XTick',1:num_regions,'XTickLabels',{af_labels.name},'XtickLabelRotation',45,'TickLength',[0 0])
axis equal
axis tight

autoArrangeFigures
%% Compare the cell types per fish and region

close all

% define the type of extraction
stimTypeNum = 4;
% get the number of fish in the pre group
pre_number = size(average_str.pre.meta,1);

% allocate memory for the average correlations (type,dir,on/off,pre/match)
type_counts = struct([]);
type_counts(1).post = zeros(size(average_str.post.meta,1),num_regions,3,8,2,2);
type_counts(1).postcontrol = zeros(size(average_str.postcontrol.meta,1),num_regions,3,8,2,2);

% set a fish counter for post and control
match_counter = struct([]);
match_counter(1).post = 1;
match_counter(1).postcontrol = 1;
% for all the fish in the pre group
for fish = 1:pre_number
    
    fprintf(strjoin({'Current fish:',num2str(fish),'out of',num2str(pre_number),'\r\n'},' '))
    % get the fish name
    fish_name = average_str.pre.meta(fish).fish_name;
    % see if there's a match in post or postcontrol
    match_idx = find(strcmp(fish_name,{average_str.post.meta.fish_name}));
    match_group = 'post';
    % if it's not in post, check postcontrol
    if isempty(match_idx)
        match_idx = find(strcmp(fish_name,{average_str.postcontrol.meta.fish_name}));
        match_group = 'postcontrol';
    end
    % if still no match, skip
    if isempty(match_idx)
        continue
    end
    
    % for all the celltypes
    for celltype = 1:3
%         corr_fig = figure;
%         histo_fig = figure;
%         % set up a counter for the matches
%         match_counter = struct([]);
        
        
        % get the fish averages per region
        pre_averages = average_str.pre.region_ave_perfish{fish,celltype};
        
        % get the matching averages
        match_averages = average_str.(match_group).region_ave_perfish{match_idx,celltype};
        
        % get the counts too
        pre_counts = average_str.pre.region_counts;
        match_counts = average_str.(match_group).region_counts;
        
        % compare the averages
        % for all the regions
        for region = 1:length(af_labels)
            % get the current region name
            current_region = af_labels(region).name;
            
            % check if the region is in the fish
            if contains(current_region,fields(pre_averages)) == 0 ||...
                contains(current_region,fields(match_averages)) == 0
                continue
            end
            % extract the region of interest
            average1 = pre_averages.(current_region);
            average2 = match_averages.(current_region);
            
            average1(isnan(average1)) = 0;
            average2(isnan(average2)) = 0;
%             if isempty(average1) || isempty(average2)
%                 continue
%             end
            
%             % reshape them to 2D
%             average1 = sort_traces(reshape(average1,size(average1,1),[]));
%             average2 = sort_traces(reshape(average2,size(average2,1),[]));
% 
%             % get the correlation matrix between them and plot
%             [correlation,pval] = corr(average1',average2','rows','pairwise');

            % run the svd on the averages
            [~,pre_svd] = fish_svd(average1,stimTypeNum);
            [~,match_svd] = fish_svd(average2,stimTypeNum);

            % separate the cells based on pref dir/ori and on/off
            [~,pref_dir] = max(pre_svd(:,51:58),[],2);
            on_off = (pre_svd(:,225)<pre_svd(:,226))+1;
            [~,match_dir] = max(match_svd(:,51:58),[],2);
            match_on_off = (match_svd(:,225)<match_svd(:,226))+1;
            % collapse into a matrix with only on/off and dir
            temp_pre = zeros(8,2);
            temp_match = zeros(8,2);
            
            % for all the clusters
            for clu = 1:size(pref_dir,1)
                temp_pre(pref_dir(clu),on_off(clu)) = ...
                    temp_pre(pref_dir(clu),on_off(clu)) + pre_counts{fish,celltype}(region,clu);
            end
            for clu = 1:size(match_dir,1)
                temp_match(match_dir(clu),match_on_off(clu)) = ...
                    temp_match(match_dir(clu),match_on_off(clu)) + match_counts{match_idx,celltype}(region,clu);
            end            
            
            % store the counts for each type
            type_counts.(match_group)(match_counter.(match_group),region,celltype,:,:,1) = temp_pre;
            type_counts.(match_group)(match_counter.(match_group),region,celltype,:,:,2) = temp_match;
            
        end
%         % save the correlation figure
%         figure(corr_fig)
%         sgtitle(celltype_labels{celltype})
%         colormap(viridis)
%         % define the path and save
%         set(gcf,'Color','w')
%         file_path = fullfile(figure_path,'PerFish',strjoin({'correlationPerRegionPerFish','pre',...
%             match_group,celltype_labels{celltype},fish_name,'.png'},'_'));
%         export_fig(file_path,'-r600')
        

    end
    % increase the fish counter
    match_counter.(match_group) = match_counter.(match_group) + 1;
    
    
end
%% Plot the results per fish
close all

% for both groups of fish
for groups = 1:2
    switch groups
        case 1
            match_group = 'post';
        case 2
            match_group = 'postcontrol';
    end
    % for all the fish
    for fish = 1:size(type_counts.(match_group),1)
        close all
        % get the fish name
        fish_name = average_str.pre.meta(fish).fish_name;
        % for cell types
        for celltype = 1:3
            figure
            % get the color scaling
            all_regions = type_counts.(match_group)(fish,:,celltype,:,:,:);
            c_max = log(max(all_regions(:)));
            % for all the regions
            for region = 1:num_regions
                % get the current fish
                pre = log(squeeze(type_counts.(match_group)(fish,region,celltype,:,:,1)));
                match = log(squeeze(type_counts.(match_group)(fish,region,celltype,:,:,2)));
                % plot the results
                subplot(2,num_regions,region)
                imagesc(pre)
                xlabel(strjoin({'pre',af_labels(region).name},' '),...
                    'Interpreter','None','Rotation',45,'HorizontalAlignment','right')
                set(gca,'YTick',[],'XTick',[])
                % set the color scaling
                set(gca,'CLim',[0 c_max])
                axis tight
                
                subplot(2,num_regions,region+num_regions)
                imagesc(match)
                xlabel(strjoin({match_group,af_labels(region).name},' '),...
                    'Interpreter','None','Rotation',45,'HorizontalAlignment','right')
                set(gca,'YTick',[],'XTick',[])
                axis tight
                
                % set the color scaling
                set(gca,'CLim',[0 c_max])
            end
            % set the title
            sgtitle(celltype_labels{celltype})
            
            
            % save the figure
            colormap(viridis)
            % define the path and save
            set(gcf,'Color','w')
            file_path = fullfile(figure_path,'PerDir',strjoin({'clusterCounts',...
                match_group,celltype_labels{celltype},fish_name,'.png'},'_'));
            export_fig(file_path,'-r300')
        end


    end
end
autoArrangeFigures
%% Get the deltas per fish and then average

close all

% allocate memory to save the deltas
delta_cell = cell(2,1);
% for both groups of fish
for groups = 1:2
    switch groups
        case 1
            match_group = 'post';
        case 2
            match_group = 'postcontrol';
    end
    % calculate the delta for the group
    delta_counts = type_counts.(match_group)(:,:,:,:,:,2) - ...
        type_counts.(match_group)(:,:,:,:,:,1);
    % save the deltas
    delta_cell{groups} = delta_counts;
    % calculate the average across animals
    delta_average = squeeze(nanmean(delta_counts,1));

    % for cell types
    for celltype = 1:3
        figure
        % get the scaling
        c_lim = [min(delta_average(:,celltype,:,:),[],'all'),max(delta_average(:,celltype,:,:),[],'all')];
        % for all the regions
        for region = 1:num_regions
            subplot(round(sqrt(num_regions)),ceil(sqrt(num_regions)),region)
            imagesc(squeeze(delta_average(region,celltype,:,:))')
            set(gca,'CLim',c_lim)
            xlabel(strjoin({match_group,af_labels(region).name},' '),...
                'Interpreter','None')
            set(gca,'YTick',[],'XTick',[])
            axis tight 
            % if it's the last one, add the colorbar
            if region == num_regions
                colorbar
            end
        end
        % set the title
        sgtitle(celltype_labels{celltype})
        % save the figure
        colormap(viridis)
        % define the path and save
        set(gcf,'Color','w')
        file_path = fullfile(figure_path,strjoin({'deltaCounts',...
            match_group,celltype_labels{celltype},'.png'},'_'));
        export_fig(file_path,'-r300')
    end
end
autoArrangeFigures
%% Test the difference between the control and ablated groups
close all

% define the alpha value
alpha = 0.05;

% allocate memory for the results
test_results = zeros(num_regions,3,8,2);

% for all the regions
for region = 1:num_regions
    % for all the celltypes
    for celltype = 1:3
        % for all directions
        for direction = 1:8
            % for on/off
            for onoff = 1:2
                
                % run the test across animals and save
                test_results(region,celltype,direction,onoff) = ...
                    ranksum(squeeze(delta_cell{1}(:,region,celltype,direction,onoff)), ...
                    squeeze(delta_cell{2}(:,region,celltype,direction,onoff)));
            end
        end
    end
end

% plot the results
% for all the celltypes
for celltype = 1:3
    figure
    % for all the regions
    for region = 1:num_regions
        subplot(round(sqrt(num_regions)),ceil(sqrt(num_regions)),region)
        % get the current results
        current_results = (squeeze(test_results(region,celltype,:,:)).*numel(1))<alpha;
        imagesc(current_results')
%         set(gca,'CLim',c_lim)
        xlabel(strjoin({match_group,af_labels(region).name},' '),...
            'Interpreter','None')
        set(gca,'YTick',[],'XTick',[])
        axis tight 
        % if it's the last one, add the colorbar
        if region == num_regions
            colorbar
        end
    end
    % set the title
    sgtitle(celltype_labels{celltype})
    % save the figure
    colormap(viridis)
end
figure
histogram(log(test_results(:)))