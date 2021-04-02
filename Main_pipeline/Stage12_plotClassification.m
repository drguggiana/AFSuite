%% Plot predictions
clearvars
close all
Paths
%% Load the target file

target_path = uipickfiles('FilterSpec',classification_path);

file_contents = load(target_path{1});

confusion_cell = file_contents.confusion_cell;
params = file_contents.params;

% load the labels
labels = load(constants_path,'labels');
labels = labels.labels;

af_labels = labels.af;
% add the all regions label
af_labels(14).name = 'All';
celltype_labels = labels.celltype;


clear file_contents
%% Calculate the average confusion matrices for each region

% get the number of fish
num_fish = size(confusion_cell,1);
% get the number of celltypes
num_celltype = size(confusion_cell,2);
% get the number of regions
num_regions = size(confusion_cell,3);

% get the number of categories
switch params.subsample_flag
    case 0 
        num_cat = 26;
    case 1 % stimulus categories
        num_cat = 4;
    case 2 % stimulus and time
        num_cat = 16;
        
        
end

% allocate memory for the averages across fish for each region and celltype
average_matrices = zeros(num_regions,num_celltype,num_cat,num_cat);
sem_matrices = zeros(num_regions,num_celltype,num_cat,num_cat);
% for all regions
for region = 1:num_regions
    % for all the celltypes
    for celltype = 1:num_celltype
        % allocate a matrix to get the average
        average_conf = zeros(num_cat,num_cat,num_fish);
        % for all the fish
        for fish = 1:num_fish
            % see if there's a nan
            if ~iscell(confusion_cell{fish,celltype,region})
                average_conf(:,:,fish) = NaN;
                continue
            end
            % get the reps
            current_matrix = nanmean(cat(3,confusion_cell{fish,celltype,region}{:,2}),3);
            average_conf(:,:,fish) = current_matrix;
        end
        average_matrices(region,celltype,:,:) = nanmean(average_conf,3);
        sem_matrices(region,celltype,:,:) = nanstd(average_conf,0,3)./sqrt(num_fish);
    end
end
%% Plot confusion matrices

close all

% for all celltypes
for celltype = 1:num_celltype
    figure
    % for all the regions
    for region = 1:num_regions
        subplot(round(sqrt(num_regions)),ceil(sqrt(num_regions)),region)
        % get the matrix
        current_matrix = squeeze(average_matrices(region,celltype,:,:));
        % normalize to rows
        current_matrix = current_matrix./sum(current_matrix,2);
        imagesc(current_matrix)
        axis square
        set(gca,'CLim',[0 1])
        set(gca,'XTick',1:4,'XTickLabel',{'ON','OMR','BE','WE'},...
            'XTickLabelRotation',45,'TickLabelInterpreter','None',...
            'YTick',[])
        ylabel(af_labels(region).name,'Interpreter','None')

        
    end
    set(gcf,'Color','w')
    sgtitle(celltype_labels{celltype})
end
autoArrangeFigures
%% Plot per stimulus performance

close all

% for all celltypes
for celltype = 1:num_celltype
    figure
    % for all the regions
    for region = 1:num_regions
        subplot(round(sqrt(num_regions)),ceil(sqrt(num_regions)),region)
        % get the matrix
        ref_matrix = squeeze(average_matrices(region,celltype,:,:));
        current_sem = squeeze(sem_matrices(region,celltype,:,:));
        % normalize to rows
        current_matrix = ref_matrix./sum(ref_matrix,2);
        current_sem = current_sem./sum(ref_matrix,2);
        
        
%         plot(1:num_cat,diag(current_matrix),'ko')
        errorbar(1:num_cat,diag(current_matrix),diag(current_sem),'ko')
%         axis square
%         set(gca,'CLim',[0 1])
        set(gca,'XTick',1:4,'XTickLabel',{'ON','OMR','BE','WE'},...
            'XTickLabelRotation',45,'TickLabelInterpreter','None',...
            'XLim',[0 num_cat+1],'YLim',[0 1])
        ylabel(af_labels(region).name,'Interpreter','None')
        hold on
        plot(get(gca,'XLim'),[0.25 0.25],'r--')
        
    end
    set(gcf,'Color','w')
    sgtitle(celltype_labels{celltype})
end
autoArrangeFigures
%% Plot per stimulus performance split by stimuli

close all

% for all celltypes
for celltype = 1:num_celltype
    figure
    % for all the stimuli
    for stim = 1:num_cat
        subplot(round(sqrt(num_cat)),ceil(sqrt(num_cat)),stim)
        % get the matrix
        ref_matrix = squeeze(average_matrices(:,celltype,stim,:));
        current_sem = squeeze(sem_matrices(:,celltype,stim,:));
        % normalize to rows
        current_matrix = ref_matrix(:,stim)./sum(ref_matrix,2);
        current_sem = current_sem(:,stim)./sum(ref_matrix,2);
        
        
%         plot(1:num_cat,diag(current_matrix),'ko')
        errorbar(1:num_regions,current_matrix,current_sem,'ko')
%         axis square
%         set(gca,'CLim',[0 1])
        set(gca,'XTick',1:num_regions,'XTickLabel',{af_labels.name},...
            'XTickLabelRotation',45,'TickLabelInterpreter','None',...
            'XLim',[0 num_regions+1],'YLim',[0 1])
        ylabel(num2str(stim),'Interpreter','None')
        hold on
        plot(get(gca,'XLim'),[0.25 0.25],'r--')
        
    end
    set(gcf,'Color','w')
    sgtitle(celltype_labels{celltype})
end
autoArrangeFigures