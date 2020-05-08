function [idx_clu,clu_num] = cluster_cleanup(idx_clu,clu_num,conc_trace,snr_mat,constants)
%% Use correlation to merge clusters that are too similar
% define the correlation threshold
% corr_threshold = 0.8;
corr_threshold = constants.correlation;
pval_threshold = constants.pval;
close all
%allocate memory for the averages
clu_ave = zeros(clu_num,size(conc_trace,2));
%and for the trace number
clu_number = zeros(clu_num,1);
%for all the clusters
for clu = 1:clu_num
    %calculate the cluster average
    clu_ave(clu,:) = mean(conc_trace(idx_clu==clu,:),1);
    %and store the number of traces going into each average
    clu_number(clu) = sum(idx_clu==clu);
end

% calculate the correlation matrix between clusters
[cluster_corr,pval] = corr(clu_ave');
% keep only the significant correlations
cluster_corr(pval>pval_threshold) = NaN;

% figure
% plot(clu_ave([11 3],:)')
figure
imagesc(cluster_corr)
figure
temp_corr = tril(cluster_corr,1);
histogram(temp_corr(:),'Normalization','cdf')

% get the matrix with correlation values
corr_submatrix = triu(cluster_corr,1);
% allocate memory for the new idx
new_idx = idx_clu;
% go through the clusters
for clu = 1:clu_num
    % check if this cluster has already been eliminated
    if all(isnan(corr_submatrix(clu,:)))
        continue
    end
    % get the corresponding row
    curr_corr = corr_submatrix(clu,:);
    % nan the matching cluster
    curr_corr(clu) = NaN;
    % get the clusters that match with this one above threshold
    matching_clusters = find(curr_corr>corr_threshold);
    % for all the matching clusters
    for match = 1:length(matching_clusters)
        new_idx(new_idx==matching_clusters(match)) = clu;
    end
end
% renumber the clusters consecutively
% get the new cluster numbers
new_numbers = unique(new_idx);
% get the new clu_num
new_clu_num = length(new_numbers);
% renumber
for clu = 1:new_clu_num
    new_idx(new_idx==new_numbers(clu)) = clu;
end
% redefine the cluster number and idx
idx_clu = new_idx;
clu_num = new_clu_num;
%% Exclude clusters with snr below a threshold and with too few traces

% % define the thresholds
% num_thres = constants.trace_number;
% stim_thres = constants.stimuli;
% perc_constant

%load the snr for the dataset
% snr_mat = load(name_cell{files},'snr_mat');
% snr_mat = snr_mat.snr_mat;

[idx_clu,clu_num] = cluster_snr(snr_mat,clu_num,idx_clu,constants);