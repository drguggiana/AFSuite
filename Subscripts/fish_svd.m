function [trace_num, seed_output, stim_length, stim_num3, dsi_clas_final]=fish_svd(seed_form,stimTypeNum,varargin)
% function to calculate the svd of the stimuli involved in the
% experiment(i.e. a lot of hardcoded boundaries between stimuli and such)
% the function outputs the svd'd traces and their tuning curves, plus the
% lengths of the stimuli and the number

% define verbosity
if length(varargin)>=1
    verbose = varargin{1};
else
    verbose = 0;
end

% get the classification if already provided
if length(varargin) >= 2
    dsi_clas_final = varargin{2};
    class_flag = 1;
else
    class_flag = 0;
end

% define the number of bootstraps
boot_num = 1000;
%define the set of angles
angle_vec = [135 180 225 270 315 360 45 90];
angle_rad = deg2rad(angle_vec);

%get the number of angles
ang_num = length(angle_vec);
%create the set of exponentials
exp_dsi = zeros(ang_num,1);
exp_osi = zeros(ang_num,1);
%for all the angles
for ang = 1:ang_num
    %calculate the exponential
    exp_dsi(ang) = exp(1i.*angle_rad(ang));
    %calculate the osi exponential
    exp_osi(ang) = exp(2.*1i.*angle_rad(ang));
end
 
%allocate memory for the output
dsiosi_cell = cell(stimTypeNum,4);
% %get the number of reps
% rep_num = size(reps_all,4);
%get the number of traces
trace_num = size(seed_form,1);

% permute seed_form to avoid using squeeze
seed_form = permute(seed_form,[2,3,1]);
% allocate memory for the ROI classification
dsi_clas = zeros(trace_num,stimTypeNum-1);

%allocate memory to store the length of the trace
trace_length = zeros(stimTypeNum,2);

%for the three stimulus types
for stim_type = 1:stimTypeNum
    
    % if verbosity is on
    if verbose
        fprintf(strcat('SVD Current stim:',num2str(stim_type),'\r\n'))
    end
    %based on the type of stimulus, define the stim numbers and the time
    %intervals to consider (the eclipses keep running after 60)
    switch stim_type
        case 1
            target_stim = 2:9;
            time_int = 21:70;
        case 2
            target_stim = 11:18;
            time_int = 21:70;
        case 3
            target_stim = 19:26;
            time_int = 21:70;
        case 4
            target_stim = [1,10];
            time_int = 21:70;
    end
    
    %extract the fraction of the data to use
    seed_split = reshape(seed_form(time_int,target_stim,:),length(time_int)*length(target_stim),[]);
    
    % calculate the zscore across the non-ROI dimension and reshape back
    seed_split = reshape(zscore(seed_split,1,1),length(time_int),length(target_stim),[]);
    
    %allocate memory for the singular values
    svd_mat = zeros(trace_num,length(time_int));
    tc_mat = zeros(trace_num,length(target_stim));
    real_dsiosi = zeros(trace_num,2);
    boot_dsiosi = zeros(trace_num,2);
    
    % only get the permutations if the classification is not given
    if class_flag == 0
        % allocate memory for the shuffles
        shuff_ang = zeros(boot_num,8);
        % generate the shuffles
        for ii=1:boot_num
            shuff_ang(ii,:) = randperm(8);
        end
    end
    %for every seed
    for seeds = 1:trace_num
        
        %get the vector of values
        [U,~,V] = svd(seed_split(:,:,seeds));
        svd_mat(seeds,:) = U(:,1);
        tc_mat(seeds,:) = V(:,1);
        % if it's not the on off, calculate DSI, OSI
        if stim_type~=4
            % calculate the DSI and OSI
            real_dsiosi(seeds,1) = norm(sum(exp_dsi.*V(:,1)),2);
            real_dsiosi(seeds,2) = norm(sum(exp_osi.*V(:,1)),2);
            
            % only shuffle if the class is not provided
            if class_flag == 0
                % allocate memory for the bootstrap result
                boot_res = zeros(boot_num,2);
                % allocate memory for executing the shuffle
                tempV = repmat(V(:,1),1,boot_num);

                % shuffle the angles
                tempV = tempV(shuff_ang');
                % calculate the shuffle result
                boot_res(:,1) = vecnorm(sum(tempV.*exp_dsi,1),2,1);
                boot_res(:,2) = vecnorm(sum(tempV.*exp_osi,1),2,1);
                % calculate the bootstrapped DSI and OSI
                boot_dsiosi(seeds,1) = sum(boot_res(:,1)>=real_dsiosi(seeds,1))./boot_num;
                boot_dsiosi(seeds,2) = sum(boot_res(:,2)>=real_dsiosi(seeds,2))./boot_num;
            end
        end
    end
    
    
    %save the matrix
    dsiosi_cell{stim_type,1} = horzcat(svd_mat,tc_mat);
%     dsiosi_cell{stim_type,2} = tc_mat;

    %also store the length of the interval
    trace_length(stim_type,1) = size(svd_mat,2)+size(tc_mat,2);
%     trace_length(stim_type,2) = size(tc_mat,2);
    % if it's not on off, store the classification of the ROI per stim
    if stim_type~=4
        trace_length(stim_type,2) = size(real_dsiosi,2);
        dsiosi_cell{stim_type,2} = real_dsiosi;
        dsiosi_cell{stim_type,3} = boot_dsiosi;
        dsi_clas(boot_dsiosi(:,1)<0.05&boot_dsiosi(:,2)>0.05,stim_type,:) = 1;
        dsi_clas(boot_dsiosi(:,1)>0.05&boot_dsiosi(:,2)<0.05,stim_type,:) = 2;
    end
end
    
% assemble the feature output
seed_output = cat(2,dsiosi_cell{:,1:2});
% seed_output = dsiosi_cell;

% classify the ROI regardless of stimulus if not already provided
if class_flag == 0
    dsi_clas_final = zeros(trace_num,1);
    dsi_clas_final(any(dsi_clas==2,2))=2;
    dsi_clas_final(any(dsi_clas==1,2))=1;
end

%define the stimuli lengths
% stim_length = sum(trace_length,2)';
stim_length = vertcat(trace_length(:,1),sum(trace_length(:,2)))';
%get the new number of stimuli
stim_num3 = length(stim_length);