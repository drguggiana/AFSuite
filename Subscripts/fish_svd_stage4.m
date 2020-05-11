function [rep_num, trace_num, reps_all, stim_length, stim_num3, dsi_clas_final]=fish_svd_stage4(reps_all,stimTypeNum,stim_num2,time_num)


% function to calculate the svd of the stimuli involved in the
% experiment(i.e. a lot of hardcoded boundaries between stimuli and such)
% the function outputs the svd'd traces and their tuning curves, plus the
% lengths of the stimuli and the number

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
%get the number of reps
rep_num = size(reps_all,4);
%get the number of traces
trace_num = size(reps_all,1);
%allocate memory to store the independent rep data
rep_tempcell = cell(rep_num,1);
dsi_clas = zeros(trace_num,stimTypeNum-1,rep_num);
% temp = cell(rep_num,1);
%for all the reps
for reps = 1:rep_num
    %reshape the data matrix for svd
    seed_form = reshape(reps_all(:,:,:,reps),trace_num,time_num,stim_num2);
    % seed_form = reshape(conc_trace,trace_num,time_num,stim_num2);
    %allocate memory to store the length of the trace
    trace_length = zeros(stimTypeNum,2);
    
    %for the three stimulus types
    for stim_type = 1:stimTypeNum
        
        fprintf(strcat('SVD Current stim:',num2str(stim_type),'\r\n'))
        %based on the type of stimulus, define the stim numbers and the time
        %intervals to consider (the eclipses keep running after 60)
        switch stim_type
            case 1
                omr_stim = 2:9;
                time_int = 21:70;
            case 2
                omr_stim = 11:18;
                time_int = 21:70;
            case 3
                omr_stim = 19:26;
                time_int = 21:70;
            case 4
                omr_stim = [1,10];
                time_int = 21:70;
        end
        
        %extract the fraction of the data to use
        seed_split = seed_form(:,time_int,omr_stim);
        
        %allocate memory for the singular values
        svd_mat = zeros(trace_num,length(time_int));
        tc_mat = zeros(trace_num,length(omr_stim));
        dsi_mat = zeros(trace_num,2);
        osi_mat = zeros(trace_num,2);
        
        %for every cluster
        for seeds = 1:trace_num
            
            %get the vector of values
            tar_vec = squeeze(seed_split(seeds,:,:));
            [U,~,V] = svd(zscore(tar_vec,1,'all'));
            svd_mat(seeds,:) = U(:,1);
            tc_mat(seeds,:) = V(:,1);
            
            if stim_type~=4
                dsi_mat(seeds,1) = norm(sum(exp_dsi.*V(:,1)),2);
                dsi_mat(seeds,2) = norm(sum(exp_osi.*V(:,1)),2);
                boot_res = zeros(boot_num,2);
                tempV =  repmat(V(:,1),1,boot_num);
                shuff_ang =zeros(boot_num,8);
                for ii=1:boot_num
                    shuff_ang(ii,:) = randperm(8);
                end
                tempV = tempV(shuff_ang');
                boot_res(:,1) = vecnorm(sum(tempV.*exp_dsi,1),2,1);
                boot_res(:,2) = vecnorm(sum(tempV.*exp_osi,1),2,1);
                
                osi_mat(seeds,1) = sum(boot_res(:,1)>=dsi_mat(seeds,1))./boot_num;
                osi_mat(seeds,2) = sum(boot_res(:,2)>=dsi_mat(seeds,2))./boot_num;
            end
        end
        
        
        %save the matrix
        dsiosi_cell{stim_type,1} = svd_mat;
        dsiosi_cell{stim_type,2} = tc_mat;
        dsiosi_cell{stim_type,3} = dsi_mat;
        dsiosi_cell{stim_type,4} = osi_mat;
        %also store the length of the interval
        trace_length(stim_type,1) = size(svd_mat,2);
        trace_length(stim_type,2) = size(tc_mat,2);
        trace_length(stim_type,3) = size(dsi_mat,2);
        
        if stim_type~=4
            dsi_clas(osi_mat(:,1)<0.05&osi_mat(:,2)>0.05,stim_type,:) = 1;
            dsi_clas(osi_mat(:,1)>0.05&osi_mat(:,2)<0.05,stim_type,:) = 2;
        end
    end
    
    %      rep_tempcell{reps} = [horzcat(dsiosi_cell{:,1}),horzcat(dsiosi_cell{:,2})];
%     rep_tempcell{reps} = [cat(2, dsiosi_cell{1,:}),cat(2, dsiosi_cell{2,:}),cat(2, dsiosi_cell{3,:}),cat(2, dsiosi_cell{4,:})];

    rep_tempcell{reps} = horzcat([cat(2, dsiosi_cell{1,1:3}),cat(2, dsiosi_cell{2,1:3}),cat(2, dsiosi_cell{3,1:3}),cat(2, dsiosi_cell{4,1:3})],...
        [cat(2, dsiosi_cell{1,4}),cat(2, dsiosi_cell{2,4}),cat(2, dsiosi_cell{3,4}),cat(2, dsiosi_cell{4,4})]);
end
 
%assemble the new matrix of SVD'd data 
reps_all = squeeze(cat(3,rep_tempcell{:}));
% % reps_all = abs(squeeze(cat(3,rep_tempcell{:})));
% dsi_clas(reps_all(:,end-1)<0.05&reps_all(:,end)>0.05,:) = 1;
% dsi_clas(reps_all(:,end-1)>0.05&reps_all(:,end)<0.05,:) = 2;
dsi_clas_final = zeros(trace_num,1);
dsi_clas_final(any(dsi_clas==2,2))=2;
dsi_clas_final(any(dsi_clas==1,2))=1;

%define the stimuli lengths
stim_length = sum(trace_length,2)';
% stim_length = trace_length(:)';
%get the new number of stimuli
stim_num3 = length(stim_length);
