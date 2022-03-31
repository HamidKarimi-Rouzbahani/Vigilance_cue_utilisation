%% Sample response processing
% developed by Hamid Karimi-Rouzbahani on 31/March/2022
clc
clear all;
rng('default') % do not play with this

Sets_of_subjects=20; % how many groups of subjects did you set at Line 56 of "TestVigCueUtilisation" script?
subjects=[2] ; % subjects you want the include in analysis
Testing_blocks=[1:30]; % blocks you want to include in analysis

dirs=dir();
% or Determine where the data is stored on PC
% dirs=dir('C:\');

% Regenerating task parameters
% Task parameters
percentage_target_cond=[0.5 0.09]; % Frequency of targets across conditions
% the length of this vector also determines how many (target frequency)
% conditions you will have
Bias_in_target_side= 0.7; % (range: 0-1) >0.5 means more target on the left vs right; 0.5 means balanced
Blocks_per_condition=15; % Number of blocks per target frequency condition

Num_of_conditions=length(percentage_target_cond);
% Generating counter-balanced blocks in terms of frequency conditions, cued
% colours and cue bias
conds_perm=repmat(perms([1:Num_of_conditions]),[Sets_of_subjects 1]);
Block_condition=nan(size(conds_perm,1),Num_of_conditions*Blocks_per_condition);
Cued_color_in_block=ones(size(conds_perm,1),Num_of_conditions*Blocks_per_condition);
for Subj=1:size(conds_perm,1)
    for cond=1:Num_of_conditions
        blk_inds_of_currnt_cond=(cond-1)*Blocks_per_condition+1:cond*Blocks_per_condition;
        Block_condition(Subj,blk_inds_of_currnt_cond)=repmat(conds_perm(Subj,cond),[1 Blocks_per_condition]);
        if length(blk_inds_of_currnt_cond)>1
            indices_to_make_zero=randsample(blk_inds_of_currnt_cond,round(size(blk_inds_of_currnt_cond,2)./2));
            Cued_color_in_block(Subj,indices_to_make_zero)=0;
        else
            Cued_color_in_block(Subj,blk_inds_of_currnt_cond)=randi([0 1],1);
        end
    end
end

% Counter-balance for Cue direction as well
Bias_in_target_side_subj(1:size(Cued_color_in_block,1),1:size(Cued_color_in_block,2))=Bias_in_target_side;
Bias_in_target_side_subj(size(Cued_color_in_block,1)+1:2*size(Cued_color_in_block,1),1:size(Cued_color_in_block,2))=1-Bias_in_target_side;
Block_condition=repmat(Block_condition,[2 1]);
Cued_color_in_block=repmat(Cued_color_in_block,[2 1]);

rand_inds=randperm(size(Bias_in_target_side_subj,1));
Bias_in_target_side_subj = Bias_in_target_side_subj(rand_inds,:);
Block_condition = Block_condition(rand_inds,:);
Cued_color_in_block = Cued_color_in_block(rand_inds,:);

%% Data preparation
Num_blks=length(Testing_blocks);
accuracy_att=nan(Num_blks,length(subjects));
TPR_att=nan(Num_blks,length(subjects));
FPR_att=nan(Num_blks,length(subjects));
TNR_att=nan(Num_blks,length(subjects));
FNR_att=nan(Num_blks,length(subjects));
correct_reaction_times_att=zeros(Num_blks,max(subjects));

for Subj=[subjects]
    for blk=Testing_blocks
        percentage_target=percentage_target_cond(1,Block_condition(Subj,blk));
        Condition_string=['Freq_',sprintf('%.2f', percentage_target)];
        for i=3:size(dirs,1)
            if strcmp(dirs(i).name,['Subj_',num2str(Subj),'_Blk_',num2str(blk),'_',Condition_string,'_test_CueUtilisation.mat'])
                load(dirs(i).name);
                break;
            end
        end
        Targ_Freq_Condition_blk=str2double(dirs(i).name(end-27:end-24));
        
        mean_sampling_time=1./60;
        for dot_num=1:Num_moving_dots*Trials_per_block
            tr=ceil(dot_num./Num_moving_dots);
            dot_in_trial=dot_num-(tr-1).*Num_moving_dots;
            
            if ~isempty(find(key_pressed1(dot_in_trial,:,tr),1))
                
                key_press_sample=find(key_pressed1(dot_in_trial,:,tr), 1, 'first');
                if isnan(distance_traj1(dot_num,key_press_sample))
                    distance_traj1(dot_num,key_press_sample)=3000;
                end
                dist_relative_to_boundary(dot_in_trial,tr)=distance_traj1(dot_num,key_press_sample)-hitting_border_distance;
            else
                dist_relative_to_boundary(dot_in_trial,tr)=nan;
            end
            distance_change_per_sample(dot_in_trial,tr)=(distance_traj1(dot_num,appearance_time(dot_in_trial,tr)+10)-distance_traj1(dot_num,appearance_time(dot_in_trial,tr)+20))./(11);
            
            if ~isempty(find(key_pressed2(dot_in_trial,:,tr),1))
                
                key_press_sample2=find(key_pressed2(dot_in_trial,:,tr), 1, 'first' );
                if isnan(distance_traj2(dot_num,key_press_sample2))
                    distance_traj2(dot_num,key_press_sample)=3000;
                end
                dist_relative_to_boundary2(dot_in_trial,tr)=distance_traj2(dot_num,key_press_sample2)-hitting_border_distance;
            else
                dist_relative_to_boundary2(dot_in_trial,tr)=nan;
            end
            distance_change_per_sample2(dot_in_trial,tr)=(distance_traj2(dot_num,appearance_time2(dot_in_trial,tr)+10)-distance_traj2(dot_num,appearance_time2(dot_in_trial,tr)+20))./(11);
        end
        
        
        distance_change_per_sample(distance_change_per_sample<0)=mean(distance_change_per_sample(distance_change_per_sample>0));
        distance_change_per_sample2(distance_change_per_sample2<0)=mean(distance_change_per_sample2(distance_change_per_sample2>0));
        
        reaction_times=((-dist_relative_to_boundary)./distance_change_per_sample).*mean_sampling_time;
        reaction_times2=((-dist_relative_to_boundary2)./distance_change_per_sample2).*mean_sampling_time;
        %% Behavioural Performance
        
        % attended
        tp_att=0;
        tn_att=0;
        fp_F_att=0;
        fp_S_att=0;
        fp_T_att=0;
        fn_att=0;
        
        g=0;
        for dot_num=1:Num_moving_dots*Trials_per_block
            tr=ceil(dot_num./Num_moving_dots);
            dot_in_trial=dot_num-(tr-1).*Num_moving_dots;
            
            if sum(dot_in_trial==top_events(:,tr))==1 && dot_color(dot_in_trial,tr)==Cued_color_in_block(Subj,blk)
                g=g+1;
                if isnan(reaction_times(dot_in_trial,tr)) && (top_events(tr)~=top_targets(tr))
                    tn_att=tn_att+1;    % number of non-target events with no resp;
                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)~=top_targets(tr) && (reaction_times(dot_in_trial,tr)<0)
                    fp_F_att=fp_F_att+1;    % number of non-target events with fast resp;
                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)~=top_targets(tr) && (reaction_times(dot_in_trial,tr)>=0)
                    fp_S_att=fp_S_att+1;    % number of non-target events with Slow resp;
                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr) && (reaction_times(dot_in_trial,tr)<0)
                    fp_T_att=fp_T_att+1;    % number of target events with Too early resp;
                elseif isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr)
                    fn_att=fn_att+1;    % number of target events with no resp;
                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr) && reaction_times(dot_in_trial,tr)>0
                    tp_att=tp_att+1;    % number of target events with resp;
                    correct_reaction_times_att(blk,Subj)=correct_reaction_times_att(blk,Subj)+reaction_times(dot_in_trial,tr);
                end
            end
            
            if sum(dot_in_trial==top_events2(:,tr))==1 && dot_color2(dot_in_trial,tr)==Cued_color_in_block(Subj,blk)
                g=g+1;
                if isnan(reaction_times2(dot_in_trial,tr)) && (top_events2(tr)~=top_targets2(tr))
                    tn_att=tn_att+1;
                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)~=top_targets2(tr) && (reaction_times2(dot_in_trial,tr)<0)
                    fp_F_att=fp_F_att+1;
                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)~=top_targets2(tr) && (reaction_times2(dot_in_trial,tr)>=0)
                    fp_S_att=fp_S_att+1;
                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr) && (reaction_times2(dot_in_trial,tr)<0)% || reaction_times2(dot_in_trial,tr)>time_to_touch_the_obstacle2(dot_in_trial,tr))
                    fp_T_att=fp_T_att+1;
                elseif isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr)
                    fn_att=fn_att+1;
                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr) && reaction_times2(dot_in_trial,tr)>0 %&& reaction_times2(dot_in_trial,tr)<time_to_touch_the_obstacle2(dot_in_trial,tr)
                    tp_att=tp_att+1;
                    correct_reaction_times_att(blk,Subj)=correct_reaction_times_att(blk,Subj)+reaction_times2(dot_in_trial,tr);
                end
            end
        end
        
        fp_att=fp_F_att+fp_S_att+fp_T_att;
        accuracy_att(blk,Subj)=(tp_att+tn_att)./(sum(top_events>0)+sum(top_events2>0));
        TPR_att(blk,Subj)=(tp_att)./(tp_att+fn_att);
        TNR_att(blk,Subj)=(tn_att)./(tn_att+fp_att);
        FPR_att(blk,Subj)=(fp_att)./(fp_att+tn_att);
        FNR_att(blk,Subj)=(fn_att)./(tp_att+fn_att);
        correct_reaction_times_att(blk,Subj)=correct_reaction_times_att(blk,Subj)./tp_att;
        % Removed the unattended dots for simplicity of the data
        
        for cond=1:length(percentage_target_cond)
            if Targ_Freq_Condition_blk==percentage_target_cond(cond)
                % accuracy
                Data{cond,1}(blk,Subj)=(tp_att+tn_att)./(sum(top_events>0)+sum(top_events2>0));
                
                % Hit rate
                Data{cond,2}(blk,Subj)=(tp_att)./(tp_att+fn_att);
                
                % True negative rate
                Data{cond,3}(blk,Subj)=(tn_att)./(tn_att+fp_att);
                
                % False alarm
                Data{cond,4}(blk,Subj)=(fp_att)./(fp_att+tn_att);
                
                % Miss
                Data{cond,5}(blk,Subj)=(fn_att)./(tp_att+fn_att);
                
                % Dprime
                Data{cond,6}(blk,Subj)=Data{cond,2}(blk,Subj)-Data{cond,4}(blk,Subj);
                
                % Reaction time
                Data{cond,7}(blk,Subj)=correct_reaction_times_att(blk,Subj);
            else
                Data{cond,1}(blk,Subj)=nan;
                Data{cond,2}(blk,Subj)=nan;
                Data{cond,3}(blk,Subj)=nan;
                Data{cond,4}(blk,Subj)=nan;
                Data{cond,5}(blk,Subj)=nan;
                Data{cond,6}(blk,Subj)=nan;
                Data{cond,7}(blk,Subj)=nan;
            end
        end
    end
end

%% Saving data as Excel file for analysis
for Subj=1:max(subjects)
    for cond=1:size(Data,1)
        if ~ismember(Subj,subjects)
            Data{cond,1}(:,Subj)=nan;
            Data{cond,2}(:,Subj)=nan;
            Data{cond,3}(:,Subj)=nan;
            Data{cond,4}(:,Subj)=nan;
            Data{cond,5}(:,Subj)=nan;
            Data{cond,6}(:,Subj)=nan;
            Data{cond,7}(:,Subj)=nan;
        end
    end
    cond=1;
    Hit_rate_condition=Data{cond,2}(:,Subj); % Hit rate in condition
    Reaction_time_condition=Data{cond,7}(:,Subj); % reaction time in condition
    T = table(Hit_rate_condition,Reaction_time_condition);
    T.Properties.VariableNames = {['Hit_rate_target_freq_',num2str(percentage_target_cond(cond))] ['RT_target_freq_',num2str(percentage_target_cond(cond))]};
    Ttotal=T;
    for cond=2:size(Data,1)
        
        Hit_rate_condition=Data{cond,2}(:,Subj); % Hit rate in condition
        Reaction_time_condition=Data{cond,7}(:,Subj); % reaction time in condition
        T = table(Hit_rate_condition,Reaction_time_condition);
        T.Properties.VariableNames = {['Hit_rate_target_freq_',num2str(percentage_target_cond(cond))] ['RT_target_freq_',num2str(percentage_target_cond(cond))]};
        Ttotal=[Ttotal T];
    end
    Bias_target_side=Bias_in_target_side_subj(Subj,:)'; % Bias inb target side: >0.5 means more targets from left
    T = table(Bias_target_side);
    T.Properties.VariableNames = {'More_targets_on_left(>0.5)_or_right(<0.5)'};
    Ttotal=[Ttotal T];
    
    filename = ['MoM_data.xlsx']; % Change the name to anything you prefer
    writetable(Ttotal,filename,'Sheet',['Subj_' num2str(Subj)])
    [Subj]
end

%% Plotting some results
RT=0; % 1 for reaction time and 0 for hit rate
if RT==0
    dataA1=(1-Data{1,2}(:,subjects))*100;
    dataB1=(1-Data{2,2}(:,subjects))*100;
else
    dataA1=Data{1,7}*1000;
    dataB1=Data{2,7}*1000;
end

MeanA=nanmean(dataA1,2);
MeanA=MeanA(~isnan(MeanA));
MeanM=nanmean(dataB1,2);
MeanM=MeanM(~isnan(MeanM));

if length(subjects)==1
    figure;
    Shad1=plot([1:length(MeanA)],MeanA,'linewidth',3);
    hold on;
    Shad2=plot([1:length(MeanM)],MeanM,'linewidth',3);
    xlabel('Block #')
    if RT==0
        ylabel({'Percentage of misses (%)'})
    else
        ylabel('Reaction time (ms)')
    end
    legend([Shad1,Shad2],{'Active','Monitoring'},'location','northwest','edgecolor','none')
else
    StdA=nanstd(dataA1');
    StdA=StdA(~isnan(StdA));
    StdM=nanstd(dataB1');
    StdM=StdM(~isnan(StdM));
    
    if exist('shadedErrorBar')~=0
        figure;
        Shad1=shadedErrorBar([1:15],MeanA,((StdA)./sqrt(length(subjects)))'*1.96,{'color',[0.1 0.1 0.9],'LineWidth',3},1);
        hold on;
        Shad2=shadedErrorBar([1:15],MeanM,((StdM)./sqrt(length(subjects)))'*1.96,{'color',[0.9 0.1 0.1],'LineWidth',3},1);
        xlabel('Block #')
        
        if RT==0
            ylabel({'Percentage of misses (%)'})
        else
            ylabel('Reaction time (ms)')
        end
        legend([Shad1.mainLine,Shad2.mainLine],{'Active','Monitoring'},'location','northwest','edgecolor','none')
        
    else
        figure;
        Shad1=plot([1:length(MeanA)],MeanA,'linewidth',3);
        hold on;
        Shad2=plot([1:length(MeanM)],MeanM,'linewidth',3);
        xlabel('Block #')
        if RT==0
            ylabel({'Percentage of misses (%)'})
        else
            ylabel('Reaction time (ms)')
        end
        legend([Shad1,Shad2],{'Active','Monitoring'},'location','northwest','edgecolor','none')
    end
end