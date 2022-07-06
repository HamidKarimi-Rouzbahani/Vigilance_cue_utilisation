%% Sample response processing
% developed by Hamid Karimi-Rouzbahani on 15/June/2022
% Modified by Hamid Karimi-Rouzbahani on 6/July/2022

clc
clear all;

subjects=[4] ; % subjects you want the include in analysis
percentages=[0.5 0.09]; % Frequency of targets across conditions
chunks = [1:6];
Testing_blocks=[1:5]; % blocks you want to include in analysis
for p=1:length(percentages)
    for i=1:7
        Data{p,i}=nan(length(chunks)*length(Testing_blocks),max(subjects));
    end
end
dirs=dir();
% or Determine where the data is stored on PC
% dirs=dir('C:\');
Num_moving_dots=2;
Trials_per_block=32;

%% Data preparation
for Subj=subjects
    for chunk=chunks
        for blk=Testing_blocks
            for cndss=percentages
                %                 cndss=percentages(cns);
                %                 [percentages]
                correct_reaction_times_att=0;
                %                 processs=0;
                for i=3:size(dirs,1)
                    if  strcmp(dirs(i).name(end-3:end),'.mat') && strcmp(dirs(i).name(1:26+length(num2str(Subj))+length(num2str(chunk))+length(num2str(blk))+1),['Subj_',num2str(Subj),'_Blk_',num2str(blk),'_Chunk_',num2str(chunk),'_Freq_',sprintf('%.2f',cndss)])
                        load(dirs(i).name);
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
                                    correct_reaction_times_att=correct_reaction_times_att+reaction_times(dot_in_trial,tr);
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
                                    correct_reaction_times_att=correct_reaction_times_att+reaction_times2(dot_in_trial,tr);
                                end
                            end
                        end
                        
                        fp_att=fp_F_att+fp_S_att+fp_T_att;
                        correct_reaction_times_att=correct_reaction_times_att./tp_att;
                        % Removed the unattended dots for simplicity of the data
                        
                        [~,cond]=ismember(cndss,percentages);
                        [chunk blk cndss cond]
%                         if chunk>round(length(chunks)/2)
%                             chunkc=chunk-round(length(chunks)/2);
%                         else
                            chunkc=chunk;                            
%                         end
                        if Targ_Freq_Condition_blk==cndss
                            % Accuracy
                            Data{cond,1}((chunkc-1)*length(Testing_blocks)+blk,Subj)=(tp_att+tn_att)./(sum(top_events>0)+sum(top_events2>0));
                            
                            % Hit rate
                            Data{cond,2}((chunkc-1)*length(Testing_blocks)+blk,Subj)=(tp_att)./(tp_att+fn_att);
                            
                            % True negative rate
                            Data{cond,3}((chunkc-1)*length(Testing_blocks)+blk,Subj)=(tn_att)./(tn_att+fp_att);
                            
                            % False alarm
                            Data{cond,4}((chunkc-1)*length(Testing_blocks)+blk,Subj)=(fp_att)./(fp_att+tn_att);
                            
                            % Miss
                            Data{cond,5}((chunkc-1)*length(Testing_blocks)+blk,Subj)=(fn_att)./(tp_att+fn_att);
                            
                            % Dprime
                            Data{cond,6}((chunkc-1)*length(Testing_blocks)+blk,Subj)=Data{cond,2}(blk,Subj)-Data{cond,4}(blk,Subj);
                            
                            % Reaction time
                            Data{cond,7}((chunkc-1)*length(Testing_blocks)+blk,Subj)=correct_reaction_times_att;
                        else
                            Data{cond,1}((chunkc-1)*length(Testing_blocks)+blk,Subj)=nan;
                            Data{cond,2}((chunkc-1)*length(Testing_blocks)+blk,Subj)=nan;
                            Data{cond,3}((chunkc-1)*length(Testing_blocks)+blk,Subj)=nan;
                            Data{cond,4}((chunkc-1)*length(Testing_blocks)+blk,Subj)=nan;
                            Data{cond,5}((chunkc-1)*length(Testing_blocks)+blk,Subj)=nan;
                            Data{cond,6}((chunkc-1)*length(Testing_blocks)+blk,Subj)=nan;
                            Data{cond,7}((chunkc-1)*length(Testing_blocks)+blk,Subj)=nan;
                        end
                    end
                end
            end
        end
    end
end
%% Saving data as Excel file for analysis
for Subj=subjects
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
    Mean_Hit_rate=nanmean(Hit_rate_condition);
    Reaction_time_condition=Data{cond,7}(:,Subj); % reaction time in condition
    Mean_Reaction_time=nanmean(Reaction_time_condition);
    
    HR=[Hit_rate_condition;nan(5,1);Mean_Hit_rate];
    RT=[Reaction_time_condition;nan(5,1);Mean_Reaction_time];
    T = table(HR,RT);
    T.Properties.VariableNames = {['Hit_rate_target_freq_',num2str(percentages(cond)*100)] ['RT_target_freq_',num2str(percentages(cond)*100)]};
    Ttotal=T;
    Data_csv_total=[HR RT];
    
    for cond=2:size(Data,1)
        
        Hit_rate_condition=Data{cond,2}(:,Subj); % Hit rate in condition
        Mean_Hit_rate=nanmean(Hit_rate_condition);
        Reaction_time_condition=Data{cond,7}(:,Subj); % reaction time in condition
        Mean_Reaction_time=nanmean(Reaction_time_condition);
        HR=[Hit_rate_condition;nan(5,1);Mean_Hit_rate];
        RT=[Reaction_time_condition;nan(5,1);Mean_Reaction_time];
        T = table(HR,RT);
        T.Properties.VariableNames = {['Hit_rate_target_freq_',num2str(percentages(cond)*100)] ['RT_target_freq_',num2str(percentages(cond)*100)]};
        Ttotal=[Ttotal T];
        Data_csv_total=horzcat(Data_csv_total,[HR RT]);
    end
    
    filename = ['MoM_data_CU_manually_counter_balanced.xlsx']; % Change the name to anything you prefer
    writetable(Ttotal,filename,'Sheet',['Subj_' num2str(Subj)])
       
    [Subj]
end

%% Plotting some results
RT=0; % 1 for reaction time and 0 for hit rate
if RT==0
    dataA1=(Data{1,2}(:,subjects))*100;
    dataB1=(Data{2,2}(:,subjects))*100;
else
    dataA1=Data{1,7}(:,subjects)*1000;
    dataB1=Data{2,7}(:,subjects)*1000;
end

Mean1=nanmean(dataA1,2);
Mean1=Mean1(~isnan(Mean1));

Mean2=nanmean(dataB1,2);
Mean2=Mean2(~isnan(Mean2));

figure;
Shad1=plot([1:length(Mean1)],Mean1,'linewidth',3);
hold on;
Shad2=plot([1:length(Mean2)],Mean2,'linewidth',3);
xlabel('Block #')
if RT==0
    ylabel({'Percentage of hits (%)'})
else
    ylabel('Reaction time (ms)')
end
legend([Shad1,Shad2],{['Target Freq. = ',num2str(percentages(1))],['Target Freq. = ',num2str(percentages(2))]},'location','northwest','edgecolor','none')
