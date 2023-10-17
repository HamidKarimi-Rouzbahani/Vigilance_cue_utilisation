%% Sample response processing
% developed by Hamid Karimi-Rouzbahani on 15/June/2022
% Modified by Hamid Karimi-Rouzbahani on 6/July/2022
% Modified by Hamid Karimi-Rouzbahani on 4/September/2023
% Modified by Hamid Karimi-Rouzbahani on 17/September/2023 to produce FA as
% well
% Modified by Hamid Karimi-Rouzbahani on 24/September/2023 to separate
% performance on cued and uncued sides

% Modified by Hamid Karimi-Rouzbahani on 17/October/2023 to fix some errors
% with data collection: see the last lines in this section

clc
clear all;
addpath(genpath('F:\RESEARCH\Hamid\Multicentre dataset\Scripts\bayesFactor-master'))
subjects=[1:92]; % subjects you want the include in analysis
percentages=[0.5 0.09]; % Frequency of targets across conditions
chunks = [1:6]; % number of chunks
Testing_blocks=[1:5]; % blocks you want to include in analysis
for p=1:length(percentages)
    for i=1:6 % performance metrics including hit rate, FA and RT
        Data{p,i}=nan(length(chunks)*length(Testing_blocks),max(subjects));
    end
end
Num_moving_dots=2;
Trials_per_block=32;
load('Bias_in_target_side_subj.mat')
Bias_in_target_side_subj(40,:)=Bias_in_target_side_subj(41,:); % correct fot subject #40 who was labeled as #41
Bias_in_target_side_subj(89,:)=Bias_in_target_side_subj(87,:); % correct fot subject #89 who was labeled as #87
%% Data preparation
for Subj=subjects
    address=[['F:\RESEARCH\Hamid\Anina\Zaid\OneDrive_2023-09-03\All data\REXP00',sprintf( '%02d', Subj )]];
    dirs=dir(address);
    for chunk=chunks
        for blk=Testing_blocks
            for cndss=percentages
                correct_reaction_times_att_left=0;
                correct_reaction_times_att_right=0;
                for i=3:size(dirs,1)
                    if  strcmp(dirs(i).name(end-3:end),'.mat') && strcmp(dirs(i).name(1:26+length(num2str(Subj))+length(num2str(chunk))+length(num2str(blk))+1),['Subj_',num2str(Subj),'_Blk_',num2str(blk),'_Chunk_',num2str(chunk),'_Freq_',sprintf('%.2f',cndss)])
                        load([address,'\',dirs(i).name]);
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
                        tp_att_left=0;
                        tn_att_left=0;
                        fp_F_att_left=0;
                        fp_S_att_left=0;
                        fp_T_att_left=0;
                        fn_att_left=0;

                        tp_att_right=0;
                        tn_att_right=0;
                        fp_F_att_right=0;
                        fp_S_att_right=0;
                        fp_T_att_right=0;
                        fn_att_right=0;

                        g=0;
                        f=0;
                        for dot_num=1:Num_moving_dots*Trials_per_block
                            tr=ceil(dot_num./Num_moving_dots);
                            dot_in_trial=dot_num-(tr-1).*Num_moving_dots;

                            if sum(dot_in_trial==top_events(:,tr))==1 && dot_color(dot_in_trial,tr)==Cued_color_in_block(Subj,blk)
                                if isnan(reaction_times(dot_in_trial,tr)) && (top_events(tr)~=top_targets(tr))
                                    tn_att_left=tn_att_left+1;    % number of non-target events with no resp;
                                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)~=top_targets(tr) && (reaction_times(dot_in_trial,tr)<0)
                                    fp_F_att_left=fp_F_att_left+1;    % number of non-target events with fast resp;
                                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)~=top_targets(tr) && (reaction_times(dot_in_trial,tr)>=0)
                                    fp_S_att_left=fp_S_att_left+1;    % number of non-target events with Slow resp;
                                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr) && (reaction_times(dot_in_trial,tr)<0)
                                    fp_T_att_left=fp_T_att_left+1;    % number of target events with Too early resp;
                                    f=f+1;
                                elseif isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr)
                                    fn_att_left=fn_att_left+1;    % number of target events with no resp;
                                    f=f+1;
                                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr) && reaction_times(dot_in_trial,tr)>0
                                    tp_att_left=tp_att_left+1;    % number of target events with resp;
                                    f=f+1;
                                    correct_reaction_times_att_left=correct_reaction_times_att_left+reaction_times(dot_in_trial,tr);
                                end
                            end
                            
                            if sum(dot_in_trial==top_events2(:,tr))==1 && dot_color2(dot_in_trial,tr)==Cued_color_in_block(Subj,blk)
                                if isnan(reaction_times2(dot_in_trial,tr)) && (top_events2(tr)~=top_targets2(tr))
                                    tn_att_right=tn_att_right+1;
                                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)~=top_targets2(tr) && (reaction_times2(dot_in_trial,tr)<0)
                                    fp_F_att_right=fp_F_att_right+1;
                                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)~=top_targets2(tr) && (reaction_times2(dot_in_trial,tr)>=0)
                                    fp_S_att_right=fp_S_att_right+1;
                                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr) && (reaction_times2(dot_in_trial,tr)<0)% || reaction_times2(dot_in_trial,tr)>time_to_touch_the_obstacle2(dot_in_trial,tr))
                                    fp_T_att_right=fp_T_att_right+1;
                                    g=g+1;
                                elseif isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr)
                                    g=g+1;
                                    fn_att_right=fn_att_right+1;
                                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr) && reaction_times2(dot_in_trial,tr)>0 %&& reaction_times2(dot_in_trial,tr)<time_to_touch_the_obstacle2(dot_in_trial,tr)
                                    tp_att_right=tp_att_right+1;
                                    g=g+1;
                                    correct_reaction_times_att_right=correct_reaction_times_att_right+reaction_times2(dot_in_trial,tr);
                                end
                            end
                        end

                        fp_att_left=fp_F_att_left+fp_S_att_left+fp_T_att_left;
                        correct_reaction_times_att_left=correct_reaction_times_att_left./tp_att_left;

                        fp_att_right=fp_F_att_right+fp_S_att_right+fp_T_att_right;
                        correct_reaction_times_att_right=correct_reaction_times_att_right./tp_att_right;
                        
                        [~,cond]=ismember(cndss,percentages);

                        chunkc=chunk;
                        blk_all=(chunkc-1)*length(Testing_blocks)+blk;
                        [Subj chunk blk blk_all]

                        if Targ_Freq_Condition_blk==cndss
                            if Bias_in_target_side_subj(Subj,blk)>0.5
                                % Hit rate left
                                Data{cond,1}(blk_all,Subj)=(tp_att_left)./(tp_att_left+fn_att_left);

                                % Hit rate right
                                Data{cond,2}(blk_all,Subj)=(tp_att_right)./(tp_att_right+fn_att_right);

                                % False alarm left
                                Data{cond,3}(blk_all,Subj)=(fp_att_left)./(fp_att_left+tn_att_left);

                                % False alarm right
                                Data{cond,4}(blk_all,Subj)=(fp_att_right)./(fp_att_right+tn_att_right);

                                % Reaction time left
                                Data{cond,5}(blk_all,Subj)=correct_reaction_times_att_left;

                                % Reaction time right
                                Data{cond,6}(blk_all,Subj)=correct_reaction_times_att_right;                            
                            else
                                % Hit rate right
                                Data{cond,1}(blk_all,Subj)=(tp_att_right)./(tp_att_right+fn_att_right);

                                % Hit rate left
                                Data{cond,2}(blk_all,Subj)=(tp_att_left)./(tp_att_left+fn_att_left);

                                % False alarm right
                                Data{cond,3}(blk_all,Subj)=(fp_att_right)./(fp_att_right+tn_att_right);

                                % False alarm left
                                Data{cond,4}(blk_all,Subj)=(fp_att_left)./(fp_att_left+tn_att_left);

                                % Reaction time right
                                Data{cond,5}(blk_all,Subj)=correct_reaction_times_att_right;

                                % Reaction time left
                                Data{cond,6}(blk_all,Subj)=correct_reaction_times_att_left;
                            end
                        else
                            Data{cond,1}(blk_all,Subj)=nan;
                            Data{cond,2}(blk_all,Subj)=nan;
                            Data{cond,3}(blk_all,Subj)=nan;
                            Data{cond,4}(blk_all,Subj)=nan;
                            Data{cond,5}(blk_all,Subj)=nan;
                            Data{cond,6}(blk_all,Subj)=nan;
                        end
                    end
                end
            end
        end
    end    
end
ccc
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
        end
    end
    cond=1;
    Hit_rate_conditionB=Data{cond,1}(:,Subj); % Hit rate in cued
    Mean_Hit_rate_allB=nanmean(Hit_rate_conditionB);
    Hit_rate_conditionU=Data{cond,2}(:,Subj); % Hit rate in uncued
    Mean_Hit_rate_allU=nanmean(Hit_rate_conditionU);

    FA_rate_conditionB=Data{cond,3}(:,Subj); % False alarm rate in cued
    Mean_FA_rate_allB=nanmean(FA_rate_conditionB);
    FA_rate_conditionU=Data{cond,4}(:,Subj); % False alarm rate in uncued
    Mean_FA_rate_allU=nanmean(FA_rate_conditionU);

    Reaction_time_conditionB=Data{cond,5}(:,Subj); % Reaction time in cued
    Mean_Reaction_time_allB=nanmean(Reaction_time_conditionB);
    Reaction_time_conditionU=Data{cond,6}(:,Subj); % Reaction time in uncued
    Mean_Reaction_time_allU=nanmean(Reaction_time_conditionU);

    HRB=[Hit_rate_conditionB;nan(5,1);Mean_Hit_rate_allB];
    FAB=[FA_rate_conditionB;nan(5,1);Mean_FA_rate_allB];
    RTB=[Reaction_time_conditionB;nan(5,1);Mean_Reaction_time_allB];

    HRU=[Hit_rate_conditionU;nan(5,1);Mean_Hit_rate_allU];
    FAU=[FA_rate_conditionU;nan(5,1);Mean_FA_rate_allU];
    RTU=[Reaction_time_conditionU;nan(5,1);Mean_Reaction_time_allU];

    T = table(HRB,HRU,RTB,RTU,FAB,FAU);
    T.Properties.VariableNames = {['Hit_targ_freq_',num2str(percentages(cond)*100),'_cued'] ['Hit_targ_freq_',num2str(percentages(cond)*100),'_uncued'] ['RT_targ_freq_',num2str(percentages(cond)*100),'_cued'] ['RT_targ_freq_',num2str(percentages(cond)*100),'_uncued'] ['FA_targ_freq_',num2str(percentages(cond)*100),'_cued'] ['FA_targ_freq_',num2str(percentages(cond)*100),'_uncued'] };
    Ttotal=T;

    for cond=2:size(Data,1)

        Hit_rate_conditionB=Data{cond,1}(:,Subj); % Hit rate in cued
        Mean_Hit_rate_allB=nanmean(Hit_rate_conditionB);
        Hit_rate_conditionU=Data{cond,2}(:,Subj); % Hit rate in uncued
        Mean_Hit_rate_allU=nanmean(Hit_rate_conditionU);

        FA_rate_conditionB=Data{cond,3}(:,Subj); % False alarm rate in cued
        Mean_FA_rate_allB=nanmean(FA_rate_conditionB);
        FA_rate_conditionU=Data{cond,4}(:,Subj); % False alarm rate in uncued
        Mean_FA_rate_allU=nanmean(FA_rate_conditionU);

        Reaction_time_conditionB=Data{cond,5}(:,Subj); % Reaction time in cued
        Mean_Reaction_time_allB=nanmean(Reaction_time_conditionB);
        Reaction_time_conditionU=Data{cond,6}(:,Subj); % Reaction time in uncued
        Mean_Reaction_time_allU=nanmean(Reaction_time_conditionU);

        HRB=[Hit_rate_conditionB;nan(5,1);Mean_Hit_rate_allB];
        FAB=[FA_rate_conditionB;nan(5,1);Mean_FA_rate_allB];
        RTB=[Reaction_time_conditionB;nan(5,1);Mean_Reaction_time_allB];

        HRU=[Hit_rate_conditionU;nan(5,1);Mean_Hit_rate_allU];
        FAU=[FA_rate_conditionU;nan(5,1);Mean_FA_rate_allU];
        RTU=[Reaction_time_conditionU;nan(5,1);Mean_Reaction_time_allU];
        T = table(HRB,HRU,RTB,RTU,FAB,FAU);
        T.Properties.VariableNames = {['Hit_targ_freq_',num2str(percentages(cond)*100),'_cued'] ['Hit_targ_freq_',num2str(percentages(cond)*100),'_uncued'] ['RT_targ_freq_',num2str(percentages(cond)*100),'_cued'] ['RT_targ_freq_',num2str(percentages(cond)*100),'_uncued'] ['FA_targ_freq_',num2str(percentages(cond)*100),'_cued'] ['FA_targ_freq_',num2str(percentages(cond)*100),'_uncued'] };
        Ttotal=[Ttotal T];

    end

    Tlab=table(num2str([nan(1,35) 'A']'));
    Tlab.Properties.VariableNames = {['Average']};
    Ttotal2=[Ttotal Tlab];

    filename = ['MoM_data_CU_manually_counter_balanced_cue_side_separated.xlsx']; % Change the name to anything you prefer
    writetable(Ttotal2,filename,'Sheet',['Subj_' num2str(Subj)])

    [Subj]
end
ccc
%% Plotting the results
plott=3; % 1= Hit Rate; 2= False Alarm; 3=RT
gca = axes('Position',[0.22 0.25 0.775 0.72]);
if plott==1
    miny=20;
    maxy=100;
    yticks=[0:20:100];
    data_to_plot1=(Data{1,1}(:,subjects))*100;
    data_to_plot2=(Data{2,1}(:,subjects))*100;
    data_to_plot3=(Data{1,2}(:,subjects))*100;
    data_to_plot4=(Data{2,2}(:,subjects))*100;
elseif plott==2
    miny=0;
    maxy=40;
    yticks=[0:20:40];
    data_to_plot1=(Data{1,3}(:,subjects))*100;
    data_to_plot2=(Data{2,3}(:,subjects))*100;
    data_to_plot3=(Data{1,4}(:,subjects))*100;
    data_to_plot4=(Data{2,4}(:,subjects))*100;    
elseif plott==3
    miny=200;
    maxy=400;
    yticks=[200:50:400];
    data_to_plot1=(Data{1,5}(:,subjects))*1000;
    data_to_plot2=(Data{2,5}(:,subjects))*1000;
    data_to_plot3=(Data{1,6}(:,subjects))*1000;
    data_to_plot4=(Data{2,6}(:,subjects))*1000;    
end


MeanAc=(nanmean([data_to_plot1(1:15,:) data_to_plot1(16:30,:)],2));
StdAc=nanstd([data_to_plot1(1:15,:) data_to_plot1(16:30,:)]');
MeanMc=(nanmean([data_to_plot2(1:15,:) data_to_plot2(16:30,:)],2));
StdMc=nanstd([data_to_plot2(1:15,:) data_to_plot2(16:30,:)]');

MeanAu=(nanmean([data_to_plot3(1:15,:) data_to_plot3(16:30,:)],2));
StdAu=nanstd([data_to_plot3(1:15,:) data_to_plot3(16:30,:)]');
MeanMu=(nanmean([data_to_plot4(1:15,:) data_to_plot4(16:30,:)],2));
StdMu=nanstd([data_to_plot4(1:15,:) data_to_plot4(16:30,:)]');

plots(1)=shadedErrorBar([1:15],MeanAc',((StdAc)./sqrt(size(data_to_plot1,2)))*1.96,'lineprops',{'b','markerfacecolor','b','LineWidth',2},'transparent',1);
hold on;
plots(2)=shadedErrorBar([1:15],MeanMc',((StdMc)./sqrt(size(data_to_plot2,2)))*1.96,'lineprops',{'r','markerfacecolor','r','LineWidth',2},'transparent',1);

plots(3)=shadedErrorBar([1:15],MeanAu',((StdAc)./sqrt(size(data_to_plot3,2)))*1.96,'lineprops',{'--b','markerfacecolor','b','LineWidth',2},'transparent',1);
hold on;
plots(4)=shadedErrorBar([1:15],MeanMu',((StdMc)./sqrt(size(data_to_plot4,2)))*1.96,'lineprops',{'--r','markerfacecolor','r','LineWidth',2},'transparent',1);

xticks={'','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15',''};
xlim([0 16])
ylim([miny maxy])
if plott==1
    ylabel({'Hit rate (%)'})
elseif plott==2
    ylabel({'False alarm (%)'})
elseif plott==3
    ylabel({'Reaction time (ms)'})
end
set(gca,'FontSize',20,'LineWidth',3,'XTick',...
    [0:16],'XTickLabel',{'','','','','','','','','','','','','','','','',''},...
    'YTick',yticks,'YTickLabel',num2str(yticks'),'ycolor','k','xcolor','k');
box off;
legend ([plots(1).mainLine,plots(2).mainLine,plots(3).mainLine,plots(4).mainLine],{'Active cued','Monitoring cued','Active uncued','Monitoring uncued'},'location','southwest','edgecolor','none')
