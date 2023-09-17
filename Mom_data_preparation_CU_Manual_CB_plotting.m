%% Sample response processing
% developed by Hamid Karimi-Rouzbahani on 15/June/2022
% Modified by Hamid Karimi-Rouzbahani on 6/July/2022
% Modified by Hamid Karimi-Rouzbahani on 4/September/2023
% Modified by Hamid Karimi-Rouzbahani on 17/September/2023 to produce FA as
% well

clc
clear all;
addpath(genpath('F:\RESEARCH\Hamid\Multicentre dataset\Scripts\bayesFactor-master'))
subjects=[1:85 87 89:92]; % subjects you want the include in analysis
percentages=[0.5 0.09]; % Frequency of targets across conditions
chunks = [1:6]; % number of chunks
Testing_blocks=[1:5]; % blocks you want to include in analysis
for p=1:length(percentages)
    for i=1:5 % performance metrics including hit rate and RT
        Data{p,i}=nan(length(chunks)*length(Testing_blocks),max(subjects));
    end
end
Num_moving_dots=2;
Trials_per_block=32;

%% Data preparation
for Subj=subjects
    address=[['F:\RESEARCH\Hamid\Anina Zaid\OneDrive_2023-09-03\All data\REXP00',sprintf( '%02d', Subj )]];
    dirs=dir(address);
    for chunk=chunks
        for blk=Testing_blocks
            for cndss=percentages
                correct_reaction_times_att=0;
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
                        [Subj chunk blk]

                        chunkc=chunk;
                        if Targ_Freq_Condition_blk==cndss
                            % Hit rate
                            Data{cond,1}((chunkc-1)*length(Testing_blocks)+blk,Subj)=(tp_att)./(tp_att+fn_att);

                            % True negative rate
                            Data{cond,2}((chunkc-1)*length(Testing_blocks)+blk,Subj)=(tn_att)./(tn_att+fp_att);

                            % False alarm
                            Data{cond,3}((chunkc-1)*length(Testing_blocks)+blk,Subj)=(fp_att)./(fp_att+tn_att);

                            % Miss rate
                            Data{cond,4}((chunkc-1)*length(Testing_blocks)+blk,Subj)=(fn_att)./(tp_att+fn_att);

                            % Reaction time
                            Data{cond,5}((chunkc-1)*length(Testing_blocks)+blk,Subj)=correct_reaction_times_att;
                        else
                            Data{cond,1}((chunkc-1)*length(Testing_blocks)+blk,Subj)=nan;
                            Data{cond,2}((chunkc-1)*length(Testing_blocks)+blk,Subj)=nan;
                            Data{cond,3}((chunkc-1)*length(Testing_blocks)+blk,Subj)=nan;
                            Data{cond,4}((chunkc-1)*length(Testing_blocks)+blk,Subj)=nan;
                            Data{cond,5}((chunkc-1)*length(Testing_blocks)+blk,Subj)=nan;
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
        end
    end
    cond=1;
    Hit_rate_condition=Data{cond,1}(:,Subj); % Hit rate in condition
    Mean_Hit_rate=nanmean(Hit_rate_condition);
    FA_rate_condition=Data{cond,3}(:,Subj); % False alarm rate in condition
    Mean_FA_rate=nanmean(FA_rate_condition);
    Reaction_time_condition=Data{cond,5}(:,Subj); % reaction time in condition
    Mean_Reaction_time=nanmean(Reaction_time_condition);

    HR=[Hit_rate_condition;nan(5,1);Mean_Hit_rate];
    FA=[FA_rate_condition;nan(5,1);Mean_FA_rate];
    RT=[Reaction_time_condition;nan(5,1);Mean_Reaction_time];

    T = table(HR,FA,RT);
    T.Properties.VariableNames = {['Hit_rate_target_freq_',num2str(percentages(cond)*100)] ['FA_rate_target_freq_',num2str(percentages(cond)*100)] ['RT_target_freq_',num2str(percentages(cond)*100)]};
    Ttotal=T;
    Data_csv_total=[HR FA RT];

    for cond=2:size(Data,1)

        Hit_rate_condition=Data{cond,1}(:,Subj); % Hit rate in condition
        Mean_Hit_rate=nanmean(Hit_rate_condition);
        FA_rate_condition=Data{cond,3}(:,Subj); % FA rate in condition
        Mean_FA_rate=nanmean(FA_rate_condition);
        Reaction_time_condition=Data{cond,5}(:,Subj); % reaction time in condition
        Mean_Reaction_time=nanmean(Reaction_time_condition);
        HR=[Hit_rate_condition;nan(5,1);Mean_Hit_rate];
        FA=[FA_rate_condition;nan(5,1);Mean_FA_rate];
        RT=[Reaction_time_condition;nan(5,1);Mean_Reaction_time];
        T = table(HR,FA,RT);
        T.Properties.VariableNames = {['Hit_rate_target_freq_',num2str(percentages(cond)*100)] ['FA_rate_target_freq_',num2str(percentages(cond)*100)] ['RT_target_freq_',num2str(percentages(cond)*100)]};
        Ttotal=[Ttotal T];
        Data_csv_total=horzcat(Data_csv_total,[HR FA RT]);
    end

    filename = ['MoM_data_CU_manually_counter_balanced.xlsx']; % Change the name to anything you prefer
    writetable(Ttotal,filename,'Sheet',['Subj_' num2str(Subj)])

    [Subj]
end

%% Plotting the results

plott=3; % 1= Hit Rate; 3= False Alarm; 5=RT
gca = axes('Position',[0.22 0.25 0.775 0.72]);
if plott==1
    data_to_plot1=(Data{1,plott}(:,subjects))*100;
    data_to_plot2=(Data{2,plott}(:,subjects))*100;
    miny=20;
    maxy=100;
    yticks=[0:20:100];
elseif plott==3
    data_to_plot1=(Data{1,plott}(:,subjects))*100;
    data_to_plot2=(Data{2,plott}(:,subjects))*100;
    miny=0;
    maxy=40;
    yticks=[0:20:40];
elseif plott==5
    data_to_plot1=Data{1,5}(:,subjects)*1000;
    data_to_plot2=Data{2,5}(:,subjects)*1000;
    miny=200;
    maxy=400;
    yticks=[200:50:400];
end
MeanA=(nanmean([data_to_plot1(1:15,:) data_to_plot1(16:30,:)],2));
StdA=nanstd([data_to_plot1(1:15,:) data_to_plot1(16:30,:)]');
MeanM=(nanmean([data_to_plot2(1:15,:) data_to_plot2(16:30,:)],2));
StdM=nanstd([data_to_plot2(1:15,:) data_to_plot2(16:30,:)]');
plots(1)=shadedErrorBar([1:15],MeanA',((StdA)./sqrt(size(data_to_plot1,2)))*1.96,'lineprops',{'b','markerfacecolor','b','LineWidth',2},'transparent',1);
hold on;
plots(2)=shadedErrorBar([1:15],MeanM',((StdM)./sqrt(size(data_to_plot2,2)))*1.96,'lineprops',{'r','markerfacecolor','r','LineWidth',2},'transparent',1);

xticks={'','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15',''};
xlim([0 16])
ylim([miny maxy])
if plott==1
    ylabel({'Hit rate (%)'})
elseif plott==3
    ylabel({'False alarm (%)'})
elseif plott==5
    ylabel({'Reaction time (ms)'})
end
set(gca,'FontSize',20,'LineWidth',3,'XTick',...
    [0:16],'XTickLabel',{'','','','','','','','','','','','','','','','',''},...
    'YTick',yticks,'YTickLabel',num2str(yticks'),'ycolor','k','xcolor','k');
box off;
if plott==1
    legend ([plots(1).mainLine,plots(2).mainLine],{'Active','Monitoring'},'location','southwest','edgecolor','none')
end

% set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
pdf_paper_size         = [20 20];
fig.PaperSize       = pdf_paper_size;
pdf_print_resolution   = '-r300';
pdf_file_name=['Beh_data_',num2str(plott)];
% print(['\' pdf_file_name '.pdf'], '-dpdf', pdf_print_resolution)

% Bayes
for blk=1:15
    dt1=data_to_plot1(blk,:);
    dt1=dt1(~isnan(dt1));
    dt2=data_to_plot2(15+blk,:);
    dt2=dt2(~isnan(dt2));
    t(blk)=bf.ttest2(dt1,dt2);
    [ttest(blk),pp(blk)]=ttest2(dt1,dt2);
end

Effects=t';
p=1;
Threshold=6;
Neutral_color=[0.3 0.3 0.3];

maxy_sig=16;
miny_sig=-16;
figure;

gca = axes('Position',[0.22 0.77 0.775 0.15]);

colors={'k'};

line([0 16],[0 0],'Color','k','linewidth',1);
hold on;
for dots=1:size(Effects,1)
    if Effects(dots,p)>Threshold
        plot([dots],log10(Effects(dots,p)),'LineStyle','none','marker','o','MarkerFaceColor',colors{p},'Color',colors{p},'linewidth',1,'markersize',6);
    elseif Effects(dots,p)<Threshold && Effects(dots,p)>1./Threshold
        plot([dots],log10(Effects(dots,p)),'LineStyle','none','marker','o','MarkerFaceColor',[1 1 1],'Color',Neutral_color,'linewidth',1,'markersize',6);
    elseif Effects(dots,p)<1./Threshold
        plot([dots],log10(Effects(dots,p)),'LineStyle','none','marker','o','MarkerFaceColor',Neutral_color,'Color',Neutral_color,'linewidth',1,'markersize',6);
    end
    hold on;
end


set(gca,'FontSize',20,'LineWidth',3,'XTick',...
    [0:16],'XTickLabel',xticks,'YTick',...
    [(miny_sig) 0 (maxy_sig)],'YTickLabel',...
    {['{',num2str(miny_sig),'}'],'{0}',...
    ['{',num2str(maxy_sig),'}']},'YMinorTick','on',...
    'XTickLabelRotation',60,'ycolor','k','xcolor','k');
ylabel('BF (Log_{10})')
box off;
xlim ([0 16]);
ylim([miny_sig maxy_sig])
xlabel('Block #')
title ('Active vs. Monitoring')

%     % set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
pdf_paper_size         = [20 20];
fig.PaperSize       = pdf_paper_size;
pdf_print_resolution   = '-r300';
pdf_file_name=['Beh_data_Bayes_',num2str(plott)];
% print(['\' pdf_file_name '.pdf'], '-dpdf', pdf_print_resolution)
