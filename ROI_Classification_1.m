AligData = 'AlineCaSig_h35_20161127';
BehData ='ZL_h35_Rig2P_20161127_Virables';
DffFile = 'dff_ZL_h35_20161127_fied1_d150_3x';
load(BehData)
load(DffFile,'F_num_T','F_num_A','MinActFrame','MinActFrame');

isPlot = 1;
SaveFig = 0; 
SavePPT = 0;
SaveVar = 0;
SortWithReward = 1;
sigAlpha = 0.05;
load(AligData);
%%
AligTonData = AlineCaSigData.CaSigAlineTone;
AligActData = AlineCaSigData.CaSigAlineAct;
Tone_frequency = Data_extract.Tone_frequency(TestTrialNum);
Action_choice = Data_extract.Action_choice(TestTrialNum);
Probe_Ind = Data_extract.Probe_index(TestTrialNum);
RewardTime = Data_extract.Reward_time(TestTrialNum);
Freq = unique(Tone_frequency);
Boundary = Freq(end)/2;
LowTrain_temp = Tone_frequency == Freq(1) & Action_choice == 0; 
HighTrain_temp = Tone_frequency == Freq(end) & Action_choice == 1; 
LowProbe_temp = Tone_frequency < Boundary & Probe_Ind & Action_choice == 0;
HighProbe_temp = Tone_frequency > Boundary & Probe_Ind & Action_choice == 1;

LowTrain_temp_Wro = Tone_frequency == Freq(1) & Action_choice == 1; 
HighTrain_temp_Wro = Tone_frequency == Freq(end) & Action_choice == 0; 
LowProbe_temp_Wro = Tone_frequency < Boundary & Probe_Ind & Action_choice == 1;
HighProbe_temp_Wro = Tone_frequency > Boundary & Probe_Ind & Action_choice == 0;

RewTime = F_num_A(TestTrialNum,1)-F_num_T(TestTrialNum,1)+MinActFrame;
LowTrain = LowTrain_temp(TestTrialNum) ; 
HighTrain = HighTrain_temp(TestTrialNum) ; 
LowProbe = LowProbe_temp(TestTrialNum) ;
HighProbe = HighProbe_temp(TestTrialNum); 
LowTrain_Wro = LowTrain_temp_Wro(TestTrialNum) ; 
HighTrain_Wro = HighTrain_temp_Wro(TestTrialNum) ; 
LowProbe_Wro = LowProbe_temp_Wro(TestTrialNum) ;
HighProbe_Wro = HighProbe_temp_Wro(TestTrialNum); 
AllTypes = {LowTrain,LowProbe,HighProbe,HighTrain};
AllTypes = {LowTrain_Wro,LowProbe_Wro,HighProbe_Wro,HighTrain_Wro};
FrameLength = AlineCaSigData.length_frames_T;
TitleName = {'Train Low' 'Probe Low' 'Probe High' 'Train High'};
TitleName_Wro = {'Train Low Wro' 'Probe Low Wro' 'Probe High Wro' 'Train High Wro'};
FormerHaft = TestTrialNum(1:round(length(TestTrialNum)/2)) ; 
LaterHaft = TestTrialNum(round(length(TestTrialNum)/2)+1:end); 
%% First Lick Time
FirstLick = FirstLickTime_fun(Data_extract);
FirstLick_test = FirstLick(TestTrialNum);
FirstLick_test_frame = round(FirstLick_test/FrameTime) + double(MinOnsFrame);
FirstLick_test_b1 = FirstLick_test(FormerHaft);
FirstLick_test_b2 = FirstLick_test(LaterHaft);
for kk = 1:4
    WholeFirstLickTime{kk} = FirstLick_test(AllTypes{kk});
    WholeFirstLickTime_mean(kk) = mean(WholeFirstLickTime{kk});
    WholeFirstLickTime_SEM(kk) = std(WholeFirstLickTime{kk})/(length(WholeFirstLickTime{kk}))^0.5;
    WholeFirstLickTime_FirstHaf{kk} = FirstLick_test_b1(AllTypes{kk}(FormerHaft));
    WholeFirstLickTime_FirstHaf_mean(kk) = mean(WholeFirstLickTime_FirstHaf{kk});
    WholeFirstLickTime_FirstHaf_SEM(kk) = std(WholeFirstLickTime_FirstHaf{kk})/(length(WholeFirstLickTime_FirstHaf{kk}))^0.5;    
    WholeFirstLickTime_SecondHaf{kk} = FirstLick_test_b2(AllTypes{kk}(LaterHaft));
    WholeFirstLickTime_SecondHaf_mean(kk) = mean(WholeFirstLickTime_SecondHaf{kk});
    WholeFirstLickTime_SecondHaf_SEM(kk) = std(WholeFirstLickTime_SecondHaf{kk})/(length(WholeFirstLickTime_SecondHaf{kk}))^0.5;      
end
figure;hold on; set(gcf,'position',[2900 290 250 150]);
line([0 5],[800 800],'color',[.7 .7 .7],'linestyle','--');
errorbar(1:4,WholeFirstLickTime_FirstHaf_mean,WholeFirstLickTime_FirstHaf_SEM,'b.','markersize',15,'linestyle','-');
errorbar(1:4,WholeFirstLickTime_SecondHaf_mean,WholeFirstLickTime_SecondHaf_SEM,'r.','markersize',15,'linestyle','-');
errorbar(1:4,WholeFirstLickTime_mean,WholeFirstLickTime_SEM,'k.','markersize',15,'linestyle','-');
xlim([0.5 4.5]);
ylim([0 1500]);
title('First Lick Time');
%%
ColLimit = 80;

for ROI  = 32%1:size(testData,2)
    close all;    

    SelectData_Tone = squeeze(AligTonData(:,ROI,:))*100;
    SelectData_Act = squeeze(AligActData(:,ROI,:))*100;
    testData1 = SelectData_Tone(FormerHaft,:);
    testData2 = SelectData_Tone(LaterHaft,:);
%     ClimLow = prctile(reshape(SelectData_Tone,1,[]),30);
    ClimLow = 0;
    ClimHigh = prctile(reshape(SelectData_Tone,1,[]),ColLimit);
    ClimHighTrace = prctile(reshape(SelectData_Tone,1,[]),90);
%     if isPlot
    fig = figure; set(fig,'position',[1900 750 1000 350],'color','w');
%     end
    for i =1:4
        tempData_Tone = SelectData_Tone(AllTypes{i},:);
        tempData_Act = SelectData_Act(AllTypes{i},:);
%         [ToneSig(ROI,i),~] = ttest(median(tempData(:,1:MinOnsFrame),2),median(tempData(:,MinOnsFrame:MinOnsFrame+round(1000/FrameTime)),2),'alpha',0.01);
        [~,ToneSig(ROI,i)] = ranksum(mean(tempData_Tone(:,1:MinOnsFrame),2),mean(tempData_Tone(:,MinOnsFrame:MinOnsFrame+round(1000/FrameTime)),2),'alpha',sigAlpha);
        ToneOverBase(ROI,i) = mean(mean(tempData_Tone(:,MinOnsFrame:MinOnsFrame+round(800/FrameTime)),2)-mean(tempData_Tone(:,1:MinOnsFrame),2));
%         [ChoicSig(ROI,i),~] = ttest(median(tempData(:,1:MinOnsFrame),2),median(DataAlig2Act(:,MinActFrame-round(800/FrameTime):MinActFrame),2),'alpha',0.01);
        [~,ChoicSig(ROI,i),~] = ranksum(mean(tempData_Tone(:,1:MinOnsFrame),2),mean(tempData_Act(:,MinActFrame-round(800/FrameTime):MinActFrame),2),'alpha',sigAlpha);
        ChoiceOverBase(ROI,i) = mean(mean(tempData_Act(:,MinActFrame-round(800/FrameTime):MinActFrame),2)-mean(tempData_Tone(:,1:MinOnsFrame),2));
        [~,RewardSig(ROI,i),~] = ranksum(mean(tempData_Act(:,MinActFrame-round(800/FrameTime):MinActFrame),2),mean(tempData_Act(:,MinActFrame:MinActFrame+45),2),'alpha',sigAlpha);
        RewardOverAct(ROI,i) = mean(mean(tempData_Act(:,MinActFrame:MinActFrame+45),2)-mean(tempData_Act(:,MinActFrame-round(800/FrameTime):MinActFrame),2));
%         [LateVsBaseSig(ROI,i),~] = ttest(median(DataAlig2Act(:,MinActFrame+45:end),2),median(tempData(:,1:MinOnsFrame),2),'alpha',0.01);
        [~,LateVsBaseSig(ROI,i),~] = ranksum(mean(tempData_Tone(:,1:MinOnsFrame),2),mean(tempData_Act(:,MinActFrame+45:end),2),'alpha',sigAlpha);
        LateOverBase(ROI,i) = mean(mean(tempData_Act(:,MinActFrame+45:end),2)-mean(tempData_Tone(:,1:MinOnsFrame),2));
        [~,LateVsRewardSig(ROI,i),~] = ranksum(mean(tempData_Act(:,MinActFrame:MinActFrame+45),2),mean(tempData_Act(:,MinActFrame+50:end),2),'alpha',sigAlpha);
        LateOverReward(ROI,i) = mean(mean(tempData_Act(:,MinActFrame+30:end),2))-mean(mean(tempData_Act(:,MinActFrame:MinActFrame+15),2));
        subplot(1,4,i);hold on;
        [RewardTime Inds] = sort(RewTime(AllTypes{i}));
        ChoiceTime = FirstLick_test_frame(AllTypes{i});
        imagesc(tempData_Tone(Inds,:));
        plot(RewardTime,[1:sum(AllTypes{i})],'r.','markersize',8);
        plot(ChoiceTime(Inds),[1:sum(AllTypes{i})],'marker','.','color','c','markersize',8,'linestyle','none');
        set(gca,'clim',[ClimLow ClimHigh],'fontsize',15,'fontweight','bold','xtick',[double(MinOnsFrame):1000/FrameTime:FrameLength],'xticklabel',[0:1:10]);
        ylim([1 sum(AllTypes{i})]);
        xlim([0 FrameLength]);
        line([MinOnsFrame MinOnsFrame],[0 sum(AllTypes{i})],'color','w','linewidth',2);
        if i == 1
           ylabel('Number of Trials');
           xlabel('Time(s)')
    %        text(-20,sum(AllTypes{i})+5,'Stim Onset');
        end
        if i ==4
           ColBar = colorbar('position',[0.92 0.1 0.03 0.5]); 
           text(260,sum(AllTypes{i}*0.7),'\DeltaF/F0(%)','fontweight','bold');
        end
        title(TitleName{i});
    end
    fig = figure; set(fig,'position',[1900 250 1000 200],'color','w');  % Wrong trial color plot
    for i = 1:4
        tempDataWro_Tone = SelectData_Tone(AllTypes_Wro{i},:);
        subplot(1,4,i);hold on;
        [RewardTime Inds] = sort(RewTime(AllTypes_Wro{i}));
        ChoiceTime = FirstLick_test_frame(AllTypes_Wro{i});
        imagesc(tempDataWro_Tone(Inds,:));
        plot(RewardTime,[1:sum(AllTypes_Wro{i})],'r.','markersize',8);
        plot(ChoiceTime(Inds),[1:sum(AllTypes_Wro{i})],'marker','.','color','c','markersize',8,'linestyle','none');
        set(gca,'clim',[ClimLow ClimHigh],'fontsize',15,'fontweight','bold','xtick',[double(MinOnsFrame):1000/FrameTime:FrameLength],'xticklabel',[0:1:10]);
        ylim([1 sum(AllTypes_Wro{i})]);
        xlim([0 FrameLength]);
        line([MinOnsFrame MinOnsFrame],[0 sum(AllTypes_Wro{i})],'color','w','linewidth',2);
        if i == 1
           ylabel('Number of Trials');
           xlabel('Time(s)')
    %        text(-20,sum(AllTypes{i})+5,'Stim Onset');
        end
        if i ==4
           ColBar = colorbar('position',[0.92 0.1 0.03 0.5]); 
           text(260,sum(AllTypes_Wro{i}*0.7),'\DeltaF/F0(%)','fontweight','bold');
        end
        title(TitleName_Wro{i});        
    end
    fig = figure; set(fig,'position',[2900 750 1000 350],'color','w');   % Unsort Correct trials color plot
    for i =1:4
        subplot(1,4,i);hold on;
        tempData_Tone = SelectData_Tone(AllTypes{i},:);
        ChoiceTime = FirstLick_test_frame(AllTypes{i});
        imagesc(tempData_Tone);
        plot(RewTime(AllTypes{i}),[1:sum(AllTypes{i})],'r.','markersize',8);
        plot(ChoiceTime,[1:sum(AllTypes{i})],'c.','markersize',8);
        set(gca,'clim',[ClimLow ClimHigh],'fontsize',15,'fontweight','bold','xtick',[double(MinOnsFrame):1000/FrameTime:FrameLength],'xticklabel',[0:1:10]);
        ylim([1 sum(AllTypes{i})]);
        xlim([0 FrameLength]);
        line([MinOnsFrame MinOnsFrame],[0 sum(AllTypes{i})],'color','w','linewidth',2);        
        if i == 1
           ylabel('Number of Trials');
           xlabel('Time(s)')
        end
        if i ==4
           ColBar = colorbar('position',[0.92 0.1 0.03 0.5]); 
           text(260,sum(AllTypes{i}*0.7),'\DeltaF/F0(%)','fontweight','bold');
        end
        title(TitleName{i});        
    end
% 
%         figTab = figure; % create a table
%         set(figTab,'position',[2000 200 950 100]);
%         Tone_Sig=ToneSig(ROI,:); ToneSubBase=ToneOverBase(ROI,:); Choic_Sig=ChoicSig(ROI,:); ChoiceSubBase = ChoiceOverBase(ROI,:);
%         Reward_Sig=RewardSig(ROI,:); RewardSubAct=RewardOverAct(ROI,:); LatVsBas_Sig = LateVsBaseSig(ROI,:); LatVsRew_Sig=LateVsRewardSig(ROI,:);
%         LateSubReward=LateOverReward(ROI,:); LateSubBase=LateOverBase(ROI,:);
%         T = table(Tone_Sig',ToneSubBase',Choic_Sig',ChoiceSubBase',Reward_Sig',...
%             RewardSubAct',LatVsBas_Sig',LateSubBase',LatVsRew_Sig',LateSubReward',...
%             'VariableNames',{'Tone_Sig' 'ToneSubBase' 'Choic_Sig' 'ChoiceSubBase','Reward_Sig','RewardSubAct','LatVsBas_Sig','LateSubBase','LatVsRew_Sig','LateSubReward'},...
%             'RowNames',TitleName');
%         t_h = uitable(figTab,'Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
%         'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
%         t_h.Position(3) = 1;
%         t_h.Position(4) = 1;

        
        fig = figure; % plot Align to tone trace
        set(fig,'position',[1900 470 1000 220],'color','w');
        for i =1:4
            subplot(1,4,i);hold on;
%             axis square;
            temp_T = SelectData_Tone(AllTypes{i},:);
            temp_A = SelectData_Act(AllTypes{i},:);
            SEM_T = std(temp_T)/(sum(AllTypes{i}))^0.5;
            SEM_A = std(temp_A)/(sum(AllTypes{i}))^0.5;
            Act_Trace = mean(temp_A);
            line([1 FrameLength],[0 0],'linestyle','--','color',[.7 .7 .7],'linewidth',1);
            plot(mean(temp_T),'color','k');
            plot([1-MinActFrame+MinOnsFrame+round(800/FrameTime):1-MinActFrame+MinOnsFrame+round(800/FrameTime)+length(Act_Trace)-1],Act_Trace,'color','b');
            patch([1:FrameLength fliplr(1:FrameLength)],[mean(temp_T)-SEM_T fliplr(mean(temp_T)+SEM_T)],'k','edgecolor','none');
            patch([1:length(Act_Trace) fliplr(1:length(Act_Trace))],[Act_Trace-SEM_A fliplr(Act_Trace+SEM_A)],'b','edgecolor','none');
            line([MinOnsFrame MinOnsFrame],[ClimHighTrace*-0.2 ClimHighTrace*0.2],'color',[.7 .7 .7],'linewidth',1);
            line([MinOnsFrame+round(300/FrameTime)+1 round(MinOnsFrame+300/FrameTime)+1],[ClimHighTrace*-0.2 ClimHighTrace*0.2],'color','r','linewidth',1);
            line([MinOnsFrame+round(800/FrameTime)+1 round(MinOnsFrame+800/FrameTime)+1],[ClimHighTrace*-0.2 ClimHighTrace*0.2],'color','g','linewidth',1);            
            set(gca,'fontsize',15,'fontweight','bold','xtick',[double(MinOnsFrame):1000/FrameTime:FrameLength],'xticklabel',[0:1:10]);
            ylim([-50 ClimHighTrace]);
            xlim([0 FrameLength]);
            alpha 0.3

            if i == 1
               ylabel(['\DeltaF/F0']);
               xlabel('Time(s)')
%                plot(MinOnsFrame+15,ClimHigh-20,'marker','<','color','k','markersize',5,'markerfacecolor','k');
%                text(double(MinOnsFrame)+25,ClimHigh-18,'Stim Onset','fontweight','bold');
            end
            if i > 1
               set(gca,'ycolor','w','ytick',[],'yticklabel',[]); 
            end
            title(TitleName{i});
        end
        fig = figure; % plot Align to tone trace for Wrong trials
        set(fig,'position',[1900 -20 1000 220],'color','w');
        for i =1:4
            subplot(1,4,i);hold on;
            temp_T = SelectData_Tone(AllTypes_Wro{i},:);
            temp_A = SelectData_Act(AllTypes_Wro{i},:);
            SEM_T = std(temp_T)/(sum(AllTypes_Wro{i}))^0.5;
            SEM_A = std(temp_A)/(sum(AllTypes_Wro{i}))^0.5;
            Act_Trace = mean(temp_A);
            line([1 FrameLength],[0 0],'linestyle','--','color',[.7 .7 .7],'linewidth',1);
            plot(mean(temp_T),'color','k');
            plot([1-MinActFrame+MinOnsFrame+round(800/FrameTime):1-MinActFrame+MinOnsFrame+round(800/FrameTime)+length(Act_Trace)-1],Act_Trace,'color','b');
            patch([1:FrameLength fliplr(1:FrameLength)],[mean(temp_T)-SEM_T fliplr(mean(temp_T)+SEM_T)],'k','edgecolor','none');
            patch([1:length(Act_Trace) fliplr(1:length(Act_Trace))],[Act_Trace-SEM_A fliplr(Act_Trace+SEM_A)],'b','edgecolor','none');
            line([MinOnsFrame MinOnsFrame],[ClimHighTrace*-0.2 ClimHighTrace*0.2],'color',[.7 .7 .7],'linewidth',1);
            line([MinOnsFrame+round(300/FrameTime)+1 round(MinOnsFrame+300/FrameTime)+1],[ClimHighTrace*-0.2 ClimHighTrace*0.2],'color','r','linewidth',1);
            line([MinOnsFrame+round(800/FrameTime)+1 round(MinOnsFrame+800/FrameTime)+1],[ClimHighTrace*-0.2 ClimHighTrace*0.2],'color','g','linewidth',1);            
            set(gca,'fontsize',15,'fontweight','bold','xtick',[double(MinOnsFrame):1000/FrameTime:FrameLength],'xticklabel',[0:1:10]);
            ylim([-50 ClimHighTrace]);
            xlim([0 FrameLength]);
            alpha 0.3
            if i == 1
               ylabel(['\DeltaF/F0']);
               xlabel('Time(s)')
            end
            if i > 1
               set(gca,'ycolor','w','ytick',[],'yticklabel',[]); 
            end
            title(TitleName{i});
        end    
        
        fig = figure; hold on; % plot Align to Reward trace
        set(fig,'position',[2900 470 1000 220],'color','w');
        for i =1:4
            subplot(1,4,i);hold on;
            axis square;
            line([MinOnsFrame MinOnsFrame],[ClimHighTrace*-0.2 ClimHighTrace*0.2],'color',[.7 .7 .7],'linewidth',1);
            line([MinOnsFrame+round(300/FrameTime)+1 round(MinOnsFrame+300/FrameTime)+1],[ClimHighTrace*-0.2 ClimHighTrace*0.2],'color','r','linewidth',1);
            line([MinOnsFrame+round(800/FrameTime)+1 round(MinOnsFrame+800/FrameTime)+1],[ClimHighTrace*-0.2 ClimHighTrace*0.2],'color','g','linewidth',1);            
            line([1 FrameLength],[0 0],'linestyle','--','color',[.7 .7 .7],'linewidth',1);
            temp1 = testData1(AllTypes{i}(FormerHaft),:);
            temp2 = testData2(AllTypes{i}(LaterHaft),:);
            SEM1 = std(temp1)/(sum(AllTypes{i}(FormerHaft)))^0.5;
            SEM2 = std(temp2)/(sum(AllTypes{i}(LaterHaft)))^0.5;
            plot(mean(temp1),'color','k');
            plot(mean(temp2),'color','r');

            patch([1:FrameLength fliplr(1:FrameLength)],[mean(temp1)-SEM1 fliplr(mean(temp1)+SEM1)],'k','edgecolor','none');
            patch([1:FrameLength fliplr(1:FrameLength)],[mean(temp2)-SEM2 fliplr(mean(temp2)+SEM2)],'r','edgecolor','none');

            set(gca,'fontsize',15,'fontweight','bold','xtick',[double(MinOnsFrame):1000/FrameTime:FrameLength],'xticklabel',[0:1:10]);
            ylim([-50 ClimHighTrace]);
            xlim([0 FrameLength]);
            if i == 1
               ylabel(['\DeltaF/F0']);
               xlabel('Time(s)');
               plot(MinOnsFrame+15,ClimHigh-20,'marker','<','color','k','markersize',5,'markerfacecolor','k');
               text(double(MinOnsFrame)+25,ClimHigh-18,'Stim Onset','fontweight','bold');
            end
            if i > 1
               set(gca,'ycolor','w','ytick',[],'yticklabel',[]);
            end
            title(TitleName{i});
            alpha 0.3
        end   
end






















