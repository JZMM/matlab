% if is opto, there are 6 groups, Low and High, Train low and Train probe,
% Train high and train probbe
% if is not opto, there are only 2 groups, Low and High

% the session is for statistis
%
%%  plot the miss
files = dir('*.mat');
for f = 14%:length(files)
    clearvars TestTrialNum
    fn = files(f).name;
    load(fn);
%     if exist('TestTrialNum')
%         continue;
%     end
    Miss_Ind = Data_extract.Miss_Ind;
    figure;
    plot(Miss_Ind,'ko','linestyle','none');  
    disp(['The filename is ' fn '. The maximum trial number is ' num2str(length(Miss_Ind))]);
    endN = input('Input the end trial number:');
    if endN == 0
        TestTrialNum = [1:length(Miss_Ind)];
    else
        TestTrialNum = [1:endN];
    end
    save(fn,'TestTrialNum','-append');
end
%%  this session is to plot trace seperate the correct and wrong
close all

savefig = 0;     % ==============================
isOpto = 0;      % ============================== in one session, has either opto and non-opto
StarT = 0;
EndT = 7;
Xlimit = EndT*1000;
titleName = {'Low','High','Low Opto','Low NonOpto','High Opto','High NonOpto'};

if isOpto FieldsName = titleName; else FieldsName = titleName(1:2);end
if savefig
    if ~exist('LickTrace_Seper')
        mkdir('LickTrace_Seper');
    end
    cd LickTrace_Seper;
    SaveP = pwd;
    cd ..;
end

bin_size = 0.001;
kx=-1:bin_size:1;
krn=normpdf(kx,0,0.2);  % kernel is a half-gaussian with sd of 50 ms.
krn(kx<0)=[];
krn = krn/max(krn);
LengOfGau = length(krn);

files = dir('*.mat');
for f = 1%:length(files)
    fn = files(f).name;
    load(fn);     
    if strcmp('ZL',fn(1:2))
        saveName = [fn(4:7) fn(14:21)];   
    else
        disp(fn);
        saveName = input('Input a saveName name:');
    end
    AnimName{f} = fn(4:6);
%
    Min_Onset_Time = double(Data_extract.Min_Onset_Time);
    Left_lick_time = Data_extract.Left_lick_time(TestTrialNum);
    Right_lick_time = Data_extract.Right_lick_time(TestTrialNum);
    Tone_onset_time = double(Data_extract.Tone_onset_time(TestTrialNum));
    Tone_frequency = Data_extract.Tone_frequency(TestTrialNum);
    Action_choice = Data_extract.Action_choice(TestTrialNum);
    Miss_Ind = Data_extract.Miss_Ind(TestTrialNum);
    Frequencies = Data_extract.Frequencies;
    Answer_time = Data_extract.Answer_time(TestTrialNum);
    Probe_index = Data_extract.Probe_index(TestTrialNum);
    Opto_trial_index = Data_extract.Opto_trial_index(TestTrialNum);

    Boundary = Frequencies(end)/2;
    
    Low_Ind = Tone_frequency < Boundary & Action_choice == 0; 
    High_Ind = Tone_frequency > Boundary & Action_choice == 1; 
    if isOpto
        LowOpto_Ind = Low_Ind & Opto_trial_index; 
        LowNonOpto_Ind = Low_Ind & ~Opto_trial_index; 
        HighOpto_Ind = High_Ind & Opto_trial_index; 
        HighNonOpto_Ind = High_Ind & ~Opto_trial_index;  
        Inds = [Low_Ind;High_Ind;LowOpto_Ind;LowNonOpto_Ind;HighOpto_Ind;HighNonOpto_Ind];
    else
        Inds = [Low_Ind;High_Ind];
    end
    for tr = 1:length(Left_lick_time)   % transform to gaussian
        L_Lick_temp = Left_lick_time{tr}-Tone_onset_time(tr)+Min_Onset_Time;
        R_Lick_temp = Right_lick_time{tr}-Tone_onset_time(tr)+Min_Onset_Time;
        L_Lick_temp(L_Lick_temp>EndT*1000) = [];
        R_Lick_temp(R_Lick_temp>EndT*1000) = [];
        L_Lick_temp(L_Lick_temp<=0) = [];
        R_Lick_temp(R_Lick_temp<=0) = [];
        L_Lick(tr,:) = zeros(1,length([StarT:bin_size:EndT+1]));   % initiate lick trace
        for i = 1:length(L_Lick_temp)
            L_Lick(tr,L_Lick_temp(i):1:L_Lick_temp(i)+LengOfGau-1) = L_Lick(tr,L_Lick_temp(i):1:L_Lick_temp(i)+LengOfGau-1) + krn;
        end
        R_Lick(tr,:) = zeros(1,length([StarT:bin_size:EndT+1]));   % initiate lick trace
        for i = 1:length(R_Lick_temp)
            R_Lick(tr,R_Lick_temp(i):1:R_Lick_temp(i)+LengOfGau-1) = R_Lick(tr,R_Lick_temp(i):1:R_Lick_temp(i)+LengOfGau-1) + krn;
        end    
    end
    for L = 1:size(Inds,1)
        LLicks(L,:) = mean(L_Lick(Inds(L,:),:));
        Lstd(L,:) = std(L_Lick(Inds(L,:),:))/(size(L_Lick(Inds(L,:),:),1))^0.5;
        RLicks(L,:) =  mean(R_Lick(Inds(L,:),:));
        Rstd(L,:) = std(R_Lick(Inds(L,:),:))/(size(R_Lick(Inds(L,:),:),1))^0.5;
    end
    for L = 1:size(Inds,1)
        if size(Inds,1) > 2
            if L == 1 | L == 3 | L == 4
                temp1 = LLicks(L,:);
                [~,IndMax] = max(RLicks(2,:));
                temp1(1:IndMax) = temp1(1:IndMax)-LLicks(2,1:IndMax);
                temp1(temp1<0) = 0;
                LLicks(L,:) = temp1;
            else
                temp2 = RLicks(L,:);
                [~,IndMax] = max(LLicks(1,:));
                temp2(1:IndMax) = temp2(1:IndMax)-RLicks(1,1:IndMax);                
                temp2(temp1<0) = 0;
                RLicks(L,:) = temp2;                               
            end
        else
            if L == 1 
                temp1 = LLicks(L,:);
                [~,IndMax] = max(RLicks(2,:));
                temp1(1:IndMax) = temp1(1:IndMax)-LLicks(2,1:IndMax);
                temp1(temp1<0) = 0;
                LLicks(L,:) = temp1;
            else
                temp2 = RLicks(L,:);
                [~,IndMax] = max(LLicks(1,:));
                temp2(1:IndMax) = temp2(1:IndMax)-RLicks(1,1:IndMax);  
                temp2(temp1<0) = 0;
                RLicks(L,:) = temp2;                               
            end                            
        end       
    end
    fig1 = figure;set(fig1,'position',[2000 500 1000 300],'color','w');
    MaxV = max([max(max(LLicks)) max(max(RLicks))]);
    pks = zeros(size(Inds,1),2);
    locs = zeros(size(Inds,1),2);
    w = zeros(size(Inds,1),2);
    for j = 1:size(Inds,1)
        subplot(1,size(Inds,1),j);hold on;axis square;        
        findpeaks(LLicks(j,:),'MinPeakProminence',1,'Annotate','extents');
        [pks,locs,w] = findpeaks(LLicks(j,:),'MinPeakProminence',0.5,'Annotate','extents');
        if isempty(pks)
            temp1(j,1) = NaN; PeaksInd(j,1) = NaN; Widt(j,1) = NaN;
            HafInd(j,1) = NaN;
            OneRisingTime(j,1) = NaN;
            TwoRisingTime(j,1) = NaN;
        else
            temp1(j,1) = pks; PeaksInd(j,1) = locs; Widt(j,1) = w; 
            [~,HafInd(j,1)] = min(abs(LLicks(j,1:locs)-pks/2));
            [~,OneRisingTime(j,1)] = min(abs(LLicks(j,1:locs)-std(LLicks(j,:))));
            [~,TwoRisingTime(j,1)] = min(abs(LLicks(j,1:locs)-2*std(LLicks(j,:))));
        end
        line([0 2000],[2*std(LLicks(j,:)) 2*std(LLicks(j,:))]);
        line([0 2000],[std(LLicks(j,:)) std(LLicks(j,:))]);
        legend('off');
        findpeaks(RLicks(j,:),'MinPeakProminence',1,'Annotate','extents');
        [pks,locs,w] = findpeaks(RLicks(j,:),'MinPeakProminence',0.7,'Annotate','extents');
        if isempty(pks)
            temp1(j,2) = NaN; PeaksInd(j,2) = NaN; Widt(j,2) = NaN;
            HafInd(j,2) = NaN;
            OneRisingTime(j,2) = NaN;
            TwoRisingTime(j,2) = NaN;
        else
            temp1(j,2) = pks; PeaksInd(j,2) = locs; Widt(j,2) = w; 
            [~,HafInd(j,2)] = min(abs(RLicks(j,1:locs)-pks/2));
            [~,OneRisingTime(j,2)] = min(abs(RLicks(j,1:locs)-std(RLicks(j,:))));
            [~,TwoRisingTime(j,2)] = min(abs(RLicks(j,1:locs)-2*std(RLicks(j,:))));
        end
        line([0 2000],[2*std(RLicks(j,:)) 2*std(RLicks(j,:))]);
        line([0 2000],[std(RLicks(j,:)) std(RLicks(j,:))]);
        legend('off');        
%         line([ind1(j) ind1(j)],[0 2],'color','k');
%         line([ind2(j) ind2(j)],[0 2],'color','r');
%         line(LefHafInd(j,:),LefHafPeak(j,:)+Max1(j)/2,'color','k');
%         line(RighHafInd(j,:),RighHafPeak(j,:)+Max2(j)/2,'color','r');
        
%         plot(LLicks(j,:),'color','k','linewidth',1);
        patch([1:size(LLicks,2) fliplr(1:size(LLicks,2))],[LLicks(j,:)+Lstd(j,:) fliplr(LLicks(j,:)-Lstd(j,:))],'k','edgecolor','none')
%         plot(RLicks(j,:),'color','r','linewidth',1);
        patch([1:size(LLicks,2) fliplr(1:size(LLicks,2))],[RLicks(j,:)+Rstd(j,:) fliplr(RLicks(j,:)-Rstd(j,:))],'r','edgecolor','none')

       alpha 0.3
       line([Min_Onset_Time Min_Onset_Time],[0 0.1],'color','b','linewidth',3);
       line([Min_Onset_Time+600 Min_Onset_Time+600],[0 0.1],'color','c','linewidth',3);

       set(gca,'fontsize',15,'fontweight','bold','xtick',[Min_Onset_Time:1000:8000],'xticklabel',[0:1:8],'ytick',[0 1],'yticklabel',[0 1]);
       ylim([0 MaxV*1.1]);
       xlim([1 Xlimit]);
       if j ==1, xlabel('Time(s)'),ylabel('Lick Counts');end
       title(titleName{j});
    end
    t = suptitle(saveName);set(t,'fontsize',15,'fontweight','bold','interpreter','none');
    if savefig
        saveas(fig1,[SaveP '/' saveName '_LickTrace_corr.fig']);
        saveas(fig1,[SaveP '/' saveName '_LickTrace_corr.png']);
    end
    PeakValue{f} = temp1;
    PeakInd{f} = PeaksInd-Min_Onset_Time;
    HafPeakWidt{f} = Widt;
    HafPeakTime{f} = HafInd-Min_Onset_Time;
    RiseTime_1std{f} = OneRisingTime-Min_Onset_Time;
    RiseTime_2std{f} = TwoRisingTime-Min_Onset_Time;
end
%%
save('VgatChR-1_1stOpto_Lick_Stat','AnimName','PeakValue','PeakInd','HafPeakWidt','HafPeakTime','RiseTime_1std','RiseTime_2std','FieldsName');
%%  plotting  for opto
load('VgatChR-1_1stOpto_Lick_Stat');
testData1 = HafPeakWidt; % PeakValue PeakInd RiseTime_1std RiseTime_2std HafPeakWidt HafPeakTime
load('VgatChR-2_3rdOpto_Lick_Stat');
testData2 = HafPeakWidt; % PeakValue PeakInd RiseTime_1std RiseTime_2std HafPeakWidt HafPeakTime
for i = 1:length(testData1)
    temp1(i,:) = [testData1{i}(3:4,1)' testData1{i}(5:6,2)'];
end
for i = 1:length(testData2)
    temp2(i,:) = [testData2{i}(3:4,1)' testData2{i}(5:6,2)'];
end
fig1 = figure;hold on;set(fig1,'position',[2000 200 400 300],'color','w');
axis square;
set(gca,'fontsize',15,'fontweight','bold','xtick',[1.5 3.5],'xticklabel',...
    {'Low Freqs','High Freqs'},'xticklabelrotation',0);


xlim([0.5 4.5]);
% first plot
patch([1.7 2.3 2.3 1.7],[0 0 5000 5000],'c','edgecolor','none');
patch([3.7 4.3 4.3 3.7],[0 0 5000 5000],'c','edgecolor','none');
plot(temp1(:,1:2)','color','k','linewidth',3);
% plot(ones(length(temp1(:,1)),1),temp1(:,1),'ko');
% plot(2*ones(length(temp1(:,1)),1),temp1(:,2),'bo');

plot(repmat([3;4],1,length(temp1(:,1))),temp1(:,3:4)','color','k','linewidth',3);
% plot(3*ones(length(temp1(:,1)),1),temp1(:,3),'ko');
% plot(4*ones(length(temp1(:,1)),1),temp1(:,4),'bo');

% second plot
plot(temp2(:,1:2)','color',[.7 .7 .7],'linewidth',3);
% plot(ones(length(temp2(:,1)),1),temp2(:,1),'marker','o','color',[.7 .7 .7],'linestyle','none');
% plot(2*ones(length(temp2(:,1)),1),temp2(:,2),'co');

plot(repmat([3;4],1,length(temp2(:,1))),temp2(:,3:4)','color',[.7 .7 .7],'linewidth',3);
% plot(3*ones(length(temp2(:,1)),1),temp2(:,3),'marker','o','color',[.7 .7 .7],'linestyle','none');
% plot(4*ones(length(temp2(:,1)),1),temp2(:,4),'co');
alpha 0.5

%%
ylim([200 2200]);
title('Lick Haft Height Width(VgtChR)');
ylabel('Time/ms');
%%   t-test
[h1,p1] = ttest([temp1(:,1);temp2(:,1)],[temp1(:,2);temp2(:,2)])
[h2,p2] = ttest([temp1(:,3);temp2(:,3)],[temp1(:,4);temp2(:,4)])
% [h,p] = ttest(temp1(:,1),temp1(:,2))
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  plotting  for between sessions
load('VgatChrim_2ndNonOpto_Lick_Stat');
testData1 = PeakValue; % PeakValue PeakInd RiseTime_1std RiseTime_2std HafPeakWidt HafPeakTime
load('VgatChrim_1stOpto_Lick_Stat');
testData2 = PeakValue; % PeakValue PeakInd RiseTime_1std RiseTime_2std HafPeakWidt HafPeakTime
for i = 1:length(testData1)
    temp1(i,:) = [testData1{i}(1,1) testData1{i}(2,2)];
end
for i = 1:length(testData2)
    temp2(i,:) = [testData2{i}(1,1) testData2{i}(2,2)];
end
fig1 = figure;hold on;set(fig1,'position',[2000 200 400 300],'color','w');
axis square;
set(gca,'fontsize',15,'fontweight','bold','xtick',[1.5 3.5],'xticklabel',...
    {'Low Freqs','High Freqs'},'xticklabelrotation',0);


xlim([0.5 4.5]);
% first plot
patch([1.7 2.3 2.3 1.7],[0 0 5000 5000],'c','edgecolor','none');
patch([3.7 4.3 4.3 3.7],[0 0 5000 5000],'c','edgecolor','none');
plot([temp1(:,1) temp2(:,1)]','color','k','linewidth',3);
% plot([temp1(:,1) temp2(:,1)]','color',[.7 .7 .7],'linewidth',3);
% plot(ones(length(temp1(:,1)),1),temp1(:,1),'ko');
% plot(2*ones(length(temp1(:,1)),1),temp1(:,2),'bo');

plot(repmat([3;4],1,length(temp1(:,1))),[temp1(:,2) temp2(:,2)]','color','k','linewidth',3);
% plot(repmat([3;4],1,length(temp1(:,1))),[temp1(:,2) temp2(:,2)]','color',[.7 .7 .7],'linewidth',3);
% plot(3*ones(length(temp1(:,1)),1),temp1(:,3),'ko');
% plot(4*ones(length(temp1(:,1)),1),temp1(:,4),'bo');

alpha 0.5

%%
ylim([0.5 2.5]);
title('Lick Haft Height Width(VgtChR)');
ylabel('Time/ms');
%%   t-test
[h1,p1] = ttest(temp1(:,1),temp2(:,1))
[h2,p2] = ttest(temp1(:,2),temp2(:,2))
% [h,p] = ttest(temp1(:,1),temp1(:,2))
%%
[h1,p1] = ttest([temp1(:,1);data1(:,1)],[temp2(:,1);data2(:,1)])
[h2,p2] = ttest([temp1(:,2);data1(:,2)],[temp2(:,2);data2(:,2)])





