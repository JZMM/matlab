% for probe session, include 4 groups, Low train, Low probe, High probe,
% High train

load('ZL_h37_Rig2P_20161113_Virables');

LimitTime = 5000;


Min_Onset_Time = Data_extract.Min_Onset_Time;
Left_lick_time = Data_extract.Left_lick_time(TestTrialNum);
Right_lick_time = Data_extract.Right_lick_time(TestTrialNum);
Tone_onset_time = double(Data_extract.Tone_onset_time(TestTrialNum));
Tone_frequency = Data_extract.Tone_frequency(TestTrialNum);
Action_choice = Data_extract.Action_choice(TestTrialNum);
Miss_Ind = Data_extract.Miss_Ind(TestTrialNum);
Frequencies = Data_extract.Frequencies;
Answer_time = Data_extract.Answer_time(TestTrialNum);
Probe_index = Data_extract.Probe_index(TestTrialNum);

Boundary = Frequencies(end)/2;
LowTrain_temp = Tone_frequency == Frequencies(1) & Action_choice == 0; 
HighTrain_temp = Tone_frequency == Frequencies(end) & Action_choice == 1; 
LowProbe_temp = Tone_frequency < Boundary & Probe_index & Action_choice == 0;
HighProbe_temp = Tone_frequency > Boundary & Probe_index & Action_choice == 1;
%%
LowTrain_LLick = [];
LowTrain_RLick = [];
LowProbe_LLick = [];
LowProbe_RLick = [];
HighProbe_LLick = [];
HighProbe_RLick = [];
HighTrain_LLick = [];
HighTrain_RLick = [];
for i = TestTrialNum
    if LowTrain_temp(i)
        LowTrain_LLick = [LowTrain_LLick Left_lick_time{i}-Tone_onset_time(i)];
        LowTrain_RLick = [LowTrain_RLick Right_lick_time{i}-Tone_onset_time(i)];
    end
    if LowProbe_temp(i)
        LowProbe_LLick = [LowProbe_LLick Left_lick_time{i}-Tone_onset_time(i)];
        LowProbe_RLick = [LowProbe_RLick Right_lick_time{i}-Tone_onset_time(i)];        
    end
    if HighProbe_temp(i)
        HighProbe_LLick = [HighProbe_LLick Left_lick_time{i}-Tone_onset_time(i)];
        HighProbe_RLick = [HighProbe_RLick Right_lick_time{i}-Tone_onset_time(i)];    
    end
    if HighTrain_temp(i)
        HighTrain_LLick = [HighTrain_LLick Left_lick_time{i}-Tone_onset_time(i)];
        HighTrain_RLick = [HighTrain_RLick Right_lick_time{i}-Tone_onset_time(i)];       
    end
end

LowTrain_LLick(LowTrain_LLick > LimitTime) = [];
LowTrain_RLick(LowTrain_RLick > LimitTime) = [];
LowProbe_LLick(LowProbe_LLick > LimitTime) = [];
LowProbe_RLick(LowProbe_RLick > LimitTime) = [];
HighProbe_LLick(HighProbe_LLick > LimitTime) = [];
HighProbe_RLick(HighProbe_RLick > LimitTime) = [];
HighTrain_LLick(HighTrain_LLick > LimitTime) = [];
HighTrain_RLick(HighTrain_RLick > LimitTime) = [];
%%
load('AlineCaSig_h36_20161127_1st_of_3','AlineCaSigData','FrameTime');
load('ZL_h36_Rig2P_20161127-fix_Virables');
saveName = 'h37_20161115';
is2P = 1;

titleName = {'Train Low','Probe Low','Probe High','Train High'};
Xlimit = AlineCaSigData.length_frames_T*FrameTime;

StarT = 0;
EndT = 7;

bin_size = 0.001;
kx=-1:bin_size:1;
krn=normpdf(kx,0,0.2);  % kernel is a half-gaussian with sd of 50 ms.
krn(kx<0)=[];
krn = krn/max(krn);
LengOfGau = length(krn);
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

Boundary = Frequencies(end)/2;
LowTrain_Ind = Tone_frequency == Frequencies(1) & Action_choice == 0; 
HighTrain_Ind = Tone_frequency == Frequencies(end) & Action_choice == 1; 
LowProbe_Ind = Tone_frequency < Boundary & Probe_index & Action_choice == 0;
HighProbe_Ind = Tone_frequency > Boundary & Probe_index & Action_choice == 1;
Inds = [LowTrain_Ind;LowProbe_Ind;HighProbe_Ind;HighTrain_Ind];


for tr = 1:length(Left_lick_time)
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

for L = 1:4
    LLicks(L,:) =  mean(L_Lick(Inds(L,:),:));
    Lstd(L,:) = std(L_Lick(Inds(L,:),:))/(size(L_Lick(Inds(L,:),:),1))^0.5;
    RLicks(L,:) =  mean(R_Lick(Inds(L,:),:));
    Rstd(L,:) = std(R_Lick(Inds(L,:),:))/(size(R_Lick(Inds(L,:),:),1))^0.5;
end
fig = figure;set(fig,'position',[2000 500 1000 300],'color','w');
MaxV = max([max(max(LLicks)) max(max(RLicks))]);
for j = 1:4
   subplot(1,4,j);hold on;axis square;
   if j <3
        plot(LLicks(j,:),'color','k','linewidth',1);
        patch([1:size(LLicks,2) fliplr(1:size(LLicks,2))],[LLicks(j,:)+Lstd(j,:) fliplr(LLicks(j,:)-Lstd(j,:))],'k','edgecolor','none')
   else
       plot(RLicks(j,:),'color','r','linewidth',1);
       patch([1:size(LLicks,2) fliplr(1:size(LLicks,2))],[RLicks(j,:)+Rstd(j,:) fliplr(RLicks(j,:)-Rstd(j,:))],'r','edgecolor','none')
   end
   alpha 0.3
   line([Min_Onset_Time Min_Onset_Time],[0 0.1],'color','b','linewidth',3);
   if is2P
       line([Min_Onset_Time+800 Min_Onset_Time+800],[0 0.1],'color','c','linewidth',3);
   else
       line([Min_Onset_Time+600 Min_Onset_Time+600],[0 0.1],'color','c','linewidth',3);
   end
   set(gca,'fontsize',15,'fontweight','bold','xtick',[Min_Onset_Time:2000:8000],'xticklabel',[0:2:8],'ytick',[0 1],'yticklabel',[0 1]);
   ylim([0 MaxV*1.1]);
   xlim([1 Xlimit]);
   if j ==1, xlabel('Time(s)'),ylabel('Lick Rate/s');end
   title(titleName{j});
end
t = suptitle(saveName);set(t,'fontsize',15,'fontweight','bold','interpreter','none');
fig = figure;set(fig,'position',[2000 500 1000 300],'color','w');
for j = 1:4
   subplot(1,4,j);hold on;axis square;
   if j <3
        plot(RLicks(j,:),'color','k','linewidth',1);
        patch([1:size(LLicks,2) fliplr(1:size(LLicks,2))],[RLicks(j,:)+Lstd(j,:) fliplr(RLicks(j,:)-Lstd(j,:))],'k','edgecolor','none')
   else
       plot(LLicks(j,:),'color','r','linewidth',1);
       patch([1:size(LLicks,2) fliplr(1:size(LLicks,2))],[LLicks(j,:)+Rstd(j,:) fliplr(LLicks(j,:)-Rstd(j,:))],'r','edgecolor','none')
   end
   alpha 0.3
   line([Min_Onset_Time Min_Onset_Time],[0 0.1],'color','b','linewidth',3);
   if is2P
       line([Min_Onset_Time+800 Min_Onset_Time+800],[0 0.1],'color','c','linewidth',3);
   else
       line([Min_Onset_Time+600 Min_Onset_Time+600],[0 0.1],'color','c','linewidth',3);
   end
   set(gca,'fontsize',15,'fontweight','bold','xtick',[Min_Onset_Time:2000:8000],'xticklabel',[0:2:8],'ytick',[0 1],'yticklabel',[0 1]);
   ylim([0 MaxV*1.1]);
   xlim([1 Xlimit]);
   if j ==1, xlabel('Time(s)'),ylabel('Lick Rate/s');end
   title(titleName{j});
end
t = suptitle(saveName);set(t,'fontsize',15,'fontweight','bold','interpreter','none');
%%
figure;hold on;
bar(tempLick);
figure;hold on;
for i =1:length(tempLick)
    Licks(i,:) = tempLick(i) + krn;
end

 plot(Licks(1,:));






