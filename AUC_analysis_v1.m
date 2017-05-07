%% load a AUC of ROC data

fn1 = 'h37_20161113_2nd_of_3_AUC';
fn2 = 'h37_20161114_3rd_of_3_AUC';
fn3 = 'h37_2nd3rd';
load(fn3);
load(fn1);
AUCofNovLef = AUC_Out.ToneROCout_LowCorrTri_NovUnnovTyp.AUC;
SIGofNovLef = AUC_Out.ToneROCout_LowCorrTri_NovUnnovTyp.SigValue;
AUCofNovRig = AUC_Out.ToneROCout_HighCorrTri_NovUnnovTyp.AUC;
SIGofNovRIg = AUC_Out.ToneROCout_HighCorrTri_NovUnnovTyp.SigValue;
% make all value >= 0.5
% is an option
AUCofNovLef(AUCofNovLef<0.5)=1-AUCofNovLef(AUCofNovLef<0.5);
AUCofNovRig(AUCofNovRig<0.5)=1-AUCofNovRig(AUCofNovRig<0.5);

FrameTime = AUC_Out.FrameTime;
MinOnsFrame = AUC_Out.MinOnsFrame;

            % get differents of AUC  use abs and max
% Dff_AUC_Lef = max(abs(AUCofNovLef(:,TargFrames)-squeeze(SIGofNovLef(:,TargFrames,1))),[],2);
% Dff_AUC_Rig = max(abs(AUCofNovRig(:,TargFrames)-squeeze(SIGofNovRIg(:,TargFrames,1))),[],2);
% AUC_Lef = max(abs(AUCofNovLef(:,TargFrames)),[],2);
% AUC_Rig = max(abs(AUCofNovRig(:,TargFrames)),[],2);
            % get differents of AUC  don't make to abs and max
% Dff_AUC_Lef = max(AUCofNovLef(:,TargFrames)-squeeze(SIGofNovLef(:,TargFrames,1)),[],2);
% Dff_AUC_Rig = max(AUCofNovRig(:,TargFrames)-squeeze(SIGofNovRIg(:,TargFrames,1)),[],2);
% AUC_Lef = max(AUCofNovLef(:,TargFrames),[],2);
% AUC_Rig = max(AUCofNovRig(:,TargFrames),[],2);
            % get difference of AUC  don't make to abs and mean
PreLefDif = AUCofNovLef'-SIGofNovLef(:,1);
PreRigDif = AUCofNovRig'-SIGofNovRIg(:,1);
PreLef = AUCofNovLef';
PreRig = AUCofNovRig';

load(fn2);
AUCofNovLef = AUC_Out.ToneROCout_LowCorrTri_NovUnnovTyp.AUC;
SIGofNovLef = AUC_Out.ToneROCout_LowCorrTri_NovUnnovTyp.SigValue;
AUCofNovRig = AUC_Out.ToneROCout_HighCorrTri_NovUnnovTyp.AUC;
SIGofNovRIg = AUC_Out.ToneROCout_HighCorrTri_NovUnnovTyp.SigValue;
% make all value >= 0.5
% is an option
AUCofNovLef(AUCofNovLef<0.5)=1-AUCofNovLef(AUCofNovLef<0.5);
AUCofNovRig(AUCofNovRig<0.5)=1-AUCofNovRig(AUCofNovRig<0.5);
% Assign to second group
ProLefDif = AUCofNovLef'-SIGofNovLef(:,1);
ProRigDif = AUCofNovRig'-SIGofNovRIg(:,1);
ProLef = AUCofNovLef';
ProRig = AUCofNovRig';

LefDiffTimeOrd = [PreLefDif ProLefDif];
RigDiffTimeOrd = [PreRigDif ProRigDif];
LefAucTimeOrd = [PreLef ProLef];
RigAucTimeOrd = [PreRig ProRig];

LefAucSigTimeOrd(:,1) = PreLefDif >0; % significant than threshold for Current day and Former day
LefAucSigTimeOrd(:,2) = ProLefDif >0;   % the first colum is Former, second is Current
LefAucSigTimeOrd(:,3) = ProLefDif >0 & PreLefDif >0;
RigAucSigTimeOrd(:,1) = PreRigDif >0; 
RigAucSigTimeOrd(:,2) = ProRigDif >0;
RigAucSigTimeOrd(:,3) = ProRigDif >0 & PreRigDif >0;

if LefImprov>=0
    LefAucSigPerfOrd(:,1) = PreLefDif >0; % significant than threshold for Improve day and Impaire day; 
    LefAucSigPerfOrd(:,2) = ProLefDif >0;   % the first colum is impair, second is improve
    LefDiffPerOrd = [PreLefDif ProLefDif];
    LefAucPerOrd = [PreLef ProLef];
else
    LefAucSigPerfOrd(:,2) = PreLefDif >0; 
    LefAucSigPerfOrd(:,1) = ProLefDif >0;
    LefDiffPerOrd = [ProLefDif PreLefDif];
    LefAucPerOrd = [ProLef PreLef];
end
LefAucSigPerfOrd(:,3) = ProLefDif >0 & PreLefDif >0;
if RigImprov>=0
    RigAucSigPerfOrd(:,1) = PreRigDif >0; 
    RigAucSigPerfOrd(:,2) = ProRigDif >0;
    RigDiffPerOrd = [PreRigDif ProRigDif];
    RigAucPerOrd = [PreRig ProRig];
else
    RigAucSigPerfOrd(:,2) = PreRigDif >0; 
    RigAucSigPerfOrd(:,1) = ProRigDif >0;
    RigDiffPerOrd = [ProRigDif PreRigDif];
    RigAucPerOrd = [ProRig PreRig];
end
RigAucSigPerfOrd(:,3) = ProRigDif >0 & PreRigDif >0;

%% for actual auc time order
figure; hold on;axis square; xlim([0.3 1]); ylim([0.3 1]); line([0 1],[0 1],'linestyle',':','color','k','linewidth',1);
line([-1 1],[0.5 0.5],'linestyle',':','color','k','linewidth',1); line([0.5 0.5],[-1 1],'linestyle',':','color','k','linewidth',1);
plot(LefAucTimeOrd(:,1),LefAucTimeOrd(:,2),'k.','markersize',20); plot(LefAucTimeOrd(:,1),LefAucTimeOrd(:,2),'w.','markersize',5);
plot(RigAucTimeOrd(:,1),RigAucTimeOrd(:,2),'r.','markersize',20); plot(RigAucTimeOrd(:,1),RigAucTimeOrd(:,2),'w.','markersize',5);
title('Actual auc time order');
% for actual auc Perform order
figure; hold on;axis square; xlim([0.3 1]); ylim([0.3 1]); line([0 1],[0 1],'linestyle',':','color','k','linewidth',1);
line([-1 1],[0.5 0.5],'linestyle',':','color','k','linewidth',1); line([0.5 0.5],[-1 1],'linestyle',':','color','k','linewidth',1);
plot(LefAucPerOrd(:,1),LefAucPerOrd(:,2),'k.','markersize',20); plot(LefAucPerOrd(:,1),LefAucPerOrd(:,2),'w.','markersize',5);
plot(RigAucPerOrd(:,1),RigAucPerOrd(:,2),'r.','markersize',20); plot(RigAucPerOrd(:,1),RigAucPerOrd(:,2),'w.','markersize',5);
title('Actual auc Perform order');
% for diff auc time order
figure; hold on;axis square; xlim([-0.2 0.5]); ylim([-0.2 0.5]); line([-1 1],[-1 1],'linestyle',':','color','k','linewidth',1);
line([-1 1],[0 0],'linestyle',':','color','k','linewidth',1); line([0 0],[-1 1],'linestyle',':','color','k','linewidth',1);
plot(LefDiffTimeOrd(:,1),LefDiffTimeOrd(:,2),'k.','markersize',20); plot(LefDiffTimeOrd(:,1),LefDiffTimeOrd(:,2),'w.','markersize',5);
plot(RigDiffTimeOrd(:,1),RigDiffTimeOrd(:,2),'r.','markersize',20); plot(RigDiffTimeOrd(:,1),RigDiffTimeOrd(:,2),'w.','markersize',5);
title('Diff auc time order');
% for diff auc Perform order
figure; hold on;axis square; xlim([-0.2 0.5]); ylim([-0.2 0.5]); line([-1 1],[-1 1],'linestyle',':','color','k','linewidth',1);
line([-1 1],[0 0],'linestyle',':','color','k','linewidth',1); line([0 0],[-1 1],'linestyle',':','color','k','linewidth',1);
plot(LefDiffPerOrd(:,1),LefDiffPerOrd(:,2),'k.','markersize',20); plot(LefDiffPerOrd(:,1),LefDiffPerOrd(:,2),'w.','markersize',5);
plot(RigDiffPerOrd(:,1),RigDiffPerOrd(:,2),'r.','markersize',20); plot(RigDiffPerOrd(:,1),RigDiffPerOrd(:,2),'w.','markersize',5);
title('Diff auc Perform order');
%%
[h1,p1] = ttest(PreLefDif,ProLef)
[h2,p2] = ttest(PreRigDif,ProRig)
%%
if exist('AUCvsScore_UnAbs.mat')
    load('AUCvsScore_UnAbs.mat');
    Num = length(AUCvsScore);
    n =Num+1;
else
    n=1;
end


AUCvsScore(n).name = fn3;
AUCvsScore(n).RoiNum = size(LefDiffTimeOrd,1);
AUCvsScore(n).LeftImpro = LefImprov;
AUCvsScore(n).RightImpro = RigImprov;

AUCvsScore(n).LefDiffTimeOrd = LefDiffTimeOrd;
AUCvsScore(n).RigDiffTimeOrd = RigDiffTimeOrd;
AUCvsScore(n).LefAucTimeOrd = LefAucTimeOrd;
AUCvsScore(n).RigAucTimeOrd = RigAucTimeOrd;
AUCvsScore(n).LefAucSigTimeOrd = LefAucSigTimeOrd;
AUCvsScore(n).RigAucSigTimeOrd = RigAucSigTimeOrd;

AUCvsScore(n).LefDiffPerOrd = LefDiffPerOrd;
AUCvsScore(n).RigDiffPerOrd = RigDiffPerOrd;
AUCvsScore(n).LefAucPerOrd = LefAucPerOrd;
AUCvsScore(n).RigAucPerOrd = RigAucPerOrd;
AUCvsScore(n).LefAucSigPerfOrd = LefAucSigPerfOrd;
AUCvsScore(n).RigAucSigPerfOrd = RigAucSigPerfOrd;

save('AUCvsScore_UnAbs.mat','AUCvsScore')
%%  use diff
% load('AUCvsScore_UnAbs.mat');

% %%
figure;hold on;axis square; xlim([-0.2 0.4]); ylim([-0.2 0.4]);
line([-1 1],[-1 1],'linestyle','-','color','k','linewidth',2); 
line([-1 1],[0 0],'linestyle','-.','color',[.7 .7 .7],'linewidth',2); line([0 0],[-1 1],'linestyle','-.','color',[.7 .7 .7],'linewidth',2);
scatter(L_DiffTimeOrd(:,1),L_DiffTimeOrd(:,2),80,'k','filled'); alpha 0.1;      title('AuRoc Sig-Index(Low Freqs)');
set(gca,'fontsize',15,'fontweight','bold'); xlabel('Impaired Sessions'); ylabel('Improved Sessions');

Ind = L_AucSigTimeOrd(:,1) | L_AucSigTimeOrd(:,2);
figure;hold on;axis square;
plot([1 2],L_DiffTimeOrd(Ind,1:2));
% figure;hold on;axis square; xlim([0 0.4]); ylim([0 0.4]);
% line([0 1],[0 1],'linestyle',':','color','k','linewidth',1); scatter(tempRigCurr,tempRigForm,80,'k','filled'); alpha 0.1;      title('UnSort Right');
% figure;hold on;axis square; xlim([0 0.4]); ylim([0 0.4]);
% line([0 1],[0 1],'linestyle',':','color','k','linewidth',1); scatter(tempLefImpair,tempLefImprov,80,'k','filled'); alpha 0.1;      title('Sort Left');
% figure;hold on;axis square; xlim([0 0.4]); ylim([0 0.4]);
% line([0 1],[0 1],'linestyle',':','color','k','linewidth',1); scatter(tempRigImpair,tempRigImprov,80,'k','filled'); alpha 0.1;      title('Sort Right');

% figure;hold on;axis square; xlim([0 0.4]); ylim([0 0.4]);
% line([0 1],[0 1],'linestyle',':','color','k','linewidth',3);
% for pl = 1:length(tempLefImpair)
%     plot(tempLefImpair(pl),tempLefImprov(pl),'wo','markersize',8);
%     scatter(tempLefImpair(pl),tempLefImprov(pl),80,'k','filled');
% %     plot(tempLefImpair(pl),tempLefImprov(pl),'.','markersize',25,'color',[.5 .5 .5]); 
% end
% alpha 0.4;
% title('Sort Left');
% figure;hold on;axis square; xlim([0 0.4]); ylim([0 0.4]);
% line([0 1],[0 1],'linestyle',':','color','k','linewidth',1); scatter(tempRigImpair,tempRigImprov,80,'k','filled'); alpha 0.2;      title('Sort Right');
%%  use raw AUC
load('AUCvsScore_UnAbs_mean.mat');
tempLefImprov =[];
tempLefImpair =[];
tempRigImprov =[];
tempRigImpair =[];
for i = 1:length(AUCvsScore)
    if AUCvsScore(i).LeftImpro>=0
        tempLefImprov = [tempLefImprov;AUCvsScore(i).ProLef];
        tempLefImpair = [tempLefImpair;AUCvsScore(i).PreLef];
    else
        tempLefImprov = [tempLefImprov;AUCvsScore(i).PreLef];
        tempLefImpair = [tempLefImpair;AUCvsScore(i).ProLef];
    end
    if AUCvsScore(i).RightImpro>=0
        tempRigImprov = [tempRigImprov;AUCvsScore(i).ProRig];
        tempRigImpair = [tempRigImpair;AUCvsScore(i).PreRig];
    else
        tempRigImprov = [tempRigImprov;AUCvsScore(i).PreRig];
        tempRigImpair = [tempRigImpair;AUCvsScore(i).ProRig];
    end    
end
%%
figure;hold on;axis square;
xlim([0.5 0.9]);
ylim([0.5 0.9]);
line([0 1],[0 1],'linestyle',':','color','k','linewidth',1);
scatter(tempLefImpair,tempLefImprov,80,'k','filled');
% plot(PreRigDif,ProRig,'w.','markersize',5);
alpha 0.1;
%%
for j = 1:length(AUCvsScore)
    LefImp(j) = AUCvsScore(j).LeftImpro;
    RigImp(j) = AUCvsScore(j).RightImpro;
end    
    VarIndLeftMean_Perf = mean(abs(LefImp));
    VarIndLeftSEM_Perf = std(abs(LefImp))/(length(AUCvsScore))^0.5;
    VarIndRightMean_Perf = mean(abs(RigImp));
    VarIndRightSEM_Perf = std(abs(RigImp))/(length(AUCvsScore))^0.5;    

    VarIndLeftMean_Time = mean(LefImp);
    VarIndLeftSEM_Time = std(LefImp)/(length(AUCvsScore))^0.5;
    VarIndRightMean_Time = mean(RigImp);
    VarIndRightSEM_Time = std(RigImp)/(length(AUCvsScore))^0.5;     

    
    

%% transform struct to cell array
LefDiffTimeOrd = arrayfun(@(x) x.LefDiffTimeOrd,AUCvsScore,'UniformOutput',0);
RigDiffTimeOrd = arrayfun(@(x) x.RigDiffTimeOrd,AUCvsScore,'UniformOutput',0);
LefAucTimeOrd = arrayfun(@(x) x.LefAucTimeOrd,AUCvsScore,'UniformOutput',0);
RigAucTimeOrd = arrayfun(@(x) x.RigAucTimeOrd,AUCvsScore,'UniformOutput',0);
LefAucSigTimeOrd = arrayfun(@(x) x.LefAucSigTimeOrd,AUCvsScore,'UniformOutput',0);
RigAucSigTimeOrd = arrayfun(@(x) x.RigAucSigTimeOrd,AUCvsScore,'UniformOutput',0);
LefDiffPerOrd = arrayfun(@(x) x.LefDiffPerOrd,AUCvsScore,'UniformOutput',0);
RigDiffPerOrd = arrayfun(@(x) x.RigDiffPerOrd,AUCvsScore,'UniformOutput',0);
LefAucPerOrd = arrayfun(@(x) x.LefAucPerOrd,AUCvsScore,'UniformOutput',0);
RigAucPerOrd = arrayfun(@(x) x.RigAucPerOrd,AUCvsScore,'UniformOutput',0);
LefAucSigPerfOrd = arrayfun(@(x) x.LefAucSigPerfOrd,AUCvsScore,'UniformOutput',0);
RigAucSigPerfOrd = arrayfun(@(x) x.RigAucSigPerfOrd,AUCvsScore,'UniformOutput',0);
%%     transform cell to matix array
L_AucSigPerfOrd = [];
R_AucSigPerfOrd = [];
L_AucSigTimeOrd = [];
R_AucSigTimeOrd = [];
L_DiffPerfOrd = [];
R_DiffPerfOrd = [];
L_DiffTimeOrd = [];
R_DiffTimeOrd = [];
L_AucPerfOrd = [];
R_AucPerfOrd = [];
L_AucTimeOrd = [];
R_AucTimeOrd = [];
for i =1:length(AUCvsScore)
L_AucSigPerfOrd = [L_AucSigPerfOrd; LefAucSigPerfOrd{i}];
R_AucSigPerfOrd = [R_AucSigPerfOrd; RigAucSigPerfOrd{i}];
L_AucSigTimeOrd = [L_AucSigTimeOrd; LefAucSigTimeOrd{i}];
R_AucSigTimeOrd = [R_AucSigTimeOrd; RigAucSigTimeOrd{i}];
L_DiffPerfOrd = [L_DiffPerfOrd; LefDiffPerOrd{i}];
R_DiffPerfOrd = [R_DiffPerfOrd; RigDiffPerOrd{i}];
L_DiffTimeOrd = [L_DiffTimeOrd; LefDiffTimeOrd{i}];
R_DiffTimeOrd = [R_DiffTimeOrd; RigDiffTimeOrd{i}];
L_AucPerfOrd = [L_AucPerfOrd; LefAucPerOrd{i}];
R_AucPerfOrd = [R_AucPerfOrd; RigAucPerOrd{i}];
L_AucTimeOrd = [L_AucTimeOrd; LefAucTimeOrd{i}];
R_AucTimeOrd = [R_AucTimeOrd; RigAucTimeOrd{i}];
end
%%  plotting 
Ind_L_T_1 = L_AucSigTimeOrd(:,1);
Ind_R_T_1 = R_AucSigTimeOrd(:,1);
Ind_L_T_2 = L_AucSigTimeOrd(:,2);
Ind_R_T_2 = R_AucSigTimeOrd(:,2);
Ind_L_T_3 = L_AucSigTimeOrd(:,3);
Ind_R_T_3 = R_AucSigTimeOrd(:,3);


temp_L_T_1 = find(Ind_L_T_3==1);
temp_R_T_1 = find(Ind_R_T_3==1);

% for diff auc time order
figure; hold on;axis square; xlim([-0.2 0.5]); ylim([-0.2 0.5]); line([-1 1],[-1 1],'linestyle',':','color','k','linewidth',1);
line([-1 1],[0 0],'linestyle',':','color','k','linewidth',1); line([0 0],[-1 1],'linestyle',':','color','k','linewidth',1);
plot(L_DiffTimeOrd(temp_L_T_1,1),L_DiffTimeOrd(temp_L_T_1,2),'k.','markersize',20); plot(L_DiffTimeOrd(temp_L_T_1,1),L_DiffTimeOrd(temp_L_T_1,2),'w.','markersize',5);
plot(R_DiffTimeOrd(temp_R_T_1,1),R_DiffTimeOrd(temp_R_T_1,2),'r.','markersize',20); plot(R_DiffTimeOrd(temp_R_T_1,1),R_DiffTimeOrd(temp_R_T_1,2),'w.','markersize',5);
title('Diff auc time order');
%%
Ind_L_P_1 = L_AucSigPerfOrd(:,1);
Ind_R_P_1 = R_AucSigPerfOrd(:,1);
Ind_L_P_2 = L_AucSigPerfOrd(:,2);
Ind_R_P_2 = R_AucSigPerfOrd(:,2);
Ind_L_P_3 = L_AucSigPerfOrd(:,3);
Ind_R_P_3 = R_AucSigPerfOrd(:,3);


temp_L_P = find(Ind_L_P_2==1);
temp_R_P = find(Ind_R_P_2==1);

% for diff auc Perform order
figure; hold on;axis square; xlim([-0.2 0.5]); ylim([-0.2 0.5]); line([-1 1],[-1 1],'linestyle',':','color','k','linewidth',1);
line([-1 1],[0 0],'linestyle',':','color','k','linewidth',1); line([0 0],[-1 1],'linestyle',':','color','k','linewidth',1);
plot(L_DiffPerfOrd(temp_L_P,1),L_DiffPerfOrd(temp_L_P,2),'k.','markersize',20); plot(L_DiffPerfOrd(temp_L_P,1),L_DiffPerfOrd(temp_L_P,2),'w.','markersize',5);
plot(R_DiffPerfOrd(temp_R_P,1),R_DiffPerfOrd(temp_R_P,2),'r.','markersize',20); plot(R_DiffPerfOrd(temp_R_P,1),R_DiffPerfOrd(temp_R_P,2),'w.','markersize',5);
title('Diff auc Perform order');












