%%  11111111111111111111111111111111111111111    1st block 
files = dir('*.mat');
Bound = 16000;
for i = 1:length(files)
    fn = files(i).name;
    load(fn);
    PerfDataAnimName{i} = fn([4:7 20:21]);
    Action_choice = Data_extract.Action_choice(TestTrialNum);
    Tones = Data_extract.Tone_frequency(TestTrialNum);
    Freqs = unique(Tones);
    MisInd = Data_extract.Miss_Ind(TestTrialNum);
    TrainLow = Tones == Freqs(1) & ~MisInd;
    ProbeLow = Tones > Freqs(1) & ~MisInd & Tones < Bound;
    ProbeHigh = Tones < Freqs(end) & ~MisInd & Tones > Bound;
    TrainHigh = Tones == Freqs(end) & ~MisInd;
    
    TrainLowCorr = TrainLow & Action_choice==0;
    ProbeLowCorr = ProbeLow & Action_choice==0;
    ProbeHighCorr = ProbeHigh & Action_choice==1;
    TrainHighCorr = TrainHigh & Action_choice==1;   
    CorrRate(i,:) = [sum(TrainLowCorr)/sum(TrainLow) sum(ProbeLowCorr)/sum(ProbeLow) ...
        sum(ProbeHighCorr)/sum(ProbeHigh) sum(TrainHighCorr)/sum(TrainHigh)];
    RightChoiceScore = [1-CorrRate(i,1:2) CorrRate(i,3:4)];
%     NormScore(1) = (RightChoiceScore(2)-RightChoiceScore(1)-(max(RightChoiceScore)-RightChoiceScore(2)))/(max(RightChoiceScore)-RightChoiceScore(1));
%     NormScore(2) = (RightChoiceScore(3)-RightChoiceScore(1)-(max(RightChoiceScore)-RightChoiceScore(3)))/(max(RightChoiceScore)-RightChoiceScore(1));
    NormScore(i,:) = [(RightChoiceScore(2)-RightChoiceScore(1)-(RightChoiceScore(end)-RightChoiceScore(2)))/(RightChoiceScore(end)-RightChoiceScore(1)) ...
        (RightChoiceScore(3)-RightChoiceScore(1)-(RightChoiceScore(end)-RightChoiceScore(3)))/(RightChoiceScore(end)-RightChoiceScore(1))];

%     RightChoiceScore(1:2) = 1-RightChoiceScore(1:2);
    
end
%%  2222222222222222222222222222222222222222222222       block

files = dir('*.mat');
Bound = 16000;
for i = 1:length(files)
    fn = files(i).name;
    load(fn);
    AnimName{i} = fn([4:7 18:21]);
    Block1 = [1:round(2*length(TestTrialNum)/3)];
    Block2 = [round(length(TestTrialNum)/3):length(TestTrialNum)];
    Action_choice = Data_extract.Action_choice(TestTrialNum);
    Tones = Data_extract.Tone_frequency(TestTrialNum);
    Freqs = unique(Tones);
    MisInd = Data_extract.Miss_Ind(TestTrialNum);
    TrainLow = Tones == Freqs(1) & ~MisInd;
    ProbeLow = Tones > Freqs(1) & ~MisInd & Tones < Bound;
    ProbeHigh = Tones < Freqs(end) & ~MisInd & Tones > Bound;
    TrainHigh = Tones == Freqs(end) & ~MisInd;
    
    TrainLowCorr = TrainLow & Action_choice==0;
    ProbeLowCorr = ProbeLow & Action_choice==0;
    ProbeHighCorr = ProbeHigh & Action_choice==1;
    TrainHighCorr = TrainHigh & Action_choice==1;   
    Score_1stBlock(i,:) = [sum(TrainLowCorr(Block1))/sum(TrainLow(Block1)) sum(ProbeLowCorr(Block1))/sum(ProbeLow(Block1)) ...
        sum(ProbeHighCorr(Block1))/sum(ProbeHigh(Block1)) sum(TrainHighCorr(Block1))/sum(TrainHigh(Block1))];
    Score_2ndBlock(i,:) = [sum(TrainLowCorr(Block2))/sum(TrainLow(Block2)) sum(ProbeLowCorr(Block2))/sum(ProbeLow(Block2)) ...
        sum(ProbeHighCorr(Block2))/sum(ProbeHigh(Block2)) sum(TrainHighCorr(Block2))/sum(TrainHigh(Block2))];  

    Score_1st_Nor(i,:) = Score_1stBlock(i,:);
    Score_2nd_Nor(i,:) = Score_2ndBlock(i,:);
    Score_1st_Nor(i,1:2) = 1-Score_1st_Nor(i,1:2);
    Score_2nd_Nor(i,1:2) = 1-Score_2nd_Nor(i,1:2);
    Score_1st_Nor(i,:) = (Score_1st_Nor(i,:)-Score_1st_Nor(i,1))/(Score_1st_Nor(i,4)-Score_1st_Nor(i,1));
    Score_2nd_Nor(i,:) = (Score_2nd_Nor(i,:)-Score_2nd_Nor(i,1))/(Score_2nd_Nor(i,4)-Score_2nd_Nor(i,1));
    Score_1st_Nor(i,1:2) = 1-Score_1st_Nor(i,1:2);
    Score_2nd_Nor(i,1:2) = 1-Score_2nd_Nor(i,1:2);
    ScoreImpr_1st2ndBlock(i,:) = Score_2nd_Nor(i,:) - Score_1st_Nor(i,:);
%     save(fn,'Score_1stBlock','Score_2ndBlock','Score_1st_Nor','Score_2nd_Nor','ScoreImpr_1st2ndBlock','-append');
end
% save('Performance_2P.mat','Score_1stBlock','Score_2ndBlock','Score_1st_Nor','Score_2nd_Nor','ScoreImpr_1st2ndBlock','AnimName');
%%
files = dir('*.mat');
Bound = 16000;
for i = 1:length(files)
    fn = files(i).name;
    load(fn);
    Form_scor(i,:) = Score_1st_Nor(2:3);
    Latt_scor(i,:) = Score_2nd_Nor(2:3);
end
figure; hold on;  xlim([0.5 4.5]);
AxisForm([0.2 0.2 0.6 0.7],20,{0,[1 2 3 4],{'1stBlock','2ndBlock','1stBlock','2ndBlock'},20});
plot([1 2],[Form_scor(:,1) Latt_scor(:,1)],'color','k','linewidth',2);
plot([3 4],[Form_scor(:,2) Latt_scor(:,2)],'color','r','linewidth',2);
Mean1 = mean([Form_scor Latt_scor]);
SEMs = std([Form_scor Latt_scor])/(size(Form_scor,1))^0.5;
line([0.7 0.8],[Mean1(1) Mean1(1)],'color','k','linewidth',3);
line([2.2 2.3],[Mean1(3) Mean1(3)],'color','k','linewidth',3);
line([2.7 2.8],[Mean1(2) Mean1(2)],'color','r','linewidth',3);
line([4.2 4.3],[Mean1(4) Mean1(4)],'color','r','linewidth',3);
line([0.75 0.75],[Mean1(1)+SEMs(1) Mean1(1)-SEMs(1)],'color','k');
line([2.25 2.25],[Mean1(3)+SEMs(3) Mean1(3)-SEMs(3)],'color','k');
line([2.75 2.75],[Mean1(2)+SEMs(2) Mean1(2)-SEMs(2)],'color','r');
line([4.25 4.25],[Mean1(4)+SEMs(4) Mean1(4)-SEMs(4)],'color','r');
%%
files = dir('*.mat');
Bound = 16000;
for i = 1:length(files)
    fn = files(i).name;
    load(fn);
    L_scor(i,:) = ScoreImpr_1st2ndBlock(2);
    R_scor(i,:) = ScoreImpr_1st2ndBlock(3);
end
%%
L_Frac = [];
R_Frac = [];
files = dir('*.mat');
for h = 1:length(files)
    fn = files(h).name;
    load(fn);
    AUCofNovLef = AUC_Out.ToneROCout_LowCorrTri_NovUnnovTyp.AUC;
    SIGofNovLef = AUC_Out.ToneROCout_LowCorrTri_NovUnnovTyp.SigValue;
    AUCofNovRig = AUC_Out.ToneROCout_HighCorrTri_NovUnnovTyp.AUC;
    SIGofNovRIg = AUC_Out.ToneROCout_HighCorrTri_NovUnnovTyp.SigValue;

    PreLefDifAUC = AUCofNovLef'-SIGofNovLef(:,1);
    PreRigDifAUC = AUCofNovRig'-SIGofNovRIg(:,1);

    PreLefDifAUC_Frac(h) = sum(PreLefDifAUC>0)/size(PreLefDifAUC,1);
    PreRigDifAUC_Frac(h) = sum(PreRigDifAUC>0)/size(PreRigDifAUC,1);     
end
%%
figure; hold on;  
% xlim([0.5 4.5]);
AxisForm([0.2 0.2 0.6 0.7],20,{0,[1 2 3 4],{'1stBlock','2ndBlock','1stBlock','2ndBlock'},20});
plot(PreLefDifAUC_Frac,L_scor,'k.','markersize',10);
plot(PreRigDifAUC_Frac,R_scor,'r.','markersize',10);
figure; hold on;  
% xlim([0.5 4.5]);
AxisForm([0.2 0.2 0.6 0.7],20,{0,[1 2 3 4],{'1stBlock','2ndBlock','1stBlock','2ndBlock'},20});
plot(PreLefDifAUC_Frac,R_scor,'k.','markersize',10);
plot(PreRigDifAUC_Frac,L_scor,'r.','markersize',10);
figure; hold on;  
AxisForm([0.2 0.2 0.6 0.7],20,{0,[1 2 3 4],{'1stBlock','2ndBlock','1stBlock','2ndBlock'},20});
plot(-PreLefDifAUC_Frac+PreRigDifAUC_Frac,L_scor-R_scor,'k.','markersize',10);

