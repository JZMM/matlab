%% tone data
F = dir('*.mat');

for i = 1:length(F)
 fn = F(i).name;
 load(fn);
 AnimName{i} = fn([1:4 14:16]);
 T_L_data = AUC_Out.ToneROCout_LowCorrTri_NovUnnovTyp.AUC;
 T_R_data = AUC_Out.ToneROCout_HighCorrTri_NovUnnovTyp.AUC;
 T_L_label = AUC_Out.ToneROCout_LowCorrTri_NovUnnovTyp.ProOverPre;
 T_R_label = AUC_Out.ToneROCout_HighCorrTri_NovUnnovTyp.ProOverPre;
 T_L_sig = AUC_Out.ToneROCout_LowCorrTri_NovUnnovTyp.SigValue;
 T_R_sig = AUC_Out.ToneROCout_HighCorrTri_NovUnnovTyp.SigValue;
 
 T_L_sigInd = T_L_data - T_L_sig(:,1)'>0;
 T_R_sifInd = T_R_data - T_R_sig(:,1)'>0; 
 T_L_data(~T_L_label) = 1 - T_L_data(~T_L_label);
 T_R_data(~T_R_label) = 1 - T_R_data(~T_R_label);
 T_L_data = 2*(T_L_data - 0.5);
 T_R_data = 2*(T_R_data - 0.5);


L_NonSigData{i} = T_L_data(~T_L_sigInd);
R_NonSigData{i} = T_R_data(~T_R_sifInd);
L_Prob{i} = T_L_data(T_L_sigInd & T_L_data>0);
L_Train{i} = T_L_data(T_L_sigInd & T_L_data<0);
R_Prob{i} = T_R_data(T_R_sifInd & T_R_data>0);
R_Train{i} = T_R_data(T_R_sifInd & T_R_data<0);
NeuronNum(i) = length(T_L_data);
end
% NonSig.Low = L_NonSigData;
% NonSig.High = R_NonSigData;
% ProbOverTrain.Low = L_Prob;
% ProbOverTrain.High = R_Prob;
% TrainbOverPro.Low = L_Train;
% TrainbOverPro.High = R_Train;
%%                                                       
binSize = 10;
figure; hold on; set(gcf,'position',[2000 100 400 400]); xlim([-0.5 0.5]);
set(gca,'position',[0.2 0.2 0.2 0.4],'fontsize',15,'fontweight','bold');
[h1,c1] = hist(L_NonSigData,binSize);
[h2,c2] = hist(L_Prob,binSize/2);
[h3,c3] = hist(L_Train,binSize/2);
h1 = h1/length(L_Data)*100;
h2 = h2/length(L_Data)*100;
h3 = h3/length(L_Data)*100;
bar(c1,h1,'facecolor',[.5 .5 .5]);
bar(c2,h2,'facecolor','b');
bar(c3,h3,'facecolor','c');

figure; hold on; set(gcf,'position',[2000 100 400 400]);xlim([-0.5 0.5]);
set(gca,'position',[0.2 0.2 0.2 0.4],'fontsize',15,'fontweight','bold');
[H1,C1] = hist(R_NonSigData,binSize);
[H2,C2] = hist(R_Prob,binSize/2);
[H3,C3] = hist(R_Train,binSize/2);
H1 = H1/length(R_Data)*100;
H2 = H2/length(R_Data)*100;
H3 = H3/length(R_Data)*100;
bar(C1,H1,'facecolor',[.5 .5 .5]);
bar(C2,H2,'facecolor','b');
bar(C3,H3,'facecolor','c');
%%
L_pre(2,:) = [sum(h2) sum(h3)];
R_pre(2,:) = [sum(H2) sum(H3)];
%%
L_pro(2,:) = [sum(h2) sum(h3)];
R_pro(2,:) = [sum(H2) sum(H3)];
%% act data
F = dir('*.mat');
for i = 1:length(F)
 fn = F(i).name;
 load(fn);
 A_L_data = AUC_Out.ActROCout_LowCorrTri_NovUnnovTyp.AUC;
 A_R_data = AUC_Out.ActROCout_HighCorrTri_NovUnnovTyp.AUC;
 A_L_label = AUC_Out.ActROCout_LowCorrTri_NovUnnovTyp.ProOverPre;
 A_R_label = AUC_Out.ActROCout_HighCorrTri_NovUnnovTyp.ProOverPre;
 A_L_sig = AUC_Out.ActROCout_LowCorrTri_NovUnnovTyp.SigValue;
 A_R_sig = AUC_Out.ActROCout_HighCorrTri_NovUnnovTyp.SigValue;
 
 A_L_All(i,:) = [sum(A_L_data-A_L_sig(:,1)'>0)/length(A_L_data) sum(A_L_data-A_L_sig(:,1)'>0 & A_L_label)/length(A_L_data) sum(A_L_data-A_L_sig(:,1)'>0 & ~A_L_label)/length(A_L_data)];
 A_R_All(i,:) = [sum(A_R_data-A_R_sig(:,1)'>0)/length(A_R_data) sum(A_R_data-A_R_sig(:,1)'>0 & A_R_label)/length(A_R_data) sum(A_R_data-A_R_sig(:,1)'>0 & ~A_R_label)/length(A_R_data)];

end

%%
%% tone data
F = dir('*.mat');

for i = 1:length(F)
 fn = F(i).name;
 load(fn);
 T_L_data = AUC_Out.ToneROCout_LowCorrTri_NovUnnovTyp.AUC;
 T_R_data = AUC_Out.ToneROCout_HighCorrTri_NovUnnovTyp.AUC;
 T_L_label = AUC_Out.ToneROCout_LowCorrTri_NovUnnovTyp.ProOverPre;
 T_R_label = AUC_Out.ToneROCout_HighCorrTri_NovUnnovTyp.ProOverPre;
 T_L_sig = AUC_Out.ToneROCout_LowCorrTri_NovUnnovTyp.SigValue;
 T_R_sig = AUC_Out.ToneROCout_HighCorrTri_NovUnnovTyp.SigValue;
 
 T_L_sigInd = T_L_data - T_L_sig(:,1)'>0;
 T_R_sifInd = T_R_data - T_R_sig(:,1)'>0; 
 T_L_data(~T_L_label) = 1 - T_L_data(~T_L_label);
 T_R_data(~T_R_label) = 1 - T_R_data(~T_R_label);
 T_L_data = 2*(T_L_data - 0.5);
 T_R_data = 2*(T_R_data - 0.5);

 
L_NonSigData = T_L_data(~T_L_sigInd);
R_NonSigData = T_R_data(~T_R_sifInd);
L_Prob = T_L_data(T_L_sigInd & T_L_data>0);
L_Train = T_L_data(T_L_sigInd & T_L_data<0);
R_Prob = T_R_data(T_R_sifInd & T_R_data>0);
R_Train = T_R_data(T_R_sifInd & T_R_data<0);
 
end
%%

Grp1_1 = [1 4 6 9 12];
Grp1_2 = [2 5 7 10 13];
Grp2_1 = [2 7 10 13 15];
Grp2_2 = [3 8 11 14 16];
for i = 1:length(Grp1_1)
    Num_Lprob_1st2nd_1(i) = length(L_Prob{Grp1_1(i)});
    Num_Lprob_1st2nd_2(i) = length(L_Prob{Grp1_2(i)});
    Num_Rprob_1st2nd_1(i) = length(R_Prob{Grp1_1(i)});
    Num_Rprob_1st2nd_2(i) = length(R_Prob{Grp1_2(i)});    
    
    NumFrac_Lprob_1st2nd_1(i) = Num_Lprob_1st2nd_1(i)/NeuronNum(Grp1_1(i));
    NumFrac_Lprob_1st2nd_2(i) = Num_Lprob_1st2nd_2(i)/NeuronNum(Grp1_2(i));
    NumFrac_Rprob_1st2nd_1(i) = Num_Rprob_1st2nd_1(i)/NeuronNum(Grp1_1(i));
    NumFrac_Rprob_1st2nd_2(i) = Num_Rprob_1st2nd_2(i)/NeuronNum(Grp1_2(i)); 
    Perf_1st2nd_1(i,:) = NormScore(Grp1_1(i));
    Perf_1st2nd_2(i,:) = NormScore(Grp1_2(i));
end
for j = 1:length(Grp2_1)
    Num_Lprob_2nd3rd_1(j) = length(L_Prob{Grp2_1(j)});
    Num_Lprob_2nd3rd_2(j) = length(L_Prob{Grp2_2(j)});
    Num_Rprob_2nd3rd_1(j) = length(R_Prob{Grp2_1(j)});
    Num_Rprob_2nd3rd_2(j) = length(R_Prob{Grp2_2(j)});    
    
    NumFrac_Lprob_2nd3rd_1(j) = Num_Lprob_2nd3rd_1(j)/NeuronNum(Grp2_1(j));
    NumFrac_Lprob_2nd3rd_2(j) = Num_Lprob_2nd3rd_2(j)/NeuronNum(Grp2_2(j));
    NumFrac_Rprob_2nd3rd_1(j) = Num_Rprob_2nd3rd_1(j)/NeuronNum(Grp2_1(j));
    NumFrac_Rprob_2nd3rd_2(j) = Num_Rprob_2nd3rd_2(j)/NeuronNum(Grp2_2(j));    
    
    Perf_2nd3rd_1(j,:) = NormScore(Grp2_1(j));
    Perf_2nd3rd_2(j,:) = NormScore(Grp2_2(j));    
end

LefFracDiff = [NumFrac_Lprob_1st2nd_1 NumFrac_Lprob_2nd3rd_1]-[NumFrac_Lprob_1st2nd_2 NumFrac_Lprob_2nd3rd_2];
RigFracDiff = [NumFrac_Rprob_1st2nd_1 NumFrac_Rprob_2nd3rd_1]-[NumFrac_Rprob_1st2nd_2 NumFrac_Rprob_2nd3rd_2];

PerfDiff = [Perf_1st2nd_1;Perf_2nd3rd_1]-[Perf_1st2nd_2;Perf_2nd3rd_2];

figure;hold on;
plot(LefFracDiff',PerfDiff,'ko','linestyle','none');
% figure;hold on;
plot(RigFracDiff',PerfDiff,'ro','linestyle','none');

%%
for i = 1:length(L_Prob)
    Num_Lprob(i) = length(L_Prob{i});
    Num_Rprob(i) = length(R_Prob{i});   
    NumFrac_Lprob(i) = Num_Lprob(i)/NeuronNum(i);
    NumFrac_Rprob(i) = Num_Rprob(i)/NeuronNum(i);
    
    Num_Ltrain(i) = length(L_Train{i});
    Num_Rtrain(i) = length(R_Train{i});   
    NumFrac_Ltrain(i) = Num_Ltrain(i)/NeuronNum(i);
    NumFrac_Rtrain(i) = Num_Rtrain(i)/NeuronNum(i);   
end
figure;hold on;
plot(NumFrac_Lprob',ScoreImpr_1st2ndBlock(:,2),'ko','linestyle','none');
plot(NumFrac_Rprob',ScoreImpr_1st2ndBlock(:,3),'ro','linestyle','none');
figure;hold on;
xlim([0.5 4.5]);
plot([1 2],[Score_1stBlock(:,2) Score_2ndBlock(:,2)],'color','k');
plot([3 4],[Score_1stBlock(:,3) Score_2ndBlock(:,3)],'color','r');
for j = 1:length(L_Prob)
    if NumFrac_Lprob(j)>0
    plot(1,Score_1stBlock(j,2),'ko','markersize',NumFrac_Lprob(j)*100);
    end
    if NumFrac_Ltrain(j)>0
        plot(2,Score_2ndBlock(j,2),'ko','markersize',NumFrac_Ltrain(j)*100);
    end
    if NumFrac_Rprob(j)>0
    plot(3,Score_1stBlock(j,3),'ro','markersize',NumFrac_Rprob(j)*100);
    end
    if NumFrac_Rtrain(j)>0
        plot(4,Score_2ndBlock(j,3),'ro','markersize',NumFrac_Rtrain(j)*100);
    end
end
%%
figure;hold on; AxisForm([],15);
a=polyfit(Score_2ndBlock(:,2)-Score_1stBlock(:,2),Score_2ndBlock(:,3)-Score_1stBlock(:,3),1);
x1 =linspace(-1,1);
y1=polyval(a,x1);
plot(x1,y1,'linewidth',3)
plot(Score_2ndBlock(:,2)-Score_1stBlock(:,2),Score_2ndBlock(:,3)-Score_1stBlock(:,3),'k.','markersize',25,'linestyle','none');
plot(Score_2ndBlock(:,2)-Score_1stBlock(:,2),Score_2ndBlock(:,3)-Score_1stBlock(:,3),'w.','markersize',6,'linestyle','none');


xlabel('Left Improvement');
ylabel('Right Improvement');
title('Probe Trials Perf');
%%
Grp1_1 = [1 4 6 9 12];
Grp1_2 = [2 5 7 10 13];
Grp2_1 = [2 7 10 13 15];
Grp2_2 = [3 8 11 14 16];

temp1 = NumFrac_Lprob;  % NumFrac_Lprob  NumFrac_Ltrain
temp2 = NumFrac_Rprob; % NumFrac_Rprob  NumFrac_Rtrain
figure;hold on;
xlim([0 5]);
plot([1 2],[Score_2ndBlock(Grp1_1,2)-Score_1stBlock(Grp1_1,2) Score_2ndBlock(Grp1_2,2)-Score_1stBlock(Grp1_2,2)],'k.','linestyle','-')
plot([2 3],[Score_2ndBlock(Grp2_1,2)-Score_1stBlock(Grp2_1,2) Score_2ndBlock(Grp2_2,2)-Score_1stBlock(Grp2_2,2)],'ko','linestyle','-')
L_FracDiff = [temp1(Grp1_2)-temp1(Grp1_1) temp1(Grp2_2)-temp1(Grp2_1)]';
L_PerfDiff = [Score_2ndBlock(Grp1_2,2)-Score_1stBlock(Grp1_2,2)-(Score_2ndBlock(Grp1_1,2)-Score_1stBlock(Grp1_1,2));...
    Score_2ndBlock(Grp2_2,2)-Score_1stBlock(Grp2_2,2)-(Score_2ndBlock(Grp2_1,2)-Score_1stBlock(Grp2_1,2))];

R_FracDiff = [temp2(Grp1_2)-temp2(Grp1_1) temp2(Grp2_2)-temp2(Grp2_1)]';
R_PerfDiff = [Score_2ndBlock(Grp1_2,3)-Score_1stBlock(Grp1_2,3)-(Score_2ndBlock(Grp1_1,3)-Score_1stBlock(Grp1_1,3));...
    Score_2ndBlock(Grp2_2,3)-Score_1stBlock(Grp2_2,3)-(Score_2ndBlock(Grp2_1,3)-Score_1stBlock(Grp2_1,3))];

% [h,p]=ttest2([Score_2ndBlock(Ind1,2)-Score_1stBlock(Ind1,2);Score_1stBlock(Ind1,3)-Score_2ndBlock(Ind1,3)],...
%     [Score_2ndBlock(Ind2,2)-Score_1stBlock(Ind2,2);Score_1stBlock(Ind2,3)-Score_2ndBlock(Ind2,3)])
[h,p]=ttest2([Score_2nd_Nor(Ind1,2)-Score_1st_Nor(Ind1,2);Score_1st_Nor(Ind1,3)-Score_2nd_Nor(Ind1,3)],...
    [Score_2nd_Nor(Ind2,2)-Score_1st_Nor(Ind2,2);Score_1st_Nor(Ind2,3)-Score_2nd_Nor(Ind2,3)])