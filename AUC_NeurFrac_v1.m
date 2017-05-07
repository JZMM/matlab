%% tone data
F = dir('*.mat');
L_ind = [];
R_ind = [];
L_Data = [];
R_Data = [];
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

 
 L_Data = [L_Data T_L_data];
 R_Data = [R_Data T_R_data]; 
 L_ind = [L_ind T_L_sigInd];
 R_ind = [R_ind T_R_sifInd]; 
 
end
L_NonSigData = L_Data(~L_ind);
R_NonSigData = R_Data(~R_ind);
L_Prob = L_Data(L_ind & L_Data>0);
L_Train = L_Data(L_ind & L_Data<0);
R_Prob = R_Data(R_ind & R_Data>0);
R_Train = R_Data(R_ind & R_Data<0);
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





