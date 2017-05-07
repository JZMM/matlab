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
% AUCofNovLef(AUCofNovLef<0.5)=1-AUCofNovLef(AUCofNovLef<0.5);
% AUCofNovRig(AUCofNovRig<0.5)=1-AUCofNovRig(AUCofNovRig<0.5);

PreLefDifAUC = AUCofNovLef'-SIGofNovLef(:,1);
PreRigDifAUC = AUCofNovRig'-SIGofNovRIg(:,1);
PreLefAUC = AUCofNovLef';
PreRigAUC = AUCofNovRig';

PreLefDifAUC_Frac = sum(PreLefDifAUC>0)/size(PreLefDifAUC,1);
PreRigDifAUC_Frac = sum(PreRigDifAUC>0)/size(PreRigDifAUC,1);

load(fn2);
AUCofNovLef = AUC_Out.ToneROCout_LowCorrTri_NovUnnovTyp.AUC;
SIGofNovLef = AUC_Out.ToneROCout_LowCorrTri_NovUnnovTyp.SigValue;
AUCofNovRig = AUC_Out.ToneROCout_HighCorrTri_NovUnnovTyp.AUC;
SIGofNovRIg = AUC_Out.ToneROCout_HighCorrTri_NovUnnovTyp.SigValue;
% make all value >= 0.5
% is an option
% AUCofNovLef(AUCofNovLef<0.5)=1-AUCofNovLef(AUCofNovLef<0.5);
% AUCofNovRig(AUCofNovRig<0.5)=1-AUCofNovRig(AUCofNovRig<0.5);

% Assign to second group
ProLefDifAUC = AUCofNovLef'-SIGofNovLef(:,1);
ProRigDifAUC = AUCofNovRig'-SIGofNovRIg(:,1);
ProLefAUC = AUCofNovLef';
ProRigAUC = AUCofNovRig';

ProLefDifAUC_Frac = sum(ProLefDifAUC>0)/size(ProLefDifAUC,1);
ProRigDifAUC_Frac = sum(ProRigDifAUC>0)/size(ProRigDifAUC,1);

LefDifAUC_Frac = [PreLefDifAUC_Frac ProLefDifAUC_Frac];
RigDifAUC_Frac = [PreRigDifAUC_Frac ProRigDifAUC_Frac];
save(fn3,'LefDifAUC_Frac','RigDifAUC_Frac','-append');
%%
files = dir('*.mat');
for i =1:length(files)
    fn = files(i).name;
    load(fn);
%     ScoreLef(i,1) = NormScore_1(1);  % use score
%     ScoreLef(i,2) = NormScore_2(1);
%     ScoreRig(i,1) = NormScore_1(2);
%     ScoreRig(i,2) = NormScore_2(2);

    ScoreLef(i,:) = LefDiff;  % use Diff
    ScoreRig(i,:) = RigDiff;

    
    LefAUC_Frac(i,:) = LefDifAUC_Frac;
    RigAUC_Frac(i,:) = RigDifAUC_Frac;
end
%%   1111111111111111111
figure; hold on; axis square; set(gca,'fontsize',15,'fontweight','bold');
set(gcf,'position',[2000 400 400 400],'color','w');
title('AUC_diff vs PerfDiff','interpreter','none');
xlabel('AUC Diff');
ylabel('Perf Diff');
for j =1:size(LefAUC_Frac,1)
    line(LefAUC_Frac(j,:),ScoreLef(j,:));
    plot(LefAUC_Frac(j,1),ScoreLef(j,1),'ko');
    plot(LefAUC_Frac(j,2),ScoreLef(j,2),'co','markersize',10);
end
%%   2222222222222
figure; hold on; axis square; set(gca,'fontsize',15,'fontweight','bold','position',[0.2 0.2 0.5 0.5]);
set(gcf,'position',[2000 400 400 400],'color','w');
title('AUC_diff vs PerfDiff','interpreter','none');
xlabel('-    Neurons Fraction    +');
ylabel('Perf Improve Index');
line([-1 1],[0 0],'color',[.7 .7 .7],'linewidth',2,'linestyle',':');
for j =1:size(LefAUC_Frac,1)
%     plot(LefAUC_Frac(j,2)-LefAUC_Frac(j,1),ScoreLef(j,1)-ScoreLef(j,2),'w.','markersize',35);
%     PL1 = plot(LefAUC_Frac(j,2)-LefAUC_Frac(j,1),ScoreLef(j,1)-ScoreLef(j,2),'k.','markersize',30);
%     plot(RigAUC_Frac(j,2)-RigAUC_Frac(j,1),ScoreRig(j,2)-ScoreRig(j,1),'w.','markersize',35);
%     PL2 = plot(RigAUC_Frac(j,2)-RigAUC_Frac(j,1),ScoreRig(j,2)-ScoreRig(j,1),'r.','markersize',30);
    plot(LefAUC_Frac(j,2)-LefAUC_Frac(j,1),ScoreLef(j,1)-ScoreLef(j,2),'w.','markersize',35);
    PL1 = plot(LefAUC_Frac(j,2)-LefAUC_Frac(j,1),ScoreLef(j,1)-ScoreLef(j,2),'k.','markersize',30);
    plot(RigAUC_Frac(j,2)-RigAUC_Frac(j,1),ScoreRig(j,2)-ScoreRig(j,1),'w.','markersize',35);
    PL2 = plot(RigAUC_Frac(j,2)-RigAUC_Frac(j,1),ScoreRig(j,2)-ScoreRig(j,1),'r.','markersize',30);    
end
% line([-1 1],[-1 1]);
xlim([-0.5 0.5]);
p_Lef = polyfit(LefAUC_Frac(:,2)-LefAUC_Frac(:,1),ScoreLef(:,1)-ScoreLef(:,2),1);
p_Rig = polyfit(RigAUC_Frac(:,2)-RigAUC_Frac(:,1),ScoreRig(:,2)-ScoreRig(:,1),1);
yL = polyval(p_Lef,[-1:0.1:1]);
yR = polyval(p_Rig,[-1:0.1:1]);
plot([-1:0.1:1],yL,'color','k','linewidth',2);
plot([-1:0.1:1],yR,'color','r','linewidth',2);
[r_L,p_L]=corrcoef(LefAUC_Frac(:,2)-LefAUC_Frac(:,1),ScoreLef(:,1)-ScoreLef(:,2));
[r_R,p_R]=corrcoef(RigAUC_Frac(:,2)-RigAUC_Frac(:,1),ScoreRig(:,2)-ScoreRig(:,1));
[r_All,p_All]=corrcoef([LefAUC_Frac(:,2)-LefAUC_Frac(:,1);RigAUC_Frac(:,2)-RigAUC_Frac(:,1)],[ScoreLef(:,2)-ScoreLef(:,1);ScoreRig(:,2)-ScoreRig(:,1)]);
legend([PL1 PL2],{'Low' ,'High'});
text(0,0,['r_L = ' num2str(r_L(1,2))],'fontsize',10,'fontweight','bold');
text(0,-0.2,['r_H = ' num2str(r_R(1,2))],'fontsize',10,'fontweight','bold','color','r');
%%
figure; hold on; axis square; set(gca,'fontsize',15,'fontweight','bold');
set(gcf,'position',[2000 400 400 400],'color','w');
title('AUC_diff vs PerfDiff','interpreter','none');
xlabel('-    Neurons Fraction    +');
ylabel('Perf Impaired');
% for j =1:size(LefAUC_Frac,1)
%     if ScoreLef(j,2)-ScoreLef(j,1) >0
%         line(LefAUC_Frac(j,:),ScoreLef(j,:));
%         plot(LefAUC_Frac(j,1),ScoreLef(j,1),'ko');
%         plot(LefAUC_Frac(j,2),ScoreLef(j,2),'co','markersize',10);
%     else
%         line(fliplr(LefAUC_Frac(j,:)),fliplr(ScoreLef(j,:)));
%         plot(LefAUC_Frac(j,2),ScoreLef(j,2),'ko');
%         plot(LefAUC_Frac(j,1),ScoreLef(j,1),'co','markersize',10);        
%     end
% end
% LefAUC_Frac(:,1) = LefAUC_Frac(:,1)-LefAUC_Frac(:,1);
% LefAUC_Frac(:,2) = LefAUC_Frac(:,2)-LefAUC_Frac(:,1);
% ScoreLef(:,1) = ScoreLef(:,1)-ScoreLef(:,1);
% ScoreLef(:,2) = ScoreLef(:,2)-ScoreLef(:,1);
for j =1:size(LefAUC_Frac,1)
    if ScoreLef(j,2)-ScoreLef(j,1) >0
        LefAUC_Frac(j,:) = LefAUC_Frac(j,:)-[LefAUC_Frac(j,1) LefAUC_Frac(j,1)];
        ScoreLef(j,:) = ScoreLef(j,:)-[ScoreLef(j,1) ScoreLef(j,1)];
        line(LefAUC_Frac(j,:),ScoreLef(j,:));
        plot(LefAUC_Frac(j,1),ScoreLef(j,1),'k.','markersize',30);
        plot(LefAUC_Frac(j,2),ScoreLef(j,2),'c.','markersize',20);
    else
        LefAUC_Frac(j,:) = LefAUC_Frac(j,:)-[LefAUC_Frac(j,2) LefAUC_Frac(j,2)];
        ScoreLef(j,:) = ScoreLef(j,:)-[ScoreLef(j,2) ScoreLef(j,2)];        
        line(fliplr(LefAUC_Frac(j,:)),fliplr(ScoreLef(j,:)));
        plot(LefAUC_Frac(j,2),ScoreLef(j,2),'k.','markersize',30);
        plot(LefAUC_Frac(j,1),ScoreLef(j,1),'c.','markersize',20);        
    end
end
line([0 mean(LefAUC_Frac(:,2))],[0,mean(ScoreLef(:,2))],'linewidth',4,'color','r');
plot(mean(LefAUC_Frac(:,2)),mean(ScoreLef(:,2)),'r.','markersize',30);
line([0 1],[0 0]);
xlim([-0.4 0.4]);
%%
figure; hold on; axis square; set(gca,'fontsize',15,'fontweight','bold','position',[0.2 0.2 0.5 0.6]);
set(gcf,'position',[2000 400 400 400],'color','w');
title('AUC_diff vs PerfDiff','interpreter','none');
% xlabel('AUC Diff');
% ylabel('Perf Diff');
for j =1:size(LefAUC_Frac,1)
%     line(LefAUC_Frac(j,:),ScoreLef(j,:));
    plot([1 2],LefAUC_Frac(j,:),'ko');
    line([1 2],LefAUC_Frac(j,:),'color','k');
    plot([4 5],ScoreLef(j,:),'co');
    line([4 5],ScoreLef(j,:),'color','c');
%     plot(LefAUC_Frac(j,:),'co','markersize',10);
end
%%  plot Fraction pre vs pro
figure; hold on; set(gca,'fontsize',20,'fontweight','bold','position',[0.2 0.2 0.5 0.7],'xtick',[1 2],'xticklabel',{'Former','Latter'},'xticklabelrotation',30);
xlim([0.5 2.5]);
ylim([-0.02 0.6])
set(gcf,'position',[2000 400 400 400],'color','w');
plot([1 2],LefAUC_Frac,'color','k','linewidth',2);
plot([1 2],RigAUC_Frac,'color','r','linewidth',2);
% plot(0.75,mean([LefAUC_Frac(:,1);RigAUC_Frac(:,1)]),'k.','markersize',30);
Mean1 = mean([LefAUC_Frac(:,1);RigAUC_Frac(:,1)]);
Mean2 = mean([LefAUC_Frac(:,2);RigAUC_Frac(:,2)]);
line([0.7 0.8],[Mean1 Mean1],'color','k','linewidth',3);
line([2.2 2.3],[Mean2 Mean2],'color','k','linewidth',3);
SEMs =[std([LefAUC_Frac(:,1);RigAUC_Frac(:,1)])/(length([LefAUC_Frac(:,1);RigAUC_Frac(:,1)]))^0.5 ...
    std([LefAUC_Frac(:,2);RigAUC_Frac(:,2)])/(length([LefAUC_Frac(:,2);RigAUC_Frac(:,2)]))^0.5];
line([0.75 0.75],[Mean1+SEMs(1) Mean1-SEMs(1)],'color','k');
line([2.25 2.25],[Mean2+SEMs(2) Mean2-SEMs(2)],'color','k');
ylabel('Neruons Fraction');
%%  plot Fraction Improve vs impair
figure; hold on; set(gca,'fontsize',20,'fontweight','bold','position',[0.2 0.2 0.5 0.7],'xtick',[1 2],'xticklabel',{'Former','Latter'},'xticklabelrotation',30);
xlim([0.5 2.5]);
ylim([-0.02 0.6])
set(gcf,'position',[2000 400 400 400],'color','w');
Ind_1 = ScoreLef(:,2)-ScoreLef(:,1);
Ind_2 = ScoreRig(:,2)-ScoreRig(:,1);
for tt = 1:size(ScoreLef)
    if Ind_1(tt)<0
        Lef(tt,:) = fliplr(LefAUC_Frac(tt,:));
    else
        Lef(tt,:) = LefAUC_Frac(tt,:);
    end
    if Ind_2(tt)<0
        Rig(tt,:) = fliplr(RigAUC_Frac(tt,:));
    else
        Rig(tt,:) = RigAUC_Frac(tt,:);
    end    
end
plot([1 2],Lef,'color','k','linewidth',2);
plot([1 2],Rig,'color','r','linewidth',2);
% plot(0.75,mean([LefAUC_Frac(:,1);RigAUC_Frac(:,1)]),'k.','markersize',30);
Mean1 = mean([Lef(:,1);Rig(:,1)]);
Mean2 = mean([Lef(:,2);Rig(:,2)]);
line([0.7 0.8],[Mean1 Mean1],'color','k','linewidth',3);
line([2.2 2.3],[Mean2 Mean2],'color','k','linewidth',3);
SEMs =[std([Lef(:,1);Rig(:,1)])/(length([Lef(:,1);Rig(:,1)]))^0.5 ...
    std([Lef(:,2);Rig(:,2)])/(length([Lef(:,2);Rig(:,2)]))^0.5];
line([0.75 0.75],[Mean1+SEMs(1) Mean1-SEMs(1)],'color','k');
line([2.25 2.25],[Mean2+SEMs(2) Mean2-SEMs(2)],'color','k');
ylabel('Neruons Fraction');
%% plot performance diff or Normalize score
figure; hold on; 
AxisForm([0.2 0.2 0.5 0.7],20,{0,[1 2],{'Former','Latter'},30});

xlim([0.5 2.5]);
% ylim([-0.02 0.6])
set(gcf,'position',[2000 400 400 400],'color','w');

plot([1 2],ScoreLef,'color','k','linewidth',2);
plot([1 2],ScoreRig,'color','r','linewidth',2);
Mean1 = mean(ScoreLef);
Mean2 = mean(ScoreRig);
SEM1 = std(ScoreLef)/(size(ScoreLef,1))^0.5;
SEM2 = std(ScoreRig)/(size(ScoreRig,1))^0.5;
line([0.7 0.8],[Mean1(1) Mean1(1)],'color','k','linewidth',3);
line([2.2 2.3],[Mean1(2) Mean1(2)],'color','k','linewidth',3);
line([0.7 0.8],[Mean2(1) Mean2(1)],'color','r','linewidth',3);
line([2.2 2.3],[Mean2(2) Mean2(2)],'color','r','linewidth',3);
line([0.75 0.75],[Mean1(1)+SEM1(1) Mean1(1)-SEM1(1)],'color','k');
line([2.25 2.25],[Mean1(2)+SEM1(2) Mean1(2)-SEM1(2)],'color','k');
line([0.75 0.75],[Mean2(1)+SEM2(1) Mean2(1)-SEM2(1)],'color','r');
line([2.25 2.25],[Mean2(2)+SEM2(2) Mean2(2)-SEM2(2)],'color','r');



