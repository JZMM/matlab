%%  find the correlationship of probe performance and neuron diff-fraction of left and right
%
%
%
% the 1st block is to find the correlation of left and right performace
load('All_beh_from_other_0.5');
block1st_1 = Score_1st_Nor;  % Score_1stBlock  Score_1st_Nor
block2nd_1 = Score_2nd_Nor;  % Score_2ndBlock  Score_2nd_Nor
load('Performance_2P_0.5');
block1st_2 = Score_1st_Nor; 
block2nd_2 = Score_2nd_Nor;    

block1st = [block1st_1;block1st_2];
block2nd = [block2nd_1;block2nd_2];
%
figure;hold on; axis square; xlim([-0.4 0.4]); ylim([-0.4 0.4]);
set(gcf,'position',[2500 200 400 300]);
line([-1 1],[0 0],'linestyle',':','color',[.5 .5 .5]);
line([0 0],[-1 1],'linestyle',':','color',[.5 .5 .5]);
set(gca,'fontsize',15,'fontweight','bold');
a = polyfit(block2nd(:,2)-block1st(:,2),block2nd(:,3)-block1st(:,3),1);
x1 = linspace(-1,1);
y1 = polyval(a,x1);
plot(x1,y1,'linewidth',3,'color','k');
[r,p] = corrcoef(block2nd(:,2)-block1st(:,2),block2nd(:,3)-block1st(:,3))
text(0,0,['r = ' num2str(r(1,2))],'fontsize',12,'fontweight','bold');
text(0,0.1,['p = ' num2str(p(1,2))],'fontsize',12,'fontweight','bold');

% plot(block2nd(:,2)-block1st(:,2),block2nd(:,3)-block1st(:,3),'k.','markersize',25)
% plot(block2nd(:,2)-block1st(:,2),block2nd(:,3)-block1st(:,3),'w.','markersize',6)
plot(block2nd_1(:,2)-block1st_1(:,2),block2nd_1(:,3)-block1st_1(:,3),'k.','markersize',25);
plot(block2nd_1(:,2)-block1st_1(:,2),block2nd_1(:,3)-block1st_1(:,3),'w.','markersize',6);
plot(block2nd_2(:,2)-block1st_2(:,2),block2nd_2(:,3)-block1st_2(:,3),'r.','markersize',25);
plot(block2nd_2(:,2)-block1st_2(:,2),block2nd_2(:,3)-block1st_2(:,3),'w.','markersize',6);
%%  this block is find correlation of neuron fraction with performance
fn1 = 'All_beh_from_other_0.66.mat';
fn2 = 'Performance_2P_0.66.mat';
UseNor = 1;
PlotAll = 0;
Conv = 0;

Name1 = fn2(end-6:end-4);
load(fn1);
if UseNor
block1st_1 = Score_1st_Nor;  % Score_1stBlock  Score_1st_Nor
block2nd_1 = Score_2nd_Nor;  % Score_2ndBlock  Score_2nd_Nor
Name2 = 'Nor_';
else
block1st_1 = Score_1stBlock; 
block2nd_1 = Score_2ndBlock;  
Name2 = 'unNor_';
end
load(fn2);
if UseNor
block1st_2 = Score_1st_Nor; 
block2nd_2 = Score_2nd_Nor;   
else
block1st_2 = Score_1stBlock; 
block2nd_2 = Score_2ndBlock;       
end

block1st = [block1st_1;block1st_2];
block2nd = [block2nd_1;block2nd_2];
[r,p] = corrcoef(block2nd(:,2)-block1st(:,2),block2nd(:,3)-block1st(:,3))
load('AUC_Frac.mat');
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
temp1 = NumFrac_Lprob;  % NumFrac_Lprob  NumFrac_Ltrain
temp2 = NumFrac_Rprob; % NumFrac_Rprob  NumFrac_Rtrain
if UseNor
    Scor1 = Score_1st_Nor; 
    Scor2 = Score_2nd_Nor; 
else
    Scor1 = Score_1stBlock;
    Scor2 = Score_2ndBlock; 
end
FracDiff = temp1 - temp2;
Ind1 = find(FracDiff < median(FracDiff));
Ind2 = find(FracDiff > median(FracDiff));

LefDiff = Scor2(:,2)-Scor1(:,2);
RigDiff_1 = Scor2(:,3)-Scor1(:,3);
if Conv
    RigDiff = RigDiff_1*r(1,2);
else
    RigDiff = -RigDiff_1;
end
aa = repmat([0.1 0.2 0.3],1,10);

[Lef_L_R,~] = sort(LefDiff(Ind1));
[Rig_L_R,~]= sort(RigDiff(Ind1));
[Lef_R_L,~]= sort(LefDiff(Ind2));
[Rig_R_L,~]= sort(RigDiff(Ind2));

figure;hold on;
ylabel('Left-ward Improv');
line([-1 5],[0 0],'linestyle','--','color',[.5 .5 .5]);
set(gcf,'position',[2500 200 400 300]);set(gca,'fontsize',15,'fontweight','bold','xtick',[1 2],'xticklabel',{'Lef<Rig','Lef>Rig'},'position',[0.2 0.2 0.6 0.7])
xlim([0.5 2.5]);
for m = 1:length(Ind1)
    Pl1 = plot(0.8+aa(m),Lef_L_R(m),'ko');
    if PlotAll
    Pl2 = plot(0.8+aa(m),Rig_L_R(m),'rv');
    end
end
for n = 1:length(Ind2)
    plot(1.8+aa(n),Lef_R_L(n),'ko');
    if PlotAll
    plot(1.8+aa(n),Rig_R_L(n),'rv');
    end
end
if PlotAll
    if Conv
        Name3 = '_Conv';
    else
        Name3 = '';
    end
	mu1 = mean([LefDiff(Ind1);RigDiff(Ind1)]);
    sem1 = std([LefDiff(Ind1);RigDiff(Ind1)])/16^0.5;
    mu2 = mean([LefDiff(Ind2);RigDiff(Ind2)]);
    sem2 = std([LefDiff(Ind2);RigDiff(Ind2)])/16^0.5;
else
    Name3 = '';
    mu1 = mean([LefDiff(Ind1)])  ;  
    sem1 = std([LefDiff(Ind1)])/8^0.5;
    mu2 = mean([LefDiff(Ind2)]);
    sem2 = std([LefDiff(Ind2)])/8^0.5;
end
line([0.6 0.7],[mu1 mu1],'linewidth',3,'color','k');
line([0.65 0.65],[mu1-sem1 mu1+sem1],'linewidth',1,'color','k');
line([1.6 1.7],[mu2 mu2],'linewidth',3,'color','k');
line([1.65 1.65],[mu2-sem2 mu2+sem2],'linewidth',1,'color','k');
if PlotAll
[p1,h]=ranksum([LefDiff(Ind1);RigDiff(Ind1)],[LefDiff(Ind2);RigDiff(Ind2)]);
[h,p2]=ttest2([LefDiff(Ind1);RigDiff(Ind1)],[LefDiff(Ind2);RigDiff(Ind2)]);
[a1,b1]=normfit([LefDiff(Ind1);RigDiff(Ind1)]);
[a2,b2]=normfit([LefDiff(Ind2);RigDiff(Ind2)]);
P1 = normcdf([LefDiff(Ind1);RigDiff(Ind1)],a1,b1);
P2 = normcdf([LefDiff(Ind2);RigDiff(Ind2)],a2,b2);
hh = [kstest([LefDiff(Ind1);RigDiff(Ind1)],[[LefDiff(Ind1);RigDiff(Ind1)],P1],0.05) kstest([LefDiff(Ind2);RigDiff(Ind2)],[[LefDiff(Ind2);RigDiff(Ind2)],P2],0.05)];
legend([Pl1 Pl2],{'Lef Perf','Rig Perf'});
else
[p1,h]=ranksum([LefDiff(Ind1)],[LefDiff(Ind2)]);
[h,p2]=ttest2([LefDiff(Ind1)],[LefDiff(Ind2)]);
hh = [kstest([LefDiff(Ind1)]) kstest([LefDiff(Ind2)])];
end
text(1.65,-0.17,['p(rank) = ' num2str(p1)],'fontsize',12,'fontweight','bold');
if sum(hh)==0
text(1.65,-0.15,['p(tt) = ' num2str(p2)],'fontsize',12,'fontweight','bold');
end
title([Name2 Name1 Name3],'interpreter','none');
