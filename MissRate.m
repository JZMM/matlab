%  This session code calculate the miss rate of manipulation for opto and
%  non opto in same session
files = dir('*.mat');
for i =1:length(files)
    fn = files(i).name;
    load(fn);
    if strcmp(fn(1:2),'ZL')
        AnimName{i} = fn(4:6);
    else
        AnimName{i} = input('Input an animal name:');
    end
    Miss_Ind = Data_extract.Miss_Ind;
    Opto_trial_index = Data_extract.Opto_trial_index;
    OptoMisRate(i) = sum(Opto_trial_index  & Miss_Ind)/sum(Opto_trial_index); 
    UnOptoMisRate(i) = sum(~Opto_trial_index  & Miss_Ind)/sum(~Opto_trial_index);       
end
%%
% cd 'G:\Behavior\ChR2_1\Results\MisRate';
save('VgChR_3rdOpto_MisRate_OptoUnOpto','AnimName','OptoMisRate','UnOptoMisRate');
%%
%  This session code calculate the miss rate of manipulation for whole
%  session
files = dir('*.mat');
for i =1:length(files)
    fn = files(i).name;
    load(fn);
    if strcmp(fn(1:2),'ZL')
        AnimName{i} = fn(4:6);
    else
        AnimName{i} = input('Input an animal name:');
    end
    Miss_Ind = Data_extract.Miss_Ind;
    MisRate(i) = sum(Miss_Ind(1:end-30))/(length(Miss_Ind)-30);       
end
%%
save('VgChR_1stUnOpto_MisRate_session','AnimName','MisRate');


%% this session use to plot opto and non-opto in one session
fn1 = 'VgChR_3rdOpto_MisRate_OptoUnOpto';
fn2 = 'VgChR1_1stOpto_MisRate_OptoUnOpto';
load(fn1);
OptoMisRate1 = OptoMisRate;
UnOptoMisRate1 = UnOptoMisRate;
load(fn2);
OptoMisRate2 = OptoMisRate;
UnOptoMisRate2 = UnOptoMisRate;


Mar = {'o','o'};
Col = {'k',[.5 .5 .5]};
MarSize = 9;
ExtraSize = 2;

figure;hold on;
axis square;
set(gcf,'color','w','position',[2000 200 300 300]);
title('Miss Rate');
xlabel('Non-Opto');
ylabel('Opto');
set(gca,'fontsize',18,'fontweight','bold'); % axis po: [0.3 0.3 0.4 0.6] for 2, [0.3 0.3 0.4 0.6] for 3

xlim([0 .14]);
ylim([0 .14]);
LineP = line([0 1],[0 1],'linestyle','-','color',[.8 .8 .8],'linewidth',1);

PL1 = plot(UnOptoMisRate1,OptoMisRate1,'wo','linewidth',1,'markersize',8,'markerfacecolor','k');
PL2 = plot(UnOptoMisRate2,OptoMisRate2,'wo','linewidth',1,'markersize',8,'markerfacecolor',[.5 .5 .5]);

%% this session use to plot opto and non-opto for different session

Msize = 40;
fn1 = 'C57_1stUnMusci_MisRate_session';
fn2 = 'C57_3rdMusci_MisRate_session';

fn3 = 'C57_2ndCB_MisRate_session';
fn4 = 'C57_1stMusci_MisRate_session';

fn5 = 'VgChR_1stUnOpto_MisRate_session';
fn6 = 'VgChR_3rdOpto_MisRate_session';

fn7 = 'VgChR1_2ndUnOpto_MisRate_session';
fn8 = 'VgChR1_1stOpto_MisRate_session';

fn9 = 'VgChrim_2ndUnOpto_MisRate_session';
fn10 = 'VgChrim_1stOpto_MisRate_session';
load(fn1);
MisRate1 = MisRate;
load(fn2);
MisRate2 = MisRate;
load(fn3);
MisRate3 = MisRate;
load(fn4);
MisRate4 = MisRate;
load(fn5);
MisRate5 = MisRate;
load(fn6);
MisRate6 = MisRate;
load(fn7);
MisRate7 = MisRate;
load(fn8);
MisRate8 = MisRate;
load(fn9);
MisRate9 = MisRate;
load(fn10);
MisRate10 = MisRate;


Mar = {'o','o','d','d','*'};
Col = {'k',[.5 .5 .5],'k',[.5 .5 .5],[.5 .5 .5]};
MarSize = 9;
ExtraSize = 2;

figure;hold on;
axis square;
set(gcf,'color','w','position',[2000 200 300 300]);
title('Miss Rate');
xlabel('Control');
ylabel('Manipulation');
set(gca,'fontsize',18,'fontweight','bold'); % axis po: [0.3 0.3 0.4 0.6] for 2, [0.3 0.3 0.4 0.6] for 3

xlim([0 .2]);
ylim([0 .4]);
LineP = line([0 1],[0 1],'linestyle','-','color',[.8 .8 .8],'linewidth',1);

% PL1 = plot(MisRate1,MisRate2,'w.','linewidth',1,'markersize',Msize,'markerfacecolor','k');
% PL2 = plot(MisRate3,MisRate4,'w.','linewidth',1,'markersize',Msize,'markerfacecolor',[.5 .5 .5]);
% PL3 = plot(MisRate5,MisRate6,'w.','linewidth',1,'markersize',Msize,'markerfacecolor','k');
% PL4 = plot(MisRate7,MisRate8,'w.','linewidth',1,'markersize',Msize,'markerfacecolor',[.5 .5 .5]);
% PL5 = plot(MisRate9,MisRate10,'w.','linewidth',1,'markersize',Msize,'markerfacecolor',[.5 .5 .5]);
PL1 = scatter(MisRate1,MisRate2,Msize,'k','o','filled');
PL2 = scatter(MisRate3,MisRate4,Msize,[.5 .5 .5],'o','filled');
PL3 = scatter(MisRate5,MisRate6,Msize,'b','d','filled');
PL4 = scatter(MisRate7,MisRate8,Msize,'c','d','filled');
PL5 = scatter(MisRate9,MisRate10,Msize,[1 .25 .3],'o','filled');
alpha .5










