%% % This code use 9/10(or ..) of trials to train the model and test the left 1/10 trials for history trials
% And plotting
load('AlineCaSig_h36_20161127_1st_of_3.mat');
load('ZL_h36_Rig2P_20161127-fix_Virables');
TitleN = 'h36_20161127';

TestFrac = 20; % choose propotion of trials are used to test the model(10 means 1/10)

trainIndName = '~Miss_Ind & ~Probe_Ind'; % select trial type to train
% trainIndName = '~Probe_Ind'; % select trial type to train
% trainIndName = 'Action_choice ~= Trial_type & ~Miss_Ind'; % select trial type to train

testIndName = '~Miss_Ind'; % select trial type to test
% testIndName = '~Miss_Ind | Miss_Ind'; % select trial type to test

RawDataName = 'AlineCaSigData.CaSigAlineAct';  % Select aline data   'SpikeAlignedTone'   'SpikeAlignedAct' 'AlineCaSigData.CaSigAlineTone'
                                                                        %  'AlineCaSigData.CaSigAlineAct'
RawLabelName = 'Data_extract.Trial_type'; % Select outcome label   'Data_extract.Action_choice'    'Data_extract.Trial_type'   
TrainWind = [0 2]; % Select which period's activity to train 
% TestWind = [-0.5 0;...
%             0 0.5;...
%             0.5 1.5;...
%             1.5 2.5;...
%             2.5 3.5;...
%             3.5 10]; % Select which period's activity to test 
TestWind = [0 2]; % Select which period's activity to test         
isSpiOrCaS = 0;     % 1: Spike data     0: CaSig data
isBoot = 1;    % 1: bootstrap the data
IsShuffle = 1; % 1: shuffle the label then train  0: not shuffle label
SaveFig = 1; % 1: save 0: don't save
% if IsShuffle
%     RawLabelName = ['Data_extract.' RawLabelName];
% end
niTer = [300 20]; % niTer(1) is interact time;  niTer(2) is n-trials to mean

if strcmp(RawDataName,'SpikeAlignedTone') | strcmp(RawDataName,'AlineCaSigData.CaSigAlineTone')
    VarNote = ['Tone' num2str(TrainWind(1)) '-' num2str(TrainWind(2))]; %##############################
else
    VarNote = ['Act' num2str(TrainWind(1)) '-' num2str(TrainWind(2))];
end

trainInd = eval(trainIndName);
testInd = eval(testIndName); 
trainSeq = find(trainInd==1);
testSeq = find(testInd==1);
RawData = eval(RawDataName);  
if IsShuffle
    RawLabel = Vshuffle(eval(RawLabelName)); 
else
    RawLabel = eval(RawLabelName); 
end

% making directory
CurrPath = pwd;
if ~exist('Decoding_Plot')
    mkdir('Decoding_Plot');
end

if SaveFig
    cd('Decoding_Plot');
    if isSpiOrCaS
        if ~exist('WithSpike')
            mkdir('WithSpike')
        end
        cd('WithSpike');
    else
        if ~exist('WithCaS')
            mkdir('WithCaS')
        end
        cd('WithCaS');    
    end
    if ~exist('CurrentTrialDecoding')
        mkdir('CurrentTrialDecoding')
    end
    cd('CurrentTrialDecoding');
    if ~exist(VarNote)
        mkdir(VarNote)
    end
    cd(VarNote);
    if strcmp(RawLabelName,'Data_extract.Action_choice')
        FoldN = 'ActLable';
        if ~exist(FoldN)
            mkdir(FoldN)
        end
    elseif strcmp(RawLabelName,'Data_extract.Trial_type')
        FoldN = 'ToneLable';
        if ~exist(FoldN)
            mkdir(FoldN)
        end    
    end
    cd(FoldN);
    if IsShuffle
        mkdir('shuffLable')
        cd('shuffLable');
    end
    SavePath = pwd;
    cd(CurrPath)
end

Tone_frequency = Data_extract.Tone_frequency;

if length(TrainWind) < 2
    TimeWindFrame =[round(TrainWind(1)*1000/FrameTime)+1+double(MinOnsFrame) size(RawData,3)];
else
    TimeWindFrame =[round(TrainWind(1)*1000/FrameTime)+1 floor(TrainWind(2)*1000/FrameTime)]+double(MinOnsFrame);
end
if TimeWindFrame(2) > size(RawData,3)
    TimeWindFrame(2) = size(RawData,3);
    warning('The input end time is out of range!!!! ');
end
ShuffInd = Vshuffle(trainSeq);  % shuffle the order of trial sequence
for i = 1:TestFrac
    if i < TestFrac
        testBlock = [1:round(length(trainSeq)/TestFrac)]+round(length(trainSeq)/TestFrac)*(i-1);
    else
        testBlock = [1+round(length(trainSeq)/TestFrac)*(i-1):length(trainSeq)];
    end
   traiInd{i} = setdiff(ShuffInd,ShuffInd(testBlock));
   tesInt{i} = setdiff(testSeq,traiInd{i});
end
% get model
Fre =  unique(Tone_frequency);

for j = 1:TestFrac
    [ mdl{j},CVmodel{j},TrainErro{j} ] = TbTdecoding_v1( RawData,traiInd{j},RawLabel,TimeWindFrame,isSpiOrCaS,isBoot,niTer );
end

% test time window
for ii = 1:size(TestWind,1)
if TestWind(ii,2)*1000/FrameTime > size(RawData,3)
    testFrame =[round(TestWind(ii,1)*1000/FrameTime)+1+double(MinOnsFrame) size(RawData,3)];
    warning('The input end time is out of range!!!! ');
else
    testFrame =[round(TestWind(ii,1)*1000/FrameTime)+1 floor(TestWind(ii,2)*1000/FrameTime)]+double(MinOnsFrame);
end   
%predict
for m = 1:TestFrac
    Tone = Tone_frequency(tesInt{m});
    if isSpiOrCaS       
        TestData = sum(RawData(tesInt{m},:,testFrame(1):testFrame(2)),3);
    else
        TestData = max(RawData(tesInt{m},:,testFrame(1):testFrame(2)),[],3);
%         TestData = mean(RawData(tesInt{m},:,testFrame(1):testFrame(2)),3);
    end
    [yp{m},~] = predict(mdl{m},TestData);
    for F =1:8
        FreqInd = Tone == Fre(F);
        OutCome{ii}(m,F) = sum(yp{m}(FreqInd));
        nTrials{ii}(m,F) = sum(FreqInd);
    end    
end

% Calculate score
PredictScor(ii,:) = sum(OutCome{ii})./sum(nTrials{ii});
PreSEM(ii,:) = std(OutCome{ii}./nTrials{ii})/TestFrac^0.5;

end
PredictResults.mdl = mdl;
PredictResults.yp = yp;
PredictResults.OutCome = OutCome;
PredictResults.nTrials = nTrials;
PredictResults.PredictScor = PredictScor;
PredictResults.PreSEM = PreSEM;
PredictResults.Scores = Scores;

PredictResults.Fre = Fre;
PredictResults.TestFrac = TestFrac;
PredictResults.trainIndName = trainIndName;
PredictResults.testIndName = testIndName;
PredictResults.RawDataName = RawDataName;
PredictResults.RawLabelName = RawLabelName;
PredictResults.TrainWind = TrainWind;
PredictResults.TestWind = TestWind;
PredictResults.isSpiOrCaS = isSpiOrCaS;
PredictResults.isShuff = isBoot;
PredictResults.niTer = niTer;
PredictResults.TitleN = TitleN;


%% save the results
save ([SavePath '\PredictResults.mat'], 'PredictResults');

%%

figure;hold on
axis square;
x = log2(double(PredictResults.Fre));
set(gcf,'color','w','position',[2200 500 560 420]);
set(gca,'fontsize',16,'fontweight','bold','xtick',x,'xticklabel',round(PredictResults.Fre/100)/10,'position',[0.05 0.1738 0.7750 0.7307]);
if IsShuffle
    if isSpiOrCaS
        title([PredictResults.TitleN ' Decoding-Spike(' VarNote 's)Shuffle'],'interpreter','none');
    else
        title([PredictResults.TitleN ' Decoding-CaS(' VarNote 's)Shuffle'],'interpreter','none');
    end    
else    
    if isSpiOrCaS
        title([PredictResults.TitleN ' Decoding-Spike(' VarNote 's)'],'interpreter','none');
    else
        title([PredictResults.TitleN ' Decoding-CaS(' VarNote 's)'],'interpreter','none');
    end
end

xlabel('Frequecy(kHz)');
ylabel('Right Choice Probility');

xlim([x(1)-0.01 x(end)+0.01]);
ylim([-0.01 1.01]);

AllScores = [PredictResults.Scores;PredictResults.PredictScor];
Col = ['k','g','r','b','c','m','y']; 
%%
PlotNum = 2;

b0 = [0; 1;12;4];      % Initial guess of the parameters. Try different values to get good fit!
y = AllScores(PlotNum,:);
if PlotNum>1
    useWeights = 1;
    nTrials = PredictResults.nTrials{PlotNum-1}(1,:);
    nTrials(1) = sum(PredictResults.nTrials{PlotNum-1}(:,1));
    nTrials(end) = sum(PredictResults.nTrials{PlotNum-1}(:,end));
else
    useWeights = 0;
end

% The logistic function
modelfun = @(b,x) (b(1)+ b(2)./(1+exp(-(x - b(3))./b(4))));
% The 4 parameters to fit.
% b(1), the minimum ?P(Choice Right)?; 
% b(2), the maximum ?P(Choice Right)?.
% b(3), the inflection point of the sigmoid; 
% 1/b(4), the slope of the sigmoid; 

opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
if useWeights == 1 
    [b,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,modelfun,b0,'Weights',nTrials);
else
    [b,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,modelfun,b0,opts);
end
X{PlotNum} = linspace(min(x), max(x), 100);
% y1 = modelfun(b,x1);

% Plot confidence intervel
[Y{PlotNum}, delta] = nlpredci(modelfun, X{PlotNum}, b, R, 'Covar',...
    CovB, 'MSE', MSE);
lower = Y{PlotNum} - delta;
upper = Y{PlotNum} + delta;

%
h_data(PlotNum) = scatter(x,y,80,Col(PlotNum),'filled');

h_fit(PlotNum) = plot(X{PlotNum}, Y{PlotNum},'color',Col(PlotNum),'linewidth',3);
if PlotNum>1
    errorbar(x,y,PredictResults.PreSEM(PlotNum-1,:),'linestyle','none','linewidth',1);
end
%%
alpha .5;
% legend(h_fit(2:end),'-0.5-0s','0-0.5s','0.5-1.5s','1.5-2.5s','2.5-3.5s','3.5-end');
legend(h_fit(1:end),'Beh','Pred');
legend('boxoff');
%  
cd(SavePath);

PsychoVar.X = X;
PsychoVar.Y = Y;
save([PredictResults.TitleN ' PsychoVar(' num2str(PredictResults.TrainWind(1)) '-' num2str(PredictResults.TrainWind(2)) 's).mat'],'PsychoVar');
saveas(gcf,[PredictResults.TitleN ' Neural Psycho(' VarNote 's).fig']);
%% for normalize
NarmBlock = 1;
b0 = [0; 1; 13;2];      % Initial guess of the parameters. Try different values to get good fit!
useWeights = 1;

ActualScor = PredictResults.Scores;
NormScore = PredictResults.PredictScor(NarmBlock,:);
for i=1:8
   Scorfix(i)= ActualScor(1)+(NormScore(i)-NormScore(1))/(NormScore(end)-NormScore(1))*(ActualScor(end)-ActualScor(1));
end

nTrials = PredictResults.nTrials{NarmBlock}(1,:);
nTrials(1) = sum(PredictResults.nTrials{NarmBlock}(:,1));
nTrials(end) = sum(PredictResults.nTrials{NarmBlock}(:,end));
x = log2(double(PredictResults.Fre));
modelfun = @(b,x) (b(1)+ b(2)./(1+exp(-(x - b(3))./b(4))));
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
if useWeights == 1 
    [b,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,Scorfix,modelfun,b0,'Weights',nTrials);
else
    [b,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,Scorfix,modelfun,b0,opts);
end
X_nor = linspace(min(x), max(x), 100);
% y1 = modelfun(b,x1);

% Plot confidence intervel
[Y_nor, delta] = nlpredci(modelfun, X_nor, b, R, 'Covar',...
    CovB, 'MSE', MSE);
lower = Y_nor - delta;
upper = Y_nor + delta;

%
h_data(PlotNum) = scatter(x,Scorfix,80,[.8 .8 .8],'filled');

h_fit(PlotNum) = plot(X_nor, Y_nor,'color',[.8 .8 .8],'linewidth',3);



