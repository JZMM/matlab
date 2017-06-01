%% This code use for calsium-spike transforming, the noise argurment for  fast_oopsi()
% is std of noise

% std: std of data - smooth(data,5)

load('dff_ZL_h37_20161113_field1_d167_3x_Same_ROI_for_2nd_of_3_consecutive_days');
rawData = dff_prctileF0;
% temp_data = zeros(size(AlineCaSigData.CaSigAlineTone));
for nROI = 1:nROIs
    temp_data(:,nROI,:) = rawData{nROI};   % dff_subtr_preSoundmean   dff_prctileF0 
end                                                       % dff_meanF0               fraw_corr_trials   

%%
FoldStd = 1;
SmoBinSize = 5;

V.est_lam = 1;
V.est_sig = 1;
V.fast_iter_max = 10;

[nTrials,nROIs,nF] = size(temp_data);
V.Ncells = 1;
% V.T = nF; 
% V.T = size(temp_data,1) * size(temp_data,3);
V.Npixels = 1;
V.dt = FrameTime/1000;
% P.lam = 10;
PrcValue = 0.5;     % remove percentile of outliers

nSpikes = zeros(nTrials,nROIs,nF);
ROIstd = zeros(nROIs,1);
% for nROI = 1:nROIs
%     cROIdata = squeeze(temp_data(:,nROI,:)); 
%     Noise{nROI} = [];   
%     Noise{nROI} = reshape(cROIdata',1,[]) - smooth(reshape(cROIdata',1,[]),SmoBinSize)';
%     cStd = 3*std(Noise{nROI});
%     if isnan(cStd)
%         error('Input data cannot contains nan value.');
%     end
%     P.sig = cStd;    
%     ROIstd(nROI) = cStd;
% %     parfor tr = 1:nTrials
% %         nTrace = cROIdata(tr,:);
% %         nsTrace = smooth(nTrace,11,'rloess');
% % %         nsTrace = zscore(nsTrace);
% %         [n_best,~,~,~]=fast_oopsi(nsTrace,V,P);
% % %         n_best(1:4) = 0;
% % %         n_best(n_best < 0.05) = 0;
% %         nSpikes(tr,nROI,:) = n_best;
% %     end  
%     Trace = reshape(cROIdata',1,[]);
%     SmoTrace = smooth(Trace,11,'rloess');
%     [n_best,~,~,~]=fast_oopsi(SmoTrace',V,P);
%     nSpikes(:,nROI,:) = reshape(n_best,fliplr(size(cROIdata)))';
% end
Trace_connect = zeros(nROIs,TrialN*250);
SmoTrace = zeros(nROIs,TrialN*250);
ROIstd = zeros(1,nROIs);
tic
parfor nROI = 1:nROIs

    Trace_connect(nROI,:) = reshape((squeeze(temp_data(:,nROI,:)))',1,[]);
    SmoTrace(nROI,:) = (smooth(Trace_connect(nROI,:),SmoBinSize,'rloess'))';
    cStd = FoldStd*std(Trace_connect(nROI,:) - SmoTrace(nROI,:));
    ROIstd(nROI) = cStd;
    if isnan(cStd)
        error('Input data cannot contains nan value.');
    end 
%     parfor tr = 1:nTrials
%         nTrace = cROIdata(tr,:);
%         nsTrace = smooth(nTrace,11,'rloess');
% %         nsTrace = zscore(nsTrace);
%         [n_best,~,~,~]=fast_oopsi(nsTrace,V,P);
% %         n_best(1:4) = 0;
% %         n_best(n_best < 0.05) = 0;
%         nSpikes(tr,nROI,:) = n_best;
%     end  

end
toc
%
for m = 1:nROIs
    P.sig =  ROIstd(m);
    PP{m} = P;
%     V.Npixels = sum(sum(SavedCaTrials.ROIinfo.ROImask{m}));
    VV{m} = V;
end
%
% tic
% for nROI = 1:nROIs
%     CaSdata = squeeze(temp_data(:,nROI,:));
%     parfor tr = 1:TrialN
%         nTrace = CaSdata(tr,:);
% %         nsTrace = zscore(nsTrace);
%         [n_best,~,~,~]=fast_oopsi(nTrace,VV{nROI},PP{nROI});
% %         n_best(1:4) = 0;
% %         n_best(n_best < 0.05) = 0;
%         nSpikes(tr,nROI,:) = n_best;        
%     end
% end
% toc
%
tic
parfor n = 1:nROIs 
    [n_best,~,~,~]=fast_oopsi(Trace_connect(n,:),VV{n},PP{n});  %use un-smooth data
%     [n_best,~,~,~]=fast_oopsi(SmoTrace(n,:),V,PP{n});  %use smooth data
    nSpikes_connect(n,:) = n_best;
          
end
toc

for nn = 1:nROIs
    nSpikes(:,nn,:) = reshape(nSpikes_connect(nn,:),fliplr([TrialN 250]))'; 
end
%% save
Spikes = nSpikes;
% save h36_20161127_1st_of_3_spikes(1xstdOfNoise_unReshape_unSmomth_withoutPixel_unSubstract_10inter).mat Spikes Miss_Ind Action_choice Probe_Ind Trial_type FrameTime Tone_frequency OnsetF...
%     F_num_T F_num_A MinOnsFrame MinActFrame ; 
save h36_20161127_1st_of_3_spikes(1xstdOfNoise_unReshape_unSmomth_withoutPixel_unSubstract_10inter).mat Spikes Miss_Ind Action_choice Probe_Ind Trial_type FrameTime Tone_frequency OnsetF...
    F_num_T F_num_A MinOnsFrame MinActFrame Trace_connect nSpikes_connect; 
%%  Plot and check
rawData = dff_prctileF0;
[Trials,ROIs,Frames] = size(Spikes);
nROI = 40;
tr = 31;
tempSpike_N = Spikes(tr,nROI,:);
tempCaS = squeeze(rawData{nROI}(tr,:));
close all;
Spike_data1 = squeeze(tempSpike_noi);
Spike_data2 = smooth(tempSpike_noi,5);
figure;hold on;set(gcf,'position',[2000 100 800 400]);
[hax1,plot1,plot2] = plotyy([1:Frames],Spike_data1,[1:Frames],tempCaS);
ylabel(hax1(1),'Spikes');
ylabel(hax1(2),'CaSig trce');
ylim(hax1(1),[0 0.7]);
ylim(hax1(2),[min(tempCaS)-(max(tempCaS)-min(tempCaS))/20 ...
    max(tempCaS)+(max(tempCaS)-min(tempCaS))/20]);
xlim(hax1(1),[1 Frames]);
xlim(hax1(2),[1 Frames]);
title('UnSmooth');
figure;hold on;set(gcf,'position',[2000 500 800 400]);
[hax2,plot3,plot3] = plotyy([1:Frames],Spike_data2,[1:Frames],smooth(tempCaS,5));
ylabel(hax2(1),'Spikes');
ylabel(hax2(2),'CaSig trce');
ylim(hax2(1),[0 0.7]);
ylim(hax2(2),[min(tempCaS)-(max(tempCaS)-min(tempCaS))/20 ...
    max(tempCaS)+(max(tempCaS)-min(tempCaS))/20]);
title('Smooth');
xlim(hax2(1),[1 Frames]);
xlim(hax2(2),[1 Frames]);

figure;hold on;set(gcf,'position',[2800 500 800 400]);
[hax2,plot3,plot3] = plotyy([1:Frames],Spike_data1,[1:Frames],smooth(tempCaS,5));
ylabel(hax2(1),'Spikes');
ylabel(hax2(2),'CaSig trce');
ylim(hax2(1),[0 0.7]);
ylim(hax2(2),[min(tempCaS)-(max(tempCaS)-min(tempCaS))/20 ...
    max(tempCaS)+(max(tempCaS)-min(tempCaS))/20]);
% plot1 = [.7 .8 .8];
% plot2.color = 'r';
%% plot and check
rawData = dff_prctileF0;
% rawData = dff_subtr_preSoundmean;


load('h36_20161127_1st_of_3_spikes(3xstdOfNoise_reshape)');
Spike_3Noi_re = Spikes;
load('h36_20161127_1st_of_3_spikes(2xstdOfNoise)');
Spike_2Noi = Spikes;
load('h36_20161127_1st_of_3_spikes(3xstdOfNoise)');
Spike_3Noi = Spikes;
%%
nROI = 42;
tr = 127;
[Trials,ROIs,Frames] = size(Spike_3Noi_re);
tempSpike_noi = squeeze(Spike_3Noi_re(tr,nROI,:));

tempSpike_2noi = squeeze(Spike_2Noi(tr,nROI,:));


tempSpike_pre = squeeze(Spike_3Noi(tr,nROI,:));
CaS = reshape(rawData{nROI},1,[]);
tempCaS = squeeze(rawData{nROI}(tr,:));

close all;

figure;hold on;set(gcf,'position',[2000 100 800 400]);
[hax1,plot1,plot2] = plotyy([1:Frames],tempSpike_noi,[1:Frames],tempCaS);
ylim(hax1(1),[0 0.7]);
ylim(hax1(2),[min(CaS)-(max(CaS)-min(CaS))/20 ...
    max(CaS)+(max(CaS)-min(CaS))/20]);

figure;hold on;set(gcf,'position',[2000 500 800 400]);
[hax2,plot2,plot3] = plotyy([1:Frames],tempSpike_2noi,[1:Frames],tempCaS);
ylim(hax2(1),[0 0.7]);
ylim(hax2(2),[min(CaS)-(max(CaS)-min(CaS))/20 ...
    max(CaS)+(max(CaS)-min(CaS))/20]);
figure;hold on;set(gcf,'position',[2800 500 800 400]);
[hax3,plot3,plot4] = plotyy([1:Frames],tempSpike_pre,[1:Frames],tempCaS);
ylim(hax3(1),[0 0.7]);
ylim(hax3(2),[min(CaS)-(max(CaS)-min(CaS))/20 ...
    max(CaS)+(max(CaS)-min(CaS))/20]);
%%
ylabel(hax1(1),'Spikes');
ylabel(hax1(2),'CaSig trce');
ylim(hax1(1),[0 0.7]);
ylim(hax1(2),[min(tempCaS)-(max(tempCaS)-min(tempCaS))/20 ...
    max(tempCaS)+(max(tempCaS)-min(tempCaS))/20]);
xlim(hax1(1),[1 Frames]);
xlim(hax1(2),[1 Frames]);
title('UnSmooth');
figure;hold on;set(gcf,'position',[2000 500 800 400]);
[hax2,plot3,plot3] = plotyy([1:Frames],Spike_data2,[1:Frames],smooth(tempCaS,5));
ylabel(hax2(1),'Spikes');
ylabel(hax2(2),'CaSig trce');
ylim(hax2(1),[0 0.7]);
ylim(hax2(2),[min(tempCaS)-(max(tempCaS)-min(tempCaS))/20 ...
    max(tempCaS)+(max(tempCaS)-min(tempCaS))/20]);
title('Smooth');
xlim(hax2(1),[1 Frames]);
xlim(hax2(2),[1 Frames]);

figure;hold on;set(gcf,'position',[2800 500 800 400]);
[hax2,plot3,plot3] = plotyy([1:Frames],Spike_data1,[1:Frames],smooth(tempCaS,5));
ylabel(hax2(1),'Spikes');
ylabel(hax2(2),'CaSig trce');
ylim(hax2(1),[0 0.7]);
ylim(hax2(2),[min(tempCaS)-(max(tempCaS)-min(tempCaS))/20 ...
    max(tempCaS)+(max(tempCaS)-min(tempCaS))/20]);
%% to check cdf and mode of pre-tone
saveName = 'h37 20161119 preTone Corr'; % remember to remove '_'  !!!!!!!!!!
% load(dff_ZL_h36_20161127_fied1_d150_3x);
isOpen  = exportToPPTX();
if ~isempty(isOpen),
    % If PowerPoint already started, then close first and then open a new one
    exportToPPTX('close');
end
if ~exist([saveName '.pptx'])
    exportToPPTX('new','Dimensions',[16 9]);
else
    exportToPPTX('open',[saveName '.pptx']);
end   
for nROI = 1:size(dff_prctileF0,2)
    cROIdata = squeeze(dff_prctileF0{nROI}); 
    tempSpike_noi = [];
    NormTemp = [];
    for tr = 1:size(cROIdata,1)
        temp = cROIdata(tr,[1:round(OnsetT(tr)/FrameTime)]);
        tempSpike_noi = [tempSpike_noi temp];
        NormTemp = [NormTemp temp-mean(temp)];
    end
    

    a= prctile(NormTemp,0.5);
    b= prctile(NormTemp,99.5);
    temp_Norm = NormTemp;
    temp_Norm(temp_Norm>b) = [];
    temp_Norm(temp_Norm<a) = [];
    STD1 = std(tempSpike_noi);
    STD2 = std(NormTemp);
    STD3 = std(temp_Norm);
    subplot(2,2,1);cdfplot(tempSpike_noi);
    text(mean(tempSpike_noi),0.5,['std=' num2str(STD1)]);
    subplot(2,2,2);hist(tempSpike_noi,1000);
    subplot(2,2,3);cdfplot(NormTemp);
    line([a a],[0 1],'color','c');
    line([b b],[0 1],'color','c');
    text(mean(NormTemp),0.5,['std=' num2str(STD2)]);
    text(mean(NormTemp),0.3,['std=' num2str(STD3)]);
    subplot(2,2,4);hist(NormTemp,1000);
    line([a a],[0 100]);
    line([b b],[0 100]);
    suptitle([saveName ' ROI-' num2str(nROI)]);
    exportToPPTX('addslide');
    exportToPPTX('addpicture',gcf);
    close;
end
exportToPPTX('saveandclose',saveName);