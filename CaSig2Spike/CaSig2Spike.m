%% This code use for calsium-spike transforming, the noise argurment for  fast_oopsi()
% is std of noise

% std: std of data - smooth(data,5)
fn = 'dff_ZL_h37_20161113_field1_d167_3x_Same_ROI_for_2nd_of_3_consecutive_days';
DataType = 0;  % 0:dff_prctileF0  1:dff_meanF0   2:dff_subtr_preSoundmean  fraw_corr_trials

load(fn);
if DataType==0
    rawData = dff_prctileF0;
elseif DataType==1
    rawData = dff_meanF0;
elseif DataType==2
    rawData = dff_subtr_preSoundmean;
end
if ~isempty(gcp)
    delete(gcp);
end
ppool = parpool('local',6);
% temp_data = zeros(size(AlineCaSigData.CaSigAlineTone));
if size(rawData,3)==1
    for nROI = 1:nROIs
        temp_data(:,nROI,:) = rawData{nROI};   
    end                                                                        
end













