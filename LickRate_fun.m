function [ output_args ] = LickRate_fun( ff, timeWindow , AlignType)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here
%  timeWindow  ms   0: Tone onset
LLickTime = ff.Left_lick_time;
RLickTime = ff.Right_lick_time;
OnsetT = ff.Tone_onset_time;
RewT = ff.Reward_time;
FirstLickT = FirstLickTime_fun(ff);
Act = ff.Action_choice;
FirstLickTime = zeros(1,length(Act));

for i = 1:length(Act)
   if  Act(i) == 1
       temp = RLickTime{i};
       temp(temp<OnsetT(i)) = [];
       FirstLickTime(i) = min(temp)-OnsetT(i);
   elseif Act(i) == 0
       temp = LLickTime{i};
       temp(temp<OnsetT(i)) = [];
       FirstLickTime(i) = min(temp)-OnsetT(i);    
   elseif Act(i) == 2
       FirstLickTime(i) = NaN;
   end
end
end



end

