function [ FirstLickTime ] = FirstLickTime_fun( ff )
% UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here
LLickTime = ff.Left_lick_time;
RLickTime = ff.Right_lick_time;
OnsetT = ff.Tone_onset_time;
Act = ff.Action_choice;
FirstLickTime = zeros(1,length(Act));
% Test = zeros(1,length(Act));
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
%    if FirstLickTime(i)>3000
%        Test(i) = 1;
%    end
end
end
