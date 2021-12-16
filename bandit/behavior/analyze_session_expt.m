function analyze_session_expt(trialData,trials,tlabel,save_path)
% % analyze_session_expt %
%PURPOSE:   Analyze a single session of bandit behavior
%           --- for measurements that are only valid for an experiment
%           --- e.g., response times are not relevant for a computer stimulation
%AUTHORS:   H Atilgan and AC Kwan 191203
%
%INPUT ARGUMENTS
%   trialData:    the 'trialData' structure for the session
%   trials:       the 'trials' structure for the session
%   tlabel:       the title that will appear on some of the output figures
%   save_path:    path for saving the plots
%
%OUTPUT ARGUMENTS
%

if ~exist(save_path,'dir')
    mkdir(save_path);
end

%% plot lick rates
trialType={{'left','reward'},{'left','noreward'},{'right','reward'},{'right','noreward'}};
lick_trType=bandit_get_lickrate_byTrialType(trialData,trials,trialType,edges);

plot_lickrate_byTrialType(lick_trType);
print(gcf,'-dpng',fullfile(save_path,'lickrates_byTrialType'));    
saveas(gcf, fullfile(save_path,'lickrates_byTrialType'), 'fig');

%% plot response times
valLabel='Response time (s)';
trialType={'go','left','right'};
edges=[0:0.02:1];
respTime_trType=get_val_byTrialType(trialData.rt,trials,trialType,edges,valLabel);

plot_val_byTrialType(respTime_trType);
print(gcf,'-dpng',fullfile(save_path,'rt_byTrialType'));    
saveas(gcf, fullfile(save_path,'rt_byTrialType'), 'fig');

%% plot ITI
valLabel='Inter-trial interval (s)';
trialType={'go','reward','noreward'};
edges=[0:0.1:20];
iti_trType=get_val_byTrialType(trialData.iti,trials,trialType,edges,valLabel);

plot_val_byTrialType(iti_trType);
print(gcf,'-dpng',fullfile(save_path,'iti_byTrialType'));    
saveas(gcf, fullfile(save_path,'iti_byTrialType'), 'fig');

end
