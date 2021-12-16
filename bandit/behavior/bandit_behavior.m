function bandit_behavior(dataIndex,save_path)
% % bandit_behavior %
%PURPOSE:   Analyze bandit behavior averaged across animals
%AUTHORS:   H Atilgan and AC Kwan 191204
%
%INPUT ARGUMENTS
%   dataIndex:    a database index table for the sessions to analyze
%   save_path:    path for saving the plots
%
%OUTPUT ARGUMENTS
%

%%
if ~exist(save_path,'dir')
    mkdir(save_path);
end

%% go through each animal
animalList = unique(dataIndex.Animal);

disp('-----------------------------------------------------------');
disp(['--- Analyzing - summary of ', int2str(numel(animalList)) ' animals']);
disp('-----------------------------------------------------------');

for j = 1:numel(animalList)
    
    %which session belong to this one animal
    currAnimalSessions = ismember(dataIndex.Animal,animalList(j));
    
    %concatenate the sessions for this one animal
    [trialData, trials, nRules] = merge_sessions(dataIndex(currAnimalSessions,:));
    stats = value_getTrialStats(trials, nRules);
    stats = value_getTrialStatsMore(stats);
    
    %% plot choice behavior - around switches left to right
    trials_back=10;  % set number of previous trials
    
    sw_output{j}=choice_switch(stats,trials_back);
    
    %% plot choice behavior - around switch high-probability side to low-probability side
    sw_hrside_output{j}=choice_switch_hrside(stats,trials_back);

    %% plot choice behavior - around switches left to right, as a function of the statistics of the block preceding the switch
    L1_ranges=[10 20;10 20;10 20;10 20]; %consider only subset of blocks within the range, for trials to criterion
    L2_ranges=[1 3;4 7;8 15;16 30];      %consider only subset of blocks within the range, for random added number of trials
    sw_random_output{j}=choice_switch_random(stats,trials_back,L1_ranges,L2_ranges);

    %% plot lick rates
    trialType={'reward','noreward'};
    edges=[-2:0.02:5];   % edges to plot the lick rate histogram
    lick_trType{j}=get_lickrate_byTrialType(trialData,trials,trialType,edges);
    
    %% plot response times
    valLabel='Response time (s)';
    trialType={'go','left','right'};
    edges=[-0.2:0.02:2];
    respTime_trType{j}=get_val_byTrialType(trialData.rt,trials,trialType,edges,valLabel);
    
    %% plot ITI
    valLabel='Inter-trial interval (s)';
    trialType={'go','reward','noreward'};
    edges=[0:0.1:30];
    iti_trType{j}=get_val_byTrialType(trialData.iti,trials,trialType,edges,valLabel);
    
    %% plot logistic regression analysis
    num_regressor = 10;
    [lreg{j}, ~, ~, ~]=logreg_RCUC(stats,num_regressor);

end

tlabel = ['Summary of ' int2str(numel(animalList)) ' animals'];

plot_switch(sw_output,tlabel,stats.rule_labels);
print(gcf,'-dpng',fullfile(save_path,'switches'));
saveas(gcf, fullfile(save_path,'switches'), 'fig');

plot_switch_hrside(sw_hrside_output,tlabel);
print(gcf,'-dpng',fullfile(save_path,'switches_hrside'));
saveas(gcf, fullfile(save_path,'switches_hrside'), 'fig');

plot_switch_random(sw_random_output,tlabel,stats.rule_labels);
print(gcf,'-dpng',fullfile(save_path,'switches_random'));
saveas(gcf, fullfile(save_path,'switches_random'), 'fig');

plot_lickrate_byTrialType(lick_trType);
print(gcf,'-dpng',fullfile(save_path,'lickrates_byTrialType'));
saveas(gcf, fullfile(save_path,'lickrates_byTrialType'), 'fig');

plot_val_byTrialType(respTime_trType);
print(gcf,'-dpng',fullfile(save_path,'rt_byTrialType'));
saveas(gcf, fullfile(save_path,'rt_byTrialType'), 'fig');

plot_val_byTrialType(iti_trType);
print(gcf,'-dpng',fullfile(save_path,'iti_byTrialType'));
saveas(gcf, fullfile(save_path,'iti_byTrialType'), 'fig');

plot_logreg(lreg,tlabel);
print(gcf,'-dpng',fullfile(save_path,'logreg'));    
saveas(gcf, fullfile(save_path,'logreg'), 'fig');
