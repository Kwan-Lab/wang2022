function bandit_behaviorPerSession(dataIndex,save_path)
% % bandit_behaviorPerSession %
%PURPOSE:   Analyze bandit behavior averaged across sessions
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

%% go through each session

disp('-----------------------------------------------------------');
disp(['--- Analyzing - summary of ', int2str(size(dataIndex,1)) ' sessions']);
disp('-----------------------------------------------------------');

for j = 1:size(dataIndex,1)
    
    load(fullfile(dataIndex.BehPath{j},[dataIndex.LogFileName{j}(1:end-4),'_beh.mat']));
   
    trials = value_getTrialMasks(trialData);
    stats = value_getTrialStats(trials, sessionData.nRules);
    stats = value_getTrialStatsMore(stats);
    
    %% plot choice behavior - around switches left to right
    trials_back=10;  % set number of previous trials
    
    sw_output{j}=choice_switch(stats,trials_back);
    
    %% plot choice behavior - around switch high-probability side to low-probability side
    sw_hrside_output{j}=choice_switch_hrside(stats,trials_back);

    %% plot choice behavior - around switches left to right, as a function of the statistics of the block preceding the switch
    L1_ranges=[10 20;10 20;10 20;10 20]; %consider only subset of blocks within the range, for trials to criterion
    L2_ranges=[0 4;5 9;10 14;15 30];      %consider only subset of blocks within the range, for random added number of trials
    sw_random_output{j}=choice_switch_random(stats,trials_back,L1_ranges,L2_ranges);

    sw_hrside_random_output{j}=choice_switch_hrside_random(stats,trials_back,L1_ranges,L2_ranges);

    %% plot tendency to predict upcoming reversal
    %store the x and y variables to be made into a scatter plot
    sw_randomVsChoice{j}.dat=[stats.blockTrialRandomAdded stats.blockPreSwitchWorseChoiceAtSwitch]; 
    sw_randomVsChoice{j}.range{1}=[0 30];
    sw_randomVsChoice{j}.range{2}=[0 0.4];
    sw_randomVsChoice{j}.label{1}={'L_2 (random number of trials added)'};
    sw_randomVsChoice{j}.label{2}={'Fraction of trials';'selecting initial worse option'};
    
    %% plot lick rates
    trialType={{'left','reward'},{'left','noreward'},{'right','reward'},{'right','noreward'}};
    edges=[-0.5:0.02:3];   % edges to plot the lick rate histogram
    lick_trType{j}=get_lickrate_byTrialType(trialData,trials,trialType,edges);
    
    %% plot response times
    valLabel='Response time (s)';
    trialType={'go','left','right'};
    edges=[0:0.01:1];
    respTime_trType{j}=get_val_byTrialType(trialData.rt,trials,trialType,edges,valLabel);
    
    %% plot ITI
    valLabel='Inter-trial interval (s)';
    trialType={'go','reward','noreward'};
    edges=[0:0.1:30];
    iti_trType{j}=get_val_byTrialType(trialData.iti,trials,trialType,edges,valLabel);

    %% plot logistic regression analysis - rewarded/unrewarded choices
    num_regressor = 10;
    [lreg{j}, ~, ~, ~]=logreg_RCUC(stats,num_regressor);

    %% plot logistic regression analysis - rewarded/unrewarded left/right choices
    % this analysis is ill-conditioned -- too many regressors for too few trials in some sessions
%    num_regressor = 10;
%    [lreg_LR{j}, ~, ~, ~]=logreg_RCUC_LR(stats,num_regressor);
end

%% plotting
tlabel = ['Summary of ' int2str(size(dataIndex,1)) ' sessions'];

plot_switch(sw_output,tlabel,stats.rule_labels);
print(gcf,'-dpng',fullfile(save_path,'switches'));
saveas(gcf, fullfile(save_path,'switches'), 'fig');

plot_switch_hrside(sw_hrside_output,tlabel);
print(gcf,'-dpng',fullfile(save_path,'switches_hrside'));
saveas(gcf, fullfile(save_path,'switches_hrside'), 'fig');

plot_switch_random(sw_random_output,tlabel,stats.rule_labels);
print(gcf,'-dpng',fullfile(save_path,'switches_random'));
saveas(gcf, fullfile(save_path,'switches_random'), 'fig');
plot_switch_hrside_random(sw_hrside_random_output,tlabel);
print(gcf,'-dpng',fullfile(save_path,'switches_hrside_random'));
saveas(gcf, fullfile(save_path,'switches_hrside_random'), 'fig');
plot_switch_random_distn(sw_random_output,tlabel);
print(gcf,'-dpng',fullfile(save_path,'switches_random_distn'));
saveas(gcf, fullfile(save_path,'switches_random_distn'), 'fig');

plot_binxaveragey(sw_randomVsChoice,tlabel);
print(gcf,'-dpng',fullfile(save_path,'switches_randomVsChoice'));
saveas(gcf, fullfile(save_path,'switches_randomVsChoice'), 'fig');

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

%plot_logreg(lreg_LR,tlabel);
%print(gcf,'-dpng',fullfile(save_path,'logreg_LR'));    
%saveas(gcf, fullfile(save_path,'logreg_LR'), 'fig');
