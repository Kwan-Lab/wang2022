function analyze_session(stats,tlabel,save_path)
% % analyze_session %
%PURPOSE:   Analyze a single session of bandit behavior
%AUTHORS:   H Atilgan and AC Kwan 191203
%
%INPUT ARGUMENTS
%   stats:        the 'stats' structure for the session
%   tlabel:       the title that will appear on some of the output figures
%   save_path:    path for saving the plots
%
%OUTPUT ARGUMENTS
%

if ~exist(save_path,'dir')
    mkdir(save_path);
end

%% analysis of behavioral performance
% plot choice behavior - whole sessions
n_plot = 100*ceil(numel(stats.c)/100); %plot up to the nearest 100 trials
if n_plot > 1000   %for computer simulations, #trials can be too high to plot effectively
    n_plot = 1000;
end

plot_session_task(stats,n_plot,tlabel);
print(gcf,'-dpng',fullfile(save_path,'session'));
saveas(gcf, fullfile(save_path,'session'), 'fig');

%% plot choice behavior - around switches left to right
trials_back=10;  % set number of previous trials
sw_output=choice_switch(stats,trials_back);

plot_switch(sw_output,tlabel,stats.rule_labels);
print(gcf,'-dpng',fullfile(save_path,'switches'));    
saveas(gcf, fullfile(save_path,'switches'), 'fig');

%% plot choice behaviour - around switch high-probability side to low-probability side
sw_hrside_output=choice_switch_hrside(stats,trials_back);

plot_switch_hrside(sw_hrside_output,tlabel);
print(gcf,'-dpng',fullfile(save_path,'switches_hrside'));    
saveas(gcf, fullfile(save_path,'switches_hrside'), 'fig');

%% plot choice behavior - around switches left to right, as a function of the statistics of the block preceding the switch
L1_ranges=[10 20;10 20;10 20;10 20]; %consider only subset of blocks within the range, for trials to criterion
L2_ranges=[0 4;5 9;10 14;15 30];      %consider only subset of blocks within the range, for random added number of trials
sw_random_output=choice_switch_random(stats,trials_back,L1_ranges,L2_ranges);

plot_switch_random(sw_random_output,tlabel,stats.rule_labels);
print(gcf,'-dpng',fullfile(save_path,'switches_random'));
saveas(gcf, fullfile(save_path,'switches_random'), 'fig');

%% plot tendency to predict upcoming reversal
%store the x and y variables to be made into a scatter plot
sw_randomVsChoice.dat=[stats.blockTrialRandomAdded stats.blockPreSwitchWorseChoiceAtSwitch];
sw_randomVsChoice.range{1}=[0 30];
sw_randomVsChoice.range{2}=[0 0.4];
sw_randomVsChoice.label{1}={'L_2 (random number of trials added)'};
sw_randomVsChoice.label{2}={'Fraction of trials';'selecting initial worse option'};

plot_binxaveragey(sw_randomVsChoice,tlabel);
print(gcf,'-dpng',fullfile(save_path,'switches_randomVsChoice'));
saveas(gcf, fullfile(save_path,'switches_randomVsChoice'), 'fig');

%% plot logistic regression analysis
% num_regressor = 10;
% [lreg, ~, ~, ~]=bandit_logreg_RCUC(stats,num_regressor);
% plot_logreg(lreg,tlabel);
% print(gcf,'-dpng',fullfile(save_path,'logreg'));    
% saveas(gcf, fullfile(save_path,'logreg'), 'fig');

end
