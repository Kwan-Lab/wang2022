function matchingPennies_behaviorAll(dataIndex,save_path)
% % bandit_behaviorPerAnimal %
%PURPOSE:   Analyze bandit behavior averaged across animals
%AUTHORS:   H Atilgan and AC Kwan 191204
%MODIFIED: H WAng
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
cd(save_path);
%% go through each animal

    
    %concatenate the sessions for this one animal
pennies_summary(dataIndex);


close all;
%plot_logreg(lreg_LR,tlabel);
%print(gcf,'-dpng',fullfile(save_path,'logreg_LR'));    
%saveas(gcf, fullfile(save_path,'logreg_LR'), 'fig');
