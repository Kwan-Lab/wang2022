function savebehfigpath = bandit_predictSession(dataIndex,exampleLogName,exampleAnimal,model_path,model_type)
% % bandit_predictSession %
%PURPOSE:   Take parameters extracted from model fitting to predict
%           behavior using experimental data
%AUTHORS:   AC Kwan 191212
%
%INPUT ARGUMENTS
%   dataIndex:      a database index table for the sessions to analyze
%   exampleLogName: name of the logfile to compare animal vs model
%   exampleAnimal:  name of the logfile to compare animal vs model
%   model_path:     path containing the .mat file for the models
%   model_type:     the name of the model, e.g. 'Q_RPE' or 'WSLS',
%                   corresponding to the .m file names
%
%OUTPUT ARGUMENTS
%   savebehfigpath: the path related to analysis figure for that session

%% load the animal data
idx = find(ismember(dataIndex.LogFileName,exampleLogName)==1);
load(fullfile(dataIndex.BehPath{idx},[dataIndex.LogFileName{idx}(1:end-4),'_beh.mat']));
   
trials = value_getTrialMasks(trialData);
stats = value_getTrialStats(trials, sessionData.nRules);
stats = value_getTrialStatsMore(stats);

%% use the animal data and fitted model parameters to simulate agent updating

load(fullfile(model_path,[model_type '.mat']));
if exist('session','var')   %if the model fitting was done per session, then find the parameters for that specific session
    idx_sim = find(ismember(session,exampleLogName)==1);
else  %else the model fitting was done per animal
    idx_sim = find(ismember(animal,exampleAnimal)==1);
end

player_fit.label=['algo_' model_type];
player_fit.params=fitpar{idx_sim};
    
stats_sim=predictAgent(player_fit,stats);

%% plot the animal vs simulated behaviors

n_plot = 100*ceil(numel(stats.c)/100); %plot up to the nearest 100 trials
if n_plot > 1000   %for computer simulations, #trials can be too high to plot effectively
    n_plot = 1000;
end

% Create a subfolder to save the images for this session
% folder named using year/month/day of file
yr=num2str(sessionData.dateTime{1}(9:10));
mo=num2str(sessionData.dateTime{1}(1:2));
day=num2str(sessionData.dateTime{1}(4:5));
savebehfigpath = fullfile(dataIndex.BehPath{idx},[yr mo day]);

if ~exist(savebehfigpath,'dir')
    mkdir(savebehfigpath);
end

tlabel=['Animal: ' exampleLogName];
tlabel2=['Model: ' player_fit.label];
plot_session_sim(stats,stats_sim,n_plot,tlabel,tlabel2);

end

