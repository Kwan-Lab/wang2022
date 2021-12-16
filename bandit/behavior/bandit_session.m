function bandit_session(BehPath,LogFileName,savematpath)
% % bandit_session %
%PURPOSE:   Preparing to analyze a single session of mouse behavior
%AUTHORS:   H Atilgan and AC Kwan 191203
%
%INPUT ARGUMENTS
%   BehPath:        path for the location of the analysis folder containing
%                   the behavioral .mat file
%   LogFileName:    name of the logfile
%
%OUTPUT ARGUMENTS
%

%% load the behavioral data
setup_figprop;
disp('-----------------------------------------------------------');
disp('--- Analyzing a single behavioral session ');
disp('-----------------------------------------------------------');
disp(['Loading ' LogFileName]);

load(fullfile(BehPath, ['bandit_',LogFileName(end-29:end-4),'_beh.mat']));

% Get trial information
trials = value_getTrialMasks(trialData);
stats = value_getTrialStats(trials,sessionData.nRules);
stats = value_getTrialStatsMore(stats);

% What to put as title for some of the figures generated
tlabel=strcat('Subject=',sessionData.subject{1},', Time=',sessionData.dateTime(1),'-',sessionData.dateTime(2));

% Create a subfolder to save the images for this session
% folder named using year/month/day of file
yr=num2str(sessionData.dateTime{1}(9:10));
mo=num2str(sessionData.dateTime{1}(1:2));
day=num2str(sessionData.dateTime{1}(4:5));
savebehfigpath = fullfile(BehPath,[yr mo day]);

    
    time_range=[-2 6];
%    plot_session_beh_vert(trialData,trials,[],tlabel,time_range);

    % plot choice behavior - whole sessions

if ~exist(savebehfigpath)
    mkdir(savebehfigpath)
end
cd(savebehfigpath)
stats.playerlabel{1} = 1;
bandit_plot_session_game(stats,sessionData.nTrials,tlabel);
    
analyze_session(stats,tlabel,savebehfigpath);

%analyze_session_expt(trialData,trials,tlabel,savebehfigpath);

%save([savematpath,'_beh_cut.mat'],...
%             'trialData','sessionData','logfileData','trials',...
%             'lregRCUC_output','lregCRInt_output',...
%             'lick_trType','iti_trType','respTime_trType',...
%             'fitpar','bic','nlike',...
%             'entro',...
%             'stats');

    close all;
    clearvars -except i dirs expData;
end