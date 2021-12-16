function MP_pupilMLR_all(dataIndex,modelpath, savefigpath)
% % bandit_behaviorPerAnimal %
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
if ~exist(savefigpath,'dir')
    mkdir(savefigpath);
end

% animalInd = 0 means all animals
pennies_pupilMLR_perAnimal(dataIndex, modelpath, savefigpath, 0);


close all;