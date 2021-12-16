function MP_pupilMLR_perAnimal(dataIndex,modelpath, savefigpath)
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

%% go through each animal
animalFolder = unique(dataIndex.LogFilePath);
for ii = 1:length(animalFolder)
    animalList{ii} = animalFolder{ii}(29:end);
end


disp('-----------------------------------------------------------');
disp(['--- Analyzing - summary of ', int2str(numel(animalList)) ' animals']);
disp('-----------------------------------------------------------');

for j = 1:numel(animalList)
    
    %which session belong to this one animal
    currAnimalSessions = contains(dataIndex.LogFilePath,animalList(j));
    
    %concatenate the sessions for this one animal
    pennies_pupilMLR_perAnimal(dataIndex(currAnimalSessions,:), modelpath, savefigpath, j);
end

close all;