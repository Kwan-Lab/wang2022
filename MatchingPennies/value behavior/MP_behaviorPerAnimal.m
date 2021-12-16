function MP_behaviorPerAnimal(dataIndex,save_path)
% % MP_behaviorPerAnimal %
%PURPOSE:   Analyze matching pennies behavior averaged within animals
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
animalFolder = unique(dataIndex.LogFilePath);
for ii = 1:length(animalFolder)
    Ind = strfind(animalFolder{ii},filesep);
    startInd = Ind(end);
    animalList{ii} = animalFolder{ii}(startInd+1:end);
end

disp('-----------------------------------------------------------');
disp(['--- Analyzing - summary of ', int2str(numel(animalList)) ' animals']);
disp('-----------------------------------------------------------');

for j = 1:numel(animalList)
    
    %which session belong to this one animal
    currAnimalSessions = ismember(dataIndex.Animal,animalList(j));
    
    %concatenate the sessions for this one animal
    pennies_summary(dataIndex(currAnimalSessions,:));
end

close all;

