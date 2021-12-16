function MP_pupilRL_MLR_all(dataIndex, savefigpath);


%%
if ~exist(savefigpath,'dir')
    mkdir(savefigpath);
end
%% go through each animal
animalFolder = unique(dataIndex.LogFilePath);
for ii = 1:length(animalFolder)
    animalList{ii} = animalFolder{ii}(29:end);
end

for j = 1:numel(animalList)
    
    %which session belong to this one animal
    currAnimalSessions = contains(dataIndex.LogFilePath,animalList(j));
    savefig = fullfile(savefigpath(1:end-18),animalList{j});
% animalInd = 0 means all animals
    if ~exist(savefig,'dir')
        mkdir(savefig);
    end
    MP_pupilRL_acrossSessions(dataIndex(currAnimalSessions,:),savefig);
    %pennies_pupilRL_perAnimal(dataIndex, modelpath, savefigpath, j);
end


close all;
end