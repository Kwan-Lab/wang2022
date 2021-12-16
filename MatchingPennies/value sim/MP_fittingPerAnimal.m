function MP_fittingPerAnimal(dataIndex,save_path)
% %MP_fittingPerAnimal %
%PURPOSE:   Fit various learning algorithms to experimental data, on a
%           per-animal basis by merging sessions from same animal
%AUTHORS:   H Atilgan and AC Kwan 191212
%MODIFIED:  H Wang

%INPUT ARGUMENTS
%   dataIndex:    a database index table for the sessions to analyze
%   save_path:    path for saving the analysis
%
%OUTPUT ARGUMENTS
% some subject has multiple training periods, fit them separately

model{1}.name = 'WSLS';   % text label to refer to the model
model{1}.fun = 'funWSLS'; % the corresponding .m code for the model
model{1}.initpar=0.8;     % initial parameter: [prob_WSLS]
model{1}.lb=0;            % upper bound of parameters
model{1}.ub=1;            % lower bound of parameters

model{2}.name = 'Q_RPE';            % text label to refer to the model
model{2}.fun = 'funQ_RPE';      % the corresponding .m code for the model
model{2}.initpar=[0.2 5];       % initial [alpha beta]
model{2}.lb=[0 0];              % upper bound of parameters
model{2}.ub=[1 inf];            % lower bound of parameters

model{3}.name = 'DQ_RPE';           % text label to refer to the model
model{3}.fun = 'funDQ_RPE';     % the corresponding .m code for the model
model{3}.initpar=[0.5 5 0.2];   % initial [alpha_reward beta alpha_noreward]
model{3}.lb=[0 0 0];            % upper bound of parameters
model{3}.ub=[1 inf 1];          % lower bound of parameters

model{4}.name = 'FQ_RPE';      % text label to refer to the model
model{4}.fun = 'funFQ_RPE';    % the corresponding .m code for the model
model{4}.initpar=[0.5 5];      % initial [alpha_reward beta]
model{4}.lb=[0 0];             % upper bound of parameters
model{4}.ub=[1 inf];           % lower bound of parameters

model{5}.name = 'FQ_RPE_CK'; % with a choice autocorrelation term
model{5}.fun = 'funFQ_RPE_CK';
model{5}.initpar = [0.1 1 0 0]; % initial [alpha beta alpha_K beta_K]
model{5}.lb = [0 0 0 0];
model{5}.ub = [1 inf 1 inf];


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
disp(['--- Fitting models - summary of ', int2str(numel(animalList)) ' animals']);
disp('-----------------------------------------------------------');

for k=1:numel(model)
    
    disp(['Considering model "' model{k}.name '"']);

    % Does the analysis file exist?
    fn = dir(fullfile(save_path, [model{k}.fun(4:end) '_model.mat']));
    if size(fn,1)>0
        answer = questdlg(['There is already .mat file for model fitting for ' model{k}.fun '. Run the analysis again? (may take time)'], ...
            'Run model fitting?', ...
            'Yes','No','Yes');
        if strcmp(answer,'Yes')
            runFit = true;
        else
            runFit = false;
        end
    else
        runFit = true;
    end
    
    if (runFit)
        for j = 1:numel(animalList)
            
            disp(['Processing animal # ' int2str(j) '...']);
            disp(['   ' num2str(sum(contains(dataIndex.LogFilePath,animalList(j)))) ' sessions associated with this animal']);
            
            %which session belong to this one animal
            currAnimalSessions = contains(dataIndex.LogFilePath,animalList(j));
            
            %concatenate the sessions for this one animal
            stats_all = MP_merge_sessions(dataIndex(currAnimalSessions,:));
            
            stats_fit.c = stats_all.c(:,1);
            stats_fit.r =stats_all.r;
            if isfield(model{k},'lb')
                [fitpar{j}, ~, bic{j}, ~]=fit_fun(stats_fit,model{k}.fun,model{k}.initpar,model{k}.lb,model{k}.ub);
            else
                [fitpar{j}, ~, bic{j}, ~]=fit_fun(stats_fit,model{k}.fun,model{k}.initpar);
            end
            
            %get the latent variable, save them in sessions
            if strcmp(model{k}.name, 'FQ_RPE')
                player1.label=['algo_',model{k}.name];   
                player1.params.a=fitpar{j}(1);
                player1.params.b=fitpar{j}(2);

                
                stats_sim=predictAgent(player1,stats_all);
                stats_sim.sLength = stats_all.sessionLength;
                stats_sim.alpha = fitpar{j}(1);
                stats_sim.beta = fitpar{j}(2);
                % save stats_sim
                save(fullfile(save_path,[animalList{j},'_FQRPE.mat']),'stats_sim');
            elseif strcmp(model{k}.name, 'FQ_RPE_CK')
                player1.label=['algo_',model{k}.name];   
                player1.params.a=fitpar{j}(1);
                player1.params.b=fitpar{j}(2);
                player1.params.ac = fitpar{j}(3);
                player1.params.bc = fitpar{j}(4);
                
                stats_sim=predictAgent(player1,stats_all);
                stats_sim.sLength = stats_all.sessionLength;
                stats_sim.alpha = fitpar{j}(1);
                stats_sim.beta = fitpar{j}(2);
                stats_sim.alphac = fitpar{j}(3);
                stats_sim.betac = fitpar{j}(4);
                
                % save stats_sim
                save(fullfile(save_path,[animalList{j},'_FQRPECK.mat']),'stats_sim');
            end
            animal{j} = animalList{j};
        end
        
        %save behavioral .mat file
        save(fullfile(save_path, [model{k}.fun(4:end) '_model.mat']),...
            'animal','fitpar','bic');

        fitparMat = cell2mat(fitpar');
        disp('Median fitted parameters:');
        nanmedian(fitparMat,1)
    end
    
end

end

