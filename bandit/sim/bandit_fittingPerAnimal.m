function bandit_fittingPerAnimal(dataIndex,save_path)
% % bandit_fittingPerAnimal %
%PURPOSE:   Fit various learning algorithms to experimental data, on a
%           per-animal basis by merging sessions from same animal
%AUTHORS:   H Atilgan and AC Kwan 191212
%
%INPUT ARGUMENTS
%   dataIndex:    a database index table for the sessions to analyze
%   save_path:    path for saving the analysis
%
%OUTPUT ARGUMENTS
%

model{1}.name = 'WSLS';   % text label to refer to the model
model{1}.fun = 'funWSLS'; % the corresponding .m code for the model
model{1}.initpar=0.8;     % initial parameter: [prob_WSLS]
model{1}.lb=0.0001;            % upper bound of parameters
model{1}.ub=1;            % lower bound of parameters

model{2}.name = 'Q_RPE';            % text label to refer to the model
model{2}.fun = 'funQ_RPE';      % the corresponding .m code for the model
model{2}.initpar=[0.2 5];       % initial [alpha beta]
model{2}.lb=[0.0001 0.0001];              % upper bound of parameters
model{2}.ub=[1 100];            % lower bound of parameters

model{3}.name = 'DQ_RPE';           % text label to refer to the model
model{3}.fun = 'funDQ_RPE';     % the corresponding .m code for the model
model{3}.initpar=[0.5 5 0.2];   % initial [alpha_reward beta alpha_noreward]
model{3}.lb=[0.0001 0.0001 0.0001];            % upper bound of parameters
model{3}.ub=[1 100 1];          % lower bound of parameters

model{4}.name = 'FQ_RPE';      % text label to refer to the model
model{4}.fun = 'funFQ_RPE';    % the corresponding .m code for the model
model{4}.initpar=[0.5 5];      % initial [alpha_reward beta]
model{4}.lb=[0.0001 0.0001];             % upper bound of parameters
model{4}.ub=[1 100];           % lower bound of parameters

model{5}.name = 'FQ_RPE_CK'; % with a choice autocorrelation term
model{5}.fun = 'funFQ_RPE_CK';
model{5}.initpar = [0.1 0 0.1 0.1]; % initial [alpha beta tau phi]
model{5}.lb = [0 0 0 0];
model{5}.ub = [1 inf 1 inf];
% model{5}.name = 'bayesian_baseline';      % text label to refer to the model
% model{5}.fun = 'funbayesian_baseline';    % the corresponding .m code for the model
% model{5}.initpar=0.05;   % initial [alpha_reward beta lambda]
% model{5}.lb=0.0001;          % upper bound of parameters
% model{5}.ub=1;          % lower bound of parameters
% 
% model{6}.name = 'bayesian_block';      % text label to refer to the model
% model{6}.fun = 'funbayesian_block';    % the corresponding .m code for the model
% model{6}.initpar=0.4;   % initial [alpha_reward beta lambda]
% model{6}.lb=0.0001;          % upper bound of parameters
% model{6}.ub=1;          % lower bound of parameters
% 
% model{7}.name = 'bayesian_twoparams';      % text label to refer to the model
% model{7}.fun = 'funbayesian_twoparams';    % the corresponding .m code for the model
% model{7}.initpar=[0.05 0.05];   % initial [alpha_reward beta lambda]
% model{7}.lb=[0.0001 0.0001];          % upper bound of parameters
% model{7}.ub=[1 1];          % lower bound of parameters

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
    fn = dir(fullfile(save_path, [model{k}.fun(4:end) '.mat']));
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
            disp(['   ' num2str(sum(strcmp(dataIndex.Animal,animalList(j)))) ' sessions associated with this animal']);
            
            %which session belong to this one animal
            currAnimalSessions = ismember(dataIndex.Animal,animalList(j));
            
            %concatenate the sessions for this one animal
            [trialData, trials, nRules] = merge_sessions(dataIndex(currAnimalSessions,:));
            
            %stats = value_getTrialStats(trials, nRules);
            %stats = value_getTrialStatsMore(stats);
            stats_all = bandit_merge_sessions(dataIndex(currAnimalSessions,:));
            
            if isfield(model{k},'lb')
                [fitpar{j}, ~, bic{j}, ~]=fit_fun(stats_all,model{k}.fun,model{k}.initpar,model{k}.lb,model{k}.ub);
            else
                [fitpar{j}, ~, bic{j}, ~]=fit_fun(stats_all,model{k}.fun,model{k}.initpar);
            end
            
             % get the latent variable, save them in sessions
            if strcmp(model{k}.name,'FQ_RPE')
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

